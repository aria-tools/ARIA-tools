# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & Brett Buzzanga & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import glob
import logging
import shutil
import osgeo
import numpy as np
import matplotlib.pyplot as plt

import ARIAtools.util.shp
import ARIAtools.util.vrt

LOGGER = logging.getLogger(__name__)

def prep_mask(
        product_dict, maskfilename, bbox_file, prods_TOTbbox, proj,
        amp_thresh=None, arrres=None, workdir='./', outputFormat='ENVI',
        num_threads='2', multilooking=None, rankedResampling=False):
    """
    Function to load and export mask file.
    If "Download" flag, GSHHS water mask will be donwloaded on the fly.
    If full resolution NLCD landcover data is given (NLCD...img) it gets cropped
    """
    _world_watermask = [f' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/'
                        f'gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L{i}.shp'
                        for i in range(1, 5)]

    _world_watermask.append(' /vsizip/vsicurl/https://osmdata.openstreetmap.de'
                            '/download/land-polygons-complete-4326.zip/'
                            'land-polygons-complete-4326/land_polygons.shp')

    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir, 'mask')
    os.makedirs(workdir, exist_ok=True)

    # Get bounds of user bbox_file
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds

    # File must be physically extracted, cannot proceed with VRT format.
    # Defaulting to ENVI format.
    if outputFormat == 'VRT':
        outputFormat = 'ENVI'

    # Download mask
    if maskfilename.lower() == 'download':
        LOGGER.info("Downloading water mask...")
        maskfilename = os.path.join(workdir, 'watermask' + '.msk')
        os.environ['CPL_ZIP_ENCODING'] = 'UTF-8'

        # Make coastlines/islands union VRT
        ARIAtools.util.vrt.renderOGRVRT(
            os.path.join(workdir, 'watermsk_shorelines.vrt'),
            _world_watermask[::2])

        # Make lakes/ponds union VRT
        ARIAtools.util.vrt.renderOGRVRT(
            os.path.join(workdir, 'watermsk_lakes.vrt'),
            _world_watermask[1::2])

        # Initiate water-mask with coastlines/islands union VRT
        # save uncropped mask
        osgeo.gdal.Rasterize(
            os.path.join(workdir, 'watermask_uncropped.msk'),
            os.path.join(workdir, 'watermsk_shorelines.vrt'),
           format=outputFormat, outputBounds=bounds,
           outputType=osgeo.gdal.GDT_Byte, xRes=arrres[0], yRes=arrres[1],
           targetAlignedPixels=True, burnValues=[1], layers='merged')

        osgeo.gdal.Translate(
            os.path.join(workdir, 'watermask_uncropped.msk.vrt'),
            os.path.join(workdir, 'watermask_uncropped.msk'), format="VRT")

        # save cropped mask
        osgeo.gdal.Warp(
            maskfilename, os.path.join(workdir, 'watermask_uncropped.msk.vrt'),
            format=outputFormat, outputBounds=bounds,
            outputType=osgeo.gdal.GDT_Byte, xRes=arrres[0], yRes=arrres[1],
            targetAlignedPixels=True, multithread=True,
            options=['NUM_THREADS=%s' % (num_threads)])

        update_file = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
        update_file.SetProjection(proj)
        update_file.GetRasterBand(1).SetNoDataValue(0.)
        osgeo.gdal.Translate(f'{maskfilename}.vrt', maskfilename, format='VRT')

        # Take inverse of lakes/ponds union because of opposite designation
        # (1 for water, 0 for land) as desired (0 for water, 1 for land)
        lake_masks = osgeo.gdal.Rasterize(
            '', os.path.join(workdir, 'watermsk_lakes.vrt'), format='MEM',
            outputBounds=bounds, outputType=osgeo.gdal.GDT_Byte,
            xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
            burnValues=[1], layers='merged', inverse=True)

        lake_masks.SetProjection(proj)
        lake_masks = lake_masks.ReadAsArray()

        # Update water-mask with lakes/ponds union
        mask_file = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
        mask_array = lake_masks * osgeo.gdal.Open(maskfilename).ReadAsArray()
        mask_file.GetRasterBand(1).WriteArray(mask_array)

    # Use NLCD Mask
    elif os.path.basename(maskfilename).lower().startswith('nlcd'):
        LOGGER.info("Accessing and cropping the NLCD mask.")
        maskfilename = NLCDMasker(os.path.dirname(workdir), product_dict)(
            proj, bounds, arrres, outputFormat)

    # User specified mask
    else:
        # Path to local version of user specified mask
        user_mask = maskfilename  # for clarity
        user_mask_n = os.path.basename(os.path.splitext(user_mask)[0])
        local_mask = os.path.join(workdir, f'{user_mask_n}.msk')
        local_mask_unc = os.path.join(workdir, f'{user_mask_n}_uncropped.msk')
        if user_mask == local_mask:
            LOGGER.warning(
                'The mask you specified already exists in %s, '
                'using the existing one...' % os.path.dirname(local_mask))
            ds = osgeo.gdal.Open(local_mask)

        else:
            # move the mask to the local directory and built a VRT for it
            osgeo.gdal.UseExceptions()

            # shutil.copy(user_mask, local_mask_unc)
            ds = osgeo.gdal.BuildVRT(
                f'{local_mask_unc}.vrt', user_mask, outputBounds=bounds)
            assert ds is not None, f'GDAL could not open user mask: {user_mask}'

        # crop the user mask and write
        osgeo.gdal.Warp(
            f'{local_mask}', ds, format=outputFormat,
            cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrres[0],
            yRes=arrres[1], targetAlignedPixels=True, multithread=True,
            options=[f'NUM_THREADS={num_threads}'])

        # set projection of the local mask
        mask_file = osgeo.gdal.Open(local_mask, osgeo.gdal.GA_Update)
        mask_file.SetProjection(proj)

        # create vrt for local cropped mask
        osgeo.gdal.Translate(f'{local_mask}.vrt', local_mask, format='VRT')

        # assign new local mask for amp thresh
        maskfilename = local_mask

    # Make average amplitude mask
    if amp_thresh is not None:
        amp_file = ARIAtools.util.vrt.rasterAverage(
            os.path.join(workdir, 'avgamplitude'), product_dict, bounds,
            prods_TOTbbox, arrres, outputFormat=outputFormat, thresh=amp_thresh)

        # Update mask with average amplitude
        mask_file = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
        mask_arr = mask_file.ReadAsArray()
        mask_file.GetRasterBand(1).WriteArray(mask_arr * amp_file)

    # crop/expand mask to DEM size?
    osgeo.gdal.Warp(
        maskfilename, maskfilename, format=outputFormat,
        cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrres[0],
        yRes=arrres[1], targetAlignedPixels=True, multithread=True,
        options=[f'NUM_THREADS={num_threads} -overwrite'])

    mask = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
    mask.SetProjection(proj)
    mask.SetDescription(maskfilename)
    mask_array = mask.ReadAsArray()
    mask_array[mask_array != 1] = 0
    mask.GetRasterBand(1).WriteArray(mask_array)

    # Update VRT
    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    osgeo.gdal.Translate(
        maskfilename + '.vrt', maskfilename, options=translate_options)

    # Apply multilooking, if specified
    if multilooking is not None:
        ARIAtools.util.vrt.resampleRaster(
            maskfilename, multilooking, bounds, prods_TOTbbox, rankedResampling,
            outputFormat=outputFormat, num_threads=num_threads)

    # return filename of mask
    return maskfilename


def make_mask(ds_crop, lc):
    # values outside crop; 0 is really just extra check
    lc.extend([0, 255])

    arr = ds_cros.path.ReadAsArray()
    for lclass in lc:
        arr = np.where(arr == lclass, np.nan, arr)

    arr = np.where(np.isnan(arr), np.nan, 1)
    ds = arr2ds(ds_crop, arr)
    return ds


def resamp(src, proj, bounds, arrres, view=False):
    """
    Resample a dataset from src dims using result of merged_productbbox
    """
    if isinstance(src, str) and os.path.exists(src):
        src = osgeo.gdal.Open(src)

    # TODO this doesn't make any sense, path is not a variable here
    path = path if os.path.exists(src.GetDescription()) else ''

    # compute geotranform
    pix_width, pix_height = arrres
    pix_width *= np.sign(bounds[2] - bounds[0])
    pix_height *= np.sign(bounds[1] - bounds[3])
    width = int((bounds[2] - bounds[0]) / pix_width)
    height = int((bounds[1] - bounds[3]) / pix_height)
    gt = [bounds[0], pix_width, 0, bounds[3], 0, pix_height]

    if path is not None:
        dst = osgeo.gdal.GetDriverByName('ENVI').Create(path, width, height, 1,
                                                  osgeo.gdalconst.GDT_Float32)

    else:
        LOGGER.info('width, height', width, height)
        dst = osgeo.gdal.GetDriverByName('MEM').Create(
            '', width, height, 1, osgeo.gdalconst.GDT_Int16)

    dst.SetGeoTransform(gt)
    dst.SetProjection(proj)
    osgeo.gdal.ReprojectImage(
        src, dst, src.GetProjection(), proj,
        osgeo.gdalconst.GRA_NearestNeighbour)
    return dst


def crop_ds(path_raster, path_poly, path_write=None):
    """
    Crop path_raster (or dataset opened with gdal) to a polygon on disk
    """
    assert os.path.exists(path_raster) or path_raster.startswith(
        '/vsi'), 'Invalid path to NLCD mask'
    path_raster = osgeo.gdal.BuildVRT('', path_raster)
    ds = osgeo.gdal.Warp(
        path_write, path_raster, format='VRT', cutlineDSName=path_poly,
        cropToCutline=True, dstNodata=0)

    if path_write is not None:
        ds_vrt = osgeo.gdal.BuildVRT('{}.vrt'.format(path_write), ds)
    return ds


def arr2ds(ds_orig, arr, noData=0):
    """
    Create a new dataset identical to ds_orig except with values of 'arr'
    Performed in memory
    """
    if isinstance(ds_orig, str):
        ds_orig = osgeo.gdal.Open(ds_orig)

    arr = np.where(np.isnan(arr), noData, arr)

    if arr.ndim == 3:
        ds_out = osgeo.gdal.GetDriverByName('MEM').Create(
            '', ds_orig.RasterXSize, ds_orig.RasterYSize, arr.shape[0],
            osgeo.gdal.GDT_Float32)

        for i in range(arr.shape[0]):
            ds_out.GetRasterBand(i + 1).WriteArray(arr[i, :, :])
            ds_out.GetRasterBand(i + 1).SetNoDataValue(noData)
    else:
        ds_out = osgeo.gdal.GetDriverByName('MEM').Create(
            '', ds_orig.RasterXSize, ds_orig.RasterYSize, 1,
            osgeo.gdal.GDT_Float32)
        ds_out.GetRasterBand(1).WriteArray(arr)
        ds_out.GetRasterBand(1).SetNoDataValue(noData)

    ds_out.SetGeoTransform(ds_orig.GetGeoTransform())
    ds_out.SetProjection(ds_orig.GetProjection())
    return ds_out


class NLCDMasker(object):
    def __init__(self, path_aria, product_dict, lc=None):
        self.path_aria = path_aria
        self.product_dict = product_dict
        self.path_bbox = os.path.join(
            self.path_aria, 'productBoundingBox', 'productBoundingBox.json')
        if lc is None:
            self.lc = [11, 12, 90, 95]
        else:
            self.lc = lc
        osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')

    def __call__(self, proj, bounds, arrres, outputFormat='ENVI', test=False):
        """
        view=True to plot the mask; test=True to apply mask to dummy data
        """
        # File must be physically extracted, cannot proceed with VRT format.
        # Defaulting to ENVI format.
        if outputFormat == 'VRT':
            outputFormat = 'ENVI'

        # crop the raw nlcd in memory
        PATH_NLCD = (
            '/vsizip/vsicurl/https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019_'
            'land_cover_l48_20210604.zip/nlcd_2019_land_cover_l48_20210604.img')

        LOGGER.info('Cropping raw NLCD...')
        ds_crop = crop_ds(PATH_NLCD, self.path_bbox)

        # resample the mask to the size of the products (mimic ariaExtract)
        LOGGER.info('Resampling NLCD to selected bbox...')
        ds_resamp = resamp(ds_crop, proj, bounds, arrres)

        # mask the landcover classes in lc
        LOGGER.info('Masking water and wetlands...')
        ds_mask = make_mask(ds_resamp, self.lc)

        # write mask to disk
        path_mask = os.path.join(self.path_aria, 'mask')
        os.mkdirs(path_mask) if not os.path.exists(path_mask) else ''
        dst = os.path.join(path_mask, 'NLCD_cros.path.msk')
        ds = osgeo.gdal.Translate(
            dst, ds_mask, format=outputFormat, outputType=osgeo.gdal.GDT_Byte)
        ds1 = osgeo.gdal.BuildVRT(dst + '.vrt', ds)

        # save a view of the mask
        arr = ds.ReadAsArray()
        plt.imshow(arr, cmap=plt.cm.Greys_r, interpolation='nearest')
        plt.colorbar()
        plt.title(f'Resampled mask\nDims: {arr.shape}')
        plt.savefig(f'{os.path.splitext(dst)[0]}.png')

        if test:
            # create temporary reference file to access geolocation info
            warp_options = osgeo.gdal.WarpOptions(
                format="MEM", outputBounds=bounds)
            self.path_temp = osgeo.gdal.Warp(
                '', self.product_dict[0][1], xRes=arrres[0], yRes=arrres[1],
                targetAlignedPixels=True, options=warp_options)
            self.__test__(ds_mask)
            self.path_temp = None

        ds.FlushCache()
        return dst

    def __test__(self, ds_maskre):
        # apply mask to dummy data FOR TESTING
        ds_dummy = self._dummy_data(self.path_temp)
        ds_masked = self._apply_mask(ds_maskre, ds_dummy)

        arr = ds_masked.ReadAsArray()
        arr = np.where(arr > 9990, np.nan, arr)
        plt.imshow(arr, cmap='jet_r', interpolation='nearest')
        plt.colorbar()
        plt.title('Mask applied to dummy data')
        plt.show()

    def _dummy_data(self, ds2match):
        """ Create raster of dummy data using the dem (for sizing); For test """
        if isinstance(ds2match, str) and os.path.exists(ds2match):
            arr = osgeo.gdal.Open(ds2match).ReadAsArray()
        else:
            arr = ds2match.ReadAsArray()

        # arr = np.random.rand(*arr.shape)
        arr = np.ones(arr.shape) * 10
        ds = arr2ds(ds2match, arr)
        arr1 = ds.ReadAsArray()
        return ds

    def _apply_mask(self, ds_mask, ds_2mask):
        """ Apply mask to test viewing """
        arr = ds_mask.ReadAsArray()
        arr = np.where(arr == 0, np.nan, 1)
        masked = ds_2mask.ReadAsArray() * arr
        ds_masked = arr2ds(ds_mask, masked)
        return ds_masked
