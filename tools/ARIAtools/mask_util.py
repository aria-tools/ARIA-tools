#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & Brett Buzzanga & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS RESERVED
# United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import os.path as op
import glob
import numpy as np
import logging
import shutil
from osgeo import gdal, ogr, gdalconst

from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.vrtmanager import rasterAverage
from ARIAtools.logger import logger

log = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def prep_mask(product_dict, maskfilename, bbox_file, prods_TOTbbox, proj,
                    amp_thresh=None, arrshape=None, workdir='./',
                    outputFormat='ENVI', num_threads='2'):
    """Function to load and export mask file.
    If "Download" flag, GSHHS water mask will be donwloaded on the fly.
    If full resolution NLCD landcover data is given (NLCD...img) it gets cropped
    """

    # Import functions
    from ARIAtools.vrtmanager import renderOGRVRT
    _world_watermask = [f' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/'\
                        f'gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L{i}.shp'\
                        for i in range(1, 5)]
    _world_watermask.append( ' /vsizip/vsicurl/https://osmdata.openstreetmap.de'\
                             '/download/land-polygons-complete-4326.zip/'\
                             'land-polygons-complete-4326/land_polygons.shp')
    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir,'mask')
    os.makedirs(workdir, exist_ok=True)

    # Get bounds of user bbox_file
    bounds = open_shapefile(bbox_file, 0, 0).bounds

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
        outputFormat='ENVI'

    # Download mask
    if maskfilename.lower()=='download':
        log.info("***Downloading water mask... ***")
        maskfilename=os.path.join(workdir,'watermask'+'.msk')
        os.environ['CPL_ZIP_ENCODING'] = 'UTF-8'
        ###Make coastlines/islands union VRT
        renderOGRVRT(os.path.join(workdir,'watermsk_shorelines.vrt'), _world_watermask[::2])
        ###Make lakes/ponds union VRT
        renderOGRVRT(os.path.join(workdir,'watermsk_lakes.vrt'), _world_watermask[1::2])

        ###Initiate water-mask with coastlines/islands union VRT
        # save uncropped mask
        gdal.Rasterize(os.path.join(workdir,'watermask_uncropped.msk'), os.path.join(workdir,'watermsk_shorelines.vrt'),
                                    format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1],
                                    height=arrshape[0], burnValues=[1], layers='merged')
        gdal.Translate(os.path.join(workdir,'watermask_uncropped.msk.vrt'), os.path.join(workdir,'watermask_uncropped.msk'), format="VRT")
        # save cropped mask
        gdal.Warp(maskfilename, os.path.join(workdir,'watermask_uncropped.msk.vrt'), format=outputFormat, outputBounds=bounds,
                         outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)])

        update_file=gdal.Open(maskfilename,gdal.GA_Update)
        update_file.SetProjection(proj)
        update_file.GetRasterBand(1).SetNoDataValue(0.); del update_file
        gdal.Translate(f'{maskfilename}.vrt', maskfilename, format='VRT')

        ###Must take inverse of lakes/ponds union because of opposite designation (1 for water, 0 for land) as desired (0 for water, 1 for land)
        lake_masks=gdal.Rasterize('', os.path.join(workdir,'watermsk_lakes.vrt'),
                        format='MEM', outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1],
                        height=arrshape[0], burnValues=[1], layers='merged', inverse=True)

        lake_masks.SetProjection(proj)
        lake_masks=lake_masks.ReadAsArray()

        ## Update water-mask with lakes/ponds union
        mask_file = gdal.Open(maskfilename, gdal.GA_Update)
        mask_file.GetRasterBand(1).WriteArray(lake_masks*gdal.Open(maskfilename).ReadAsArray())
        #Delete temp files
        del lake_masks, mask_file

    ## Use NLCD Mask
    elif os.path.basename(maskfilename).lower().startswith('nlcd'):
            log.info("***Accessing and cropping the NLCD mask...***")
            maskfilename = NLCDMasker(os.path.dirname(workdir))(
                                            proj, bounds, arrshape, outputFormat)

    ## User specified mask
    else:
        # Path to local version of user specified mask
        user_mask   = maskfilename # for clarity
        user_mask_n = os.path.basename(os.path.splitext(user_mask)[0])
        local_mask  = os.path.join(workdir, f'{user_mask_n}.msk')
        local_mask_unc = os.path.join(workdir, f'{user_mask_n}_uncropped.msk')
        if user_mask == local_mask:
            log.warning('The mask you specified already exists in %s, '\
                    'using the existing one...', os.path.dirname(local_mask))

        else:
            # move the mask to the local directory and built a VRT for it
            gdal.UseExceptions()
            # shutil.copy(user_mask, local_mask_unc)
            ds = gdal.BuildVRT(f'{local_mask_unc}.vrt', user_mask,
                                                    outputBounds=bounds)

            assert ds is not None, f'GDAL could not open user mask: {user_mask}'

            # crop the user mask and write
            gdal.Warp(f'{local_mask}', ds,
                            format=outputFormat, cutlineDSName=prods_TOTbbox,
                            outputBounds=bounds, width=arrshape[1],
                            height=arrshape[0], multithread=True,
                            options=[f'NUM_THREADS={num_threads}'])

            # set projection of the local mask
            mask_file = gdal.Open(local_mask, gdal.GA_Update)
            mask_file.SetProjection(proj); del mask_file

            # create vrt for local cropped mask
            gdal.Translate(f'{local_mask}.vrt', local_mask, format='VRT')
        # assign new local mask for amp thresh
        maskfilename = local_mask

    ## Make average amplitude mask
    if amp_thresh is not None:
        amp_file = rasterAverage(os.path.join(workdir, 'avgamplitude'),
                        product_dict, bounds, prods_TOTbbox,
                        outputFormat=outputFormat, thresh=amp_thresh)

        # Update mask with average amplitude
        mask_file = gdal.Open(maskfilename, gdal.GA_Update)
        mask_arr  = mask_file.ReadAsArray()
        mask_file.GetRasterBand(1).WriteArray(mask_arr*amp_file)

        # Delete temp files
        del mask_file, amp_file

    # crop/expand mask to DEM size?
    print (maskfilename)
    breakpoint()
    mask = gdal.Warp('', maskfilename, format='MEM',
                    cutlineDSName=prods_TOTbbox, outputBounds=bounds,
                    width=arrshape[1], height=arrshape[0], multithread=True,
                    options=[f'NUM_THREADS={num_threads}'])

    mask.SetProjection(proj)
    mask.SetDescription(maskfilename)

    try:
        ## remove extra files
        os.remove(path_shorelines) if os.path.exists(path_shorelines) else ''
        os.remove(path_lakes) if os.path.exists(path_lakes) else ''
    except:
        pass

    return mask


def make_mask(ds_crop, lc):
    # values outside crop; 0 is really just extra check
    lc.extend([0,255])

    arr = ds_crop.ReadAsArray()
    for lclass in lc:
        arr = np.where(arr == lclass, np.nan, arr)

    arr = np.where(np.isnan(arr), np.nan, 1)
    ds  = arr2ds(ds_crop, arr)
    return ds


def resamp(src, proj, bounds, arrshape, view=False):
    """ Resample a dataset from src dims using result of merged_productbbox """
    if isinstance(src, str) and op.exists(src):
        src = gdal.Open(src, gdal.GA_ReadOnly)
    path = path if op.exists(src.GetDescription()) else ''

    ## compute geotranform
    height, width = arrshape
    pix_width = (bounds[2]-bounds[0]) / width
    pix_height = (bounds[1] - bounds[3]) / height
    gt = [bounds[0], pix_width, 0, bounds[3], 0, pix_height]

    if path:
        dst = gdal.GetDriverByName('ENVI').Create(path, width, height, 1,
                                                        gdalconst.GDT_Float32)

    else:
        dst = gdal.GetDriverByName('MEM').Create(path, width, height, 1,
                                                        gdalconst.GDT_Int16)

    dst.SetGeoTransform(gt)
    dst.SetProjection(proj)
    gdal.ReprojectImage(src, dst, src.GetProjection(), proj,
                                                gdalconst.GRA_NearestNeighbour)
    del src
    return dst


def crop_ds(path_raster, path_poly, path_write=''):
    """ Crop path_raster (or dataset opened with gdal) to a polygon on disk """
    assert op.exists(path_raster) or path_raster.startswith('/vsi'), 'Invalid path to NLCD mask'
    path_raster = gdal.BuildVRT('', path_raster)
    ds          = gdal.Warp(path_write, path_raster, format='VRT', cutlineDSName=path_poly,
                                             cropToCutline=True, dstNodata=0)

    if path_write: ds_vrt = gdal.BuildVRT('{}.vrt'.format(path_write), ds)
    return ds


def arr2ds(ds_orig, arr, noData=0):
    """Create a new dataset identical to ds_orig except with values of 'arr'

    Performed in memory
    """
    ds_orig = gdal.Open(ds_orig, gdal.GA_ReadOnly) if isinstance(ds_orig, str) else ds_orig
    arr     = np.where(np.isnan(arr), noData, arr)

    if arr.ndim == 3:
        ds_out = gdal.GetDriverByName('MEM').Create('', ds_orig.RasterXSize,
                        ds_orig.RasterYSize, arr.shape[0], gdal.GDT_Float32)

        for i in range(arr.shape[0]):
            ds_out.GetRasterBand(i+1).WriteArray(arr[i, :, :])
            ds_out.GetRasterBand(i+1).SetNoDataValue(noData)
    else:
        ds_out = gdal.GetDriverByName('MEM').Create('', ds_orig.RasterXSize,
                    ds_orig.RasterYSize, 1, gdal.GDT_Float32)
        ds_out.GetRasterBand(1).WriteArray(arr)
        ds_out.GetRasterBand(1).SetNoDataValue(noData)

    ds_out.SetGeoTransform(ds_orig.GetGeoTransform())
    ds_out.SetProjection(ds_orig.GetProjection())
    del ds_orig
    return ds_out


class NLCDMasker(object):
    def __init__(self, path_aria, lc=[11, 12, 90, 95]):
        self.path_aria = path_aria
        self.path_bbox = op.join(self.path_aria, 'productBoundingBox',
                                                    'productBoundingBox.json')
        self.path_dem  = op.join(self.path_aria, 'DEM', 'SRTM_3arcsec.dem')
        self.lc        = lc # landcover classes to mask
        gdal.PushErrorHandler('CPLQuietErrorHandler')

    def __call__(self, proj, bounds, arrshape, outputFormat='ENVI', test=False):
        """ view=True to plot the mask; test=True to apply mask to dummy data """
        import matplotlib.pyplot as plt

        # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
        if outputFormat=='VRT':
            outputFormat='ENVI'

        ## crop the raw nlcd in memory
        path_nlcd = '/vsizip/vsicurl/https://s3-us-west-2.amazonaws.com/mrlc/nlcd_2019'\
                    '_land_cover_l48_20210604.zip/nlcd_2019_land_cover_l48_20210604.img'

        log.info('Cropping raw NLCD...')
        ds_crop   = crop_ds(path_nlcd, self.path_bbox)

        ## resample the mask to the size of the products (mimic ariaExtract)
        log.info('Resampling NLCD to selected bbox...')
        ds_resamp = resamp(ds_crop, proj, bounds, arrshape)

        ## mask the landcover classes in lc
        log.info('Masking water and wetlands...')
        ds_mask   = make_mask(ds_resamp, self.lc)

        ## write mask to disk
        path_mask = op.join(self.path_aria, 'mask')
        os.mkdirs(path_mask) if not op.exists(path_mask) else ''
        dst = op.join(path_mask, 'NLCD_crop.msk')
        ds  = gdal.Translate(dst, ds_mask, options=gdal.TranslateOptions(format=outputFormat, outputType=gdal.GDT_Byte))
        gdal.BuildVRT(dst + '.vrt' ,ds)

        ## save a view of the mask
        arr = ds.ReadAsArray()
        plt.imshow(arr, cmap=plt.cm.Greys_r, interpolation='nearest')
        plt.colorbar();
        plt.title(f'Resampled mask\nDims: {arr.shape}')
        plt.savefig(op.join(path_mask, 'NLCD_crop.msk.png'))

        if test: self.__test__(ds_mask)

        return dst

    def __test__(self, ds_maskre):
        ## apply mask to dummy data FOR TESTING
        import matplotlib.pyplot as plt
        ds_dummy  = self._dummy_data(self.path_dem)
        ds_masked = self._apply_mask(ds_maskre, ds_dummy)

        arr = ds_masked.ReadAsArray()
        arr = np.where(arr>9990, np.nan, arr)
        plt.imshow(arr, cmap='jet_r', interpolation='nearest'); plt.colorbar()
        plt.title('Mask applied to dummy data'); plt.show()

    def _dummy_data(self, ds2match):
        """ Create raster of dummy data using the dem (for sizing); For test """
        if isinstance(ds2match, str) and op.exists(ds2match):
            arr = gdal.Open(ds2match, gdal.GA_ReadOnly).ReadAsArray()
        else:
            arr = ds2match.ReadAsArray()

        # arr = np.random.rand(*arr.shape)
        arr = np.ones(arr.shape) * 10
        ds  = arr2ds(ds2match, arr)
        arr1 = ds.ReadAsArray()
        return ds

    def _apply_mask(self, ds_mask, ds_2mask):
        """ Apply mask to test viewing """
        arr = ds_mask.ReadAsArray()
        arr = np.where(arr == 0, np.nan, 1)
        masked    = ds_2mask.ReadAsArray() * arr
        ds_masked = arr2ds(ds_mask, masked)
        return ds_masked
