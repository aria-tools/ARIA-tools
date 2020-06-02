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
from osgeo import gdal, ogr, gdalconst
from ARIAtools.shapefile_util import open_shapefile

def prep_mask(product_dict, maskfilename, bbox_file, prods_TOTbbox, proj, amp_thresh=None, arrshape=None, workdir='./', outputFormat='ENVI', num_threads='2'):
    '''
        Function to load and export mask file.
        If "Download" flag is specified, GSHHS water mask will be donwloaded on the fly.
        If the full resolution NLCD landcover data is given (NLCD...img) it will be cropped to match product
    '''

    # Import functions
    from ARIAtools.vrtmanager import renderOGRVRT

    _world_watermask = [' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L1.shp',' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L2.shp',' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L3.shp', ' /vsizip/vsicurl/http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip/GSHHS_shp/f/GSHHS_f_L4.shp',' /vsizip/vsicurl/https://osmdata.openstreetmap.de/download/land-polygons-complete-4326.zip/land-polygons-complete-4326/land_polygons.shp']

    # If specified DEM subdirectory exists, delete contents
    workdir=os.path.join(workdir,'mask')
    if os.path.exists(workdir) and os.path.abspath(maskfilename)!=os.path.abspath(os.path.join(workdir,os.path.basename(maskfilename).split('.')[0].split('uncropped')[0]+'.msk')) and os.path.abspath(maskfilename)!=os.path.abspath(os.path.join(workdir,os.path.basename(maskfilename).split('.')[0]+'.msk')) or maskfilename.lower()=='download':
        for i in glob.glob(os.path.join(workdir,'*.*')): os.remove(i)
    if not os.path.exists(workdir):
        os.mkdir(workdir)


    # Get bounds of user bbox_file
    bounds=open_shapefile(bbox_file, 0, 0).bounds

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
        outputFormat='ENVI'

    # Download mask
    if maskfilename.lower()=='download':
        print("***Downloading water mask... ***")
        maskfilename=os.path.join(workdir,'watermask'+'.msk')
        os.environ['CPL_ZIP_ENCODING'] = 'UTF-8'
        ###Make coastlines/islands union VRT
        renderOGRVRT(os.path.join(workdir,'watermsk_shorelines.vrt'), _world_watermask[::2])

        ###Make lakes/ponds union VRT
        renderOGRVRT(os.path.join(workdir,'watermsk_lakes.vrt'), _world_watermask[1::2])

        ###Initiate water-mask with coastlines/islands union VRT
        # save uncropped mask
        gdal.Rasterize(os.path.join(workdir,'watermask_uncropped.msk'), os.path.join(workdir,'watermsk_shorelines.vrt'), options=gdal.RasterizeOptions(format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], burnValues=[1], layers='merged'))
        gdal.Translate(os.path.join(workdir,'watermask_uncropped.msk.vrt'), os.path.join(workdir,'watermask_uncropped.msk'), options=gdal.TranslateOptions(format="VRT"))
        # save cropped mask
        gdal.Warp(maskfilename, os.path.join(workdir,'watermask_uncropped.msk.vrt'), options=gdal.WarpOptions(format=outputFormat, outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
        update_file=gdal.Open(maskfilename,gdal.GA_Update)
        update_file.SetProjection(proj)
        update_file.GetRasterBand(1).SetNoDataValue(0.) ; del update_file
        gdal.Translate(maskfilename+'.vrt', maskfilename, options=gdal.TranslateOptions(format="VRT"))

        ###Must take inverse of lakes/ponds union because of opposite designation (1 for water, 0 for land) as desired (0 for water, 1 for land)
        lake_masks=gdal.Rasterize('', os.path.join(workdir,'watermsk_lakes.vrt'), options=gdal.RasterizeOptions(format='MEM', outputBounds=bounds, outputType=gdal.GDT_Byte, width=arrshape[1], height=arrshape[0], burnValues=[1], layers='merged', inverse=True))
        lake_masks.SetProjection(proj)
        lake_masks=lake_masks.ReadAsArray()

        if amp_thresh:
            ###Make average amplitude mask
            # Delete existing average amplitude file
            for i in glob.glob(os.path.join(workdir,'avgamplitude*')): os.remove(i)

            # building the VRT
            gdal.BuildVRT(os.path.join(workdir,'avgamplitude') +'.vrt', product_dict[0], options=gdal.BuildVRTOptions(resolution='highest', resampleAlg='average'))
            # taking average of all scenes and generating raster
            gdal.Warp(os.path.join(workdir,'avgamplitude'), os.path.join(workdir,'avgamplitude')+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, resampleAlg='average', width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
            # Update VRT
            gdal.Translate(os.path.join(workdir,'avgamplitude')+'.vrt', os.path.join(workdir,'avgamplitude'), options=gdal.TranslateOptions(format="VRT"))
            print('Saved average amplitude here: '+ os.path.join(workdir,'avgamplitude'))

            # Using specified amplitude threshold, pass amplitude mask
            amp_file=gdal.Open(os.path.join(workdir,'avgamplitude')).ReadAsArray()
            amp_file = np.where(amp_file < float(amp_thresh), 0, 1)
        else:
            amp_file=np.ones((lake_masks.shape[0],lake_masks.shape[1]))

        ###Update water-mask with lakes/ponds union and average amplitude
        mask_file = gdal.Open(maskfilename,gdal.GA_Update)
        mask_file.GetRasterBand(1).WriteArray(lake_masks*gdal.Open(maskfilename).ReadAsArray()*amp_file)
        #Delete temp files
        del lake_masks, amp_file, mask_file
        os.remove(os.path.join(workdir,'watermsk_shorelines.vrt')); os.remove(os.path.join(workdir,'watermsk_lakes.vrt'))

    if os.path.basename(maskfilename).lower().startswith('nlcd'):
        print("***Downloading and cropping the NLCD mask...***")
        maskfilename = NLCDMasker(os.path.dirname(workdir))(proj, bounds, arrshape) ## write mask to disk

    # Load mask
    try:
        # Check if uncropped/cropped maskfiles exist in 'mask' subdirectory
        if not os.path.exists(os.path.join(workdir,os.path.basename(maskfilename).split('.')[0]+'.msk')):
            # save uncropped masfile
            gdal.BuildVRT(os.path.join(workdir,os.path.basename(maskfilename).split('.')[0]+'_uncropped.msk.vrt'), maskfilename, options=gdal.BuildVRTOptions(outputBounds=bounds))
            # update maskfilename
            maskfilename=os.path.join(workdir,os.path.basename(maskfilename).split('.')[0].split('uncropped')[0]+'.msk')
            # save cropped maskfile
            gdal.Warp(maskfilename, os.path.join(workdir,os.path.basename(maskfilename).split('.')[0]+'_uncropped.msk.vrt'), options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
            mask_file = gdal.Open(maskfilename,gdal.GA_Update)
            mask_file.SetProjection(proj) ; del mask_file
            gdal.Translate(maskfilename+'.vrt', maskfilename, options=gdal.TranslateOptions(format="VRT"))
        #pass cropped DEM
        mask=gdal.Warp('', maskfilename, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
        mask.SetProjection(proj)
        mask.SetDescription(maskfilename)
    except:
        raise Exception('Failed to open user mask')

    return mask

class NLCDMasker(object):
    def __init__(self, path_aria, lc=[11, 12, 90, 95]):
        self.path_aria = path_aria
        self.path_bbox = op.join(self.path_aria, 'productBoundingBox', 'productBoundingBox.shp')
        self.path_dem  = op.join(self.path_aria, 'DEM', 'SRTM_3arcsec.dem')
        self.lc        = lc # landcover classes to mask
        gdal.PushErrorHandler('CPLQuietErrorHandler')

    def __call__(self, proj, bounds, arrshape, test=False):
        """ view=True to plot the mask; test=True to apply mask to dummy data """
        import matplotlib.pyplot as plt
        ## crop the raw nlcd in memory
        path_nlcd = '/vsizip/vsicurl/https://s3-us-west-2.amazonaws.com/mrlc/NLCD_2016'\
                    '_Land_Cover_L48_20190424.zip/NLCD_2016_Land_Cover_L48_20190424.img'

        ds_crop   = crop_ds(path_nlcd, self.path_bbox)

        ## mask the landcover classes in lc
        ds_mask = make_mask(ds_crop, self.lc)

        ## resample the mask to the size of the products (mimic ariaExtract)
        ds_maskre = resamp(ds_mask, proj, bounds, arrshape)

        ## write mask to disk
        path_mask = op.join(self.path_aria, 'mask')
        os.mkdirs(path_mask) if not op.exists(path_mask) else ''
        dst = op.join(path_mask, 'NLCD_crop.msk')
        ds  = gdal.Translate(dst, ds_maskre, format='ENVI', outputType=gdal.GDT_UInt16, noData=0)
        gdal.BuildVRT(dst + '.vrt' ,ds)

        ## save a view of the mask
        arr = ds.ReadAsArray()
        plt.imshow(arr, cmap=plt.cm.Greys_r, interpolation='nearest')
        plt.colorbar(); plt.title('Resampled mask\nDims: {}'.format(arr.shape))
        plt.savefig(op.join(path_mask, 'NLCD_crop.msk.png'))

        if test: self.__test__(ds_maskre)

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
        masked = ds_2mask.ReadAsArray() * arr
        ds_masked = arr2ds(ds_mask, masked)
        return ds_masked

def crop_ds(path_raster, path_poly, path_write=''):
    """ Crop to a raster (or path to one) using a polygon stored on disk """
    if isinstance(path_raster, str) and (op.exists(path_raster) or path_raster.startswith('/vsi')):
        ds  = gdal.Open(path_raster)
    else:
        raise FileNotFoundError(path_raster)
    ds  = gdal.Warp(path_write, ds, format='VRT', cutlineDSName=path_poly,
                                             cropToCutline=True, dstNodata=9999)
    if path_write: ds_vrt = gdal.BuildVRT('{}.vrt'.format(path_write), ds)
    return ds

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
    """ Resample a dataset from src dimensions using outputs from merged_productbbox """
    path = src.GetDescription()
    path = path if op.exists(path) else ''
    if isinstance(src, str) and op.exists(src):
        src = gdal.Open(src, gdal.GA_ReadOnly)

    ## compute geotranform
    height, width = arrshape
    pix_width = (bounds[2]-bounds[0]) / width
    pix_height = (bounds[1] - bounds[3]) / height
    gt = [bounds[0], pix_width, 0, bounds[3], 0, pix_height]

    if path:
        dst = gdal.GetDriverByName('ENVI').Create(path, width, height, 1, gdalconst.GDT_Float32)
    else:
        dst = gdal.GetDriverByName('MEM').Create('', width, height, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform(gt)
    dst.SetProjection(proj)
    gdal.ReprojectImage(src, dst, src.GetProjection(), proj, gdalconst.GRA_NearestNeighbour)
    del src
    return dst

def arr2ds(ds_orig, arr, noData=9999):
    """
    Create a new dataset identical to ds_orig except with values of 'arr'
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

## deprecated as it depends on DEM
def resamp_DEP(src, match, view=False):
    """ Resample a dataset from src dimensions to match """
    path = src.GetDescription()
    path = path if op.exists(path) else ''
    if isinstance(src, str) and op.exists(src):
        src = gdal.Open(src, gdal.GA_ReadOnly)
    if isinstance(match, str) and op.exists(match):
        match = gdal.Open(match, gdal.GA_ReadOnly)

    src_proj     = src.GetProjection()
    src_geotrans = src.GetGeoTransform()
    # match_nodata   = match.GetRasterBand(1).GetNoDataValue()

    match_proj     = match.GetProjection()
    match_geotrans = match.GetGeoTransform()
    width  = match.RasterXSize
    height = match.RasterYSize

    if path:
        dst = gdal.GetDriverByName('ENVI').Create(path, width, height, 1, gdalconst.GDT_Float32)
    else:
        dst = gdal.GetDriverByName('MEM').Create('', width, height, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform(match_geotrans)
    dst.SetProjection(match_proj)
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_NearestNeighbour)
    del src, match
    return dst
