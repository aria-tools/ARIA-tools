#!/usr/bin/env python3
import os
import os.path as op
import argparse
import numpy as np
import matplotlib.pyplot as plt

from osgeo import gdal, gdalconst

def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Command line interface to create a mask using the NLCD land cover dataset\nNLCD can be downloaded from: '\
                                                 'https://www.mrlc.gov/data?f%5B0%5D=category%3Aland%20cover&f%5B1%5D=year%3A2016'\
                                                 '\nMust specify the directory to ARIA results (e.g. from ariaTSsetup.py or ariaExtract.py) and NLCD product',
                                     epilog='Examples of use:\n\t ariaNLCDmask.py ARIAresults ./NLCD_2016_Land_Cover_L48_20190424.img --view',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('path_aria', type=str, help='Path to file to ARIA results (i.e., where DEM and productBoundingbox directories are located)')
    parser.add_argument('path_nlcd', type=str, help='Path to file to NLCD product')
    parser.add_argument('-t', '--test', dest='test', action='store_true', help='Run test on dummy data')
    parser.add_argument('-v', '--view', dest='view', action='store_true', help='Show plots of steps; for debugging')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    # iargs  = 'test'
    if len(os.sys.argv) < 3:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)
    # print(inps); quit()
    return inps

def crop_ds(path_raster, path_poly):
    """ Crop to a raster (or path to one) using the a polygon stored on disk """
    if isinstance(path_raster, str) and op.exists(path_raster):
        ds  = gdal.Open(path_raster)
    ds  = gdal.Warp('', ds, format='VRT', cutlineDSName=path_poly, cropToCutline=True, dstNodata=9999)
    # dst = op.join(self.path_data, 'cropped_NLCD')
    # ds  = gdal.Warp(dst, ds, format='ENVI', cutlineDSName=self.path_prod, cropToCutline=True, dstNodata=9999)
    # ds_vrt = gdal.BuildVRT('{}.vrt'.format(dst), ds)
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

def resamp(src, bounds, proj, arrshape, view=False):
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
    # dst.SetProjection(proj)
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
    print (match_geotrans); quit()
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

class NLCDMasker(object):
    def __init__(self, path_aria):
        self.path_aria = path_aria
        self.path_bbox = op.join(self.path_aria, 'productBoundingBox', 'productBoundingBox.shp')
        self.path_dem  = op.join(self.path_aria, 'DEM', 'SRTM_3arcsec.dem')

    def __call__(self, path_nlcd, lc=[11,12, 90, 95], view=False, test=False):
        ## crop the raw nlcd in memory
        ds_crop = crop_ds(path_nlcd, self.path_bbox)

        ## mask the landcover classes in lc
        ds_mask = make_mask(ds_crop, lc)

        ## resample the mask to the size of the products (mimic ariaExtract)
        proj     = 'GEOGCS["unknown",DATUM["unnamed",SPHEROID["Spheroid",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST]]'
        bounds   = (-76.655, 36.75, -75.928, 37.225)
        arrshape = [570, 872]
        ds_maskre = resamp(ds_mask, bounds, proj, arrshape)
        # ds_maskre = resamp_DEP(ds_mask, self.path_dem)

        ## write mask to disk
        path_mask = op.join(self.path_aria, 'mask')
        os.mkdirs(path_mask) if not op.exists(path_mask) else ''
        dst = op.join(path_mask, 'NLCD_crop.msk')
        ds  = gdal.Translate(dst, ds_maskre, format='ENVI', outputType=gdal.GDT_UInt16, noData=0)
        gdal.BuildVRT(dst + '.vrt' ,ds)

        ## optionally view the mask
        if view:
            arr = ds.ReadAsArray()
            plt.imshow(arr, cmap=plt.cm.Greys_r); plt.colorbar()
            plt.title('Resampled mask\nDims: {}'.format(arr.shape)); plt.show()

        if test: self.__test__(ds_maskre)

        return dst

    def __test__(self, ds_maskre):
        ## apply mask to dummy data FOR TESTING
        ds_dummy  = self._dummy_data(self.path_dem)
        ds_masked = self._apply_mask(ds_maskre, ds_dummy)

        arr = ds_masked.ReadAsArray()
        arr = np.where(arr>9990, np.nan, arr)
        plt.imshow(arr, cmap='jet_r'); plt.colorbar()
        plt.title('Mask applied to dummy data'); plt.show()

    def _dummy_data(self, ds2match):
        """ Create a raster of dummy data using the dem (for sizing) """
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
        arr = ds_mask.ReadAsArray()
        arr = np.where(arr == 0, np.nan, 1)
        masked = ds_2mask.ReadAsArray() * arr
        ds_masked = arr2ds(ds_mask, masked)
        return ds_masked

if __name__ == '__main__':
    # inps = cmdLineParse()
    # path_aria = op.abspath(inps.path_aria)
    # path_nlcd = op.abspath(inps.path_nlcd)
    # gdal.SetConfigOption('AWS_NO_SIGN_REQUEST', 'YES')
    # gdal.SetConfigOption('AWS_SECRET_ACCESS_KEY', '')
    # ds = gdal.Info('/vsizip/vsicurl/https://s3-us-west-2.amazonaws.com/mrlc/NLCD_2016_Land_Cover_L48_20190424.zip/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img')
    # NLCDMasker(path_aria)(path_nlcd, view=inps.view, test=inps.test)
    path_aria = '/Volumes/BB_1TB/data/VLM/Sentinel1/track_004/Results/SP-HR/TEMP'
    path_nlcd = '/Volumes/BB_1TB/data/mask/NLCD_2016_Land_Cover_L48_20190424/NLCD_2016_Land_Cover_L48_20190424.img'

    NLCDMasker(path_aria)(path_nlcd, view=True, test=False)
