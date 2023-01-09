#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, Brett Buzzanga, David Bekaert,
#         Marin Govorcin, Zhang Yunjun
# Copyright 2022, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import datetime as dt
import logging
import os
import shutil

import h5py

from netCDF4 import Dataset

import numpy as np

from osgeo import gdal

import pysolid

import xarray as xr

log = logging.getLogger(__name__)

# hardcode correction layer variables for GUNW
LYR_GROUP = 'science/grids/corrections/external/tides'
LYR_NAMES = ['solidEarthTide']
DIM_NAMES = ['heightsMeta', 'latitudeMeta', 'longitudeMeta']


def create_parser(iargs=None):
    """ Parser to read command line arguments """
    p = argparse.ArgumentParser(description='Update ARIA products '
                                'with correction layer(s).')
    p.add_argument('-f', '--file', type=str, help='ARIA file')
    p.add_argument('-o', '--out_dir', default=os.getcwd(), type=str,
                   help='Directory to store results (Default=./).')
    p.add_argument('-of', '--output_format', default='ISCE', type=str,
                   help='Specify GDAL compatible format for '
                   'intermediate outputs (Default=ISCE).')
    p.add_argument('-v', '--gunw_version', default='1c', type=str,
                   help='Specify GUNW version for file update '
                   '(Default=1c).')
    p.add_argument('-ir', '--interp_result', action='store_true',
                   help='Intersect with DEM and interpolate layer '
                   'and deposit into separate subdirectory (Default=False).')
    p.add_argument('-kt', '--keep_tmp', action='store_false',
                   help='Do not remove tmp directory with '
                   'intermediate outputs (Default=False).')

    args = p.parse_args(args=iargs)
    args.argv = iargs if iargs else os.sys.argv[1:]

    return args


def create_gdalrast(inname, outname, output_format):
    """ Initiate gdal raster + VRT """
    # Initiate raster
    gdal.Warp(outname, inname + '.vrt',
              options=gdal.WarpOptions(format=output_format))
    # Update VRT
    gdal.Translate(outname + '.vrt', outname,
                   options=gdal.TranslateOptions(format="VRT"))

    return


def update_gunw(path_gunw: str, ds_ifg):
    """ Update the path_gunw file using input correction layers """
    # first need to delete the variable; only can seem to with h5
    with h5py.File(path_gunw, 'a') as h5:
        for k in LYR_GROUP.split():
            h5 = h5[k]
        for i in LYR_NAMES:
            del h5[i]
        for k in 'crs'.split():
            if k in h5.keys():
                del h5[k]

    with Dataset(path_gunw, mode='a') as ds:
        ds_grp = ds[LYR_GROUP]
        for dim in DIM_NAMES:
            # dimension may already exist if updating
            try:
                ds_grp.createDimension(dim, len(ds_ifg.coords[dim]))
                # necessary for transform
                v = ds_grp.createVariable(dim, np.float32, dim)
                v[:] = ds_ifg[dim]
                v.setncatts(ds_ifg[dim].attrs)
            except:
                pass

        for name in LYR_NAMES:
            da = ds_ifg[name]
            nodata = da.encoding['_FillValue']
            chunksize = da.encoding['chunksizes']
            v = ds_grp.createVariable(name, np.float32, DIM_NAMES,
                                      chunksizes=chunksize,
                                      fill_value=nodata)
            v[:] = da.data
            v.setncatts(da.attrs)

        # add the projection
        v_proj = ds_grp.createVariable('crs', 'i')
        v_proj.setncatts(ds_ifg['crs'].attrs)

    log.info('Updated %s group in: %s',
             os.path.basename(LYR_GROUP), path_gunw)

    return


class setGUNW(object):
    def __init__(self, f: str, out_name: str, out_dir: str,
                 output_format: str, gunw_version: str, interp_result: bool,
                 keep_tmp: bool):
        self.path_gunw = f
        self.gunw_prefix = 'NETCDF:"' + self.path_gunw + '"'
        self.out_name = out_name
        self.out_dir = os.path.join(out_dir, self.out_name)
        self.tmp_dir = os.path.join(self.out_dir, 'tmp')
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
            os.makedirs(self.tmp_dir)
        self.output_format = output_format
        self.gunw_version = gunw_version
        self.interp_result = interp_result
        self.keep_tmp = keep_tmp
        # get geom filenames
        self.inc_angle_fn = self.gunw_prefix + \
            ':/science/grids/imagingGeometry/incidenceAngle'
        self.az_angle_fn = self.gunw_prefix + \
            ':/science/grids/imagingGeometry/azimuthAngle'
        # hardcode
        self.gdal_fmt = 'float32'
        self.nodata = 0.

    def __call__(self):
        self.dates = self.get_dates()
        self.ref_time = self.get_time()
        self.wavelength = self.get_wavelength()
        self.name = f'{self.dates[0]}_{self.dates[1]}'

        # access metadata
        self.atr, self.geotrans, self.proj, \
            self.hgt_field, self.zMeta, self.latMeta, \
            self.lonMeta = self.get_spatialarrs()

        # screen geom files
        self.inc_angle = self.screen_geomfiles(self.inc_angle_fn)
        self.az_angle = self.screen_geomfiles(self.az_angle_fn)

        # estimate differential SET ENU
        self.tide_e_ref, self.tide_n_ref, self.tide_u_ref = \
            self.get_setENU(self.dates[0])
        self.tide_e_sec, self.tide_n_sec, self.tide_u_sec = \
            self.get_setENU(self.dates[1])
        self.tide_e = np.subtract(self.tide_e_sec, self.tide_e_ref)
        self.tide_n = np.subtract(self.tide_n_sec, self.tide_n_ref)
        self.tide_u = np.subtract(self.tide_u_sec, self.tide_u_ref)

        # get SET LOS estimate for each height level
        self.tide_los = np.zeros(self.inc_angle.shape)
        for i in range(len(self.zMeta)):
            self.tide_los[i] = self.get_setLOS(self.inc_angle[i],
                                               self.az_angle[i])

        # Save file
        self.outname = os.path.join(self.out_dir, self.name)
        # building the VRT
        gdal.BuildVRT(self.outname + '.vrt', self.inc_angle_fn)
        gdal.Open(self.outname + '.vrt').SetMetadataItem(
                  self.hgt_field, gdal.Open(self.inc_angle_fn
                                            ).GetMetadataItem(self.hgt_field))
        # initiate raster
        create_gdalrast(self.outname, self.outname, self.output_format)
        # Update VRT with new raster
        update_file = gdal.Open(self.outname, gdal.GA_Update)
        for i in range(len(self.zMeta)):
            update_file.GetRasterBand(i+1).WriteArray(self.tide_los[i])
        del update_file

        # interpolate SET (bypass if you wish to have the downsampled cube)
        if self.interp_result:
            self.interpolate_set()

        # remove tmp directory
        if self.keep_tmp:
            shutil.rmtree(self.tmp_dir)

    def get_time(self):
        """ Get the center time of the secondary date from the filename """
        tt = os.path.basename(self.path_gunw).split('-')[7]
        return f'{tt[:2]}:{tt[2:4]}:{tt[4:]}'

    def get_wavelength(self):
        group = 'science/radarMetaData'
        with xr.open_dataset(self.path_gunw, group=group) as ds:
            wavelength = ds['wavelength'].item()
        return wavelength

    def get_dates(self):
        """ Get the ref/sec date from the filename """
        ref, sec = os.path.basename(self.path_gunw).split('-')[6].split('_')
        return [ref, sec]

    def screen_geomfiles(self, fname):
        """ Metadata layer quality check, correction applied if necessary """
        from ARIAtools.extractProduct import metadata_qualitycheck

        # create tmp file
        fname_key = fname.split('/')[-1]
        fname_tmp = os.path.join(self.tmp_dir, fname_key)
        if not os.path.exists(fname_tmp):
            os.makedirs(fname_tmp)
        fname_tmp += '/' + self.name
        gdal.BuildVRT(fname_tmp + '.vrt', fname)
        gdal.Open(fname_tmp + '.vrt').SetMetadataItem(
                  self.hgt_field, gdal.Open(fname
                                            ).GetMetadataItem(self.hgt_field))

        # initialize data array
        data_arr = gdal.Open(fname_tmp + '.vrt')
        # check and screen layer, if necessary
        verbose = False
        data_arr = metadata_qualitycheck(data_arr,
                                         fname_key,
                                         fname_tmp,
                                         verbose).data_array

        return data_arr.ReadAsArray()

    def get_spatialarrs(self):
        """ Get lat/lon/height arrays + attributes for metadata layers """
        # use inc angle file
        data_arr = gdal.Open(self.inc_angle_fn)
        # get attributes
        geotrans = data_arr.GetGeoTransform()
        proj = data_arr.GetProjection()
        y_siz, x_siz = data_arr.RasterYSize, data_arr.RasterXSize
        atr = {
            'LENGTH': y_siz,
            'WIDTH': x_siz,
            'X_FIRST': geotrans[0],
            'Y_FIRST': geotrans[3],
            'X_STEP':  geotrans[1],
            'Y_STEP': geotrans[-1],
        }

        # get heights
        zdim = data_arr.GetMetadataItem('NETCDF_DIM_EXTRA')[1:-1]
        hgt_field = f'NETCDF_DIM_{zdim}_VALUES'
        zMeta = np.array(data_arr.GetMetadataItem(
                         hgt_field)[1:-1].split(','), dtype='float32')
        # get lats and lons
        latMeta = np.linspace(geotrans[3],
                              geotrans[3] + (geotrans[5] * y_siz), y_siz)
        lonMeta = np.linspace(geotrans[0],
                              geotrans[0] + (geotrans[1] * x_siz), x_siz)

        return atr, geotrans, proj, hgt_field, zMeta, latMeta, lonMeta

    def get_setENU(self, scene):
        """ Estimate SET in ENU """
        # source:
        # https://github.com/insarlab/
        # PySolid/blob/main/docs/plot_grid_SET.ipynb

        # call
        dt_obj = dt.datetime.strptime(
            scene + '-' + self.ref_time, "%Y%m%d-%H:%M:%S"
        )
        tide_e, tide_n, tide_u = pysolid.calc_solid_earth_tides_grid(
            dt_obj, self.atr,
            display=False,
            verbose=False,
        )

        return tide_e, tide_n, tide_u

    def get_setLOS(self, inc_angle, az_angle):
        """ Estimate SET in LOS """
        # project ENU to radar line-of-sight (LOS)
        # with positive for motion towards satellite
        inc_angle_rad = np.deg2rad(inc_angle)
        az_angle_rad = np.deg2rad(az_angle-90)
        tide_los = (self.tide_e * np.sin(inc_angle_rad) * np.sin(az_angle_rad)
                    * -1 + self.tide_n * np.sin(inc_angle_rad)
                    * np.cos(az_angle_rad) + self.tide_u
                    * np.cos(inc_angle_rad))

        # Convert m to rad
        tide_los = np.divide(tide_los, float(self.wavelength) / (4*np.pi))

        return tide_los

    def interpolate_set(self):
        """ Interpolate SET and intersect with DEM """
        from ARIAtools.ARIAProduct import ARIA_standardproduct
        from ARIAtools.extractProduct import merged_productbbox, prep_dem, \
            finalize_metadata
        from ARIAtools.shapefile_util import open_shapefile

        nc_version = self.gunw_version
        # hardcode
        num_threads = 10
        url_version = None
        verbose = True
        croptounion = False
        minimumOverlap = 0.0081
        hgt_field = 'NETCDF_DIM_heightsMeta_VALUES'
        mask = None
        # run ARIA class
        standardproduct_info = ARIA_standardproduct(self.path_gunw,
                                                    bbox=None,
                                                    workdir=self.tmp_dir,
                                                    num_threads=num_threads,
                                                    url_version=url_version,
                                                    nc_version=nc_version,
                                                    verbose=verbose)
        # get product bounding box
        standardproduct_info.products[0], standardproduct_info.products[1], \
            standardproduct_info.bbox_file, prods_TOTbbox, \
            prods_TOTbbox_metadatalyr, arrshape, proj = merged_productbbox(
                standardproduct_info.products[0],
                standardproduct_info.products[1],
                os.path.join(self.tmp_dir, 'productBoundingBox'),
                standardproduct_info.bbox_file, croptounion,
                num_threads=num_threads,
                minimumOverlap=minimumOverlap, verbose=verbose)

        # get DEM and lat/lon arrays
        _, dem, lat, lon = prep_dem('download',
                                    standardproduct_info.bbox_file,
                                    prods_TOTbbox, prods_TOTbbox_metadatalyr,
                                    proj, arrshape=arrshape,
                                    workdir=self.tmp_dir,
                                    outputFormat=self.output_format,
                                    num_threads=num_threads)

        # initiate raster
        outname_interp = self.outname + '_interp'
        create_gdalrast(self.outname, outname_interp, self.output_format)

        # Interpolate/intersect with DEM before cropping
        bounds = open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds
        dem_bounds = [dem.GetGeoTransform()[0], dem.GetGeoTransform()[3] +
                      (dem.GetGeoTransform()[-1]*dem.RasterYSize),
                      dem.GetGeoTransform()[0] +
                      (dem.GetGeoTransform()[1]*dem.RasterXSize),
                      dem.GetGeoTransform()[3]]
        finalize_metadata(outname_interp, bounds, dem_bounds,
                          prods_TOTbbox, dem, lat, lon, hgt_field,
                          mask, self.output_format, verbose=verbose)


if __name__ == '__main__':
    args = create_parser()
    for out_name in LYR_NAMES:
        set_gunwobj = setGUNW(args.file, out_name, args.out_dir,
                              args.output_format, args.gunw_version,
                              args.interp_result, args.keep_tmp)
    set_gunwobj()
