# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Extract and organize specified layer(s).
If no layer is specified, extract product bounding box shapefile(s)
"""
import os
import sys

# Force HDF5 and GDAL environment overrides before C-libraries initialize
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['GDAL_PAM_ENABLED'] = 'NO'

import glob
import time
import copy
import json
import shutil
import logging
import datetime
import tarfile
import subprocess
import threading

import dask
import rioxarray
import rasterio
import osgeo
import osgeo_utils.gdal_calc
import pyproj
import numpy as np
import scipy.interpolate
import shapely.geometry
import h5py

import ARIAtools.product
import ARIAtools.util.ionosphere
import ARIAtools.util.interp
import ARIAtools.util.vrt
import ARIAtools.util.shp
import ARIAtools.util.misc
import ARIAtools.util.seq_stitch

from ARIAtools.constants import ARIA_PX_SIZES

LOGGER = logging.getLogger(__name__)
METADATA_STAGING_DRIVER = 'GTiff'
# metadata layer quality check, correction applied if necessary
# only apply to geometry layers and prods derived from older ISCE versions
GEOM_LYRS = ['bPerpendicular', 'bParallel', 'incidenceAngle',
             'lookAngle', 'azimuthAngle']

class MetadataQualityCheck:
    """
    Metadata quality control function.
    Artifacts recognized based off of covariance of cross-profiles.
    Bug-fix varies based off of layer of interest.
    Verbose mode generates a series of quality control plots with
    these profiles.
    """

    def __init__(self, data_array, prod_key, outname, verbose=None):
        # Pass inputs
        self.data_array = data_array
        self.prod_key = prod_key
        self.outname = outname
        self.verbose = verbose
        self.data_array_band = data_array.GetRasterBand(1).ReadAsArray()

        # mask by nodata value
        no_data_value = self.data_array.GetRasterBand(1).GetNoDataValue()
        self.data_array_band = np.ma.masked_where(
            self.data_array_band == no_data_value, self.data_array_band)

        # Run class
        self.__run__()

    def __truncateArray__(self, data_array_band, Xmask, Ymask):
        # Mask columns/rows which are entirely made up of 0s
        # first must crop all columns with no valid values
        nancols = np.all(data_array_band.mask == True, axis=0)
        data_array_band = data_array_band[:, ~nancols]
        Xmask = Xmask[:, ~nancols]
        Ymask = Ymask[:, ~nancols]
        # first must crop all rows with no valid values
        nanrows = np.all(data_array_band.mask == True, axis=1)
        data_array_band = data_array_band[~nanrows]
        Xmask = Xmask[~nanrows]
        Ymask = Ymask[~nanrows]

        return data_array_band, Xmask, Ymask

    def __getCovar__(self, prof_direc, profprefix=''):
        from scipy.stats import linregress
        # Mask columns/rows which are entirely made up of 0s
        if (self.data_array_band.mask.size != 1 and
                True in self.data_array_band.mask):
            Xmask, Ymask = np.meshgrid(
                np.arange(0, self.data_array_band.shape[1], 1),
                np.arange(0, self.data_array_band.shape[0], 1))
            self.data_array_band, Xmask, Ymask = self.__truncateArray__(
                self.data_array_band, Xmask, Ymask)

        # append prefix for plot names
        prof_direc = profprefix + prof_direc

        # Cycle between range and azimuth profiles
        rsquaredarr = []
        std_errarr = []

        # iterate through transpose of matrix if looking in azimuth
        data_array_band = (
            self.data_array_band.T if 'azimuth' in prof_direc else
            self.data_array_band)
        for i in enumerate(data_array_band):
            mid_line = i[1]
            xarr = np.array(range(len(mid_line)))
            # remove masked values from slice
            if mid_line.mask.size != 1:
                if True in mid_line.mask:
                    xarr = xarr[~mid_line.mask]
                    mid_line = mid_line[~mid_line.mask]

            # chunk array to better isolate artifacts
            chunk_size = 4

            for j in range(0, len(mid_line.tolist()), chunk_size):
                chunk = mid_line.tolist()[j:j + chunk_size]
                xarr_chunk = xarr[j:j + chunk_size]
                # make sure each iteration contains at least minimum number of
                # elements
                if j == range(0, len(mid_line.tolist()), chunk_size)[-2] and \
                        len(mid_line.tolist()) % chunk_size != 0:
                    chunk = mid_line.tolist()[j:]
                    xarr_chunk = xarr[j:]
                # linear regression and get covariance
                slope, bias, rsquared, p_value, std_err = linregress(
                    xarr_chunk, chunk)
                rsquaredarr.append(abs(rsquared)**2)
                std_errarr.append(std_err)
                # terminate early if last iteration would have small chunk size
                if len(chunk) > chunk_size:
                    break

            # exit loop/make plots in verbose mode if R^2 and standard error
            # anomalous, or if on last iteration
            if (min(rsquaredarr) < 0.9 and max(std_errarr) > 0.01) or \
                    (i[0] == (len(data_array_band) - 1)):
                if self.verbose:
                    # Make quality-control plots
                    import matplotlib.pyplot as plt
                    ax0 = plt.figure().add_subplot(111)
                    ax0.scatter(xarr, mid_line, c='k', s=7)
                    refline = np.linspace(min(xarr), max(xarr), 100)
                    ax0.plot(refline, (refline * slope) + bias,
                             linestyle='solid', color='red')
                    ax0.set_ylabel('%s array' % (self.prod_key))
                    ax0.set_xlabel('distance')
                    ax0.set_title('Profile along %s' % (prof_direc))
                    ax0.annotate(
                        'R\u00b2 = %f\nStd error= %f' % (
                            min(rsquaredarr), max(std_errarr)),
                        (0, 1), xytext=(4, -4), xycoords='axes fraction',
                        textcoords='offset points', fontweight='bold',
                        ha='left', va='top')

                    if min(rsquaredarr) < 0.9 and max(std_errarr) > 0.01:
                        ax0.annotate(
                            'WARNING: R\u00b2 and standard error\nsuggest '
                            'artifact exists', (1, 1), xytext=(4, -4),
                            xycoords='axes fraction',
                            textcoords='offset points', fontweight='bold',
                            ha='right', va='top')

                    plt.margins(0)
                    plt.tight_layout()
                    plt.savefig(os.path.join(
                        os.path.dirname(os.path.dirname(self.outname)),
                        'metadatalyr_plots', self.prod_key,
                        os.path.basename(self.outname) + '_%s.png' % (
                            prof_direc)))
                    plt.close()
                break

        return rsquaredarr, std_errarr

    def __run__(self):
        from scipy.linalg import lstsq

        # Get R^2/standard error across range
        rsquaredarr_rng, std_errarr_rng = self.__getCovar__('range')
        # Get R^2/standard error across azimuth
        rsquaredarr_az, std_errarr_az = self.__getCovar__('azimuth')

        # filter out normal values from arrays
        rsquaredarr = [0.97]
        std_errarr = [0.0015]

        if min(rsquaredarr_rng) < 0.97 and max(std_errarr_rng) > 0.0015:
            rsquaredarr.append(min(rsquaredarr_rng))
            std_errarr.append(max(std_errarr_rng))

        if min(rsquaredarr_az) < 0.97 and max(std_errarr_az) > 0.0015:
            rsquaredarr.append(min(rsquaredarr_az))
            std_errarr.append(max(std_errarr_az))

        # if R^2 and standard error anomalous, fix array
        if min(rsquaredarr) < 0.97 and max(std_errarr) > 0.0015:
            # Cycle through each band
            for i in range(1, 5):

                self.data_array_band = self.data_array.GetRasterBand(
                    i).ReadAsArray()

                # mask by nodata value
                no_data_value = self.data_array.GetRasterBand(
                    i).GetNoDataValue()
                self.data_array_band = np.ma.masked_where(
                    self.data_array_band == no_data_value,
                    self.data_array_band)
                negs_percent = ((self.data_array_band < 0).sum()
                                / self.data_array_band.size) * 100

                # Unique bug-fix for bPerp layers with sign-flips
                if ((self.prod_key == 'bPerpendicular' and
                    min(rsquaredarr) < 0.8 and max(std_errarr) > 0.1) and
                        (negs_percent != 100 or negs_percent != 0)):

                    # Circumvent Bperp sign-flip bug by comparing percentage of
                    # positive and negative values
                    self.data_array_band = abs(self.data_array_band)
                    if negs_percent > 50:
                        self.data_array_band *= -1
                else:

                    # regular grid covering the domain of the data
                    X, Y = np.meshgrid(
                        np.arange(0, self.data_array_band.shape[1], 1),
                        np.arange(0, self.data_array_band.shape[0], 1))

                    Xmask, Ymask = np.meshgrid(
                        np.arange(0, self.data_array_band.shape[1], 1),
                        np.arange(0, self.data_array_band.shape[0], 1))

                    # best-fit linear plane: for very large artifacts, must
                    # mask array for outliers to get best fit
                    if min(rsquaredarr) < 0.85 and max(std_errarr) > 0.0015:
                        maj_percent = ((self.data_array_band <
                                        self.data_array_band.mean()).sum()
                                       / self.data_array_band.size) * 100

                        # mask all values above mean
                        if maj_percent > 50:
                            self.data_array_band = np.ma.masked_where(
                                self.data_array_band > self.data_array_band.mean(),
                                self.data_array_band)

                        # mask all values below mean
                        else:
                            self.data_array_band = np.ma.masked_where(
                                self.data_array_band < self.data_array_band.mean(),
                                self.data_array_band)

                    # Mask columns/rows which are entirely made up of 0s
                    if (self.data_array_band.mask.size != 1 and
                            True in self.data_array_band.mask):
                        self.data_array_band, Xmask, Ymask = \
                            self.__truncateArray__(
                                self.data_array_band, Xmask, Ymask)

                    # truncated grid covering the domain of the data
                    Xmask = Xmask[~self.data_array_band.mask]
                    Ymask = Ymask[~self.data_array_band.mask]

                    self.data_array_band = self.data_array_band[
                        ~self.data_array_band.mask]

                    XX = Xmask.flatten()
                    YY = Ymask.flatten()
                    A = np.c_[XX, YY, np.ones(len(XX))]
                    C, _, _, _ = lstsq(A, self.data_array_band.data.flatten())

                    # evaluate it on grid
                    self.data_array_band = C[0] * X + C[1] * Y + C[2]

                    # mask by nodata value
                    no_data_value = self.data_array.GetRasterBand(
                        i).GetNoDataValue()
                    self.data_array_band = np.ma.masked_where(
                        self.data_array_band == no_data_value,
                        self.data_array_band)
                    np.ma.set_fill_value(self.data_array_band, no_data_value)

                # update band
                self.data_array.GetRasterBand(i).WriteArray(
                    self.data_array_band.filled())

                # Pass warning and get R^2/standard error across range/azimuth
                # (only do for first band)
                if i == 1:
                    # make sure appropriate unit is passed to print statement
                    lyrunit = "\N{DEGREE SIGN}"
                    if (self.prod_key == 'bPerpendicular' or
                            self.prod_key == 'bParallel'):
                        lyrunit = 'm'

                    LOGGER.warning((
                        "%s layer for IFG %s has R\u00b2 of %.4f and standard "
                        "error of %.4f%s, automated correction applied") % (
                        self.prod_key, os.path.basename(self.outname),
                        min(rsquaredarr), max(std_errarr), lyrunit))

                    rsquaredarr_rng, std_errarr_rng = self.__getCovar__(
                        'range', profprefix='corrected')

                    rsquaredarr_az, std_errarr_az = self.__getCovar__(
                        'azimuth', profprefix='corrected')

        self.data_array_band = None
        
        safe_return_array = self.data_array
        self.data_array = None
        
        return safe_return_array


def crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs):
    """
    Manage cropping of existing, extracted layers
    """
    LOGGER.debug('Cropping %s - %s', ifg_tag, lyrname)

    tmp_crop = outname + '_crop.tmp.tif'
    if os.path.exists(tmp_crop):
        os.remove(tmp_crop)

    crop_kwargs = gdal_warp_kwargs.copy()
    crop_kwargs['format'] = 'GTiff'
    crop_kwargs['multithread'] = False
    warp_options = osgeo.gdal.WarpOptions(**crop_kwargs)

    ds = osgeo.gdal.Warp(
        tmp_crop, outname + '.vrt', options=warp_options
    )
    ds = None

    out_crop = outname + '_crop'
    for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        f_to_rm = f"{out_crop}{ext}"
        if os.path.exists(f_to_rm):
            os.remove(f_to_rm)

    ds_env = osgeo.gdal.Translate(out_crop, tmp_crop, format='ENVI')
    if ds_env is not None:
        ds_env.FlushCache()
    ds_env = None

    if os.path.exists(tmp_crop):
        os.remove(tmp_crop)

    for crop_name in glob.glob(outname + '_crop*'):
        fname = os.path.basename(crop_name).replace('_crop', '')
        fname = os.path.join(os.path.dirname(crop_name), fname)
        os.rename(crop_name, fname)

    # Update VRT safely via explicit dataset handle
    ds_src = osgeo.gdal.Open(outname, osgeo.gdal.GA_ReadOnly)
    if ds_src is not None:
        ds_trans = osgeo.gdal.Translate(outname + '.vrt', ds_src, format='VRT')
        if ds_trans is not None:
            ds_trans.FlushCache()
        ds_trans = None
        ds_src = None

    return


def merged_productbbox(
        metadata_dict, product_dict, workdir='./', bbox_file=None,
        croptounion=False, num_threads='2', minimumOverlap=0.0081,
        verbose=None, runlog=None):
    """
    Extract/merge productBoundingBox layers for each pair.
    Also update dict, report common track bbox
    (default is to take common intersection, but user may specify union),
    report common track union to accurately interpolate metadata fields,
    and expected shape for DEM.
    """
    os.makedirs(workdir, exist_ok=True)

    is_nisar_file = False
    track_fileext = product_dict[0]['unwrappedPhase'][0].split('"')[1]
    if track_fileext.endswith('.h5'):
        is_nisar_file = True

    lyr_proj = int(metadata_dict[0]['projection'][0])
    if bbox_file is not None:
        user_bbox = ARIAtools.util.shp.open_shp(bbox_file)
        overlap_area = ARIAtools.util.shp.shp_area(user_bbox, lyr_proj)
        if overlap_area < minimumOverlap:
            raise Exception(f"User bound box {bbox_file} has an area of only "
                            f"{overlap_area}km\u00b2, below specified "
                            f"minimum threshold area "
                            f"{minimumOverlap}km\u00b2")

    prods_TOTbbox = os.path.join(workdir, 'productBoundingBox.json')
    prods_TOTbbox_metadatalyr = os.path.join(
        workdir, 'productBoundingBox_croptounion_formetadatalyr.json')
    if os.path.exists(prods_TOTbbox) \
        and os.path.exists(prods_TOTbbox_metadatalyr):
        exist_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        exist_metadatalyr = \
            ARIAtools.util.shp.open_shp(prods_TOTbbox_metadatalyr)

        if runlog:
            log_data = runlog.load()
            run_time = log_data['run_times'][-1]
        else:
            run_time = datetime.datetime.now().strftime('%Y%m%d-%H%M%S')

        copy_ext = f"{run_time}.json"
        bbox_copyname = prods_TOTbbox.replace('.json', copy_ext)
        shutil.copyfile(prods_TOTbbox, bbox_copyname)

        metadatalyr_copyname = \
            prods_TOTbbox_metadatalyr.replace('.json', copy_ext)
        shutil.copyfile(prods_TOTbbox_metadatalyr, metadatalyr_copyname)

        LOGGER.debug(
            'Copying existing productBoundingBox to %s', bbox_copyname)
        LOGGER.debug(
            'Copying existing metadatalyr to %s', metadatalyr_copyname)

    else:
        exist_bbox = None

    for scene in product_dict:

        pair_name = scene["pair_name"][0]
        outname = os.path.join(workdir, pair_name + '.json')
        if os.path.exists(outname):
            os.remove(outname)

        for prods_bbox in scene["productBoundingBox"]:
            if os.path.exists(outname):
                union_bbox = ARIAtools.util.shp.open_shp(outname)
                prods_bbox = prods_bbox.union(union_bbox)
            ARIAtools.util.shp.save_shp(
                outname, prods_bbox, lyr_proj)
        scene["productBoundingBox"] = [outname]

    prods_TOTbbox = os.path.join(workdir, 'productBoundingBox.json')

    sceneareas = [
        ARIAtools.util.shp.open_shp(i['productBoundingBox'][0]).area
        for i in product_dict]
    ind_max_area = sceneareas.index(max(sceneareas))
    product_bbox = product_dict[ind_max_area]['productBoundingBox'][0]
    ARIAtools.util.shp.save_shp(
        prods_TOTbbox_metadatalyr, ARIAtools.util.shp.open_shp(product_bbox),
        lyr_proj)

    if bbox_file is not None:
        ARIAtools.util.shp.save_shp(
            prods_TOTbbox, ARIAtools.util.shp.open_shp(bbox_file),
            lyr_proj)

    else:
        ARIAtools.util.shp.save_shp(
            prods_TOTbbox, ARIAtools.util.shp.open_shp(product_bbox),
            lyr_proj)

    rejected_scenes = []
    for scene in product_dict:
        scene_obj = scene['productBoundingBox'][0]
        prods_bbox = ARIAtools.util.shp.open_shp(scene_obj)
        total_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        total_bbox_metadatalyr = ARIAtools.util.shp.open_shp(
            prods_TOTbbox_metadatalyr)

        if croptounion:
            total_bbox = total_bbox.union(prods_bbox)
            total_bbox_metadatalyr = total_bbox_metadatalyr.union(prods_bbox)

            ARIAtools.util.shp.save_shp(
                prods_TOTbbox, total_bbox, lyr_proj)
            ARIAtools.util.shp.save_shp(
                prods_TOTbbox_metadatalyr, total_bbox_metadatalyr,
                lyr_proj)

        else:
            prods_bbox = prods_bbox.intersection(total_bbox)

            if prods_bbox.geom_type == 'MultiPolygon':
                LOGGER.debug(
                    f'Rejected scene {scene_obj} is type MultiPolygon')
                rejected_scenes.append(product_dict.index(scene))
                os.remove(scene_obj)
                continue

            if prods_bbox.bounds == () or prods_bbox.is_empty:
                LOGGER.debug(f'Rejected scene {scene_obj} '
                             f'has no common overlap with bbox')
                rejected_scenes.append(product_dict.index(scene))
                os.remove(scene_obj)

            else:
                overlap_area = ARIAtools.util.shp.shp_area(
                    prods_bbox, lyr_proj)

                if overlap_area < minimumOverlap:
                    LOGGER.debug(f'Rejected scene {scene_obj} has only '
                                 f'{overlap_area}km\u00b2 overlap with bbox')
                    rejected_scenes.append(product_dict.index(scene))
                    os.remove(scene_obj)

                else:
                    ARIAtools.util.shp.save_shp(
                        prods_TOTbbox, prods_bbox, lyr_proj)

                    total_bbox_metadatalyr = total_bbox_metadatalyr.union(
                        ARIAtools.util.shp.open_shp(
                            scene['productBoundingBox'][0]))
                    ARIAtools.util.shp.save_shp(
                        prods_TOTbbox_metadatalyr, total_bbox_metadatalyr,
                        lyr_proj)

    if rejected_scenes != []:
        LOGGER.warning(("%d out of %d interferograms rejected for not "
                        "meeting specified spatial thresholds"),
                        len(rejected_scenes), len(product_dict))
    metadata_dict = [
        i for j, i in enumerate(metadata_dict) if j not in rejected_scenes]
    product_dict = [
        i for j, i in enumerate(product_dict) if j not in rejected_scenes]
    if product_dict == []:
        raise Exception(
            'No common track overlap, footprints cannot be generated.')

    if bbox_file is not None:
        user_bbox = ARIAtools.util.shp.open_shp(bbox_file)
        total_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        user_bbox = user_bbox.intersection(total_bbox)
        ARIAtools.util.shp.save_shp(
            prods_TOTbbox, user_bbox, lyr_proj)

    else:
        bbox_file = prods_TOTbbox

    if exist_bbox and runlog:
        exist_area = ARIAtools.util.shp.shp_area(exist_bbox, lyr_proj)

        new_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        new_area = ARIAtools.util.shp.shp_area(new_bbox, lyr_proj)
        area_ratio = new_area / exist_area

        olap_bbox = exist_bbox.intersection(new_bbox)
        olap_area = ARIAtools.util.shp.shp_area(olap_bbox, lyr_proj)

        delta_area = np.abs(olap_area - exist_area)
        delta_area = np.round(delta_area * 1E7) * 1E-7

        olap_ratio = olap_area / exist_area
        olap_ratio = np.round(olap_ratio * 1E7) * 1E-7

        LOGGER.debug(
            'Area difference (|prev - new|): %f km\u00b2', delta_area)
        LOGGER.debug('Area ratio (new/prev): %f', area_ratio)
        LOGGER.debug('Overlap ratio (new/prev): %f', olap_ratio)

        if (delta_area != 0.0) or (olap_ratio != 1.0):
            LOGGER.warning('Product bbox changed in size from previous '
                           'run %f vs %f', new_area, exist_area)

        if shapely.equals(new_bbox, exist_bbox):
            update_mode = 'skip'
        elif olap_ratio < 1.0:
            update_mode = 'crop_only'
        else:
            update_mode = 'full_extract'
        runlog.update('update_mode', update_mode)

        LOGGER.info('Update mode: %s', update_mode)

    OG_bounds = list(
        ARIAtools.util.shp.open_shp(bbox_file).bounds)
    gdal_warp_kwargs = {
        'format': 'MEM', 'multithread': False, 'dstSRS': f'EPSG:{lyr_proj}'}
    ds_vrt = osgeo.gdal.BuildVRT('', [product_dict[0]['unwrappedPhase'][0]])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', ds_vrt, options=warp_options)
        arrres = [abs(ds.GetGeoTransform()[1]),
                  abs(ds.GetGeoTransform()[-1])]
        ds = None
        ds_vrt = None

    for i, res in enumerate(arrres):
        res_ndx = np.argmin([
            np.abs(res - px_size) for px_size in ARIA_PX_SIZES
        ])
        arrres[i] = ARIA_PX_SIZES[res_ndx]

    gdal_warp_kwargs['outputBounds'] = OG_bounds
    gdal_warp_kwargs['xRes'] = arrres[0]
    gdal_warp_kwargs['yRes'] = arrres[1]
    gdal_warp_kwargs['targetAlignedPixels'] = True
    ds_vrt = osgeo.gdal.BuildVRT('', [product_dict[0]['unwrappedPhase'][0]])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', ds_vrt, options=warp_options)

        arrshape = [ds.RasterYSize, ds.RasterXSize]
        ds_gt = ds.GetGeoTransform()
        new_bounds = [ds_gt[0], ds_gt[3] + (ds_gt[-1] * arrshape[0]),
                      ds_gt[0] + (ds_gt[1] * arrshape[1]), ds_gt[3]]

        if OG_bounds != new_bounds:
            user_bbox = shapely.geometry.Polygon(np.column_stack((
                np.array([new_bounds[0], new_bounds[2], new_bounds[2],
                          new_bounds[0], new_bounds[0]]),
                np.array([new_bounds[1], new_bounds[1], new_bounds[3],
                          new_bounds[3], new_bounds[1]]))))

            bbox_file = os.path.join(
                os.path.dirname(workdir), 'user_bbox.json')
            ARIAtools.util.shp.save_shp(
                bbox_file, user_bbox, lyr_proj)
            total_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
            user_bbox = user_bbox.intersection(total_bbox)
            ARIAtools.util.shp.save_shp(
                prods_TOTbbox, user_bbox, lyr_proj)

        proj = ds.GetProjection()
        ds = None
        ds_vrt = None

    if runlog is None:
        update_mode = 'full_extract'
    else:
        log_data = runlog.load()
        update_mode = log_data['update_mode']

        if ('arrres' in log_data.keys()) \
                and (arrres != log_data['arrres']):
            runlog.update('update_mode', 'full_extract')
            LOGGER.warning('arrres has changed. '
                           'Setting update mode to full_extract.')

        if ('lyr_proj' in log_data.keys()) \
                and (lyr_proj != log_data['lyr_proj']):
            runlog.update('update_mode', 'full_extract')
            LOGGER.warning('lyr_proj has changed. '
                           'Setting update mode to full_extract.')

        runlog.update('croptounion', croptounion)
        runlog.update('prods_TOTbbox', prods_TOTbbox)
        runlog.update('prods_TOTbbox_metadatalyr', prods_TOTbbox_metadatalyr)
        runlog.update('arrres', arrres)
        runlog.update('lyr_proj', lyr_proj)
        runlog.update('metadata_dict', metadata_dict)
        runlog.update('product_dict', product_dict)
        runlog.update('bbox_file', bbox_file)
        runlog.update('is_nisar_file', is_nisar_file)

    return (metadata_dict, product_dict, bbox_file, prods_TOTbbox,
            prods_TOTbbox_metadatalyr, arrres, proj, update_mode,
            is_nisar_file)


def create_raster_from_gunw(fname, data_lis, proj, driver, hgt_field=None,
    sign_multiplier=1, dem=None):
    """Create a local raster from a GUNW subdataset.

    NISAR HDF5 layers use GeoTIFF as their on-disk backing raster.  ARIA-tools
    accesses these extensionless files through GDAL/VRT, so the backing driver
    does not need to be ENVI.  Avoiding ENVI is necessary on systems where its
    scanline writer leaves partial 3-D correction rasters.
    """

    def raster_is_complete(path, expected_bands=None,
                           expected_ny=None, expected_nx=None):
        """Confirm the raster opens and its first/last scanlines are readable."""
        check_ds = osgeo.gdal.Open(path, osgeo.gdal.GA_ReadOnly)
        if check_ds is None:
            return False
        try:
            if expected_bands is not None and \
                    check_ds.RasterCount != expected_bands:
                return False
            if expected_ny is not None and check_ds.RasterYSize != expected_ny:
                return False
            if expected_nx is not None and check_ds.RasterXSize != expected_nx:
                return False
            if check_ds.RasterCount < 1 or check_ds.RasterYSize < 1 or \
                    check_ds.RasterXSize < 1:
                return False

            for band_no in range(1, check_ds.RasterCount + 1):
                check_band = check_ds.GetRasterBand(band_no)
                first = check_band.ReadRaster(
                    0, 0, check_ds.RasterXSize, 1)
                last = check_band.ReadRaster(
                    0, check_ds.RasterYSize - 1,
                    check_ds.RasterXSize, 1)
                if first is None or last is None:
                    return False
            return True
        except Exception:
            return False
        finally:
            check_ds = None

    # Metadata lists may contain one subdataset per spatial frame. Exporting
    # only data_lis[0] leaves the rest of the DEM-intersected product empty.
    # Build and validate each frame cube first, then mosaic all bands on a
    # common grid before finalize_metadata performs the height interpolation.
    if len(data_lis) > 1:
        frame_names = [
            f'{fname}_frame_{index:03d}_stage'
            for index in range(len(data_lis))]
        frame_vrts = []
        cleanup_suffixes = ['', '.vrt', '.aux.xml', '.xml', '.hdr']

        for suffix in cleanup_suffixes:
            final_path = fname + suffix
            if os.path.isfile(final_path):
                os.remove(final_path)

        try:
            expected_bands = None
            common_height_meta = None
            frame_info = []

            for index, (frame_name, frame_source) in enumerate(
                    zip(frame_names, data_lis)):
                create_raster_from_gunw(
                    frame_name, [frame_source], proj,
                    METADATA_STAGING_DRIVER, hgt_field,
                    sign_multiplier=sign_multiplier, dem=dem)
                frame_vrt = frame_name + '.vrt'
                if not raster_is_complete(frame_vrt):
                    raise RuntimeError(
                        f'Metadata frame {index + 1}/{len(data_lis)} did not '
                        f'produce a complete raster: {frame_source}')

                ds_frame = osgeo.gdal.Open(
                    frame_vrt, osgeo.gdal.GA_ReadOnly)
                if ds_frame is None:
                    raise RuntimeError(
                        f'Could not reopen metadata frame: {frame_vrt}')

                band_count = ds_frame.RasterCount
                if expected_bands is None:
                    expected_bands = band_count
                elif band_count != expected_bands:
                    ds_frame = None
                    raise RuntimeError(
                        f'Metadata frame band mismatch for {fname}: '
                        f'expected {expected_bands}, frame {index + 1} has '
                        f'{band_count}')

                frame_gt = ds_frame.GetGeoTransform()
                frame_proj = ds_frame.GetProjection()
                x_edge_2 = (frame_gt[0]
                            + frame_gt[1] * ds_frame.RasterXSize)
                y_edge_2 = (frame_gt[3]
                            + frame_gt[5] * ds_frame.RasterYSize)
                frame_bounds = (
                    min(frame_gt[0], x_edge_2),
                    min(frame_gt[3], y_edge_2),
                    max(frame_gt[0], x_edge_2),
                    max(frame_gt[3], y_edge_2))
                frame_info.append({
                    'projection': frame_proj,
                    'geotransform': frame_gt,
                    'bounds': frame_bounds})
                ds_frame = None

                if hgt_field is not None:
                    frame_height_meta = ARIAtools.util.vrt.get_hgt_meta(
                        frame_vrt, hgt_field)
                    if not frame_height_meta:
                        raise RuntimeError(
                            f'Missing height metadata for frame '
                            f'{index + 1}: {frame_source}')
                    if common_height_meta is None:
                        common_height_meta = frame_height_meta
                    else:
                        reference_heights = np.array(
                            common_height_meta[1:-1].split(','),
                            dtype='float64')
                        frame_heights = np.array(
                            frame_height_meta[1:-1].split(','),
                            dtype='float64')
                        if (reference_heights.shape != frame_heights.shape
                                or not np.allclose(
                                    reference_heights, frame_heights,
                                    rtol=0.0, atol=1e-4)):
                            raise RuntimeError(
                                f'Inconsistent height levels across metadata '
                                f'frames for {fname}')

                frame_vrts.append(frame_vrt)

            first_projection = frame_info[0]['projection']
            same_projection = bool(first_projection)
            if same_projection:
                first_crs = pyproj.CRS.from_wkt(first_projection)
                same_projection = all(
                    info['projection']
                    and first_crs.equals(
                        pyproj.CRS.from_wkt(info['projection']))
                    for info in frame_info[1:])

            warp_kwargs = {
                'format': METADATA_STAGING_DRIVER,
                'outputType': osgeo.gdal.GDT_Float32,
                # Let GDAL honor each frame's intrinsic NoData value; NISAR
                # revisions may encode it differently between datasets.
                'dstNodata': np.nan,
                'multithread': False,
                'creationOptions': [
                    'TILED=YES', 'COMPRESS=DEFLATE',
                    'PREDICTOR=3', 'BIGTIFF=IF_SAFER']}

            if same_projection:
                first_gt = frame_info[0]['geotransform']
                x_res = abs(first_gt[1])
                y_res = abs(first_gt[5])
                for info in frame_info[1:]:
                    frame_gt = info['geotransform']
                    if (not np.isclose(abs(frame_gt[1]), x_res)
                            or not np.isclose(abs(frame_gt[5]), y_res)):
                        raise RuntimeError(
                            f'Metadata frame resolution mismatch for {fname}')
                warp_kwargs.update({
                    'dstSRS': first_projection,
                    'outputBounds': (
                        min(info['bounds'][0] for info in frame_info),
                        min(info['bounds'][1] for info in frame_info),
                        max(info['bounds'][2] for info in frame_info),
                        max(info['bounds'][3] for info in frame_info)),
                    'xRes': x_res, 'yRes': y_res,
                    'targetAlignedPixels': True})
            else:
                # Frames crossing a projection-zone boundary are mosaicked on
                # the DEM grid CRS, which is also finalize_metadata's target.
                mosaic_projection = (
                    dem.GetProjection() if dem is not None else proj)
                if isinstance(mosaic_projection, int):
                    mosaic_projection = f'EPSG:{mosaic_projection}'
                warp_kwargs['dstSRS'] = mosaic_projection

            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": "1"}):
                ds_mosaic = osgeo.gdal.Warp(
                    fname, frame_vrts,
                    options=osgeo.gdal.WarpOptions(**warp_kwargs))
            if ds_mosaic is None:
                raise RuntimeError(
                    f'Failed to mosaic {len(frame_vrts)} metadata frames '
                    f'for {fname}')
            ds_mosaic.FlushCache()
            mosaic_gt = ds_mosaic.GetGeoTransform()
            mosaic_proj = ds_mosaic.GetProjection()
            mosaic_width = ds_mosaic.RasterXSize
            mosaic_height = ds_mosaic.RasterYSize
            ds_mosaic = None

            if not raster_is_complete(
                    fname, expected_bands=expected_bands,
                    expected_ny=mosaic_height, expected_nx=mosaic_width):
                raise RuntimeError(
                    f'Mosaicked metadata raster is incomplete: {fname}')

            ds_mosaic_src = osgeo.gdal.Open(fname, osgeo.gdal.GA_ReadOnly)
            vrt_options = osgeo.gdal.BuildVRTOptions(
                outputSRS=mosaic_proj or proj)
            ds_mosaic_vrt = osgeo.gdal.BuildVRT(
                fname + '.vrt', [ds_mosaic_src], options=vrt_options)
            if ds_mosaic_vrt is None:
                ds_mosaic_src = None
                raise RuntimeError(
                    f'Could not create VRT for metadata mosaic: {fname}')
            if hgt_field is not None and common_height_meta is not None:
                ds_mosaic_vrt.SetMetadataItem(
                    hgt_field, common_height_meta)
            ds_mosaic_vrt.SetMetadataItem(
                'ARIA_FRAME_COUNT', str(len(data_lis)))
            ds_mosaic_vrt.FlushCache()
            ds_mosaic_vrt = None
            ds_mosaic_src = None

            x_edge_2 = mosaic_gt[0] + mosaic_gt[1] * mosaic_width
            y_edge_2 = mosaic_gt[3] + mosaic_gt[5] * mosaic_height
            LOGGER.info(
                'Mosaicked %d metadata frames for %s: %d bands, %dx%d, '
                'bounds=(%.3f, %.3f, %.3f, %.3f)',
                len(frame_vrts), fname, expected_bands,
                mosaic_width, mosaic_height,
                min(mosaic_gt[0], x_edge_2),
                min(mosaic_gt[3], y_edge_2),
                max(mosaic_gt[0], x_edge_2),
                max(mosaic_gt[3], y_edge_2))
            return
        except Exception:
            for suffix in cleanup_suffixes:
                failed_path = fname + suffix
                if os.path.isfile(failed_path):
                    os.remove(failed_path)
            raise
        finally:
            for frame_name in frame_names:
                for suffix in cleanup_suffixes:
                    frame_path = frame_name + suffix
                    if os.path.isfile(frame_path):
                        os.remove(frame_path)

    # Clean up any leftover temp or corrupted files
    for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr', '.tmp.tif', '_warp.tif', '_ref.vrt']:
        f_to_rm = f"{fname}{ext}"
        if os.path.exists(f_to_rm):
            try:
                os.remove(f_to_rm)
            except OSError:
                pass

    src_path = data_lis[0]
    is_h5 = '.h5' in src_path
    translated_ok = False

    if is_h5:
        export_stage = 'parsing the HDF5 subdataset path'
        try:
            parts = src_path.split('":')
            raw_path = parts[0].replace('NETCDF:"', '').replace('HDF5:"', '')
            ds_path = parts[1] if len(parts) > 1 else ''

            if os.path.exists(raw_path):
                # GDAL can read the subdataset header even on systems where
                # translating all of its scanlines fails.  Preserve its grid
                # and metadata on the directly written output.
                export_stage = 'reading the GDAL subdataset header'
                src_ds = osgeo.gdal.Open(src_path, osgeo.gdal.GA_ReadOnly)
                src_gt = src_ds.GetGeoTransform() if src_ds else None
                src_proj = src_ds.GetProjection() if src_ds else None
                src_meta = src_ds.GetMetadata() if src_ds else {}
                src_band_count = src_ds.RasterCount if src_ds else None
                src_ny = src_ds.RasterYSize if src_ds else None
                src_nx = src_ds.RasterXSize if src_ds else None
                src_nodata = None
                if src_ds and src_ds.RasterCount:
                    src_nodata = src_ds.GetRasterBand(1).GetNoDataValue()
                src_ds = None

                export_stage = 'opening the HDF5 file'
                with h5py.File(raw_path, 'r') as h5f:
                    if ds_path not in h5f:
                        raise KeyError(
                            f'HDF5 dataset {ds_path!r} not found in {raw_path}')

                    h5_ds = h5f[ds_path]
                    if h5_ds.ndim not in (2, 3):
                        raise ValueError(
                            f'Unsupported shape {h5_ds.shape} for {ds_path}')

                    if h5_ds.ndim == 2:
                        bands, ny, nx = 1, h5_ds.shape[0], h5_ds.shape[1]
                        band_axis = None
                    else:
                        # GDAL exposes a 3-D NetCDF/HDF5 variable as one
                        # raster band per value of its non-spatial dimension.
                        # The NISAR radarGrid cubes normally use axis 0, but
                        # infer it from the GDAL header when possible.
                        gdal_band_count = src_band_count or h5_ds.shape[0]
                        gdal_ny = src_ny or h5_ds.shape[-2]
                        gdal_nx = src_nx or h5_ds.shape[-1]

                        candidates = [
                            axis for axis, size in enumerate(h5_ds.shape)
                            if size == gdal_band_count
                            and tuple(h5_ds.shape[i] for i in range(3)
                                      if i != axis) == (gdal_ny, gdal_nx)
                        ]
                        band_axis = candidates[0] if candidates else 0
                        bands = h5_ds.shape[band_axis]
                        spatial_shape = tuple(
                            h5_ds.shape[i] for i in range(3)
                            if i != band_axis)
                        ny, nx = spatial_shape

                    # Do not rely on the geotransform reported by GDAL for
                    # NISAR radarGrid HDF5 subdatasets. Some GDAL builds
                    # expose these arrays with an identity transform even
                    # though sibling datasets contain projected pixel centers.
                    grid_group = os.path.dirname(ds_path.rstrip('/'))
                    x_path = f'{grid_group}/xCoordinates'
                    y_path = f'{grid_group}/yCoordinates'
                    projection_path = f'{grid_group}/projection'
                    coordinate_georef = False
                    grid_epsg = None

                    if x_path in h5f and y_path in h5f:
                        x_coords = np.asarray(h5f[x_path][()]).squeeze()
                        y_coords = np.asarray(h5f[y_path][()]).squeeze()
                        if (x_coords.ndim == 1 and y_coords.ndim == 1
                                and len(x_coords) == nx
                                and len(y_coords) == ny
                                and nx > 1 and ny > 1):
                            dx = float(np.median(np.diff(x_coords)))
                            dy = float(np.median(np.diff(y_coords)))
                            if (np.isfinite(dx) and np.isfinite(dy)
                                    and dx != 0.0 and dy != 0.0):
                                # Coordinate arrays describe pixel centers;
                                # GDAL geotransforms begin at an outer corner.
                                src_gt = (
                                    float(x_coords[0]) - dx / 2.0, dx, 0.0,
                                    float(y_coords[0]) - dy / 2.0, 0.0, dy)
                                coordinate_georef = True
                        else:
                            LOGGER.warning(
                                'Ignoring incompatible NISAR coordinate '
                                'arrays for %s: x=%s, y=%s, raster=%dx%d',
                                ds_path, x_coords.shape, y_coords.shape,
                                nx, ny)

                    if projection_path in h5f:
                        try:
                            projection_value = np.asarray(
                                h5f[projection_path][()]).squeeze()
                            if projection_value.size == 1:
                                grid_epsg = int(projection_value.item())
                                src_proj = pyproj.CRS.from_epsg(
                                    grid_epsg).to_wkt()
                        except (TypeError, ValueError,
                                pyproj.exceptions.CRSError):
                            LOGGER.warning(
                                'Could not interpret NISAR projection at %s',
                                projection_path)

                    if coordinate_georef:
                        x_edge_2 = src_gt[0] + src_gt[1] * nx
                        y_edge_2 = src_gt[3] + src_gt[5] * ny
                        LOGGER.info(
                            'Using NISAR coordinate georeferencing for %s: '
                            'EPSG=%s, bounds=(%.3f, %.3f, %.3f, %.3f), '
                            'resolution=(%.6g, %.6g)',
                            ds_path, grid_epsg,
                            min(src_gt[0], x_edge_2),
                            min(src_gt[3], y_edge_2),
                            max(src_gt[0], x_edge_2),
                            max(src_gt[3], y_edge_2),
                            abs(src_gt[1]), abs(src_gt[5]))
                    else:
                        LOGGER.warning(
                            'NISAR coordinate georeferencing was unavailable '
                            'for %s; retaining the GDAL subdataset transform %s',
                            ds_path, src_gt)

                    # Direct multi-band writes through ENVI/ISCE can fail
                    # partway through a scanline. This raster is an internal
                    # cube consumed by finalize_metadata, not the published
                    # product, so always use a robust GeoTIFF backing store.
                    storage_driver = METADATA_STAGING_DRIVER
                    out_driver = osgeo.gdal.GetDriverByName(storage_driver)
                    if out_driver is None:
                        raise RuntimeError(
                            f'GDAL driver {storage_driver!r} not found')

                    create_options = []
                    if storage_driver == 'GTiff':
                        create_options = [
                            'TILED=YES', 'COMPRESS=DEFLATE',
                            'PREDICTOR=3', 'BIGTIFF=IF_SAFER']
                    export_stage = f'creating the {storage_driver} raster'
                    out_ds = out_driver.Create(
                        fname, nx, ny, bands, osgeo.gdal.GDT_Float32,
                        options=create_options)
                    if out_ds is None:
                        raise RuntimeError(
                            f'GDAL could not create {fname} with '
                            f'{storage_driver}')

                    if src_gt:
                        out_ds.SetGeoTransform(src_gt)
                    if src_proj or proj:
                        out_ds.SetProjection(src_proj or proj)
                    if src_meta:
                        out_ds.SetMetadata(src_meta)

                    try:
                        for b in range(bands):
                            export_stage = (
                                f'reading HDF5 band {b + 1} of {bands}')
                            if band_axis is None:
                                band_data = np.asarray(h5_ds[:, :],
                                                       dtype=np.float32)
                            else:
                                index = [slice(None)] * 3
                                index[band_axis] = b
                                band_data = np.asarray(
                                    h5_ds[tuple(index)], dtype=np.float32)

                            if sign_multiplier == -1:
                                band_data *= -1.0

                            export_stage = (
                                f'writing {storage_driver} band '
                                f'{b + 1} of {bands}')
                            out_band = out_ds.GetRasterBand(b + 1)
                            out_band.WriteArray(band_data)
                            out_band.SetNoDataValue(
                                src_nodata if src_nodata is not None
                                else np.nan)
                            out_band.FlushCache()
                    except Exception:
                        out_ds = None
                        raise

                    out_ds.FlushCache()
                    out_ds = None

                    export_stage = 'validating the completed raster'
                    translated_ok = raster_is_complete(
                        fname, expected_bands=bands,
                        expected_ny=ny, expected_nx=nx)
                    if not translated_ok:
                        raise RuntimeError(
                            f'Raster validation failed after writing {fname}')

                    LOGGER.info(
                        'Direct HDF5 export completed: %s '
                        '(%s, %d bands, %dx%d)',
                        fname, storage_driver, bands, nx, ny)
        except Exception as e:
            LOGGER.warning(
                'Direct HDF5 export failed while %s for %s: %s',
                export_stage, src_path, e)
            for ext in ['', '.aux.xml', '.xml', '.hdr']:
                partial = f'{fname}{ext}'
                if os.path.exists(partial):
                    try:
                        os.remove(partial)
                    except OSError:
                        pass

    if not translated_ok:
        try:
            # Any cube carrying a height dimension will be consumed by
            # finalize_metadata. Keep it in the staging driver regardless of
            # whether its source container is HDF5 or NetCDF.
            fallback_driver = (
                METADATA_STAGING_DRIVER
                if is_h5 or hgt_field is not None else driver)
            ds_trans = osgeo.gdal.Translate(
                fname, src_path, format=fallback_driver)
            if ds_trans is None:
                raise RuntimeError(
                    f'GDAL Translate returned no dataset for {src_path}')
            ds_trans.FlushCache()
            ds_trans = None
            translated_ok = raster_is_complete(fname)
            if not translated_ok:
                raise RuntimeError(
                    f'Raster validation failed after translating {fname}')
        except Exception as e:
            LOGGER.warning(f"GDAL Translate failed for {src_path}: {e}")
            for ext in ['', '.aux.xml', '.xml', '.hdr']:
                partial = f'{fname}{ext}'
                if os.path.exists(partial):
                    try:
                        os.remove(partial)
                    except OSError:
                        pass
            return

    if not os.path.exists(fname) or os.path.getsize(fname) < 100:
        LOGGER.warning(f"Output raster {fname} was not created or is empty.")
        return

    # Create VRT dataset
    ds_src = osgeo.gdal.Open(fname, osgeo.gdal.GA_ReadOnly)
    if ds_src is None:
        LOGGER.warning('Completed raster could not be reopened: %s', fname)
        return

    # outputSRS assigns a CRS; it does not transform coordinates. Preserve the
    # projected CRS read from NISAR rather than relabelling the grid.
    raster_proj = ds_src.GetProjection()
    buildvrt_options = osgeo.gdal.BuildVRTOptions(
        outputSRS=raster_proj or proj)
    ds_vrt = osgeo.gdal.BuildVRT(
        fname + '.vrt', [ds_src], options=buildvrt_options
    )
    if ds_vrt is None:
        LOGGER.warning('Could not create VRT for completed raster: %s', fname)
        ds_src = None
        return

    ds_vrt.SetMetadataItem('ARIA_FRAME_COUNT', '1')

    ds_vrt.FlushCache()
    ds_vrt = None
    ds_src = None

    if not raster_is_complete(fname + '.vrt'):
        LOGGER.warning('VRT validation failed for completed raster: %s', fname)
        try:
            os.remove(fname + '.vrt')
        except OSError:
            pass
        return

    # Attach height metadata to VRT if present
    if hgt_field is not None and os.path.exists(fname + '.vrt'):
        hgt_meta = ARIAtools.util.vrt.get_hgt_meta(src_path, hgt_field)
        if not hgt_meta and is_h5:
            try:
                parts = src_path.split('":')
                raw_path = parts[0].replace('NETCDF:"', '').replace('HDF5:"', '')
                if os.path.exists(raw_path):
                    with h5py.File(raw_path, 'r') as h5f:
                        h_path = '/science/LSAR/GUNW/metadata/radarGrid/heightAboveEllipsoid'
                        if h_path in h5f:
                            h_vals = h5f[h_path][:]
                            hgt_meta = '{' + ','.join(str(h) for h in h_vals) + '}'
            except Exception:
                pass

        if hgt_meta:
            ds_meta_update = osgeo.gdal.Open(
                fname + '.vrt', osgeo.gdal.GA_Update
            )
            if ds_meta_update is not None:
                ds_meta_update.SetMetadataItem(hgt_field, hgt_meta)
                ds_meta_update.FlushCache()
            ds_meta_update = None

    return


def validate_metadata_staging_raster(
        path, context, expected_bands=None,
        expected_width=None, expected_height=None):
    """Validate a GeoTIFF metadata intermediate before finalization."""
    ds_check = osgeo.gdal.Open(path, osgeo.gdal.GA_ReadOnly)
    if ds_check is None:
        raise RuntimeError(f'Could not open {context} staging raster: {path}')
    try:
        actual_driver = ds_check.GetDriver().ShortName
        if actual_driver.upper() != METADATA_STAGING_DRIVER.upper():
            raise RuntimeError(
                f'{context} staging raster uses {actual_driver}, expected '
                f'{METADATA_STAGING_DRIVER}: {path}')
        if (expected_bands is not None
                and ds_check.RasterCount != expected_bands):
            raise RuntimeError(
                f'{context} staging raster has {ds_check.RasterCount} bands; '
                f'expected {expected_bands}: {path}')
        if (expected_width is not None
                and ds_check.RasterXSize != expected_width):
            raise RuntimeError(
                f'{context} staging raster width is {ds_check.RasterXSize}; '
                f'expected {expected_width}: {path}')
        if (expected_height is not None
                and ds_check.RasterYSize != expected_height):
            raise RuntimeError(
                f'{context} staging raster height is {ds_check.RasterYSize}; '
                f'expected {expected_height}: {path}')
        for band_index in range(1, ds_check.RasterCount + 1):
            band = ds_check.GetRasterBand(band_index)
            first_row = band.ReadRaster(0, 0, ds_check.RasterXSize, 1)
            last_row = band.ReadRaster(
                0, ds_check.RasterYSize - 1, ds_check.RasterXSize, 1)
            if first_row is None or last_row is None:
                raise RuntimeError(
                    f'{context} staging raster band {band_index} is '
                    f'incomplete: {path}')
    finally:
        ds_check = None


def prep_metadatalayers(
        outname, metadata_arr, dem, layer, layers, is_nisar_file=False,
        proj='4326', driver='ENVI', model_name=None, sign_multiplier=1):
    """Wrapper to prep metadata layer for extraction"""
    if dem is None:
        raise Exception('No DEM input specified. '
                        'Cannot extract 3D imaging geometry '
                        'layers without DEM to intersect with.')

    ifg = os.path.basename(outname)
    out_dir = os.path.dirname(outname)
    ref_outname = copy.deepcopy(outname)

    if metadata_arr[0].split('/')[-1] == 'ionosphere':
        ds_vrt = osgeo.gdal.BuildVRT(outname + '.vrt', metadata_arr)
        if ds_vrt is not None:
            ds_vrt.FlushCache()
        ds_vrt= None
        return [0], None, outname

    if 'tropo' in layer:
        if not is_nisar_file:
            out_dir = os.path.join(out_dir, model_name)
        outname = os.path.join(out_dir, ifg)

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    ds_meta = osgeo.gdal.Open(metadata_arr[0])
    zdim = ds_meta.GetMetadataItem('NETCDF_DIM_EXTRA')[1:-1]
    ds_meta = None
    hgt_field = f'NETCDF_DIM_{zdim}_VALUES'

    # A run made with the older azimuth path may have left a VRT whose UTM
    # coordinates were merely labelled EPSG:4326. Detect that impossible
    # combination so a retry regenerates the derived raster instead of
    # reusing the poisoned cache entry.
    existing_vrt = outname + '.vrt'
    if is_nisar_file and os.path.exists(existing_vrt):
        ds_existing = osgeo.gdal.Open(
            existing_vrt, osgeo.gdal.GA_ReadOnly)
        invalid_georef = False
        if ds_existing is not None:
            existing_gt = ds_existing.GetGeoTransform()
            existing_proj = ds_existing.GetProjection()
            try:
                existing_crs = (pyproj.CRS.from_wkt(existing_proj)
                                if existing_proj else None)
                invalid_georef = bool(
                    existing_crs is not None
                    and existing_crs.is_geographic
                    and (abs(existing_gt[0]) > 360.0
                         or abs(existing_gt[3]) > 90.0))
            except pyproj.exceptions.CRSError:
                invalid_georef = True
        ds_existing = None

        if invalid_georef:
            LOGGER.warning(
                'Removing cached NISAR raster with inconsistent CRS and '
                'coordinates: %s', outname)
            for suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                stale_path = outname + suffix
                if os.path.isfile(stale_path):
                    os.remove(stale_path)

    if (
        not os.path.exists(outname + ".vrt")
        and len(
            {
                ARIAtools.util.vrt.get_hgt_meta(i, hgt_field)
                for i in metadata_arr
            }
        )
        != 1
    ):
        heights = [
            ARIAtools.util.vrt.get_hgt_meta(i, hgt_field)
            for i in metadata_arr
        ]
        raise Exception(
            "Inconsistent heights for metadata layer(s) "
            f"{metadata_arr}; corresponding heights: {heights}"
        )

    if 'tropo' in layer or layer == 'solidEarthTide':
        if not is_nisar_file:
            date_dir = os.path.join(out_dir, 'dates')
            if not os.path.exists(date_dir):
                os.mkdir(date_dir)

            ref_outname = os.path.join(date_dir, ifg.split('_')[0])
            sec_outname = os.path.join(date_dir, ifg.split('_')[1])
            ref_str = 'reference/' + layer
            sec_str = 'secondary/' + layer
            
            if ref_str in metadata_arr[0]:
                sec_metadata_arr = [i.replace(ref_str, sec_str) for i in metadata_arr]
            elif 'reference' in metadata_arr[0]:
                sec_metadata_arr = [i.replace('reference', 'secondary') for i in metadata_arr]
            else:
                sec_metadata_arr = copy.deepcopy(metadata_arr)

            tup_outputs = [
                (ref_outname, metadata_arr), (sec_outname, sec_metadata_arr)]

        else:
            ref_outname = os.path.join(out_dir, ifg)
            sec_outname = ref_outname
            tup_outputs = [(ref_outname, metadata_arr)]

        for i in tup_outputs:
            for j in glob.glob(i[0] + '*'):
                if os.path.isfile(j):
                    os.remove(j)
            create_raster_from_gunw(i[0], i[1], proj, driver, hgt_field,
                sign_multiplier, dem=dem)

        if not is_nisar_file:
            generate_diff(
                ref_outname, sec_outname, outname, layer, layer, False,
                hgt_field, proj, driver, dem=dem)

        if layer in layers:
            for i in [ref_outname, sec_outname]:
                if not os.path.exists(i):
                    create_raster_from_gunw(i, [i], proj, driver, hgt_field)

    else:
        if not os.path.exists(outname + '.vrt'):
            if is_nisar_file:
                if layer == 'azimuthAngle':
                    losx_arr = copy.deepcopy(metadata_arr)
                    losy_arr = [
                        path.replace('losUnitVectorX', 'losUnitVectorY')
                        for path in metadata_arr
                    ]

                    losx_name = os.path.join(out_dir, f'temp_{ifg}_losx_arr')
                    losy_name = os.path.join(out_dir, f'temp_{ifg}_losy_arr')

                    create_raster_from_gunw(
                        losx_name, losx_arr, proj, driver, hgt_field,
                        dem=dem
                    )
                    create_raster_from_gunw(
                        losy_name, losy_arr, proj, driver, hgt_field,
                        dem=dem
                    )

                    ds_temp = osgeo.gdal.Open(losx_name + '.vrt')
                    if ds_temp is None:
                        raise RuntimeError(
                            f'Could not open temporary LOS-X raster '
                            f'{losx_name}.vrt')
                    src_nodata = ds_temp.GetRasterBand(1).GetNoDataValue()
                    azimuth_gt = ds_temp.GetGeoTransform()
                    azimuth_proj = ds_temp.GetProjection()
                    azimuth_width = ds_temp.RasterXSize
                    azimuth_height = ds_temp.RasterYSize
                    azimuth_bands = ds_temp.RasterCount
                    ds_temp = None

                    if src_nodata is not None and not np.isnan(src_nodata):
                        calc_cmd = (
                            f"numpy.where("
                            f"(A=={src_nodata})|(B=={src_nodata}), "
                            f"numpy.nan, "
                            f"numpy.degrees(numpy.arctan2(-B, -A)))"
                        )
                    else:
                        calc_cmd = "numpy.degrees(numpy.arctan2(-B, -A))"

                    # Keep NISAR derived rasters GeoTIFF-backed as well. More
                    # importantly, preserve the projected LOS grid on the
                    # calculated azimuth raster; assigning ``proj`` here can
                    # relabel UTM coordinates as longitude/latitude.
                    calc_driver = METADATA_STAGING_DRIVER
                    # A failed ISCE-backed calculation from an earlier run
                    # may leave a truncated raster and XML sidecar. Remove
                    # those before recreating the derived cube in GeoTIFF.
                    for calc_suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                        calc_path = outname + calc_suffix
                        if os.path.isfile(calc_path):
                            os.remove(calc_path)
                    osgeo_utils.gdal_calc.Calc(
                        A=losx_name + '.vrt',
                        B=losy_name + '.vrt',
                        outfile=outname,
                        calc=calc_cmd,
                        format=calc_driver,
                        allBands="A",
                        quiet=True,
                        overwrite=True
                    )

                    ds_update = osgeo.gdal.Open(
                        outname, osgeo.gdal.GA_Update
                    )
                    if ds_update is None:
                        raise RuntimeError(
                            f'Azimuth calculation did not create {outname}')
                    if (ds_update.RasterXSize != azimuth_width
                            or ds_update.RasterYSize != azimuth_height
                            or ds_update.RasterCount != azimuth_bands):
                        raise RuntimeError(
                            f'Unexpected azimuth raster dimensions for '
                            f'{outname}: got {ds_update.RasterCount} bands '
                            f'of {ds_update.RasterXSize}x'
                            f'{ds_update.RasterYSize}; expected '
                            f'{azimuth_bands} bands of '
                            f'{azimuth_width}x{azimuth_height}')
                    ds_update.SetGeoTransform(azimuth_gt)
                    if azimuth_proj:
                        ds_update.SetProjection(azimuth_proj)
                    for b in range(1, ds_update.RasterCount + 1):
                        ds_update.GetRasterBand(b).SetNoDataValue(np.nan)
                    ds_update.FlushCache()
                    ds_update = None
                    validate_metadata_staging_raster(
                        outname, 'azimuthAngle',
                        expected_bands=azimuth_bands,
                        expected_width=azimuth_width,
                        expected_height=azimuth_height)
                    try:
                        azimuth_crs_label = (pyproj.CRS.from_wkt(
                            azimuth_proj).to_string()
                            if azimuth_proj else 'unset')
                    except pyproj.exceptions.CRSError:
                        azimuth_crs_label = 'unrecognized'
                    LOGGER.info(
                        'Derived NISAR azimuthAngle on source grid: '
                        'CRS=%s, transform=%s, bands=%d, size=%dx%d',
                        azimuth_crs_label, azimuth_gt, azimuth_bands,
                        azimuth_width, azimuth_height)

                    buildvrt_options = osgeo.gdal.BuildVRTOptions(
                        outputSRS=azimuth_proj or proj
                    )
                    ds_vrt = osgeo.gdal.BuildVRT(
                        outname + '.vrt',
                        [outname],
                        options=buildvrt_options
                    )
                    if ds_vrt is not None:
                        ds_vrt.SetMetadataItem(
                            'ARIA_FRAME_COUNT', str(len(metadata_arr)))
                        ds_vrt.FlushCache()
                    ds_vrt = None

                    if hgt_field is not None:
                        ds_meta = osgeo.gdal.Open(losx_name + '.vrt')
                        hgt_meta = ds_meta.GetMetadataItem(hgt_field)
                        ds_meta = None

                        ds_vrt = osgeo.gdal.Open(outname + '.vrt', osgeo.gdal.GA_Update)
                        if ds_vrt is not None:
                            ds_vrt.SetMetadataItem(hgt_field, hgt_meta)
                            ds_vrt.FlushCache()
                        ds_vrt = None

                    for f in [losx_name, losy_name]:
                        for junk_file in glob.glob(f"{f}*"):
                            try:
                                os.remove(junk_file)
                            except OSError:
                                pass
                else:
                    create_raster_from_gunw(outname, metadata_arr,
                        proj, driver, hgt_field, sign_multiplier,
                        dem=dem)
            else:
                ds_vrt = osgeo.gdal.BuildVRT(outname + '.vrt', metadata_arr)
                
                ds_src = osgeo.gdal.Open(metadata_arr[0])
                hgt_val = ds_src.GetMetadataItem(hgt_field)
                ds_src = None

                ds_vrt.SetMetadataItem(hgt_field, hgt_val)
                ds_vrt.FlushCache()
                ds_vrt = None

    return hgt_field, ref_outname


def generate_diff(ref_outname, sec_outname, outname, key, OG_key, tropo_total,
                  hgt_field, proj, driver, sign_multiplier=1, dem=None):
    """ Compute differential from reference and secondary scenes (Multi-dim safe) """

    output_dir = os.path.dirname(outname)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    subset_vrts = []
    subset_height_meta = None
    ref_vrt_path = ref_outname + '.vrt'
    sec_vrt_path = sec_outname + '.vrt'

    if (dem is not None and hgt_field
            and not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET')):
        try:
            heightsMeta_str = ARIAtools.util.vrt.get_hgt_meta(
                ref_vrt_path, hgt_field)
            if heightsMeta_str:
                heightsMeta = np.array(
                    heightsMeta_str[1:-1].split(','), dtype='float32')
                dem_min, dem_max = (
                    ARIAtools.util.interp._compute_dem_range(dem))
                band_indices = (
                    ARIAtools.util.interp._get_height_subset_indices(
                        heightsMeta, dem_min, dem_max, pad=0))

                if len(band_indices) < len(heightsMeta):
                    band_list = [int(i + 1) for i in band_indices]
                    selected_heights = heightsMeta[band_indices]
                    subset_height_meta = '{' + ','.join(
                        str(float(value)) for value in selected_heights) + '}'
                    translate_opts = osgeo.gdal.TranslateOptions(
                        format='VRT', bandList=band_list)

                    for src_path, label in [
                            (ref_vrt_path, 'ref'), (sec_vrt_path, 'sec')]:
                        sub_vrt = outname + f'_{label}_hsubset.vrt'
                        ds_sub = osgeo.gdal.Translate(
                            sub_vrt, src_path, options=translate_opts)
                        if ds_sub is not None:
                            ds_sub.SetMetadataItem(
                                hgt_field, subset_height_meta)
                            ds_sub.FlushCache()
                        ds_sub = None
                        subset_vrts.append(sub_vrt)

                    ref_vrt_path = subset_vrts[0]
                    sec_vrt_path = subset_vrts[1]

                    LOGGER.debug(
                        'Subsetting %s inputs from %d to %d height '
                        'levels (DEM range: %.1f to %.1f)',
                        key, len(heightsMeta), len(band_indices),
                        dem_min, dem_max)
        except Exception:
            pass

    if not os.path.exists(sec_vrt_path) or not os.path.exists(ref_vrt_path):
        LOGGER.warning(f"Missing input VRTs for generate_diff: {ref_vrt_path} or {sec_vrt_path}")
        return

    ds_frame_marker = osgeo.gdal.Open(
        ref_vrt_path, osgeo.gdal.GA_ReadOnly)
    source_frame_count = (
        ds_frame_marker.GetMetadataItem('ARIA_FRAME_COUNT')
        if ds_frame_marker is not None else None)
    ds_frame_marker = None

    with rioxarray.open_rasterio(sec_vrt_path, masked=True) as da_sec:
        sec_attrs = da_sec.attrs
        sec_crs = da_sec.rio.crs
        sec_nodata = da_sec.rio.nodata
        
        with rioxarray.open_rasterio(ref_vrt_path, masked=True) as da_ref:
            arr_ref = da_ref.data
            arr_sec = da_sec.data

            if tropo_total:
                arr_total = arr_sec + arr_ref
            else:
                arr_total = arr_sec - arr_ref

            if sign_multiplier == -1:
                arr_total = arr_total * -1

            da_total = da_sec.copy()
            da_total.data = arr_total

            da_total.name = key
            og_da_attrs = sec_attrs
            da_attrs = {}
            for k in og_da_attrs:
                new_k = k.replace(OG_key, key)
                new_v = og_da_attrs[k]
                if isinstance(new_v, str):
                    new_v = new_v.replace(OG_key, key)
                da_attrs[new_k] = new_v
            if subset_height_meta is not None:
                da_attrs[hgt_field] = subset_height_meta
            da_total = da_total.assign_attrs(da_attrs)
            
            if sec_crs:
                da_total.rio.write_crs(sec_crs, inplace=True)
            if sec_nodata is not None:
                da_total.rio.write_nodata(sec_nodata, inplace=True)

            if "_FillValue" in da_total.attrs:
                del da_total.attrs["_FillValue"]

            # This is still a multi-height intermediate. Keep it GeoTIFF-
            # backed so ENVI/ISCE scanline writers are not used until after
            # the cube has been interpolated down to the final 2-D raster.
            # Preserve the source grid CRS rather than assigning ``proj``.
            diff_bands = int(da_total.sizes.get('band', 1))
            diff_width = da_total.rio.width
            diff_height = da_total.rio.height
            for diff_suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                diff_path = outname + diff_suffix
                if os.path.isfile(diff_path):
                    os.remove(diff_path)
            da_total.rio.to_raster(
                outname, driver=METADATA_STAGING_DRIVER, crs=sec_crs,
                tiled=True, compress='DEFLATE', predictor=3,
                BIGTIFF='IF_SAFER')

    validate_metadata_staging_raster(
        outname, key, expected_bands=diff_bands,
        expected_width=diff_width, expected_height=diff_height)

    if os.path.exists(outname):
        ds_src = osgeo.gdal.Open(outname, osgeo.gdal.GA_ReadOnly)
        if ds_src is not None:
            intermediate_proj = ds_src.GetProjection()
            buildvrt_options = osgeo.gdal.BuildVRTOptions(
                outputSRS=intermediate_proj or proj)
            ds_vrt = osgeo.gdal.BuildVRT(
                f'{outname}.vrt', [ds_src], options=buildvrt_options
            )

            if hgt_field in da_attrs:
                if not isinstance(da_attrs[hgt_field], (list, tuple)):
                     if isinstance(da_attrs[hgt_field], np.ndarray):
                         da_attrs[hgt_field] = da_attrs[hgt_field].tolist()

            if ds_vrt is not None:
                ds_vrt.SetMetadata(da_attrs)
                if source_frame_count is not None:
                    ds_vrt.SetMetadataItem(
                        'ARIA_FRAME_COUNT', source_frame_count)
                ds_vrt.FlushCache()
            ds_vrt = None
            ds_src = None

    for v in subset_vrts:
        if os.path.exists(v):
            os.remove(v)

    return


def extract_bperp_dict(products, num_threads):
    """Extracts bPerpendicular mean over frames for each product in products"""

    def read_and_average_bperp(frame):
        """Helper function to read and average baseline"""
        ds = osgeo.gdal.Open(frame, osgeo.gdal.GA_ReadOnly)
        arr = ds.ReadAsArray().astype(float)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        ds = None 
        
        if nodata is not None and not np.isnan(nodata):
            arr = np.where(arr == nodata, np.nan, arr)
        
        res = np.nanmean(arr)
        return res

    bperp_dict = {}
    for product in products:
        mean_bperp_by_frames = []
        for frame in product['bPerpendicular']:
            mean_bperp_by_frames.append(read_and_average_bperp(frame))

        LOGGER.debug('Pair name: %s, bPerpendicular %s' % (
            product['pair_name'][0], mean_bperp_by_frames))
        bperp_dict[product['pair_name'][0]] = float(
            np.mean(mean_bperp_by_frames))
            
    return bperp_dict


def track_existing_outputs(workdir, layers, valid_layers, ignore_names=[]):
    """Track existing layer outputs to assist dedup"""
    extracted_lyrnames = []
    for d in os.listdir(workdir):
        subdir_path = os.path.join(workdir, d)
        if os.path.isdir(subdir_path) and d not in ignore_names and \
            d not in layers and d in valid_layers:
            extracted_lyrnames.append(d)

    layers.extend(extracted_lyrnames)
    return layers


def track_correction_outputs(all_workdirs):
    """Track correction layer outputs to assist dedup"""
    existing_outputs = []
    for i in all_workdirs:
        if os.path.exists(i):
            existing_outputs.extend(
                glob.glob(os.path.join(i, '*/*[0-9].vrt')))
            existing_outputs.extend(
                glob.glob(os.path.join(i, '*/dates/*[0-9].vrt')))
            existing_outputs.extend(
                glob.glob(os.path.join(i, '*[0-9].vrt')))
            existing_outputs.extend(
                glob.glob(os.path.join(i, 'dates/*[0-9].vrt')))

    existing_outputs = list(set(existing_outputs))
    return existing_outputs


def ensure_requested_output_driver(
        path, requested_driver, expected_frame_count=None):
    """Remove cached rasters with a stale driver or frame coverage."""
    if not (os.path.exists(path) and os.path.exists(path + '.vrt')):
        return False

    ds_existing = osgeo.gdal.Open(path, osgeo.gdal.GA_ReadOnly)
    actual_driver = None
    if ds_existing is not None and ds_existing.GetDriver() is not None:
        actual_driver = ds_existing.GetDriver().ShortName
    ds_existing = None

    recorded_frame_count = None
    if expected_frame_count is not None:
        ds_existing_vrt = osgeo.gdal.Open(
            path + '.vrt', osgeo.gdal.GA_ReadOnly)
        if ds_existing_vrt is not None:
            recorded_frame_count = ds_existing_vrt.GetMetadataItem(
                'ARIA_FRAME_COUNT')
        ds_existing_vrt = None

    frame_count_matches = (
        expected_frame_count is None
        or recorded_frame_count == str(expected_frame_count))

    if (actual_driver is not None
            and actual_driver.upper() == requested_driver.upper()
            and frame_count_matches):
        return True

    LOGGER.info(
        'Regenerating %s because cached driver/frame count '
        '(%s, %s) does not match requested values (%s, %s)',
        path, actual_driver, recorded_frame_count,
        requested_driver, expected_frame_count)
    for suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        stale_path = path + suffix
        if os.path.isfile(stale_path):
            os.remove(stale_path)
    return False


def handle_epoch_layers(
        layers, product_dict, update_mode, gdal_warp_kwargs, proj, lyr_path,
        user_lyrs, map_lyrs, key, sec_key, ref_key, tropo_total, workdir,
        bounds, arrres, dem_bounds, prods_TOTbbox, dem, lat, lon, mask,
        outputFormat, verbose, multilooking, rankedResampling, num_threads,
        is_nisar_file):
    """
    Manage reference/secondary components for correction layers.
    Specifically record reference/secondary components within a `dates` subdir
    and deposit the differential fields in the level above.
    """
    LOGGER.debug('handle_epoch_layers %s' % key)
    if key == 'troposphereTotal':
        layers.append(key)
        sec_workdir = os.path.join(os.path.dirname(workdir),
                                   sec_key)
        ref_workdir = os.path.join(os.path.dirname(workdir),
                                   ref_key)

    else:
        sec_workdir = copy.deepcopy(workdir)
        ref_workdir = copy.deepcopy(workdir)

    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    all_workdirs = [workdir, sec_workdir, ref_workdir]
    all_workdirs = list(set(all_workdirs))
    existing_outputs = track_correction_outputs(all_workdirs)

    if update_mode == 'crop_only' and existing_outputs != []:
        for outname in existing_outputs:
            ifg_tag = os.path.basename(outname).split('.vrt')[0]
            crop_only_manager(outname[:-4], key, ifg_tag, gdal_warp_kwargs)

        return existing_outputs

    for i in all_workdirs:
        if not os.path.exists(i):
            os.mkdir(i)

    # NISAR SET currently agrees with the independent MintPy SET solution,
    # so preserve its existing sign behavior.
    #
    # Do NOT flip NISAR wet/hydrostatic tropo here for this diagnostic;
    # prep_nisar.py ingests the native NISAR tropo phase without this flip.
    if is_nisar_file and key == 'solidEarthTide':
        sign_multiplier = -1
    else:
        sign_multiplier = 1

    if (dem is not None
            and not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET')):
        try:
            sample_vrt = None
            for d in all_workdirs:
                vrts = glob.glob(os.path.join(d, '*.vrt'))
                if vrts:
                    sample_vrt = vrts[0]
                    break
            if sample_vrt is None:
                ds_tmp = osgeo.gdal.Open(product_dict[0][0][0])
                zdim = ds_tmp.GetMetadataItem('NETCDF_DIM_EXTRA')[1:-1]
                hgt_field_tmp = f'NETCDF_DIM_{zdim}_VALUES'
                ds_tmp = None
                hgt_str = ARIAtools.util.vrt.get_hgt_meta(
                    product_dict[0][0][0], hgt_field_tmp)
            else:
                hgt_field_tmp = [k for k in
                    osgeo.gdal.Open(sample_vrt).GetMetadata()
                    if 'DIM_' in k and 'VALUES' in k]
                hgt_field_tmp = hgt_field_tmp[0] if hgt_field_tmp else None
                hgt_str = (ARIAtools.util.vrt.get_hgt_meta(
                    sample_vrt, hgt_field_tmp) if hgt_field_tmp else None)
            if hgt_str:
                hts = np.array(hgt_str[1:-1].split(','), dtype='float32')
                dem_min, dem_max = (
                    ARIAtools.util.interp._compute_dem_range(dem))
                idx = ARIAtools.util.interp._get_height_subset_indices(
                    hts, dem_min, dem_max, pad=0)
                if len(idx) < len(hts):
                    LOGGER.info(
                        'Height subsetting %s: %d → %d levels '
                        '(DEM range: %.0f to %.0f m)',
                        key, len(hts), len(idx), dem_min, dem_max)
                else:
                    LOGGER.info(
                        'Using all %d height levels for %s '
                        '(DEM range: %.0f to %.0f m)',
                        len(hts), key, dem_min, dem_max)
        except Exception:
            pass

    all_outputs = []
    prog_bar = ARIAtools.util.misc.ProgressBar(
        maxValue=len(product_dict[0]), prefix=f'Exporting {key}: '
    )
    for i in enumerate(product_dict[0]):
        ifg = product_dict[1][i[0]][0]
        outname = os.path.abspath(os.path.join(workdir, ifg))

        model_name = None
        if 'tropo' in key:
            out_dir = os.path.dirname(outname)
            if not is_nisar_file:
                model_name = i[1][0].split('/')[-3]
                out_dir = os.path.join(out_dir, model_name)
                outname = os.path.join(out_dir, ifg)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

        expected_frame_count = (
            len(i[1]) if is_nisar_file and isinstance(i[1], list) else None)
        if ensure_requested_output_driver(
                outname, outputFormat, expected_frame_count):
            continue

        if ref_key in user_lyrs or tropo_total:
            ref_outname = os.path.abspath(os.path.join(ref_workdir, ifg))
            hgt_field, ref_outname = prep_metadatalayers(
                ref_outname, i[1], dem, ref_key, layers, is_nisar_file, proj,
                outputFormat, model_name, sign_multiplier)

        if model_name is not None:
            all_outputs.append(os.path.join(workdir, model_name))
            all_outputs.append(os.path.join(ref_workdir, model_name))
            all_outputs.append(os.path.join(sec_workdir, model_name))
        else:
            all_outputs.append(workdir)
            all_outputs.append(ref_workdir)
            all_outputs.append(sec_workdir)

        if 'tropo' in key:
            sec_outname = os.path.abspath(os.path.join(sec_workdir, ifg))

            if not is_nisar_file:
                wet_path = os.path.join(
                    lyr_path, model_name, 'reference', ref_key)
                dry_path = os.path.join(
                    lyr_path, model_name, 'reference', sec_key)
            else:
                wet_path = os.path.join(lyr_path, map_lyrs[0])
                dry_path = os.path.join(lyr_path, map_lyrs[1])

            sec_comp = [
                j.replace(wet_path, dry_path) for j in i[1]]

            if sec_key in user_lyrs or tropo_total:
                hgt_field, sec_outname = prep_metadatalayers(
                    sec_outname, sec_comp, dem, sec_key, layers, is_nisar_file,
                    proj, outputFormat, model_name, sign_multiplier)

            if tropo_total:
                model_dir = os.path.abspath(workdir)

                if not is_nisar_file:
                    model_dir = os.path.join(model_dir, model_name)
                    ref_diff = ref_outname
                    sec_diff = sec_outname
                    outname_diff = os.path.join(model_dir, 'dates',
                                                os.path.basename(ref_diff))
                    if not os.path.exists(outname_diff):
                        generate_diff(
                            ref_diff, sec_diff, outname_diff, key, sec_key,
                            tropo_total, hgt_field, proj, outputFormat,
                            dem=dem)
                    ref_diff = os.path.join(os.path.dirname(ref_outname),
                                            ifg.split('_')[1])
                    sec_diff = os.path.join(os.path.dirname(sec_outname),
                                            ifg.split('_')[1])
                    outname_diff = os.path.join(model_dir, 'dates',
                                                os.path.basename(ref_diff))
                    if not os.path.exists(outname_diff):
                        generate_diff(
                            ref_diff, sec_diff, outname_diff, key, sec_key,
                            tropo_total, hgt_field, proj, outputFormat,
                            dem=dem)

                    ref_diff = os.path.join(model_dir, 'dates', ifg.split('_')[0])
                    sec_diff = os.path.join(model_dir, 'dates', ifg.split('_')[1])
                    outname = os.path.join(model_dir, ifg)
                    generate_diff(
                        ref_diff, sec_diff, outname, key, sec_key, False,
                        hgt_field, proj, outputFormat, dem=dem)

                else:
                    ref_diff = os.path.join(ref_workdir, ifg)
                    sec_diff = os.path.join(sec_workdir, ifg)
                    outname = os.path.join(model_dir, ifg)
                    generate_diff(
                        ref_diff, sec_diff, outname, key, sec_key, tropo_total,
                        hgt_field, proj, outputFormat, dem=dem)

        else:
            sec_outname = os.path.dirname(ref_outname)
            sec_outname = os.path.abspath(
                os.path.join(sec_outname, ifg.split('_')[0]))

        prog_bar.update(i[0] + 1)
    prog_bar.close()

    prod_ver_list = i[1]
    for i in all_workdirs:
        key_name = os.path.basename(i)
        if os.path.exists(i):
            if key_name not in layers or len(os.listdir(i)) == 0:
                shutil.rmtree(i)

    all_outputs = list(set(all_outputs))
    for i in enumerate(all_outputs):
        if os.path.exists(i[1]):
            record_epochs = []
            record_epochs.extend(
                glob.glob(os.path.join(i[1], '*[0-9].vrt')))
            record_epochs.extend(
                glob.glob(os.path.join(i[1], 'dates/*[0-9].vrt')))

            for j in enumerate(record_epochs):
                if not os.path.exists(j[1]):
                    continue
                ds_count = osgeo.gdal.Open(j[1], osgeo.gdal.GA_ReadOnly)
                if ds_count is None:
                    continue
                band_count = ds_count.RasterCount
                ds_count = None
                if band_count == 1:
                    if j[0] == 0 and os.path.exists(j[1][:-4] + '.vrt'):
                        ref_wid, ref_hgt, ref_geotrans, _, _ = \
                            ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                        ref_arr = [ref_wid, ref_hgt, ref_geotrans, j[1][:-4]]

                    continue

                finalize_metadata(
                    j[1][:-4], bounds, arrres, dem_bounds, prods_TOTbbox,
                    dem, lat, lon, hgt_field, prod_ver_list, is_nisar_file,
                    outputFormat, verbose=verbose)

                if mask is not None and os.path.exists(j[1][:-4] + '.vrt'):
                    ds_vrt_read = osgeo.gdal.Open(
                        j[1][:-4] + '.vrt', osgeo.gdal.GA_ReadOnly
                    )
                    if ds_vrt_read is not None:
                        vrt_arr = ds_vrt_read.ReadAsArray()
                        ds_vrt_read = None
                        mask_arr = mask.ReadAsArray() * vrt_arr

                        update_file = osgeo.gdal.Open(
                            j[1][:-4], osgeo.gdal.GA_Update
                        )
                        if update_file is not None:
                            update_file.GetRasterBand(1).WriteArray(mask_arr)
                            update_file.FlushCache()
                        update_file = None
                        mask_arr = None

                if j[0] == 0 and os.path.exists(j[1][:-4] + '.vrt'):
                    ref_wid, ref_hgt, ref_geotrans, _, _ = \
                        ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                    ref_arr = [ref_wid, ref_hgt, ref_geotrans, j[1][:-4]]

                elif os.path.exists(j[1][:-4] + '.vrt'):
                    prod_wid, prod_hgt, prod_geotrans, _, _ = \
                        ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                    prod_arr = [prod_wid, prod_hgt, prod_geotrans, j[1][:-4]]
                    ARIAtools.util.vrt.dim_check(ref_arr, prod_arr)
                prev_outname = j[1][:-4]

    existing_outputs = track_correction_outputs(all_workdirs)

    return existing_outputs


def export_product_worker_helper(args):
    """Calls export_product_worker with * expanded args"""
    return export_product_worker(*args)


def export_product_worker(
        ii, ilayer, product, proj, full_product_dict_file, layers, workdir,
        bounds, prods_TOTbbox, demfile, demfile_expanded, maskfile,
        outputFormat, outputFormatPhys, layer, outDir,
        arrres, num_threads, multilooking, verbose, is_nisar_file,
        range_correction, rankedResampling, update_mode):
    """
    Worker function for export_products for parallel execution with
    multiprocessing package.
    """
    ARIAtools.product._configure_gdal_virtual_access()

    gdal_warp_kwargs = {
        'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
        'outputBounds': bounds, 'xRes': arrres[0], 'yRes': arrres[1],
        'targetAlignedPixels': True, 'multithread': False, 'dstSRS': proj}
    warp_options = osgeo.gdal.WarpOptions(
        **gdal_warp_kwargs
    )

    sign_multiplier = -1 if (is_nisar_file and layer in [
        'bPerpendicular', 'bParallel']) else 1

    mask = None if maskfile is None else osgeo.gdal.Open(maskfile)
    dem = None if demfile is None else osgeo.gdal.Open(demfile)
    dem_expanded = (
        None
        if demfile_expanded is None else osgeo.gdal.Open(demfile_expanded))

    if dem_expanded is not None:
        gt = dem_expanded.GetGeoTransform()
        xs, ys = dem_expanded.RasterXSize, dem_expanded.RasterYSize

        lat = np.linspace(gt[3], gt[3] + (gt[5] * (ys - 1)), ys)
        lat = np.repeat(lat[:, np.newaxis], xs, axis=1)
        lon = np.linspace(gt[0], gt[0] + (gt[1] * (xs - 1)), xs)
        lon = np.repeat(lon[:, np.newaxis], ys, axis=1).T

        dem_bounds = [
            gt[0], gt[3] + (gt[-1] * dem_expanded.RasterYSize),
            gt[0] + (gt[1] * dem_expanded.RasterXSize), gt[3]]

    with open(full_product_dict_file, 'r') as ifp:
        full_product_dict = json.load(ifp)

    product_dict = [[j[layers[ilayer]] for j in full_product_dict],
                    [j["pair_name"] for j in full_product_dict]]

    ifg_tag = product_dict[1][ii][0]
    outname = os.path.abspath(os.path.join(workdir, ifg_tag))

    is_metadata_product = (
        any(':/science/grids/imagingGeometry' in s for s in product)
        or any(':/science/LSAR/GUNW/metadata/radarGrid' in s
               for s in product))
    expected_frame_count = (
        len(product) if is_nisar_file and is_metadata_product else None)
    existing_output_matches = ensure_requested_output_driver(
        outname, outputFormatPhys, expected_frame_count)

    if update_mode != 'crop_only' and existing_output_matches:
        LOGGER.debug('Skipping %s - %s', ifg_tag,
                     {os.path.dirname(outname).split('/')[-1]})

    elif update_mode == 'crop_only' \
            and os.path.exists(outname) \
            and os.path.exists(outname + '.vrt'):
        lyrname = os.path.dirname(outname).split('/')[-1]
        crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs)
        if os.path.dirname(outname).split('/')[-1] == 'unwrappedPhase':
            lyrname = 'connectedComponents'
            path_parts = outname.split('/')

            if path_parts[-2] == 'unwrappedPhase':
                path_parts[-2] = lyrname

            outname = '/'.join(path_parts)
            crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs)

    else:
        LOGGER.debug('Extracting %s - %s', ifg_tag,
                     {os.path.dirname(outname).split('/')[-1]})

        if is_metadata_product:
            hgt_field, outname = prep_metadatalayers(
                outname, product, dem_expanded, layer, layers,
                is_nisar_file, proj, outputFormatPhys,
                sign_multiplier=sign_multiplier)

            finalize_metadata(
                outname, bounds, arrres, dem_bounds, prods_TOTbbox,
                dem_expanded, lat, lon, hgt_field, product, is_nisar_file,
                outputFormatPhys, verbose=verbose)

        elif layer != 'unwrappedPhase' and layer != 'connectedComponents':

            if is_nisar_file:

                if layer == 'amplitude':
                    amp_ds_list = []
                    prod_list = product if isinstance(product, list) else [product]
                    
                    mem_driver = osgeo.gdal.GetDriverByName('MEM')
                    for prod_frame in prod_list:
                        ds_in = osgeo.gdal.Open(prod_frame)
                        complex_data = ds_in.GetRasterBand(1).ReadAsArray()
                        
                        ds_amp = mem_driver.Create(
                            '', ds_in.RasterXSize, ds_in.RasterYSize, 1, osgeo.gdal.GDT_Float32
                        )
                        ds_amp.SetProjection(ds_in.GetProjection())
                        ds_amp.SetGeoTransform(ds_in.GetGeoTransform())
                        
                        amp_band = ds_amp.GetRasterBand(1)
                        amp_arr = np.abs(complex_data)
                        
                        amp_arr[np.isnan(amp_arr)] = 0
                        
                        amp_band.WriteArray(amp_arr)
                        amp_band.SetNoDataValue(0)
                        
                        amp_ds_list.append(ds_amp)
                        ds_in = None

                    amp_kwargs = gdal_warp_kwargs.copy()
                    amp_kwargs['format'] = 'GTiff'
                    amp_kwargs['multithread'] = False
                    if ('dstSRS' in amp_kwargs and
                            isinstance(amp_kwargs['dstSRS'], int)):
                        amp_kwargs['dstSRS'] = f"EPSG:{amp_kwargs['dstSRS']}"

                    amp_warp_opts = osgeo.gdal.WarpOptions(
                        outputType=osgeo.gdal.GDT_Float32,
                        srcNodata=0,
                        dstNodata=np.nan,
                        **amp_kwargs
                    )

                    tmp_amp = outname + '.tmp.tif'
                    if os.path.exists(tmp_amp):
                        os.remove(tmp_amp)

                    ds_amp_warp = osgeo.gdal.Warp(
                        tmp_amp, amp_ds_list, options=amp_warp_opts
                    )
                    ds_amp_warp = None
                    amp_ds_list = None

                    for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                        f_to_rm = f"{outname}{ext}"
                        if os.path.exists(f_to_rm):
                            os.remove(f_to_rm)

                    osgeo.gdal.Translate(outname, tmp_amp, format=outputFormatPhys)
                    if os.path.exists(tmp_amp):
                        os.remove(tmp_amp)

                else:
                    if isinstance(product, list) and len(product) > 1:
                        tmp_mosaic = str(outname) + "_uncropped.vrt"
                        
                        tmp_vrts = []
                        for idx, p in enumerate(product):
                            t_vrt = f"{outname}_{idx}_tmp.vrt"
                            ds_tmp = osgeo.gdal.Warp(
                                t_vrt, p, format="VRT", dstSRS=proj
                            )
                            if ds_tmp is not None:
                                ds_tmp.FlushCache()
                            ds_tmp = None
                            tmp_vrts.append(t_vrt)
                            
                        ds_mosaic = osgeo.gdal.BuildVRT(tmp_mosaic, tmp_vrts)
                        if ds_mosaic is not None:
                            ds_mosaic.FlushCache()
                        ds_mosaic = None
                        warp_inputs = tmp_mosaic
                    else:
                        warp_inputs = (
                            product[0] if isinstance(product, list)
                            else product
                        )

                    if outputFormat == 'VRT':
                        ds = osgeo.gdal.Warp(
                            outname + '.vrt', warp_inputs, options=warp_options
                        )
                        if ds is not None:
                            ds.FlushCache()
                        ds = None
                    else:
                        tmp_tif = outname + '.tmp.tif'
                        if os.path.exists(tmp_tif):
                            os.remove(tmp_tif)

                        warp_kwargs_gtiff = gdal_warp_kwargs.copy()
                        warp_kwargs_gtiff['format'] = 'GTiff'
                        warp_kwargs_gtiff['multithread'] = False
                        warp_opts_gtiff = osgeo.gdal.WarpOptions(**warp_kwargs_gtiff)

                        ds = osgeo.gdal.Warp(
                            tmp_tif, warp_inputs, options=warp_opts_gtiff
                        )
                        ds = None

                        for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                            f_to_rm = f"{outname}{ext}"
                            if os.path.exists(f_to_rm):
                                os.remove(f_to_rm)

                        osgeo.gdal.Translate(outname, tmp_tif, format=outputFormatPhys)
                        if os.path.exists(tmp_tif):
                            os.remove(tmp_tif)

            else:
                with osgeo.gdal.config_options(
                    {"GDAL_NUM_THREADS": num_threads}
                ):

                    if outputFormat == 'VRT':
                        ds_vrt = osgeo.gdal.BuildVRT(
                            outname + "_uncropped.vrt", product
                        )
                        if ds_vrt is not None:
                            ds_vrt.FlushCache()
                        ds_vrt = None
                        ds = osgeo.gdal.Warp(
                            outname + '.vrt',
                            outname + '_uncropped.vrt',
                            options=warp_options
                        )
                        if ds is not None:
                            ds.FlushCache()
                        ds = None
                    else:
                        ds_vrt = osgeo.gdal.BuildVRT(outname + '.vrt', product)
                        if ds_vrt is not None:
                            ds_vrt.FlushCache()
                        ds_vrt = None

                        tmp_tif = outname + '.tmp.tif'
                        if os.path.exists(tmp_tif):
                            os.remove(tmp_tif)

                        warp_kwargs_gtiff = gdal_warp_kwargs.copy()
                        warp_kwargs_gtiff['format'] = 'GTiff'
                        warp_kwargs_gtiff['multithread'] = False
                        warp_opts_gtiff = osgeo.gdal.WarpOptions(**warp_kwargs_gtiff)

                        ds = osgeo.gdal.Warp(
                            tmp_tif,
                            outname + '.vrt',
                            options=warp_opts_gtiff
                        )
                        ds = None

                        for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
                            f_to_rm = f"{outname}{ext}"
                            if os.path.exists(f_to_rm):
                                os.remove(f_to_rm)

                        osgeo.gdal.Translate(outname, tmp_tif, format=outputFormatPhys)
                        if os.path.exists(tmp_tif):
                            os.remove(tmp_tif)

                        ds_trans = osgeo.gdal.Translate(
                            outname + '.vrt',
                            outname,
                            options=osgeo.gdal.TranslateOptions(
                                format="VRT"
                            )
                        )
                        if ds_trans is not None:
                            ds_trans.FlushCache()
                        ds_trans = None

            if os.path.exists(outname):
                ds_src = osgeo.gdal.Open(outname, osgeo.gdal.GA_ReadOnly)
                if ds_src is not None:
                    ds_trans = osgeo.gdal.Translate(
                        outname + '.vrt', ds_src, format="VRT"
                    )
                    if ds_trans is not None:
                        ds_trans.FlushCache()
                    ds_trans = None
                    ds_src = None

            if outputFormat != 'VRT':
                tmp_mosaic = str(outname) + "_uncropped.vrt"
                if os.path.exists(tmp_mosaic):
                    os.remove(tmp_mosaic)
                
                if isinstance(product, list) and len(product) > 1:
                    for idx in range(len(product)):
                        t_vrt = f"{outname}_{idx}_tmp.vrt"
                        if os.path.exists(t_vrt):
                            os.remove(t_vrt)

        else:
            conn_files = full_product_dict[ii]['connectedComponents']
            prod_bbox_files = full_product_dict[ii][
                'productBoundingBoxFrames']
            outFileConnComp = os.path.join(
                outDir, 'connectedComponents', ifg_tag)

            outFilePhs = os.path.join(outDir, 'unwrappedPhase', ifg_tag)

            phs_files = full_product_dict[ii]['unwrappedPhase']

            ARIAtools.util.seq_stitch.product_stitch_sequential(
                phs_files, conn_files, arrres=arrres, epsg=proj,
                bounds=bounds, clip_json=prods_TOTbbox, output_unw=outFilePhs,
                output_conn=outFileConnComp,
                output_format=outputFormatPhys,
                is_nisar_file=is_nisar_file,
                range_correction=range_correction, save_fig=False,
                overwrite=True)

            if multilooking is not None:
                ARIAtools.util.vrt.resampleRaster(
                    outFilePhs, multilooking, bounds, prods_TOTbbox,
                    rankedResampling, outputFormat=outputFormatPhys,
                    num_threads=num_threads)

            if mask is not None:
                for j in [outFileConnComp, outFilePhs]:
                    ds_vrt_read = osgeo.gdal.Open(
                        j + '.vrt', osgeo.gdal.GA_ReadOnly
                    )
                    vrt_arr = ds_vrt_read.ReadAsArray()
                    ds_vrt_read = None
                    mask_arr = mask.ReadAsArray() * vrt_arr

                    update_file = osgeo.gdal.Open(
                        j, osgeo.gdal.GA_Update
                    )
                    if update_file is not None:
                        update_file.GetRasterBand(1).WriteArray(mask_arr)
                        update_file.FlushCache()
                    update_file = None
                    mask_arr = None

        if layer != 'unwrappedPhase' and layer != 'connectedComponents':

            if multilooking is not None:
                ARIAtools.util.vrt.resampleRaster(
                    outname, multilooking, bounds, prods_TOTbbox,
                    rankedResampling, outputFormat=outputFormatPhys,
                    num_threads=num_threads)

            if mask is not None:
                ds_vrt_read = osgeo.gdal.Open(
                    outname + '.vrt', osgeo.gdal.GA_ReadOnly
                )
                vrt_arr = ds_vrt_read.ReadAsArray()
                ds_vrt_read = None
                mask_arr = mask.ReadAsArray() * vrt_arr

                update_file = osgeo.gdal.Open(
                    outname, osgeo.gdal.GA_Update
                )
                if update_file is not None:
                    update_file.GetRasterBand(1).WriteArray(mask_arr)
                    update_file.FlushCache()
                update_file = None
                mask_arr = None

    prod_wid, prod_hgt, prod_geotrans, _, _ = \
        ARIAtools.util.vrt.get_basic_attrs(outname + '.vrt')
    prev_outname = os.path.abspath(os.path.join(workdir, ifg_tag))
    prod_arr = [
        prod_wid, prod_hgt, prod_geotrans, os.path.join(workdir, ifg_tag)]

    return ii, ilayer, prev_outname, prod_arr


def export_products(
        full_product_dict, proj, bbox_file, prods_TOTbbox, layers, arrres,
        iono_filter, is_nisar_file, rankedResampling=False, demfile=None,
        demfile_expanded=None, lat=None, lon=None, maskfile=None, outDir='./',
        outputFormat='VRT', verbose=None, num_threads='2', multilooking=None,
        tropo_total=False, model_names=[], multiproc_method='single',
        runlog=None):
    """
    Export layer and 2D meta-data layers (at the product resolution).
    The function finalize_metadata is called to derive the 2D metadata layer.
    Dem/lat/lon arrays must be passed for this process.
    The keys specify which layer to extract from the dictionary.
    All products are cropped by the bounds from the input bbox_file,
    and clipped to the track extent denoted by the input prods_TOTbbox.
    Optionally, a user may pass a mask-file.
    """

    start_time = time.time()
    LOGGER.debug('export_products, layers: {}'.format(layers))

    if not layers and not tropo_total:
        return

    ref_wid = None
    ref_hgt = None
    ref_geotrans = None
    ref_arr = None

    mask = None if maskfile is None else osgeo.gdal.Open(maskfile)
    dem = None if demfile is None else osgeo.gdal.Open(demfile)
    dem_expanded = (
        None if demfile_expanded is None
        else osgeo.gdal.Open(demfile_expanded))

    srs = osgeo.osr.SpatialReference()
    srs.ImportFromWkt(proj)
    srs.AutoIdentifyEPSG()
    epsg_code = f'EPSG:{int(srs.GetAuthorityCode(None))}'
    srs = None
    lyr_input_dict = {
        'layers': layers, 'prods_TOTbbox': prods_TOTbbox,
        'proj': epsg_code, 'dem': dem_expanded, 'lat': lat, 'lon': lon,
        'mask': mask, 'verbose': verbose, 'multilooking': multilooking,
        'rankedResampling': rankedResampling, 'num_threads': num_threads,
        'is_nisar_file': is_nisar_file}

    range_correction = True
    track_fileext = full_product_dict[0]['unwrappedPhase'][0]
    if is_nisar_file:
        range_correction = False
        model_names = ['']
    else:
        model_names = [f'_{i}' for i in model_names]
    lyr_input_dict['is_nisar_file'] = is_nisar_file

    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds
    lyr_input_dict['bounds'] = bounds
    lyr_input_dict['arrres'] = arrres
    if dem_expanded is not None:
        dem_gt = dem_expanded.GetGeoTransform()
        dem_bounds = [
            dem_gt[0], dem_gt[3] + (dem_gt[-1] * dem_expanded.RasterYSize),
            dem_gt[0] + (dem_gt[1] * dem_expanded.RasterXSize), dem_gt[3]]
        lyr_input_dict['dem_bounds'] = dem_bounds

    if (outputFormat == 'VRT' and mask is not None) or \
            (outputFormat == 'VRT' and multilooking is not None):
        outputFormat = 'ENVI'

    outputFormatPhys = 'ENVI'
    if outputFormat != 'VRT':
        outputFormatPhys = outputFormat
    lyr_input_dict['outputFormat'] = outputFormatPhys

    if runlog is None:
        update_mode = 'full_extract'
    else:
        log_data = runlog.load()
        update_mode = log_data['update_mode']
        if 'update_mode' in log_data.keys():
            update_mode = log_data['update_mode']

        prev_maskfile = log_data['maskfilename'] if 'maskfilename' \
            in log_data.keys() else None

        if maskfile != prev_maskfile and prev_maskfile is not None:
            update_mode = 'full_extract'
            LOGGER.warning(
                'Mask file has changed. Setting update mode to full_extract.')

        runlog.update('maskfilename', maskfile)

        prev_demfile = log_data['demfile'] if 'demfile' \
            in log_data.keys() else None

        if demfile != prev_demfile and prev_demfile is not None:
            update_mode = 'full_extract'
            LOGGER.warning(
                'DEM file has changed. Setting update mode to full_extract.')

        runlog.update('demfile', demfile)
        runlog.update('update_mode', update_mode)

    extracted_files = []

    gdal_warp_kwargs = {
        'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
        'outputBounds': bounds, 'xRes': arrres[0], 'yRes': arrres[1],
        'targetAlignedPixels': True, 'multithread': False, 'dstSRS': epsg_code}

    lyr_input_dict['update_mode'] = update_mode
    lyr_input_dict['gdal_warp_kwargs'] = gdal_warp_kwargs

    tropo_lyrs = ['troposphereWet', 'troposphereHydrostatic']
    user_lyrs = list(set(layers).intersection(tropo_lyrs))
    if tropo_total or user_lyrs != []:
        if is_nisar_file:
            lyr_prefix = '/science/LSAR/GUNW/metadata/radarGrid/'
        else:
            lyr_prefix = '/science/grids/corrections/external/troposphere/'
        key = 'troposphereTotal'
        wet_key = 'troposphereWet'
        dry_key = 'troposphereHydrostatic'
        workdir = os.path.join(outDir, key)
        lyr_input_dict['lyr_path'] = lyr_prefix
        lyr_input_dict['user_lyrs'] = user_lyrs
        lyr_input_dict['key'] = key
        lyr_input_dict['sec_key'] = dry_key
        lyr_input_dict['ref_key'] = wet_key
        lyr_input_dict['tropo_total'] = tropo_total
        lyr_input_dict['workdir'] = workdir

        for i in model_names:
            model = wet_key + f'{i}'
            tropo_lyrs.append(model)
            tropo_lyrs.append(dry_key + f'{i}')
            product_dict = [
                [j[model] for j in full_product_dict if model in j.keys()],
                [j["pair_name"]
                    for j in full_product_dict if model in j.keys()]
            ]
            product_dict_dry = [
                j[dry_key + f'{i}'] for j in full_product_dict if
                dry_key + f'{i}' in j.keys()
            ]

            map_lyrs = [
                product_dict[0][0][0].split('/')[-1],
                product_dict_dry[0][0].split('/')[-1]
            ]
            lyr_input_dict['map_lyrs'] = map_lyrs

            lyr_input_dict['product_dict'] = product_dict

            extracted_files.extend(handle_epoch_layers(**lyr_input_dict))

            tag = i.split('_')[-1]
            prev_outname = os.path.abspath(
                os.path.join(workdir,
                             tag,
                             product_dict[1][0][0])
            )
            if os.path.exists(prev_outname + '.vrt'):
                prev_outname_check = copy.deepcopy(prev_outname)

        if 'prev_outname_check' in locals():
            ref_wid, ref_hgt, ref_geotrans, _, _ = \
                ARIAtools.util.vrt.get_basic_attrs(prev_outname_check + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans, prev_outname]

    tropo_lyrs = list(set(tropo_lyrs))
    ext_corr_lyrs = tropo_lyrs + ['solidEarthTide', 'troposphereTotal']
    if 'solidEarthTide' in layers:
        if is_nisar_file:
            lyr_prefix = '/science/LSAR/GUNW/metadata/radarGrid/'
        else:
            lyr_prefix = '/science/grids/corrections/external/tides/solidEarth/'
        key = 'solidEarthTide'
        ref_key = key
        sec_key = key
        product_dict = [
            [j[key] for j in full_product_dict if key in j.keys()],
            [j["pair_name"] for j in full_product_dict if key in j.keys()]]

        map_lyrs = [product_dict[0][0][0].split('/')[-1]]
        lyr_input_dict['map_lyrs'] = map_lyrs

        workdir = os.path.join(outDir, key)
        prev_outname = copy.deepcopy(workdir)

        lyr_input_dict['product_dict'] = product_dict
        lyr_input_dict['lyr_path'] = lyr_prefix
        lyr_input_dict['user_lyrs'] = ['solidEarthTide']
        lyr_input_dict['key'] = key
        lyr_input_dict['sec_key'] = sec_key
        lyr_input_dict['ref_key'] = ref_key
        lyr_input_dict['tropo_total'] = False
        lyr_input_dict['workdir'] = workdir

        extracted_files.extend(handle_epoch_layers(**lyr_input_dict))

        prev_outname = os.path.abspath(os.path.join(workdir,
                                       product_dict[1][0][0]))
        if os.path.exists(prev_outname + '.vrt'):
            ref_wid, ref_hgt, ref_geotrans, \
                _, _ = ARIAtools.util.vrt.get_basic_attrs(prev_outname + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans,
                       prev_outname]

    ext_corr_lyrs += ['ionosphere']
    if 'ionosphere' in layers:
        lyr_prefix = '/science/grids/corrections/derived/ionosphere/ionosphere'
        key = 'ionosphere'
        product_dict = \
            [[j[key] for j in full_product_dict if key in j.keys()],
             [j["pair_name"] for j in full_product_dict if key in j.keys()]]

        workdir = os.path.join(outDir, key)
        prev_outname = copy.deepcopy(workdir)

        if multilooking is not None:
            iono_arrres = [arrres[0] * multilooking, arrres[1] * multilooking]
        else:
            iono_arrres = arrres

        lyr_input_dict = dict(input_iono_files=None,
                              arrres=iono_arrres,
                              epsg=epsg_code,
                              output_iono=None,
                              output_format=outputFormat,
                              bounds=bounds,
                              clip_json=prods_TOTbbox,
                              mask_file=mask,
                              iono_filter=iono_filter,
                              is_nisar_file=is_nisar_file,
                              verbose=verbose,
                              overwrite=True)

        prog_bar = ARIAtools.util.misc.ProgressBar(
            maxValue=len(product_dict[0]), prefix='Exporting ionosphere: '
        )
        for i, layer in enumerate(product_dict[0]):
            outname = os.path.abspath(
                os.path.join(
                    workdir,
                    product_dict[1][i][0]))
            lyr_input_dict['input_iono_files'] = layer
            lyr_input_dict['output_iono'] = outname

            if os.path.exists(outname) and update_mode == 'crop_only':
                crop_only_manager(outname, 'ionosphere',
                    product_dict[1][i][0], gdal_warp_kwargs)

            if not os.path.exists(outname):
                ARIAtools.util.ionosphere.export_ionosphere(**lyr_input_dict)

            extracted_files.append(outname)

            if os.path.exists(outname + '.vrt'):
                prev_outname_check = copy.deepcopy(outname)
                
            prog_bar.update(i + 1)
        prog_bar.close()

        if 'prev_outname_check' in locals():
            ref_wid, ref_hgt, ref_geotrans, _, _ = \
                ARIAtools.util.vrt.get_basic_attrs(prev_outname_check + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans, prev_outname]

    if runlog is not None:
        runlog.update('extracted_files', extracted_files)

    layers = [i for i in layers if i not in ext_corr_lyrs]

    full_product_dict_file = os.path.join(outDir, 'full_product_dict.json')
    with open(full_product_dict_file, 'w') as ofp:
        json.dump(full_product_dict, ofp)

    for ilayer, layer in enumerate(layers):

        product_dict = [[j[layer] for j in full_product_dict],
                        [j["pair_name"] for j in full_product_dict]]

        workdir = os.path.join(outDir, layer)
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        if (dem_expanded is not None
                and not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET')
                and any(':/science/grids/imagingGeometry' in s
                        or ':/science/LSAR/GUNW/metadata/radarGrid' in s
                        for s in product_dict[0][0])
                and layer not in ['ionosphere']):
            try:
                ds_tmp = osgeo.gdal.Open(product_dict[0][0][0])
                zdim = ds_tmp.GetMetadataItem('NETCDF_DIM_EXTRA')
                ds_tmp = None
                if zdim:
                    hgt_field_tmp = f'NETCDF_DIM_{zdim[1:-1]}_VALUES'
                    hgt_str = ARIAtools.util.vrt.get_hgt_meta(
                        product_dict[0][0][0], hgt_field_tmp)
                    if hgt_str:
                        hts = np.array(
                            hgt_str[1:-1].split(','), dtype='float32')
                        d_min, d_max = (
                            ARIAtools.util.interp
                               ._compute_dem_range(
                                   dem_expanded))
                        idx = (ARIAtools.util.interp
                               ._get_height_subset_indices(
                                   hts, d_min, d_max, pad=0))
                        if len(idx) < len(hts):
                            LOGGER.info(
                                'Height subsetting %s: %d → %d levels '
                                '(DEM range: %.0f to %.0f m)',
                                layer, len(hts), len(idx), d_min, d_max)
                        else:
                            LOGGER.info(
                                'Using all %d height levels for %s '
                                '(DEM range: %.0f to %.0f m)',
                                len(hts), layer, d_min, d_max)
            except Exception:
                pass

        mp_args = []
        for ii, product in enumerate(product_dict[0]):
            ifg_tag = product_dict[1][ii][0]
            outname = os.path.abspath(os.path.join(workdir, ifg_tag))
            extracted_files.append(outname)
            if layer == 'unwrappedPhase':
                extracted_files.append(outname.replace(
                    'unwrappedPhase', 'connectedComponents'))

            mp_args.append((
                ii, ilayer, product, epsg_code, full_product_dict_file, layers,
                workdir, bounds, prods_TOTbbox, demfile,
                demfile_expanded, maskfile, outputFormat, outputFormatPhys,
                layer, outDir, arrres, num_threads,
                multilooking, verbose, is_nisar_file, range_correction,
                rankedResampling, update_mode))

        if int(num_threads) == 1 or multiproc_method in ['single', 'threads']:
            
            prog_bar = ARIAtools.util.misc.ProgressBar(
                maxValue=len(mp_args), prefix=f'Exporting {layer}: '
            )

            if multiproc_method == 'single':
                outputs = []
                for i, arg in enumerate(mp_args):
                    outputs.append(export_product_worker_helper(arg))
                    prog_bar.update(i + 1)
                    sys.stdout.flush()
                prog_bar.close()
                
            else:
                LOGGER.debug('Running %d total jobs with threads', len(mp_args))

                lock = threading.Lock()
                completed = 0

                def update_progress(result):
                    nonlocal completed
                    with lock:
                        completed += 1
                        prog_bar.update(completed)
                    return result

                jobs = []
                for arg in mp_args:
                    job = dask.delayed(
                        lambda x: update_progress(export_product_worker(*x))
                    )(arg)
                    jobs.append(job)

                outputs = dask.compute(
                    jobs, num_workers=int(num_threads), scheduler='threads'
                )[0]
                prog_bar.close()

            for ii_out, ilayer_out, outname, prod_arr in outputs:
                if ref_arr is None:
                    ref_arr = copy.deepcopy(prod_arr)
                else:
                    ARIAtools.util.vrt.dim_check(ref_arr, prod_arr)
                prev_outname = outname

        elif multiproc_method == 'gnu_parallel':
            export_workers_temp_dir = os.path.join(outDir, 'export_workers')
            if os.path.isdir(export_workers_temp_dir):
                shutil.rmtree(export_workers_temp_dir)
            os.mkdir(export_workers_temp_dir)

            for ii_arg, args in enumerate(mp_args):
                this_json_file = os.path.join(
                    outDir, 'export_workers',
                    'export_product_args_%2.2d.json' % ii_arg)
                with open(this_json_file, 'w') as ofp:
                    json.dump(args, ofp)

            LOGGER.debug('Running %d total jobs in parallel' % len(mp_args))
            
            prog_bar = ARIAtools.util.misc.ProgressBar(
                maxValue=len(mp_args), prefix=f'Exporting {layer}: '
            )

            proc = subprocess.Popen((
                'find %s/export_workers -name "export_product_args_*.json" | '
                'parallel -j %d export_product.py {}') % (
                    outDir, int(num_threads)), shell=True)

            while proc.poll() is None:
                num_done = len(glob.glob(
                    os.path.join(export_workers_temp_dir, 'outputs_*.json')
                ))
                prog_bar.update(num_done)
                time.sleep(1.0)

            num_done = len(glob.glob(
                os.path.join(export_workers_temp_dir, 'outputs_*.json')
            ))
            prog_bar.update(num_done)
            prog_bar.close()

            output_files = glob.glob(os.path.join(
                export_workers_temp_dir, 'outputs_*.json'))

            if len(output_files) > 0:
                if ref_arr is None:
                    with open(output_files[0]) as ifp:
                        output_dict = json.load(ifp)
                        ref_arr = copy.deepcopy(output_dict['prod_arr'])

                for output_file in output_files:
                    with open(output_file) as ifp:
                        output_dict = json.load(ifp)
                    ARIAtools.util.vrt.dim_check(ref_arr, output_dict['prod_arr'])
                    prev_outname = output_dict['outname']

    end_time = time.time()
    LOGGER.debug(
        "export_product_worker took %f seconds" % (end_time - start_time))

    if runlog is not None:
        runlog.update('extracted_files', extracted_files)

    plots_subdir = os.path.abspath(
        os.path.join(outDir, 'metadatalyr_plots'))
    if os.path.exists(plots_subdir) and len(os.listdir(plots_subdir)) == 0:
        shutil.rmtree(plots_subdir)

    retval = [None]*4 if ref_arr is None else ref_arr

    return retval


def finalize_metadata(outname, bbox_bounds, arrres, dem_bounds, prods_TOTbbox,
                      dem, lat, lon, hgt_field, prod_list, is_nisar_file=False,
                      outputFormat='ENVI', verbose=None, num_threads='2'):
    """Interpolate and extract 2D metadata layer.
    2D metadata layer is derived by interpolating and then intersecting
    3D layers with a DEM.
    Lat/lon arrays must also be passed for this process.
    """
    ref_geotrans = dem.GetGeoTransform()
    dem_arrres = [abs(ref_geotrans[1]), abs(ref_geotrans[-1])]

    NOHGT_LYRS = ['ionosphere']
    metadatalyr_name = outname.split('/')[-2]
    needs_height_interp = metadatalyr_name not in NOHGT_LYRS

    tmp_name = outname + '.vrt'
    warp_src = tmp_name
    heightsMeta = None
    subset_vrt = None
    ds_frame_marker = osgeo.gdal.Open(tmp_name, osgeo.gdal.GA_ReadOnly)
    source_frame_count = (
        ds_frame_marker.GetMetadataItem('ARIA_FRAME_COUNT')
        if ds_frame_marker is not None else None)
    ds_frame_marker = None

    if needs_height_interp:
        heightsMeta_str = ARIAtools.util.vrt.get_hgt_meta(
            tmp_name, hgt_field)
        if not heightsMeta_str:
            raise RuntimeError(
                f'Missing height metadata {hgt_field!r} in {tmp_name}')
        heightsMeta = np.array(
            heightsMeta_str[1:-1].split(','), dtype='float32')

        dem_min, dem_max = ARIAtools.util.interp._compute_dem_range(dem)

        if not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET'):
            band_indices = ARIAtools.util.interp._get_height_subset_indices(
                heightsMeta, dem_min, dem_max, pad=0)
        else:
            band_indices = np.arange(len(heightsMeta))

        if len(band_indices) < len(heightsMeta):
            LOGGER.debug(
                'Subsetting 3D cube from %d to %d height levels '
                '(DEM range: %.1f to %.1f)',
                len(heightsMeta), len(band_indices), dem_min, dem_max)

            band_list = [int(i + 1) for i in band_indices]
            subset_vrt = outname + '_hsubset.vrt'
            translate_opts = osgeo.gdal.TranslateOptions(
                format='VRT', bandList=band_list)
            ds_sub = osgeo.gdal.Translate(
                subset_vrt, tmp_name, options=translate_opts)
            if ds_sub is not None:
                ds_sub.FlushCache()
            ds_sub = None
            warp_src = subset_vrt
            heightsMeta = heightsMeta[band_indices]

    ds_src = osgeo.gdal.Open(warp_src, osgeo.gdal.GA_ReadOnly)
    if ds_src is None:
        raise RuntimeError(f'Could not open metadata cube {warp_src}')
    src_gt = ds_src.GetGeoTransform()
    src_proj_before_warp = ds_src.GetProjection()
    src_x_edge_2 = src_gt[0] + src_gt[1] * ds_src.RasterXSize
    src_y_edge_2 = src_gt[3] + src_gt[5] * ds_src.RasterYSize
    src_bounds = (
        min(src_gt[0], src_x_edge_2),
        min(src_gt[3], src_y_edge_2),
        max(src_gt[0], src_x_edge_2),
        max(src_gt[3], src_y_edge_2))
    src_xres = abs(src_gt[1])
    src_yres = abs(src_gt[5])
    try:
        src_crs_label = (pyproj.CRS.from_wkt(
            src_proj_before_warp).to_string()
            if src_proj_before_warp else 'unset')
    except pyproj.exceptions.CRSError:
        src_crs_label = 'unrecognized'
    try:
        dem_crs_label = (pyproj.CRS.from_wkt(
            dem.GetProjection()).to_string()
            if dem.GetProjection() else 'unset')
    except pyproj.exceptions.CRSError:
        dem_crs_label = 'unrecognized'
    LOGGER.info(
        'Metadata cube before warp for %s: CRS=%s, '
        'bounds=(%.3f, %.3f, %.3f, %.3f); DEM CRS=%s, '
        'requested bounds=(%.3f, %.3f, %.3f, %.3f)',
        metadatalyr_name, src_crs_label, *src_bounds,
        dem_crs_label, *dem_bounds)
    ds_src = None
    pad_cells = 2
    padded_bounds = [
        dem_bounds[0] - pad_cells * dem_arrres[0],
        dem_bounds[1] - pad_cells * dem_arrres[1],
        dem_bounds[2] + pad_cells * dem_arrres[0],
        dem_bounds[3] + pad_cells * dem_arrres[1],
    ]
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_kwargs = {
            'format': 'MEM', 'outputBounds': padded_bounds,
            'multithread': False}
        if dem.GetProjection():
            # outputBounds are expressed in the DEM CRS.
            warp_kwargs['dstSRS'] = dem.GetProjection()
        warp_options = osgeo.gdal.WarpOptions(**warp_kwargs)
        ds_warp = osgeo.gdal.Warp('', warp_src, options=warp_options)
        if ds_warp is None:
            raise RuntimeError(
                f'GDAL could not warp metadata cube {warp_src} from '
                f'{src_crs_label} to {dem_crs_label}')
        data_array_nodata = ds_warp.GetRasterBand(1).GetNoDataValue()
        data_array = ds_warp.ReadAsArray().astype('float32')
        gt_mem = ds_warp.GetGeoTransform()
        cube_proj = ds_warp.GetProjection()
        x_size = ds_warp.RasterXSize
        y_size = ds_warp.RasterYSize
        ds_warp = None

    if subset_vrt is not None and os.path.exists(subset_vrt):
        os.remove(subset_vrt)

    if data_array.ndim == 2:
        data_array = data_array[np.newaxis, ...]

    if data_array_nodata is not None:
        if np.isnan(data_array_nodata):
            data_array[~np.isfinite(data_array)] = np.nan
        else:
            data_array[data_array == data_array_nodata] = np.nan

    valid_cube_samples = np.isfinite(data_array)
    if not np.any(valid_cube_samples):
        raise RuntimeError(
            f'Warped metadata cube contains no valid samples for {outname}. '
            f'Source CRS={src_crs_label}, source bounds={src_bounds}; '
            f'DEM CRS={dem_crs_label}, requested bounds={tuple(dem_bounds)}. '
            'This indicates incorrect cube georeferencing or no spatial '
            'overlap with the DEM.')

    version_check = []
    for i in prod_list:
        if not is_nisar_file:
            v_num = i.split(':')[-2].split('/')[-1]
            v_num = v_num.split('.nc')[0][-5:]
        else:
            basename = os.path.basename(i.split('"')[1])
            v_num = basename.split('_')[-1][:-3]
            v_num = '.'.join(v_num)
        version_check.append(v_num)
    version_check = min(version_check)

    if ((metadatalyr_name in GEOM_LYRS and version_check < '2_0_4')
            and not is_nisar_file):
        plots_subdir = os.path.abspath(os.path.join(outname, '../..',
                                       'metadatalyr_plots', metadatalyr_name))
        if not os.path.exists(plots_subdir):
            os.makedirs(plots_subdir)

        data_array = MetadataQualityCheck(
            data_array,
            os.path.basename(os.path.dirname(outname)),
            outname,
            verbose).data_array

    if needs_height_interp:
        if data_array.shape[0] != len(heightsMeta):
            raise RuntimeError(
                f'Height-band mismatch for {outname}: raster has '
                f'{data_array.shape[0]} bands but metadata has '
                f'{len(heightsMeta)} heights')

        tmp_name = outname + '_interpolated.tif'

        latitudeMeta = np.linspace(
            gt_mem[3] + gt_mem[5] / 2.0,
            gt_mem[3] + (gt_mem[5] * (y_size - 0.5)),
            y_size, dtype='float32')

        longitudeMeta = np.linspace(
            gt_mem[0] + gt_mem[1] / 2.0,
            gt_mem[0] + (gt_mem[1] * (x_size - 0.5)),
            x_size, dtype='float32')

        dem_band = dem.GetRasterBand(1)
        nodata = dem_band.GetNoDataValue()
        dem_elevation = dem_band.ReadAsArray().astype('float32')
        if nodata is not None:
            if np.isnan(nodata):
                dem_elevation[~np.isfinite(dem_elevation)] = np.nan
            else:
                dem_elevation[dem_elevation == nodata] = np.nan

        # The grids normally share a CRS.  If they do not, transform the DEM
        # x/y coordinates into the metadata cube's CRS before interpolation.
        dem_proj = dem.GetProjection()
        same_crs = True
        if dem_proj and cube_proj:
            same_crs = pyproj.CRS.from_wkt(dem_proj).equals(
                pyproj.CRS.from_wkt(cube_proj))

        if same_crs:
            interp_y, interp_x, interp_z = lat, lon, dem_elevation
        else:
            transformer = pyproj.Transformer.from_crs(
                pyproj.CRS.from_wkt(dem_proj),
                pyproj.CRS.from_wkt(cube_proj), always_xy=True)
            interp_x, interp_y, interp_z = transformer.transform(
                lon, lat, dem_elevation)

        pnts = np.stack(
            (interp_y, interp_x, interp_z), axis=-1).astype('float32')

        # RegularGridInterpolator requires monotonic axes.  Normalize every
        # axis to ascending order and reorder the cube to match.
        if latitudeMeta[0] > latitudeMeta[-1]:
            latitudeMeta = latitudeMeta[::-1]
            data_array = data_array[:, ::-1, :]
        if longitudeMeta[0] > longitudeMeta[-1]:
            longitudeMeta = longitudeMeta[::-1]
            data_array = data_array[:, :, ::-1]
        height_order = np.argsort(heightsMeta)
        heightsMeta = heightsMeta[height_order]
        data_array = data_array[height_order, :, :]
        if len(np.unique(heightsMeta)) != len(heightsMeta):
            raise RuntimeError(
                f'Duplicate height levels found while finalizing {outname}')

        valid_dem_points = (
            np.isfinite(interp_y) & np.isfinite(interp_x)
            & np.isfinite(interp_z))
        if not np.any(valid_dem_points):
            raise RuntimeError(
                f'DEM contains no valid interpolation points for {outname}')

        point_y = interp_y[valid_dem_points]
        point_x = interp_x[valid_dem_points]
        point_z = interp_z[valid_dem_points]
        cube_ranges = (
            float(latitudeMeta[0]), float(latitudeMeta[-1]),
            float(longitudeMeta[0]), float(longitudeMeta[-1]),
            float(heightsMeta[0]), float(heightsMeta[-1]))
        point_ranges = (
            float(np.nanmin(point_y)), float(np.nanmax(point_y)),
            float(np.nanmin(point_x)), float(np.nanmax(point_x)),
            float(np.nanmin(point_z)), float(np.nanmax(point_z)))
        LOGGER.info(
            'Interpolation coordinates for %s: cube y=(%.3f, %.3f), '
            'x=(%.3f, %.3f), height=(%.3f, %.3f); DEM/query '
            'y=(%.3f, %.3f), x=(%.3f, %.3f), height=(%.3f, %.3f); '
            'finite cube samples=%d/%d',
            metadatalyr_name, *cube_ranges, *point_ranges,
            int(valid_cube_samples.sum()), valid_cube_samples.size)

        ranges_overlap = (
            point_ranges[1] >= cube_ranges[0]
            and point_ranges[0] <= cube_ranges[1]
            and point_ranges[3] >= cube_ranges[2]
            and point_ranges[2] <= cube_ranges[3]
            and point_ranges[5] >= cube_ranges[4]
            and point_ranges[4] <= cube_ranges[5])
        if not ranges_overlap:
            raise RuntimeError(
                f'DEM/query coordinates do not overlap the metadata cube '
                f'for {outname}. Cube ranges (ymin, ymax, xmin, xmax, '
                f'hmin, hmax)={cube_ranges}; query ranges={point_ranges}')

        interper = scipy.interpolate.RegularGridInterpolator(
            (latitudeMeta, longitudeMeta, heightsMeta),
            data_array.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        out_interpolated = interper(pnts)

        valid_interpolated = np.isfinite(out_interpolated)
        if not np.any(valid_interpolated):
            raise RuntimeError(
                f'Height interpolation produced no valid pixels for '
                f'{outname}; check the DEM and cube CRS/bounds')

        LOGGER.info(
            'Finalized %s interpolation: %d/%d valid pixels, '
            'range %.6g to %.6g',
            metadatalyr_name, int(valid_interpolated.sum()),
            out_interpolated.size,
            float(np.nanmin(out_interpolated)),
            float(np.nanmax(out_interpolated)))

        ARIAtools.util.vrt.renderVRT(
            tmp_name, out_interpolated, geotrans=dem.GetGeoTransform(),
            drivername='GTiff',
            gdal_fmt='float32',
            proj=dem.GetProjection(), nodata=np.nan)
        out_interpolated = None

    dem_crop = outname + '_demcrop.tif'
    if os.path.exists(dem_crop):
        os.remove(dem_crop)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': 'GTiff', 'cutlineDSName': prods_TOTbbox,
            'outputBounds': dem_bounds, 'srcNodata': np.nan,
            'dstNodata': np.nan,
            'xRes': dem_arrres[0], 'yRes': dem_arrres[1],
            'targetAlignedPixels': True, 'multithread': False,
            'creationOptions': [
                'TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=3',
                'BIGTIFF=IF_SAFER']}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds_crop1 = osgeo.gdal.Warp(
            dem_crop, tmp_name, options=warp_options
        )
        if ds_crop1 is None:
            raise RuntimeError(
                f'Failed to crop interpolated metadata to DEM bounds: '
                f'{outname}')
        ds_crop1 = None
    validate_metadata_staging_raster(
        dem_crop, f'{metadatalyr_name} DEM crop', expected_bands=1)

    if needs_height_interp:
        for tmp_suffix in ['', '.vrt', '.aux.xml']:
            tmp_path = tmp_name + tmp_suffix
            if os.path.exists(tmp_path):
                os.remove(tmp_path)

    for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        f_to_rm = f"{outname}{ext}"
        if os.path.exists(f_to_rm):
            os.remove(f_to_rm)

    final_driver = outputFormat
    if osgeo.gdal.GetDriverByName(final_driver) is None:
        raise RuntimeError(
            f'Requested GDAL output driver is unavailable: {final_driver}')

    # GDAL's ISCE writer cannot reliably serve as a gdal.Warp destination
    # (IReadBlock/scanline failures occur even for an otherwise valid source).
    # Perform all spatial operations in GeoTIFF, then translate the completed
    # single-band raster to the user-requested physical format.
    final_stage = (outname if final_driver.upper() == 'GTIFF'
                   else outname + '_final_stage.tif')
    for stage_suffix in ['', '.vrt', '.aux.xml']:
        stage_path = final_stage + stage_suffix
        if final_stage != outname and os.path.exists(stage_path):
            os.remove(stage_path)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': 'GTiff', 'cutlineDSName': prods_TOTbbox,
            'outputBounds': bbox_bounds, 'srcNodata': np.nan,
            'dstNodata': np.nan,
            'xRes': arrres[0], 'yRes': arrres[1], 'targetAlignedPixels': True,
            'multithread': False,
            'creationOptions': [
                'TILED=YES', 'COMPRESS=DEFLATE', 'PREDICTOR=3',
                'BIGTIFF=IF_SAFER']}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds_crop2 = osgeo.gdal.Warp(
            final_stage, dem_crop, options=warp_options
        )
        if ds_crop2 is None:
            raise RuntimeError(
                f'Failed to create finalized metadata staging raster: '
                f'{final_stage}')
        ds_crop2.FlushCache()
        ds_crop2 = None
    validate_metadata_staging_raster(
        final_stage, f'{metadatalyr_name} final', expected_bands=1)

    if final_driver.upper() != 'GTIFF':
        translate_options = osgeo.gdal.TranslateOptions(format=final_driver)
        ds_final = osgeo.gdal.Translate(
            outname, final_stage, options=translate_options)
        if ds_final is None:
            raise RuntimeError(
                f'Failed to translate finalized metadata raster to '
                f'{final_driver}: {outname}')
        ds_final.FlushCache()
        ds_final = None

        for stage_suffix in ['', '.vrt', '.aux.xml']:
            stage_path = final_stage + stage_suffix
            if os.path.exists(stage_path):
                os.remove(stage_path)

    for crop_suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        crop_path = dem_crop + crop_suffix
        if os.path.exists(crop_path):
            os.remove(crop_path)

    # Reject all-NoData outputs instead of silently publishing blank layers.
    ds_target = osgeo.gdal.Open(outname, osgeo.gdal.GA_ReadOnly)
    if ds_target is None:
        raise RuntimeError(f'Could not reopen finalized raster {outname}')

    actual_final_driver = ds_target.GetDriver().ShortName
    if actual_final_driver.upper() != final_driver.upper():
        ds_target = None
        raise RuntimeError(
            f'Finalized raster driver mismatch for {outname}: requested '
            f'{final_driver}, created {actual_final_driver}')

    final_array = ds_target.GetRasterBand(1).ReadAsArray()
    final_valid = np.isfinite(final_array)
    if not np.any(final_valid):
        ds_target = None
        raise RuntimeError(
            f'Finalized metadata raster contains no valid pixels: {outname}')

    LOGGER.info(
        'Final metadata raster %s (%s): %d/%d valid pixels, '
        'range %.6g to %.6g',
        outname, actual_final_driver, int(final_valid.sum()), final_array.size,
        float(np.nanmin(final_array)), float(np.nanmax(final_array)))
    final_array = None

    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    vrt_ds = osgeo.gdal.Translate(
        outname + '.vrt', ds_target, options=translate_options)
    if vrt_ds is None:
        ds_target = None
        raise RuntimeError(f'Could not create VRT for {outname}')
    if source_frame_count is not None:
        vrt_ds.SetMetadataItem('ARIA_FRAME_COUNT', source_frame_count)
    vrt_ds.FlushCache()
    vrt_ds = None
    ds_target = None

    data_array = None
    dem = None
    lat = None
    lon = None

    return


def transformPoints(lats: np.ndarray, lons: np.ndarray, hgts: np.ndarray,
                    old_proj: pyproj.CRS, new_proj: pyproj.CRS) -> np.ndarray:
    '''
    Transform lat/lon/hgt data to an array of points in a new
    projection
    '''
    transformer = pyproj.Transformer.from_crs(old_proj, new_proj)

    if not isinstance(new_proj, pyproj.CRS):
        new_proj = pyproj.CRS.from_epsg(new_proj.lstrip('EPSG:'))
    if not isinstance(old_proj, pyproj.CRS):
        old_proj = pyproj.CRS.from_epsg(old_proj.lstrip('EPSG:'))

    in_flip = old_proj.axis_info[0].direction
    out_flip = new_proj.axis_info[0].direction

    if in_flip == 'east':
        res = transformer.transform(lons, lats, hgts)
    else:
        res = transformer.transform(lats, lons, hgts)

    if out_flip == 'east':
        return np.stack((res[1], res[0], res[2]), axis=-1, dtype='float32').T
    else:
        return np.stack(res, axis=-1, dtype='float32').T
