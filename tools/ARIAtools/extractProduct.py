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

import ARIAtools.product
import ARIAtools.util.ionosphere
import ARIAtools.util.interp
import ARIAtools.util.vrt
import ARIAtools.util.shp
import ARIAtools.util.misc
import ARIAtools.util.seq_stitch

from ARIAtools.constants import ARIA_PX_SIZES

LOGGER = logging.getLogger(__name__)
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
        
        # --- CLOSE THE CLASS-LEVEL POINTER ---
        # Capture the dataset to return it, then explicitly sever 
        # the class's internal link to the GDAL memory object.
        safe_return_array = self.data_array
        self.data_array = None
        
        return safe_return_array


def crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs):
    """
    Manage cropping of existing, extracted layers
    """
    LOGGER.debug('Cropping %s - %s', ifg_tag, lyrname)

    # Crop
    gdal_warp_kwargs['format'] = 'ENVI'
    warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
    ds = osgeo.gdal.Warp(
        outname + '_crop', outname + '.vrt', options=warp_options
    )
    ds = None

    for crop_name in glob.glob(outname + '_crop*'):
        fname = os.path.basename(crop_name).replace('_crop', '')
        fname = os.path.join(os.path.dirname(crop_name), fname)
        os.rename(crop_name, fname)

    # Update VRT
    ds_trans = osgeo.gdal.Translate(outname + '.vrt', outname, format='VRT')
    ds_trans = None

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
    # If specified workdir doesn't exist, create it
    os.makedirs(workdir, exist_ok=True)

    # determine if NISAR GUNW
    is_nisar_file = False
    track_fileext = product_dict[0]['unwrappedPhase'][0].split('"')[1]
    if track_fileext.endswith('.h5'):
        is_nisar_file = True

    # If specified, check if user's bounding box meets minimum threshold area
    lyr_proj = int(metadata_dict[0]['projection'][0])
    if bbox_file is not None:
        user_bbox = ARIAtools.util.shp.open_shp(bbox_file)
        overlap_area = ARIAtools.util.shp.shp_area(user_bbox, lyr_proj)
        if overlap_area < minimumOverlap:
            raise Exception(f"User bound box {bbox_file} has an area of only "
                            f"{overlap_area}km\u00b2, below specified "
                            f"minimum threshold area "
                            f"{minimumOverlap}km\u00b2")

    # Check if product bounding box exists from previous run
    prods_TOTbbox = os.path.join(workdir, 'productBoundingBox.json')
    prods_TOTbbox_metadatalyr = os.path.join(
        workdir, 'productBoundingBox_croptounion_formetadatalyr.json')
    if os.path.exists(prods_TOTbbox) \
        and os.path.exists(prods_TOTbbox_metadatalyr):
        exist_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        exist_metadatalyr = \
            ARIAtools.util.shp.open_shp(prods_TOTbbox_metadatalyr)

        # Save copy of file to disk
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

    # Extract/merge productBoundingBox layers
    for scene in product_dict:

        # Get pair name, expected in dictionary
        pair_name = scene["pair_name"][0]
        outname = os.path.join(workdir, pair_name + '.json')
        if os.path.exists(outname):
            os.remove(outname)

        # Create union of productBoundingBox layers
        for prods_bbox in scene["productBoundingBox"]:
            if os.path.exists(outname):
                union_bbox = ARIAtools.util.shp.open_shp(outname)
                prods_bbox = prods_bbox.union(union_bbox)
            ARIAtools.util.shp.save_shp(
                outname, prods_bbox, lyr_proj)
        scene["productBoundingBox"] = [outname]

    prods_TOTbbox = os.path.join(workdir, 'productBoundingBox.json')

    # Need to track bounds of max extent
    # to avoid metadata interpolation issues
    sceneareas = [
        ARIAtools.util.shp.open_shp(i['productBoundingBox'][0]).area
        for i in product_dict]
    ind_max_area = sceneareas.index(max(sceneareas))
    product_bbox = product_dict[ind_max_area]['productBoundingBox'][0]
    ARIAtools.util.shp.save_shp(
        prods_TOTbbox_metadatalyr, ARIAtools.util.shp.open_shp(product_bbox),
        lyr_proj)

    # Initiate intersection file with bbox, if bbox specified
    if bbox_file is not None:
        ARIAtools.util.shp.save_shp(
            prods_TOTbbox, ARIAtools.util.shp.open_shp(bbox_file),
            lyr_proj)

    # Initiate intersection with largest scene, if bbox NOT specified
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
        # Generate footprint for the union of all products
        if croptounion:
            # Get union
            total_bbox = total_bbox.union(prods_bbox)
            total_bbox_metadatalyr = total_bbox_metadatalyr.union(prods_bbox)

            # Save to file
            ARIAtools.util.shp.save_shp(
                prods_TOTbbox, total_bbox, lyr_proj)
            ARIAtools.util.shp.save_shp(
                prods_TOTbbox_metadatalyr, total_bbox_metadatalyr,
                lyr_proj)

        # Generate footprint for the common intersection of all products
        else:
            # Now pass track intersection for cutline
            prods_bbox = prods_bbox.intersection(total_bbox)

            # Estimate percentage of overlap with bbox
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

                # Kick out scenes below specified overlap threshold
                if overlap_area < minimumOverlap:
                    LOGGER.debug(f'Rejected scene {scene_obj} has only '
                                 f'{overlap_area}km\u00b2 overlap with bbox')
                    rejected_scenes.append(product_dict.index(scene))
                    os.remove(scene_obj)

                else:
                    ARIAtools.util.shp.save_shp(
                        prods_TOTbbox, prods_bbox, lyr_proj)

                    # Need to track bounds of max extent
                    # to avoid metadata interpolation issues
                    total_bbox_metadatalyr = total_bbox_metadatalyr.union(
                        ARIAtools.util.shp.open_shp(
                            scene['productBoundingBox'][0]))
                    ARIAtools.util.shp.save_shp(
                        prods_TOTbbox_metadatalyr, total_bbox_metadatalyr,
                        lyr_proj)

    # Remove scenes with insufficient overlap w.r.t. bbox
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

    # If bbox specified, intersect with common track intersection/union
    if bbox_file is not None:
        user_bbox = ARIAtools.util.shp.open_shp(bbox_file)
        total_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        user_bbox = user_bbox.intersection(total_bbox)
        ARIAtools.util.shp.save_shp(
            prods_TOTbbox, user_bbox, lyr_proj)

    else:
        bbox_file = prods_TOTbbox

    # Compare current bbox to existing bbox
    if exist_bbox and runlog:
        exist_area = ARIAtools.util.shp.shp_area(exist_bbox, lyr_proj)

        # Calculate overlap area
        new_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
        new_area = ARIAtools.util.shp.shp_area(new_bbox, lyr_proj)
        area_ratio = new_area / exist_area

        olap_bbox = exist_bbox.intersection(new_bbox)
        olap_area = ARIAtools.util.shp.shp_area(olap_bbox, lyr_proj)

        # Compare areas
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
            # Same bbox within machine precision
            update_mode = 'skip'
        elif olap_ratio < 1.0:
            # For smaller bbox, need to crop
            update_mode = 'crop_only'
        else:
            # If no prior products exist, or new AOI is larger
            update_mode = 'full_extract'
        runlog.update('update_mode', update_mode)

        LOGGER.info('Update mode: %s', update_mode)

    # Warp the first scene with the output-bounds defined above
    # ensure output-bounds are an integer multiple of interferometric grid
    # and adjust if necessary
    OG_bounds = list(
        ARIAtools.util.shp.open_shp(bbox_file).bounds)
    gdal_warp_kwargs = {
        'format': 'MEM', 'multithread': True, 'dstSRS': f'EPSG:{lyr_proj}'}
    ds_vrt = osgeo.gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', ds_vrt, options=warp_options)
        arrres = [abs(ds.GetGeoTransform()[1]),
                  abs(ds.GetGeoTransform()[-1])]
        ds = None
        ds_vrt = None

    # Adjust arrres to supported resolution
    for i, res in enumerate(arrres):
        res_ndx = np.argmin([
            np.abs(res - px_size) for px_size in ARIA_PX_SIZES
        ])
        arrres[i] = ARIA_PX_SIZES[res_ndx]

    # warp again with fixed transform and bounds
    gdal_warp_kwargs['outputBounds'] = OG_bounds
    gdal_warp_kwargs['xRes'] = arrres[0]
    gdal_warp_kwargs['yRes'] = arrres[1]
    gdal_warp_kwargs['targetAlignedPixels'] = True
    ds_vrt = osgeo.gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', ds_vrt, options=warp_options)

        # Get shape of full res layers
        arrshape = [ds.RasterYSize, ds.RasterXSize]
        ds_gt = ds.GetGeoTransform()
        new_bounds = [ds_gt[0], ds_gt[3] + (ds_gt[-1] * arrshape[0]),
                      ds_gt[0] + (ds_gt[1] * arrshape[1]), ds_gt[3]]

        if OG_bounds != new_bounds:
            # Use shapely to make list
            user_bbox = shapely.geometry.Polygon(np.column_stack((
                np.array([new_bounds[0], new_bounds[2], new_bounds[2],
                          new_bounds[0], new_bounds[0]]),
                np.array([new_bounds[1], new_bounds[1], new_bounds[3],
                          new_bounds[3], new_bounds[1]]))))

            # Save polygon in shapefile
            bbox_file = os.path.join(
                os.path.dirname(workdir), 'user_bbox.json')
            ARIAtools.util.shp.save_shp(
                bbox_file, user_bbox, lyr_proj)
            total_bbox = ARIAtools.util.shp.open_shp(prods_TOTbbox)
            user_bbox = user_bbox.intersection(total_bbox)
            ARIAtools.util.shp.save_shp(
                prods_TOTbbox, user_bbox, lyr_proj)

        # Get projection of full res layers
        proj = ds.GetProjection()
        ds = None
        ds_vrt = None

    # Run additional checks and update runlog if provided
    if runlog is None:
        update_mode = 'full_extract'
    else:
        log_data = runlog.load()
        update_mode = log_data['update_mode']

        # Check other parameters
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

        # Update log
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
    """Wrapper to create raster and apply projection using Rioxarray (Safe)"""

    # --- Height-based band subsetting optimisation ---
    # If a DEM is provided, subset to only the height bands that span the
    # DEM elevation range *before* the expensive remote I/O Warp.
    subset_vrts = []
    effective_data_lis = data_lis
    subsetted_heightsMeta = None

    if (dem is not None and hgt_field
            and not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET')):
        try:
            heightsMeta_str = ARIAtools.util.vrt.get_hgt_meta(
                data_lis[0], hgt_field)
            if heightsMeta_str:
                heightsMeta = np.array(
                    heightsMeta_str[1:-1].split(','), dtype='float32')
                dem_min, dem_max = (
                    ARIAtools.util.interp._compute_dem_range(dem))
                band_indices = (
                    ARIAtools.util.interp._get_height_subset_indices(
                        heightsMeta, dem_min, dem_max, pad=1))

                if len(band_indices) < len(heightsMeta):
                    band_list = [int(i + 1) for i in band_indices]
                    translate_opts = osgeo.gdal.TranslateOptions(
                        format='VRT', bandList=band_list)
                    subsetted_data = []
                    for idx, src in enumerate(data_lis):
                        sub_vrt = fname + f'_src{idx}_hsubset.vrt'
                        ds_sub = osgeo.gdal.Translate(
                            sub_vrt, src, options=translate_opts)
                        ds_sub = None
                        subset_vrts.append(sub_vrt)
                        subsetted_data.append(sub_vrt)
                    effective_data_lis = subsetted_data
                    subsetted_heightsMeta = heightsMeta[band_indices]
        except Exception:
            pass  # fall back to reading all bands

    # 1) Build a lightweight reference warp (VRT)
    ref_vrt = fname + "_ref.vrt"
    ds = osgeo.gdal.Warp(
        ref_vrt,
        effective_data_lis[0],
        format='VRT',
        dstSRS=proj,
        dstNodata=np.nan,
        multithread=True
    )
    ds = None # Close immediately

    # 2) Read the derived resolution
    ds = osgeo.gdal.Open(ref_vrt, osgeo.gdal.GA_ReadOnly)
    gt = ds.GetGeoTransform()
    xres, yres = gt[1], abs(gt[5])
    ds = None # Close immediately

    # 3) Warp + mosaic to temp Tiff
    mosaic_tif = fname + "_warp.tif"
    ds = osgeo.gdal.Warp(
        mosaic_tif,
        effective_data_lis,
        format='GTiff',
        xRes=xres, yRes=yres,
        dstSRS=proj,
        dstNodata=np.nan,
        multithread=True,
        creationOptions=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=IF_SAFER"]
    )
    ds = None # Close immediately

    # 3b) Fix stale NETCDF dimension metadata after band subsetting.
    #     gdal.Warp propagates the original height dimension metadata
    #     (e.g. 20 values) even though the data now has fewer bands.
    #     rioxarray uses this metadata to build coordinates, causing a
    #     dimension mismatch.  Update it to match the actual band count.
    if subsetted_heightsMeta is not None:
        ds_fix = osgeo.gdal.Open(mosaic_tif, osgeo.gdal.GA_Update)
        if ds_fix is not None:
            meta = ds_fix.GetMetadata()
            for key in list(meta.keys()):
                if key.startswith('NETCDF_DIM_') and key.endswith('_VALUES'):
                    dim_name = key[len('NETCDF_DIM_'):-len('_VALUES')]
                    new_vals = '{' + ','.join(
                        str(h) for h in subsetted_heightsMeta) + '}'
                    ds_fix.SetMetadataItem(key, new_vals)
                    def_key = f'NETCDF_DIM_{dim_name}_DEF'
                    if def_key in meta:
                        old_def = meta[def_key]
                        dtype_str = old_def.strip('{}[]').split(',')[-1]
                        ds_fix.SetMetadataItem(
                            def_key,
                            '{' + str(len(subsetted_heightsMeta)) +
                            ',' + dtype_str + '}')
            ds_fix.FlushCache()
            ds_fix = None

    # 4) Open with rioxarray (Context Manager prevents locking)
    with rioxarray.open_rasterio(mosaic_tif, masked=True) as da:
        da = da.rio.write_nodata(np.nan, encoded=True)
        
        # --- FIX: Remove _FillValue from attrs ---
        if "_FillValue" in da.attrs:
            del da.attrs["_FillValue"]
        # -----------------------------------------
        
        # Flip the sign for NISAR convention
        if sign_multiplier == -1:
            da = da * -1
            
        # Enforce threading during the write
        with rasterio.Env(GDAL_NUM_THREADS='ALL_CPUS'): 
            da.rio.to_raster(fname, driver=driver, crs=proj)

    # 5) Clean up (Now safe because 'da' is closed)
    if os.path.exists(mosaic_tif):
        os.remove(mosaic_tif)
    if os.path.exists(ref_vrt):
        os.remove(ref_vrt)

    # 6) Create VRT file
    buildvrt_options = osgeo.gdal.BuildVRTOptions(outputSRS=proj)
    ds_vrt = osgeo.gdal.BuildVRT(
        fname + '.vrt', fname, options=buildvrt_options
    )
    ds_vrt = None # Close immediately

    # 7) Add height info
    if hgt_field is not None:
        if subsetted_heightsMeta is not None:
            # Write the subsetted height values
            hgt_meta = '{' + ','.join(
                str(h) for h in subsetted_heightsMeta) + '}'
        else:
            hgt_meta = ARIAtools.util.vrt.get_hgt_meta(
                data_lis[0], hgt_field)
        
        ds_meta_update = osgeo.gdal.Open(
            fname + '.vrt', osgeo.gdal.GA_Update
        )
        ds_meta_update.SetMetadataItem(hgt_field, hgt_meta)
        ds_meta_update = None # Close immediately

    # Clean up temporary subset VRTs
    for v in subset_vrts:
        if os.path.exists(v):
            os.remove(v)

    return


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

    # ionosphere layer, heights do not exist to exit
    if metadata_arr[0].split('/')[-1] == 'ionosphere':
        ds_vrt = osgeo.gdal.BuildVRT(outname + '.vrt', metadata_arr)
        ds_vrt= None
        return [0], None, outname

    # capture model if tropo product
    if 'tropo' in layer:
        if not is_nisar_file:
            out_dir = os.path.join(out_dir, model_name)
        outname = os.path.join(out_dir, ifg)

        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    # Get height values
    ds_meta = osgeo.gdal.Open(metadata_arr[0])
    zdim = ds_meta.GetMetadataItem('NETCDF_DIM_EXTRA')[1:-1]
    ds_meta = None # Close
    hgt_field = f'NETCDF_DIM_{zdim}_VALUES'

    # Check if height layers are consistent
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
        # get ref and sec paths
        if not is_nisar_file:
            date_dir = os.path.join(out_dir, 'dates')
            if not os.path.exists(date_dir):
                os.mkdir(date_dir)

            ref_outname = os.path.join(date_dir, ifg.split('_')[0])
            sec_outname = os.path.join(date_dir, ifg.split('_')[1])
            ref_str = 'reference/' + layer
            sec_str = 'secondary/' + layer
            sec_metadata_arr = [
                i[:-len(ref_str)] + sec_str for i in metadata_arr]
            tup_outputs = [
                (ref_outname, metadata_arr), (sec_outname, sec_metadata_arr)]

        else:
            ref_outname = os.path.join(out_dir, ifg)
            sec_outname = ref_outname
            tup_outputs = [(ref_outname, metadata_arr)]

        # write ref and sec files
        for i in tup_outputs:

            # delete temporary files to circumvent potential inconsistent dims
            for j in glob.glob(i[0] + '*'):
                if os.path.isfile(j):
                    os.remove(j)
            create_raster_from_gunw(i[0], i[1], proj, driver, hgt_field,
                sign_multiplier, dem=dem)

        if not is_nisar_file:
            # compute differential
            generate_diff(
                ref_outname, sec_outname, outname, layer, layer, False,
                hgt_field, proj, driver, dem=dem)

        # write raster to file if it does not exist
        if layer in layers:
            for i in [ref_outname, sec_outname]:
                if not os.path.exists(i):
                    create_raster_from_gunw(i, [i], proj, driver, hgt_field)

    else:
        if not os.path.exists(outname + '.vrt'):
            if is_nisar_file:
                # Need to compute azimuthAngle from
                # losUnitVectorX and losUnitVectorY
                if layer == 'azimuthAngle':
                    losx_arr = copy.deepcopy(metadata_arr)
                    losy_arr = [
                        path.replace('losUnitVectorX', 'losUnitVectorY')
                        for path in metadata_arr
                    ]

                    # Initiate temp los dimension files
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

                    # Get NoData value from input
                    ds_temp = osgeo.gdal.Open(losx_name + '.vrt')
                    src_nodata = ds_temp.GetRasterBand(1).GetNoDataValue()
                    ds_temp = None

                    # Construct the Calc String safely
                    # If src_nodata exists and is NOT nan, we must mask it
                    # manually.
                    if src_nodata is not None and not np.isnan(src_nodata):
                        # Logic: If (A == nodata) OR (B == nodata), return
                        # NaN. Else calculate the angle.
                        calc_cmd = (
                            f"numpy.where("
                            f"(A=={src_nodata})|(B=={src_nodata}), "
                            f"numpy.nan, "
                            f"numpy.degrees(numpy.arctan2(-B, -A)))"
                        )
                    else:
                        # If input is already NaN (or None), the math
                        # handles it naturally.
                        calc_cmd = "numpy.degrees(numpy.arctan2(-B, -A))"

                    # Compute azimuthAngle from the outputs above
                    # We map path_x to 'A' and path_y to 'B'
                    osgeo_utils.gdal_calc.Calc(
                        A=losx_name + '.vrt',
                        B=losy_name + '.vrt',
                        outfile=outname,
                        calc=calc_cmd,
                        format=driver,
                        allBands="A",  # processes every band
                        quiet=True
                    )

                    # Manually enforce NoData = NaN on the output header
                    ds_update = osgeo.gdal.Open(
                        outname, osgeo.gdal.GA_Update
                    )
                    if ds_update:
                        for b in range(1, ds_update.RasterCount + 1):
                            ds_update.GetRasterBand(b).SetNoDataValue(np.nan)
                        ds_update = None

                    # Create VRT file
                    buildvrt_options = osgeo.gdal.BuildVRTOptions(
                        outputSRS=proj
                    )
                    ds_vrt = osgeo.gdal.BuildVRT(
                        outname + '.vrt',
                        outname,
                        options=buildvrt_options
                    )
                    ds_vrt = None

                    # Add height info
                    if hgt_field is not None:
                        # Fetch height from the LOCAL file
                        ds_meta = osgeo.gdal.Open(losx_name + '.vrt')
                        hgt_meta = ds_meta.GetMetadataItem(hgt_field)
                        ds_meta = None  # Close file

                        # Also added GA_Update so GDAL does not fail
                        # silently when saving metadata)
                        ds_vrt = osgeo.gdal.Open(outname + '.vrt', osgeo.gdal.GA_Update)
                        ds_vrt.SetMetadataItem(hgt_field, hgt_meta)
                        ds_vrt = None  # Close file

                    # Cleanup Input Files
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
                
                # Get metadata from source safely
                ds_src = osgeo.gdal.Open(metadata_arr[0])
                hgt_val = ds_src.GetMetadataItem(hgt_field)
                ds_src = None # Close source

                # Set metadata on VRT and THEN close
                ds_vrt.SetMetadataItem(hgt_field, hgt_val)
                ds_vrt = None # Close VRT

    return hgt_field, ref_outname


def generate_diff(ref_outname, sec_outname, outname, key, OG_key, tropo_total,
                  hgt_field, proj, driver, sign_multiplier=1, dem=None):
    """ Compute differential from reference and secondary scenes (Multi-dim safe) """

    # if specified workdir doesn't exist, create it
    output_dir = os.path.dirname(outname)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # --- Height-based band subsetting optimisation ---
    # Subset to only the height bands spanning the DEM range.
    subset_vrts = []  # track temp VRTs for cleanup
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
                        heightsMeta, dem_min, dem_max, pad=1))

                if len(band_indices) < len(heightsMeta):
                    band_list = [int(i + 1) for i in band_indices]
                    translate_opts = osgeo.gdal.TranslateOptions(
                        format='VRT', bandList=band_list)

                    for src_path, label in [
                            (ref_vrt_path, 'ref'), (sec_vrt_path, 'sec')]:
                        sub_vrt = outname + f'_{label}_hsubset.vrt'
                        ds_sub = osgeo.gdal.Translate(
                            sub_vrt, src_path, options=translate_opts)
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
            pass  # fall back to reading all bands

    # 1. Open Inputs with Context Managers (Closes files automatically)
    with rioxarray.open_rasterio(sec_vrt_path, masked=True) as da_sec:
        # Copy attributes and crs while file is open
        sec_attrs = da_sec.attrs
        sec_crs = da_sec.rio.crs
        sec_nodata = da_sec.rio.nodata
        
        # Open Reference inside the first block or separately
        with rioxarray.open_rasterio(ref_vrt_path, masked=True) as da_ref:
            arr_ref = da_ref.data
            arr_sec = da_sec.data # Read data while open

            # 2. Math (Preserves dimensions)
            if tropo_total:
                arr_total = arr_sec + arr_ref
            else:
                arr_total = arr_sec - arr_ref

            if sign_multiplier == -1:
                arr_total = arr_total * -1

            # 3. Create Output DataArray
            # We copy da_sec to preserve coordinates/dims/attrs
            da_total = da_sec.copy()
            da_total.data = arr_total

            # Update attributes
            da_total.name = key
            og_da_attrs = sec_attrs
            da_attrs = {}
            for k in og_da_attrs:
                new_k = k.replace(OG_key, key)
                new_v = og_da_attrs[k]
                if isinstance(new_v, str):
                    new_v = new_v.replace(OG_key, key)
                da_attrs[new_k] = new_v
            da_total = da_total.assign_attrs(da_attrs)
            
            # Ensure CRS/Nodata is carried over
            if sec_crs:
                da_total.rio.write_crs(sec_crs, inplace=True)
            if sec_nodata is not None:
                da_total.rio.write_nodata(sec_nodata, inplace=True)

            # --- FIX: Remove _FillValue from attrs to prevent conflict ---
            if "_FillValue" in da_total.attrs:
                del da_total.attrs["_FillValue"]
            # -------------------------------------------------------------

            # 4. Write to disk
            # Using rasterio.Env to ensure threading settings are respected
            with rasterio.Env(GDAL_NUM_THREADS='ALL_CPUS'):
                da_total.rio.to_raster(outname, driver=driver, crs=proj)

    # 5. Build VRT (Pure GDAL)
    buildvrt_options = osgeo.gdal.BuildVRTOptions(outputSRS=proj)
    ds_vrt = osgeo.gdal.BuildVRT(
        f'{outname}.vrt', outname, options=buildvrt_options
    )

    # Fix numpy array attributes for VRT metadata
    if hgt_field in da_attrs:
        if not isinstance(da_attrs[hgt_field], (list, tuple)):
             if isinstance(da_attrs[hgt_field], np.ndarray):
                 da_attrs[hgt_field] = da_attrs[hgt_field].tolist()
             else:
                 # Fallback if it's a scalar or something else
                 pass
                 
    ds_vrt.SetMetadata(da_attrs)
    ds_vrt = None

    # Clean up temporary subset VRTs
    for v in subset_vrts:
        if os.path.exists(v):
            os.remove(v)

    return


def extract_bperp_dict(products, num_threads):
    """Extracts bPerpendicular mean over frames for each product in products"""

    def read_and_average_bperp(frame):
        """Helper function for dask multiprocessing"""

        # Re-authenticate GDAL for the isolated worker process
        ARIAtools.product._configure_gdal_virtual_access()

        # 1. Open explicitly
        ds = osgeo.gdal.Open(frame, osgeo.gdal.GA_ReadOnly)
        
        # 2. Read data and get nodata value
        arr = ds.ReadAsArray().astype(float)
        nodata = ds.GetRasterBand(1).GetNoDataValue()
        
        # 3. CRITICAL: Close the file explicitly
        ds = None 
        
        # 4. Replace nodata with NaN (if nodata is not already NaN)
        if nodata is not None and not np.isnan(nodata):
            arr = np.where(arr == nodata, np.nan, arr)
        
        # 5. Take mean ignoring NaN values
        res = np.nanmean(arr)
        
        return res

    bperp_dict = {}
    for product in products:
        jobs = []
        for frame in product['bPerpendicular']:
            jobs.append(dask.delayed(read_and_average_bperp)(frame))

        mean_bperp_by_frames = dask.compute(
            jobs, num_workers=int(num_threads), scheduler='processes')[0]

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
    # Depending on type, set sec/ref output dirs
    if key == 'troposphereTotal':
        layers.append(key)
        sec_workdir = os.path.join(os.path.dirname(workdir),
                                   sec_key)
        ref_workdir = os.path.join(os.path.dirname(workdir),
                                   ref_key)

    else:
        sec_workdir = copy.deepcopy(workdir)
        ref_workdir = copy.deepcopy(workdir)

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    all_workdirs = [workdir, sec_workdir, ref_workdir]
    all_workdirs = list(set(all_workdirs))
    existing_outputs = track_correction_outputs(all_workdirs)

    # update existing outputs, if necessary
    if update_mode == 'crop_only' and existing_outputs != []:
        for outname in existing_outputs:
            ifg_tag = os.path.basename(outname).split('.vrt')[0]
            crop_only_manager(outname[:-4], key, ifg_tag, gdal_warp_kwargs)

        return existing_outputs

    # If specified workdirs do not exist, create them
    for i in all_workdirs:
        if not os.path.exists(i):
            os.mkdir(i)

    # Flip sign for external corrections if NISAR
    sign_multiplier = -1 if (is_nisar_file and key in [
        'solidEarthTide', 'troposphereWet', 
        'troposphereHydrostatic', 'troposphereTotal']) else 1

    # Log height subsetting info once for this layer group
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
                # VRTs not yet created; use first product to peek at heights
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
                    hts, dem_min, dem_max, pad=1)
                if len(idx) < len(hts):
                    LOGGER.info(
                        'Height subsetting %s: %d → %d levels '
                        '(DEM range: %.0f to %.0f m)',
                        key, len(hts), len(idx), dem_min, dem_max)
        except Exception:
            pass

    # Iterate through all IFGs
    all_outputs = []
    prog_bar = ARIAtools.util.misc.ProgressBar(
        maxValue=len(product_dict[0]), prefix=f'Exporting {key}: '
    )
    for i in enumerate(product_dict[0]):
        ifg = product_dict[1][i[0]][0]
        outname = os.path.abspath(os.path.join(workdir, ifg))

        # capture model if tropo product
        model_name = None
        if 'tropo' in key:
            out_dir = os.path.dirname(outname)
            if not is_nisar_file:
                model_name = i[1][0].split('/')[-3]
                out_dir = os.path.join(out_dir, model_name)
                outname = os.path.join(out_dir, ifg)
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

        # skip if product exists
        if os.path.exists(outname):
            continue

        # create temp files for ref/sec components
        if ref_key in user_lyrs or tropo_total:
            ref_outname = os.path.abspath(os.path.join(ref_workdir, ifg))
            hgt_field, ref_outname = prep_metadatalayers(
                ref_outname, i[1], dem, ref_key, layers, is_nisar_file, proj,
                outputFormat, model_name)

        # record output directories
        if model_name is not None:
            all_outputs.append(os.path.join(workdir, model_name))
            all_outputs.append(os.path.join(ref_workdir, model_name))
            all_outputs.append(os.path.join(sec_workdir, model_name))
        else:
            all_outputs.append(workdir)
            all_outputs.append(ref_workdir)
            all_outputs.append(sec_workdir)

        # capture if tropo and separate distinct wet and hydro layers
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
                    proj, outputFormat, model_name)

            # if specified, compute total delay
            if tropo_total:
                model_dir = os.path.abspath(workdir)

                if not is_nisar_file:
                    model_dir = os.path.join(model_dir, model_name)
                    # compute reference diff
                    ref_diff = ref_outname
                    sec_diff = sec_outname
                    outname_diff = os.path.join(model_dir, 'dates',
                                                os.path.basename(ref_diff))
                    if not os.path.exists(outname_diff):
                        generate_diff(
                            ref_diff, sec_diff, outname_diff, key, sec_key,
                            tropo_total, hgt_field, proj, outputFormat,
                            dem=dem)
                    # compute secondary diff
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

                    # compute total diff
                    ref_diff = os.path.join(ref_workdir, model_name, ifg)
                    sec_diff = os.path.join(sec_workdir, model_name, ifg)

                else:
                    # compute total diff
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

    # delete temporary files if layers not requested
    prod_ver_list = i[1]
    for i in all_workdirs:
        key_name = os.path.basename(i)
        if os.path.exists(i):
            if key_name not in layers or len(os.listdir(i)) == 0:
                shutil.rmtree(i)

    # interpolate and intersect epochs for user requested layers
    all_outputs = list(set(all_outputs))
    for i in enumerate(all_outputs):
        if os.path.exists(i[1]):
            record_epochs = []
            record_epochs.extend(
                glob.glob(os.path.join(i[1], '*[0-9].vrt')))
            record_epochs.extend(
                glob.glob(os.path.join(i[1], 'dates/*[0-9].vrt')))

            for j in enumerate(record_epochs):
                # dedup check for interpolating only new files
                ds_count = osgeo.gdal.Open(j[1], osgeo.gdal.GA_ReadOnly)
                band_count = ds_count.RasterCount
                ds_count = None
                if band_count == 1:
                    # Track consistency of dimensions
                    if j[0] == 0:
                        ref_wid, ref_hgt, ref_geotrans, _, _ = \
                            ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                        ref_arr = [ref_wid, ref_hgt, ref_geotrans, j[1][:-4]]

                    continue

                # Interpolate/intersect with DEM before cropping
                finalize_metadata(
                    j[1][:-4], bounds, arrres, dem_bounds, prods_TOTbbox,
                    dem, lat, lon, hgt_field, prod_ver_list, is_nisar_file,
                    outputFormat, verbose=verbose)

                # Apply mask (if specified)
                if mask is not None:
                    # Load mask
                    ds_vrt_read = osgeo.gdal.Open(
                        j[1][:-4] + '.vrt', osgeo.gdal.GA_ReadOnly
                    )
                    vrt_arr = ds_vrt_read.ReadAsArray()
                    ds_vrt_read = None
                    mask_arr = mask.ReadAsArray() * vrt_arr

                    # Initiate file update with mask
                    update_file = osgeo.gdal.Open(
                        j[1][:-4], osgeo.gdal.GA_Update
                    )
                    update_file.GetRasterBand(1).WriteArray(mask_arr)

                    # Clear variables
                    update_file = None
                    mask_arr = None

                # Track consistency of dimensions
                if j[0] == 0:
                    ref_wid, ref_hgt, ref_geotrans, _, _ = \
                        ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                    ref_arr = [ref_wid, ref_hgt, ref_geotrans, j[1][:-4]]

                else:
                    prod_wid, prod_hgt, prod_geotrans, _, _ = \
                        ARIAtools.util.vrt.get_basic_attrs(j[1][:-4])
                    prod_arr = [prod_wid, prod_hgt, prod_geotrans, j[1][:-4]]
                    ARIAtools.util.vrt.dim_check(ref_arr, prod_arr)
                prev_outname = j[1][:-4]

    # pass final list of outputs
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
    # Re-authenticate GDAL for the isolated worker process
    ARIAtools.product._configure_gdal_virtual_access()

    # Initialize warp dict
    gdal_warp_kwargs = {
        'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
        'outputBounds': bounds, 'xRes': arrres[0], 'yRes': arrres[1],
        'targetAlignedPixels': True, 'multithread': True, 'dstSRS': proj}
    warp_options = osgeo.gdal.WarpOptions(
        **gdal_warp_kwargs
    )

    # Flip sign for baselines if NISAR
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

    # Load product dict file
    with open(full_product_dict_file, 'r') as ifp:
        full_product_dict = json.load(ifp)

    product_dict = [[j[layers[ilayer]] for j in full_product_dict],
                    [j["pair_name"] for j in full_product_dict]]

    ifg_tag = product_dict[1][ii][0]
    outname = os.path.abspath(os.path.join(workdir, ifg_tag))

    if update_mode != 'crop_only' \
            and os.path.exists(outname) \
            and os.path.exists(outname + '.vrt'):
        LOGGER.debug('Skipping %s - %s', ifg_tag,
                     {os.path.dirname(outname).split('/')[-1]})

    elif update_mode == 'crop_only' \
            and os.path.exists(outname) \
            and os.path.exists(outname + '.vrt'):
        lyrname = os.path.dirname(outname).split('/')[-1]
        crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs)
        # make sure to update conn comp file(s)
        if os.path.dirname(outname).split('/')[-1] == 'unwrappedPhase':
            lyrname = 'connectedComponents'
            # Split the path into components
            path_parts = outname.split('/')

            # Replace "unwrappedPhase" only at the second-to-last index
            if path_parts[-2] == 'unwrappedPhase':
                path_parts[-2] = lyrname

            # Rejoin the path
            outname = '/'.join(path_parts)
            crop_only_manager(outname, lyrname, ifg_tag, gdal_warp_kwargs)

    else:
        LOGGER.debug('Extracting %s - %s', ifg_tag,
                     {os.path.dirname(outname).split('/')[-1]})

        # Extract/crop metadata layers
        if (any(':/science/grids/imagingGeometry' in s for s in product) or
            any(':/science/LSAR/GUNW/metadata/radarGrid' in s
                for s in product)):
            # make VRT pointing to metadata layers in standard product
            hgt_field, outname = prep_metadatalayers(
                outname, product, dem_expanded, layer, layers,
                is_nisar_file, proj, sign_multiplier=sign_multiplier)

            # Interpolate/intersect with DEM before cropping
            finalize_metadata(
                outname, bounds, arrres, dem_bounds, prods_TOTbbox,
                dem_expanded, lat, lon, hgt_field, product, is_nisar_file,
                outputFormatPhys, verbose=verbose)

        # Extract/crop full res layers, except for "unw" and "conn_comp"
        # which requires advanced stitching
        elif layer != 'unwrappedPhase' and layer != 'connectedComponents':

            if is_nisar_file:

                if layer == 'amplitude':
                    amp_ds_list = []
                    # Ensure product is a list
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
                        
                        # 1. Standardize any weird NaNs back to 0 
                        # so GDAL's C++ engine can safely identify the transparent edge padding
                        amp_arr[np.isnan(amp_arr)] = 0
                        
                        amp_band.WriteArray(amp_arr)
                        amp_band.SetNoDataValue(0)
                        
                        amp_ds_list.append(ds_amp)
                        ds_in = None

                    amp_kwargs = gdal_warp_kwargs.copy()
                    amp_kwargs['format'] = outputFormatPhys
                    if ('dstSRS' in amp_kwargs and
                            isinstance(amp_kwargs['dstSRS'], int)):
                        amp_kwargs['dstSRS'] = f"EPSG:{amp_kwargs['dstSRS']}"

                    # 2. srcNodata=0 forces the overlapping blank edges to be completely transparent.
                    # 3. dstNodata=np.nan converts the final stitched background safely back to NaN!
                    amp_warp_opts = osgeo.gdal.WarpOptions(
                        outputType=osgeo.gdal.GDT_Float32,
                        srcNodata=0,
                        dstNodata=np.nan,
                        **amp_kwargs
                    )

                    # Warp directly from the MEM datasets
                    ds_amp_warp = osgeo.gdal.Warp(
                        outname, amp_ds_list, options=amp_warp_opts
                    )
                    
                    # Cleanup
                    ds_amp_warp = None
                    amp_ds_list = None

                else:
                    # Standard NISAR layer options
                    # 1. If multiple frames are passed, build a VRT mosaic
                    if isinstance(product, list) and len(product) > 1:
                        tmp_mosaic = str(outname) + "_uncropped.vrt"
                        
                        # Reproject heterogeneous UTM zones safely via VRTs
                        tmp_vrts = []
                        for idx, p in enumerate(product):
                            t_vrt = f"{outname}_{idx}_tmp.vrt"
                            osgeo.gdal.Warp(
                                t_vrt, p, format="VRT", dstSRS=proj
                            )
                            tmp_vrts.append(t_vrt)
                            
                        osgeo.gdal.BuildVRT(tmp_mosaic, tmp_vrts)
                        warp_inputs = tmp_mosaic
                    else:
                        warp_inputs = (
                            product[0] if isinstance(product, list)
                            else product
                        )

                    # 2. Safely warp the single mosaic/file
                    if outputFormat == 'VRT':
                        ds = osgeo.gdal.Warp(
                            outname + '.vrt', warp_inputs, options=warp_options
                        )
                        ds = None
                    else:
                        ds = osgeo.gdal.Warp(
                            outname, warp_inputs, options=warp_options
                        )
                        ds = None

            else:
                # Legacy handling
                with osgeo.gdal.config_options(
                    {"GDAL_NUM_THREADS": num_threads}
                ):

                    if outputFormat == 'VRT':
                        ds_vrt = osgeo.gdal.BuildVRT(
                            outname + "_uncropped.vrt", product
                        )
                        ds_vrt = None
                        ds = osgeo.gdal.Warp(
                            outname + '.vrt',
                            outname + '_uncropped.vrt',
                            options=warp_options
                        )
                        ds = None
                    else:
                        ds_vrt = osgeo.gdal.BuildVRT(outname + '.vrt', product)
                        ds_vrt = None
                        ds = osgeo.gdal.Warp(
                            outname,
                            outname + '.vrt',
                            options=warp_options
                        )
                        ds = None
                        ds_trans = osgeo.gdal.Translate(
                            outname + '.vrt',
                            outname,
                            options=osgeo.gdal.TranslateOptions(
                                format="VRT"
                            )
                        )
                        ds_trans = None

            # Create VRT pointing to physical file if a physical file was written
            if os.path.exists(outname):
                ds_trans = osgeo.gdal.Translate(
                    outname + '.vrt', outname, format="VRT"
                )
                ds_trans = None

            # VRT formats require the source mosaic to remain on disk.
            # Only delete the uncropped mosaic if data was physically extracted.
            if outputFormat != 'VRT':
                tmp_mosaic = str(outname) + "_uncropped.vrt"
                if os.path.exists(tmp_mosaic):
                    os.remove(tmp_mosaic)
                
                # Clean up intermediate heterogeneous projection VRTs!
                if isinstance(product, list) and len(product) > 1:
                    for idx in range(len(product)):
                        t_vrt = f"{outname}_{idx}_tmp.vrt"
                        if os.path.exists(t_vrt):
                            os.remove(t_vrt)

        # Extract/crop phs and conn_comp layers
        else:
            # get connected component input files
            conn_files = full_product_dict[ii]['connectedComponents']
            prod_bbox_files = full_product_dict[ii][
                'productBoundingBoxFrames']
            outFileConnComp = os.path.join(
                outDir, 'connectedComponents', ifg_tag)

            # Check if phs phase and conn_comp files are already generated
            outFilePhs = os.path.join(outDir, 'unwrappedPhase', ifg_tag)
            # if (not os.path.exists(outFilePhs) or
            #         not os.path.exists(outFileConnComp)):

            phs_files = full_product_dict[ii]['unwrappedPhase']

            # stitching
            ARIAtools.util.seq_stitch.product_stitch_sequential(
                phs_files, conn_files, arrres=arrres, epsg=proj,
                bounds=bounds, clip_json=prods_TOTbbox, output_unw=outFilePhs,
                output_conn=outFileConnComp,
                output_format=outputFormatPhys,
                is_nisar_file=is_nisar_file,
                range_correction=range_correction, save_fig=False,
                overwrite=True)

            # If necessary, resample phs/conn_comp file
            if multilooking is not None:
                ARIAtools.util.vrt.resampleRaster(
                    outFilePhs, multilooking, bounds, prods_TOTbbox,
                    rankedResampling, outputFormat=outputFormatPhys,
                    num_threads=num_threads)

            # Apply mask (if specified)
            if mask is not None:
                for j in [outFileConnComp, outFilePhs]:
                    # Load mask
                    ds_vrt_read = osgeo.gdal.Open(
                        j + '.vrt', osgeo.gdal.GA_ReadOnly
                    )
                    vrt_arr = ds_vrt_read.ReadAsArray()
                    ds_vrt_read = None
                    mask_arr = mask.ReadAsArray() * vrt_arr

                    # Initiate file update with mask
                    update_file = osgeo.gdal.Open(
                        j, osgeo.gdal.GA_Update
                    )
                    update_file.GetRasterBand(1).WriteArray(mask_arr)

                    # Clear variables
                    update_file = None
                    mask_arr = None

        if layer != 'unwrappedPhase' and layer != 'connectedComponents':

            # If necessary, resample raster
            if multilooking is not None:
                ARIAtools.util.vrt.resampleRaster(
                    outname, multilooking, bounds, prods_TOTbbox,
                    rankedResampling, outputFormat=outputFormatPhys,
                    num_threads=num_threads)

            # Apply mask (if specified)
            if mask is not None:
                # Load mask
                ds_vrt_read = osgeo.gdal.Open(
                    outname + '.vrt', osgeo.gdal.GA_ReadOnly
                )
                vrt_arr = ds_vrt_read.ReadAsArray()
                ds_vrt_read = None
                mask_arr = mask.ReadAsArray() * vrt_arr

                # Initiate file update with mask
                update_file = osgeo.gdal.Open(
                    outname, osgeo.gdal.GA_Update
                )
                update_file.GetRasterBand(1).WriteArray(mask_arr)

                # Clear variables
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
        return  # only bbox

    # initiate tracker of output dimensions
    ref_wid = None
    ref_hgt = None
    ref_geotrans = None
    ref_arr = None

    mask = None if maskfile is None else osgeo.gdal.Open(maskfile)
    dem = None if demfile is None else osgeo.gdal.Open(demfile)
    dem_expanded = (
        None if demfile_expanded is None
        else osgeo.gdal.Open(demfile_expanded))

    # create dictionary of all inputs needed for correction lyr extraction
    # Get the authority code (EPSG code)
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

    # track if product stack is NISAR GUNW or not
    range_correction = True
    track_fileext = full_product_dict[0]['unwrappedPhase'][0]
    if is_nisar_file:
        range_correction = False
        model_names = ['']
    else:
        model_names = [f'_{i}' for i in model_names]
    lyr_input_dict['is_nisar_file'] = is_nisar_file

    # get bounds
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds
    lyr_input_dict['bounds'] = bounds
    lyr_input_dict['arrres'] = arrres
    if dem_expanded is not None:
        dem_gt = dem_expanded.GetGeoTransform()
        dem_bounds = [
            dem_gt[0], dem_gt[3] + (dem_gt[-1] * dem_expanded.RasterYSize),
            dem_gt[0] + (dem_gt[1] * dem_expanded.RasterXSize), dem_gt[3]]
        lyr_input_dict['dem_bounds'] = dem_bounds

    # Mask specified, so file must be physically extracted,
    # cannot proceed with VRT format. Defaulting to ENVI format.
    if (outputFormat == 'VRT' and mask is not None) or \
            (outputFormat == 'VRT' and multilooking is not None):
        outputFormat = 'ENVI'

    # Set output format layers that must always be physically extracted
    outputFormatPhys = 'ENVI'
    if outputFormat != 'VRT':
        outputFormatPhys = outputFormat
    lyr_input_dict['outputFormat'] = outputFormatPhys

    # Recall update mode and conduct final checks for extraction
    if runlog is None:
        update_mode = 'full_extract'
    else:
        log_data = runlog.load()
        update_mode = log_data['update_mode']
        if 'update_mode' in log_data.keys():
            update_mode = log_data['update_mode']

        # Check water mask
        prev_maskfile = log_data['maskfilename'] if 'maskfilename' \
            in log_data.keys() else None

        if maskfile != prev_maskfile and prev_maskfile is not None:
            update_mode = 'full_extract'
            LOGGER.warning(
                'Mask file has changed. Setting update mode to full_extract.')

        runlog.update('maskfilename', maskfile)

        # Check DEM
        prev_demfile = log_data['demfile'] if 'demfile' \
            in log_data.keys() else None

        if demfile != prev_demfile and prev_demfile is not None:
            update_mode = 'full_extract'
            LOGGER.warning(
                'DEM file has changed. Setting update mode to full_extract.')

        runlog.update('demfile', demfile)

        # Update final mode
        runlog.update('update_mode', update_mode)

    # track extracted layers
    extracted_files = []

    # Initialize warp dict
    gdal_warp_kwargs = {
        'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
        'outputBounds': bounds, 'xRes': arrres[0], 'yRes': arrres[1],
        'targetAlignedPixels': True, 'multithread': True, 'dstSRS': epsg_code}

    # track if files need to be updated
    lyr_input_dict['update_mode'] = update_mode
    lyr_input_dict['gdal_warp_kwargs'] = gdal_warp_kwargs

    # If specified, extract tropo layers
    tropo_lyrs = ['troposphereWet', 'troposphereHydrostatic']
    user_lyrs = list(set(layers).intersection(tropo_lyrs))
    if tropo_total or user_lyrs != []:
        # set input keys
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

        # loop through valid models
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

            # get unique layer names from path
            map_lyrs = [
                product_dict[0][0][0].split('/')[-1],
                product_dict_dry[0][0].split('/')[-1]
            ]
            lyr_input_dict['map_lyrs'] = map_lyrs

            # set iterative keys
            lyr_input_dict['product_dict'] = product_dict

            # extract layers
            extracted_files.extend(handle_epoch_layers(**lyr_input_dict))

            # remove leading underscore from model name to get subdir name
            tag = i.split('_')[-1]
            # track valid files
            prev_outname = os.path.abspath(
                os.path.join(workdir,
                             tag,
                             product_dict[1][0][0])
            )
            if os.path.exists(prev_outname + '.vrt'):
                prev_outname_check = copy.deepcopy(prev_outname)

        # track consistency of dimensions
        if 'prev_outname_check' in locals():
            ref_wid, ref_hgt, ref_geotrans, _, _ = \
                ARIAtools.util.vrt.get_basic_attrs(prev_outname_check + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans, prev_outname]

    # If specified, extract solid earth tides
    tropo_lyrs = list(set(tropo_lyrs))
    ext_corr_lyrs = tropo_lyrs + ['solidEarthTide', 'troposphereTotal']
    if 'solidEarthTide' in layers:
        lyr_prefix = '/science/grids/corrections/external/tides/solidEarth/'
        key = 'solidEarthTide'
        ref_key = key
        sec_key = key
        product_dict = [
            [j[key] for j in full_product_dict if key in j.keys()],
            [j["pair_name"] for j in full_product_dict if key in j.keys()]]

        # get unique layer names from path
        map_lyrs = [product_dict[0][0][0].split('/')[-1]]
        lyr_input_dict['map_lyrs'] = map_lyrs

        workdir = os.path.join(outDir, key)
        prev_outname = copy.deepcopy(workdir)

        # set input keys
        lyr_input_dict['product_dict'] = product_dict
        lyr_input_dict['lyr_path'] = lyr_prefix
        lyr_input_dict['user_lyrs'] = ['solidEarthTide']
        lyr_input_dict['key'] = key
        lyr_input_dict['sec_key'] = sec_key
        lyr_input_dict['ref_key'] = ref_key
        lyr_input_dict['tropo_total'] = False
        lyr_input_dict['workdir'] = workdir

        # extract layers
        extracted_files.extend(handle_epoch_layers(**lyr_input_dict))

        # Track consistency of dimensions
        prev_outname = os.path.abspath(os.path.join(workdir,
                                       product_dict[1][0][0]))
        ref_wid, ref_hgt, ref_geotrans, \
            _, _ = ARIAtools.util.vrt.get_basic_attrs(prev_outname + '.vrt')
        ref_arr = [ref_wid, ref_hgt, ref_geotrans,
                   prev_outname]

    # If specified, extract ionosphere long wavelength
    ext_corr_lyrs += ['ionosphere']
    if 'ionosphere' in layers:
        lyr_prefix = '/science/grids/corrections/derived/ionosphere/ionosphere'
        key = 'ionosphere'
        product_dict = \
            [[j[key] for j in full_product_dict if key in j.keys()],
             [j["pair_name"] for j in full_product_dict if key in j.keys()]]

        workdir = os.path.join(outDir, key)
        prev_outname = copy.deepcopy(workdir)

        # Set output res
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

            # if file exists and needs to be cropped, avoid iono routine
            if os.path.exists(outname) and update_mode == 'crop_only':
                crop_only_manager(outname, 'ionosphere',
                    product_dict[1][i][0], gdal_warp_kwargs)

            # only extract if file does not exist
            if not os.path.exists(outname):
                ARIAtools.util.ionosphere.export_ionosphere(**lyr_input_dict)

            # track output
            extracted_files.append(outname)

            # track valid files
            if os.path.exists(outname + '.vrt'):
                prev_outname_check = copy.deepcopy(outname)
                
            prog_bar.update(i + 1)
        prog_bar.close()

        # track consistency of dimensions
        if 'prev_outname_check' in locals():
            ref_wid, ref_hgt, ref_geotrans, _, _ = \
                ARIAtools.util.vrt.get_basic_attrs(prev_outname_check + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans, prev_outname]

    # Update runlog if provided
    if runlog is not None:
        runlog.update('extracted_files', extracted_files)

    # Loop through other user expected layers
    layers = [i for i in layers if i not in ext_corr_lyrs]

    full_product_dict_file = os.path.join(outDir, 'full_product_dict.json')
    with open(full_product_dict_file, 'w') as ofp:
        json.dump(full_product_dict, ofp)

    for ilayer, layer in enumerate(layers):

        product_dict = [[j[layer] for j in full_product_dict],
                        [j["pair_name"] for j in full_product_dict]]

        # If specified workdir doesn't exist, create it
        workdir = os.path.join(outDir, layer)
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        # Log height subsetting info once per geometry layer
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
                            ARIAtools.util.interp._compute_dem_range(
                                dem_expanded))
                        idx = (ARIAtools.util.interp
                               ._get_height_subset_indices(
                                   hts, d_min, d_max, pad=1))
                        if len(idx) < len(hts):
                            LOGGER.info(
                                'Height subsetting %s: %d \u2192 %d levels '
                                '(DEM range: %.0f to %.0f m)',
                                layer, len(hts), len(idx), d_min, d_max)
            except Exception:
                pass

        mp_args = []
        # Iterate through all IFGs
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
            
            # Initialize the custom ARIA progress bar
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

                # Set up a thread-safe counter for Dask
                lock = threading.Lock()
                completed = 0

                def update_progress(result):
                    nonlocal completed
                    with lock:
                        completed += 1
                        prog_bar.update(completed)
                    return result

                # Create jobs wrapped with our thread-safe progress updater
                jobs = []
                for arg in mp_args:
                    job = dask.delayed(
                        lambda x: update_progress(export_product_worker(*x))
                    )(arg)
                    jobs.append(job)

                # Compute all jobs
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

            # Run the export worker jobs with GNU parallel in the background
            proc = subprocess.Popen((
                'find %s/export_workers -name "export_product_args_*.json" | '
                'parallel -j %d export_product.py {}') % (
                    outDir, int(num_threads)), shell=True)

            # Poll the directory for completed JSON files to update progress
            while proc.poll() is None:
                num_done = len(glob.glob(
                    os.path.join(export_workers_temp_dir, 'outputs_*.json')
                ))
                prog_bar.update(num_done)
                time.sleep(1.0)

            # Catch the final update immediately after the process finishes
            num_done = len(glob.glob(
                os.path.join(export_workers_temp_dir, 'outputs_*.json')
            ))
            prog_bar.update(num_done)
            prog_bar.close()

            # load in output files and verify dimensions
            output_files = glob.glob(os.path.join(
                export_workers_temp_dir, 'outputs_*.json'))

            if len(output_files) > 0:
                # Remove hardcoded 0_0 so it grabs the correct layer file
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

    # Update runlog if provided
    if runlog is not None:
        runlog.update('extracted_files', extracted_files)

    # delete directory for quality control plots if empty
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

    # Check if this layer needs height-based DEM intersection
    NOHGT_LYRS = ['ionosphere']
    metadatalyr_name = outname.split('/')[-2]
    needs_height_interp = metadatalyr_name not in NOHGT_LYRS

    # --- Height-based band subsetting optimisation ---
    # Only load the vertical layers that span the DEM elevation range
    # instead of the entire 3D cube.  This reduces I/O, memory, and
    # interpolation cost.
    tmp_name = outname + '.vrt'
    warp_src = tmp_name
    heightsMeta = None
    subset_vrt = None

    if needs_height_interp:
        # Get height levels from VRT metadata (no data loading)
        heightsMeta_str = ARIAtools.util.vrt.get_hgt_meta(
            tmp_name, hgt_field)
        heightsMeta = np.array(
            heightsMeta_str[1:-1].split(','), dtype='float32')

        # Get DEM elevation range (nodata-aware)
        dem_min, dem_max = ARIAtools.util.interp._compute_dem_range(dem)

        # Compute which height bands are needed
        # Set ARIA_DISABLE_HEIGHT_SUBSET=1 to bypass for benchmarking
        if not os.environ.get('ARIA_DISABLE_HEIGHT_SUBSET'):
            band_indices = ARIAtools.util.interp._get_height_subset_indices(
                heightsMeta, dem_min, dem_max, pad=1)
        else:
            band_indices = np.arange(len(heightsMeta))

        if len(band_indices) < len(heightsMeta):
            LOGGER.debug(
                'Subsetting 3D cube from %d to %d height levels '
                '(DEM range: %.1f to %.1f)',
                len(heightsMeta), len(band_indices), dem_min, dem_max)

            # Create band-selected VRT to avoid loading unnecessary bands
            band_list = [int(i + 1) for i in band_indices]  # GDAL 1-based
            subset_vrt = outname + '_hsubset.vrt'
            translate_opts = osgeo.gdal.TranslateOptions(
                format='VRT', bandList=band_list)
            ds_sub = osgeo.gdal.Translate(
                subset_vrt, tmp_name, options=translate_opts)
            ds_sub = None
            warp_src = subset_vrt
            heightsMeta = heightsMeta[band_indices]

    # load layered metadata array (possibly band-subsetted)
    # Spatially crop to DEM extent, padded by 2 native grid cells of the
    # source 3D cube to avoid interpolation artefacts at the edges.
    # RegularGridInterpolator (linear) needs ≥1 cell; we use 2 for safety.
    ds_src = osgeo.gdal.Open(warp_src, osgeo.gdal.GA_ReadOnly)
    src_gt = ds_src.GetGeoTransform()
    src_xres = abs(src_gt[1])
    src_yres = abs(src_gt[5])
    ds_src = None
    pad_cells = 2
    padded_bounds = [
        dem_bounds[0] - pad_cells * src_xres,   # xmin
        dem_bounds[1] - pad_cells * src_yres,   # ymin
        dem_bounds[2] + pad_cells * src_xres,   # xmax
        dem_bounds[3] + pad_cells * src_yres,   # ymax
    ]
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(
            format="MEM", outputBounds=padded_bounds)
        ds_warp = osgeo.gdal.Warp('', warp_src, options=warp_options)
        data_array_nodata = ds_warp.GetRasterBand(1).GetNoDataValue()
        data_array = ds_warp.ReadAsArray().astype('float32')
        gt_mem = ds_warp.GetGeoTransform()
        x_size = ds_warp.RasterXSize
        y_size = ds_warp.RasterYSize
        ds_warp = None  # CLOSE

    # Clean up subset VRT if created
    if subset_vrt is not None and os.path.exists(subset_vrt):
        os.remove(subset_vrt)

    # Ensure data_array is 3-D even when only one band was loaded
    if data_array.ndim == 2:
        data_array = data_array[np.newaxis, ...]

    # get minimum version
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
        # create directory for quality control plots
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
        tmp_name = outname + '_temp'

        # heightsMeta already extracted above

        latitudeMeta = np.linspace(
            gt_mem[3], gt_mem[3] + (gt_mem[5] * (y_size - 1)),
            y_size, dtype='float32')

        longitudeMeta = np.linspace(
            gt_mem[0], gt_mem[0] + (gt_mem[1] * (x_size - 1)),
            x_size, dtype='float32')

        # --- SAFE RIOXARRAY BLOCK ---
        # Using 'with' ensures the handle to the DEM file is dropped
        # immediately after reading.
        with rioxarray.open_rasterio(
            dem.GetDescription(), band_as_variable=True, masked=True
        ) as rds:
            da_dem = rds['band_1']
            
            # interpolate the DEM to the GUNW lat/lon
            nodata = dem.GetRasterBand(1).GetNoDataValue()
            
            # Force close check: Ensure we don't hold the file open during compute
            da_dem1 = da_dem.interp(
                x=lon[0, :], y=lat[:, 0]
            ).fillna(nodata)

        # hack to get an stack of coordinates for the interpolator
        pnts = transformPoints(
            lat, lon, da_dem1.data, 'EPSG:4326', 'EPSG:4326')

        # set up the interpolator with the (subsetted) GUNW cube
        interper = scipy.interpolate.RegularGridInterpolator(
            (latitudeMeta, longitudeMeta, heightsMeta),
            data_array.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        # interpolate cube to DEM points
        out_interpolated = interper(pnts.transpose(2, 1, 0))

        # Save file (Using GDAL to ensure clean write)
        ARIAtools.util.vrt.renderVRT(
            tmp_name, out_interpolated, geotrans=dem.GetGeoTransform(),
            drivername=outputFormat,
            gdal_fmt='float32',
            proj=dem.GetProjection(), nodata=nodata)
        out_interpolated = None

    # Since metadata layer extends at least one grid node
    # outside of the expected track bounds,
    # it must be cut to conform with these bounds.
    # Crop to track extents
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
            'outputBounds': dem_bounds, 'dstNodata': data_array_nodata,
            'xRes': dem_arrres[0], 'yRes': dem_arrres[1],
            'targetAlignedPixels': True, 'multithread': True}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp(
            tmp_name + '_temp', tmp_name, options=warp_options
        )
        ds = None

    # Adjust shape
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
            'outputBounds': bbox_bounds, 'dstNodata': data_array_nodata,
            'xRes': arrres[0], 'yRes': arrres[1], 'targetAlignedPixels': True,
            'multithread': True}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp(
            outname, tmp_name + '_temp', options=warp_options
        )
        ds = None

    # remove temp files
    for i in glob.glob(outname + '*_temp*'):
        os.remove(i)

    # Update VRT
    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    vrt_ds = osgeo.gdal.Translate(
        outname + '.vrt', outname, options=translate_options)
    vrt_ds = None

    data_array = None

    # --- ADD THIS: Safely destroy the incoming DEM dataset object ---
    # This prevents anonymous gdal.Open() calls from the parent wrapper
    # from surviving past the end of this function and crashing the GC.
    dem = None
    lat = None
    lon = None

    return


def transformPoints(lats: np.ndarray, lons: np.ndarray, hgts: np.ndarray,
                    old_proj: pyproj.CRS, new_proj: pyproj.CRS) -> np.ndarray:
    '''
    Transform lat/lon/hgt data to an array of points in a new
    projection
    Args:
        lats: ndarray - WGS-84 latitude (EPSG: 4326)
        lons: ndarray - ditto for longitude
        hgts: ndarray - Ellipsoidal height in meters
        old_proj: pyproj.CRS - original projection of the points
        new_proj: pyproj.CRS - new projection in which to return the points
    Returns:
        ndarray: array of query points in weather model coordinate system (YX)
    '''
    transformer = pyproj.Transformer.from_crs(old_proj, new_proj)

    # Flags for flipping inputs or outputs
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
