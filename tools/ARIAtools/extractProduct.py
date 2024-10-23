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
import tqdm

import dask
import rioxarray
import rasterio
import osgeo
import pyproj
import numpy as np
import scipy.interpolate
import shapely.geometry

import ARIAtools.util.ionosphere
import ARIAtools.util.vrt
import ARIAtools.util.shp
import ARIAtools.util.misc
import ARIAtools.util.seq_stitch

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
        return self.data_array


def merged_productbbox(
        metadata_dict, product_dict, workdir='./', bbox_file=None,
        croptounion=False, num_threads='2', minimumOverlap=0.0081,
        verbose=None):
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
            raise Exception(f'User bound box {bbox_file} has an area of only '
                            f'{overlap_area}km\u00b2, below specified '
                            f'minimum threshold area '
                            f'{minimumOverlap}km\u00b2')

    # Extract/merge productBoundingBox layers
    for scene in product_dict:

        # Get pair name, expected in dictionary
        pair_name = scene["pair_name"][0]
        outname = os.path.join(workdir, pair_name + '.json')

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
    prods_TOTbbox_metadatalyr = os.path.join(
        workdir, 'productBoundingBox_croptounion_formetadatalyr.json')
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
        LOGGER.info(("%d out of %d interferograms rejected for not meeting "
                     "specified spatial thresholds"),
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

    # Warp the first scene with the output-bounds defined above
    # ensure output-bounds are an integer multiple of interferometric grid
    # and adjust if necessary
    OG_bounds = list(
        ARIAtools.util.shp.open_shp(bbox_file).bounds)
    gdal_warp_kwargs = {
        'format': 'MEM', 'multithread': True, 'dstSRS': f'EPSG:{lyr_proj}'}
    vrt = osgeo.gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', vrt, options=warp_options)
        arrres = [abs(ds.GetGeoTransform()[1]),
                  abs(ds.GetGeoTransform()[-1])]
        ds = None

    # warp again with fixed transform and bounds
    gdal_warp_kwargs['outputBounds'] = OG_bounds
    gdal_warp_kwargs['xRes'] = arrres[0]
    gdal_warp_kwargs['yRes'] = arrres[1]
    gdal_warp_kwargs['targetAlignedPixels'] = True
    vrt = osgeo.gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0])
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        ds = osgeo.gdal.Warp('', vrt, options=warp_options)

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

    return (metadata_dict, product_dict, bbox_file, prods_TOTbbox,
            prods_TOTbbox_metadatalyr, arrres, proj, is_nisar_file)


def create_raster_from_gunw(fname, data_lis, proj, driver, hgt_field=None):
    """Wrapper to create raster and apply projection"""
    # open original raster
    osgeo.gdal.BuildVRT(fname + '_temp.vrt', data_lis)
    da = rioxarray.open_rasterio(fname + '_temp.vrt', masked=True)

    # Reproject the raster to the desired projection
    reproj_da = da.rio.reproject(f'EPSG:{proj}',
                                 resampling=rasterio.enums.Resampling.nearest,
                                 nodata=0)
    reproj_da.rio.to_raster(fname, driver=driver, crs=f'EPSG:{proj}')
    os.remove(fname + '_temp.vrt')
    da.close()
    reproj_da.close()
    buildvrt_options = osgeo.gdal.BuildVRTOptions(outputSRS=f'EPSG:{proj}')
    osgeo.gdal.BuildVRT(fname + '.vrt', fname, options=buildvrt_options)

    if hgt_field is not None:
        # write height layers
        hgt_meta = osgeo.gdal.Open(data_lis[0]).GetMetadataItem(hgt_field)
        osgeo.gdal.Open(
            fname + '.vrt').SetMetadataItem(hgt_field, hgt_meta)

    return


def prep_metadatalayers(
        outname, metadata_arr, dem, layer, layers, is_nisar_file=False,
        proj='4326', driver='ENVI', model_name=None):
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
        osgeo.gdal.BuildVRT(outname + '.vrt', metadata_arr)
        return [0], None, outname

    # capture model if tropo product
    if 'tropo' in layer:
        if not is_nisar_file:
            out_dir = os.path.join(out_dir, model_name)
        outname = os.path.join(out_dir, ifg)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    # Get height values
    zdim = osgeo.gdal.Open(metadata_arr[0]).GetMetadataItem(
        'NETCDF_DIM_EXTRA')[1:-1]
    hgt_field = f'NETCDF_DIM_{zdim}_VALUES'

    # Check if height layers are consistent
    if not os.path.exists(outname + '.vrt') and \
        not len(set([osgeo.gdal.Open(i).GetMetadataItem(hgt_field)
                     for i in metadata_arr])) == 1:
        raise Exception(
            'Inconsistent heights for metadata layer(s) ', metadata_arr,
            ' corresponding heights: ', [osgeo.gdal.Open(i).GetMetadataItem(
                hgt_field) for i in metadata_arr])

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
            create_raster_from_gunw(i[0], i[1], proj, driver, hgt_field)

        if not is_nisar_file:
            # compute differential
            generate_diff(
                ref_outname, sec_outname, outname, layer, layer, False,
                hgt_field, proj, driver)

        # write raster to file if it does not exist
        if layer in layers:
            for i in [ref_outname, sec_outname]:
                if not os.path.exists(i):
                    create_raster_from_gunw(i, [i], proj, driver, hgt_field)

    else:
        if not os.path.exists(outname + '.vrt'):
            osgeo.gdal.BuildVRT(outname + '.vrt', metadata_arr)

            # write height layers
            osgeo.gdal.Open(outname + '.vrt').SetMetadataItem(
                hgt_field, osgeo.gdal.Open(
                    metadata_arr[0]).GetMetadataItem(hgt_field))

    return hgt_field, ref_outname


def generate_diff(ref_outname, sec_outname, outname, key, OG_key, tropo_total,
                  hgt_field, proj, driver):
    """ Compute differential from reference and secondary scenes """

    # if specified workdir doesn't exist, create it
    output_dir = os.path.dirname(outname)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # open intermediate files
    with rioxarray.open_rasterio(sec_outname + '.vrt', masked=True) as da_sec:
        arr_sec = da_sec.data
    with rioxarray.open_rasterio(ref_outname + '.vrt', masked=True) as da_ref:
        arr_ref = da_ref.data

    # make the total arr
    # if computing total delay, must add dry and wet
    if tropo_total:
        arr_total = arr_sec + arr_ref
    else:
        arr_total = arr_sec - arr_ref
    da_total = da_sec.copy()
    da_total.data = arr_total

    # update attributes
    da_total.name = key
    og_da_attrs = da_total.attrs
    da_attrs = {}
    for key in og_da_attrs:
        new_key = key.replace(OG_key, key)
        new_val = og_da_attrs[key]
        if isinstance(new_val, str):
            new_val = new_val.replace(OG_key, key)
        da_attrs[new_key] = new_val
    da_total = da_total.assign_attrs(da_attrs)

    # write initial array to file
    da_total.rio.to_raster(outname, driver=driver, crs=f'EPSG:{proj}')
    buildvrt_options = osgeo.gdal.BuildVRTOptions(outputSRS=f'EPSG:{proj}')
    ds = osgeo.gdal.BuildVRT(
        f'{outname}.vrt', outname, options=buildvrt_options)
    # fix if numpy array not set properly
    if not isinstance(da_attrs[hgt_field], np.ndarray):
        da_attrs[hgt_field] = np.array(da_attrs[hgt_field])
    da_attrs[hgt_field] = da_attrs[hgt_field].tolist()
    ds.SetMetadata(da_attrs)
    ds = None

    return


def extract_bperp_dict(products, num_threads):
    """Extracts bPerpendicular mean over frames for each product in products"""
    def read_and_average_bperp(frame):
        """Helper function for dask multiprocessing"""
        return osgeo.gdal.Open(frame).ReadAsArray().mean()

    bperp_dict = {}
    for product in products:
        jobs = []
        for frame in product['bPerpendicular']:
            jobs.append(dask.delayed(read_and_average_bperp)(frame))

        mean_bperp_by_frames = dask.compute(
            jobs, num_workers=int(num_threads), scheduler='threads')[0]

        LOGGER.debug('Pair name: %s, bPerpendicular %s' % (
            product['pair_name'][0], mean_bperp_by_frames))
        bperp_dict[product['pair_name'][0]] = float(
            np.mean(mean_bperp_by_frames))
    return bperp_dict


def handle_epoch_layers(
        layers, product_dict, proj, lyr_path, user_lyrs, map_lyrs, key,
        sec_key, ref_key, tropo_total, workdir, bounds, arrres,
        dem_bounds, prods_TOTbbox, dem, lat, lon, mask, outputFormat, verbose,
        multilooking, rankedResampling, num_threads, is_nisar_file):
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

    # If specified workdir doesn't exist, create it
    all_workdirs = [workdir, sec_workdir, ref_workdir]
    all_workdirs = list(set(all_workdirs))
    existing_outputs = []
    for i in all_workdirs:
        if not os.path.exists(i):
            os.mkdir(i)

        # for dedup, record previous outputs to avoid reprocessing
        existing_outputs.extend(
            glob.glob(os.path.join(i, '*/*[0-9].vrt')))
        existing_outputs.extend(
            glob.glob(os.path.join(i, '*/dates/*[0-9].vrt')))
        existing_outputs.extend(
            glob.glob(os.path.join(i, '*[0-9].vrt')))
        existing_outputs.extend(
            glob.glob(os.path.join(i, 'dates/*[0-9].vrt')))
    existing_outputs = list(set(existing_outputs))

    # Iterate through all IFGs
    all_outputs = []
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
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)

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
                            tropo_total, hgt_field, proj, outputFormat)
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
                            tropo_total, hgt_field, proj, outputFormat)

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
                    hgt_field, proj, outputFormat)

        else:
            sec_outname = os.path.dirname(ref_outname)
            sec_outname = os.path.abspath(
                os.path.join(sec_outname, ifg.split('_')[0]))

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

            # dedup check for interpolating only new files
            record_epochs = [
                j for j in record_epochs if j not in existing_outputs]

            for j in enumerate(record_epochs):
                # Interpolate/intersect with DEM before cropping
                finalize_metadata(
                    j[1][:-4], bounds, arrres, dem_bounds, prods_TOTbbox,
                    dem, lat, lon, hgt_field, prod_ver_list, is_nisar_file,
                    outputFormat, verbose=verbose)

                # Apply mask (if specified)
                if mask is not None:
                    update_file = osgeo.gdal.Open(
                        j[1][:-4], osgeo.gdal.GA_Update)
                    mask_arr = mask.ReadAsArray() * \
                        osgeo.gdal.Open(j[1][:-4] + '.vrt').ReadAsArray()
                    update_file.GetRasterBand(1).WriteArray(mask_arr)
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

    return


def export_product_worker_helper(args):
    """Calls export_product_worker with * expanded args"""
    return export_product_worker(*args)


def export_product_worker(
        ii, ilayer, product, proj, full_product_dict_file, layers, workdir,
        bounds, prods_TOTbbox, demfile, demfile_expanded, maskfile,
        outputFormat, outputFormatPhys, layer, outDir,
        arrres, epsg_code, num_threads, multilooking, verbose, is_nisar_file,
        range_correction, rankedResampling):
    """
    Worker function for export_products for parallel execution with
    multiprocessing package.
    """
    # Initialize warp dict
    gdal_warp_kwargs = {
        'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
        'outputBounds': bounds, 'xRes': arrres[0], 'yRes': arrres[1],
        'targetAlignedPixels': True, 'multithread': True, 'dstSRS': proj}

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

    # Extract/crop metadata layers
    if (any(':/science/grids/imagingGeometry' in s for s in product) or
        any(':/science/LSAR/GUNW/metadata/radarGrid/' in s for s in product)):
        # make VRT pointing to metadata layers in standard product
        hgt_field, outname = prep_metadatalayers(
            outname, product, dem_expanded, layer, layers,
            is_nisar_file, proj)

        # Interpolate/intersect with DEM before cropping
        finalize_metadata(
            outname, bounds, arrres, dem_bounds, prods_TOTbbox, dem_expanded,
            lat, lon, hgt_field, product, is_nisar_file, outputFormatPhys,
            verbose=verbose)

    # Extract/crop full res layers, except for "unw" and "conn_comp"
    # which requires advanced stitching
    elif layer != 'unwrappedPhase' and layer != 'connectedComponents':
        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
            if outputFormat == 'VRT':
                # building the virtual vrt
                osgeo.gdal.BuildVRT(outname + "_uncropped" + '.vrt', product)

                # building the cropped vrt
                osgeo.gdal.Warp(
                    outname + '.vrt', outname + '_uncropped.vrt',
                    options=warp_options)
            else:
                # building the VRT
                osgeo.gdal.BuildVRT(outname + '.vrt', product)
                osgeo.gdal.Warp(
                    outname, outname + '.vrt', options=warp_options)

                # Update VRT
                osgeo.gdal.Translate(
                    outname + '.vrt', outname,
                    options=osgeo.gdal.TranslateOptions(format="VRT"))

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
        if (not os.path.exists(outFilePhs) or
                not os.path.exists(outFileConnComp)):

            phs_files = full_product_dict[ii]['unwrappedPhase']

            # stitching
            ARIAtools.util.seq_stitch.product_stitch_sequential(
                phs_files, conn_files, arrres=arrres, epsg=epsg_code,
                bounds=bounds, clip_json=prods_TOTbbox, output_unw=outFilePhs,
                output_conn=outFileConnComp,
                output_format=outputFormatPhys,
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
                    update_file = osgeo.gdal.Open(
                        j, osgeo.gdal.GA_Update)
                    mask_arr = mask.ReadAsArray() * \
                        osgeo.gdal.Open(j + '.vrt').ReadAsArray()
                    update_file.GetRasterBand(1).WriteArray(mask_arr)

    if layer != 'unwrappedPhase' and layer != 'connectedComponents':

        # If necessary, resample raster
        if multilooking is not None:
            ARIAtools.util.vrt.resampleRaster(
                outname, multilooking, bounds, prods_TOTbbox,
                rankedResampling, outputFormat=outputFormatPhys,
                num_threads=num_threads)

        # Apply mask (if specified)
        if mask is not None:
            update_file = osgeo.gdal.Open(
                outname, osgeo.gdal.GA_Update)
            mask_arr = mask.ReadAsArray() * \
                osgeo.gdal.Open(outname + '.vrt').ReadAsArray()
            update_file.GetRasterBand(1).WriteArray(mask_arr)
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
        is_nisar_file, rankedResampling=False, demfile=None,
        demfile_expanded=None, lat=None, lon=None, maskfile=None, outDir='./',
        outputFormat='VRT', verbose=None, num_threads='2', multilooking=None,
        tropo_total=False, model_names=[], multiproc_method='single'):
    """
    Export layer and 2D meta-data layers (at the product resolution).
    The function finalize_metadata is called to derive the 2D metadata layer.
    Dem/lat/lon arrays must be passed for this process.
    The keys specify which layer to extract from the dictionary.
    All products are cropped by the bounds from the input bbox_file,
    and clipped to the track extent denoted by the input prods_TOTbbox.
    Optionally, a user may pass a mask-file.
    """
    LOGGER.debug('export_products, layers: {}'.format(layers))
    # Remove these directories to avoid state dependent VRT parsing bug
    # with pre-existing files (TODO to fix this)
    for layer in layers:
        if layer in GEOM_LYRS:
            target = os.path.join(outDir, layer)
            if os.path.isdir(target):
                LOGGER.warning('Deleting %s to avoid VRT header bug!' % target)
                shutil.rmtree(target)

    if not layers and not tropo_total:
        return  # only bbox

    # initiate tracker of output dimensions
    ref_wid = None
    ref_hgt = None
    ref_geotrans = None

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
    epsg_code = srs.GetAuthorityCode(None)
    srs = None
    lyr_input_dict = {
        'layers': layers, 'prods_TOTbbox': prods_TOTbbox,
        'proj': int(epsg_code), 'dem': dem_expanded, 'lat': lat, 'lon': lon,
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
            prev_outname = os.path.abspath(os.path.join(workdir, i))
            lyr_input_dict['product_dict'] = product_dict

            # extract layers
            handle_epoch_layers(**lyr_input_dict)

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
        handle_epoch_layers(**lyr_input_dict)

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
                              verbose=verbose,
                              overwrite=True)

        for i, layer in enumerate(product_dict[0]):
            outname = os.path.abspath(
                os.path.join(
                    workdir,
                    product_dict[1][i][0]))
            lyr_input_dict['input_iono_files'] = layer
            lyr_input_dict['output_iono'] = outname
            ARIAtools.util.ionosphere.export_ionosphere(**lyr_input_dict)

            # track valid files
            if os.path.exists(outname + '.vrt'):
                prev_outname_check = copy.deepcopy(outname)

        # track consistency of dimensions
        if 'prev_outname_check' in locals():
            ref_wid, ref_hgt, ref_geotrans, _, _ = \
                ARIAtools.util.vrt.get_basic_attrs(prev_outname_check + '.vrt')
            ref_arr = [ref_wid, ref_hgt, ref_geotrans, prev_outname]

    # Loop through other user expected layers
    layers = [i for i in layers if i not in ext_corr_lyrs]

    full_product_dict_file = os.path.join(outDir, 'full_product_dict.json')
    with open(full_product_dict_file, 'w') as ofp:
        json.dump(full_product_dict, ofp)

    mp_args = []
    for ilayer, layer in enumerate(layers):

        product_dict = [[j[layer] for j in full_product_dict],
                        [j["pair_name"] for j in full_product_dict]]

        # If specified workdir doesn't exist, create it
        workdir = os.path.join(outDir, layer)
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        # Iterate through all IFGs
        # TODO can we wrap this into funtion and run it
        # with multiprocessing, to gain speed up
        for ii, product in enumerate(product_dict[0]):
            mp_args.append((
                ii, ilayer, product, proj, full_product_dict_file, layers,
                workdir, bounds, prods_TOTbbox, demfile,
                demfile_expanded, maskfile, outputFormat, outputFormatPhys,
                layer, outDir, arrres, epsg_code, num_threads,
                multilooking, verbose, is_nisar_file, range_correction,
                rankedResampling))

    start_time = time.time()
    if int(num_threads) == 1 or multiproc_method in ['single', 'threads']:
        if multiproc_method == 'single':
            outputs = []
            for arg in tqdm.tqdm(mp_args, total=len(mp_args), desc='Exporting'):
                outputs.append(export_product_worker_helper(arg))
                sys.stdout.flush()
        else:
            LOGGER.debug('Running %d total jobs with threads' % len(mp_args))

            # Create a progress bar
            with tqdm.tqdm(total=len(mp_args), desc='Exporting') as pbar:
                # Create jobs with a wrapper function that includes progress update
                jobs = []
                for arg in mp_args:
                    # Create a delayed job that includes both work and progress update
                    job = dask.delayed(lambda x: (export_product_worker(*x),
                                                  pbar.update(1))[0])(arg)
                    jobs.append(job)

                # Compute all jobs
                outputs = dask.compute(jobs, num_workers=int(num_threads),
                                       scheduler='threads')[0]

        for ii, ilayer, outname, prod_arr in outputs:
            if ii == 0 and ilayer == 0:
                ref_arr = copy.deepcopy(prod_arr)
            else:
                ARIAtools.util.vrt.dim_check(ref_arr, prod_arr)
            prev_outname = outname

    elif multiproc_method == 'gnu_parallel':
        export_workers_temp_dir = os.path.join(outDir, 'export_workers')
        if os.path.isdir(export_workers_temp_dir):
            shutil.rmtree(export_workers_temp_dir)
        os.mkdir(export_workers_temp_dir)

        for ii, args in enumerate(mp_args):
            this_json_file = os.path.join(
                outDir, 'export_workers',
                'export_product_args_%2.2d.json' % ii)
            with open(this_json_file, 'w') as ofp:
                json.dump(args, ofp)

        LOGGER.debug('Running %d total jobs in parallel' % len(mp_args))
        # Run the export worker jobs with GNU parallel
        subprocess.call((
            'find %s/export_workers -name "export_product_args_*.json" | '
            'parallel -j %d export_product.py {}') % (
                outDir, int(num_threads)), shell=True)
        end_time = time.time()

        # load in output files and verify dimensions
        output_files = glob.glob(os.path.join(
            export_workers_temp_dir, 'outputs_*.json'))

        if len(output_files) > 0:
            with open(os.path.join(
                    export_workers_temp_dir, 'outputs_0_0.json')) as ifp:
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

    # delete directory for quality control plots if empty
    plots_subdir = os.path.abspath(
        os.path.join(outDir, 'metadatalyr_plots'))
    if os.path.exists(plots_subdir) and len(os.listdir(plots_subdir)) == 0:
        shutil.rmtree(plots_subdir)

    try:
        retval = ref_arr
    except UnboundLocalError:
        retval = [None, None, None, None]
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

    # load layered metadata array
    tmp_name = outname + '.vrt'
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(format="MEM")
        data_array = osgeo.gdal.Warp('', tmp_name, options=warp_options)

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

    metadatalyr_name = outname.split('/')[-2]
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

    # only perform DEM intersection for rasters with valid height levels
    NOHGT_LYRS = ['ionosphere']
    if metadatalyr_name not in NOHGT_LYRS:
        tmp_name = outname + '_temp'

        # Define lat/lon/height arrays for metadata layers
        heightsMeta = np.array(
            osgeo.gdal.Open(outname + '.vrt').GetMetadataItem(
                hgt_field)[1:-1].split(','), dtype='float32')

        latitudeMeta = np.linspace(
            data_array.GetGeoTransform()[3],
            data_array.GetGeoTransform()[3] +
            (data_array.GetGeoTransform()[5] * (data_array.RasterYSize - 1)),
            data_array.RasterYSize, dtype='float32')

        longitudeMeta = np.linspace(
            data_array.GetGeoTransform()[0],
            data_array.GetGeoTransform()[0] +
            (data_array.GetGeoTransform()[1] * (data_array.RasterXSize - 1)),
            data_array.RasterXSize, dtype='float32')

        da_dem = rioxarray.open_rasterio(
            dem.GetDescription(), band_as_variable=True,
            masked=True)['band_1']

        # interpolate the DEM to the GUNW lat/lon
        nodata = dem.GetRasterBand(1).GetNoDataValue()
        da_dem1 = da_dem.interp(
            x=lon[0, :], y=lat[:, 0]).fillna(nodata)

        # hack to get an stack of coordinates for the interpolator
        # to interpolate in the right shape
        pnts = transformPoints(
            lat, lon, da_dem1.data, 'EPSG:4326', 'EPSG:4326')

        # set up the interpolator with the GUNW cube
        data_array_inp = data_array.ReadAsArray().astype('float32')
        interper = scipy.interpolate.RegularGridInterpolator(
            (latitudeMeta, longitudeMeta, heightsMeta),
            data_array_inp.transpose(1, 2, 0),
            fill_value=np.nan, bounds_error=False)

        # interpolate cube to DEM points
        out_interpolated = interper(pnts.transpose(2, 1, 0))

        # Save file
        ARIAtools.util.vrt.renderVRT(
            tmp_name, out_interpolated, geotrans=dem.GetGeoTransform(),
            drivername=outputFormat,
            gdal_fmt=data_array.ReadAsArray().dtype.name,
            proj=dem.GetProjection(), nodata=nodata)
        out_interpolated = None

    # Since metadata layer extends at least one grid node
    # outside of the expected track bounds,
    # it must be cut to conform with these bounds.
    # Crop to track extents
    data_array_nodata = data_array.GetRasterBand(1).GetNoDataValue()
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
            'outputBounds': dem_bounds, 'dstNodata': data_array_nodata,
            'xRes': dem_arrres[0], 'yRes': dem_arrres[1],
            'targetAlignedPixels': True, 'multithread': True}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        osgeo.gdal.Warp(tmp_name + '_temp', tmp_name, options=warp_options)

    # Adjust shape
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
            'outputBounds': bbox_bounds, 'dstNodata': data_array_nodata,
            'xRes': arrres[0], 'yRes': arrres[1], 'targetAlignedPixels': True,
            'multithread': True}
        warp_options = osgeo.gdal.WarpOptions(**gdal_warp_kwargs)
        osgeo.gdal.Warp(outname, tmp_name + '_temp', options=warp_options)

    # remove temp files
    for i in glob.glob(outname + '*_temp*'):
        os.remove(i)

    # Update VRT
    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    osgeo.gdal.Translate(
        outname + '.vrt', outname, options=translate_options)

    data_array = None


def gacos_correction(full_product_dict, gacos_products, bbox_file,
                     prods_TOTbbox, outDir='./', outputFormat='VRT',
                     verbose=None, num_threads='2'):
    """Perform tropospheric corrections.
    Must provide valid path to GACOS products.
    All products are cropped by the bounds from the input bbox_file,
    and clipped to the track extent denoted by the input prods_TOTbbox.
    """
    # File must be physically extracted, cannot proceed with VRT format.
    # Defaulting to ENVI format.
    outputFormat = 'ENVI' if outputFormat == 'VRT' else outputFormat

    user_bbox = ARIAtools.util.shp.open_shp(bbox_file)
    bounds = user_bbox.bounds

    product_dict = [[j['unwrappedPhase'] for j in full_product_dict[1]],
                    [j['lookAngle'] for j in full_product_dict[1]],
                    [j["pair_name"] for j in full_product_dict[1]]]
    metadata_dict = [[j['azimuthZeroDopplerMidTime'] for j in
                      full_product_dict[0]],
                     [j['wavelength'] for j in full_product_dict[0]]]
    workdir = os.path.join(outDir, 'gacos_corrections')

    # If specified workdir doesn't exist, create it
    os.makedirs(workdir, exist_ok=True)

    # Get list of all dates for which standard products exist
    date_list = []
    for i in product_dict[2]:
        date_list.append(i[0][:8])
        date_list.append(i[0][9:])
    date_list = list(set(date_list))  # trim to unique dates only

    # Determine if file input is single file, a list, or wildcard.
    # If list of files
    if len([str(val) for val in gacos_products.split(',')]) > 1:
        gacos_products = [str(i) for i in gacos_products.split(',')]
        # Sort and parse tropo products
        tropo_sublist = []
        for i in gacos_products:
            # If wildcard
            if '*' in i:
                tropo_sublist.append(
                    glob.glob(os.path.expanduser(os.path.expandvars(i))))
            else:
                tropo_sublist.append([i])
        # Finalize list of tropo products
        gacos_products = []
        for sublist in tropo_sublist:
            for item in sublist:
                gacos_products.append(os.path.abspath(item))

    # If single file or wildcard
    else:
        # If single file
        if os.path.isfile(gacos_products):
            gacos_products = [gacos_products]
        # If wildcard
        else:
            gacos_products = glob.glob(os.path.expanduser(
                                       os.path.expandvars(gacos_products)))
        # Convert relative paths to absolute paths
        gacos_products = [os.path.abspath(i) for i in gacos_products]
    if len(gacos_products) == 0:
        raise Exception('No file match found')

    # Extract tarfiles
    # Setup dictionary to track for products that are to be merged
    tropo_date_dict = {}
    for i in date_list:
        tropo_date_dict[i] = []
        tropo_date_dict[i + "_UTC"] = []
    for i in enumerate(gacos_products):
        if not os.path.isdir(i[1]) and not i[1].endswith('.ztd') \
                and not i[1].endswith('.tif'):
            untar_dir = os.path.join(os.path.abspath(os.path.join(
                i[1], os.pardir)),
                os.path.basename(i[1]).split('.')[0] + '_extracted')
            if not tarfile.is_tarfile(i[1]):
                raise Exception('Cannot extract %s because it is not a '
                                'valid tarfile. Resolve this '
                                'and relaunch' % (i[1]))
            LOGGER.info('Extracting GACOS tarfile %s to %s.',
                        os.path.basename(i[1]), untar_dir)
            tarfile.open(i[1]).extractall(path=untar_dir)
            gacos_products[i[0]] = untar_dir
        # Loop through each GACOS product file, differentiating between direct
        # input list of GACOS products vs parent directory
        if i[1].endswith('.ztd') or i[1].endswith('.tif'):
            ztd_list = [i[1]]
        else:
            ztd_list = glob.glob(os.path.join(gacos_products[i[0]], '*.ztd')) \
                + glob.glob(os.path.join(gacos_products[i[0]], '*.tif'))
            # prioritize older .ztd over .tif duplicates if the former exists
            # older .ztd files contain UTC time info
            ztd_list = [i for i in ztd_list if not (
                i.endswith('.tif') and i.split('.tif')[0] in ztd_list)]
        for k in ztd_list:
            # Only check files corresponding to standard product dates
            ztd_basename = os.path.basename(k)
            ztd_basename = ztd_basename.split('.tif')[0].split('.ztd')[0]
            if ztd_basename in date_list:
                tropo_date_dict[ztd_basename].append(k)
                if os.path.exists(k + '.rsc'):
                    tropo_date_dict[ztd_basename + "_UTC"].append(
                        ztd_basename[:4] + '-'
                        + ztd_basename[4:6]
                        + '-' + ztd_basename[6:8] + '-'
                        + open(k + '.rsc', 'r').readlines()[-1].split()[1])
                else:
                    tropo_date_dict[ztd_basename + "_UTC"].append(
                        ztd_basename[:4] + '-'
                        + ztd_basename[4:6]
                        + '-' + ztd_basename[6:8] + '-'
                        + 'NoneUTC')
                # make corresponding VRT file, if it doesn't exist
                if not os.path.exists(k + '.vrt'):
                    # parse metadata from RSC if it exists
                    if os.path.exists(k + '.rsc'):
                        tropo_rsc_dict = {}
                        for line in open(k + '.rsc', 'r').readlines():
                            tropo_rsc_dict[line.split()[0]] = line.split()[1]
                        gacos_prod = np.fromfile(k, dtype='float32').reshape(
                            int(tropo_rsc_dict['FILE_LENGTH']),
                            int(tropo_rsc_dict['WIDTH']))
                    else:
                        tropo_rsc_dict = ARIAtools.util.vrt.tifGacos(k)
                        gacos_prod = osgeo.gdal.Open(k).ReadAsArray()
                    # Save as GDAL file, using proj from first unwrappedPhase
                    # file
                    geotrans = (float(tropo_rsc_dict['X_FIRST']),
                                float(tropo_rsc_dict['X_STEP']),
                                0.0, float(tropo_rsc_dict['Y_FIRST']),
                                0.0, float(tropo_rsc_dict['Y_STEP']))
                    proj = osgeo.gdal.Open(os.path.join(
                        outDir, 'unwrappedPhase',
                        product_dict[2][0][0])).GetProjection()
                    ARIAtools.util.vrt.renderVRT(
                        k, gacos_prod, geotrans=geotrans,
                        drivername=outputFormat, gdal_fmt='float32',
                        proj=proj, nodata=0.)
                    gacos_prod = None
                    LOGGER.debug('GACOS product %s successfully converted to '
                                 'GDAL-readable raster', k)
                # make corresponding RSC file, if it doesn't exist
                if not os.path.exists(k + '.rsc'):
                    ARIAtools.util.vrt.rscGacos(
                        os.path.join(k + '.vrt'), os.path.join(k + '.rsc'),
                        tropo_date_dict)

    # If multiple GACOS directories, merge products.
    gacos_products = list(set([os.path.dirname(i)
                          if (i.endswith('.ztd') or i.endswith('.tif'))
                          else i for i in gacos_products]))
    if len(gacos_products) > 1:
        gacos_products = os.path.join(outDir, 'merged_GACOS')
        LOGGER.info('Stitching/storing GACOS products in %s.', gacos_products)
        # If specified merged directory doesn't exist, create it
        if not os.path.exists(os.path.join(outDir, 'merged_GACOS')):
            os.mkdir(os.path.join(outDir, 'merged_GACOS'))

        for i in tropo_date_dict:
            # only make rsc/vrt files if valid product
            if 'UTC' not in i and tropo_date_dict[i] != []:
                outname = os.path.join(outDir, 'merged_GACOS', i + '.ztd.vrt')
                # building the VRT
                osgeo.gdal.BuildVRT(outname, tropo_date_dict[i])
                geotrans = osgeo.gdal.Open(outname).GetGeoTransform()
                # Create merged rsc file
                ARIAtools.util.vrt.rscGacos(
                    outname, merged_rsc, tropo_date_dict)
    else:
        gacos_products = gacos_products[0]

    # Estimate percentage of overlap with tropospheric product
    for i in glob.glob(os.path.join(gacos_products, '*.vrt')):
        # create shapefile
        geotrans = osgeo.gdal.Open(i).GetGeoTransform()
        bbox = [geotrans[3]
                + (osgeo.gdal.Open(i).ReadAsArray().shape[0] * geotrans[-1]),
                geotrans[3],
                geotrans[0],
                geotrans[0]
                + (osgeo.gdal.Open(i).ReadAsArray().shape[1] * geotrans[1])
                ]
        bbox = shapely.geometry.Polygon(np.column_stack((
            np.array([bbox[2], bbox[3], bbox[3], bbox[2], bbox[2]]),
            np.array([bbox[0], bbox[0], bbox[1], bbox[1], bbox[0]]))))
        ARIAtools.util.shp.save_shp(i + '.json', bbox)
        per_overlap = ((user_bbox.intersection(
            ARIAtools.util.shp.open_shp(i + '.json')).area)
            / (user_bbox.area)) * 100
        if per_overlap != 100. and per_overlap != 0.:
            LOGGER.warning('Common track extent only has %d overlap with'
                           'tropospheric product %s\n', per_overlap, i[0])
        if per_overlap == 0.:
            raise Exception('No spatial overlap between tropospheric '
                            'product %s and defined bounding box. '
                            'Resolve conflict and relaunch', i[1])

    # Iterate through all IFGs and apply corrections
    missing_products = []
    for i in range(len(product_dict[0])):
        ifg = product_dict[2][i][0]
        outname = os.path.join(workdir, ifg)
        gacos_epochs_dir = os.path.join(workdir, 'dates')
        ref_outname = os.path.join(gacos_epochs_dir, f'{ifg[:8]}')
        sec_outname = os.path.join(gacos_epochs_dir, f'{ifg[9:]}')
        outname = os.path.join(workdir, ifg)
        unwname = os.path.join(outDir, 'unwrappedPhase', ifg)
        if i == 0:
            meta = osgeo.gdal.Info(unwname, format='json')
            geoT = meta['geoTransform']
            proj = meta['coordinateSystem']['wkt']
            arrres = [abs(geoT[1]), abs(geoT[-1])]

        tropo_reference = os.path.join(gacos_products, f'{ifg[:8]}.ztd.vrt')
        tropo_secondary = os.path.join(gacos_products, f'{ifg[9:]}.ztd.vrt')

        # if .ztd products don't exist, check if .tif exists
        if not os.path.exists(tropo_reference):
            tropo_reference = os.path.join(
                gacos_products, f'{ifg[:8]}.ztd.tif.vrt')
        if not os.path.exists(tropo_secondary):
            tropo_secondary = os.path.join(
                gacos_products, f'{ifg[9:]}.ztd.tif.vrt')

        # skip if corrected already generated and does not need to be updated
        if os.path.exists(outname):
            # get unwrappedPhase geotrans and productbounding box
            unw_prodcheck = osgeo.gdal.Open(unwname)
            unw_geotrans = unw_prodcheck.GetGeoTransform()
            unw_prodcheck = np.isfinite(unw_prodcheck.ReadAsArray())
            tropo_prodcheck = osgeo.gdal.Open(outname)
            output_geotrans = tropo_prodcheck.GetGeoTransform()
            tropo_prodcheck = np.isfinite(tropo_prodcheck.ReadAsArray())
            if unw_geotrans == output_geotrans and np.array_equal(
                    unw_prodcheck, tropo_prodcheck):
                continue
            unw_prodcheck = None
            tropo_prodcheck = None

        if os.path.exists(tropo_reference) and os.path.exists(tropo_secondary):
            # Check if tropo products are temporally consistent with IFG
            for j in [tropo_reference, tropo_secondary]:
                # Get ARIA product times
                aria_rsc_dict = {}
                aria_rsc_dict['azimuthZeroDopplerMidTime'] = [
                    datetime.datetime.strptime(
                        os.path.basename(j)[: 4] + '-' + os.path.
                        basename(j)[4: 6] + '-' + os.path.basename(j)[6: 8] +
                        '-' + m[11:],
                        "%Y-%m-%d-%H:%M:%S.%f") for m in metadata_dict[0][0]]
                # Get tropo product UTC times
                tropo_rsc_dict = {}
                tropo_rsc_dict['TIME_OF_DAY'] = open(
                    j[:-4] + '.rsc',
                    'r').readlines()[-1].split()[1].split('UTC')[:-1]
                # If new TIF product, UTC times not available
                if 'None' in tropo_rsc_dict['TIME_OF_DAY'][0]:
                    tropo_rsc_dict['TIME_OF_DAY'] = [
                        max(aria_rsc_dict['azimuthZeroDopplerMidTime'])]
                # If stitched tropo product, must account for date change (if
                # applicable)
                elif '-' in tropo_rsc_dict['TIME_OF_DAY'][0]:
                    tropo_rsc_dict['TIME_OF_DAY'] = [
                        datetime.datetime.strptime(
                            m[: 10] + '-' + m[11:].split('.')[0] + '-' +
                            str(
                                round(
                                    float('0.' + m[11:].split('.')[1]) * 60)),
                            "%Y-%m-%d-%H-%M")
                        for m in tropo_rsc_dict['TIME_OF_DAY']]
                else:
                    tropo_rsc_dict['TIME_OF_DAY'] = [
                        datetime.datetime.strptime(
                            os.path.basename(j)[: 4] + '-' + os.path.
                            basename(j)[4: 6] + '-' + os.path.basename(j)
                            [6: 8] + '-' +
                            tropo_rsc_dict['TIME_OF_DAY'][0].split('.')[0] +
                            '-' +
                            str(
                                round(
                                    float(
                                        '0.' +
                                        tropo_rsc_dict['TIME_OF_DAY'][0].
                                        split('.')[-1]) * 60)),
                            "%Y-%m-%d-%H-%M")]

                # Check and report if tropospheric product falls outside of
                # standard product range
                latest_start = max(aria_rsc_dict['azimuthZeroDopplerMidTime']
                                   + [min(tropo_rsc_dict['TIME_OF_DAY'])])
                earliest_end = min(aria_rsc_dict['azimuthZeroDopplerMidTime']
                                   + [max(tropo_rsc_dict['TIME_OF_DAY'])])
                delta = (earliest_end - latest_start).total_seconds() + 1
                if delta < 0:
                    LOGGER.warning('tropospheric product was generated %f '
                                   'secs outside of acquisition interval for '
                                   'scene %s in IFG %s',
                                   abs(delta), os.path.basename(j)[:8],
                                   product_dict[2][i][0])

            # Open corresponding tropo products and pass the difference
            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
                gdal_warp_kwargs = {'format': outputFormat,
                                    'cutlineDSName': prods_TOTbbox,
                                    'outputBounds': bounds,
                                    'xRes': arrres[0],
                                    'yRes': arrres[1],
                                    'targetAlignedPixels': True,
                                    'multithread': True}
                tropo_reference = osgeo.gdal.Warp(
                    '', tropo_reference, options=osgeo.gdal.WarpOptions(
                        **gdal_warp_kwargs)).ReadAsArray()
                tropo_secondary = osgeo.gdal.Warp(
                    '', tropo_secondary, options=osgeo.gdal.WarpOptions(
                        **gdal_warp_kwargs)).ReadAsArray()
                tropo_product = np.subtract(tropo_secondary, tropo_reference)

            # Convert troposphere from m to rad
            scale = float(metadata_dict[1][i][0]) / (4 * np.pi)
            tropo_product /= scale

            # Account for incAngle
            # if in TS mode, only 1 incfile would be generated, so check for
            # this
            path_inc = os.path.join(outDir, 'incidenceAngle', ifg)
            if os.path.exists(path_inc):
                da = rasterio.open(path_inc)
            else:
                da = rasterio.open(
                    path_inc.replace(
                        ifg, product_dict[2][0][0]))
            inc_arr = da.read().squeeze()
            inc_arr = np.where(np.isclose(inc_arr, da.nodata), np.nan, inc_arr)
            cos_inc = np.cos(np.deg2rad(inc_arr))

            tropo_product /= cos_inc

            # Save differential field to file
            tropo_product = np.where(
                np.isnan(tropo_product), 0., tropo_product)
            ARIAtools.util.vrt.renderVRT(
                outname, tropo_product, geotrans=geoT, drivername=outputFormat,
                gdal_fmt='float32', proj=proj, nodata=0.)

            # check if reference and secondary scenes are written to file
            if not os.path.exists(ref_outname):
                tropo_reference /= scale
                tropo_reference /= cos_inc
                tropo_reference = np.where(
                    np.isnan(tropo_reference), 0., tropo_reference)
                ARIAtools.util.vrt.renderVRT(
                    ref_outname, tropo_reference, geotrans=geoT,
                    drivername=outputFormat, gdal_fmt='float32', proj=proj,
                    nodata=0.)

            if not os.path.exists(sec_outname):
                tropo_secondary /= scale
                tropo_secondary /= cos_inc
                tropo_secondary = np.where(
                    np.isnan(tropo_secondary), 0., tropo_secondary)
                ARIAtools.util.vrt.renderVRT(
                    sec_outname, tropo_secondary, geotrans=geoT,
                    drivername=outputFormat, gdal_fmt='float32', proj=proj,
                    nodata=0.)

            tropo_product = None
            tropo_reference = None
            tropo_secondary = None
            da = None
            inc_arr = None
            cos_inc = None

            # Track consistency of dimensions
            if i[0] == 0:
                ref_wid, ref_hgt, ref_geotrans, \
                    _, _ = ARIAtools.util.vrt.get_basic_attrs(outname)
                ref_arr = [ref_wid, ref_hgt, ref_geotrans,
                           os.path.join(workdir, ifg)]
            else:
                prod_wid, prod_hgt, prod_geotrans, \
                    _, _ = ARIAtools.util.vrt.get_basic_attrs(outname)
                prod_arr = [prod_wid, prod_hgt, prod_geotrans,
                            os.path.join(workdir, ifg)]
                ARIAtools.util.vrt.dim_check(ref_arr, prod_arr)
            prev_outname = os.path.join(workdir, ifg)

        else:
            LOGGER.warning('Must skip IFG %s, because the tropospheric '
                           'products corresponding to the reference and/or '
                           'secondary products are not found in the '
                           'specified folder %s',
                           ifg, gacos_products)
            for j in [tropo_reference, tropo_secondary]:
                if not os.path.exists(j) and j not in missing_products:
                    missing_products.append(j)
    # Print list of dates missing tropospheric corrections
    if len(missing_products) > 0:
        missing_products = [os.path.basename(i)[:8] for i in missing_products]
        LOGGER.debug(
            "Tropo products for the following dates are missing:",
            missing_products)


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
