# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & Brett Buzzanga & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import glob
import logging
import os
import shutil
from time import sleep

import affine
import copy
import osgeo.gdal
import pyproj
import rasterio
import tile_mate

import ARIAtools.util.shp
import ARIAtools.util.vrt


LOGGER = logging.getLogger(__name__)


def _remove_raster_files(path):
    """Remove a raster and common GDAL sidecars."""
    for suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        candidate = path + suffix
        if os.path.isfile(candidate):
            os.remove(candidate)


def _validate_raster(path, context, expected_driver=None):
    """Reject missing, truncated, or incorrectly formatted mask rasters."""
    ds = osgeo.gdal.Open(path, osgeo.gdal.GA_ReadOnly)
    if ds is None:
        raise RuntimeError(f'Could not open {context}: {path}')
    try:
        if ds.RasterCount < 1 or ds.RasterXSize < 1 or ds.RasterYSize < 1:
            raise RuntimeError(f'{context} has invalid dimensions: {path}')
        actual_driver = ds.GetDriver().ShortName
        if (expected_driver is not None
                and actual_driver.upper() != expected_driver.upper()):
            raise RuntimeError(
                f'{context} driver mismatch: expected {expected_driver}, '
                f'created {actual_driver}: {path}')
        band = ds.GetRasterBand(1)
        if (band.ReadRaster(0, 0, ds.RasterXSize, 1) is None
                or band.ReadRaster(
                    0, ds.RasterYSize - 1, ds.RasterXSize, 1) is None):
            raise RuntimeError(f'{context} is incomplete: {path}')
    finally:
        ds = None


def _translate_mask_stage(stage_path, output_path, output_format, context):
    """Translate a verified GeoTIFF mask stage to the requested driver."""
    _validate_raster(stage_path, f'{context} GeoTIFF stage', 'GTiff')
    _remove_raster_files(output_path)
    ds_out = osgeo.gdal.Translate(
        output_path, stage_path, format=output_format)
    if ds_out is None:
        raise RuntimeError(
            f'Failed to translate {context} to {output_format}: '
            f'{output_path}')
    ds_out.FlushCache()
    ds_out = None
    _validate_raster(output_path, context, output_format)


def prep_mask(
        product_dict, maskfilename, bbox_file, prods_TOTbbox, proj,
        amp_thresh=None, arrres=None, workdir='./', outputFormat='ENVI',
        num_threads='2', multilooking=None, rankedResampling=False,
        runlog=None):
    """
    Function to load and export mask file with tile_mate
    """
    LOGGER.debug("prep_mask")

    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir, 'mask')
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)

    # Get bounds of user bbox_file
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds

    # File must be physically extracted, cannot proceed with VRT format
    # Defaulting to ENVI format
    if outputFormat == 'VRT':
        outputFormat = 'ENVI'

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    # Retrieve update mode
    update_mode = 'full_extract'
    if runlog is not None:
        log_data = runlog.load()
        if 'update_mode' in log_data.keys():
            update_mode = log_data['update_mode']

    # set temp directory
    temp_workdir = os.path.join(workdir, 'tmp_dir')

    # delete temporary directory
    if os.path.exists(temp_workdir):
        shutil.rmtree(temp_workdir)

    # Detect if mask is ESA WorldCover download or pre-existing ESA mask
    is_esa_mask = 'esa_world_cover' in os.path.basename(maskfilename).lower()

    if maskfilename.lower() == 'download' or \
            maskfilename.lower() in tile_mate.stitcher.DATASET_SHORTNAMES or \
            is_esa_mask:
        # if download specified or auto-detected esa mask, default to esa world cover mask
        if maskfilename.lower() == 'download' or is_esa_mask:
            maskfilename = 'esa_world_cover_2021'
        lyr_name = copy.deepcopy(maskfilename)
        LOGGER.info('Preparing water mask: %s', lyr_name)

        # set file names
        uncropped_maskfilename = os.path.join(workdir,
                                              f'{maskfilename}_uncropped.tif')
        maskfilename = os.path.join(workdir, f'{maskfilename}.msk')
        ref_file = os.path.join(workdir, 'tmp_referencefile.tif')

        # Check if uncropped mask exists and is valid
        uncropped_valid = (
            os.path.exists(uncropped_maskfilename)
            and os.path.getsize(uncropped_maskfilename) > 1000
        )

        if uncropped_valid and update_mode != 'full_extract':
            LOGGER.warning(
                '%s has already been downloaded. Skipping download.',
                uncropped_maskfilename)
        else:
            # download mask
            dat_arr, dat_prof = tile_mate.get_raster_from_tiles(
                bounds, tile_shortname=lyr_name)

            # fill permanent water body
            if lyr_name in ('esa_world_cover_2020', 'esa_world_cover_2021'):
                dat_arr[dat_arr == 80] = 0
                dat_arr[dat_arr != 0] = 1

            # assign datatype and set resampling mode
            dat_arr = dat_arr.astype('byte')
            f_dtype = 'uint8'
            resampling_mode = rasterio.warp.Resampling.nearest

            # get output parameters from temp file
            crs = pyproj.CRS.from_wkt(proj)
            if os.path.exists(ref_file):
                os.remove(ref_file)

            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
                ds_ref = osgeo.gdal.Warp(
                    ref_file, product_dict[0], format='GTiff',
                    outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                    targetAlignedPixels=True, multithread=False)
                if ds_ref is None:
                    raise RuntimeError(
                        f'Failed to create mask reference raster: {ref_file}')
                ds_ref.FlushCache()
                ds_ref = None
            _validate_raster(ref_file, 'mask reference raster', 'GTiff')

            with rasterio.open(ref_file) as src:
                reference_gt = src.transform
                resize_col = src.width
                resize_row = src.height

            # remove temporary file
            for j in glob.glob(ref_file + '*'):
                if os.path.isfile(j):
                    os.remove(j)

            # save uncropped raster to file
            with rasterio.open(uncropped_maskfilename, 'w', driver='GTiff',
                               height=resize_row, width=resize_col, count=1,
                               dtype=f_dtype, crs=crs,
                               transform=affine.Affine(*reference_gt)) as dst:
                rasterio.warp.reproject(
                    source=dat_arr, destination=rasterio.band(dst, 1),
                    src_transform=dat_prof['transform'],
                    src_crs=dat_prof['crs'], dst_transform=reference_gt,
                    dst_crs=crs, resampling=resampling_mode)
            _validate_raster(
                uncropped_maskfilename, 'uncropped downloaded mask', 'GTiff')

        # save cropped mask with precise spacing
        tmp_crop_mask = f'{maskfilename}.tmp.tif'
        if os.path.exists(tmp_crop_mask):
            os.remove(tmp_crop_mask)

        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            ds_crop = osgeo.gdal.Warp(
                tmp_crop_mask, uncropped_maskfilename, format='GTiff',
                outputBounds=bounds, outputType=osgeo.gdal.GDT_Byte,
                xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
                multithread=False, options=['-overwrite'])
            if ds_crop is None:
                raise RuntimeError(
                    f'Failed to create downloaded mask staging raster: '
                    f'{tmp_crop_mask}')
            ds_crop.FlushCache()
            ds_crop = None

        # Keep subsequent spatial operations on the verified GeoTIFF stage.
        # The requested physical format is written only once, at the end.
        mask_processing_source = tmp_crop_mask
        update_file = osgeo.gdal.Open(
            mask_processing_source, osgeo.gdal.GA_Update)
        if update_file is None:
            raise RuntimeError(
                f'Could not update downloaded mask stage: '
                f'{mask_processing_source}')
        update_file.SetProjection(proj)
        update_file.GetRasterBand(1).SetNoDataValue(0.)
        update_file.FlushCache()
        update_file = None

    # User specified mask
    else:
        LOGGER.info("Using user specified mask %s" % maskfilename)
        user_mask = os.path.abspath(maskfilename)
        user_mask_n = os.path.basename(os.path.splitext(user_mask)[0])
        local_mask = os.path.join(workdir, f'{user_mask_n}.msk')
        local_mask_unc = os.path.join(workdir, f'{user_mask_n}_uncropped.msk')

        if user_mask == local_mask:
            LOGGER.debug(
                'The mask you specified already exists in %s, '
                'using the existing one...' % os.path.dirname(local_mask))

            os.makedirs(temp_workdir, exist_ok=True)
            local_mask_noext = os.path.join(workdir, '%s.' % (user_mask_n))
            for j in glob.glob(local_mask_noext + '*'):
                if not j.startswith(temp_workdir):
                    shutil.move(j, temp_workdir)

            temp_local_mask = os.path.join(
                temp_workdir, '%s.msk' % (user_mask_n))

            if not os.path.exists(temp_local_mask) or os.path.getsize(temp_local_mask) < 100:
                raise RuntimeError(
                    f"Corrupted local mask found at {temp_local_mask}. "
                    "Please remove TS_extract_masked/mask directory and re-run."
                )

            ds = osgeo.gdal.Open(temp_local_mask)

        else:
            osgeo.gdal.UseExceptions()

            ds = osgeo.gdal.BuildVRT(
                f'{local_mask_unc}.vrt', [user_mask], outputBounds=bounds)
            assert ds is not None, f'Could not open user mask: {user_mask}'

        tmp_user_crop = f'{local_mask}.tmp.tif'
        if os.path.exists(tmp_user_crop):
            os.remove(tmp_user_crop)

        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            ds_crop = osgeo.gdal.Warp(
                tmp_user_crop, ds, format='GTiff',
                cutlineDSName=prods_TOTbbox, outputBounds=bounds,
                xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
                multithread=False, options=['-overwrite'])
            if ds_crop is None:
                raise RuntimeError(
                    f'Failed to create user-mask staging raster: '
                    f'{tmp_user_crop}')
            ds_crop.FlushCache()
            ds_crop = None

        mask_processing_source = tmp_user_crop
        mask_file = osgeo.gdal.Open(
            mask_processing_source, osgeo.gdal.GA_Update)
        if mask_file is None:
            raise RuntimeError(
                f'Could not update user-mask stage: '
                f'{mask_processing_source}')
        mask_file.SetProjection(proj)
        mask_file.FlushCache()
        mask_file = None
        maskfilename = local_mask

    # Make average amplitude mask
    if amp_thresh is not None:
        amp_file = ARIAtools.util.vrt.rasterAverage(
            os.path.join(workdir, 'avgamplitude'), product_dict, bounds,
            prods_TOTbbox, arrres, outputFormat=outputFormat,
            thresh=amp_thresh)

        mask_file = osgeo.gdal.Open(
            mask_processing_source, osgeo.gdal.GA_Update)
        if mask_file is None:
            raise RuntimeError(
                f'Could not apply amplitude threshold to mask stage: '
                f'{mask_processing_source}')
        mask_arr = mask_file.ReadAsArray()
        mask_file.GetRasterBand(1).WriteArray(mask_arr * amp_file)
        mask_file.FlushCache()
        mask_file = None

    # Crop/expand mask to DEM size
    tmp_mask = f'{maskfilename}.final.tmp.tif'
    if os.path.exists(tmp_mask):
        os.remove(tmp_mask)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        ds_final_stage = osgeo.gdal.Warp(
            tmp_mask, mask_processing_source, format='GTiff',
            cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrres[0],
            yRes=arrres[1], targetAlignedPixels=True, multithread=False,
            options=['-overwrite'])
        if ds_final_stage is None:
            raise RuntimeError(
                f'Failed to create final mask staging raster: {tmp_mask}')
        ds_final_stage.FlushCache()
        ds_final_stage = None

    for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        f_to_rm = f'{maskfilename}{ext}'
        if os.path.exists(f_to_rm):
            os.remove(f_to_rm)

    _translate_mask_stage(
        tmp_mask, maskfilename, outputFormat, 'final mask')
    if os.path.exists(tmp_mask):
        os.remove(tmp_mask)
    _remove_raster_files(mask_processing_source)

    mask = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
    if mask is None:
        raise RuntimeError(f'Could not reopen final mask: {maskfilename}')
    mask.SetProjection(proj)
    mask.SetDescription(maskfilename)
    mask_array = mask.ReadAsArray()
    mask_array[mask_array != 1] = 0
    mask.GetRasterBand(1).WriteArray(mask_array)
    mask.FlushCache()
    mask = None
    _validate_raster(maskfilename, 'final normalized mask', outputFormat)

    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    ds_vrt = osgeo.gdal.Translate(
        maskfilename + '.vrt', maskfilename, options=translate_options)
    if ds_vrt is None:
        raise RuntimeError(f'Could not create final mask VRT: {maskfilename}')
    ds_vrt.FlushCache()
    ds_vrt = None

    return maskfilename
