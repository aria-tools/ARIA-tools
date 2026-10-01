# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Digital Elevation Model utilities
"""
import os
import shutil
import logging
import numpy as np

import osgeo
import dem_stitcher

import ARIAtools.util.shp

LOGGER = logging.getLogger(__name__)


def _remove_raster_files(path):
    """Remove a raster and the sidecars used by supported GDAL drivers."""
    for suffix in ['', '.vrt', '.aux.xml', '.xml', '.hdr']:
        candidate = path + suffix
        if os.path.isfile(candidate):
            os.remove(candidate)


def _validate_raster(path, context, expected_driver=None):
    """Ensure a raster is complete enough to safely use downstream."""
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


def _translate_staged_raster(stage_path, output_path, output_format, context):
    """Translate a completed GeoTIFF stage to the requested final driver."""
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


def prep_dem(demfilename, bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr,
             proj, arrres=None, workdir='./',
             outputFormat='ENVI', num_threads='2', dem_name: str = 'glo_90',
             multilooking=None, rankedResampling=False, runlog=None):
    """
    Function to load and export DEM, lat, lon arrays.
    If "Download" flag is specified, DEM will be downloaded on the fly.
    """
    LOGGER.debug('prep_dem')
    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir, 'DEM')
    aria_dem = os.path.join(workdir, f'{dem_name}.dem')
    os.makedirs(workdir, exist_ok=True)

    # bounds of user bbox
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds

    # File must be physically extracted, cannot proceed with VRT format.
    # Defaulting to ENVI format.
    if outputFormat == 'VRT':
        outputFormat = 'ENVI'

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    if demfilename.lower() == 'download':
        if dem_name not in dem_stitcher.datasets.DATASETS:
            raise ValueError(
                '%s must be in %s' % (
                    dem_name, ', '.join(dem_stitcher.datasets.DATASETS)))

        LOGGER.info('Downloading DEM: %s', dem_name)
        demfilename = download_dem(
            aria_dem, prods_TOTbbox_metadatalyr, num_threads, dem_name, runlog)

    # checks for user specified DEM, ensure it's georeferenced
    else:
        LOGGER.info("Using user specified DEM %s" % demfilename)
        demfilename = os.path.abspath(demfilename)
        assert os.path.exists(demfilename), (
            f'Cannot open DEM at: {demfilename}')

        ds_u = osgeo.gdal.Open(demfilename)
        epsg = osgeo.osr.SpatialReference(
            wkt=ds_u.GetProjection()).GetAttrValue('AUTHORITY', 1)
        assert epsg is not None, (
            f'No projection information in DEM: {demfilename}')

    # write cropped DEM
    if demfilename == os.path.abspath(aria_dem):
        LOGGER.warning('The DEM you specified already exists in %s, '
                       'using the existing one...', os.path.dirname(aria_dem))
        ds_aria = osgeo.gdal.Open(aria_dem)
        dem_processing_stage = aria_dem + '.processing.tmp.tif'
        _remove_raster_files(dem_processing_stage)
        ds_processing = osgeo.gdal.Translate(
            dem_processing_stage, ds_aria, format='GTiff')
        if ds_processing is None:
            raise RuntimeError(
                f'Failed to create DEM processing stage: '
                f'{dem_processing_stage}')
        ds_processing.FlushCache()
        ds_processing = None
        _validate_raster(
            dem_processing_stage, 'DEM processing stage', 'GTiff')
        ds_aria = None

    else:
        # Clean up pre-existing output files
        for ext in ['', '.vrt', '.aux.xml', '.xml', '.hdr', '.tmp.tif']:
            f_to_rm = f'{aria_dem}{ext}'
            if os.path.exists(f_to_rm):
                os.remove(f_to_rm)
    
        tmp_tif = f'{aria_dem}.tmp.tif'
    
        # Warp to intermediate GeoTIFF to avoid ENVI block I/O read errors
        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            gdal_warp_kwargs = {
                'format': 'GTiff', 'cutlineDSName': prods_TOTbbox,
                'outputBounds': bounds, 'outputType': osgeo.gdal.GDT_Int16,
                'dstNodata': -32768,
                'xRes': arrres[0], 'yRes': arrres[1],
                'targetAlignedPixels': True, 'multithread': False,
                'options': ['-overwrite']}
            ds_stage = osgeo.gdal.Warp(
                tmp_tif, demfilename,
                options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs))
            if ds_stage is None:
                raise RuntimeError(
                    f'Failed to create cropped DEM staging raster: {tmp_tif}')
            ds_stage.FlushCache()
            ds_stage = None
    
        # Convert the completed intermediate to the requested physical format.
        _translate_staged_raster(
            tmp_tif, aria_dem, outputFormat, 'cropped DEM')
        dem_processing_stage = tmp_tif
    
        update_file = osgeo.gdal.Open(aria_dem, osgeo.gdal.GA_Update)
        if update_file is None:
            raise RuntimeError(f'Could not update cropped DEM: {aria_dem}')
        update_file.SetProjection(proj)
        update_file.FlushCache()
        update_file = None
        ds_aria = osgeo.gdal.Translate(
            f'{aria_dem}.vrt', aria_dem, format='VRT')
        if ds_aria is None:
            raise RuntimeError(f'Could not create DEM VRT: {aria_dem}.vrt')
        ds_aria.FlushCache()
        ds_aria = None
        LOGGER.info(
            'Applied cutline to produce 3 arc-sec SRTM DEM: %s', aria_dem)

    # Load DEM and setup lat and lon arrays
    # pass expanded DEM for metadata field interpolation
    bounds = list(
        ARIAtools.util.shp.open_shp(prods_TOTbbox_metadatalyr).bounds)

    demfile_expanded = aria_dem.replace('.dem', '_expanded.dem')
    expanded_stage = demfile_expanded + '.tmp.tif'
    _remove_raster_files(demfile_expanded)
    _remove_raster_files(expanded_stage)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': 'GTiff', 'outputBounds': bounds,
            'dstNodata': -32768,
            'xRes': arrres[0],
            'yRes': arrres[1], 'targetAlignedPixels': True,
            'multithread': False, 'options': ['-overwrite'],
            'creationOptions': [
                'TILED=YES', 'COMPRESS=DEFLATE',
                'PREDICTOR=2', 'BIGTIFF=IF_SAFER']}
        ds_expanded_stage = osgeo.gdal.Warp(
            expanded_stage, dem_processing_stage,
            options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs))
        if ds_expanded_stage is None:
            raise RuntimeError(
                f'Failed to create expanded DEM staging raster: '
                f'{expanded_stage}')
        ds_expanded_stage.FlushCache()
        ds_expanded_stage = None

    _translate_staged_raster(
        expanded_stage, demfile_expanded, outputFormat, 'expanded DEM')
    _remove_raster_files(expanded_stage)
    _remove_raster_files(dem_processing_stage)
    ds_aria_expanded = osgeo.gdal.Open(
        demfile_expanded, osgeo.gdal.GA_ReadOnly)
    if ds_aria_expanded is None:
        raise RuntimeError(f'Could not reopen expanded DEM: {demfile_expanded}')

    # Define lat/lon arrays for fullres layers
    gt = ds_aria_expanded.GetGeoTransform()
    xs, ys = ds_aria_expanded.RasterXSize, ds_aria_expanded.RasterYSize

    lat = np.linspace(gt[3], gt[3] + (gt[5] * (ys - 1)), ys)
    lat = np.repeat(lat[:, np.newaxis], xs, axis=1)
    lon = np.linspace(gt[0], gt[0] + (gt[1] * (xs - 1)), xs)
    lon = np.repeat(lon[:, np.newaxis], ys, axis=1).T

    ds_aria_expanded = None

    return aria_dem, demfile_expanded, lat, lon


def download_dem(
        path_dem, path_prod_union, num_threads, dem_name='glo_90',
        runlog=None):
    """Download the DEM over product bbox union."""
    LOGGER.debug('download_dem')
    root = os.path.splitext(path_dem)[0]
    vrt_path = f"{root}_uncropped.vrt"

    # Check that VRT tiles exist
    tiles_exist = False
    if os.path.exists(vrt_path):
        ds = osgeo.gdal.Open(vrt_path, osgeo.gdal.GA_ReadOnly)
        tile_names = ds.GetFileList()
        tile_checks = [os.path.exists(tile_name) and os.path.getsize(tile_name) > 0 for tile_name in tile_names]
        if False not in tile_checks:
            tiles_exist = True
        ds = None  # Close dataset handle after reading tile file list

    # Retrieve update mode
    update_mode = 'full_extract'
    if runlog is not None:
        log_data = runlog.load()
        if 'update_mode' in log_data.keys():
            update_mode = log_data['update_mode']

    # Check if DEM has already been downloaded and overlaps necessary area
    if tiles_exist and update_mode != 'full_extract':
        LOGGER.warning(
            '%s has already been downloaded. Skipping download.', vrt_path)

    else:
        dirname = os.path.dirname(path_dem)
        tile_dir = os.path.join(dirname, f"{dem_name}_tiles")

        # Download DEM
        prod_shapefile = ARIAtools.util.shp.open_shp(path_prod_union)
        extent = prod_shapefile.bounds

        localize_tiles_to_gtiff = False if dem_name == 'glo_30' else True
        dem_tile_paths = dem_stitcher.get_dem_tile_paths(
            bounds=extent, dem_name=dem_name,
            localize_tiles_to_gtiff=localize_tiles_to_gtiff,
            tile_dir=tile_dir)

        ds = osgeo.gdal.BuildVRT(vrt_path, dem_tile_paths)
        ds = None  # Flush VRT header and close dataset handle

    return vrt_path
