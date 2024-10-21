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


def prep_dem(demfilename, bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr,
             proj, arrres=None, workdir='./',
             outputFormat='ENVI', num_threads='2', dem_name: str = 'glo_90',
             multilooking=None, rankedResampling=False):
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
        LOGGER.warning(
            "Cannot proceed with VRT format, using ENVI format instead")
        outputFormat = 'ENVI'

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    if demfilename.lower() == 'download':
        if dem_name not in dem_stitcher.datasets.DATASETS:
            raise ValueError(
                '%s must be in %s' % (
                    dem_name, ', '.join(dem_stitcher.datasets.DATASETS)))

        LOGGER.info("Downloading DEM...")
        demfilename = download_dem(
            aria_dem, prods_TOTbbox_metadatalyr, num_threads, dem_name)

    else:  # checks for user specified DEM, ensure it's georeferenced
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

    else:
        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            gdal_warp_kwargs = {
                'format': outputFormat, 'cutlineDSName': prods_TOTbbox,
                'outputBounds': bounds, 'outputType': osgeo.gdal.GDT_Int16,
                'xRes': arrres[0], 'yRes': arrres[1],
                'targetAlignedPixels': True, 'multithread': True}
            osgeo.gdal.Warp(
                aria_dem, demfilename,
                options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs))

        update_file = osgeo.gdal.Open(aria_dem, osgeo.gdal.GA_Update)
        update_file.SetProjection(proj)
        ds_aria = osgeo.gdal.Translate(
            f'{aria_dem}.vrt', aria_dem, format='VRT')
        LOGGER.info(
            'Applied cutline to produce 3 arc-sec SRTM DEM: %s',
            aria_dem)

    # Load DEM and setup lat and lon arrays
    # pass expanded DEM for metadata field interpolation
    bounds = list(
        ARIAtools.util.shp.open_shp(prods_TOTbbox_metadatalyr).bounds)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            'format': outputFormat, 'outputBounds': bounds, 'xRes': arrres[0],
            'yRes': arrres[1], 'targetAlignedPixels': True,
            'multithread': True, 'options': ['-overwrite']}
        demfile_expanded = aria_dem.replace('.dem', '_expanded.dem')
        ds_aria_expanded = osgeo.gdal.Warp(
            demfile_expanded, aria_dem,
            options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs))

    # Delete temporary dem-stitcher directory
    # if os.path.exists(f'{dem_name}_tiles'):
    #     shutil.rmtree(f'{dem_name}_tiles')

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
        path_dem, path_prod_union, num_threads, dem_name: str = 'glo_90'):
    """Download the DEM over product bbox union."""
    LOGGER.debug('download_dem')
    root = os.path.splitext(path_dem)[0]
    vrt_path = f'{root}_uncropped.vrt'

    # Check if DEM has already been downloaded
    if os.path.exists(os.path.abspath(vrt_path)):
        LOGGER.warning('%s has already been downloaded. Skipping download.',
                       vrt_path)

    else:
        dirname = os.path.dirname(path_dem)
        tile_dir = os.path.join(dirname, f'{dem_name}_tiles')

        # Download DEM
        prod_shapefile = ARIAtools.util.shp.open_shp(path_prod_union)
        extent = prod_shapefile.bounds

        localize_tiles_to_gtiff = False if dem_name == 'glo_30' else True
        dem_tile_paths = dem_stitcher.get_dem_tile_paths(
            bounds=extent, dem_name=dem_name,
            localize_tiles_to_gtiff=localize_tiles_to_gtiff,
            tile_dir=tile_dir)

        ds = osgeo.gdal.BuildVRT(vrt_path, dem_tile_paths)

    return vrt_path
