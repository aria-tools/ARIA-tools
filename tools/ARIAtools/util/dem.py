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

import logging
import os

import dem_stitcher
import numpy as np
import osgeo

import ARIAtools.util.shp

LOGGER = logging.getLogger(__name__)


def prep_dem(
    demfilename,
    bbox_file,
    prods_TOTbbox,
    prods_TOTbbox_metadatalyr,
    proj,
    arrres=None,
    workdir="./",
    outputFormat="ENVI",
    num_threads="2",
    dem_name: str = "glo_90",
    multilooking=None,
    rankedResampling=False,
    runlog=None,
):
    """
    Function to load and export DEM, lat, lon arrays.
    If "Download" flag is specified, DEM will be downloaded on the fly.
    """
    LOGGER.debug("prep_dem")
    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir, "DEM")
    aria_dem = os.path.join(workdir, f"{dem_name}.dem")
    os.makedirs(workdir, exist_ok=True)

    # bounds of user bbox
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds

    # File must be physically extracted, cannot proceed with VRT format.
    # Defaulting to ENVI format.
    if outputFormat == "VRT":
        outputFormat = "ENVI"

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    if demfilename.lower() == "download":
        if dem_name not in dem_stitcher.datasets.DATASETS:
            raise ValueError(
                f"{dem_name} must be in " f"{', '.join(dem_stitcher.datasets.DATASETS)}"
            )

        LOGGER.info("Downloading DEM: %s", dem_name)
        demfilename = download_dem(
            aria_dem, prods_TOTbbox_metadatalyr, num_threads, dem_name, runlog
        )

    # checks for user specified DEM, ensure it's georeferenced
    else:
        LOGGER.info("Using user specified DEM %s", demfilename)
        demfilename = os.path.abspath(demfilename)
        assert os.path.exists(demfilename), f"Cannot open DEM at: {demfilename}"

        ds_u = osgeo.gdal.Open(demfilename)
        epsg = osgeo.osr.SpatialReference(wkt=ds_u.GetProjection()).GetAttrValue(
            "AUTHORITY", 1
        )
        assert epsg is not None, f"No projection information in DEM: {demfilename}"

    # write cropped DEM
    if demfilename == os.path.abspath(aria_dem):
        LOGGER.warning(
            "The DEM you specified already exists in %s, " "using the existing one...",
            os.path.dirname(aria_dem),
        )
        existing_dem = osgeo.gdal.Open(aria_dem)
        assert existing_dem is not None, f"Could not open DEM at: {aria_dem}"
        existing_dem = None

    else:
        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            gdal_warp_kwargs = {
                "format": outputFormat,
                "cutlineDSName": prods_TOTbbox,
                "outputBounds": bounds,
                "outputType": osgeo.gdal.GDT_Int16,
                "dstNodata": -32768,
                "xRes": arrres[0],
                "yRes": arrres[1],
                "targetAlignedPixels": True,
                "multithread": True,
            }
            osgeo.gdal.Warp(
                aria_dem,
                demfilename,
                options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs),
            )

        update_file = osgeo.gdal.Open(aria_dem, osgeo.gdal.GA_Update)
        update_file.SetProjection(proj)
        vrt_ds = osgeo.gdal.Translate(f"{aria_dem}.vrt", aria_dem, format="VRT")
        assert vrt_ds is not None, f"Could not build VRT for DEM at: {aria_dem}"
        vrt_ds = None
        LOGGER.info("Applied cutline to produce 3 arc-sec SRTM DEM: %s", aria_dem)

    # Load DEM and setup lat and lon arrays
    # pass expanded DEM for metadata field interpolation
    bounds = list(ARIAtools.util.shp.open_shp(prods_TOTbbox_metadatalyr).bounds)

    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        gdal_warp_kwargs = {
            "format": outputFormat,
            "outputBounds": bounds,
            "dstNodata": -32768,
            "xRes": arrres[0],
            "yRes": arrres[1],
            "targetAlignedPixels": True,
            "multithread": True,
            "options": ["-overwrite"],
        }
        demfile_expanded = aria_dem.replace(".dem", "_expanded.dem")
        ds_aria_expanded = osgeo.gdal.Warp(
            demfile_expanded,
            aria_dem,
            options=osgeo.gdal.WarpOptions(**gdal_warp_kwargs),
        )

    # Define lat/lon arrays for fullres layers
    gt = ds_aria_expanded.GetGeoTransform()
    xs, ys = ds_aria_expanded.RasterXSize, ds_aria_expanded.RasterYSize

    lat = np.linspace(gt[3], gt[3] + (gt[5] * (ys - 1)), ys)
    lat = np.repeat(lat[:, np.newaxis], xs, axis=1)
    lon = np.linspace(gt[0], gt[0] + (gt[1] * (xs - 1)), xs)
    lon = np.repeat(lon[:, np.newaxis], ys, axis=1).T

    ds_aria_expanded = None

    return aria_dem, demfile_expanded, lat, lon


def _validate_tile_integrity(tile_path):
    """
    Validate that a DEM tile file exists, has non-zero size, and is readable.

    Parameters
    ----------
    tile_path : str
        Path to the tile file to validate

    Returns
    -------
    bool
        True if tile is valid, False if needs re-download
    """
    if not os.path.exists(tile_path):
        return False

    # Check file size is non-zero
    try:
        file_size = os.path.getsize(tile_path)
        if file_size == 0:
            LOGGER.warning("Tile %s has zero size, needs re-download", tile_path)
            return False
    except OSError as e:
        LOGGER.warning("Cannot access tile %s: %s", tile_path, e)
        return False

    # Verify file is readable by GDAL
    try:
        ds = osgeo.gdal.Open(tile_path, osgeo.gdal.GA_ReadOnly)
        if ds is None:
            LOGGER.warning(
                "Tile %s cannot be opened by GDAL, needs re-download", tile_path
            )
            return False
        ds = None
        return True
    except Exception as e:
        LOGGER.warning("Error validating tile %s: %s", tile_path, e)
        return False


def download_dem(
    path_dem, path_prod_union, num_threads, dem_name="glo_90", runlog=None
):
    """Download the DEM over product bbox union."""
    LOGGER.debug("download_dem")
    root = os.path.splitext(path_dem)[0]
    vrt_path = f"{root}_uncropped.vrt"

    # Check that VRT tiles exist and are valid
    tiles_exist = False
    if os.path.exists(vrt_path):
        try:
            ds = osgeo.gdal.Open(vrt_path, osgeo.gdal.GA_ReadOnly)
            if ds is not None:
                tile_names = ds.GetFileList()
                ds = None

                # Validate each tile
                tile_checks = [
                    _validate_tile_integrity(tile_name) for tile_name in tile_names
                ]
                if False not in tile_checks:
                    tiles_exist = True
                else:
                    invalid_count = tile_checks.count(False)
                    LOGGER.warning(
                        "%d/%d tiles are invalid or missing, will re-download",
                        invalid_count,
                        len(tile_names),
                    )
        except Exception as e:
            LOGGER.warning("Error validating VRT %s: %s", vrt_path, e)

    # Retrieve update mode
    update_mode = "full_extract"
    if runlog is not None:
        log_data = runlog.load()
        if "update_mode" in log_data.keys():
            update_mode = log_data["update_mode"]

    # Check if DEM has already been downloaded and overlaps necessary area
    if tiles_exist and update_mode != "full_extract":
        LOGGER.warning("%s has already been downloaded. Skipping download.", vrt_path)

    else:
        dirname = os.path.dirname(path_dem)
        tile_dir = os.path.join(dirname, f"{dem_name}_tiles")

        # Download DEM
        prod_shapefile = ARIAtools.util.shp.open_shp(path_prod_union)
        extent = prod_shapefile.bounds

        localize_tiles_to_gtiff = False if dem_name == "glo_30" else True
        dem_tile_paths = dem_stitcher.get_dem_tile_paths(
            bounds=extent,
            dem_name=dem_name,
            localize_tiles_to_gtiff=localize_tiles_to_gtiff,
            tile_dir=tile_dir,
        )

        ds = osgeo.gdal.BuildVRT(vrt_path, dem_tile_paths)

    return vrt_path
