# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import logging

import numpy as np
import pyproj
import shapely
import shapely.wkt
import osgeo.gdal

LOGGER = logging.getLogger(__name__)

osgeo.gdal.UseExceptions()
# Suppress warnings
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')


def open_shp(fname, lyrind=0, ftind=0):
    """Open a existing shapefile and pass the coordinates back."""
    # opening the file
    file_bbox = osgeo.ogr.Open(fname)

    # If layer name provided
    if isinstance(lyrind, str):
        file_bbox = file_bbox.GetLayerByName(lyrind).GetFeature(ftind)

    # If layer index provided
    else:
        file_bbox = file_bbox.GetLayerByIndex(lyrind).GetFeature(ftind)
    geom = file_bbox.GetGeometryRef()
    file_bbox = shapely.wkt.loads(geom.ExportToWkt())
    return file_bbox


def save_shp(fname, polygon, projection, drivername='GeoJSON'):
    """Save a polygon shapefile."""
    # Define the target spatial reference system
    target_srs = osgeo.ogr.osr.SpatialReference()
    target_srs.ImportFromEPSG(projection)

    # open file
    ds = osgeo.ogr.GetDriverByName(drivername).CreateDataSource(fname)

    # create layer
    layer = ds.CreateLayer('', target_srs, osgeo.ogr.wkbPolygon)

    # Add 1 attribute
    layer.CreateField(osgeo.ogr.FieldDefn('id', osgeo.ogr.OFTInteger))

    # Create a new feature (attribute and geometry)
    feat = osgeo.ogr.Feature(layer.GetLayerDefn())
    feat.SetField('id', 0)

    # Make a geometry, from input Shapely object
    geom = osgeo.ogr.CreateGeometryFromWkb(polygon.wkb)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)
    return


def shp_area(file_bbox, projection, bounds=False):
    """Compute km\u00b2 area of shapefile."""
    # loop through polygons
    shape_area = 0

    # pass single polygon as list
    if file_bbox.geom_type == 'Polygon':
        file_bbox = [file_bbox]

    for polyobj in file_bbox:

        # need to reproject if not already in UTM
        if projection == 4326:

            # get coords
            if bounds:
                # Pass coordinates of bounds as opposed to cutline
                # Necessary for estimating DEM/mask footprints
                WSEN = polyobj.bounds
                lon = np.array([WSEN[0], WSEN[0], WSEN[2], WSEN[2], WSEN[0]])
                lat = np.array([WSEN[1], WSEN[3], WSEN[3], WSEN[1], WSEN[1]])
            else:
                lon, lat = polyobj.exterior.coords.xy

            # use equal area projection centered on/bracketing AOI
            pa = pyproj.Proj(
                "+proj=aea +lat_1={} +lat_2={} +lat_0={} +lon_0={}".format(
                    min(lat), max(lat), (max(lat) + min(lat)) / 2,
                    (max(lon) + min(lon)) / 2))
            x, y = pa(lon, lat)
            cop = {"type": "Polygon", "coordinates": [zip(x, y)]}

            # area in km^2
            shape_area += shapely.geometry.shape(cop).area / 1e6

        else:
            # area in km^2
            shape_area = polyobj.area / 1e6
    return shape_area


def chunk_area(WSEN):
    """
    Chunk an area ~evenly pieces < 450000 km required by the SRTM server.
    """
    # need buffer for projections inconsistencies
    MAX_AREA = 400000
    W, S, E, N = WSEN

    # start with a 2 x 2 chunk
    n = 2
    area = MAX_AREA + 1
    while area > MAX_AREA:
        cols = np.linspace(W, E, n + 1)
        rows = np.linspace(S, N, n + 1)
        Wi, Si, Ei, Ni = [cols[0], rows[0], cols[1], rows[1]]
        poly = shapely.geometry.Polygon(
            [(Wi, Ni), (Wi, Si), (Ei, Si), (Ei, Ni)])
        area = shp_area(poly)
        n += 1
        if n > 100:
            LOGGER.error(
                'There was a problem chunking the DEM; check input bounds')
            raise Exception(
                "There was a problem chunking the DEM; check input bounds")
    return rows, cols


def plot_shp(fname):
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    # Extract first layer of features from shapefile using OGR
    ds = osgeo.ogr.Open(fname)
    lyr = ds.GetLayer(0)

    # Get extent and calculate buffer size
    ext = lyr.GetExtent()
    xoff = (ext[1] - ext[0]) / 50
    yoff = (ext[3] - ext[2]) / 50

    paths = []
    lyr.ResetReading()

    # Read all features in layer and store as paths
    for i, feat in enumerate(lyr):
        geom = feat.geometry()
        geom_name = geom.GetGeometryName()
        codes = []
        all_x = []
        all_y = []
        for i in range(geom.GetGeometryCount()):
            # Read ring geometry and create path
            r = geom.GetGeometryRef(i)
            if geom_name == 'MULTIPOLYGON':
                r = geom.GetGeometryRef(i)
                for j in range(r.GetGeometryCount()):
                    p = r.GetGeometryRef(j)
                    r = p
            x = [r.GetX(j) for j in range(r.GetPointCount())]
            y = [r.GetY(j) for j in range(r.GetPointCount())]

            # skip boundary between individual rings
            codes += [mpath.Path.MOVETO] + (len(x) - 1) * [mpath.Path.LINETO]
            all_x += x
            all_y += y

        path = mpath.Path(np.column_stack((all_x, all_y)), codes)
        paths.append(path)

    plt_style = 'classic'
    for test_style in ['seaborn-v0_8', 'seaborn']:
        if test_style in plt.style.available:
            plt_style = test_style

    with plt.style.context(plt_style):
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111)
        ax.set_xlim(ext[0] - xoff, ext[1] + xoff)
        ax.set_ylim(ext[2] - yoff, ext[3] + yoff)

        # Add paths as patches to axes
        for path in paths:
            patch = mpatches.PathPatch(
                path, fill=False, facecolor='blue', edgecolor='black',
                linewidth=1)
            ax.add_patch(patch)

        ax.set_xlabel('longitude', labelpad=15, fontsize=15)
        ax.set_ylabel('latitude', labelpad=15, fontsize=15)
        ax.set_title(
            os.path.basename(os.path.splitext(fname)[0]), fontsize=15)
        ax.set_aspect(1.0)
        ax.grid(False)
    plt.show()
