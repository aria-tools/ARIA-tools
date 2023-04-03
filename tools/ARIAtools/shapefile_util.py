#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import numpy as np
import logging
log = logging.getLogger(__name__)

from osgeo import gdal, ogr
gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

def open_shapefile(fname, lyrind, ftind):
    """Open a existing shapefile and pass the coordinates back."""
    # import dependencies
    from shapely.wkt import loads

    # opening the file
    file_bbox = ogr.Open(fname)

    #If layer name provided
    if isinstance(lyrind, str):
        file_bbox = file_bbox.GetLayerByName(lyrind).GetFeature(ftind)

    #If layer index provided
    else:
        file_bbox = file_bbox.GetLayerByIndex(lyrind).GetFeature(ftind)
    geom = file_bbox.GetGeometryRef()
    file_bbox = loads(geom.ExportToWkt())

    return file_bbox

def save_shapefile(fname, polygon, drivername):
    """Save a polygon shapefile."""
    # open file
    ds = ogr.GetDriverByName(drivername).CreateDataSource(fname)
    # create layer
    layer = ds.CreateLayer('', None, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger)) #Add 1 attribute
    # Create a new feature (attribute and geometry)
    feat = ogr.Feature(layer.GetLayerDefn())
    feat.SetField('id', 0)

    # Make a geometry, from input Shapely object
    geom = ogr.CreateGeometryFromWkb(polygon.wkb)
    feat.SetGeometry(geom)
    layer.CreateFeature(feat)

    # Save/close everything
    ds = layer = feat = geom = None

    return

def shapefile_area(file_bbox, bounds = False):
    """Compute km\u00b2 area of shapefile."""
    # import dependencies
    from pyproj import Proj
    from shapely.geometry import shape

    # loop through polygons
    shape_area = 0
    # pass single polygon as list
    if file_bbox.type == 'Polygon': file_bbox = [file_bbox]
    for polyobj in file_bbox:
        # get coords
        if bounds:
            # Pass coordinates of bounds as opposed to cutline
            # Necessary for estimating DEM/mask footprints
            WSEN = polyobj.bounds
            lon = np.array([WSEN[0],WSEN[0],WSEN[2],WSEN[2],WSEN[0]])
            lat = np.array([WSEN[1],WSEN[3],WSEN[3],WSEN[1],WSEN[1]])
        else:
            lon, lat = polyobj.exterior.coords.xy

        # use equal area projection centered on/bracketing AOI
        pa = Proj("+proj=aea +lat_1={} +lat_2={} +lat_0={} +lon_0={}". \
             format(min(lat), max(lat), (max(lat)+min(lat))/2, \
             (max(lon)+min(lon))/2))
        x, y = pa(lon, lat)
        cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
        shape_area += shape(cop).area/1e6  # area in km^2

    return shape_area

def chunk_area(WSEN):
    """Chunk an area ~evenly pieces < 450000 km required by the
       SRTM server."""
    from shapely.geometry import Polygon
    max_area   = 400000 # need buffer for projections inconsistencies
    W, S, E, N = WSEN
    n = 2 # start with a 2 x 2 chunk
    area = max_area + 1
    while area > max_area:
        cols = np.linspace(W, E, n + 1)
        rows = np.linspace(S, N, n + 1)
        Wi, Si, Ei, Ni = [cols[0], rows[0], cols[1], rows[1]]
        poly = Polygon([(Wi,Ni), (Wi,Si), (Ei,Si), (Ei,Ni)])
        area = shapefile_area(poly)
        n+=1
        if n>100:
            log.error('There was a problem chunking the DEM; check input '
                      'bounds')
            os.sys.exit()
    return rows, cols

def plot_shapefile(fname):
    import matplotlib.path as mpath
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    # Extract first layer of features from shapefile using OGR
    ds = ogr.Open(fname, gdal.GA_ReadOnly)
    lyr = ds.GetLayer(0)

    # Get extent and calculate buffer size
    ext = lyr.GetExtent()
    xoff = (ext[1]-ext[0])/50
    yoff = (ext[3]-ext[2])/50

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
            codes += [mpath.Path.MOVETO] + (len(x)-1)*[mpath.Path.LINETO]
            all_x += x
            all_y += y

        path = mpath.Path(np.column_stack((all_x,all_y)), codes)
        paths.append(path)

    with plt.style.context(('seaborn')):
        fig = plt.figure(figsize=(12, 9))
        ax = fig.add_subplot(111)
        ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
        ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

        # Add paths as patches to axes
        for path in paths:
            patch = mpatches.PathPatch(path, fill=False, facecolor='blue', \
                 edgecolor='black', linewidth=1)

            ax.add_patch(patch)

        ax.set_xlabel('longitude', labelpad=15, fontsize=15)
        ax.set_ylabel('latitude', labelpad=15, fontsize=15)
        ax.set_title(os.path.basename(os.path.splitext(fname)[0]), \
             fontsize=15)
        ax.set_aspect(1.0)
        ax.grid(False)
    plt.show()
