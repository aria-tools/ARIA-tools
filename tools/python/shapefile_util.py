#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import ogr
import sys
import numpy as np

from osgeo import ogr
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
plt.style.use('seaborn')

from osgeo import gdal, ogr
gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')


def open_shapefile(fname, lyrind, ftind):
    '''
        Open a existing shapefile and pass the coordinates back.
    '''

    # import dependencies
    from shapely.wkt import loads

    # opening the file
    file_bbox = ogr.Open(fname)

    #If layer name provided
    if type(lyrind) is str:
        file_bbox = file_bbox.GetLayerByName(lyrind).GetFeature(ftind)
    #If layer index provided
    else:
        file_bbox = file_bbox.GetLayerByIndex(lyrind).GetFeature(ftind)
    geom = file_bbox.GetGeometryRef()
    file_bbox = loads(geom.ExportToWkt())                       
    return file_bbox


def save_shapefile(fname, polygon, drivername):
    '''
        Save a shapefile
    '''

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



def plot_shapefile(fname):
    # Extract first layer of features from shapefile using OGR
    ds = ogr.Open(fname)
    nlay = ds.GetLayerCount()
    lyr = ds.GetLayer(0)

    # Get extent and calculate buffer size
    ext = lyr.GetExtent()
    xoff = (ext[1]-ext[0])/50
    yoff = (ext[3]-ext[2])/50

    # Prepare figure
    fig, ax = plt.subplots()
    ax.set_xlim(ext[0]-xoff,ext[1]+xoff)
    ax.set_ylim(ext[2]-yoff,ext[3]+yoff)

    paths = []
    lyr.ResetReading()

    # Read all features in layer and store as paths
    for feat in lyr:
        geom = feat.geometry()
        codes = []
        all_x = []
        all_y = []
        for i in range(geom.GetGeometryCount()):
            # Read ring geometry and create path
            r = geom.GetGeometryRef(i)
            x = [r.GetX(j) for j in range(r.GetPointCount())]
            y = [r.GetY(j) for j in range(r.GetPointCount())]
            # skip boundary between individual rings
            codes += [mpath.Path.MOVETO] + \
                         (len(x)-1)*[mpath.Path.LINETO]
            all_x += x
            all_y += y
        path = mpath.Path(np.column_stack((all_x,all_y)), codes)
        paths.append(path)

    # Add paths as patches to axes
    for path in paths:
        patch = mpatches.PathPatch(path, fill=False, facecolor='blue', edgecolor='black', linewidth=1)

        ax.add_patch(patch)

    ax.set_aspect(1.0)
    ax.grid(False)
    plt.show()
