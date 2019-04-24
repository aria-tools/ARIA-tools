#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import sys
import numpy as np

from osgeo import gdal, ogr
gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')


def open_shapefile(fname, lyrind, ftind):
    '''
        Open a existing shapefile and pass the coordinates back.
        ##SS => see commend on projection of the shapefile, expand the desciption of what the function would do based on this.
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
    file_bbox = loads(geom.ExportToWkt())                       ##SS => We should track the projection of the shapefile adn pass it up? e.g. what it the data has different projection than that of the shapefile?

    return file_bbox



def save_shapefile(fname, polygon, drivername):
    '''
        ##SS => add descritpion and limitations. Is there not a way to add the projection pf the shape file as well?
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
