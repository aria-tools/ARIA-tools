#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Emre Havazli & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
from osgeo import gdal, ogr
gdal.UseExceptions()

def createParser():
    '''
        Convert .kml files of Google Earth polygons to GeoJSON files which can be used as bounding box.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Function to convert Google Earth .kml files to GeoJSON files')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory for output file. Default is local directory where script is launched.')
    parser.add_argument('-f', '--file', dest='inFile', type=str, required=True, help='Google Earth kml polygon')
    parser.add_argument('-o', '--outfile', dest='outFile', type=str, required=True, help='Output file name')

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def main(inps=None):
    inps = cmdLineParse()

    if not os.path.exists(inps.workdir):
        print('Creating directory: {0}'.format(inps.workdir))
        os.makedirs(inps.workdir)
    else:
        print('Directory {0} already exists.'.format(inps.workdir))

    srcDS = gdal.OpenEx(inps.inFile)
    ds = gdal.VectorTranslate(os.path.join(inps.workdir,inps.outFile), srcDS, format='GeoJSON')
