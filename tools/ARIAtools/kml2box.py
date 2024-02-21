#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Emre Havazli & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
from osgeo import gdal
import logging
from ARIAtools.logger import logger

gdal.UseExceptions()

log = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def createParser():
    '''
        Convert .kml files of Google Earth polygons to GeoJSON files which can be used as bounding box.
        Output GeoJSON file is compatible with ASF Vertex product search.
    '''

    import argparse
    parser = argparse.ArgumentParser(
        description='Function to convert Google Earth .kml files to GeoJSON files')
    parser.add_argument(
        '-w',
        '--workdir',
        dest='workdir',
        default='./',
        help='Specify directory for output file. Default is local directory where script is launched.')
    parser.add_argument(
        '-f',
        '--file',
        dest='inFile',
        type=str,
        required=True,
        help='Polygon kml/kmz from Google Earth')
    parser.add_argument(
        '-o',
        '--outfile',
        dest='outFile',
        type=str,
        required=True,
        help='Output file name')

    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    return parser.parse_args(args=iargs)


def main(inps=None):
    inps = cmdLineParse()

    if not os.path.exists(inps.workdir):
        log.info('Creating directory: %s', inps.workdir)
        os.makedirs(inps.workdir)
    else:
        log.info('Directory %s already exists.', inps.workdir)

    srcDS = gdal.OpenEx(inps.inFile)
    gdal.VectorTranslate(
        os.path.join(
            inps.workdir,
            inps.outFile),
        srcDS,
        format='GeoJSON',
        dim='XY')
