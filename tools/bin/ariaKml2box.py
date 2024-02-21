#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Emre Havazli & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import logging
import argparse

from osgeo import gdal

gdal.UseExceptions()

LOGGER = logging.getLogger('ariaKml2box.py')

def createParser():
    '''
    Convert .kml files of Google Earth polygons to GeoJSON files which
    can be used as bounding box.
    Output GeoJSON file is compatible with ASF Vertex product search.
    '''
    parser = argparse.ArgumentParser(
        description='Function to convert Google Earth .kml files to GeoJSON '
                    'files')
    parser.add_argument(
        '-w', '--workdir', dest='workdir', default='./',
        help='Specify directory for output file. Default is current working '
              'directory')
    parser.add_argument(
        '-f', '--file', dest='inFile', type=str, required=True,
        help='Polygon kml/kmz from Google Earth')
    parser.add_argument(
        '-o', '--outfile', dest='outFile', type=str, required=True,
        help='Output file name')
    parser.add_argument(
        '--log-level', default='warning', help='Logger log level')
    return parser

def main():
    parser = createParser()
    args = parser.parse_args()
    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    if not os.path.exists(args.workdir):
        LOGGER.info('Creating directory: %s', args.workdir)
        os.makedirs(args.workdir)
    else:
        LOGGER.info('Directory %s already exists.', args.workdir)

    srcDS = gdal.OpenEx(args.inFile)
    outfile = os.path.join(args.workdir, args.outFile)
    gdal.VectorTranslate(outfile, srcDS, format='GeoJSON', dim='XY')

if __name__ == '__main__':
    try:
        print('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except BaseException:
        pass
    main()
