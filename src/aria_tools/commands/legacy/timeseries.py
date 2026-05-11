#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author(s): Simran Sangha, David Bekaert, & Emre Havazli
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
ARIA-tool to run time series preparation.

Extract minimum required information and files to carry-out time series
analysis. Specifically, extract unwrapped interferogram, coherence, perp
baseline, LOS file(s), and (where available) tropospheric correction layers.
"""
import os
import logging
import argparse
import osgeo.gdal
import tile_mate

import ARIAtools.util.log
import ARIAtools.constants
from aria_tools.core.workflows import (
    run_timeseries_workflow,
)

from ARIAtools.constants import ARIA_STACK_DEFAULTS, ARIA_STACK_OUTFILES

osgeo.gdal.UseExceptions()

# Suppress warnings
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')
LOGGER = logging.getLogger('ariaTSsetup.py')


def create_parser():
    """Parser to read command line arguments."""
    parser = argparse.ArgumentParser(
        description='Prepare standard GUNW products products for '
                    'time series processing.')
    parser.add_argument(
        '-f', '--file', dest='imgfile', type=str, required=True,
        help='List of Sentinel-1 GUNW or NISAR GUNW products '
             '(wildcards supported) or txt file with product urls '
             'for virtual access without downloading. For virtual '
             'processing a local metadata cache is created on '
             'first run; subsequent runs read from the cache for '
             'faster initialization.')
    parser.add_argument(
        '-w', '--workdir', dest='workdir', default='./',
        help='Specify directory to deposit all outputs. Default is local '
             'directory where script is launched.')
    parser.add_argument(
        '-l', '--layers', dest='layers', default='standard',
        help='Specify the layers to extract as a comma-separated list enclosed '
             'in single quotes. Allowed values include: "unwrappedPhase", '
             '"coherence", "amplitude", "bPerpendicular", "bParallel", '
             '"incidenceAngle", "lookAngle", "azimuthAngle", "ionosphere", '
             '"troposphereWet", "troposphereHydrostatic", "troposphereTotal", '
             '"solidEarthTide". Use "all" to extract all interferometry and '
             'geometry layers. Correction layers must be explicitly specified '
             'to be extracted. If left blank, only the bounding box will be extracted.')
    parser.add_argument(
        '-tm', '--tropo_models', dest='tropo_models', type=str, default='all',
        help='Specify the weather model(s) to extract. '
             'The default is "all", which extracts all models in the product.')
    parser.add_argument(
        '-d', '--demfile', dest='demfile', type=str, default='download',
        help='DEM file. Default is to download new DEM.')
    parser.add_argument(
        '-p', '--projection', dest='projection', default='4326', type=str,
        help='EPSG projection code for DEM. By default 4326. '
             'Specify "native" to pass most common '
             'projection from stack.')
    parser.add_argument(
        '-b', '--bbox', dest='bbox', type=str, default=None,
        help='Provide either valid shapefile or Lat/Lon Bounding SNWE. -- '
             'Example : "19 20 -99.5 -98.5"')
    parser.add_argument(
        '-m', '--mask', dest='mask', type=str, default=None,
        help='Specify either path to valid water mask, or '
             'download using one of the following '
             f'data sources: {tile_mate.stitcher.DATASET_SHORTNAMES}')
    parser.add_argument(
        '-at', '--amp_thresh', dest='amp_thresh', default=None, type=str,
        help='Amplitudes below this threshold will be masked. Specify "None" '
             'to omit amplitude mask. default: "None".')
    parser.add_argument(
        '-nt', '--num_threads', dest='num_threads', default='2', type=str,
        help='Specify number of threads for multiprocessing operation '
             'in gdal. By default "2". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument(
        '-of', '--outputFormat', dest='outputFormat', type=str, default='VRT',
        help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By '
             'default files are generated virtually except for '
             '"bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", '
             '"azimuthAngle", "unwrappedPhase" as these require either DEM '
             'intersection or corrections to be applied')
    parser.add_argument(
        '-croptounion', '--croptounion', action='store_true',
        dest='croptounion',
        help='If turned on, IFGs cropped to bounds based off of union and '
             'bbox (if specified). Program defaults to crop all IFGs '
             'to bounds based off of common intersection and bbox (if '
             'specified).')
    parser.add_argument(
        '-ml', '--multilooking', dest='multilooking', type=int, default=None,
        help='Multilooking factor is an integer multiple of standard '
             'resolution. E.g. 2 = 90m*2 = 180m')
    parser.add_argument(
        '-rr', '--rankedResampling', action='store_true',
        dest='rankedResampling',
        help='If turned on, IFGs resampled based off of the average of pixels '
             'in a given resampling window corresponding to the connected '
             'component mode (if multilooking specified). Program defaults to '
             'lanczos resampling algorithm through gdal (if multilooking '
             'specified).')
    parser.add_argument(
        '-mo', '--minimumOverlap', dest='minimumOverlap', type=float,
        default=0.0081,
        help='Minimum km\u00b2 area of overlap of scenes wrt specified '
             'bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single'
             'pixel at standard 90m resolution')
    parser.add_argument(
        '-if', '--iono_filter', action='store_true', dest='iono_filter',
        help='Enable spatial filtering and quadratic surface approximation '
             'of the NISAR ionosphere layer. Caution: This may smooth out '
             'valid short-wavelength signals. (Note: This filter is always '
             'enforced for S1 GUNWs to mitigate large, unreliable artifacts).'
    )
    parser.add_argument(
        '--version', dest='version', default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods; default: all')
    parser.add_argument(
        '--nc_version', dest='nc_version', default='1b',
        help='Specify netcdf version as str, e.g. 1c or all prods; '
             'default: 1b')
    parser.add_argument(
        '-verbose', '--verbose', action='store_true', dest='verbose',
        help="Toggle verbose mode on.")
    parser.add_argument(
        '--log-level', 
        choices=['debug', 'info', 'warning', 'error'], 
        default='info', 
        help='Logger log level. Default: info.'
    )
    return parser

def main():
    """Run time series prepation."""
    parser = create_parser()
    args = parser.parse_args()
    args.workdir = os.path.abspath(args.workdir)

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]

    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    LOGGER.info('ARIAtools version: %s' % ARIAtools.__version__)
    print('*****************************************************************')
    LOGGER.info('*** Time-series Preparation Function ***')
    print('*****************************************************************')
    run_timeseries_workflow(
        args,
        logger=LOGGER,
        stack_defaults=ARIA_STACK_DEFAULTS,
        stack_outputs=ARIA_STACK_OUTFILES,
    )


if __name__ == '__main__':
    main()
