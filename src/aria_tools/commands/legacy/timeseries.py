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
import h5py
import glob
import copy
import logging
import argparse
import datetime
import contextlib
import json

import numpy as np
import osgeo.gdal
import tile_mate

import ARIAtools.extractProduct
import ARIAtools.util.log
import ARIAtools.util.misc
import ARIAtools.constants
from aria_tools.core.workflows import (
    run_timeseries_workflow,
)

from ARIAtools.constants import ARIA_EXTERNAL_CORRECTIONS, \
    ARIA_TROPO_MODELS, ARIA_STACK_DEFAULTS, ARIA_STACK_OUTFILES, \
    ARIA_STANDARD_LAYERS, ARIA_LAYERS

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


def extract_bperp_dict_ts(domain_name, aria_prod):
    """Extract mean bperp from products."""
    os.environ['GDAL_PAM_ENABLED'] = 'NO'
    meta = {}
    for i in aria_prod:
        pair_name = i[-21:-4]

        # only get stats for unw file
        # otherwise pass a dummy value
        stat = 0
        if domain_name == 'unwrappedPhase':

            # find corresponding bPerp file
            b_perp = i.split('/')
            b_perp[-2] = 'bPerpendicular'
            b_perp = '/'.join(b_perp)

            if os.path.exists(b_perp):
                data_set = None
                try:
                    data_set = osgeo.gdal.Open(
                        b_perp, osgeo.gdal.GA_ReadOnly
                    )
                    if data_set is not None:
                        band = data_set.GetRasterBand(1)

                        # returns [min, max, mean, std]
                        try:
                            # gdal~3.5
                            stat = band.GetStatistics(True, True)[2]
                        except Exception:
                            # gdal~3.4
                            stat = band.GetStatistics(False, True)[2]
                        
                        # Release band reference
                        band = None
                finally:
                    # CRITICAL: Ensure file close happens no matter what
                    data_set = None

        meta[pair_name] = stat

    return meta


def extract_utc_time(aria_dates, aztime_list):
    """Extract UTC time from products."""
    utc_dict = {}
    utc_time = None
    for i in range(0, len(aria_dates)):
        # Grab pair name of the product
        pair_name = aria_dates[i]

        # Only iterate on utc_calculation if values in list are different
        if ([aztime_list[0]] * len(aztime_list) != aztime_list) or \
                utc_time is None:
            # Grab mid-times, append to a list and find minimum and
            # maximum mid-times
            mid_time_list = aztime_list[i]
            mid_time = []
            for j in mid_time_list:
                mid_time.append(
                    datetime.datetime.strptime(j, '%Y-%m-%dT%H:%M:%S.%f'))
            min_mid_time = min(mid_time)
            max_mid_time = max(mid_time)

            # Calculate time difference
            # between minimum start and maximum end time,
            # and add it to mean start time.
            time_delta = (max_mid_time - min_mid_time) / 2
            utc_time = (min_mid_time + time_delta).time()

        # Write calculated UTC time into a dictionary
        # with associated pair names as keys
        utc_dict[pair_name] = utc_time.strftime("%H:%M:%S.%f")
    return utc_dict


def generate_stack(aria_prod, stack_layer, output_file_name,
                   workdir='./', ref_tropokey=None, ref_dlist=None,
                   is_nisar_file=False):
    """Generate time series stack."""
    os.environ['GDAL_PAM_ENABLED'] = 'YES'
    # Set up single stack file
    stack_dir = os.path.join(workdir, 'stack')
    if not os.path.exists(stack_dir):
        LOGGER.info('Creating directory: %s' % stack_dir)
        os.makedirs(stack_dir)

    domain_name = copy.deepcopy(stack_layer)

    # Datatypes -- all layers are Float32 except ConnComponents
    if domain_name == 'connectedComponents':
        data_type = "Int16"
    else:
        data_type = "Float32"

    # make sure to search subdirectory for specific tropo models if necessary
    if domain_name in ARIA_TROPO_MODELS:
        stack_layer = f'{ref_tropokey}/' + stack_layer
        stack_dir = os.path.join(stack_dir, ref_tropokey)
        if not os.path.exists(stack_dir):
            LOGGER.info('Creating directory: %s' % stack_dir)
            os.makedirs(stack_dir)

    # handle individual epochs if external correction layer
    if (domain_name in ARIA_EXTERNAL_CORRECTIONS or
        domain_name in ARIA_TROPO_MODELS) and not is_nisar_file:
        stack_layer = f'{stack_layer}/' + 'dates'

    # get dates
    aria_dates = \
        sorted([prod['pair_name'][0] for prod in aria_prod.products[0]])

    # data are extracted as dates for tropo and SET layers for ARIA-S1-GUNW
    # no NISAR layers are extracted this way
    if not is_nisar_file and (domain_name in ARIA_EXTERNAL_CORRECTIONS
                              or domain_name in ARIA_TROPO_MODELS):
        aria_indiv_dates = []
        rejected_dates = []
        for aria_date in aria_dates:
            dates = aria_date.split('_')
            # check reference date
            dt1_fname = os.path.join(workdir, stack_layer, dates[0] + '.vrt')
            if os.path.exists(dt1_fname):
                aria_indiv_dates += [dates[0]]
            else:
                rejected_dates += [dates[0]]

            # check secondary date
            dt2_fname = os.path.join(workdir, stack_layer, dates[1] + '.vrt')
            if os.path.exists(dt2_fname):
                aria_indiv_dates += [dates[1]]
            else:
                rejected_dates += [dates[1]]

        aria_dates = sorted(list(set(aria_indiv_dates)))
        rejected_dates = sorted(list(set(rejected_dates)))

        # report rejected dates
        if rejected_dates != []:
            LOGGER.warning(
                'The following %d date(s) lack %s layers: %s',
                 len(rejected_dates), domain_name, ", ".join(rejected_dates)
            )

    # Find files
    int_list = [os.path.join(workdir, stack_layer, aria_date + '.vrt')
                for aria_date in aria_dates]
    dlist = sorted(int_list)
    LOGGER.info(
        'Number of %s files discovered: %d' % (stack_layer, len(int_list)))

    # Progress bar
    prog_bar = ARIAtools.util.misc.ProgressBar(
        maxValue=len(int_list), prefix=f'Exporting {output_file_name}: ')

    # only perform following checks if a differential layer
    # all NISAR layers are differential
    b_perp = []
    new_dlist = [os.path.basename(i).split('.vrt')[0] for i in dlist]
    if is_nisar_file or (
            domain_name not in ARIA_EXTERNAL_CORRECTIONS
            and domain_name not in ARIA_TROPO_MODELS):

        # get az times for each date
        aztime_list = []
        for i in aria_prod.products[0]:
            aztime_list.append(i['azimuthZeroDopplerMidTime'])

        # get bperp value
        b_perp_json_file = os.path.join(
            workdir, 'bPerpendicular', 'bperp.json')

        if os.path.isfile(b_perp_json_file):
            LOGGER.debug("Loading bPerpendicular from bperp.json")
            with open(b_perp_json_file) as ifp:
                b_perp = json.loads(ifp.read())
        else:
            b_perp = extract_bperp_dict_ts(domain_name, dlist)

        # Confirm 1-to-1 match between UNW and other derived products
        if ref_dlist and new_dlist != ref_dlist:
            LOGGER.warning(
                'Discrepancy between "unwrappedPhase" products (%s files) and '
                '%s products (%s files), rejecting scenes not common between '
                'both', len(ref_dlist), domain_name, len(new_dlist))

            # subset to match other datasets
            subset_ind = []
            for i in enumerate(new_dlist):
                if i[1] in ref_dlist:
                    subset_ind.append(i[0])
            new_dlist = [new_dlist[i] for i in subset_ind]
            dlist = [dlist[i] for i in subset_ind]

    else:
        # get az times for each date
        aztime_list = len(aria_dates) * \
            [aria_prod.products[0][0]['azimuthZeroDopplerMidTime']]

    # get UTC times
    utc_time = extract_utc_time(aria_dates, aztime_list)

    # get attributes from first product
    width, height, geo_trans, projection, no_data = \
        ARIAtools.util.vrt.get_basic_attrs(dlist[0])

    # setting up a subset of the stack
    ymin, ymax, xmin, xmax = [0, height, 0, width]

    xsize = xmax - xmin
    ysize = ymax - ymin

    # extraction of radar meta-data
    wvl = aria_prod.products[0][0]['wavelength'][0]
    start_range = aria_prod.products[0][0]['slantRangeStart'][0]
    end_range = aria_prod.products[0][0]['slantRangeEnd'][0]
    range_spacing = aria_prod.products[0][0]['slantRangeSpacing'][0]
    if is_nisar_file:
        orbit_direction = str.split(os.path.basename(aria_prod.files[0]), '_')[6]
        platform = 'NISAR'
    else:
        orbit_direction = str.split(os.path.basename(aria_prod.files[0]), '-')[2]
        platform = 'Sen'

    with open(os.path.join(stack_dir, output_file_name + '.vrt'), 'w') as fid:
        fid.write('''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
        <SRS>{proj}</SRS>
        <GeoTransform>{GT0},{GT1},{GT2},{GT3},{GT4},{GT5}</GeoTransform>\n
        '''.format(xsize=xsize, ysize=ysize, proj=projection, GT0=geo_trans[0],
                   GT1=geo_trans[1], GT2=geo_trans[2], GT3=geo_trans[3],
                   GT4=geo_trans[4], GT5=geo_trans[5]))

        for data in enumerate(dlist):
            didx = data[0] + 1
            dates = data[1].split('/')[-1][:-4]
            path = None
            # Update progress bar
            prog_bar.update(didx, suffix=dates)

            try:
                acq = utc_time[dates]
            except BaseException:
                LOGGER.debug(
                    'Skipping %s; it likely exists in the %s, but was not '
                    'specified in the product list', dates,
                    os.path.dirname(data[1]))
                continue

            if orbit_direction == 'D':
                orbDir = 'DESCENDING'
            elif orbit_direction == 'A':
                orbDir = 'ASCENDING'
            else:
                LOGGER.warninig('Orbit direction not recognized')
                orbDir = 'UNKNOWN'

            path = os.path.relpath(os.path.abspath(data[1]), start=stack_dir)
            outstr = f'''  <VRTRasterBand dataType="{data_type}" band="{didx}">
        <NoDataValue>{no_data}</NoDataValue>
        <SimpleSource>
            <SourceFilename relativeToVRT="1">{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}"
                DataType="{data_type}"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}"
                ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain='{domain_name}'>
            <MDI key="Dates">{dates}</MDI>
            <MDI key="Wavelength (m)">{wvl}</MDI>
            <MDI key="UTCTime (HH:MM:SS.ss)">{acq}</MDI>
            <MDI key="startRange">{start_range}</MDI>
            <MDI key="endRange">{end_range}</MDI>
            <MDI key="slantRangeSpacing">{range_spacing}</MDI>
            <MDI key="orbitDirection">{orbDir}</MDI>
            <MDI key="PLATFORM">{platform}</MDI>'''
            fid.write(outstr)
            if b_perp != []:
                bPerp = b_perp[dates]
                outstr = f'''
            <MDI key="perpendicularBaseline">{bPerp}</MDI>'''
                fid.write(outstr)
            outstr = f'''
        </Metadata>
    </VRTRasterBand>\n'''
            fid.write(outstr)
        fid.write('</VRTDataset>\n')
        prog_bar.close()
        LOGGER.info('%s stack generated' % output_file_name)

    return new_dlist


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
        generate_stack_func=generate_stack,
        stack_defaults=ARIA_STACK_DEFAULTS,
        stack_outputs=ARIA_STACK_OUTFILES,
    )


if __name__ == '__main__':
    main()
