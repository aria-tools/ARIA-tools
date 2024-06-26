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
import glob
import copy
import logging
import argparse
import datetime
import collections

import osgeo.gdal
import tile_mate

import ARIAtools.product
import ARIAtools.util.vrt
import ARIAtools.util.misc
import ARIAtools.util.mask
import ARIAtools.util.log
import ARIAtools.util.dem
import ARIAtools.extractProduct

from ARIAtools.constants import ARIA_EXTERNAL_CORRECTIONS, \
    ARIA_TROPO_MODELS, ARIA_STACK_DEFAULTS, ARIA_STACK_OUTFILES

osgeo.gdal.UseExceptions()

# Suppress warnings
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')
LOGGER = logging.getLogger('ariaTSsetup.py')


def create_parser():
    """Parser to read command line arguments."""
    parser = argparse.ArgumentParser(
        description='Prepare ARIA products for time series processing.')
    parser.add_argument(
        '-f', '--file', dest='imgfile', type=str, required=True,
        help='ARIA file')
    parser.add_argument(
        '-w', '--workdir', dest='workdir', default='./',
        help='Specify directory to deposit all outputs. Default is local '
             'directory where script is launched.')
    parser.add_argument(
        '-gp', '--gacos_products', dest='gacos_products', type=str,
        default=None,
        help='Path to director(ies) or tar file(s) containing GACOS products.')
    parser.add_argument(
        '-l', '--layers', dest='layers', default=None,
        help='Specify layers to extract as a comma deliminated list bounded '
             'by single quotes. Allowed keys are: "unwrappedPhase", '
             '"coherence", "amplitude", "bPerpendicular", "bParallel", '
             '"incidenceAngle", "lookAngle", "azimuthAngle", "ionosphere", '
             '"troposphereWet", "troposphereHydrostatic", "troposphereTotal", '
             '"solidEarthTide". If "all" specified, then all layers are '
             'extracted. If blank, will only extract bounding box.')
    parser.add_argument(
        '-tm', '--tropo_models', dest='tropo_models', type=str, default=None,
        help='Provide list of weather models you wish to extract. Refer to '
             'ARIA_TROPO_INTERNAL for list of supported models')
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
        '--log-level', default='warning', help='Logger log level')
    return parser


def extract_bperp_dict(domain_name, aria_prod):
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
                data_set = osgeo.gdal.Open(b_perp)
                band = data_set.GetRasterBand(1)

                # returns [min, max, mean, std]
                try:
                    # gdal~3.5
                    stat = band.GetStatistics(True, True)[2]

                except Exception as E:
                    # gdal~3.4
                    stat = band.GetStatistics(False, True)[2]
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
                   workdir='./', ref_tropokey=None, ref_dlist=None):
    """Generate time series stack."""
    os.environ['GDAL_PAM_ENABLED'] = 'YES'
    # Progress bar
    prog_bar = ARIAtools.util.misc.ProgressBar(
        maxValue=len(aria_prod.products[1]), print_msg='Creating stack: ')

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
        domain_name in ARIA_TROPO_MODELS):
        stack_layer = f'{stack_layer}/' + 'dates'

    # Find files
    int_list = glob.glob(
        os.path.join(workdir, stack_layer, '[0-9]*[0-9].vrt'))
    dlist = sorted(int_list)
    LOGGER.info(
        'Number of %s files discovered: %d' % (stack_layer, len(int_list)))

    # get dates
    aria_dates = [
        os.path.basename(i).split('.vrt')[0] for i in int_list]

    # only perform following checks if a differential layer
    b_perp = []
    new_dlist = [os.path.basename(i).split('.vrt')[0] for i in dlist]
    if domain_name not in ARIA_EXTERNAL_CORRECTIONS and \
            domain_name not in ARIA_TROPO_MODELS:

        # get az times for each date
        aztime_list = []
        for i in aria_prod.products[0]:
            aztime_list.append(i['azimuthZeroDopplerMidTime'])

        # get bperp value
        b_perp = extract_bperp_dict(domain_name, dlist)

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
    orbit_direction = str.split(os.path.basename(aria_prod.files[0]), '-')[2]

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
            <MDI key="orbitDirection">{orbDir}</MDI>'''
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

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]

    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    # pass number of threads for gdal multiprocessing computation
    if args.num_threads.lower() == 'all':
        args.num_threads = 'ALL_CPUS'

    LOGGER.info('ARIAtools version: %s' % ARIAtools.__version__)
    LOGGER.info('Time-series Preparation Function')
    LOGGER.info(
        'Thread count specified for gdal multiprocessing = %s' % (
            args.num_threads))

    # if user bbox was specified, file(s) not meeting imposed spatial
    # criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped
    # “radarmetadata info” and “data layer keys+paths” dictionaries for each
    # standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file']
    # (if bbox specified)
    LOGGER.info('Building ARIA product instance')
    standardproduct_info = ARIAtools.product.Product(
        args.imgfile, bbox=args.bbox, projection=args.projection,
        workdir=args.workdir, num_threads=args.num_threads,
        url_version=args.version, nc_version=args.nc_version,
        verbose=args.verbose)

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection,
    # but user may specify union), and expected shape for DEM.
    LOGGER.info('Extracting and merging product bounding boxes')
    (standardproduct_info.products[0], standardproduct_info.products[1],
     standardproduct_info.bbox_file, prods_TOTbbox,
     prods_TOTbbox_metadatalyr, arrres, proj, is_nisar_file) = \
        ARIAtools.extractProduct.merged_productbbox(
            standardproduct_info.products[0], standardproduct_info.products[1],
            os.path.join(args.workdir, 'productBoundingBox'),
            standardproduct_info.bbox_file, args.croptounion,
            num_threads=args.num_threads, minimumOverlap=args.minimumOverlap,
            verbose=args.verbose)

    # Download/Load DEM & Lat/Lon arrays, providing bbox,
    # expected DEM shape, and output dir as input.
    if args.demfile is not None:
        dem_dict = {
            'demfilename': args.demfile,
            'bbox_file': standardproduct_info.bbox_file,
            'prods_TOTbbox': prods_TOTbbox,
            'prods_TOTbbox_metadatalyr': prods_TOTbbox_metadatalyr,
            'proj': proj,
            'arrres': arrres,
            'workdir': args.workdir,
            'outputFormat': args.outputFormat,
            'num_threads': args.num_threads,
            'multilooking': args.multilooking,
            'rankedResampling': args.rankedResampling
        }

        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        LOGGER.info('Download/cropping DEM')
        demfile, demfile_expanded, lat, lon = \
            ARIAtools.util.dem.prep_dem(**dem_dict)
    else:
        demfile, demfile_expanded, lat, lon = None, None, None, None

    # Load or download mask (if specified).
    if args.mask is not None:

        # Extract amplitude layers
        amplitude_products = []
        for d in standardproduct_info.products[1]:
            # for NISAR GUNW
            if is_nisar_file:
                if 'coherence' in d:
                    for item in list(set(d['coherence'])):
                        amplitude_products.append(item)
            # for S1 GUNW
            else:
                if 'amplitude' in d:
                    for item in list(set(d['amplitude'])):
                        amplitude_products.append(item)

        mask_dict = {
            'product_dict': amplitude_products,
            'maskfilename': args.mask,
            'bbox_file': standardproduct_info.bbox_file,
            'prods_TOTbbox': prods_TOTbbox,
            'proj': proj,
            'amp_thresh': args.amp_thresh,
            'arrres': arrres,
            'workdir': args.workdir,
            'outputFormat': args.outputFormat,
            'num_threads': args.num_threads,
            'multilooking': args.multilooking,
            'rankedResampling': args.rankedResampling
        }
        LOGGER.info('Download/cropping mask')
        maskfilename = ARIAtools.util.mask.prep_mask(**mask_dict)
    else:
        maskfilename = None

    # Extract
    export_dict = {
        'proj': proj,
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'demfile': demfile,
        'demfile_expanded': demfile_expanded,
        'is_nisar_file': is_nisar_file,
        'arrres': arrres,
        'lat': lat,
        'lon': lon,
        'maskfile': maskfilename,
        'outDir': args.workdir,
        'outputFormat': args.outputFormat,
        'verbose': args.verbose,
        'num_threads': args.num_threads,
        'multilooking': args.multilooking
    }

    # export unwrappedPhase
    layers = ['unwrappedPhase', 'coherence']
    LOGGER.info(
        'Extracting unwrapped phase, coherence, and connected components for '
        'each interferogram pair')
    ref_arr_record = ARIAtools.extractProduct.export_products(
        standardproduct_info.products[1], tropo_total=False, layers=layers,
        rankedResampling=args.rankedResampling, multiproc_method='threads',
        **export_dict)

    # Remove pairing and pass combined dictionary of all layers
    extract_dict = collections.defaultdict(list)
    for d in standardproduct_info.products[1]:
        for key in standardproduct_info.products[1][0].keys():
            if key in d.keys():
                for item in list(set(d[key])):
                    extract_dict[key].append(item)

    layers = ['incidenceAngle', 'lookAngle', 'azimuthAngle']
    LOGGER.info(
        'Extracting single incidence angle, look angle and azimuth angle '
        'files valid over common interferometric grid')
    prod_arr_record = ARIAtools.extractProduct.export_products(
        [extract_dict], tropo_total=False, layers=layers,
        multiproc_method='gnu_parallel', **export_dict)

    # Track consistency of dimensions
    ARIAtools.util.vrt.dim_check(ref_arr_record, prod_arr_record)

    # Extracting other layers, if specified
    (layers, args.tropo_total, model_names) = ARIAtools.util.vrt.layerCheck(
        standardproduct_info.products[1], args.layers, args.nc_version,
        args.gacos_products, args.tropo_models, extract_or_ts='tssetup')

    if layers != [] or args.tropo_total is True:
        if layers != []:
            LOGGER.info(
                'Extracting optional, user-specified layers %s for each '
                'interferogram pair' % (layers))

        if args.tropo_total is True:
            LOGGER.info(
                'Extracting, %s for each applicable interferogram pair' % (
                    'troposphereTotal'))
        prod_arr_record = ARIAtools.extractProduct.export_products(
            standardproduct_info.products[1], tropo_total=args.tropo_total,
            model_names=model_names, layers=layers,
            multiproc_method='gnu_parallel', **export_dict)

        # Track consistency of dimensions
        ARIAtools.util.vrt.dim_check(ref_arr_record, prod_arr_record)

    # Perform GACOS-based tropospheric corrections (if specified).
    if args.gacos_products:
        ARIAtools.extractProduct.gacos_correction(
            standardproduct_info.products, args.gacos_products,
            standardproduct_info.bbox_file, prods_TOTbbox, outDir=args.workdir,
            outputFormat=args.outputFormat, verbose=args.verbose,
            num_threads=args.num_threads)

    # Generate UNW stack
    ref_dlist = generate_stack(
        standardproduct_info, 'unwrappedPhase', 'unwrapStack',
        workdir=args.workdir)

    # prepare additional stacks for other layers
    layers += ARIA_STACK_DEFAULTS
    layers.remove('unwrappedPhase')

    remove_lyrs = []
    for i in layers:
        lyr_dir = os.path.join(args.workdir, i)
        if not os.path.exists(lyr_dir):
            if i in layers:
                remove_lyrs.append(i)

    layers = [i for i in layers if i not in remove_lyrs]
    if args.tropo_total is False:
        if 'troposphereTotal' in layers:
            layers.remove('troposphereTotal')

    if args.gacos_products:
        layers += ['gacos_corrections']

    # generate other stack layers
    # generate stack default parms
    stack_dict = {'workdir': args.workdir, 'ref_dlist': ref_dlist}
    for layer in layers:
        if layer in ARIA_STACK_OUTFILES.keys():

            # iterate through model dirs if necessary
            if 'tropo' in layer:
                model_dirs = glob.glob(
                    args.workdir + f'/{layer}/*', recursive=True)
                model_dirs = [os.path.basename(i) for i in model_dirs]

                for sublyr in model_dirs:
                    generate_stack(
                        standardproduct_info, sublyr,
                        ARIA_STACK_OUTFILES[sublyr], ref_tropokey=layer,
                        **stack_dict)

            else:
                generate_stack(
                    standardproduct_info, layer, ARIA_STACK_OUTFILES[layer],
                    **stack_dict)

        else:
            msg = f'Available layers are: {ARIA_STACK_OUTFILES.keys()}'
            raise Exception(f'Selected {layer} not supported in tsSetup' + msg)


if __name__ == '__main__':
    main()
