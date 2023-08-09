#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author(s): Simran Sangha, David Bekaert, & Emre Havazli
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Extract minimum required information and files to carry-out time series analysis.
Specifically, extract unwrapped interferogram, coherence, perp baseline,
LOS file(s), and (where available) tropospheric correction layers.
"""

from ARIAtools.ARIA_global_variables import ARIA_EXTERNAL_CORRECTIONS, \
    ARIA_TROPO_MODELS, ARIA_STACK_DEFAULTS, ARIA_STACK_OUTFILES
import os
import glob
import logging
import argparse
import multiprocessing
from collections import defaultdict
from datetime import datetime
from osgeo import gdal

# Import functions
from ARIAtools import progBar
from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.mask_util import prep_mask
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.vrtmanager import resampleRaster, layerCheck, \
    get_basic_attrs, dim_check
from ARIAtools.extractProduct import merged_productbbox, prep_dem, \
    export_products, gacos_correction

gdal.UseExceptions()
# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

log = logging.getLogger(__name__)

# Import TS-related global variables


def create_parser():
    """Parser to read command line arguments."""
    parser = argparse.ArgumentParser(description='Prepare ARIA products for '
                                                 'time series processing.')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
                        required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./',
                        help='Specify directory to deposit all outputs. '
                             'Default is local directory where script is launched.')
    parser.add_argument('-gp', '--gacos_products', dest='gacos_products',
                        type=str, default=None, help='Path to director(ies) '
                                                     'or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None,
                        help='Specify layers to extract as a comma '
                             'deliminated list bounded by single quotes. '
                             'Allowed keys are: "unwrappedPhase", "coherence", '
                             '"amplitude", "bPerpendicular", "bParallel", '
                             '"incidenceAngle", "lookAngle", "azimuthAngle", '
                             '"ionosphere", "troposphereWet", '
                             '"troposphereHydrostatic", "troposphereTotal", '
                             '"solidEarthTide". '
                             'If "all" specified, then all layers are extracted. '
                             'If blank, will only extract bounding box.')
    parser.add_argument('-tm', '--tropo_models', dest='tropo_models',
                        type=str, default=None, help='Provide list of '
                                                      'weather models you wish to extract. Refer to '
                                                      'ARIA_TROPO_INTERNAL for list of supported models')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
                        default='download', help='DEM file. Default is to '
                                                 'download new DEM.')
    parser.add_argument('-p', '--projection', dest='projection',
                        default='WGS84', type=str, help='projection for DEM. '
                                                        'By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None,
                        help='Provide either valid shapefile or Lat/Lon '
                             'Bounding SNWE. -- Example : "19 20 -99.5 -98.5"')
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None,
                        help='Path to mask file or "Download". '
                             'File needs to be GDAL compatabile, contain spatial '
                             'reference information, and have invalid/valid data '
                             'represented by 0/1, respectively. '
                             'If "Download", will use GSHHS water mask. '
                             'If "NLCD", will mask classes 11, 12, 90, 95; see: '
                             'www.mrlc.gov/national-land-cover-database-nlcd-2016')
    parser.add_argument('-at', '--amp_thresh', dest='amp_thresh', default=None,
                        type=str, help='Amplitudes below this threshold '
                                       'will be masked. Specify "None" to omit '
                                       'amplitude mask. default: "None".')
    parser.add_argument('-nt', '--num_threads', dest='num_threads',
                        default='2', type=str, help='Specify number of '
                                                    'threads for multiprocessing operation in gdal. '
                                                    'By default "2". Can also specify "All" to use all '
                                                    'available threads.')
    parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',
                        type=str, default='sequential', help='Method applied to '
                                                          'stitch the unwrapped data. Allowed methods are: '
                                                          '"overlap", "2stage", and "sequential". "overlap" - '
                                                          'product overlap is minimized, "2stage" - '
                                                          'minimization is done on connected components, '
                                                          '"sequential" - sequential minimization of all '
                                                          'overlapping connected components. '
                                                          'Default is "overlap".')
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str,
                        default='VRT', help='GDAL compatible output format '
                                            '(e.g., "ENVI", "GTiff"). By default files are '
                                            'generated virtually except for "bPerpendicular", '
                                            '"bParallel", "incidenceAngle", "lookAngle", '
                                            '"azimuthAngle", "unwrappedPhase" as these require '
                                            'either DEM intersection or corrections '
                                            'to be applied')
    parser.add_argument('-croptounion', '--croptounion', action='store_true',
                        dest='croptounion', help='If turned on, IFGs cropped '
                                                 'to bounds based off of union and bbox '
                                                 '(if specified). Program defaults to crop all IFGs '
                                                 'to bounds based off of common intersection and '
                                                 'bbox (if specified).')
    parser.add_argument('-ml', '--multilooking', dest='multilooking', type=int,
                        default=None, help='Multilooking factor is an integer '
                                           'multiple of standard resolution. '
                                           'E.g. 2 = 90m*2 = 180m')
    parser.add_argument('-rr', '--rankedResampling', action='store_true',
                        dest='rankedResampling', help='If turned on, IFGs '
                                                      'resampled based off of the average of pixels in a '
                                                      'given resampling window corresponding to the '
                                                      'connected component mode (if multilooking specified).'
                                                      ' Program defaults to lanczos resampling algorithm'
                                                      ' through gdal (if multilooking specified).')
    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap',
                        type=float, default=0.0081, help='Minimum km\u00b2 '
                                                         'area of overlap of scenes wrt specified bounding box.'
                                                         ' Default 0.0081 = 0.0081km\u00b2 = area of single'
                                                         ' pixel at standard 90m resolution')
    parser.add_argument('--version', dest='version', default=None,
                        help='Specify version as str, e.g. 2_0_4 or all prods; '
                             'default: all')
    parser.add_argument('--nc_version', dest='nc_version', default='1b',
                        help='Specify netcdf version as str, '
                             'e.g. 1c or all prods;'
                             'default: 1b')
    parser.add_argument('-verbose', '--verbose', action='store_true',
                        dest='verbose', help="Toggle verbose mode on.")

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    return parser.parse_args(args=iargs)


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
                data_set = gdal.Open(b_perp)
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
                mid_time.append(datetime.strptime(j, '%Y-%m-%dT%H:%M:%S.%f'))
            min_mid_time = min(mid_time)
            max_mid_time = max(mid_time)

            # Calculate time difference between minimum start and maximum end time,
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
    from copy import deepcopy

    os.environ['GDAL_PAM_ENABLED'] = 'YES'

    # Progress bar
    prog_bar = progBar.progressBar(maxValue=len(aria_prod.products[1]),
                                   print_msg='Creating stack: ')

    # Set up single stack file
    stack_dir = os.path.join(workdir, 'stack')
    if not os.path.exists(stack_dir):
        print('Creating directory: ', stack_dir)
        os.makedirs(stack_dir)

    domain_name = deepcopy(stack_layer)
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
            print('Creating directory: ', stack_dir)
            os.makedirs(stack_dir)

    # handle individual epochs if external correction layer
    if domain_name in ARIA_EXTERNAL_CORRECTIONS or \
            domain_name in ARIA_TROPO_MODELS:
        stack_layer = f'{stack_layer}/' + 'dates'

    # Find files
    int_list = glob.glob(
        os.path.join(workdir, stack_layer, '[0-9]*[0-9].vrt'))

    print(f'Number of {stack_layer} files discovered: ', len(int_list))
    dlist = sorted(int_list)
    # get dates
    aria_dates = [os.path.basename(i).split('.vrt')[0] for i in
                  int_list]

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
            log.warning('Discrepancy between "unwrappedPhase" products '
                        '(%s files) and %s products (%s files), '
                        'rejecting scenes not common between both',
                        len(ref_dlist), domain_name, len(new_dlist))
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
    width, height, geo_trans, projection, no_data = get_basic_attrs(dlist[0])

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
            except:
                log.debug('Skipping %s; it likely exists in the %s, '
                          'but was not specified in the product list',
                          dates, os.path.dirname(data[1]))
                continue

            if orbit_direction == 'D':
                orbDir = 'DESCENDING'
            elif orbit_direction == 'A':
                orbDir = 'ASCENDING'
            else:
                print('Orbit direction not recognized')
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
        print(output_file_name, ': stack generated')

    return new_dlist


def main(inps=None):
    """Run time series prepation."""
    inps = cmd_line_parse()
    print('*****************************************************************')
    print('*** Time-series Preparation Function ***')
    print('*****************************************************************')
    # if user bbox was specified, file(s) not meeting imposed spatial
    # criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped
    # “radarmetadata info” and “data layer keys+paths” dictionaries for each
    # standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file']
    # (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile,
                                                bbox=inps.bbox,
                                                workdir=inps.workdir,
                                                num_threads=inps.num_threads,
                                                url_version=inps.version,
                                                nc_version=inps.nc_version,
                                                verbose=inps.verbose)

    # pass number of threads for gdal multiprocessing computation
    if inps.num_threads.lower() == 'all':
        print('User specified use of all %s threads for '
              'gdal multiprocessing' % (str(multiprocessing.cpu_count())))
        inps.num_threads = 'ALL_CPUS'
    print('Thread count specified for '
          'gdal multiprocessing = %s' % (inps.num_threads))

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection,
    # but user may specify union), and expected shape for DEM.
    (standardproduct_info.products[0], standardproduct_info.products[1],
     standardproduct_info.bbox_file, prods_TOTbbox,
     prods_TOTbbox_metadatalyr, arrres,
     proj) = merged_productbbox(standardproduct_info.products[0],
                                standardproduct_info.products[1],
                                os.path.join(inps.workdir,
                                             'productBoundingBox'),
                                standardproduct_info.bbox_file,
                                inps.croptounion,
                                num_threads=inps.num_threads,
                                minimumOverlap=inps.minimumOverlap,
                                verbose=inps.verbose)

    # Download/Load DEM & Lat/Lon arrays, providing bbox,
    # expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        print('Download/cropping DEM')
        # DEM parms
        dem_dict = {
        'demfilename': inps.demfile,
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'prods_TOTbbox_metadatalyr': prods_TOTbbox_metadatalyr,
        'proj': proj,
        'arrres': arrres,
        'workdir': inps.workdir,
        'outputFormat': inps.outputFormat,
        'num_threads': inps.num_threads,
        'multilooking': inps.multilooking,
        'rankedResampling': inps.rankedResampling
        }
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(**dem_dict)

    # Load or download mask (if specified).
    if inps.mask is not None:
        # Extract amplitude layers
        amplitude_products = []
        for d in standardproduct_info.products[1]:
            if 'amplitude' in d:
                for item in list(set(d['amplitude'])):
                    amplitude_products.append(item)
        # mask parms
        mask_dict = {
        'product_dict': amplitude_products,
        'maskfilename': inps.mask,
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'proj': proj,
        'amp_thresh': inps.amp_thresh,
        'arrres': arrres,
        'workdir': inps.workdir,
        'outputFormat': inps.outputFormat,
        'num_threads': inps.num_threads,
        'multilooking': inps.multilooking,
        'rankedResampling': inps.rankedResampling
        }
        inps.mask = prep_mask(**mask_dict)

    # Extract
    # aria_extract default parms
    export_dict = {
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'dem': demfile,
        'arrres': arrres,
        'lat': Latitude,
        'lon': Longitude,
        'mask': inps.mask,
        'outDir': inps.workdir,
        'outputFormat': inps.outputFormat,
        'stitchMethodType': inps.stitchMethodType,
        'verbose': inps.verbose,
        'num_threads': inps.num_threads,
        'multilooking': inps.multilooking
    }

    # export unwrappedPhase
    layers = ['unwrappedPhase', 'coherence']
    print('\nExtracting unwrapped phase, coherence, '
          'and connected components for each interferogram pair')
    ref_arr_record = export_products(standardproduct_info.products[1],
                                     tropo_total=False,
                                     layers=layers,
                                     rankedResampling=inps.rankedResampling,
                                     **export_dict)

    # Remove pairing and pass combined dictionary of all layers
    extract_dict = defaultdict(list)
    for d in standardproduct_info.products[1]:
        for key in standardproduct_info.products[1][0].keys():
            if key in d.keys():
                for item in list(set(d[key])):
                    extract_dict[key].append(item)

    layers = ['incidenceAngle', 'lookAngle', 'azimuthAngle']
    print('\nExtracting single incidence angle, look angle and azimuth angle '
          'files valid over common interferometric grid')
    prod_arr_record = export_products([extract_dict],
                                      tropo_total=False,
                                      layers=layers,
                                      **export_dict)
    # Track consistency of dimensions
    dim_check(ref_arr_record, prod_arr_record)

    layers = ['bPerpendicular']
    print('\nExtracting perpendicular baseline grids for each '
          'interferogram pair')
    prod_arr_record = export_products(standardproduct_info.products[1],
                                      tropo_total=False,
                                      layers=layers,
                                      **export_dict)
    # Track consistency of dimensions
    dim_check(ref_arr_record, prod_arr_record)

    # Extracting other layers, if specified
    layers, inps.tropo_total, \
        model_names = layerCheck(standardproduct_info.products[1],
                                 inps.layers,
                                 inps.nc_version,
                                 inps.gacos_products,
                                 inps.tropo_models,
                                 extract_or_ts='tssetup')
    if layers != [] or inps.tropo_total is True:
        if layers != []:
            print('\nExtracting optional, user-specified layers %s for each '
                  'interferogram pair' % (layers))
        if inps.tropo_total is True:
            print('\nExtracting, %s for each applicable '
                  'interferogram pair' % ('troposphereTotal'))
        prod_arr_record = export_products(standardproduct_info.products[1],
                                          tropo_total=inps.tropo_total,
                                          model_names=model_names,
                                          layers=layers,
                                          **export_dict)
        # Track consistency of dimensions
        dim_check(ref_arr_record, prod_arr_record)

    # Perform GACOS-based tropospheric corrections (if specified).
    if inps.gacos_products:
        gacos_correction(standardproduct_info.products, inps.gacos_products,
                         standardproduct_info.bbox_file, prods_TOTbbox,
                         outDir=inps.workdir, outputFormat=inps.outputFormat,
                         verbose=inps.verbose, num_threads=inps.num_threads)

    # Generate UNW stack
    ref_dlist = generate_stack(standardproduct_info, 'unwrappedPhase',
                               'unwrapStack', workdir=inps.workdir)

    # prepare additional stacks for other layers
    layers += ARIA_STACK_DEFAULTS
    layers.remove('unwrappedPhase')
    remove_lyrs = []
    for i in layers:
        lyr_dir = os.path.join(inps.workdir, i)
        if not os.path.exists(lyr_dir):
            if i in layers:
                remove_lyrs.append(i)
    layers = [i for i in layers if i not in remove_lyrs]
    if inps.tropo_total is False:
        if 'troposphereTotal' in layers:
            layers.remove('troposphereTotal')
    if inps.gacos_products:
        layers += ['gacos_corrections']

    # generate other stack layers
    # generate stack default parms
    stack_dict = {
        'workdir': inps.workdir,
        'ref_dlist': ref_dlist
    }
    for layer in layers:
        print('')
        if layer in ARIA_STACK_OUTFILES.keys():
            # iterate through model dirs if necessary
            if 'tropo' in layer:
                model_dirs = glob.glob(inps.workdir + f'/{layer}/*',
                                       recursive=True)
                model_dirs = [os.path.basename(i) for i in model_dirs]
                for sublyr in model_dirs:
                    generate_stack(standardproduct_info,
                                   sublyr,
                                   ARIA_STACK_OUTFILES[sublyr],
                                   ref_tropokey=layer,
                                   **stack_dict)
            else:
                generate_stack(standardproduct_info,
                               layer,
                               ARIA_STACK_OUTFILES[layer],
                               **stack_dict)
        else:
            msg = f'Available layers are: {ARIA_STACK_OUTFILES.keys()}'
            raise Exception(f'Selected {layer} not supported in tsSetup' + msg)
