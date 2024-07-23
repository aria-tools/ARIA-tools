#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import argparse
import logging

import tile_mate

import ARIAtools.extractProduct
import ARIAtools.util.vrt
import ARIAtools.util.dem
import ARIAtools.util.log
import ARIAtools.util.mask
import ARIAtools.product

LOGGER = logging.getLogger('ariaExtract.py')


def createParser():
    """
    Extract specified product layers. The default will export all layers.
    """
    parser = argparse.ArgumentParser(
        description='Program to extract data and meta-data layers from ARIA '
                    'standard GUNW products. Program will handle cropping/'
                    'stitching when needed. By default, the program will crop '
                    'all IFGs to bounds determined by the common intersection '
                    'and bbox (if specified)')

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
        help='Provide list ofweather models you wish to extract. Refer to '
             'ARIA_TROPO_MODELS for list of supported models')
    parser.add_argument(
        '-d', '--demfile', dest='demfile', type=str, default=None,
        help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument(
        '-p', '--projection', dest='projection', default='4326', type=str,
        help='EPSG projection code for DEM. By default 4326. '
             'Specify "native" to pass most common '
             'projection from stack.')
    parser.add_argument(
        '-b', '--bbox', dest='bbox', type=str, default=None,
        help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- "
             "Example : '19 20 -99.5 -98.5'")
    parser.add_argument(
        '-m', '--mask', dest='mask', type=str, default=None,
        help='Specify either path to valid water mask, or '
             'download using one of the following '
             f'data sources: {tile_mate.stitcher.DATASET_SHORTNAMES}')
    parser.add_argument(
        '-at', '--amp_thresh', dest='amp_thresh', default=None, type=str,
        help='Amplitude threshold below which to mask. Specify "None" to not '
             'use amplitude mask. By default "None".')
    parser.add_argument(
        '-nt', '--num_threads', dest='num_threads', default='2', type=str,
        help='Specify number of threads for multiprocessing operation in '
             'gdal. By default "2". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument(
        '-of', '--outputFormat', dest='outputFormat', type=str, default='VRT',
        help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By '
             'default files are generated virtually except for '
             '"bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", '
             '"azimuthAngle", "unwrappedPhase" as these are require either '
             'DEM intersection or corrections to be applied')
    parser.add_argument(
        '-croptounion', '--croptounion', action='store_true',
        dest='croptounion',
        help="If turned on, IFGs cropped to bounds based off of union and "
             "bbox (if specified). Program defaults to crop all IFGs to "
             "bounds based off of common intersection and bbox (if "
             "specified).")
    parser.add_argument(
        '-ml', '--multilooking', dest='multilooking', type=int, default=None,
        help='Multilooking factor is an integer multiple of standard '
             'resolution. E.g. 2 = 90m*2 = 180m')
    parser.add_argument(
        '-rr', '--rankedResampling', action='store_true',
        dest='rankedResampling',
        help="If turned on, IFGs resampled based off of the average of pixels "
             "in a given resampling window corresponding to the connected "
             "component mode (if multilooking specified). Program defaults "
             "to lanczos resampling algorithm through gdal (if multilooking "
             "specified).")
    parser.add_argument(
        '-mo', '--minimumOverlap', dest='minimumOverlap', type=float,
        default=0.0081,
        help='Minimum km\u00b2 area of overlap of scenes wrt specified '
             'bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single '
             'pixel at standard 90m resolution')
    parser.add_argument(
        '--version', dest='version', default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods; default: all')
    parser.add_argument(
        '--nc_version', dest='nc_version', default='1b',
        help='Specify netcdf version as str, '
             'e.g. 1c or all prods; default: 1b')
    parser.add_argument(
        '-verbose', '--verbose', action='store_true', dest='verbose',
        help="Toggle verbose mode on.")
    parser.add_argument(
        '--log-level', default='warning', help='Logger log level')
    return parser


def main():
    """Main workflow for extracting layers from ARIA products."""
    # Parse command line args
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)
    LOGGER.info('Extract Product Function')

    # Check whether all necessary inputs were specified.
    # some products require a DEM to extract -- if any of those are requested,
    # ensure that a valid DEM is specified
    if args.layers is not None:
        # format list of layers
        layers = list(args.layers.split(','))
        layers = [i.replace(' ', '') for i in layers]
        layers = ['all' if layer.lower() == 'all'
                  else layer for layer in layers]
        layers = ['troposphere*' if layer.startswith('troposphere')
                  else layer for layer in layers]

        # list of layers requiring DEM for extraction
        LAYERS_REQUIRING_DEM = {'all',
                                'bPerpendicular',
                                'bParallel',
                                'incidenceAngle',
                                'lookAngle',
                                'azimuthAngle',
                                'solidEarthTide',
                                'troposphere*'}

        # check that DEM is specified depending on layers requested
        if len(LAYERS_REQUIRING_DEM.intersection(layers)) > 0:
            if args.demfile is None:
                error_msg = 'A valid DEM must be specified when extracting ' \
                            'any of %s' % ', '.join(LAYERS_REQUIRING_DEM)
                LOGGER.error(error_msg)
                raise Exception(error_msg)

    # if user bbox was specified, file(s) not meeting imposed spatial criteria
    # are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped
    # “radarmetadata info” and “data layer keys+paths” dictionaries for each
    # standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if
    # bbox specified)
    standardproduct_info = ARIAtools.product.Product(
        args.imgfile, bbox=args.bbox, projection=args.projection,
        workdir=args.workdir, num_threads=args.num_threads,
        url_version=args.version, nc_version=args.nc_version,
        verbose=args.verbose)

    # Perform initial layer, product, and correction sanity checks
    args.layers, args.tropo_total, \
        model_names = ARIAtools.util.vrt.layerCheck(
            standardproduct_info.products[1], args.layers, args.nc_version,
            args.gacos_products, args.tropo_models, extract_or_ts='extract')

    # pass number of threads for gdal multiprocessing computation
    if args.num_threads.lower() == 'all':
        args.num_threads = 'ALL_CPUS'

    LOGGER.info(
        'Thread count specified for gdal multiprocessing = %s' % (
            args.num_threads))

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection,
    # but user may specify union), and expected shape for DEM.
    LOGGER.info('Extracting and merging product bounding boxes')
    (standardproduct_info.products[0], standardproduct_info.products[1],
     standardproduct_info.bbox_file, prods_TOTbbox,
     prods_TOTbbox_metadatalyr, arrres,
     proj, is_nisar_file) = ARIAtools.extractProduct.merged_productbbox(
        standardproduct_info.products[0], standardproduct_info.products[1],
        os.path.join(args.workdir, 'productBoundingBox'),
        standardproduct_info.bbox_file, args.croptounion,
        num_threads=args.num_threads, minimumOverlap=args.minimumOverlap,
        verbose=args.verbose)

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

        # mask parms
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

    # Extract
    # aria_extract default parms
    export_dict = {
        'full_product_dict': standardproduct_info.products[1],
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'proj': proj,
        'layers': args.layers,
        'is_nisar_file': is_nisar_file,
        'arrres': arrres,
        'rankedResampling': args.rankedResampling,
        'demfile': demfile,
        'demfile_expanded': demfile_expanded,
        'lat': lat,
        'lon': lon,
        'maskfile': maskfilename,
        'outDir': args.workdir,
        'outputFormat': args.outputFormat,
        'verbose': args.verbose,
        'num_threads': args.num_threads,
        'multilooking': args.multilooking,
        'tropo_total': args.tropo_total,
        'model_names': model_names
    }

    # Extract user expected layers
    LOGGER.info('Extracting products')
    arrshape = ARIAtools.extractProduct.export_products(**export_dict)

    # Perform GACOS-based tropospheric corrections (if specified).
    if args.gacos_products:
        LOGGER.info('Applying gacos_correction')
        ARIAtools.extractProduct.gacos_correction(
            standardproduct_info.products, args.gacos_products,
            standardproduct_info.bbox_file, prods_TOTbbox,
            outDir=args.workdir, outputFormat=args.outputFormat,
            verbose=args.verbose, num_threads=args.num_threads)


if __name__ == '__main__':
    main()
