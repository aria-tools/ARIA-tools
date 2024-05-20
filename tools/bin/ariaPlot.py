#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import logging
import argparse

import osgeo
import numpy as np

import ARIAtools.product
import ARIAtools.util.plot
import ARIAtools.extractProduct
import ARIAtools.util.log
import ARIAtools.util.mask

osgeo.gdal.UseExceptions()

# Suppress warnings
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')
LOGGER = logging.getLogger('ariaPlot.py')


def createParser():
    '''
    Make any of the following specified plot(s): ⊥ baseline + histogram,
    coherence + histogram + average coherence raster, ⊥ baseline &
    coherence combo, and track extents. The default is to generate all of
    these.
    '''
    parser = argparse.ArgumentParser(
        description='Function to generate various quality control and '
                    'baseline figures of the spatial-temporal network of '
                    'products.')
    parser.add_argument(
        '-f', '--file', dest='imgfile', type=str, required=True,
        help='ARIA file')
    parser.add_argument(
        '-w', '--workdir', dest='workdir', default='./',
        help='Specify directory to deposit all outputs. Default is local '
             'directory where script is launched.')
    parser.add_argument(
        '-b', '--bbox', dest='bbox', type=str, default=None,
        help='Provide either valid shapefile or Lat/Lon Bounding SNWE. -- '
             'Example : "19 20 -99.5 -98.5"')
    parser.add_argument(
        '-m', '--mask', dest='mask', type=str, default=None,
        help='Path to mask file or "Download". File needs to be GDAL '
             'compatabile, contain spatial reference information, and '
             'have invalid/valid data represented by 0/1, respectively. '
             'If "Download", will use GSHHS water mask. If "NLCD", will mask '
             'classes 11, 12, 90, 95; see: https://www.mrlc.gov/national'
             '-land-cover-database-nlcd-201://www.mrlc.gov/national-land-cover'
             '-database-nlcd-2016')
    parser.add_argument(
        '-at', '--amp_thresh', dest='amp_thresh', default=None, type=str,
        help='Amplitude threshold below which to mask. Specify "None" to '
             'not use amplitude mask. By default "None".')
    parser.add_argument(
        '-nt', '--num_threads', dest='num_threads', default='2', type=str,
        help='Specify number of threads for multiprocessing operation in '
             'gdal. By default "2". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument(
        '-of', '--outputFormat', dest='outputFormat', type=str, default='ENVI',
        help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By '
        'default files are generated with ENVI format.')
    parser.add_argument(
        '-croptounion', '--croptounion', action='store_true',
        dest='croptounion',
        help='If turned on, IFGs cropped to bounds based off of union and '
             'bbox (if specified). Program defaults to crop all IFGs to '
             'bounds based off of common intersection and bbox (if '
             'specified).')
    parser.add_argument(
        '-plottracks', '--plottracks', action='store_true', dest='plottracks',
        help='Make plot of track latitude extents vs bounding bbox/common '
             'track extent.')
    parser.add_argument(
        '-plotbperp', '--plotbperp', action='store_true', dest='plotbperp',
        help="Make a baseline plot, and a histogram of perpendicular "
             "baseline.")
    parser.add_argument(
        '-plotbperpcoh', '--plotbperpcoh', action='store_true',
        dest='plotbperpcoh',
        help='Make a baseline plot that is color-coded based on average IFG '
             'coherence.')
    parser.add_argument(
        '-plotcoh', '--plotcoh', action='store_true', dest='plotcoh',
        help='Make an average IFG coherence plot in time, and histogram of '
             'IFG average coherence.')
    parser.add_argument(
        '-makeavgoh', '--makeavgoh', action='store_true', dest='makeavgoh',
        help="Generate a 2D raster of average IFG coherence.")
    parser.add_argument(
        '-plotall', '--plotall', action='store_true', dest='plotall',
        help="Generate all above plots.")
    parser.add_argument(
        '-mo', '--minimumOverlap', dest='minimumOverlap', type=float,
        default=0.0081,
        help='Minimum km\u00b2 area of overlap of scenes wrt specified '
             'bounding box. Default 0.0081 = 0.0081km\u00b2=area of single '
             'pixel at standard 90m resolution')
    parser.add_argument(
        '--figwidth', dest='figwidth', type=str, default='standard',
        help='Width of lat extents figure in inches. Default is "standard", '
             'i.e., the 6.4-inch-wide standard figure size. Optionally, the '
             'user may define the width manually, e.g,. 8 [inches] or set '
             'the parameter to "wide" format, i.e., the width of the figure '
             'automatically scales with the number of interferograms. Other '
             'options include')
    parser.add_argument(
        '--version', dest='version', default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods; default: all')
    parser.add_argument(
        '--nc_version', dest='nc_version', default='1b',
        help='Specify netcdf version as str, e.g. 1c or all prods; default: '
             '1b')
    parser.add_argument(
        '-v', '--verbose', action='store_true', dest='verbose',
        help="Toggle verbose mode on.")
    parser.add_argument(
        '--log-level', default='info', help='Logger log level')
    return parser


def main(inps=None):
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)
    print('*****************************************************************')
    print('*** Plotting Function ***')
    print('*****************************************************************')
    # if user bbox was specified, file(s) not meeting imposed spatial criteria
    # are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped
    # “radarmetadata info” and “data layer keys+paths” dictionaries for each
    # standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if
    # bbox specified)
    standardproduct_info = ARIAtools.product.Product(
        args.imgfile, bbox=args.bbox, workdir=args.workdir,
        num_threads=args.num_threads, url_version=args.version,
        nc_version=args.nc_version, verbose=args.verbose)

    # If user requests to generate all plots.
    if args.plotall:
        LOGGER.info('"-plotall"==True. All plots will be made.')
        args.plottracks = True
        args.plotbperp = True
        args.plotcoh = True
        args.plotbperpcoh = True
        args.makeavgoh = True

    # pass number of threads for gdal multiprocessing computation
    if args.num_threads.lower() == 'all':
        import multiprocessing
        LOGGER.info(
            'User specified use of all %s threads for gdal multiprocessing',
            multiprocessing.cpu_count())
        args.num_threads = 'ALL_CPUS'
    LOGGER.info(
        'Thread count specified for gdal multiprocessing = %s',
        args.num_threads)

    if args.plottracks or args.plotcoh or args.makeavgoh or args.plotbperpcoh:
        # extract/merge productBoundingBox layers for each pair and update
        # dict, report common track bbox (default is to take common
        # intersection, but user may specify union), and expected shape for
        # DEM.

        # TODO make LHS a tuple
        standardproduct_info.products[0], standardproduct_info.products[1], \
            standardproduct_info.bbox_file, prods_TOTbbox, \
            prods_TOTbbox_metadatalyr, arrres, proj, is_nisar_file = \
                ARIAtools.extractProduct.merged_productbbox(
                    standardproduct_info.products[0],
                    standardproduct_info.products[1],
                    os.path.join(args.workdir, 'productBoundingBox'),
                    standardproduct_info.bbox_file, args.croptounion,
                    num_threads=args.num_threads,
                    minimumOverlap=args.minimumOverlap, verbose=args.verbose)

        # Load or download mask (if specified).
        if args.mask is not None:
            # TODO refactor these list comps, this is not understandable
            args.mask = ARIAtools.util.mask.prep_mask(
                [[item for sublist in [list(set(d['amplitude'])) for d in standardproduct_info.products[1] if 'amplitude' in d] for item in sublist],
                                   [item for sublist in [list(set(d['pair_name'])) for d in standardproduct_info.products[1] if 'pair_name' in d] for item in sublist]],
                                  args.mask,
                                  standardproduct_info.bbox_file,
                                  prods_TOTbbox,
                                  proj,
                                  amp_thresh=args.amp_thresh,
                                  arrres=arrres,
                                  workdir=args.workdir,
                                  outputFormat=args.outputFormat,
                                  num_threads=args.num_threads)

    # Make spatial extent plot
    if args.plottracks:
        LOGGER.info(
            "- Make plot of track latitude extents vs bounding bbox/common track extent.")
        make_plot = ARIAtools.util.plot.PlotClass([[j['productBoundingBox'] for j in standardproduct_info.products[1]],
                               [j["pair_name"] for j in standardproduct_info.products[1]]],
                              workdir=args.workdir,
                              bbox_file=standardproduct_info.bbox_file,
                              prods_TOTbbox=prods_TOTbbox,
                              arrres=arrres,
                              croptounion=args.croptounion)
        make_plot.plot_extents(figwidth=args.figwidth)

    # Make pbaseline plot
    if args.plotbperp:
        LOGGER.info("- Make baseline plot and histogram.")
        make_plot = ARIAtools.util.plot.PlotClass([[j['bPerpendicular'] for j in standardproduct_info.products[1]], [
                              j["pair_name"] for j in standardproduct_info.products[1]]], arrres=arrres, workdir=args.workdir)
        make_plot.plot_pbaselines()

    # Make average land coherence plot
    if args.plotcoh:
        LOGGER.info(
            "- Make average IFG coherence plot in time, and histogram of average IFG coherence.")
        make_plot = ARIAtools.util.plot.PlotClass([[j['coherence'] for j in standardproduct_info.products[1]],
                               [j["pair_name"] for j in standardproduct_info.products[1]]],
                              workdir=args.workdir,
                              bbox_file=standardproduct_info.bbox_file,
                              prods_TOTbbox=prods_TOTbbox,
                              arrres=arrres,
                              mask=args.mask,
                              num_threads=args.num_threads)
        make_plot.plot_coherence()

    # Generate average land coherence raster
    if args.makeavgoh:
        LOGGER.info("- Generate 2D raster of average coherence.")
        make_plot = ARIAtools.util.plot.PlotClass([[j['coherence'] for j in standardproduct_info.products[1]],
                               [j["pair_name"] for j in standardproduct_info.products[1]]],
                              workdir=args.workdir,
                              bbox_file=standardproduct_info.bbox_file,
                              prods_TOTbbox=prods_TOTbbox,
                              arrres=arrres,
                              mask=args.mask,
                              outputFormat=args.outputFormat,
                              num_threads=args.num_threads)
        make_plot.plot_avgcoherence()

    # Make pbaseline/coherence combo plot
    if args.plotbperpcoh:
        LOGGER.info(
            "- Make baseline plot that is color-coded with respect to mean IFG coherence.")
        make_plot = ARIAtools.util.plot.PlotClass([[j['bPerpendicular'] for j in standardproduct_info.products[1]],
                               [j["pair_name"]
                                   for j in standardproduct_info.products[1]],
                               [j['coherence'] for j in standardproduct_info.products[1]]],
                              workdir=args.workdir,
                              bbox_file=standardproduct_info.bbox_file,
                              prods_TOTbbox=prods_TOTbbox,
                              arrres=arrres,
                              mask=args.mask)
        make_plot.plotbperpcoh()

if __name__ == '__main__':
    main()
