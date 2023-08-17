#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import dask
from pathlib import Path
from itertools import compress
from osgeo import gdal

from dask.diagnostics import ProgressBar
from dask.distributed import progress, Client

from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.tsSetup import generate_stack
from ARIAtools.mask_util import prep_mask
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.sequential_stitching import product_stitch_sequential
from ARIAtools.extractProduct import merged_productbbox, prep_dem, prep_metadatalayers, finalize_metadata


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
                        type=str, default='overlap', help='Method applied to '
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

    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap',
                        type=float, default=0.0081, help='Minimum km\u00b2 '
                                                         'area of overlap of scenes wrt specified bounding box.'
                                                         ' Default 0.0081 = 0.0081km\u00b2 = area of single'
                                                         ' pixel at standard 90m resolution')
    parser.add_argument('--n_jobs', dest='n_jobs',
                        type=float, default=1, help='Number of parallel jobs')

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


def export_unwrappedPhase(product_dict,  bbox_file, prods_TOTbbox, arres,
                          work_dir, outputFormat='ISCE', correction_method='cycle2pi',
                          verbose=True, mask=None, multilook=None, n_jobs=1):

    # verbose printing
    def vprint(x): return print(x) if verbose == True else None

    # initalize multiprocessing
    if n_jobs > 1 or len(product_dict) > 1:
        vprint('Running GUNW unwrappedPhase and connectedComponents in parallel!')
        client = Client(threads_per_worker=1,
                        n_workers=n_jobs,
                        memory_limit='10GB')
    else:
        client = None

    # Create output directories
    unw_dir = Path(work_dir) / 'unwrappedPhase'
    unw_dir.mkdir(parents=True, exist_ok=True)

    conncomp_dir = Path(work_dir) / 'connectedComponents'
    conncomp_dir.mkdir(parents=True, exist_ok=True)

    outNames = [ix['pair_name'][0] for ix in product_dict]
    outFilePhs = [unw_dir / name for name in outNames]
    outFileConnComp = [conncomp_dir / name for name in outNames]

    zipped_out = zip(outFilePhs, outFileConnComp)
    export_list = [not (outUnw.exists() and outConnComp.exists())
                   for outUnw, outConnComp in zipped_out]

    # Get only the pairs we need to run
    product_dict = compress(product_dict, export_list)
    outNames = compress(outNames, export_list)
    outFilePhs = compress(outFilePhs, export_list)
    outFileConnComp = compress(outFileConnComp, export_list)

    zipped_jobs = zip(product_dict, outNames,
                      outFilePhs, outFileConnComp)

    # Get bounds
    bounds = open_shapefile(bbox_file, 0, 0).bounds

    # Run exporting and stitching
    export_dict = dict(
        input_unw_files=None,
        input_conncomp_files=None,
        arrres=arrres,
        output_unw=None,
        output_conn=None,
        output_format=outputFormat,
        bounds=bounds,
        clip_json=prods_TOTbbox,
        mask_file=mask,
        correction_method=correction_method,
        range_correction=True,
        verbose=False,
        save_fig=False,
        overwrite=True,
    )

    jobs = []
    for product, name, outUnw, outConnComp in zipped_jobs:
        export_dict['input_unw_files'] = product['unwrappedPhase']
        export_dict['input_conncomp_files'] = product['connectedComponents']
        export_dict['output_unw'] = outUnw
        export_dict['output_conn'] = outConnComp

        job = dask.delayed(product_stitch_sequential)(
            **export_dict, dask_key_name=name)
        jobs.append(job)
    # Run export jobs
    vprint(f'Run number of jobs: {len(jobs)}')
    out = dask.compute(*jobs)
    progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


def _gdal_export(inputfiles, outname, gdal_warp_kwargs, mask=None, outputFormat='ISCE'):
    if outputFormat == 'VRT':
        # building the virtual vrt
        gdal.BuildVRT(outname + "_uncropped" + '.vrt', inputfiles)
        # building the cropped vrt
        gdal.Warp(outname+'.vrt',
                  outname+'_uncropped.vrt',
                  options=gdal.WarpOptions(**gdal_warp_kwargs))
    else:
        # building the VRT
        gdal.BuildVRT(outname + '.vrt', inputfiles)
        gdal.Warp(outname,
                  outname+'.vrt',
                  options=gdal.WarpOptions(**gdal_warp_kwargs))

        # Update VRT
        gdal.Translate(outname+'.vrt', outname,
                       options=gdal.TranslateOptions(format="VRT"))

        # Apply mask (if specified).
        if mask is not None:
            if isinstance(mask, str):
                mask_file = gdal.Open(mask)
            else:
                # for gdal instance, from prep_mask
                mask_file = mask
            update_file = gdal.Open(outname, gdal.GA_Update)
            mask_arr = mask_file.ReadAsArray() * \
                gdal.Open(outname + '.vrt').ReadAsArray()
            update_file.GetRasterBand(1).WriteArray(mask_arr)
            del update_file, mask_arr


def exportCoherenceAmplitude(product_dict, bbox_file, prods_TOTbbox, arrres,
                             work_dir, layer='coherence', mask=None, outputFormat='ISCE',
                             n_threads=1, n_jobs=1, verbose=True):
    # Layer: coherence or amplitude
    # verbose printing
    def vprint(x): return print(x) if verbose == True else None

    # Check
    if layer not in ['coherence', 'amplitude']:
        raise ValueError(f'Selected layer: {layer} is wrong'
                         ' ,available choices: coherence, amplitude')

    # initalize multiprocessing
    if n_jobs > 1 or len(product_dict) > 1:
        vprint(f'Running GUNW {layer} in parallel!')
        client = Client(threads_per_worker=1,
                        n_workers=n_jobs,
                        memory_limit='10GB')
    else:
        client = None

    # Create output directories
    out_dir = Path(work_dir) / layer
    out_dir.mkdir(parents=True, exist_ok=True)

    outNames = [ix['pair_name'][0] for ix in product_dict]
    outFiles = [out_dir / name for name in outNames]
    export_list = [not (outFile.exists()) for outFile in outFiles]

    # Get only the pairs we need to run
    product_dict = compress(product_dict, export_list)
    outNames = compress(outNames, export_list)
    outFiles = compress(outFiles, export_list)

    zipped_jobs = zip(product_dict, outNames, outFiles)

    # Get bounds
    bounds = open_shapefile(bbox_file, 0, 0).bounds

    gdal_warp_kwargs = {'format': outputFormat,
                        'cutlineDSName': prods_TOTbbox,
                        'outputBounds': bounds,
                        'xRes': arrres[0],
                        'yRes': arrres[1],
                        'targetAlignedPixels': True,
                        'multithread': True,
                        'options': [f'NUM_THREADS={n_threads}']}

    job_dict = dict(
        inputfiles=None,
        outname=None,
        mask=mask,
        outputFormat=outputFormat,
        gdal_warp_kwargs=gdal_warp_kwargs
    )

    jobs = []
    for product, name, outFile in zipped_jobs:
        job_dict['inputfiles'] = product[layer]
        job_dict['outname'] = str(outFile)
        job = dask.delayed(_gdal_export)(**job_dict, dask_key_name=name)
        jobs.append(job)
    # Run export jobs
    vprint(f'Run number of jobs: {len(jobs)}')
    out = dask.compute(*jobs)
    progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


def _export_metadata(inputfiles, outname, demfile, aria_dem_dict,
                     bounds, prods_TOTbbox, mask=None, outputFormat='ISCE',
                     n_threads=2, verbose=True):

    layer = inputfiles[0].split('/')[-1]

    # make VRT pointing to metadata layers in standard product
    hgt_field = prep_metadatalayers(outname, inputfiles,
                                    demfile, layer, layer)[0]

    # Interpolate/intersect with DEM before cropping
    finalize_metadata(outname, bounds, aria_dem_dict['dem_bounds'],
                      prods_TOTbbox, demfile, aria_dem_dict['Latitude'],
                      aria_dem_dict['Longitude'], hgt_field, inputfiles,
                      mask, outputFormat, verbose=verbose, num_threads=n_threads)


def exportImagingGeometry(product_dict, bbox_file, prods_TOTbbox, dem, Latitude, Longitude,
                          work_dir, layer='incidenceAngle', mask=None, outputFormat='ISCE',
                          n_threads=1, n_jobs=1, verbose=True):
    # verbose printing
    def vprint(x): return print(x) if verbose == True else None

    # Check
    available_layers = ['perpendicularBaseline', 'parallelBaseline',
                        'incidenceAngle', 'lookAngle', 'azimuthAngle']
    if layer not in available_layers:
        raise ValueError(f'Selected layer: {layer} is wrong'
                         f' ,available choices: {available_layers}')

    # initalize multiprocessing
    if n_jobs > 1 or len(product_dict) > 1:
        vprint(f'Running GUNW {layer} in parallel!')
        client = Client(processes=True,
                        threads_per_worker=1,
                        n_workers=n_jobs, memory_limit='20GB')
        vprint(f'Link: {client.dashboard_link}')
    else:
        client = None

    # Create output directories
    out_dir = Path(work_dir) / layer
    out_dir.mkdir(parents=True, exist_ok=True)

    outNames = [ix['pair_name'][0] for ix in product_dict]
    outFiles = [out_dir / name for name in outNames]
    export_list = [not (outFile.exists()) for outFile in outFiles]

    # Get only the pairs we need to run
    product_dict = compress(product_dict, export_list)
    outNames = compress(outNames, export_list)
    outFiles = compress(outFiles, export_list)

    zipped_jobs = zip(product_dict, outNames, outFiles)

    # Get bounds
    bounds = open_shapefile(bbox_file, 0, 0).bounds

    # original DEM shape
    aria_dem = {}
    aria_dem['dem_bounds'] = dem.GetGeoTransform()
    aria_dem['Longitude'] = Longitude
    aria_dem['Latitude'] = Latitude

    job_dict = dict(
        inputfiles=None,
        outname=None,
        mask=mask,
        demfile=dem.GetDescription(),
        aria_dem_dict=aria_dem,
        bounds=bounds,
        prods_TOTbbox=prods_TOTbbox,
        outputFormat=outputFormat,
        verbose=verbose,
        n_threads=n_threads
    )

    jobs = []
    for product, name, outFile in zipped_jobs:
        job_dict['inputfiles'] = product[layer]
        job_dict['outname'] = str(outFile)
        job = dask.delayed(_export_metadata)(**job_dict, dask_key_name=name)
        jobs.append(job)
    # Run export jobs
    vprint(f'Run number of jobs: {len(jobs)} with {n_jobs} workers')
    out = dask.compute(*jobs)
    progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


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
    (metadata_dict, product_dict, bbox_file, prods_TOTbbox,
     prods_TOTbbox_metadatalyr, arrres, proj) = merged_productbbox(standardproduct_info.products[0],
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
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(
            inps.demfile, standardproduct_info.bbox_file,
            prods_TOTbbox, prods_TOTbbox_metadatalyr, proj,
            arrres=arrres, workdir=inps.workdir,
            outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Load or download mask (if specified).
    if inps.mask is not None:
        # Extract amplitude layers
        amplitude_products = []
        for d in standardproduct_info.products[1]:
            if 'amplitude' in d:
                for item in list(set(d['amplitude'])):
                    amplitude_products.append(item)
        inps.mask = prep_mask(amplitude_products, inps.mask,
                              standardproduct_info.bbox_file,
                              prods_TOTbbox, proj, amp_thresh=inps.amp_thresh,
                              arrres=arrres,
                              workdir=inps.workdir,
                              outputFormat=inps.outputFormat,
                              num_threads=inps.num_threads)

    # export unwrappedPhase
    layers = ['unwrappedPhase', 'coherence']
    print('\nExtracting unwrapped phase, coherence, '
          'and connected components for each interferogram pair')
    export_unwrappedPhase(product_dict,
                          bbox_file,
                          prods_TOTbbox, arrres, inps.workdir,
                          mask=inps.mask, n_jobs=inps.n_jobs)

    exportCoherenceAmplitude(product_dict,
                             bbox_file,
                             prods_TOTbbox, arrres, inps.workdir,
                             mask=inps.mask, n_threads=inps.num_threads, n_jobs=inps.n_jobs)

    # Export Imaging Geometry
    # Note sure do we need lookAngle here
    layers = ['incidenceAngle', 'azimuthAngle', 'lookAngle']
    print('\nExtracting single incidence angle, look angle and azimuth angle '
          'files valid over common interferometric grid')

    # Take a lot of RAM memory per worker, 9GB per scene
    # Dask reports leak - functions need restructuring
    # This would be around solution, not perfect but ..

    for layer in layers[:-1]:
        max_jobs = len(product_dict)
        # Hack solution to stop leaking, run dask Client in loop
        # restart cluster/Client after every iteration
        for n in range(0, max_jobs, inps.n_jobs):
            if n + n_jobs > len(product_dict):
                print('Loop:', [n, max_jobs])
                product_subset = product_dict[n:max_jobs]
            else:
                print('Loop:', [n, n + n_jobs])
                product_subset = product_dict[n:n+inps.n_jobs]

        exportImagingGeometry(product_subset,
                              bbox_file,
                              prods_TOTbbox,
                              demfile, Latitude, Longitude,
                              inps.workdir, layer=layer,
                              mask=inps.mask,
                              n_threads=inps.num_threads, n_jobs=inps.jobs)

    # MG did not test how it works on bPerp
    print('\nExtracting perpendicular baseline grids for each '
          'interferogram pair')
    exportImagingGeometry(product_dict,
                          bbox_file,
                          prods_TOTbbox,
                          demfile, Latitude, Longitude,
                          inps.workdir, layer='bPerpendicular',
                          mask=inps.mask,
                          n_threads=inps.num_threads, n_jobs=1)

    # TODO missing anxiliary products

    # Generate UNW stack
    ref_dlist = generate_stack(standardproduct_info, 'unwrappedPhase',
                               'unwrapStack', workdir=inps.workdir)
    stack_dict = {
        'workdir': inps.workdir,
        'ref_dlist': ref_dlist
    }
    generate_stack(standardproduct_info,
                   'incidenceAngle',
                   ARIA_STACK_OUTFILES['incidenceAngle'],
                   **stack_dict)

    generate_stack(standardproduct_info,
                   'azimuthAngle',
                   ARIA_STACK_OUTFILES['azimuthAngle'],
                   **stack_dict)

    generate_stack(standardproduct_info,
                   'bPerpendicular',
                   ARIA_STACK_OUTFILES['bPerpendicular'],
                   **stack_dict)
