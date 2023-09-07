#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import dask
import logging
import shutil, glob
import os, os.path as op
import time
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
from ARIAtools.ionosphere import export_ionosphere
from ARIAtools.extractProduct import merged_productbbox, prep_dem, prep_metadatalayers, finalize_metadata, handle_epoch_layers
from ARIAtools import progBar
import logging


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


def exportUnwrappedPhase(product_dict,  bbox_file, prods_TOTbbox, arres,
                          work_dir, outputFormat='ISCE', correction_method='cycle2pi',
                          mask_zero_component=False, verbose=True, mask=None, multilook=None, n_jobs=1):

    # verbose printing
    def vprint(x): return print(x) if verbose == True else None

    # initalize multiprocessing
    if n_jobs > 1 or len(product_dict) > 1:
        vprint('Running GUNW unwrappedPhase and connectedComponents in parallel!')
        client = Client(threads_per_worker=1,
                        n_workers=n_jobs,
                        memory_limit='10GB')
        vprint(f'Link: {client.dashboard_link}')
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
        arrres=arres,
        output_unw=None,
        output_conn=None,
        output_format=outputFormat,
        bounds=bounds,
        clip_json=prods_TOTbbox,
        mask_file=mask,
        correction_method=correction_method,
        range_correction=True,
        mask_zero_component=mask_zero_component,
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
    vprint(f'Run {len(jobs)} jobs with {n_jobs} workers.')
    out = dask.compute(*jobs)
    # progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


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
    vprint(f'Run {len(jobs)} jobs with {n_jobs} workers.')
    out = dask.compute(*jobs)
    # progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


def exportImagingGeometry(product_dict, bbox_file, prods_TOTbbox, dem, Latitude, Longitude,
                          work_dir, layer='incidenceAngle', mask=None,
                          outputFormat='ISCE', n_threads=1, n_jobs=1, verbose=True):

    # verbose printing
    def vprint(x): return print(x) if verbose == True else None

    # Check
    available_layers = ['bPerpendicular', 'bParallel',
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

    jobs  = []
    for product, name, outFile in zipped_jobs:
        job_dict['inputfiles'] = product[layer]
        job_dict['outname'] = str(outFile)
        job = dask.delayed(_export_metadata)(**job_dict, dask_key_name=name)
        jobs.append(job)

    # Run export jobs
    vprint(f'Run {len(jobs)} jobs with {n_jobs} workers.')
    out = dask.compute(*jobs)
    # progress(out)  # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


def exportIono(product_dict,  bbox_file, prods_TOTbbox, arres,
                work_dir, outputFormat='ISCE', verbose=True,
                mask=None, multilook=None, n_jobs=1):

    # verbose printing
    vprint = lambda x: print(x) if verbose == True else None

    # initalize multiprocessing
    if n_jobs>1 or len(product_dict)>1:
        vprint('Running GUNW ionosphere in parallel!')
        client = Client(threads_per_worker=1,
                        n_workers=n_jobs,
                        memory_limit='10GB')
        vprint(f'Link: {client.dashboard_link}')
    else:
        client = None

    # Create output directories
    iono_dir = Path(work_dir) / 'ionosphere'
    iono_dir.mkdir(parents=True, exist_ok=True)

    outNames = [ix['pair_name'][0] for ix in product_dict]
    outFileIono = [iono_dir / name for name in outNames]

    export_list = [not(outIono.exists()) for outIono in outFileIono]

    # Get only the pairs we need to run
    product_dict = compress(product_dict, export_list)
    outNames = compress(outNames, export_list)
    outFileIono = compress(outFileIono, export_list)

    zipped_jobs = zip(product_dict, outNames, outFileIono)

    # Get bounds
    bounds = open_shapefile(bbox_file, 0, 0).bounds

    # Run exporting and stitching
    export_dict = dict(
                    input_iono_files = None,
                    arrres = arrres,
                    output_iono = None,
                    output_format = outputFormat,
                    bounds = bounds,
                    clip_json = prods_TOTbbox,
                    mask_file = mask,
                    verbose = verbose,
                    overwrite = True)

    jobs = []
    for product, name, outIono in zipped_jobs:
        export_dict['input_iono_files'] = product['ionosphere']
        export_dict['output_iono'] = outIono

        job = dask.delayed(export_ionosphere)(**export_dict, dask_key_name=name)
        jobs.append(job)

    # Run export jobs
    vprint(f'Run {len(jobs)} jobs with {n_jobs} workers.')
    out = dask.compute(*jobs)
    progress(out) # need to check how to make dask progress bar with dask
    # close dask
    if client:
        client.close()


def exportTropo(product_dict, bbox_file, prods_TOTbbox, dem, Latitude, Longitude, workdir,
                   wmodel='HRRR', layer='troposphereWet',  mask=None,
                   n_threads=1, n_jobs=1, verbose=True, debug=False):

    assert layer in 'troposphereWet troposphereHydrostatic'.split(), 'Wrong layer'
    def vprint(x): return print(x) if verbose == True else None

    # Get bounds of bbox and dem
    bounds = open_shapefile(bbox_file, 0, 0).bounds
    dem_gt = dem.GetGeoTransform()
    dem_bounds = [dem_gt[0],
                dem_gt[3]+(dem_gt[-1]*dem.RasterYSize),
                dem_gt[0]+(dem_gt[1]*dem.RasterXSize),
                dem_gt[3]]

    # initalize multiprocessing
    if debug:
        client = None

    elif n_jobs > 1 or len(product_dict) > 1:
        vprint(f'Running GUNW {layer} in parallel!')
        client = Client(processes=True,
                        threads_per_worker=1,
                        n_workers=n_jobs, memory_limit='20GB')
        vprint(f'Link: {client.dashboard_link}')
    else:
        client = None

    # Check if output files exist
    out_dir = Path(workdir) / wmodel / layer
    outNames = [ix['pair_name'][0] for ix in product_dict]
    outFiles = [out_dir / name for name in outNames]
    export_list = [not (outFile.exists()) for outFile in outFiles]

    # Get only the pairs we need to run
    product_dict = compress(product_dict, export_list)
    outNames     = compress(outNames, export_list)

    job_dict = dict(
        mask=mask,
        dem=dem.GetDescription(),
        dem_bounds=dem_bounds,
        bounds=bounds,
        prods_TOTbbox=prods_TOTbbox,
        verbose=verbose,
        num_threads=n_threads,
        lon=Longitude,
        lat=Latitude,
        layer=layer,
        wmodel=wmodel,
        workdir=workdir
    )

    ref_jobs, sec_jobs  = [], []
    for i, (prod_dict, name) in enumerate(zip(product_dict, outNames)):
        # not all products have tropo layer
        try:
            job_dict['prods'] = prod_dict[f'{layer}_{wmodel.upper()}'] # all tropo ifgs
        except:
            continue

        # this will just export ref/sec and intersect ref with dem
        job_r = dask.delayed(_export_tropo)(**job_dict, dask_key_name=name)
        ref_jobs.append(job_r)

        # this will just intersect the secondary with dem
        job_s = dask.delayed(_export_tropo)(**job_dict, dask_key_name=name)
        sec_jobs.append(job_s)

        if debug:
            break

    # Run export jobs
    vprint(f'Run {len(jobs)} reference jobs with {n_jobs} workers.')
    out = dask.compute(*ref_jobs)

    vprint(f'Catching breath...')
    time.sleep(len(ref_jobs)*2)

    ## then do secondary
    vprint(f'Run {len(jobs)} secondary jobs with {n_jobs} workers.')
    out = dask.compute(*sec_jobs)

    # close dask
    if client:
        client.close()


def move_tropo_layers(path_wd, ifgs, wm='HRRR'):
    """ Move tropo layers from there temporary 'ifg' directory to one AT expects

    path_wd should contain troposphereWet
    ifgs is a list of datepairs (YYYYMMDD_YYYYMMDD)
    """

    ## move the files to one outer
    path_wd = Path(path_wd)
    paths_dry = [path_wd / 'troposphereHydrostatic' / ifg / wm for ifg in ifgs]
    paths_wet = [path_wd / 'troposphereWet' / ifg / wm for ifg in ifgs]

    for paths in [paths_dry, paths_wet]:
        for j, src in enumerate(paths):
            if j == 0:
                dst = src.parents[1]
                shutil.move(src, dst)
            else:
                for f in glob.glob(f'{src}/*'):
                    if op.isdir(f) and op.basename(f) == 'dates':
                        for f1 in glob.glob(f'{f}/*'):
                            # keep the first single date thats unpacked (theoretically always same between ifgs)
                            dst = src.parents[1] / wm / 'dates' / op.basename(f1)
                            if not op.exists(dst):
                                shutil.move(f1, dst)

                    else:
                        shutil.move(f, src.parents[1] / wm )
            shutil.rmtree(src.parents[0]) # clean up temporary ifg dir
    print ('Moved paths to correct directories')
    return


def compute_tropo_total(path_wd, wm='HRRR'):
    """ Add wet and dry for the individual dates """
    import rioxarray as xrr
    i = 0
    for root, dirs, files in os.walk(Path(path_wd) / 'troposphereHydrostatic' / wm / 'dates'):
        for f in files:
            if f.endswith('vrt'):
                continue
            path_hyd = op.join(root, f)
            path_wet = path_hyd.replace('Hydrostatic', 'wet')
            path_tot = path_hyd.replace('Hydrostatic', 'Total')
            os.makedirs(op.dirname(path_tot), exist_ok=True)

            with xrr.open_rasterio(path_hyd) as da_hyd:
                arr_hyd = da_hyd.data
            with xrr.open_rasterio(path_wet) as da_wet:
                arr_wet = da_wet.data

            arr_total = arr_hyd + arr_wet
            da_total  = da_wet.copy()
            da_total.data = arr_total
            da_total.rio.to_raster(path_tot, driver='GTiff')
            gdal.BuildVRT(f'{path_tot}.vrt', path_tot)
            i+=1

    print (f'Wrote troposphereTotal for {i} dates.')


def _export_tropo(prods, layer, wmodel, workdir, bounds, prods_TOTbbox,
                    dem, dem_bounds, lat, lon, num_threads, mask=None, verbose=False):
    ifg         = op.basename(prods[0].split(':')[1]).split('-')[6]
    ref, sec    = ifg.split('_')
    hyd_workdir = op.join(workdir, layer)
    ## double up the ifg so that files dont get reused
    ifg_outname = op.abspath(op.join(hyd_workdir, ifg, ifg))
    ref_outname = op.join(hyd_workdir, ifg, wmodel, 'dates', ref)

    # generate the differential and ref/sec date wet delay for this ifg
    # this overwrites ref and sec; if sec is ref for another, there is a problem
    # not intersected with DEM
    if not op.exists(ref_outname):
        hgt_field, model_name, hyd_outname = \
            prep_metadatalayers(ifg_outname, prods, dem,
                                layer, layer, 'GTiff')
        print ('Finalizing reference:', hyd_outname)

    ## for secondary image
    else:
        hgt_field   = 'NETCDF_DIM_heightsMeta_VALUES'
        hyd_outname = op.join(hyd_workdir, ifg, wmodel, 'dates', sec)
        print ('Finalizing secondary:', hyd_outname)

    # intersect with DEM
    finalize_metadata(hyd_outname, bounds, dem_bounds,
                    prods_TOTbbox, dem, lat, lon, hgt_field,
                    prods, mask, 'GTiff', verbose, num_threads)


def _export_metadata(inputfiles, outname, demfile, aria_dem_dict,
                     bounds, prods_TOTbbox, mask=None, outputFormat='ISCE',
                     n_threads=2, verbose=True):

    layer = inputfiles[0].split('/')[-1]

    # make VRT pointing to metadata layers in standard product
    hgt_field, model_name, ref_outname = prep_metadatalayers(outname, inputfiles,
                                    demfile, layer, [layer])

    if model_name:
        outname = op.join(op.dirname(outname), model_name, op.basename(outname))

    # Interpolate/intersect with DEM before cropping
    finalize_metadata(outname, bounds, aria_dem_dict['dem_bounds'],
                      prods_TOTbbox, demfile, aria_dem_dict['Latitude'],
                      aria_dem_dict['Longitude'], hgt_field, inputfiles,
                      mask, outputFormat, verbose=verbose, num_threads=n_threads)


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

    # NOTE we should store variable above in PICKLE so we can skip
    # preparing inputs every time if not otherwise specify

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
                             mask=inps.mask, n_threads=inps.num_threads,
                             n_jobs=inps.n_jobs)

    # Export Imaging Geometry
    # Note sure do we need lookAngle here
    layers = ['incidenceAngle', 'azimuthAngle', 'lookAngle']
    print('\nExtracting single incidence angle, look angle and azimuth angle '
          'files valid over common interferometric grid')

    # Take a lot of RAM memory per worker, 9GB per scene
    # Dask reports leak - functions need restructuring
    # Maybe something useful is here: https://github.com/dask/distributed/issues/4571

    for layer in layers:
        exportImagingGeometry(product_dict[:1],
                              bbox_file,
                              prods_TOTbbox,
                              demfile, Latitude, Longitude,
                              inps.workdir, layer=layer,
                              mask=inps.mask,
                              n_threads=inps.num_threads, n_jobs=inps.jobs)


    print('\nExtracting perpendicular baseline grids for each '
          'interferogram pair')
    max_jobs = len(product_dict)
    # Hack solution to stop leaking, run dask Client in loop
    # restart cluster/Client after every iteration
    for n in range(0, max_jobs, inps.n_jobs):
        if n + inps.n_jobs > len(product_dict):
            print('Loop:', [n, max_jobs])
            product_subset = product_dict[n:max_jobs]
        else:
            print('Loop:', [n, n + inps.n_jobs])
            product_subset = product_dict[n:n+inps.n_jobs]

        exportImagingGeometry(product_subset,
                              bbox_file,
                              prods_TOTbbox,
                              demfile, Latitude, Longitude,
                              inps.workdir, layer='bPerpendicular',
                              mask=inps.mask,
                              n_threads=inps.num_threads, n_jobs=inps.n_jobs)

    # TODO missing anxiliary products

    # TROPO

    # SET

    #IONO

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
