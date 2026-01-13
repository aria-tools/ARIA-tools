# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import glob
import scipy
import copy
import numpy as np
import logging
import decimal
import osgeo
import warnings

import ARIAtools.constants

LOGGER = logging.getLogger(__name__)

osgeo.gdal.UseExceptions()
# Suppress warnings
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')

# Save file with gdal


def renderVRT(
        fname, data_lyr, geotrans=None, drivername='ENVI',
        gdal_fmt='float32', proj=None, nodata=None, verbose=False):
    """Exports raster and renders corresponding VRT file."""
    GDAL_MAP = {
        'byte': 1, 'int16': 3, 'int32': 5, 'float32': 6, 'float64': 7,
        'cfloat32': 10, 'cfloat64': 11}

    gdalfile = osgeo.gdal.GetDriverByName(drivername).Create(
        fname, data_lyr.shape[1], data_lyr.shape[0], 1, GDAL_MAP[gdal_fmt])

    gdalfile.GetRasterBand(1).WriteArray(data_lyr)

    # If user wishes to update geotrans.
    if geotrans:
        gdalfile.SetGeoTransform(geotrans)

    # If user wishes to update projection.
    if proj:
        gdalfile.SetProjection(proj)

    translate_options_dict = {'format': 'VRT'}
    # If user wishes to set nodata val.
    if nodata is not None:
        gdalfile.GetRasterBand(1).SetNoDataValue(nodata)
        translate_options_dict['noData'] = nodata

    # Finalize VRT
    translate_options = osgeo.gdal.TranslateOptions(**translate_options_dict)
    osgeo.gdal.Translate(fname + '.vrt', gdalfile, options=translate_options)
    return


# Resample raster
def resampleRaster(
        fname, multilooking, bounds, prods_TOTbbox, rankedResampling=False,
        outputFormat='ENVI', num_threads='2'):
    """Resample rasters and update corresponding VRTs."""
    # Get datasource name (inputname)
    if os.path.exists(fname.split('.vrt')[0]):
        inputname = fname
    else:
        fname += '.vrt'
        # Explicitly close
        ds = osgeo.gdal.Open(fname, osgeo.gdal.GA_ReadOnly)
        inputname = ds.GetFileList()[-1]
        ds = None

    # Access original shape
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(
            format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
            multithread=True)
        ds = osgeo.gdal.Warp('', fname, options=warp_options)
        # Get output res
        arrres = [abs(ds.GetGeoTransform()[1]) * multilooking,
                  abs(ds.GetGeoTransform()[-1]) * multilooking]
        ds = None

    # Get geotrans/proj
    with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
        warp_options = osgeo.gdal.WarpOptions(
            format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
            xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
            resampleAlg='near', multithread=True)
        ds = osgeo.gdal.Warp('', fname, options=warp_options)
        geotrans = ds.GetGeoTransform()
        proj = ds.GetProjection()
        ds = None

    # Use pixel function to downsample connected components/unw files
    # based off of frequency of connected components in each window
    if fname.split('/')[-2] == 'connectedComponents' \
            or fname.split('/')[-2] == 'unwrappedPhase':

        # Resample unw phase based off of mode of connected components
        fnameunw = os.path.join(
            '/'.join(fname.split('/')[:-2]), 'unwrappedPhase',
            ''.join(fname.split('/')[-1]).split('.vrt')[0])

        fnameconcomp = os.path.join(
            '/'.join(fname.split('/')[:-2]), 'connectedComponents',
            ''.join(fname.split('/')[-1]).split('.vrt')[0])

        if rankedResampling:
            # open connected components/unw files
            ds_concomp = osgeo.gdal.Open(fnameconcomp)
            ds_concomp_nodata = ds_concomp.GetRasterBand(1).GetNoDataValue()
            ds_concomp = ds_concomp.ReadAsArray()
            ds_concomp = np.ma.masked_where(
                ds_concomp == ds_concomp_nodata, ds_concomp)
            np.ma.set_fill_value(ds_concomp, ds_concomp_nodata)

            ds_unw = osgeo.gdal.Open(fnameunw)
            ds_unw_nodata = ds_unw.GetRasterBand(1).GetNoDataValue()
            ds_unw = ds_unw.ReadAsArray()
            ds_unw = np.ma.masked_where(
                ds_unw == ds_unw_nodata, ds_unw)
            np.ma.set_fill_value(ds_unw, ds_unw_nodata)

            unwmap = []
            for row in range(multilooking, (ds_unw.shape[0]) + multilooking,
                             multilooking):
                unwmap_row = []
                for column in range(multilooking,
                                    (ds_unw.shape[1]) + multilooking,
                                    multilooking):
                    # get subset values
                    subset_concomp = ds_concomp[
                        row - multilooking:row, column - multilooking:column]
                    subset_unw = ds_unw[
                        row - multilooking:row, column - multilooking:column]
                    concomp_mode = scipy.stats.mode(
                        subset_concomp.flatten()).mode[0]

                    # average only phase values coinciding with concomp mode
                    subset_concomp = np.where(
                        subset_concomp != concomp_mode, 0, 1)
                    subset_unw = subset_unw * subset_concomp

                    # assign downsampled pixel values
                    unwmap_row.append(subset_unw.mean())
                unwmap.append(unwmap_row)

            # finalize unw array
            unwmap = np.array(unwmap)

            # finalize unw array shape
            indx0 = int(decimal.Decimal(
                ds_unw.shape[0] / multilooking).quantize(
                    0, decimal.ROUND_HALF_UP))
            indx1 = int(decimal.Decimal(
                ds_unw.shape[1] / multilooking).quantize(
                    0, decimal.ROUND_HALF_UP))
            unwmap = unwmap[0:indx0, 0:indx1]
            unwmap = np.ma.masked_invalid(unwmap)
            np.ma.set_fill_value(unwmap, ds_unw_nodata)

            # Clear variable
            ds_unw = None

            # unwphase
            renderVRT(
                fnameunw, unwmap.filled(), geotrans=geotrans,
                drivername=outputFormat, gdal_fmt='float32', proj=proj,
                nodata=ds_unw_nodata)

            # temp workaround for gdal bug
            try:
                # Assign and close
                ds_check = osgeo.gdal.Open(fnameunw, osgeo.gdal.GA_ReadOnly)
                ds_check = None

            except RuntimeError:
                for f in glob.glob(fnameunw + "*"):
                    os.remove(f)

                unwmap[0, 0] = unwmap[0, 0] - 1e-6
                renderVRT(
                    fnameunw, unwmap.filled(), geotrans=geotrans,
                    drivername=outputFormat, gdal_fmt='float32', proj=proj,
                    nodata=ds_unw_nodata)

            # Resample connected components
            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
                warp_options = osgeo.gdal.WarpOptions(
                    format=outputFormat, cutlineDSName=prods_TOTbbox,
                    outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                    targetAlignedPixels=True, resampleAlg='mode',
                    multithread=True,
                    options=['-overwrite'])
                osgeo.gdal.Warp(
                    fnameconcomp, fnameconcomp, options=warp_options)

                # update VRT
                vrt_options = osgeo.gdal.BuildVRTOptions(
                    options=['-overwrite'])
                osgeo.gdal.BuildVRT(
                    fnameconcomp + '.vrt', fnameconcomp, options=vrt_options)

        # Default: resample unw phase with gdal average algorithm
        else:
            ds = osgeo.gdal.Open(fnameunw, osgeo.gdal.GA_ReadOnly)
            ds_unw_nodata = ds.GetRasterBand(1).GetNoDataValue()
            ds = None

            # Resample unwphase
            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
                warp_options = osgeo.gdal.WarpOptions(
                    format=outputFormat, cutlineDSName=prods_TOTbbox,
                    outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                    targetAlignedPixels=True, resampleAlg='average',
                    multithread=True,
                    options=['-overwrite'])
                osgeo.gdal.Warp(fnameunw, fnameunw, options=warp_options)

            # update VRT
            vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
            osgeo.gdal.BuildVRT(
                fnameunw + '.vrt', fnameunw, options=vrt_options)

            # temp workaround for gdal bug
            try:
                # Assign and close
                ds_check = osgeo.gdal.Open(fnameunw, osgeo.gdal.GA_ReadOnly)
                ds_check = None

            except RuntimeError:
                unwmap = np.fromfile(fnameunw, dtype=np.float32).reshape(
                    ds.GetRasterBand(1).ReadAsArray().shape)

                for f in glob.glob(fnameunw + "*"):
                    os.remove(f)

                unwmap[0, 0] = unwmap[0, 0] - 1e-6
                renderVRT(
                    fnameunw, unwmap, geotrans=geotrans,
                    drivername=outputFormat, gdal_fmt='float32', proj=proj,
                    nodata=ds_unw_nodata)

            # Resample connected components
            with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
                warp_options = osgeo.gdal.WarpOptions(
                    format=outputFormat, cutlineDSName=prods_TOTbbox,
                    outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                    targetAlignedPixels=True, resampleAlg='near',
                    multithread=True, options=['-overwrite'])
                osgeo.gdal.Warp(
                    fnameconcomp, fnameconcomp, options=warp_options)

            # update VRT
            vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
            osgeo.gdal.BuildVRT(
                fnameconcomp + '.vrt', fnameconcomp, options=vrt_options)

    # Resample all other files with lanczos
    else:
        # Resample raster
        with osgeo.gdal.config_options({"GDAL_NUM_THREADS": num_threads}):
            warp_options = osgeo.gdal.WarpOptions(
                format=outputFormat, cutlineDSName=prods_TOTbbox,
                outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                targetAlignedPixels=True, resampleAlg='lanczos',
                multithread=True, options=['-overwrite'])
            osgeo.gdal.Warp(fname, inputname, options=warp_options)

    if outputFormat != 'VRT':
        # update VRT
        vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
        osgeo.gdal.BuildVRT(fname + '.vrt', fname, options=vrt_options)
    return


# Average rasters
def rasterAverage(
        outname, product_dict, bounds, prods_TOTbbox, arrres,
        outputFormat='ENVI', thresh=None):
    """Generate average of rasters."""
    # Make average raster
    # Delete existing average raster file
    for i in glob.glob(outname + '*'):
        os.remove(i)

    # Iterate through all layers
    for i in enumerate(product_dict):
        warp_options = osgeo.gdal.WarpOptions(
            format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
            xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True)
            
        # --- FIX START ---
        # 1. Capture the Warp result
        ds_warp = osgeo.gdal.Warp('', i[1], options=warp_options)
        
        # 2. Read data immediately
        nodata_value = ds_warp.GetRasterBand(1).GetNoDataValue()
        warp_arr = ds_warp.ReadAsArray()
        
        # 3. CRITICAL: Close the Warp dataset
        ds_warp = None 
        # --- FIX END ---

        arr_file_arr = np.ma.masked_where(
            warp_arr == nodata_value, warp_arr)

        # Iteratively update average raster file
        if os.path.exists(outname):
            # Open update file
            ds_update = osgeo.gdal.Open(outname, osgeo.gdal.GA_Update)
            band = ds_update.GetRasterBand(1)
            
            # Read, Add, Write
            current_data = band.ReadAsArray()
            band.WriteArray(arr_file_arr + current_data)
            
            # Close update file
            ds_update = None

        else:
            # If looping through first raster file, nothing to sum so just save
            # Note: We need projection/geotransform. 
            # We can re-open source i[1] briefly or cache it from ds_warp above.
            # Better approach: Cache proj/gt from ds_warp before closing it above.
            
            # (Re-opening source for metadata is safer if ds_warp was MEM)
            ds_src = osgeo.gdal.Open(i[1], osgeo.gdal.GA_ReadOnly)
            renderVRT(
                outname, arr_file_arr, geotrans=ds_src.GetGeoTransform(),
                drivername=outputFormat, gdal_fmt=arr_file_arr.dtype.name,
                proj=ds_src.GetProjection(),
                nodata=nodata_value)
            ds_src = None

    # Take average of raster sum
    ds_avg = osgeo.gdal.Open(outname, osgeo.gdal.GA_Update)
    arr_sum = ds_avg.ReadAsArray()
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        arr_mean = arr_sum / len(product_dict)

    # Mask using specified raster threshold
    if thresh:
        arr_mean = np.where(arr_mean < float(thresh), 0, 1)

    # Save updated array to file
    ds_avg.GetRasterBand(1).WriteArray(arr_mean)
    ds_avg = None  # CLOSE
    arr_mean = None

    # Load raster to pass
    ds_final = osgeo.gdal.Open(outname)
    final_arr = ds_final.ReadAsArray()
    ds_final = None
    
    return final_arr


# Perform initial layer, product, and correction sanity checks
def layerCheck(
        products, layers, nc_version, tropo_models,
        extract_or_ts):
    """Check if any conflicts between netcdf versions and expected layers."""
    # track if product stack is NISAR GUNW or not
    is_nisar_file = False
    track_fileext = products[0]['unwrappedPhase'][0]
    if len(track_fileext.split('.h5')) > 1:
        is_nisar_file = True

    # Ignore productBoundingBoxes & pair-names, they are not raster layers
    IGNORE_LAYERS = [
        'productBoundingBox', 'productBoundingBoxFrames', 'pair_name']

    # TODO Comment on tropo layers???
    RAIDER_TROPO_LAYERS = ['troposphereWet', 'troposphereHydrostatic']

    # Check all available layers in stack
    products = [list(i.keys()) for i in products]
    products = [
        [sub for sub in i if sub not in IGNORE_LAYERS] for i in products]
    all_valid_layers = list(set.union(*map(set, products)))
    all_valid_layers = list(set(all_valid_layers))

    # track tropo model names
    model_names = [i.split('_')[-1] for i in all_valid_layers if '_' in i]
    model_names = list(set(model_names))
    all_valid_layers = [i.split('_')[0] for i in all_valid_layers]
    tropo_total = False

    # If valid argument for tropo models passed, parse to list
    if isinstance(tropo_models, str):
        if tropo_models.lower() == 'all':
            LOGGER.info('All available tropo models are to be extracted')
            tropo_models = copy.deepcopy(
                ARIAtools.constants.ARIA_TROPO_INTERNAL)
        else:
            tropo_models = list(tropo_models.split(','))
            tropo_models = [i.replace(' ', '') for i in tropo_models]
        model_names = list(
            set.intersection(*map(set, [model_names, tropo_models])))
        for i in tropo_models:
            if i not in model_names:
                LOGGER.warning('%s tropo model not found in product', i)
            else:
                LOGGER.info('Generating tropo model %s', i)
    else:
        model_names = []

    # If specified, extract all layers
    if layers:
        if layers.lower() == 'all':
            LOGGER.info('All layers are to be extracted, pass all keys.')
            layers = copy.deepcopy(all_valid_layers)
            if set(RAIDER_TROPO_LAYERS).issubset(all_valid_layers):
                tropo_total = True

        # If valid argument for input layers passed, parse to list
        if isinstance(layers, str):
            layers = list(layers.split(','))
            layers = [i.replace(' ', '') for i in layers]
        if 'troposphereTotal' in layers and \
                set(RAIDER_TROPO_LAYERS).issubset(all_valid_layers) and \
                (model_names != [] or is_nisar_file):
            tropo_total = True

    # differentiate between extract and TS pipeline
    # extract pipeline
    if extract_or_ts == 'extract':
        if not layers:
            LOGGER.info(
                'No layers specified; only creating bounding box shapes')
            return [], [], []

        else:
            layers = [i.replace(' ', '') for i in layers]

    # TS pipeline
    if extract_or_ts == 'tssetup':
        if layers:
            # remove layers already generated in default TS workflow
            layers = [i for i in layers if i not in
                      ARIAtools.constants.ARIA_STANDARD_LAYERS]

        else:
            layers = []

    # pass intersection of valid layers and track invalid requests
    layer_reject = list(
        set.symmetric_difference(*map(set, [all_valid_layers, layers])))

    layer_reject = list(
        set.intersection(*map(set, [layer_reject, RAIDER_TROPO_LAYERS])))

    # only report layers which user requested
    layer_reject = list(
        set.intersection(*map(set, [layer_reject, layers])))
    layers = list(set.intersection(*map(set, [layers, all_valid_layers])))

    if layer_reject != []:
        LOGGER.warning(
            f'User-requested layers {layer_reject} cannot be extracted as '
            'they are not common to all products. Consider fixing input '
            f'"-nc_version {nc_version}" constraint to filter older product '
            'variants')

    # if specified, determine if computation of
    # total tropospheric is possible
    if tropo_total:
        if not set(RAIDER_TROPO_LAYERS).issubset(all_valid_layers):
            LOGGER.warning(
                'User-requested computation of raider-derived total '
                'troposphere "-l troposphereTotal" is not possible as tropo '
                'component layers are not common to all products.')
            tropo_total = False

        if model_names == [] and not is_nisar_file:
            LOGGER.warning(
                'Extraction of raider-derived troposphere layers is not '
                'possible as specified tropo model name(s) '
                f'"-tm {tropo_models}" is not valid.')
            tropo_total = False

    return layers, tropo_total, model_names


def get_basic_attrs(fname):
    """ Access product dimensions and nodata values """
    data_set = osgeo.gdal.Open(fname)
    width = data_set.RasterXSize
    height = data_set.RasterYSize
    geo_trans = data_set.GetGeoTransform()
    proj = data_set.GetProjection()
    no_data = data_set.GetRasterBand(1).GetNoDataValue()
    data_set = None
    return width, height, geo_trans, proj, no_data


def dim_check(ref_arr, prod_arr):
    """Check dimensions between successive products"""
    # Access respective dimensions and geotrans from inputs
    ref_wid = ref_arr[0]
    ref_hgt = ref_arr[1]
    ref_geotrans = ref_arr[2]
    prev_outname = ref_arr[3]
    prod_wid = prod_arr[0]
    prod_hgt = prod_arr[1]
    prod_geotrans = prod_arr[2]
    outname = prod_arr[3]

    if (ref_wid != prod_wid) or (ref_hgt != prod_hgt):
        raise Exception(
            f'Inconsistent product dims between products {prev_outname} and '
            f'{outname}: check respective width ({ref_wid}, {prod_wid}) '
            f'and height ({ref_hgt}, {prod_hgt}) and geotrans '
            f'({ref_geotrans}, {prod_geotrans})')
    return


# Helper to check heights safely
def get_hgt_meta(fname, field):
    ds = osgeo.gdal.Open(fname)
    val = ds.GetMetadataItem(field)
    ds = None # Close immediately

    return val
