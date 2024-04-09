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


# Make OGR VRT file
def renderOGRVRT(vrt_filename, src_datasets):
    """Generate VRT of shapefile unions."""
    VRT_HEAD = ('<OGRVRTDataSource>\n  <OGRVRTUnionLayer name="merged">\n    '
                '<FieldStrategy>Union</FieldStrategy>\n')
    VRT_TAIL = '  </OGRVRTUnionLayer>\n</OGRVRTDataSource>\n'

    with open(vrt_filename, 'w') as ofp:
        ofp.write(VRT_HEAD)
        for i in enumerate(src_datasets):
            ofp.write(
                '    <OGRVRTLayer name="Dataset%i_%s">\n' %(
                i[0], os.path.basename(i[1]).split('.shp')[0]))
            ofp.write(
                '      <SrcDataSource shared="1">%s</SrcDataSource>\n' %(i[1]))
            ofp.write(
                '      <SrcLayer>%s</SrcLayer>\n' %(
                os.path.basename(i[1]).split('.shp')[0]))
            ofp.write('    </OGRVRTLayer>\n')
        ofp.write(VRT_TAIL)
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
        inputname = osgeo.gdal.Open(fname).GetFileList()[-1]

    # Access original shape
    warp_options = osgeo.gdal.WarpOptions(
        format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
        multithread=True, options=['NUM_THREADS=%s' % (num_threads)])
    ds = osgeo.gdal.Warp('', fname, options=warp_options)

    # Get output res
    arrres = [abs(ds.GetGeoTransform()[1]) * multilooking,
              abs(ds.GetGeoTransform()[-1]) * multilooking]

    # Get geotrans/proj
    warp_options = osgeo.gdal.WarpOptions(
        format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
        xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
        resampleAlg='near', multithread=True,
        options=['NUM_THREADS=%s' % (num_threads)])
    ds = osgeo.gdal.Warp('', fname, options=warp_options)
    geotrans = ds.GetGeoTransform()
    proj = ds.GetProjection()

    # Must resample mask with nearest-neighbor
    if fname.split('/')[-2] == 'mask':

        # Resample raster
        warp_options = osgeo.gdal.WarpOptions(
            format=outputFormat, cutlineDSName=prods_TOTbbox,
            outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
            targetAlignedPixels=True, resampleAlg='near', multithread=True,
            options=['NUM_THREADS=%s' % (num_threads) + ' -overwrite'])
        osgeo.gdal.Warp(fname, inputname, options=warp_options)

    # Use pixel function to downsample connected components/unw files
    # based off of frequency of connected components in each window
    elif fname.split('/')[-2] == 'connectedComponents' \
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
            indx0 = int(decimal.Decimal(ds_unw.shape[0] / multilooking).quantize(
                0, decimal.ROUND_HALF_UP))
            indx1 = int(decimal.Decimal(ds_unw.shape[1] / multilooking).quantize(
                0, decimal.ROUND_HALF_UP))
            unwmap = unwmap[0:indx0, 0:indx1]
            unwmap = np.ma.masked_invalid(unwmap)
            np.ma.set_fill_value(unwmap, ds_unw_nodata)

            # unwphase
            renderVRT(
                fnameunw, unwmap.filled(), geotrans=geotrans,
                drivername=outputFormat, gdal_fmt='float32', proj=proj,
                nodata=ds_unw_nodata)

            # temp workaround for gdal bug
            try:
                osgeo.gdal.Open(fnameunw)

            except RuntimeError:
                for f in glob.glob(fnameunw + "*"):
                    os.remove(f)

                unwmap[0, 0] = unwmap[0, 0] - 1e-6
                renderVRT(
                    fnameunw, unwmap.filled(), geotrans=geotrans,
                    drivername=outputFormat, gdal_fmt='float32', proj=proj,
                    nodata=ds_unw_nodata)

            # Resample connected components
            warp_options = osgeo.gdal.WarpOptions(
                format=outputFormat, cutlineDSName=prods_TOTbbox,
                outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                targetAlignedPixels=True, resampleAlg='mode', multithread=True,
                options=['NUM_THREADS=%s' % (num_threads) + ' -overwrite'])
            osgeo.gdal.Warp(fnameconcomp, fnameconcomp,options=warp_options)

            # update VRT
            vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
            osgeo.gdal.BuildVRT(
                fnameconcomp + '.vrt', fnameconcomp, options=vrt_options)

        # Default: resample unw phase with gdal average algorithm
        else:
            # Resample unwphase
            ds_unw_nodata = osgeo.gdal.Open(fnameunw)
            ds_unw_nodata = ds_unw_nodata.GetRasterBand(1).GetNoDataValue()

            warp_options = osgeo.gdal.WarpOptions(
                format=outputFormat, cutlineDSName=prods_TOTbbox,
                outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                targetAlignedPixels=True, resampleAlg='average',
                multithread=True,
                options=['NUM_THREADS=%s' % (num_threads) + ' -overwrite'])
            osgeo.gdal.Warp(fnameunw, fnameunw, options=warp_options)

            # update VRT
            vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
            osgeo.gdal.BuildVRT(
                fnameunw + '.vrt', fnameunw, options=vrt_options)

            # temp workaround for gdal bug
            try:
                osgeo.gdal.Open(fnameunw)

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
            warp_options = osgeo.gdal.WarpOptions(
                format=outputFormat, cutlineDSName=prods_TOTbbox,
                outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
                targetAlignedPixels=True, resampleAlg='near', multithread=True,
                options=['NUM_THREADS=%s' % (num_threads) + ' -overwrite'])
            osgeo.gdal.Warp(fnameconcomp, fnameconcomp, options=warp_options)

            # update VRT
            vrt_options = osgeo.gdal.BuildVRTOptions(options=['-overwrite'])
            osgeo.gdal.BuildVRT(
                fnameconcomp + '.vrt', fnameconcomp, options=vrt_options)

    # Resample all other files with lanczos
    else:
        # Resample raster
        warp_options = osgeo.gdal.WarpOptions(
            format=outputFormat, cutlineDSName=prods_TOTbbox,
            outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
            targetAlignedPixels=True, resampleAlg='lanczos', multithread=True,
            options=['NUM_THREADS=%s' % (num_threads) + ' -overwrite'])
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
    """Generate average of rasters.

    Currently implemented for:
    1. amplitude under 'mask_util.py'
    2. coherence under 'plot_avgcoherence' function of productPlot"""
    # Make average raster
    # Delete existing average raster file
    for i in glob.glob(outname + '*'):
        os.remove(i)

    # Iterate through all layers
    for i in enumerate(product_dict):
        warp_options = osgeo.gdal.WarpOptions(
            format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds,
            xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True)
        arr_file = osgeo.gdal.Warp('', i[1], options=warp_options)
        nodata_value = arr_file.GetRasterBand(1).GetNoDataValue()
        arr_file_arr = np.ma.masked_where(
            arr_file.ReadAsArray() == nodata_value, arr_file.ReadAsArray())

        # Iteratively update average raster file
        if os.path.exists(outname):
            arr_file = osgeo.gdal.Open(outname, osgeo.gdal.GA_Update)
            arr_file = arr_file.GetRasterBand(1).WriteArray(
                arr_file_arr + arr_file.ReadAsArray())

        else:
            # If looping through first raster file, nothing to sum so just save
            # to file
            renderVRT(
                outname, arr_file_arr, geotrans=arr_file.GetGeoTransform(),
                drivername=outputFormat, gdal_fmt=arr_file_arr.dtype.name,
                proj=arr_file.GetProjection(),
                nodata=arr_file.GetRasterBand(1).GetNoDataValue())

    # Take average of raster sum
    arr_file = osgeo.gdal.Open(outname, osgeo.gdal.GA_Update)
    arr_mean = arr_file.ReadAsArray() / len(product_dict)

    # Mask using specified raster threshold
    if thresh:
        arr_mean = np.where(arr_mean < float(thresh), 0, 1)

    # Save updated array to file
    arr_file = arr_file.GetRasterBand(1).WriteArray(arr_mean)
    arr_file = None
    arr_mean = None

    # Load raster to pass
    arr_file = osgeo.gdal.Open(outname).ReadAsArray()
    return arr_file


# Generate GACOS rsc
def rscGacos(outvrt, merged_rsc, tropo_date_dict):
    """Generate GACOS product .rsc files"""
    in_file = osgeo.gdal.Open(outvrt)
    date = os.path.basename(outvrt[:-4]).split('.tif')[0].split('.ztd')[0]

    # Create merge rsc file
    with open(outvrt[:-4] + '.rsc', 'w') as ofp:
        ofp.write('WIDTH %s\n' % (in_file.RasterXSize))
        ofp.write('FILE_LENGTH %s\n' % (in_file.RasterYSize))
        ofp.write('XMIN %s\n' % (0))
        ofp.write('XMAX %s\n' % (in_file.RasterXSize))
        ofp.write('YMIN %s\n' % (0))
        ofp.write('YMAX %s\n' % (in_file.RasterYSize))
        ofp.write('X_FIRST %f\n' % (in_file.GetGeoTransform()[0]))
        ofp.write('Y_FIRST %f\n' % (in_file.GetGeoTransform()[3]))
        ofp.write('X_STEP %f\n' % (in_file.GetGeoTransform()[1]))
        ofp.write('Y_STEP %f\n' % (in_file.GetGeoTransform()[-1]))
        ofp.write('X_UNIT %s\n' % ('degres'))
        ofp.write('Y_UNIT %s\n' % ('degres'))
        ofp.write('Z_OFFSET %s\n' % (0))
        ofp.write('Z_SCALE %s\n' % (1))
        ofp.write('PROJECTION %s\n' % ('LATLON'))
        ofp.write('DATUM %s\n' % ('WGS84'))
        ofp.write('TIME_OF_DAY %s\n' % (''.join(
            tropo_date_dict[date + "_UTC"])))

    return


# Parse GACOS metadata from TIF
def tifGacos(intif):
    """Pass GACOS metadata as dict from TIF"""
    tropo_rsc_dict = {}
    in_file = osgeo.gdal.Open(intif)

    # Populate dict
    tropo_rsc_dict['WIDTH'] = in_file.RasterXSize
    tropo_rsc_dict['FILE_LENGTH'] = in_file.RasterYSize
    tropo_rsc_dict['XMIN'] = 0
    tropo_rsc_dict['XMAX'] = in_file.RasterXSize
    tropo_rsc_dict['YMIN'] = 0
    tropo_rsc_dict['YMAX'] = in_file.RasterYSize
    tropo_rsc_dict['X_FIRST'] = in_file.GetGeoTransform()[0]
    tropo_rsc_dict['Y_FIRST'] = in_file.GetGeoTransform()[3]
    tropo_rsc_dict['X_STEP'] = in_file.GetGeoTransform()[1]
    tropo_rsc_dict['Y_STEP'] = in_file.GetGeoTransform()[-1]
    tropo_rsc_dict['X_UNIT'] = 'degres'
    tropo_rsc_dict['Y_UNIT'] = 'degres'
    tropo_rsc_dict['Z_OFFSET'] = 0
    tropo_rsc_dict['Z_SCALE'] = 1
    tropo_rsc_dict['PROJECTION'] = 'LATLON'
    tropo_rsc_dict['DATUM'] = 'WGS84'
    tropo_rsc_dict['TIME_OF_DAY'] = 'NoneUTC'
    return tropo_rsc_dict


# Perform initial layer, product, and correction sanity checks
def layerCheck(
        products, layers, nc_version, gacos_products, tropo_models,
        extract_or_ts):
    """Check if any conflicts between netcdf versions and expected layers."""
    # track if product stack is NISAR GUNW or not
    nisar_file = False
    track_fileext = products[0]['unwrappedPhase'][0]
    if len(track_fileext.split('.h5')) > 1:
        nisar_file = True

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
                LOGGER.warning(
                    f'User-requested tropo model {i} will not be generated as '
                    'it does not exist in any of the input products')
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
                set(RAIDER_TROPO_LAYERS).issubset(all_valid_layers):
            tropo_total = True

    # differentiate between extract and TS pipeline
    # extract pipeline
    if extract_or_ts == 'extract':
        if not layers and not gacos_products:
            LOGGER.info(
                'No layers specified; only creating bounding box shapes')
            return [], [], []

        elif gacos_products:
            LOGGER.info(
                'Tropospheric corrections will be applied, making sure at '
                'least unwrappedPhase and incidenceAngle are extracted.')

            # If no input layers specified, initialize list
            if not layers:
                layers = []

            if 'incidenceAngle' not in layers:
                layers.append('incidenceAngle')

            if 'unwrappedPhase' not in layers:
                layers.append('unwrappedPhase')
        else:
            layers = [i.replace(' ', '') for i in layers]

    # TS pipeline
    TS_LAYERS_DUP = ['unwrappedPhase', 'coherence', 'incidenceAngle',
                     'lookAngle', 'azimuthAngle', 'bPerpendicular']

    ts_defaults = ['ionosphere', 'solidEarthTide']
    if extract_or_ts == 'tssetup':
        if layers:
            # remove layers already generated in default TS workflow
            layers = [i for i in layers if i not in TS_LAYERS_DUP]

        else:
            layers = []

        # add additional layers for default workflow
        ts_defaults = list(
            set.intersection(*map(set, [ts_defaults, all_valid_layers])))

        if ts_defaults != []:
            tropo_total = True
            for i in ts_defaults:
                if i not in layers:
                    layers.append(i)

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

        if model_names == [] and not nisar_file:
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
