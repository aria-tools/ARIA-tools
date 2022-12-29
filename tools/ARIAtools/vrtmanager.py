#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import glob
import numpy as np
from osgeo import gdal
import logging
gdal.UseExceptions()

# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

log = logging.getLogger(__name__)

###Save file with gdal
def renderVRT(fname, data_lyr, geotrans=None, drivername='ENVI', gdal_fmt='float32', proj=None, nodata=None, verbose=False):
    """Exports raster and renders corresponding VRT file."""
    gdalMap = { 'byte'   : 1,
                'int16'  : 3,
                'int32'    : 5,
                'float32'  : 6,
                'float64' : 7,
                'cfloat32' : 10,
                'cfloat64': 11}

    gdalfile=gdal.GetDriverByName(drivername).Create(fname, data_lyr.shape[1], data_lyr.shape[0], 1, gdalMap[gdal_fmt])
    gdalfile.GetRasterBand(1).WriteArray(data_lyr)
    if geotrans: #If user wishes to update geotrans.
        gdalfile.SetGeoTransform(geotrans)
    if proj: #If user wishes to update projection.
        gdalfile.SetProjection(proj)
    if nodata is not None: #If user wishes to set nodata val.
        gdalfile.GetRasterBand(1).SetNoDataValue(nodata)
        # Finalize VRT
        gdal.Translate(fname+'.vrt', gdalfile, options=gdal.TranslateOptions(format="VRT", noData=nodata))
    else:
        # Finalize VRT
        gdal.Translate(fname+'.vrt', gdalfile, options=gdal.TranslateOptions(format="VRT"))

    gdalfile = None

    return


###Make OGR VRT file
def renderOGRVRT(vrt_filename, src_datasets):
    """Generate VRT of shapefile unions."""
    vrt_head='<OGRVRTDataSource>\n  <OGRVRTUnionLayer name="merged">\n    <FieldStrategy>Union</FieldStrategy>\n'
    vrt_tail='  </OGRVRTUnionLayer>\n</OGRVRTDataSource>\n'

    with open(vrt_filename,'w') as vrt_write:
        vrt_write.write(vrt_head)
        for i in enumerate(src_datasets):
            vrt_write.write('    <OGRVRTLayer name="Dataset%i_%s">\n'%(i[0],os.path.basename(i[1]).split('.shp')[0]))
            vrt_write.write('      <SrcDataSource shared="1">%s</SrcDataSource>\n'%(i[1]))
            vrt_write.write('      <SrcLayer>%s</SrcLayer>\n'%(os.path.basename(i[1]).split('.shp')[0]))
            vrt_write.write('    </OGRVRTLayer>\n')
        vrt_write.write(vrt_tail)

    return


###Resample raster
def resampleRaster(fname, multilooking, bounds, prods_TOTbbox, rankedResampling=False, outputFormat='ENVI', num_threads='2'):
    """Resample rasters and update corresponding VRTs."""
    # Import functions
    from scipy import stats
    from decimal import Decimal, ROUND_HALF_UP

    # Check if physical raster exists and needs to be updated
    # Also get datasource name (inputname)
    if outputFormat=='VRT' and os.path.exists(fname.split('.vrt')[0]):
        outputFormat='ENVI'
    if os.path.exists(fname.split('.vrt')[0]):
        inputname=fname
    else:
        fname+='.vrt'
        inputname=gdal.Open(fname).GetFileList()[-1]

    # Access original shape
    ds=gdal.Warp('', fname, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
    arrshape=[abs(ds.GetGeoTransform()[1])*multilooking, abs(ds.GetGeoTransform()[-1])*multilooking] # Get output res
    #get geotrans/proj
    ds=gdal.Warp('', fname, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='near',multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
    geotrans=ds.GetGeoTransform() ; proj=ds.GetProjection()

    # Must resample mask with nearest-neighbor
    if fname.split('/')[-2]=='mask':
        # Resample raster
        gdal.Warp(fname, inputname, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='near',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))

    # Use pixel function to downsample connected components/unw files based off of frequency of connected components in each window
    elif fname.split('/')[-2]=='connectedComponents' or fname.split('/')[-2]=='unwrappedPhase':
        # Resample unw phase based off of mode of connected components
        fnameunw=os.path.join('/'.join(fname.split('/')[:-2]),'unwrappedPhase',''.join(fname.split('/')[-1]).split('.vrt')[0])
        fnameconcomp=os.path.join('/'.join(fname.split('/')[:-2]),'connectedComponents',''.join(fname.split('/')[-1]).split('.vrt')[0])
        if rankedResampling:
            #open connected components/unw files
            ds_concomp=gdal.Open(fnameconcomp)
            ds_concomp_nodata=ds_concomp.GetRasterBand(1).GetNoDataValue()
            ds_concomp=ds_concomp.ReadAsArray()
            ds_concomp=np.ma.masked_where(ds_concomp == ds_concomp_nodata, ds_concomp)
            np.ma.set_fill_value(ds_concomp, ds_concomp_nodata)

            ds_unw=gdal.Open(fnameunw)
            ds_unw_nodata=ds_unw.GetRasterBand(1).GetNoDataValue()
            ds_unw=ds_unw.ReadAsArray()
            ds_unw=np.ma.masked_where(ds_unw == ds_unw_nodata, ds_unw)
            np.ma.set_fill_value(ds_unw, ds_unw_nodata)

            unwmap=[]
            for row in range(multilooking,(ds_unw.shape[0])+multilooking,multilooking):
                unwmap_row=[]
                for column in range(multilooking,(ds_unw.shape[1])+multilooking,multilooking):
                    #get subset values
                    subset_concomp = ds_concomp[row-multilooking:row,column-multilooking:column]
                    subset_unw = ds_unw[row-multilooking:row,column-multilooking:column]
                    concomp_mode = stats.mode(subset_concomp.flatten()).mode[0]
                    #average only phase values coinciding with concomp mode
                    subset_concomp = np.where(subset_concomp != concomp_mode, 0, 1)
                    subset_unw=subset_unw*subset_concomp
                    #assign downsampled pixel values
                    unwmap_row.append(subset_unw.mean())
                unwmap.append(unwmap_row)

            #finalize unw array
            unwmap=np.array(unwmap)
            #finalize unw array shape
            unwmap=unwmap[0:int(Decimal(ds_unw.shape[0]/multilooking).quantize(0, ROUND_HALF_UP)), \
                0:int(Decimal(ds_unw.shape[1]/multilooking).quantize(0, ROUND_HALF_UP))]
            unwmap=np.ma.masked_invalid(unwmap) ; np.ma.set_fill_value(unwmap, ds_unw_nodata)
            #unwphase
            renderVRT(fnameunw, unwmap.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=ds_unw_nodata)
            #temp workaround for gdal bug
            try:
                gdal.Open(fnameunw)
            except RuntimeError:
                for f in glob.glob(fnameunw+"*"):
                    os.remove(f)
                unwmap[0,0]=unwmap[0,0]-1e-6
                renderVRT(fnameunw, unwmap.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=ds_unw_nodata)
            del unwmap
            # Resample connected components
            gdal.Warp(fnameconcomp, fnameconcomp, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='mode',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))
            # update VRT
            gdal.BuildVRT(fnameconcomp+'.vrt', fnameconcomp, options=gdal.BuildVRTOptions(options=['-overwrite']))
        # Default: resample unw phase with gdal average algorithm
        else:
            # Resample unwphase
            ds_unw_nodata=gdal.Open(fnameunw)
            ds_unw_nodata=ds_unw_nodata.GetRasterBand(1).GetNoDataValue()
            gdal.Warp(fnameunw, fnameunw, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='average',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))
            # update VRT
            gdal.BuildVRT(fnameunw+'.vrt', fnameunw, options=gdal.BuildVRTOptions(options=['-overwrite']))
            #temp workaround for gdal bug
            try:
                gdal.Open(fnameunw)
            except RuntimeError:
                unwmap=np.fromfile(fnameunw,dtype=np.float32).reshape(ds.GetRasterBand(1).ReadAsArray().shape)
                for f in glob.glob(fnameunw+"*"):
                    os.remove(f)
                unwmap[0,0]=unwmap[0,0]-1e-6
                renderVRT(fnameunw, unwmap, geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=ds_unw_nodata)
                del unwmap
            # Resample connected components
            gdal.Warp(fnameconcomp, fnameconcomp, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='near',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))
            # update VRT
            gdal.BuildVRT(fnameconcomp+'.vrt', fnameconcomp, options=gdal.BuildVRTOptions(options=['-overwrite']))

    # Resample all other files with lanczos
    else:
        # Resample raster
        gdal.Warp(fname, inputname, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='lanczos',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))

    if outputFormat!='VRT':
        # update VRT
        gdal.BuildVRT(fname+'.vrt', fname, options=gdal.BuildVRTOptions(options=['-overwrite']))

    return


###Average rasters
def rasterAverage(outname, product_dict, bounds, prods_TOTbbox, outputFormat='ENVI', thresh=None):
    """Generate average of rasters.

    Currently implemented for:
    1. amplitude under 'mask_util.py'
    2. coherence under 'plot_avgcoherence' function of productPlot"""
    ###Make average raster
    # Delete existing average raster file
    for i in glob.glob(outname+'*'): os.remove(i)

    # Iterate through all layers
    for i in enumerate(product_dict):
        arr_file = gdal.Warp('', i[1], options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds))
        arr_file_arr = np.ma.masked_where(arr_file.ReadAsArray() == arr_file.GetRasterBand(1).GetNoDataValue(), arr_file.ReadAsArray())
        ### Iteratively update average raster file
        if os.path.exists(outname):
            arr_file = gdal.Open(outname,gdal.GA_Update)
            arr_file = arr_file.GetRasterBand(1).WriteArray(arr_file_arr+arr_file.ReadAsArray())
        else:
            # If looping through first raster file, nothing to sum so just save to file
            renderVRT(outname, arr_file_arr, geotrans=arr_file.GetGeoTransform(), drivername=outputFormat, \
               gdal_fmt=arr_file_arr.dtype.name, proj=arr_file.GetProjection(), \
               nodata=arr_file.GetRasterBand(1).GetNoDataValue())
        arr_file = arr_file_arr = None

    # Take average of raster sum
    arr_file = gdal.Open(outname,gdal.GA_Update)
    arr_mean = arr_file.ReadAsArray()/len(product_dict)

    # Mask using specified raster threshold
    if thresh:
        arr_mean = np.where(arr_mean < float(thresh), 0, 1)

    # Save updated array to file
    arr_file = arr_file.GetRasterBand(1).WriteArray(arr_mean)
    arr_file = None ; arr_mean = None

    # Load raster to pass
    arr_file = gdal.Open(outname).ReadAsArray()

    return arr_file


###Generate GACOS rsc
def rscGacos(outvrt, merged_rsc, tropo_date_dict):
    """Generate GACOS product .rsc files"""
    in_file = gdal.Open(outvrt)
    date = os.path.basename(outvrt[:-4]).split('.tif')[0].split('.ztd')[0]
    # Create merge rsc file
    with open(outvrt[:-4]+'.rsc','w') as merged_rsc:
        merged_rsc.write('WIDTH %s\n'%(in_file.RasterXSize))
        merged_rsc.write('FILE_LENGTH %s\n'%(in_file.RasterYSize))
        merged_rsc.write('XMIN %s\n'%(0))
        merged_rsc.write('XMAX %s\n'%(in_file.RasterXSize))
        merged_rsc.write('YMIN %s\n'%(0))
        merged_rsc.write('YMAX %s\n'%(in_file.RasterYSize))
        merged_rsc.write('X_FIRST %f\n'%(in_file.GetGeoTransform()[0]))
        merged_rsc.write('Y_FIRST %f\n'%(in_file.GetGeoTransform()[3]))
        merged_rsc.write('X_STEP %f\n'%(in_file.GetGeoTransform()[1]))
        merged_rsc.write('Y_STEP %f\n'%(in_file.GetGeoTransform()[-1]))
        merged_rsc.write('X_UNIT %s\n'%('degres'))
        merged_rsc.write('Y_UNIT %s\n'%('degres'))
        merged_rsc.write('Z_OFFSET %s\n'%(0))
        merged_rsc.write('Z_SCALE %s\n'%(1))
        merged_rsc.write('PROJECTION %s\n'%('LATLON'))
        merged_rsc.write('DATUM %s\n'%('WGS84'))
        merged_rsc.write('TIME_OF_DAY %s\n'%(''.join( \
                                             tropo_date_dict[date+"_UTC"])))

    return


###Parse GACOS metadata from TIF
def tifGacos(intif):
    """Pass GACOS metadata as dict from TIF"""
    tropo_rsc_dict={}
    in_file = gdal.Open(intif)
    #Populate dict
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


###Perform initial layer, product, and correction sanity checks
def layerCheck(products, layers, nc_version, gacos_products, extract_or_ts):
    """Check if any conflicts between netcdf versions and expected layers."""

    # Ignore productBoundingBoxes & pair-names, they are not raster layers
    ignore_lyr = ['productBoundingBox','productBoundingBoxFrames','pair_name']
    # Check all available layers in stack
    products = [list(i.keys()) for i in products]
    products = [[sub for sub in i if sub not in ignore_lyr] for i in products]
    all_valid_layers = list(set.intersection(*map(set, products)))

    # If specified, extract all layers
    if layers.lower()=='all':
        log.info('All layers are to be extracted, pass all keys.')
        layers = all_valid_layers

    # differentiate between extract and TS pipeline
    # extract pipeline
    if extract_or_ts == 'extract':
        if not layers and not gacos_products:
            log.info('No layers specified; only creating bounding box shapes')
        elif gacos_products:
            log.info('Tropospheric corrections will be applied, making sure '
                     'at least unwrappedPhase and lookAngle are extracted.')
            # If no input layers specified, initialize list
            if not layers:
                layers = []
            # If valid argument for input layers passed, parse to list
            if isinstance(layers, str):
                layers = list(layers.split(','))
                layers = [i.replace(' ','') for i in layers]
            if 'lookAngle' not in layers:
                layers.append('lookAngle')
            if 'unwrappedPhase' not in layers:
                layers.append('unwrappedPhase')
        else:
            layers = list(layers.split(','))
            layers = [i.replace(' ','') for i in layers]

    # TS pipeline
    ts_layers_dup = ['unwrappedPhase', 'coherence', 'incidenceAngle',
                     'lookAngle', 'azimuthAngle', 'bPerpendicular']
    if extract_or_ts == 'tssetup':
        if layers:
            if isinstance(layers, str):
                layers = list(layers.split(','))
            # remove layers already generated in default TS workflow
            layers = [i for i in layers if i not in ts_layers_dup]


    # Check to see if internal conflict between tropo correction methods
    raider_tropo_layers = ['troposphereWet', 'troposphereHydrostatic', \
                           'troposphere']
    user_tropo_layers = list(set.intersection(*map(set, \
                             [layers, raider_tropo_layers])))
    if gacos_products and user_tropo_layers != []:
        raise Exception('User specified extraction of raider-derived '
                        'tropo layers %s AND gacos products with "-tp %s". '
                        'Proceed with only one.'%(user_tropo_layers,
                         gacos_products))

    # pass intersection of valid layers and track invalid requests
    layer_reject = list(set.symmetric_difference(*map(set, \
                        [all_valid_layers, layers])))
    layer_reject = list(set.intersection(*map(set, \
                        [layer_reject, raider_tropo_layers])))
    # only report layers which user requested
    layer_reject = list(set.intersection(*map(set, \
                        [layer_reject, layers])))
    layers = list(set.intersection(*map(set, [layers, all_valid_layers])))
    if layer_reject != []:
        log.warning('User-requested layers %s cannot be extracted as they '
                    'are not common to all products. Consider fixing input '
                    '"-nc_version %s" constraint to filter older product '
                    'variants \n', layer_reject, nc_version)

    return layers
