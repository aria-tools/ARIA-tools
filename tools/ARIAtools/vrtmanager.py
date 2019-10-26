#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from osgeo import gdal
gdal.UseExceptions()
# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

###Save file with gdal
def renderVRT(fname, data_lyr, geotrans=None, drivername='ENVI', gdal_fmt='float32', proj=None, nodata=None, verbose=False):
    '''
        Exports raster and renders corresponding VRT file.
    '''
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
    '''
        Generate VRT of shapefile unions.
    '''

    # Import functions
    import os

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
def resampleRaster(fname, multilooking, bounds, prods_TOTbbox, outputFormat='ENVI', num_threads='2'):
    '''
        Resample rasters and update corresponding VRTs.
    '''

    # Import functions
    import os
    from scipy import stats
    import numpy as np

    if outputFormat=='VRT':
        fname+='.vrt'

    # Access original shape
    ds=gdal.Warp('', fname, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
    arrshape=[abs(ds.GetGeoTransform()[1])*multilooking, abs(ds.GetGeoTransform()[-1])*multilooking] # Get output res

    # Must resample mask with nearest-neighbor
    if fname.split('/')[-2]=='mask':
        # Resample raster
        gdal.Warp(fname, fname, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='near',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))

    # Use pixel function to downsample connected components/unw files based off of frequency of connected components in each window
    elif fname.split('/')[-2]=='connectedComponents' or fname.split('/')[-2]=='unwrappedPhase':
        #open connected components/unw files
        fnameunw=os.path.join('/'.join(fname.split('/')[:-2]),'unwrappedPhase',''.join(fname.split('/')[-1]))
        fnameconcomp=os.path.join('/'.join(fname.split('/')[:-2]),'connectedComponents',''.join(fname.split('/')[-1]))
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

        #get geotrans/proj
        ds=gdal.Warp('', fnameconcomp, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='near',multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
        geotrans=ds.GetGeoTransform() ; proj=ds.GetProjection()

        unwmap=[]
        concompmap=[]
        for row in range(multilooking,ds_unw.shape[0]+multilooking,multilooking):
            unwmap_row=[]
            concompmap_row=[]
            for column in range(multilooking,ds_unw.shape[1]+multilooking,multilooking):
                #get subset values
                subset_concomp = ds_concomp[row-multilooking:row+multilooking,column-multilooking:column+multilooking]
                subset_unw = ds_unw[row-multilooking:row+multilooking,column-multilooking:column+multilooking]
                concomp_mode = stats.mode(subset_concomp.flatten()).mode[0]
                #average only phase values coinciding with concomp mode
                subset_concomp = np.where(subset_concomp != concomp_mode, 0, 1)
                subset_unw=subset_unw*subset_concomp
                #assign downsampled pixel values
                unwmap_row.append(subset_unw.mean())
                concompmap_row.append(concomp_mode)
            unwmap.append(unwmap_row)
            concompmap.append(concompmap_row)

        #finalize arrays
        unwmap=np.array(unwmap) 
        unwmap=np.ma.masked_invalid(unwmap) ; np.ma.set_fill_value(unwmap, ds_unw_nodata)
        concompmap=np.array(concompmap)
        concompmap=np.ma.masked_where(concompmap == ds_concomp_nodata, concompmap) ; np.ma.set_fill_value(concompmap, ds_concomp_nodata)

        #unwphase
        renderVRT(fnameunw, unwmap.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=ds_unw_nodata)
        #conn comp, mode resampling through gdal
        renderVRT(fnameconcomp, concompmap.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='int16', proj=proj, nodata=ds_concomp_nodata)

    # Resample all other files with lanczos
    else:
        # Resample raster
        gdal.Warp(fname, fname, options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrshape[0], yRes=arrshape[1], resampleAlg='lanczos',multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))

    if outputFormat!='VRT':
        # update VRT
        gdal.BuildVRT(fname+'.vrt', fname, options=gdal.BuildVRTOptions(options=['-overwrite']))

    return
