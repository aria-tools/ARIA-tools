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
