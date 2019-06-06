#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, & Emre Havazli
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import sys
import glob
import numpy as np
from datetime import datetime, time
from osgeo import gdal
gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')
# Import functions
from ARIAProduct import ARIA_standardproduct
from shapefile_util import open_shapefile
from extractProduct import merged_productbbox
from extractProduct import prep_dem
from extractProduct import export_products
from extractProduct import finalize_metadata


def createParser():
    '''
        Extract unwrapped interferogram, coherence, ⊥ baseline, and LOS file(s) for TS analysis.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Get DEM')
    parser.add_argument('-f', '--file', dest='imgfile', type=str, required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str, default='download', help='DEM file. Default is to download new DEM.')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str, help='projection for DEM. By default WGS84.')
    parser.add_argument('-bpx', '--bperpextract', action='store_true', dest='bperpextract', help="If turned on, extracts perpendicular baseline grids. A single perpendicular baseline value is calculated and included in the metadata stack VRT cube.")
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Provide valid mask file.")
    parser.add_argument('-s', '--stack', action='store_true', dest='stack', help="If turned on, creates vrt files of stacks that can be used for time series processing")
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def extractBaselineDict(aria_prod):
    bPerp = {}
    for i in aria_prod.products[1]:
        baselineName = i['bPerpendicular'][0]
        ds = gdal.Open(baselineName)
        # return [min, max, mean, std]
        stat  = ds.GetRasterBand(1).GetStatistics(True, True)
        bPerp[i['pair_name'][0]] = stat[2]
        ds = None

    return bPerp

def extractUTCtime(aria_prod):
    f = list(np.sort(aria_prod.files))
    utcDict = {}
    dates= []
    times = []
    for i in range(0,len(f)):
        if f[i].split('-')[6] in dates:
            times.append(f[i].split('-')[7])
            utcDict[dates[-1]] = times
        else:
            times = []
            dates.append(f[i].split('-')[6])
            times.append(f[i].split('-')[7])
            utcDict[f[i].split('-')[6]] = times

    UTC_time = {}
    for key, value in utcDict.items():
        tFirst = datetime.strptime(value[0],'%H%M%S')
        tLast = datetime.strptime(value[-1],'%H%M%S')
        tDelta = (tLast - tFirst)/2
        tMid = tFirst+tDelta
        UTC_time[key] = str(datetime.strftime(tMid,'%H%M%S'))

    return UTC_time

def generateStack(aria_prod,inputFiles,outputFileName,workdir='./'):

    UTC_time = extractUTCtime(aria_prod)
    bPerp = extractBaselineDict(aria_prod)

    ###Set up single stack file
    if not os.path.exists(os.path.join(workdir,'stack')):
        print('Creating directory: {0}'.format(os.path.join(workdir,'stack')))
        os.makedirs(os.path.join(workdir,'stack'))
    else:
        print('Directory {0} already exists.'.format(os.path.join(workdir,'stack')))

    if inputFiles in ['unwrappedPhase', 'unwrapped', 'unw']:
        domainName = 'unwrappedPhase'
        intList = glob.glob(os.path.join(workdir,'unwrappedPhase','[0-9]*[0-9].vrt'))
        dataType = "Float32"
        print('Number of unwrapped interferograms discovered: ', len(intList))
        Dlist = intList
    elif inputFiles in ['coherence', 'Coherence', 'coh']:
        domainName = 'Coherence'
        cohList = glob.glob(os.path.join(workdir,'coherence','[0-9]*[0-9].vrt'))
        dataType = "Float32"
        print('Number of coherence discovered: ', len(cohList))
        Dlist = cohList
    elif inputFiles in ['connectedComponents','connectedComponent','connComp']:
        domainName = 'connectedComponents'
        connCompList = glob.glob(os.path.join(workdir,'connectedComponents','[0-9]*[0-9].vrt'))
        dataType = "Int16"
        print('Number of connectedComponents discovered: ', len(connCompList))
        Dlist = connCompList
    else:
        print('Stacks can be created for unwrapped interferogram, coherence and connectedComponent VRT files')

    for ind, data in enumerate(Dlist):
        width = None
        height = None

        ds = gdal.Open(data, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        height = ds.RasterYSize
        gt  = ds.GetGeoTransform()
        projection = ds.GetProjection()
        ds = None

    # setting up a subset of the stack
    ymin, ymax, xmin, xmax = [0 , height, 0 , width]

    xsize = xmax - xmin
    ysize = ymax - ymin

    # extraction of radar meta-data
    wavelength = aria_prod.products[0][0]['wavelength'][0].data
    startRange = aria_prod.products[0][0]['slantRangeStart'][0].data
    endRange = aria_prod.products[0][0]['slantRangeEnd'][0].data
    rangeSpacing = aria_prod.products[0][0]['slantRangeSpacing'][0].data

    with open( os.path.join(workdir,'stack', (outputFileName+'.vrt')), 'w') as fid:
        fid.write( '''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
        <SRS>{proj}</SRS>">
        <GeoTransform>{GT0},{GT1},{GT2},{GT3},{GT4},{GT5}</GeoTransform>
        '''.format(xsize=xsize, ysize=ysize,
        proj=projection,
        GT0=gt[0],GT1=gt[1],GT2=gt[2],GT3=gt[3],GT4=gt[4],GT5=gt[5]))

        for ind, data in enumerate(Dlist):
            metadata = {}
            dates = data.split('/')[2][:-4]
            width = None
            height = None
            path = None

            ds = gdal.Open(data, gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize
            ds = None

            metadata['wavelength'] = wavelength
            metadata['utcTime'] = UTC_time[dates]
            metadata['bPerp'] = bPerp[dates]

            path = os.path.abspath(data)

            outstr = '''    <VRTRasterBand dataType="{dataType}" band="{index}">
        <SimpleSource>
            <SourceFilename>{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="{dataType}"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain='{domainName}'>
            <MDI key="Dates">"{dates}"</MDI>
            <MDI key="Wavelength (m)">"{wvl}"</MDI>
            <MDI key="UTCTime">"{acq}"</MDI>
            <MDI key="bPerp">"{bPerp}"</MDI>
            <MDI key="startRange">"{start_range}"</MDI>
            <MDI key="endRange">"{end_range}"</MDI>
            <MDI key="rangeSpacing">"{range_spacing}"</MDI>
        </Metadata>
    </VRTRasterBand>\n'''.format(domainName=domainName,width=width, height=height,
                                xmin=xmin, ymin=ymin,
                                xsize=xsize, ysize=ysize,
                                dates=dates, acq=metadata['utcTime'],
                                wvl=metadata['wavelength'], index=ind+1,
                                path=path, dataType=dataType, bPerp=metadata['bPerp'],
                                start_range=startRange, end_range=endRange, range_spacing=rangeSpacing)
            fid.write(outstr)

        fid.write('</VRTDataset>')
        print(outputFileName, ': stack generated')


if __name__ == '__main__':
    '''
        Main driver.
    '''
    inps = cmdLineParse()

    print("***Time-series Preparation Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)

    # Load mask (if specified).
    if inps.mask is not None:
        inps.mask=gdal.Warp('', inps.mask, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds, dstNodata=0))
        inps.mask.SetProjection(proj)
        # If no data value
        if inps.mask.GetRasterBand(1).GetNoDataValue():
            inps.mask=np.ma.masked_where(inps.mask.ReadAsArray() == inps.mask.GetRasterBand(1).GetNoDataValue(), inps.mask.ReadAsArray())
        else:
            inps.mask=inps.mask.ReadAsArray()


    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        print('Download/cropping DEM')
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, proj, arrshape=arrshape, workdir=inps.workdir)

    # Only extract layers needed for TS analysis
    layers=['unwrappedPhase','coherence']
    print('Extracting unwrapped phase, coherence, perpendicular baseline and connected components for each interferogram pair')
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir)

    layers=['incidenceAngle','lookAngle','azimuthAngle']
    print('Extracting incidence angle, look angle, and heading of the first interferogram only')
    export_products([standardproduct_info.products[1][0]], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir) ##Only the first product is written

    if inps.bperpextract==True:
        layers=['bPerpendicular']
        print('Extracting perpendicular baseline grids for each interferogram pair')
        export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir)

    if inps.stack==True:
        generateStack(standardproduct_info,'unwrappedPhase','unwrapStack')
        generateStack(standardproduct_info,'coherence','cohStack')
        generateStack(standardproduct_info,'connectedComponents','connCompStack')
