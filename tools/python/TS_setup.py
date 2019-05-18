#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
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

_world_dem = ('https://cloud.sdsc.edu/v1/AUTH_opentopography/Raster/SRTM_GL1_Ellip/SRTM_GL1_Ellip_srtm.vrt')

def createParser():
    '''
        Extract unwrapped interferogram, coherence, ⊥ baseline, and LOS file(s) for TS analysis.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Get DEM')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
            default=None, help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str,
            help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Provide valid mask file.")
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-s', '--stack', action='store_true', dest='stack', help="If turned on, creates vrt files of stacks that can be used for time series processing")
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


if __name__ == '__main__':
    '''
        Main driver.
    '''
    inps = cmdLineParse()

    print("########################################")
    print("class 'Aria_standardproduct': sort input file(s) by starting time")
    print("if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected."+'\n')
    print("Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product + path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)."+'\n')
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    print('\n'+'\n'+"########################################")
    print("fn 'merged_productbbox': extract/merge productBoundingBox layers for each pair and update dict, report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM."+'\n')
    standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)

    wavelength = standardproduct_info.products[0][0]['wavelength'][0].data

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
        print('\n'+'\n'+"########################################")
        print("fn 'prep_dem': Pass DEM-filename, loaded DEM array, and lat/lon arrays."+'\n')
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, proj, arrshape=arrshape, workdir=inps.workdir)


    # Only extract layers needed for TS analysis
    layers=['unwrappedPhase','coherence']
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir)

    # layers=['bPerpendicular','incidenceAngle','lookAngle']
    # export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir)


    if inps.stack==True:
        ##Mid-frame UTC times to dictionary
        f = list(np.sort(standardproduct_info.files))
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
            UTC_time[key] = datetime.strftime(tMid,'%H%M%S')

    ####Set up single stack file
        if not os.path.exists(os.path.join(inps.workdir,'stack')):
            print('Creating directory: {0}'.format(os.path.join(inps.workdir,'stack')))
            os.makedirs(os.path.join(inps.workdir,'stack'))
        else:
            print('Directory {0} already exists.'.format(os.path.join(inps.workdir,'stack')))

        intlist = glob.glob(os.path.join(inps.workdir,'unwrappedPhase','*.vrt'))
        print('Number of interferograms discovered: ', len(intlist))
        cohlist = glob.glob(os.path.join(inps.workdir,'coherence','*.vrt'))
        print('Number of coherence discovered: ', len(cohlist))


        ####Interferograms####
        for ind, interf in enumerate(intlist):
            width = None
            height = None
            path = None

            ds = gdal.Open(interf, gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize

            ds = None

        # setting up a subset of the stack
        ymin, ymax, xmin, xmax = [0 , height, 0 , width]
        if inps.bbox:
            ymin, ymax, xmin, xmax = inps.bbox

        xsize = xmax - xmin
        ysize = ymax - ymin

        with open( os.path.join(inps.workdir,'stack', 'stack_intrf.vrt'), 'w') as fid:
            fid.write( '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(xsize=xsize, ysize=ysize))

            for ind, interf in enumerate(intlist):
                metadata = {}
                dates = interf.split('/')[2][:-4]
                width = None
                height = None
                path = None

                ds = gdal.Open(interf, gdal.GA_ReadOnly)
                width = ds.RasterXSize
                height = ds.RasterYSize
                ds = None


                metadata['WAVELENGTH'] = wavelength
                # metadata['ACQUISITION_TIME'] = os.path.basename(os.path.dirname(slc))
                metadata['UTC Time'] = UTC_time[dates]

                path = os.path.abspath(interf)

                outstr = '''    <VRTRasterBand dataType="CFloat32" band="{index}">
            <SimpleSource>
                <SourceFilename>{path}</SourceFilename>
                <SourceBand>1</SourceBand>
                <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="CFloat32"/>
                <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
                <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
            </SimpleSource>
            <Metadata domain="interferogram">
                <MDI key="Dates">{dates}</MDI>
                <MDI key="Wavelength (m)">{wvl}</MDI>
                <MDI key="UTC Time">{acq}</MDI>
            </Metadata>
        </VRTRasterBand>\n'''.format(width=width, height=height,
                                    xmin=xmin, ymin=ymin,
                                    xsize=xsize, ysize=ysize,
                                    dates=dates, acq=metadata['UTC Time'],
                                    wvl = metadata['WAVELENGTH'], index=ind+1,
                                    path = path)
                fid.write(outstr)

            fid.write('</VRTDataset>')
            print('Interferogram stack generated')


            ####Coherence####
        for ind, coh in enumerate(cohlist):
            width = None
            height = None
            path = None

            ds = gdal.Open(coh, gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize
            ds = None

        # setting up a subset of the stack
        ymin, ymax, xmin, xmax = [0 , height, 0 , width]
        if inps.bbox:
            ymin, ymax, xmin, xmax = inps.bbox

        xsize = xmax - xmin
        ysize = ymax - ymin

        with open( os.path.join(inps.workdir,'stack', 'stack_coh.vrt'), 'w') as fid:
            fid.write( '<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">\n'.format(xsize=xsize, ysize=ysize))

            for ind, coh in enumerate(cohlist):
                metadata = {}
                dates = coh.split('/')[2][:-4]
                width = None
                height = None
                path = None

                ds = gdal.Open(coh, gdal.GA_ReadOnly)
                width = ds.RasterXSize
                height = ds.RasterYSize
                ds = None


                metadata['WAVELENGTH'] = wavelength
                # metadata['ACQUISITION_TIME'] = os.path.basename(os.path.dirname(slc))
                metadata['UTC Time'] = UTC_time[dates]

                path = os.path.abspath(coh)

                outstr = '''    <VRTRasterBand dataType="CFloat32" band="{index}">
            <SimpleSource>
                <SourceFilename>{path}</SourceFilename>
                <SourceBand>1</SourceBand>
                <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="CFloat32"/>
                <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
                <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
            </SimpleSource>
            <Metadata domain="coherence">
                <MDI key="Dates">{dates}</MDI>
                <MDI key="Wavelength (m)">{wvl}</MDI>
                <MDI key="UTC Time">{acq}</MDI>
            </Metadata>
        </VRTRasterBand>\n'''.format(width=width, height=height,
                                    xmin=xmin, ymin=ymin,
                                    xsize=xsize, ysize=ysize,
                                    dates=dates, acq=metadata['UTC Time'],
                                    wvl = metadata['WAVELENGTH'], index=ind+1,
                                    path = path)
                fid.write(outstr)

            fid.write('</VRTDataset>')
            print('Coherence stack generated')
