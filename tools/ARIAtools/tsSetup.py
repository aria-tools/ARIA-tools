#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, & Emre Havazli
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import glob
from datetime import datetime
from osgeo import gdal

gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Import functions
from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.mask_util import prep_mask
from ARIAtools.extractProduct import merged_productbbox, prep_dem, export_products, tropo_correction


def createParser():
    '''
        Extract unwrapped interferogram, coherence, ⊥ baseline, and LOS file(s) for TS analysis.
    '''

    import argparse
    parser = argparse.ArgumentParser(description='Get DEM')
    parser.add_argument('-f', '--file', dest='imgfile', type=str, required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-tp', '--tropo_products', dest='tropo_products', type=str, default=None, help='Path to director(ies) or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None, help='Specify layers to extract as a comma deliminated list bounded by single quotes. Allowed keys are: "unwrappedPhase", "coherence", "amplitude", "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", "azimuthAngle", "ionosphere". If "all" is specified, then all layers are extracted. If blank, will only extract bounding box.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str, default='download', help='DEM file. Default is to download new DEM.')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str, help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Path to mask file or 'Download'. File needs to be GDAL compatabile, contain spatial reference information, and have invalid/valid data represented by 0/1, respectively. If 'Download', will use GSHHS water mask. If 'NLCD', will mask classes 11, 12, 90, 95; see: https://www.mrlc.gov/national-land-cover-database-nlcd-201://www.mrlc.gov/national-land-cover-database-nlcd-2016")
    parser.add_argument('-at', '--amp_thresh', dest='amp_thresh', default=None, type=str, help='Amplitude threshold below which to mask. Specify "None" to not use amplitude mask. By default "None".')
    parser.add_argument('-nt', '--num_threads', dest='num_threads', default='2', type=str, help='Specify number of threads for multiprocessing operation in gdal. By default "2". Can also specify "All" to use all available threads.')
#    parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',  type=str, default='overlap', help="Method applied to stitch the unwrapped data. Either 'overlap', where product overlap is minimized, or '2stage', where minimization is done on connected components, are allowed methods. Default is 'overlap'.")
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str, default='VRT', help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By default files are generated virtually except for "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", "azimuthAngle", "unwrappedPhase" as these are require either DEM intersection or corrections to be applied')
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-bp', '--bperp', action='store_true', dest='bperp', help="If turned on, extracts perpendicular baseline grids. Default: A single perpendicular baseline value is calculated and included in the metadata of stack cubes for each pair.")
    parser.add_argument('-ml', '--multilooking', dest='multilooking', type=int, default=None, help='Multilooking factor is an integer multiple of standard resolution. E.g. 2 = 90m*2 = 180m')
    parser.add_argument('-rr', '--rankedResampling', action='store_true', dest='rankedResampling', help="If turned on, IFGs resampled based off of the average of pixels in a given resampling window corresponding to the connected component mode (if multilooking specified). Program defaults to lanczos resampling algorithm through gdal (if multilooking specified).")
    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap', type=float, default=0.0081, help='Minimum km\u00b2 area of overlap of scenes wrt specified bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single pixel at standard 90m resolution"')
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

def extractMetaDict(aria_prod,metadata):
    os.environ['GDAL_PAM_ENABLED'] = 'NO'
    meta = {}
    for i in aria_prod.products[1]:
        metaName = i[metadata][0]
        ds = gdal.Open(metaName)
        # return [min, max, mean, std]
        stat  = ds.GetRasterBand(1).GetStatistics(True, True)
        meta[i['pair_name'][0]] = stat[2]
        ds = None

    return meta

def extractUTCtime(aria_prod):
    f = aria_prod.products[1]
    utcDict = {}

    for i in range(0,len(f)):
        #Grab pair name of the product
        pairName = aria_prod.products[1][i]['pair_name']

        #Grab mid-times, append to a list and find minimum and maximum mid-times
        midTimeList = aria_prod.products[0][i]['azimuthZeroDopplerMidTime']
        midTime=[]
        for j in midTimeList:
            midTime.append(datetime.strptime(j,'%Y-%m-%dT%H:%M:%S.%f'))
        minMidTime = min(midTime)
        maxMidTime = max(midTime)

        #Calculate time difference between minimum start and maximum end time, and add it to mean start time
        #Write calculated UTC time into a dictionary with associated pair names as keys
        timeDelta = (maxMidTime - minMidTime)/2
        utcTime = (minMidTime+timeDelta).time()
        utcDict[pairName[0]] = utcTime.strftime("%H:%M:%S.%f")
    return utcDict

def generateStack(aria_prod,inputFiles,outputFileName,workdir='./', refDlist=None):

    UTC_time = extractUTCtime(aria_prod)
    bPerp = extractMetaDict(aria_prod,'bPerpendicular')
    incAng = extractMetaDict(aria_prod,'incidenceAngle')
    lookAng = extractMetaDict(aria_prod,'lookAngle')
    azimuthAng = extractMetaDict(aria_prod,'azimuthAngle')
    os.environ['GDAL_PAM_ENABLED'] = 'YES'

    ##Progress bar
    from ARIAtools import progBar
    prog_bar = progBar.progressBar(maxValue=len(aria_prod.products[1]), print_msg='Creating stack:')

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
        Dlist = sorted(intList)
    elif inputFiles in ['tropocorrected_products', 'tropocorrected', 'tropo']:
        domainName = 'unwrappedPhase'
        intList = glob.glob(os.path.join(workdir,'tropocorrected_products','[0-9]*[0-9].vrt'))
        dataType = "Float32"
        print('Number of tropocorrected unwrapped interferograms discovered: ', len(intList))
        Dlist = sorted(intList)
    elif inputFiles in ['coherence', 'Coherence', 'coh']:
        domainName = 'Coherence'
        cohList = glob.glob(os.path.join(workdir,'coherence','[0-9]*[0-9].vrt'))
        dataType = "Float32"
        print('Number of coherence discovered: ', len(cohList))
        Dlist = sorted(cohList)
    elif inputFiles in ['connectedComponents','connectedComponent','connComp']:
        domainName = 'connectedComponents'
        connCompList = glob.glob(os.path.join(workdir,'connectedComponents','[0-9]*[0-9].vrt'))
        dataType = "Int16"
        print('Number of connectedComponents discovered: ', len(connCompList))
        Dlist = sorted(connCompList)
    else:
        print('Stacks can be created for unwrapped interferogram, coherence and connectedComponent VRT files')

    # Confirm 1-to-1 match between UNW and other derived products
    newDlist = [os.path.basename(i).split('.vrt')[0] for i in Dlist]
    if refDlist and newDlist != refDlist:
        raise Exception('Discrepancy between {} number of UNW products and {} number of {} products'.format(
            len(refDlist), len(newDlist), domainName))

    for data in Dlist:
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
    wavelength = aria_prod.products[0][0]['wavelength'][0]
    startRange = aria_prod.products[0][0]['slantRangeStart'][0]
    endRange = aria_prod.products[0][0]['slantRangeEnd'][0]
    rangeSpacing = aria_prod.products[0][0]['slantRangeSpacing'][0]
    orbitDirection = str.split(os.path.basename(aria_prod.files[0]),'-')[2]

    stack_dir = os.path.join(workdir, 'stack')
    with open( os.path.join(stack_dir, outputFileName+'.vrt'), 'w') as fid:
        fid.write( '''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
        <SRS>{proj}</SRS>
        <GeoTransform>{GT0},{GT1},{GT2},{GT3},{GT4},{GT5}</GeoTransform>\n'''.format(
            xsize=xsize, ysize=ysize,
            proj=projection,
            GT0=gt[0],GT1=gt[1],GT2=gt[2],GT3=gt[3],GT4=gt[4],GT5=gt[5]))

        for data in enumerate(Dlist):
            metadata = {}
            dates = data[1].split('/')[-1][:-4]
            width = None
            height = None
            path = None
            ##Update progress bar
            prog_bar.update(data[0]+1, suffix=dates)

            ds = gdal.Open(data[1], gdal.GA_ReadOnly)
            width = ds.RasterXSize
            height = ds.RasterYSize
            ds = None

            metadata['wavelength'] = wavelength
            metadata['utcTime'] = UTC_time[dates]
            metadata['bPerp'] = bPerp[dates]
            metadata['incAng'] = incAng[dates]
            metadata['lookAng'] = lookAng[dates]
            metadata['azimuthAng'] = azimuthAng[dates]
            if orbitDirection == 'D':
                metadata['orbit_direction'] = 'DESCENDING'
            elif orbitDirection == 'A':
                metadata['orbit_direction'] = 'ASCENDING'
            else:
                print('Orbit direction not recognized')
                metadata['orbit_direction'] = 'UNKNOWN'

            path = os.path.relpath(os.path.abspath(data[1]), start=stack_dir)
            outstr = '''    <VRTRasterBand dataType="{dataType}" band="{index}">
        <SimpleSource>
            <SourceFilename relativeToVRT="1">{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}" DataType="{dataType}"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}" ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain='{domainName}'>
            <MDI key="Dates">{dates}</MDI>
            <MDI key="Wavelength (m)">{wvl}</MDI>
            <MDI key="UTCTime (HH:MM:SS.ss)">{acq}</MDI>
            <MDI key="perpendicularBaseline">{bPerp}</MDI>
            <MDI key="incidenceAngle">{incAng}</MDI>
            <MDI key="lookAngle">{lookAng}</MDI>
            <MDI key="azimuthAngle">{azimuthAng}</MDI>
            <MDI key="startRange">{start_range}</MDI>
            <MDI key="endRange">{end_range}</MDI>
            <MDI key="slantRangeSpacing">{range_spacing}</MDI>
            <MDI key="orbitDirection">{orbDir}</MDI>
        </Metadata>
    </VRTRasterBand>\n'''.format(domainName=domainName,width=width, height=height,
                                xmin=xmin, ymin=ymin,
                                xsize=xsize, ysize=ysize,
                                dates=dates, acq=metadata['utcTime'],
                                wvl=metadata['wavelength'], index=data[0]+1,
                                path=path, dataType=dataType, bPerp=metadata['bPerp'],
                                incAng=metadata['incAng'],lookAng=metadata['lookAng'],azimuthAng=metadata['azimuthAng'],
                                start_range=startRange, end_range=endRange, range_spacing=rangeSpacing, orbDir=metadata['orbit_direction'])
            fid.write(outstr)
        fid.write('</VRTDataset>')
        prog_bar.close()
        print(outputFileName, ': stack generated')

    return newDlist


def main(inps=None):
    inps = cmdLineParse()

    print("***Time-series Preparation Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    # pass number of threads for gdal multiprocessing computation
    if inps.num_threads.lower()=='all':
        import multiprocessing
        print('User specified use of all %s threads for gdal multiprocessing'%(str(multiprocessing.cpu_count())))
        inps.num_threads='ALL_CPUS'
    print('Thread count specified for gdal multiprocessing = %s'%(inps.num_threads))

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[0], standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr, arrshape, proj = merged_productbbox(standardproduct_info.products[0], standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion, num_threads=inps.num_threads, minimumOverlap=inps.minimumOverlap, verbose=inps.verbose)


    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        print('Download/cropping DEM')
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, \
            prods_TOTbbox_metadatalyr, proj, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Load or download mask (if specified).
    if inps.mask is not None:
        inps.mask = prep_mask([[item for sublist in [list(set(d['amplitude'])) for d in standardproduct_info.products[1] if 'amplitude' in d] for item in sublist], [item for sublist in [list(set(d['pair_name'])) for d in standardproduct_info.products[1] if 'pair_name' in d] for item in sublist]],inps.mask, standardproduct_info.bbox_file, prods_TOTbbox, proj, amp_thresh=inps.amp_thresh, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Extract
    layers=['unwrappedPhase','coherence']
    print('\nExtracting unwrapped phase, coherence, and connected components for each interferogram pair')
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat, stitchMethodType='overlap', verbose=inps.verbose, num_threads=inps.num_threads, multilooking=inps.multilooking)

    layers=['incidenceAngle','lookAngle','azimuthAngle']
    print('\nExtracting single incidence angle, look angle and azimuth angle files valid over common interferometric grid')
    export_products([dict(zip([k for k in set(k for d in standardproduct_info.products[1] for k in d)], [[item for sublist in [list(set(d[k])) for d in standardproduct_info.products[1] if k in d] for item in sublist] for k in set(k for d in standardproduct_info.products[1] for k in d)]))], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat, stitchMethodType='overlap', verbose=inps.verbose, num_threads=inps.num_threads, multilooking=inps.multilooking)

    if inps.bperp==True:
        layers=['bPerpendicular']
        print('\nExtracting perpendicular baseline grids for each interferogram pair')
        export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat, stitchMethodType='overlap', verbose=inps.verbose, num_threads=inps.num_threads, multilooking=inps.multilooking)

    # Extracting other layers, if specified
    if inps.layers:
        if inps.layers.lower()=='all':
            inps.layers=list(standardproduct_info.products[1][0].keys())
            # Must also remove productBoundingBoxes & pair-names because they are not raster layers
            layers=[i for i in inps.layers if i not in ['unwrappedPhase','coherence','incidenceAngle','lookAngle','azimuthAngle','bPerpendicular','productBoundingBox','productBoundingBoxFrames','pair_name']]
        else:
            inps.layers=list(inps.layers.split(','))
            layers=[i for i in inps.layers if i not in ['unwrappedPhase','coherence','incidenceAngle','lookAngle','azimuthAngle','bPerpendicular']]
            layers=[i.replace(' ','') for i in layers]

        if layers!=[]:
            print('\nExtracting optional, user-specified layers %s for each interferogram pair'%(layers))
            export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat, stitchMethodType='overlap', verbose=inps.verbose, num_threads=inps.num_threads, multilooking=inps.multilooking)

    # If necessary, resample DEM/mask AFTER they have been used to extract metadata layers and mask output layers, respectively
    if inps.multilooking is not None:
        # Import functions
        from ARIAtools.shapefile_util import open_shapefile
        from ARIAtools.vrtmanager import resampleRaster
        bounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds
        # Resample mask
        if inps.mask is not None:
            resampleRaster(inps.mask.GetDescription(), inps.multilooking, bounds, prods_TOTbbox, inps.rankedResampling, outputFormat=inps.outputFormat, num_threads=inps.num_threads)
        # Resample DEM
        if demfile is not None:
            resampleRaster(demfile.GetDescription(), inps.multilooking, bounds, prods_TOTbbox, inps.rankedResampling, outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Perform GACOS-based tropospheric corrections (if specified).
    if inps.tropo_products:
        tropo_correction(standardproduct_info.products, inps.tropo_products, standardproduct_info.bbox_file, prods_TOTbbox, outDir=inps.workdir, outputFormat=inps.outputFormat, verbose=inps.verbose, num_threads=inps.num_threads)

    # Generate Stack
    refDlist = generateStack(standardproduct_info,'coherence','cohStack',workdir=inps.workdir)
    refDlist = generateStack(standardproduct_info,'connectedComponents','connCompStack',workdir=inps.workdir, refDlist=refDlist)
    if inps.tropo_products:
        refDlist = generateStack(standardproduct_info,'tropocorrected_products','unwrapStack',workdir=inps.workdir, refDlist=refDlist)
    else:
        refDlist = generateStack(standardproduct_info,'unwrappedPhase','unwrapStack',workdir=inps.workdir, refDlist=refDlist)
