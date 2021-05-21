#!/usr/bin/env python3
"""Extract unwrapped interferogram, coherence, ⊥ baseline, and LOS file(s).
This script is intended to extract required information and files to carry-out
time series analysis on ARIA products.

Copyright 2019, by the California Institute of Technology.
ALL RIGHTS RESERVED. United States Government Sponsorship acknowledged.

Author(s): Simran Sangha, David Bekaert, & Emre Havazli
"""

import os
import glob
import logging
import argparse
import multiprocessing
from datetime import datetime
from osgeo import gdal

# Import functions
from ARIAtools import progBar
from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.mask_util import prep_mask
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.vrtmanager import resampleRaster
from ARIAtools.extractProduct import (merged_productbbox, prep_dem,
                                      export_products, tropo_correction)

gdal.UseExceptions()
# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

log = logging.getLogger(__name__)


def create_parser():
    """Parser to read command line arguments."""
    parser = argparse.ArgumentParser(description='Prepare ARIA products for '
                                     'time series processing.')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
                        required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./',
                        help='Specify directory to deposit all outputs. '
                        'Default is local directory where script is launched.')
    parser.add_argument('-tp', '--tropo_products', dest='tropo_products',
                        type=str, default=None, help='Path to director(ies) '
                        'or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None,
                        help='Specify layers to extract as a comma '
                        'deliminated list bounded by single quotes. '
                        'Allowed keys are: "unwrappedPhase", "coherence", '
                        '"amplitude", "bPerpendicular", "bParallel", '
                        '"incidenceAngle", "lookAngle", "azimuthAngle", '
                        '"ionosphere". If "all" is specified, then all layers'
                        'are extracted.'
                        'If blank, will only extract bounding box.')
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
    # parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',
    #                      type=str, default='overlap', help='Method applied '
    #                      'to stitch the unwrapped data. Either "overlap", '
    #                      'where product overlap is minimized, or "2stage", '
    #                      'where minimization is done on connected '
    #                      'components, are allowed methods. '
    #                      'default: "overlap".')
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
    parser.add_argument('-bp', '--bperp', action='store_true', dest='bperp',
                        help='If turned on, extracts perpendicular baseline '
                        'grids. Default: A single perpendicular baseline '
                        'value is calculated and included in the metadata of '
                        'stack cubes for each pair.')
    parser.add_argument('-ml', '--multilooking', dest='multilooking', type=int,
                        default=None, help='Multilooking factor is an integer '
                        'multiple of standard resolution. '
                        'E.g. 2 = 90m*2 = 180m')
    parser.add_argument('-rr', '--rankedResampling', action='store_true',
                        dest='rankedResampling', help='If turned on, IFGs '
                        'resampled based off of the average of pixels in a '
                        'given resampling window corresponding to the '
                        'connected component mode (if multilooking specified).'
                        ' Program defaults to lanczos resampling algorithm'
                        ' through gdal (if multilooking specified).')
    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap',
                        type=float, default=0.0081, help='Minimum km\u00b2 '
                        'area of overlap of scenes wrt specified bounding box.'
                        ' Default 0.0081 = 0.0081km\u00b2 = area of single'
                        ' pixel at standard 90m resolution')
    parser.add_argument('-verbose', '--verbose', action='store_true',
                        dest='verbose', help="Toggle verbose mode on.")

    return parser


def cmd_line_parse(iargs=None):
    """Command line parser."""
    parser = create_parser()
    return parser.parse_args(args=iargs)


def extract_meta_dict(aria_prod, metadata):
    """Extract metadata from products."""
    os.environ['GDAL_PAM_ENABLED'] = 'NO'
    meta = {}
    for i in aria_prod.products[1]:
        meta_name = i[metadata][0]
        data_set = gdal.Open(meta_name)
        # return [min, max, mean, std]
        stat = data_set.GetRasterBand(1).GetStatistics(True, True)
        meta[i['pair_name'][0]] = stat[2]
        data_set = None

    return meta


def extract_utc_time(aria_prod):
    """Extract UTC time from products."""
    f = aria_prod.products[1]
    utc_dict = {}

    for i in range(0, len(f)):
        # Grab pair name of the product
        pair_name = aria_prod.products[1][i]['pair_name']

        # Grab mid-times, append to a list and find minimum and
        # maximum mid-times
        mid_time_list = aria_prod.products[0][i]['azimuthZeroDopplerMidTime']
        mid_time = []
        for j in mid_time_list:
            mid_time.append(datetime.strptime(j, '%Y-%m-%dT%H:%M:%S.%f'))
        min_mid_time = min(mid_time)
        max_mid_time = max(mid_time)

        # Calculate time difference between minimum start and maximum end time,
        # and add it to mean start time.
        # Write calculated UTC time into a dictionary
        # with associated pair names as keys
        time_delta = (max_mid_time - min_mid_time)/2
        utc_time = (min_mid_time+time_delta).time()
        utc_dict[pair_name[0]] = utc_time.strftime("%H:%M:%S.%f")
    return utc_dict


def generate_stack(aria_prod, input_files, output_file_name,
                   workdir='./', ref_dlist=None):
    """Generate time series stack."""
    utc_time = extract_utc_time(aria_prod)
    b_perp = extract_meta_dict(aria_prod, 'bPerpendicular')
    inc_angle = extract_meta_dict(aria_prod, 'incidenceAngle')
    look_ang = extract_meta_dict(aria_prod, 'lookAngle')
    azimuth_ang = extract_meta_dict(aria_prod, 'azimuthAngle')
    os.environ['GDAL_PAM_ENABLED'] = 'YES'

    # Progress bar
    prog_bar = progBar.progressBar(maxValue=len(aria_prod.products[1]),
                                   print_msg='Creating stack:')

    # Set up single stack file
    if not os.path.exists(os.path.join(workdir, 'stack')):
        print('Creating directory: {0}'.format(os.path.join(workdir, 'stack')))
        os.makedirs(os.path.join(workdir, 'stack'))
    else:
        print('Directory {0} already exists.'.format(os.path.join(
              workdir, 'stack')))

    if input_files in ['unwrappedPhase', 'unwrapped', 'unw']:
        domain_name = 'unwrappedPhase'
        int_list = glob.glob(os.path.join(workdir, 'unwrappedPhase',
                             '[0-9]*[0-9].vrt'))
        data_type = "Float32"
        print('Number of unwrapped interferograms discovered: ', len(int_list))
        dlist = sorted(int_list)
    elif input_files in ['tropocorrected_products', 'tropocorrected', 'tropo']:
        domain_name = 'unwrappedPhase'
        int_list = glob.glob(os.path.join(workdir, 'tropocorrected_products',
                             '[0-9]*[0-9].vrt'))
        data_type = "Float32"
        print('Number of tropocorrected unwrapped interferograms discovered: ',
              len(int_list))
        dlist = sorted(int_list)
    elif input_files in ['coherence', 'Coherence', 'coh']:
        domain_name = 'Coherence'
        coh_list = glob.glob(os.path.join(workdir, 'coherence',
                             '[0-9]*[0-9].vrt'))
        data_type = "Float32"
        print('Number of coherence files discovered: ', len(coh_list))
        dlist = sorted(coh_list)
    elif input_files in ['connectedComponents', 'connectedComponent',
                         'connComp']:
        domain_name = 'connectedComponents'
        conn_comp_list = glob.glob(os.path.join(workdir,
                                   'connectedComponents', '[0-9]*[0-9].vrt'))
        data_type = "Int16"
        print('Number of connectedComponents discovered: ',
              len(conn_comp_list))
        dlist = sorted(conn_comp_list)
    elif input_files in ['amplitude', 'Amplitude', 'amp']:
        domain_name = 'Amplitude'
        amp_list = glob.glob(os.path.join(workdir, 'amplitude',
                             '[0-9]*[0-9].vrt'))
        data_type = "Float32"
        print('Number of amplitude files discovered: ', len(amp_list))
        dlist = sorted(amp_list)
    else:
        print('Stacks can be created for unwrapped interferogram, '
              'coherence and connectedComponent VRT files')

    # Confirm 1-to-1 match between UNW and other derived products
    new_dlist = [os.path.basename(i).split('.vrt')[0] for i in dlist]
    if ref_dlist and new_dlist != ref_dlist:
        tropo_dlist = glob.glob(os.path.join(workdir,
                                'tropocorrected_products', '[0-9]*[0-9].vrt'))
        tropo_dlist = sorted([os.path.basename(i).split(
                            '.vrt')[0] for i in tropo_dlist])
        tropo_path = os.path.join(workdir, 'tropocorrected_products')
        if os.path.exists(tropo_path) and tropo_dlist == ref_dlist:
            log.warning('Discrepancy between "tropocorrected_products" '
                        'products (%s files) and %s products (%s files),'
                        'rejecting scenes not common between both',
                        len(ref_dlist), domain_name, len(new_dlist))
        else:
            log.warning('Discrepancy between "unwrappedPhase" products '
                        '(%s files) and %s products (%s files), '
                        'rejecting scenes not common between both',
                        len(ref_dlist), domain_name, len(new_dlist))
        # subset to match other datasets
        subset_ind = [i[0] for i in enumerate(new_dlist) if i[1] in ref_dlist]
        new_dlist = [new_dlist[i] for i in subset_ind]
        dlist = [dlist[i] for i in subset_ind]

    for data in dlist:
        width = None
        height = None

        data_set = gdal.Open(data, gdal.GA_ReadOnly)
        width = data_set.RasterXSize
        height = data_set.RasterYSize
        geo_trans = data_set.GetGeoTransform()
        projection = data_set.GetProjection()
        data_set = None

    # setting up a subset of the stack
    ymin, ymax, xmin, xmax = [0, height, 0, width]

    xsize = xmax - xmin
    ysize = ymax - ymin

    # extraction of radar meta-data
    wavelength = aria_prod.products[0][0]['wavelength'][0]
    start_range = aria_prod.products[0][0]['slantRangeStart'][0]
    end_range = aria_prod.products[0][0]['slantRangeEnd'][0]
    range_spacing = aria_prod.products[0][0]['slantRangeSpacing'][0]
    orbit_direction = str.split(os.path.basename(aria_prod.files[0]), '-')[2]

    stack_dir = os.path.join(workdir, 'stack')
    with open(os.path.join(stack_dir, output_file_name+'.vrt'), 'w') as fid:
        fid.write('''<VRTDataset rasterXSize="{xsize}" rasterYSize="{ysize}">
        <SRS>{proj}</SRS>
        <GeoTransform>{GT0},{GT1},{GT2},{GT3},{GT4},{GT5}</GeoTransform>\n
        '''.format(xsize=xsize, ysize=ysize, proj=projection, GT0=geo_trans[0],
                   GT1=geo_trans[1], GT2=geo_trans[2], GT3=geo_trans[3],
                   GT4=geo_trans[4], GT5=geo_trans[5]))

        for data in enumerate(dlist):
            metadata = {}
            dates = data[1].split('/')[-1][:-4]
            width = None
            height = None
            path = None
            # Update progress bar
            prog_bar.update(data[0]+1, suffix=dates)

            data_set = gdal.Open(data[1], gdal.GA_ReadOnly)
            width = data_set.RasterXSize
            height = data_set.RasterYSize
            data_set = None

            metadata['wavelength'] = wavelength
            metadata['utcTime'] = utc_time[dates]
            metadata['bPerp'] = b_perp[dates]
            metadata['incAng'] = inc_angle[dates]
            metadata['lookAng'] = look_ang[dates]
            metadata['azimuthAng'] = azimuth_ang[dates]
            if orbit_direction == 'D':
                metadata['orbit_direction'] = 'DESCENDING'
            elif orbit_direction == 'A':
                metadata['orbit_direction'] = 'ASCENDING'
            else:
                print('Orbit direction not recognized')
                metadata['orbit_direction'] = 'UNKNOWN'

            path = os.path.abspath(data[1])
            outstr = '''  <VRTRasterBand dataType="{data_type}" band="{index}">
        <SimpleSource>
            <SourceFilename relativeToVRT="1">{path}</SourceFilename>
            <SourceBand>1</SourceBand>
            <SourceProperties RasterXSize="{width}" RasterYSize="{height}"
                DataType="{data_type}"/>
            <SrcRect xOff="{xmin}" yOff="{ymin}" xSize="{xsize}"
                ySize="{ysize}"/>
            <DstRect xOff="0" yOff="0" xSize="{xsize}" ySize="{ysize}"/>
        </SimpleSource>
        <Metadata domain='{domain_name}'>
            <MDI key="Dates">{dates}</MDI>
            <MDI key="Wavelength (m)">{wvl}</MDI>
            <MDI key="UTCTime (HH:MM:SS.ss)">{acq}</MDI>
            <MDI key="perpendicularBaseline">{b_perp}</MDI>
            <MDI key="incidenceAngle">{inc_angle}</MDI>
            <MDI key="lookAngle">{look_ang}</MDI>
            <MDI key="azimuthAngle">{azimuth_ang}</MDI>
            <MDI key="startRange">{start_range}</MDI>
            <MDI key="endRange">{end_range}</MDI>
            <MDI key="slantRangeSpacing">{range_spacing}</MDI>
            <MDI key="orbitDirection">{orbDir}</MDI>
        </Metadata>
    </VRTRasterBand>\n'''.format(domain_name=domain_name, width=width,
                                 height=height, xmin=xmin, ymin=ymin,
                                 xsize=xsize, ysize=ysize, dates=dates,
                                 acq=metadata['utcTime'],
                                 wvl=metadata['wavelength'], index=data[0]+1,
                                 path=path, data_type=data_type,
                                 b_perp=metadata['bPerp'],
                                 inc_angle=metadata['incAng'],
                                 look_ang=metadata['lookAng'],
                                 azimuth_ang=metadata['azimuthAng'],
                                 start_range=start_range, end_range=end_range,
                                 range_spacing=range_spacing,
                                 orbDir=metadata['orbit_direction'])
            fid.write(outstr)
        fid.write('</VRTDataset>')
        prog_bar.close()
        print(output_file_name, ': stack generated')

    return new_dlist


def main(inps=None):
    """Run time series prepation."""
    inps = cmd_line_parse()
    print ('*****************************************************************')
    print ('*** Time-series Preparation Function ***')
    print ('*****************************************************************')
    # if user bbox was specified, file(s) not meeting imposed spatial
    # criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped
    # “radarmetadata info” and “data layer keys+paths” dictionaries for each
    # standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file']
    # (if bbox specified)
    standardproduct_info = ARIA_standardproduct(
        inps.imgfile, bbox=inps.bbox,
        workdir=inps.workdir, verbose=inps.verbose)

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
    (standardproduct_info.products[0], standardproduct_info.products[1],
     standardproduct_info.bbox_file, prods_tot_bbox,
     prods_tot_bbox_metadatalyr, arrshape,
     proj) = merged_productbbox(standardproduct_info.products[0],
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
            prods_tot_bbox, prods_tot_bbox_metadatalyr, proj,
            arrshape=arrshape, workdir=inps.workdir,
            outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Load or download mask (if specified).
    if inps.mask is not None:
        inps.mask = prep_mask([[item for sublist in [list(set(d['amplitude']))
                              for d in standardproduct_info.products[1]
                              if 'amplitude' in d]
                              for item in sublist], [item for sublist
                              in [list(set(d['pair_name'])) for d
                                  in standardproduct_info.products[1]
                                  if 'pair_name' in d]
                              for item in sublist]], inps.mask,
                              standardproduct_info.bbox_file,
                              prods_tot_bbox, proj, amp_thresh=inps.amp_thresh,
                              arrshape=arrshape,
                              workdir=inps.workdir,
                              outputFormat=inps.outputFormat,
                              num_threads=inps.num_threads)

    # Extract
    layers = ['unwrappedPhase', 'coherence']
    print('\nExtracting unwrapped phase, coherence, '
          'and connected components for each interferogram pair')
    export_products(standardproduct_info.products[1],
                    standardproduct_info.bbox_file, prods_tot_bbox, layers,
                    dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask,
                    outDir=inps.workdir, outputFormat=inps.outputFormat,
                    stitchMethodType='overlap', verbose=inps.verbose,
                    num_threads=inps.num_threads,
                    multilooking=inps.multilooking)

    layers = ['incidenceAngle', 'lookAngle', 'azimuthAngle']
    print('\nExtracting single incidence angle, look angle and azimuth angle '
          'files valid over common interferometric grid')
    export_products([dict(zip([k for k in set(k for d in
                    standardproduct_info.products[1] for k in d)],
                    [[item for sublist in [list(set(d[k]))
                      for d in standardproduct_info.products[1]
                      if k in d] for item in sublist]
                        for k in set(k for d in
                                     standardproduct_info.products[1]
                                     for k in d)]))],
                    standardproduct_info.bbox_file, prods_tot_bbox, layers,
                    dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask,
                    outDir=inps.workdir, outputFormat=inps.outputFormat,
                    stitchMethodType='overlap', verbose=inps.verbose,
                    num_threads=inps.num_threads,
                    multilooking=inps.multilooking)
    if inps.bperp:
        layers = ['bPerpendicular']
        print('\nExtracting perpendicular baseline grids for each '
              'interferogram pair')
        export_products(standardproduct_info.products[1],
                        standardproduct_info.bbox_file, prods_tot_bbox, layers,
                        dem=demfile, lat=Latitude, lon=Longitude,
                        mask=inps.mask, outDir=inps.workdir,
                        outputFormat=inps.outputFormat,
                        stitchMethodType='overlap', verbose=inps.verbose,
                        num_threads=inps.num_threads,
                        multilooking=inps.multilooking)

    # Extracting other layers, if specified
    if inps.layers:
        if inps.layers.lower() == 'all':
            inps.layers = list(standardproduct_info.products[1][0].keys())
            # Must also remove productBoundingBoxes &
            # pair-names because they are not raster layers
            layers_list = ['unwrappedPhase', 'coherence', 'incidenceAngle',
                           'lookAngle', 'azimuthAngle', 'bPerpendicular',
                           'productBoundingBox', 'productBoundingBoxFrames',
                           'pair_name']
            layers = [i for i in inps.layers if i not in layers_list]
        else:
            inps.layers = list(inps.layers.split(','))
            layers_list = ['unwrappedPhase', 'coherence', 'incidenceAngle',
                           'lookAngle', 'azimuthAngle', 'bPerpendicular']
            layers = [i for i in inps.layers if i not in layers_list]
            layers = [i.replace(' ', '') for i in layers]

        if layers != []:
            print('\nExtracting optional, user-specified layers %s for each '
                  'interferogram pair' % (layers))
            export_products(standardproduct_info.products[1],
                            standardproduct_info.bbox_file, prods_tot_bbox,
                            layers, dem=demfile, lat=Latitude, lon=Longitude,
                            mask=inps.mask, outDir=inps.workdir,
                            outputFormat=inps.outputFormat,
                            stitchMethodType='overlap', verbose=inps.verbose,
                            num_threads=inps.num_threads,
                            multilooking=inps.multilooking)
    else:
        inps.layers = []

    # If necessary, resample DEM/mask AFTER they have been used to extract
    # metadata layers and mask output layers, respectively
    if inps.multilooking is not None:
        # Import functions
        bounds = open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds
        # Resample mask
        if inps.mask is not None:
            resampleRaster(inps.mask.GetDescription(), inps.multilooking,
                           bounds, prods_tot_bbox, inps.rankedResampling,
                           outputFormat=inps.outputFormat,
                           num_threads=inps.num_threads)
        # Resample DEM
        if demfile is not None:
            resampleRaster(demfile.GetDescription(), inps.multilooking,
                           bounds, prods_tot_bbox, inps.rankedResampling,
                           outputFormat=inps.outputFormat,
                           num_threads=inps.num_threads)

    # Perform GACOS-based tropospheric corrections (if specified).
    if inps.tropo_products:
        tropo_correction(standardproduct_info.products, inps.tropo_products,
                         standardproduct_info.bbox_file, prods_tot_bbox,
                         outDir=inps.workdir, outputFormat=inps.outputFormat,
                         verbose=inps.verbose, num_threads=inps.num_threads)

    # Generate Stack
    if inps.tropo_products:
        ref_dlist = generate_stack(standardproduct_info,
                                   'tropocorrected_products', 'unwrapStack',
                                   workdir=inps.workdir)
    else:
        ref_dlist = generate_stack(standardproduct_info, 'unwrappedPhase',
                                   'unwrapStack', workdir=inps.workdir)
    ref_dlist = generate_stack(standardproduct_info, 'coherence', 'cohStack',
                               workdir=inps.workdir, ref_dlist=ref_dlist)
    ref_dlist = generate_stack(standardproduct_info, 'connectedComponents',
                               'connCompStack', workdir=inps.workdir,
                               ref_dlist=ref_dlist)
    # If amplitude files extracted
    if 'amplitude' in inps.layers:
        ref_dlist = generate_stack(standardproduct_info, 'amplitude',
                                   'ampStack', workdir=inps.workdir,
                                   ref_dlist=ref_dlist)
