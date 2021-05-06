#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import numpy as np
import glob
from osgeo import gdal, osr
import logging
import requests
from ARIAtools.logger import logger

from ARIAtools.shapefile_util import open_shapefile, chunk_area
from ARIAtools.mask_util import prep_mask
from ARIAtools.unwrapStitching import product_stitch_overlap, product_stitch_2stage
from ARIAtools import version

gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

log = logging.getLogger(__name__)

_world_dem = "https://portal.opentopography.org/API/globaldem?demtype=SRTMGL1_E&west={}&south={}&east={}&north={}&outputFormat=GTiff"

def createParser():
    '''
       Extract specified product layers. The default will export all layers.
    '''
    import argparse
    parser = argparse.ArgumentParser(description='Program to extract data and meta-data layers from ARIA standard GUNW products. Program will handle cropping/stitching when needed. By default, the program will crop all IFGs to bounds determined by the common intersection and bbox (if specified)')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-tp', '--tropo_products', dest='tropo_products', type=str, default=None, help='Path to director(ies) or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None, help='Specify layers to extract as a comma deliminated list bounded by single quotes. Allowed keys are: "unwrappedPhase", "coherence", "amplitude", "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", "azimuthAngle", "ionosphere". If "all" is specified, then all layers are extracted. If blank, will only extract bounding box.')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
            default=None, help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str,
            help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Path to mask file or 'Download'. File needs to be GDAL compatabile, contain spatial reference information, and have invalid/valid data represented by 0/1, respectively. If 'Download', will use GSHHS water mask. If 'NLCD', will mask classes 11, 12, 90, 95; see: https://www.mrlc.gov/national-land-cover-database-nlcd-201://www.mrlc.gov/national-land-cover-database-nlcd-2016")
    parser.add_argument('-at', '--amp_thresh', dest='amp_thresh', default=None, type=str, help='Amplitude threshold below which to mask. Specify "None" to not use amplitude mask. By default "None".')
    parser.add_argument('-nt', '--num_threads', dest='num_threads', default='2', type=str, help='Specify number of threads for multiprocessing operation in gdal. By default "2". Can also specify "All" to use all available threads.')
#    parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',  type=str, default='overlap', help="Method applied to stitch the unwrapped data. Either 'overlap', where product overlap is minimized, or '2stage', where minimization is done on connected components, are allowed methods. Default is 'overlap'.")
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str, default='VRT', help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By default files are generated virtually except for "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle", "unwrappedPhase" as these are require either DEM intersection or corrections to be applied')
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-ml', '--multilooking', dest='multilooking', type=int, default=None, help='Multilooking factor is an integer multiple of standard resolution. E.g. 2 = 90m*2 = 180m')
    parser.add_argument('-rr', '--rankedResampling', action='store_true', dest='rankedResampling', help="If turned on, IFGs resampled based off of the average of pixels in a given resampling window corresponding to the connected component mode (if multilooking specified). Program defaults to lanczos resampling algorithm through gdal (if multilooking specified).")
    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap', type=float, default=0.0081, help='Minimum km\u00b2 area of overlap of scenes wrt specified bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single pixel at standard 90m resolution"')
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")

    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)

class InterpCube(object):
    '''
        Class to interpolate intersection of cube with DEM
    '''

    def __init__(self, inobj, hgtobj, latobj, lonobj):
        '''
            Init with h5py dataset.
        '''
        self.data = inobj[:]
        self.hgts = hgtobj[:]
        self.offset = None
        self.interp = []
        self.latobj = latobj[:]
        self.lonobj = lonobj[:]

        self.createInterp()

    def createInterp(self):
        '''
            Create interpolators.
        '''
        from scipy.interpolate import RectBivariateSpline
        self.offset = np.mean(self.data)
        for i in range(len(self.hgts)):
            self.interp.append( RectBivariateSpline(self.latobj, self.lonobj, self.data[i]-self.offset))

    def __call__(self, line, pix, h):
        '''
            Interpolate at a single point.
        '''
        from scipy.interpolate import interp1d

        vals  = np.array( [x(line,pix)[0,0] for x in self.interp])
        est = interp1d(self.hgts, vals, kind='cubic')
        return est(h) + self.offset

class metadata_qualitycheck:
    '''
        Metadata quality control function. Artifacts recognized based off of covariance of cross-profiles.
        Bug-fix varies based off of layer of interest.
        Verbose mode generates a series of quality control plots with these profiles.
    '''

    def __init__(self, data_array, prod_key, outname, verbose=None):
        # Pass inputs
        self.data_array = data_array
        self.prod_key = prod_key
        self.outname = outname
        self.verbose = verbose
        self.data_array_band=data_array.GetRasterBand(1).ReadAsArray()
        #mask by nodata value
        self.data_array_band=np.ma.masked_where(self.data_array_band == self.data_array.GetRasterBand(1).GetNoDataValue(), self.data_array_band)

        if self.verbose: logger.setLevel(logging.DEBUG)

        # Run class
        self.__run__()

    def __truncateArray__(self, data_array_band, Xmask, Ymask):
        # Mask columns/rows which are entirely made up of 0s
        #first must crop all columns with no valid values
        nancols=np.all(data_array_band.mask == True, axis=0)
        data_array_band=data_array_band[:,~nancols]
        Xmask=Xmask[:,~nancols]
        Ymask=Ymask[:,~nancols]
        #first must crop all rows with no valid values
        nanrows=np.all(data_array_band.mask == True, axis=1)
        data_array_band=data_array_band[~nanrows]
        Xmask=Xmask[~nanrows]
        Ymask=Ymask[~nanrows]

        return data_array_band, Xmask, Ymask

    def __getCovar__(self, prof_direc, profprefix=''):
        from scipy.stats import linregress
        # Mask columns/rows which are entirely made up of 0s
        if self.data_array_band.mask.size!=1 and True in self.data_array_band.mask:
            Xmask,Ymask = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1), np.arange(0, self.data_array_band.shape[0], 1))
            self.data_array_band, Xmask, Ymask = self.__truncateArray__(self.data_array_band, Xmask, Ymask)

        #append prefix for plot names
        prof_direc=profprefix+prof_direc

        #iterate through transpose of matrix if looking in azimuth
        arrT=''
        if 'azimuth' in prof_direc:
            arrT='.T'
        # Cycle between range and azimuth profiles
        rsquaredarr=[]
        std_errarr=[]
        for i in enumerate(eval('self.data_array_band%s'%(arrT))):
            mid_line=i[1]
            xarr=np.array(range(len(mid_line)))
            #remove masked values from slice
            if mid_line.mask.size!=1:
                if True in mid_line.mask:
                    xarr=xarr[~mid_line.mask]
                    mid_line=mid_line[~mid_line.mask]

            #chunk array to better isolate artifacts
            chunk_size= 4
            for j in range(0, len(mid_line.tolist()), chunk_size):
                chunk = mid_line.tolist()[j:j+chunk_size]
                xarr_chunk = xarr[j:j+chunk_size]
                # make sure each iteration contains at least minimum number of elements
                if j==range(0, len(mid_line.tolist()), chunk_size)[-2] and len(mid_line.tolist()) % chunk_size != 0:
                    chunk = mid_line.tolist()[j:]
                    xarr_chunk = xarr[j:]
                #linear regression and get covariance
                slope, bias, rsquared, p_value, std_err = linregress(xarr_chunk,chunk)
                rsquaredarr.append(abs(rsquared)**2)
                std_errarr.append(std_err)
                #terminate early if last iteration would have small chunk size
                if len(chunk)>chunk_size:
                    break

            #exit loop/make plots in verbose mode if R^2 and standard error anomalous, or if on last iteration
            if (min(rsquaredarr) < 0.9 and max(std_errarr) > 0.01) or (i[0]==(len(eval('self.data_array_band%s'%(arrT)))-1)):
                if self.verbose:
                    #Make quality-control plots
                    import matplotlib.pyplot as plt
                    ax0=plt.figure().add_subplot(111)
                    ax0.scatter(xarr, mid_line, c='k', s=7)
                    refline = np.linspace(min(xarr),max(xarr),100)
                    ax0.plot(refline, (refline*slope)+bias, linestyle='solid', color='red')
                    ax0.set_ylabel('%s array'%(self.prod_key))
                    ax0.set_xlabel('distance')
                    ax0.set_title('Profile along %s'%(prof_direc))
                    ax0.annotate('R\u00b2 = %f\nStd error= %f'%(min(rsquaredarr),max(std_errarr)), (0, 1), xytext=(4, -4), xycoords='axes fraction', \
                        textcoords='offset points', fontweight='bold', ha='left', va='top')
                    if min(rsquaredarr) < 0.9 and max(std_errarr) > 0.01:
                        ax0.annotate('WARNING: R\u00b2 and standard error\nsuggest artifact exists', (1, 1), xytext=(4, -4), \
                            xycoords='axes fraction', textcoords='offset points', fontweight='bold', ha='right', va='top')
                    plt.margins(0)
                    plt.tight_layout()
                    plt.savefig(os.path.join(os.path.dirname(os.path.dirname(self.outname)),'metadatalyr_plots',self.prod_key, \
                        os.path.basename(self.outname)+'_%s.eps'%(prof_direc)))
                    plt.close()
                break

        return rsquaredarr, std_errarr

    def __run__(self):
        from scipy.linalg import lstsq

        # Get R^2/standard error across range
        rsquaredarr_rng, std_errarr_rng = self.__getCovar__('range')
        # Get R^2/standard error across azimuth
        rsquaredarr_az, std_errarr_az = self.__getCovar__('azimuth')

        #filter out normal values from arrays
        rsquaredarr = [0.97] ; std_errarr=[0.0015]
        if min(rsquaredarr_rng) < 0.97 and max(std_errarr_rng) > 0.0015:
            rsquaredarr.append(min(rsquaredarr_rng))
            std_errarr.append(max(std_errarr_rng))
        if min(rsquaredarr_az) < 0.97 and max(std_errarr_az) > 0.0015:
            rsquaredarr.append(min(rsquaredarr_az))
            std_errarr.append(max(std_errarr_az))

        #if R^2 and standard error anomalous, fix array
        if min(rsquaredarr) < 0.97 and max(std_errarr) > 0.0015:
            #Cycle through each band
            for i in range(1,5):
                self.data_array_band=self.data_array.GetRasterBand(i).ReadAsArray()
                #mask by nodata value
                self.data_array_band=np.ma.masked_where(self.data_array_band == self.data_array.GetRasterBand(i).GetNoDataValue(), self.data_array_band)
                negs_percent=((self.data_array_band < 0).sum()/self.data_array_band.size)*100
                #Unique bug-fix for bPerp layers with sign-flips
                if (self.prod_key=='bPerpendicular' and min(rsquaredarr) < 0.8 and max(std_errarr) > 0.1) \
                    and (negs_percent != 100 or negs_percent != 0):
                    #Circumvent Bperp sign-flip bug by comparing percentage of positive and negative values
                    self.data_array_band=abs(self.data_array_band)
                    if negs_percent>50:
                        self.data_array_band*=-1
                else:
                    # regular grid covering the domain of the data
                    X,Y = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1), np.arange(0, self.data_array_band.shape[0], 1))
                    Xmask,Ymask = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1), np.arange(0, self.data_array_band.shape[0], 1))
                    # best-fit linear plane: for very large artifacts, must mask array for outliers to get best fit
                    if min(rsquaredarr) < 0.85 and max(std_errarr) > 0.0015:
                        maj_percent=((self.data_array_band < self.data_array_band.mean()).sum()/self.data_array_band.size)*100
                        #mask all values above mean
                        if maj_percent>50:
                            self.data_array_band = np.ma.masked_where(self.data_array_band > self.data_array_band.mean(), self.data_array_band)
                        #mask all values below mean
                        else:
                            self.data_array_band = np.ma.masked_where(self.data_array_band < self.data_array_band.mean(), self.data_array_band)
                    # Mask columns/rows which are entirely made up of 0s
                    if self.data_array_band.mask.size!=1 and True in self.data_array_band.mask:
                        self.data_array_band, Xmask, Ymask = self.__truncateArray__(self.data_array_band, Xmask, Ymask)
                    # truncated grid covering the domain of the data
                    Xmask=Xmask[~self.data_array_band.mask]
                    Ymask=Ymask[~self.data_array_band.mask]
                    self.data_array_band = self.data_array_band[~self.data_array_band.mask]
                    XX = Xmask.flatten()
                    YY = Ymask.flatten()
                    A = np.c_[XX, YY, np.ones(len(XX))]
                    C,_,_,_ = lstsq(A, self.data_array_band.data.flatten())
                    # evaluate it on grid
                    self.data_array_band = C[0]*X + C[1]*Y + C[2]
                    #mask by nodata value
                    self.data_array_band=np.ma.masked_where(self.data_array_band == self.data_array.GetRasterBand(i).GetNoDataValue(), \
                        self.data_array_band)
                    np.ma.set_fill_value(self.data_array_band, self.data_array.GetRasterBand(i).GetNoDataValue())
                #update band
                self.data_array.GetRasterBand(i).WriteArray(self.data_array_band.filled())
                # Pass warning and get R^2/standard error across range/azimuth (only do for first band)
                if i==1:
                    # make sure appropriate unit is passed to print statement
                    lyrunit = "\N{DEGREE SIGN}"
                    if self.prod_key=='bPerpendicular' or self.prod_key=='bParallel':
                        lyrunit = 'm'
                    log.warning("%s layer for IFG %s has R\u00b2 of %.4f and standard error of %.4f%s, automated correction applied",
                                self.prod_key, os.path.basename(self.outname), min(rsquaredarr), max(std_errarr), lyrunit)
                    rsquaredarr_rng, std_errarr_rng = self.__getCovar__('range', profprefix='corrected')
                    rsquaredarr_az, std_errarr_az = self.__getCovar__('azimuth', profprefix='corrected')
        del self.data_array_band

        return self.data_array

def prep_dem(demfilename, bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr, proj, arrshape=None, workdir='./', outputFormat='ENVI', num_threads='2'):
    '''
        Function to load and export DEM, lat, lon arrays.
        If "Download" flag is specified, DEM will be donwloaded on the fly.
    '''
    # If specified DEM subdirectory exists, delete contents
    workdir      = os.path.join(workdir,'DEM')
    aria_dem     = os.path.join(workdir, 'SRTM_3arcsec.dem')
    os.makedirs(workdir, exist_ok=True)

    bounds       = open_shapefile(bbox_file, 0, 0).bounds # bounds of user bbox

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    outputFormat = 'ENVI' if outputFormat == 'VRT' else outputFormat

    if demfilename.lower()=='download':
        demfilename = dl_dem(aria_dem, prods_TOTbbox_metadatalyr, num_threads)

    else: # checks for user specified DEM, ensure it's georeferenced
        demfilename = os.path.abspath(demfilename)
        assert os.path.exists(demfilename), f'Cannot open DEM at: {demfilename}'
        ds_u = gdal.Open(demfilename)
        epsg = osr.SpatialReference(wkt=ds_u.GetProjection()).GetAttrValue('AUTHORITY',1)
        assert epsg is not None, f'No projection information in DEM: {demfilename}'
        del ds_u

    # write cropped DEM
    gdal.Warp(aria_dem, demfilename, format=outputFormat,
                cutlineDSName=prods_TOTbbox, outputBounds=bounds,
                outputType=gdal.GDT_Int16, width=arrshape[1], height=arrshape[0],
                multithread=True, options=['NUM_THREADS=%s'%(num_threads)])

    update_file = gdal.Open(aria_dem, gdal.GA_Update)
    update_file.SetProjection(proj); del update_file
    ds_aria     = gdal.Translate(f'{aria_dem}.vrt', aria_dem, format='VRT')
    log.info('Applied cutline to produce 3 arc-sec SRTM DEM: %s', aria_dem)

    # Load DEM and setup lat and lon arrays
    # pass expanded DEM for metadata field interpolation
    bounds  = list(open_shapefile(prods_TOTbbox_metadatalyr, 0, 0).bounds)
    gt      = ds_aria.GetGeoTransform()
    ds_aria = gdal.Warp('', aria_dem, format='MEM', outputBounds=bounds,
                             xRes=abs(gt[1]), yRes=abs(gt[-1]), multithread=True,
                                     options=['NUM_THREADS=%s'%(num_threads)])
    ds_aria.SetProjection(proj); ds_aria.SetDescription(aria_dem)

    # Define lat/lon arrays for fullres layers
    gt, xs, ys  = ds_aria.GetGeoTransform(), ds_aria.RasterXSize, ds_aria.RasterYSize
    Latitude    = np.linspace(gt[3], gt[3]+(gt[5]*ys), ys)
    Latitude    = np.repeat(Latitude[:, np.newaxis], xs, axis=1)
    Longitude   = np.linspace(gt[0], gt[0]+(gt[1]*xs), xs)
    Longitude   = np.repeat(Longitude[:, np.newaxis], ys, axis=1).T

    return aria_dem, ds_aria, Latitude, Longitude

def dl_dem(path_dem, path_prod_union, num_threads):
    """Download the DEM over product bbox union."""

    # Import functions
    from ARIAtools.shapefile_util import shapefile_area

    root      = os.path.splitext(path_dem)[0]
    prod_shapefile = open_shapefile(path_prod_union, 0, 0)
    WSEN      = prod_shapefile.bounds
    # If area > 225000 km2, must split requests into chunks to successfully access data
    chunk = False
    if shapefile_area(prod_shapefile) > 300000:
        chunk = True
        # Increase chunking size to discretize box into smaller grids
        log.warning('User-defined bounds results in an area of %d km which exceeds the maximum download area of 450000; downloading in chunks', shapefile_area(prod_shapefile))
        rows, cols = chunk_area(WSEN)

    if chunk: # Download in chunks (if necessary)
        chunked_files, k = [], 0
        for i in range(len(rows)-1):
            for j in range(len(cols)-1):
                dst = f'{root}_{k}_uncropped.tif'
                chunked_files.append(dst)
                WSEN = [cols[j], rows[i], cols[j+1], rows[i+1]]
                r    = requests.get(_world_dem.format(*WSEN), allow_redirects=True)
                with open(dst, 'wb') as fh:
                    fh.write(r.content)
                k+=1

        # Tile chunked products together after last iteration
        dst       = f'{root}_uncropped.tif'
        gdal.Warp(dst, chunked_files, multithread=True, options=[f'NUM_THREADS={num_threads}'])
        [os.remove(i) for i in glob.glob(f'{root}_*_uncropped.tif')] # remove tmp

    else:
        dst = f'{root}_uncropped.tif'
        r   = requests.get(_world_dem.format(*WSEN), allow_redirects=True)
        with open(dst, 'wb') as fh:
            fh.write(r.content)
    del r
    return dst

def merged_productbbox(metadata_dict, product_dict, workdir='./', bbox_file=None, croptounion=False, num_threads='2', minimumOverlap=0.0081, verbose=None):
    '''
        Extract/merge productBoundingBox layers for each pair and update dict, report common track bbox (default is to take common intersection, but user may specify union), report common track union to accurately interpolate metadata fields, and expected shape for DEM.
    '''

    # Import functions
    from ARIAtools.shapefile_util import save_shapefile, shapefile_area
    from shapely.geometry import Polygon

    # If specified workdir doesn't exist, create it
    os.makedirs(workdir, exist_ok=True)

    # If specified, check if user's bounding box meets minimum threshold area
    if bbox_file is not None:
        user_bbox=open_shapefile(bbox_file, 0, 0)
        overlap_area=shapefile_area(user_bbox)
        if overlap_area<minimumOverlap:
            raise Exception("User bound box '%s' has an area of only %fkm\u00b2, below specified minimum threshold area %fkm\u00b2"%(bbox_file,overlap_area,minimumOverlap))

    # Extract/merge productBoundingBox layers
    for scene in product_dict:
        # Get pair name, expected in dictionary
        pair_name=scene["pair_name"][0]
        outname=os.path.join(workdir, pair_name+'.json')

        # Create union of productBoundingBox layers
        for frame in scene["productBoundingBox"]:
            prods_bbox=open_shapefile(frame, 'productBoundingBox', 1)
            if os.path.exists(outname):
                union_bbox=open_shapefile(outname, 0, 0)
                prods_bbox=prods_bbox.union(union_bbox)
            save_shapefile(outname, prods_bbox, 'GeoJSON')              ##SS can we track and provide the proj information of the geojson?
        scene["productBoundingBox"]=[outname]

    prods_TOTbbox=os.path.join(workdir, 'productBoundingBox.json')
    # Need to track bounds of max extent to avoid metadata interpolation issues
    prods_TOTbbox_metadatalyr=os.path.join(workdir, 'productBoundingBox_croptounion_formetadatalyr.json')
    sceneareas=[open_shapefile(i['productBoundingBox'][0], 0, 0).area for i in product_dict]
    save_shapefile(prods_TOTbbox_metadatalyr, open_shapefile(product_dict[sceneareas.index(max(sceneareas))]['productBoundingBox'][0], 0, 0), 'GeoJSON')
    # Initiate intersection file with bbox, if bbox specified
    if bbox_file is not None:
        save_shapefile(prods_TOTbbox, open_shapefile(bbox_file, 0, 0), 'GeoJSON')
    # Intiate intersection with largest scene, if bbox NOT specified
    else:
        save_shapefile(prods_TOTbbox, open_shapefile(product_dict[sceneareas.index(max(sceneareas))]['productBoundingBox'][0], 0, 0), 'GeoJSON')
    rejected_scenes=[]
    for scene in product_dict:
        prods_bbox=open_shapefile(scene['productBoundingBox'][0], 0, 0)
        total_bbox=open_shapefile(prods_TOTbbox, 0, 0)
        total_bbox_metadatalyr=open_shapefile(prods_TOTbbox_metadatalyr, 0, 0)
        # Generate footprint for the union of all products
        if croptounion:
            # Get union
            total_bbox=total_bbox.union(prods_bbox)
            total_bbox_metadatalyr=total_bbox_metadatalyr.union(prods_bbox)
            # Save to file
            save_shapefile(prods_TOTbbox, total_bbox, 'GeoJSON')
            save_shapefile(prods_TOTbbox_metadatalyr, total_bbox_metadatalyr, 'GeoJSON')
        # Generate footprint for the common intersection of all products
        else:
            # Now pass track intersection for cutline
            prods_bbox=prods_bbox.intersection(total_bbox)
            # Estimate percentage of overlap with bbox
            if prods_bbox.bounds==():
                log.debug('Rejected scene %s has no common overlap with bbox'%(scene['productBoundingBox'][0]))
                rejected_scenes.append(product_dict.index(scene))
                os.remove(scene['productBoundingBox'][0])
            else:
                overlap_area=shapefile_area(prods_bbox)
                # Kick out scenes below specified overlap threshold
                if overlap_area<minimumOverlap:
                    log.debug("Rejected scene %s has only %fkm\u00b2 overlap with bbox"%(scene['productBoundingBox'][0],overlap_area))
                    rejected_scenes.append(product_dict.index(scene))
                    os.remove(scene['productBoundingBox'][0])
                else:
                    save_shapefile(prods_TOTbbox, prods_bbox, 'GeoJSON')
                    # Need to track bounds of max extent to avoid metadata interpolation issues
                    total_bbox_metadatalyr=total_bbox_metadatalyr.union(open_shapefile(scene['productBoundingBox'][0], 0, 0))
                    save_shapefile(prods_TOTbbox_metadatalyr, total_bbox_metadatalyr, 'GeoJSON')

    # Remove scenes with insufficient overlap w.r.t. bbox
    if rejected_scenes!=[]:
        log.info("%d out of %d interferograms rejected for not meeting specified spatial thresholds", len(rejected_scenes), len(product_dict))
    metadata_dict = [i for j, i in enumerate(metadata_dict) if j not in rejected_scenes]
    product_dict = [i for j, i in enumerate(product_dict) if j not in rejected_scenes]
    if product_dict==[]:
        raise Exception('No common track overlap, footprints cannot be generated.')

    # If bbox specified, intersect with common track intersection/union
    if bbox_file is not None:
        user_bbox=open_shapefile(bbox_file, 0, 0)
        total_bbox=open_shapefile(prods_TOTbbox, 0, 0)
        user_bbox=user_bbox.intersection(total_bbox)
        save_shapefile(prods_TOTbbox, user_bbox, 'GeoJSON')
    else:
        bbox_file=prods_TOTbbox

    # Warp the first scene with the output-bounds defined above
    # ensure output-bounds are an integer multiple of interferometric grid and adjust if necessary
    OG_bounds = list(open_shapefile(bbox_file, 0, 0).bounds)
    arrres = gdal.Open(product_dict[0]['unwrappedPhase'][0])
    arrres = [abs(arrres.GetGeoTransform()[1]), abs(arrres.GetGeoTransform()[-1])]
    ds = gdal.Warp('', gdal.BuildVRT('', product_dict[0]['unwrappedPhase'][0]), options=gdal.WarpOptions(format="MEM", \
        outputBounds = OG_bounds, xRes = arrres[0], yRes = arrres[1], targetAlignedPixels = True, multithread = True, \
        options = ['NUM_THREADS = %s'%(num_threads)]))
    # Get shape of full res layers
    arrshape=[ds.RasterYSize, ds.RasterXSize]
    new_bounds = [ds.GetGeoTransform()[0], ds.GetGeoTransform()[3] + (ds.GetGeoTransform()[-1] * arrshape[0]), \
        ds.GetGeoTransform()[0] + (ds.GetGeoTransform()[1] * arrshape[1]), ds.GetGeoTransform()[3]]
    if OG_bounds != new_bounds:
        # Use shapely to make list
        user_bbox = Polygon(np.column_stack((np.array([new_bounds[0],new_bounds[2],new_bounds[2],new_bounds[0],new_bounds[0]]),
                    np.array([new_bounds[1],new_bounds[1],new_bounds[3],new_bounds[3],new_bounds[1]])))) #Pass lons/lats to create polygon
        # Save polygon in shapefile
        bbox_file = os.path.join(os.path.dirname(workdir), 'user_bbox.json')
        save_shapefile(bbox_file, user_bbox, 'GeoJSON')
        total_bbox = open_shapefile(prods_TOTbbox, 0, 0)
        user_bbox = user_bbox.intersection(total_bbox)
        save_shapefile(prods_TOTbbox, user_bbox, 'GeoJSON')
    # Get projection of full res layers
    proj=ds.GetProjection()
    del ds

    return metadata_dict, product_dict, bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr, arrshape, proj

def export_products(full_product_dict, bbox_file, prods_TOTbbox, layers, rankedResampling=False, dem=None, lat=None, lon=None, mask=None, outDir='./',outputFormat='VRT', stitchMethodType='overlap', verbose=None, num_threads='2', multilooking=None):
    """
        Export layer and 2D meta-data layers (at the product resolution).
        The function finalize_metadata is called to derive the 2D metadata layer. Dem/lat/lon arrays must be passed for this process.
        The keys specify which layer to extract from the dictionary.
        All products are cropped by the bounds from the input bbox_file, and clipped to the track extent denoted by the input prods_TOTbbox.
        Optionally, a user may pass a mask-file.
    """

    ##Import functions
    from ARIAtools.vrtmanager import resampleRaster

    if not layers: return # only bbox

    bounds=open_shapefile(bbox_file, 0, 0).bounds
    if dem is not None:
        dem_bounds=[dem.GetGeoTransform()[0],dem.GetGeoTransform()[3]+(dem.GetGeoTransform()[-1]*dem.RasterYSize),dem.GetGeoTransform()[0]+(dem.GetGeoTransform()[1]*dem.RasterXSize),dem.GetGeoTransform()[3]]
    # Loop through user expected layers
    for key in layers:
        product_dict=[[j[key] for j in full_product_dict], [j["pair_name"] for j in full_product_dict]]
        workdir=os.path.join(outDir,key)

        ##Progress bar
        from ARIAtools import progBar
        prog_bar = progBar.progressBar(maxValue=len(product_dict[0]),prefix='Generating: '+key+' - ')

        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        # Mask specified, so file must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
        if outputFormat=='VRT' and mask is not None:
           outputFormat='ENVI'

        # Iterate through all IFGs
        for i in enumerate(product_dict[0]):
            outname=os.path.abspath(os.path.join(workdir, product_dict[1][i[0]][0]))
            ##Update progress bar
            prog_bar.update(i[0]+1,suffix=product_dict[1][i[0]][0])

            # Extract/crop metadata layers
            if any(":/science/grids/imagingGeometry" in s for s in [i[1]][0]):
                #create directory for quality control plots
                if verbose and not os.path.exists(os.path.join(outDir,'metadatalyr_plots',key)):
                    os.makedirs(os.path.join(outDir,'metadatalyr_plots',key))

                #make VRT pointing to metadata layers in standard product
                gdal.BuildVRT(outname +'.vrt', [i[1]][0])

                if dem is None:
                    raise Exception('No DEM input specified. Cannot extract 3D imaging geometry layers without DEM to intersect with.')

                # Check if height layers are consistent, and if not exit with error
                if len(set([gdal.Open(i).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES') for i in [i[1]][0]]))==1:
                    gdal.Open(outname+'.vrt').SetMetadataItem('NETCDF_DIM_heightsMeta_VALUES',gdal.Open([i[1]][0][0]).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES'))
                else:
                    raise Exception('Inconsistent heights for metadata layer(s) ', [i[1]][0], ' corresponding heights: ', [gdal.Open(i).GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES') for i in [i[1]][0]])

                # Pass metadata layer VRT, with DEM filename and output name to interpolate/intersect with DEM before cropping
                finalize_metadata(outname, bounds, dem_bounds, prods_TOTbbox, dem, lat, lon, mask, outputFormat, verbose=verbose)

            # Extract/crop full res layers, except for "unw" and "conn_comp" which requires advanced stitching
            elif key!='unwrappedPhase' and key!='connectedComponents':
                if outputFormat=='VRT' and mask is None:
                    # building the virtual vrt
                    gdal.BuildVRT(outname+ "_uncropped" +'.vrt', [i[1]][0])
                    # building the cropped vrt
                    gdal.Warp(outname+'.vrt', outname+"_uncropped"+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))
                else:
                    # building the VRT
                    gdal.BuildVRT(outname +'.vrt', [i[1]][0])
                    gdal.Warp(outname, outname+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))

                    # Update VRT
                    gdal.Translate(outname+'.vrt', outname, options=gdal.TranslateOptions(format="VRT"))

                    # Apply mask (if specified).
                    if mask is not None:
                        update_file=gdal.Open(outname,gdal.GA_Update)
                        update_file.GetRasterBand(1).WriteArray(mask.ReadAsArray()*gdal.Open(outname+'.vrt').ReadAsArray())
                        del update_file

            # Extract/crop "unw" and "conn_comp" layers leveraging the two stage unwrapper
            else:
                # Check if unw phase and connected components are already generated
                if not os.path.exists(os.path.join(outDir,'unwrappedPhase',product_dict[1][i[0]][0])) or not os.path.exists(os.path.join(outDir,'connectedComponents',product_dict[1][i[0]][0])):
                    # extract the inputs required for stitching of unwrapped and connected component files
                    unw_files = full_product_dict[i[0]]['unwrappedPhase']
                    conn_files = full_product_dict[i[0]]['connectedComponents']
                    prod_bbox_files = full_product_dict[i[0]]['productBoundingBoxFrames']
                    # based on the key define the output directories
                    outFileUnw=os.path.join(outDir,'unwrappedPhase',product_dict[1][i[0]][0])
                    outFileConnComp=os.path.join(outDir,'connectedComponents',product_dict[1][i[0]][0])

                    # calling the stitching methods
                    if stitchMethodType == 'overlap':
                        product_stitch_overlap(unw_files,conn_files,prod_bbox_files,bounds,prods_TOTbbox, outFileUnw=outFileUnw,outFileConnComp= outFileConnComp, mask=mask,outputFormat = outputFormat,verbose=verbose)
                    elif stitchMethodType == '2stage':
                        product_stitch_2stage(unw_files,conn_files,bounds,prods_TOTbbox,outFileUnw=outFileUnw,outFileConnComp= outFileConnComp, mask=mask,outputFormat = outputFormat,verbose=verbose)

                    #If necessary, resample both unw/conn_comp files
                    if multilooking is not None:
                        resampleRaster(outFileConnComp, multilooking, bounds, prods_TOTbbox, rankedResampling, outputFormat=outputFormat, num_threads=num_threads)

            #If necessary, resample raster
            if multilooking is not None and key!='unwrappedPhase' and key!='connectedComponents':
                resampleRaster(outname, multilooking, bounds, prods_TOTbbox, rankedResampling, outputFormat=outputFormat, num_threads=num_threads)

        prog_bar.close()
    return

def finalize_metadata(outname, bbox_bounds, dem_bounds, prods_TOTbbox, dem, lat, lon, mask=None, outputFormat='ENVI', verbose=None, num_threads='2'):
    '''
        2D metadata layer is derived by interpolating and then intersecting 3D layers with a DEM. Lat/lon arrays must also be passed for this process.
    '''

    # import dependencies
    import scipy

    # Import functions
    from ARIAtools.vrtmanager import renderVRT

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
        outputFormat='ENVI'

    # get final shape
    arrshape=gdal.Open(dem.GetDescription()).ReadAsArray().shape
    # load layered metadata array
    data_array=gdal.Warp('', outname+'.vrt', options=gdal.WarpOptions(format="MEM", multithread=True, options=['NUM_THREADS=%s'%(num_threads)]))

    #metadata layer quality check, correction applied if necessary
    data_array = metadata_qualitycheck(data_array, os.path.basename(os.path.dirname(outname)), outname, verbose).data_array

    # Define lat/lon/height arrays for metadata layers
    heightsMeta=np.array(gdal.Open(outname+'.vrt').GetMetadataItem('NETCDF_DIM_heightsMeta_VALUES')[1:-1].split(','), dtype='float32')
    ##SS Do we need lon lat if we would be doing gdal reproject using projection and transformation? See our earlier discussions.
    latitudeMeta=np.linspace(data_array.GetGeoTransform()[3],data_array.GetGeoTransform()[3]+(data_array.GetGeoTransform()[5]*data_array.RasterYSize),data_array.RasterYSize)
    longitudeMeta=np.linspace(data_array.GetGeoTransform()[0],data_array.GetGeoTransform()[0]+(data_array.GetGeoTransform()[1]*data_array.RasterXSize),data_array.RasterXSize)

    # First, using the height/latitude/longitude arrays corresponding to the metadata layer, set-up spatial 2D interpolator. Using this, perform vertical 1D interpolation on cube, and then use result to set-up a regular-grid-interpolator. Using this, pass DEM and full-res lat/lon arrays in order to get intersection with DEM.

    # 2D interpolation
    #mask by nodata value
    interp_2d = InterpCube(np.ma.masked_where(data_array.ReadAsArray() == data_array.GetRasterBand(1).GetNoDataValue(), \
        data_array.ReadAsArray()),heightsMeta,np.flip(latitudeMeta, axis=0),longitudeMeta)
    out_interpolated=np.zeros((heightsMeta.shape[0],latitudeMeta.shape[0],longitudeMeta.shape[0]))

    # 3D interpolation
    for hgt in enumerate(heightsMeta):
        for line in enumerate(latitudeMeta):
            for pixel in enumerate(longitudeMeta):
                out_interpolated[hgt[0], line[0], pixel[0]] = interp_2d(line[1], pixel[1], hgt[1])
    out_interpolated=np.flip(out_interpolated, axis=0)
    # interpolate to interferometric grid
    interpolator = scipy.interpolate.RegularGridInterpolator((heightsMeta,np.flip(latitudeMeta, axis=0),longitudeMeta), out_interpolated, method='linear', fill_value=data_array.GetRasterBand(1).GetNoDataValue(), bounds_error=False)

    try:
        out_interpolated = interpolator(np.stack((np.flip(dem, axis=0), lat, lon), axis=-1))
    except:
        #chunk data to conserve memory
        out_interpolated = []
        # need to mask nodata
        dem_array = dem.ReadAsArray()
        dem_array = np.ma.masked_where(dem_array == dem.GetRasterBand(1).GetNoDataValue(), dem_array)
        dem_array=np.array_split(dem.ReadAsArray(), 100) ; dem_array=[x for x in dem_array if x.size > 0]
        lat=np.array_split(lat, 100) ; dem_array=[x for x in lat if x.size > 0]
        lon=np.array_split(lon, 100) ; dem_array=[x for x in lon if x.size > 0]
        for i in enumerate(dem_array):
            out_interpolated.append(interpolator(np.stack((np.flip(i[1], axis=0), lat[i[0]], lon[i[0]]), axis=-1)))
        out_interpolated=np.concatenate(out_interpolated, axis=0)
        del dem_array

    # Save file
    renderVRT(outname+'_temp', out_interpolated, geotrans=dem.GetGeoTransform(), drivername=outputFormat, gdal_fmt=data_array.ReadAsArray().dtype.name, proj=dem.GetProjection(), nodata=data_array.GetRasterBand(1).GetNoDataValue())

    # Since metadata layer extends at least one grid node outside of the expected track bounds, it must be cut to conform with these bounds.
    # Crop to track extents
    gdal.Warp(outname, outname+'_temp', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bbox_bounds, dstNodata=data_array.GetRasterBand(1).GetNoDataValue(), width=arrshape[1], height=arrshape[0], multithread=True, options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))
    #remove temp files
    for i in glob.glob(outname+'_temp*'): os.remove(i)

    # Update VRT
    gdal.Translate(outname+'.vrt', outname, options=gdal.TranslateOptions(format="VRT"))

    # Apply mask (if specified)
    if mask is not None:
        out_interpolated = gdal.Open(outname).ReadAsArray()
        out_interpolated = mask.ReadAsArray()*out_interpolated
        # Update VRT with new raster
        update_file=gdal.Open(outname,gdal.GA_Update)
        update_file.GetRasterBand(1).WriteArray(out_interpolated)
        del update_file

    del out_interpolated, interpolator, interp_2d, data_array

def tropo_correction(full_product_dict, tropo_products, bbox_file, prods_TOTbbox, outDir='./',outputFormat='VRT', verbose=None, num_threads='2'):
    """
        Perform tropospheric corrections. Must provide valid path to GACOS products.
        All products are cropped by the bounds from the input bbox_file, and clipped to the track extent denoted by the input prods_TOTbbox.
    """

    # Import functions
    from ARIAtools.vrtmanager import renderVRT
    import tarfile
    from datetime import datetime
    from ARIAtools.shapefile_util import save_shapefile
    from shapely.geometry import Polygon

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    outputFormat = 'ENVI' if outputFormat == 'VRT' else outputFormat

    user_bbox    = open_shapefile(bbox_file, 0, 0)
    bounds       = user_bbox.bounds

    product_dict=[[j['unwrappedPhase'] for j in full_product_dict[1]], [j['lookAngle'] for j in full_product_dict[1]], [j["pair_name"] for j in full_product_dict[1]]]
    metadata_dict=[[j['azimuthZeroDopplerMidTime'] for j in full_product_dict[0]], [j['wavelength'] for j in full_product_dict[0]]]
    workdir=os.path.join(outDir,'tropocorrected_products')

    # If specified workdir doesn't exist, create it
    os.makedirs(workdir, exist_ok=True)

    # Get list of all dates for which standard products exist
    date_list=[]
    for i in product_dict[2]:
        date_list.append(i[0][:8]); date_list.append(i[0][9:])
    date_list=list(set(date_list)) # trim to unique dates only

    ### Determine if file input is single file, a list, or wildcard.
    # If list of files
    if len([str(val) for val in tropo_products.split(',')])>1:
        tropo_products=[str(i) for i in tropo_products.split(',')]
        # If wildcard
        tropo_products=[os.path.abspath(item) for sublist in [glob.glob(os.path.expanduser(os.path.expandvars(i))) if '*' in i else [i] for i in tropo_products] for item in sublist]
    # If single file or wildcard
    else:
        # If single file
        if os.path.isfile(tropo_products):
            tropo_products=[tropo_products]
        # If wildcard
        else:
            tropo_products=glob.glob(os.path.expanduser(os.path.expandvars(tropo_products)))
        # Convert relative paths to absolute paths
        tropo_products=[os.path.abspath(i) for i in tropo_products]
    if len(tropo_products)==0:
        raise Exception('No file match found')

    ###Extract tarfiles
    # Setup dictionary to track for products that are to be merged
    tropo_date_dict={}
    for i in date_list: tropo_date_dict[i]=[] ; tropo_date_dict[i+"_UTC"]=[]
    for i in enumerate(tropo_products):
        if not os.path.isdir(i[1]) and not i[1].endswith('.ztd') and not i[1].endswith('.tif'):
            untar_dir=os.path.join(os.path.abspath(os.path.join(i[1], os.pardir)),os.path.basename(i[1]).split('.')[0]+'_extracted')
            if not tarfile.is_tarfile(i[1]):
                raise Exception('Cannot extract %s because it is not a valid tarfile. Resolve this and relaunch'%(i[1]))
            log.info('Extracting GACOS tarfile %s to %s.', os.path.basename(i[1]), untar_dir)
            tarfile.open(i[1]).extractall(path=untar_dir)
            tropo_products[i[0]]=untar_dir
        # Loop through each GACOS product file, differentiating between direct input list of GACOS products vs parent directory
        if i[1].endswith('.ztd') or i[1].endswith('.tif'):
            ztd_list = [i[1]]
        else:
            ztd_list = glob.glob(os.path.join(tropo_products[i[0]],'*.ztd')) + glob.glob(os.path.join(tropo_products[i[0]],'*.tif'))
            # prioritize .ztd over .tif duplicates if the former exists
            ztd_list = [i for i in ztd_list if not (i.endswith('.tif') and i.split('.tif')[0] in ztd_list)]
        for k in ztd_list:
            # Only check files corresponding to standard product dates
            if os.path.basename(k)[:-4] in date_list:
                tropo_date_dict[os.path.basename(k)[:-4]].append(k)
                tropo_date_dict[os.path.basename(k)[:-4]+"_UTC"].append(os.path.basename(k)[:4]+'-'+os.path.basename(k)[4:6]+'-'+os.path.basename(k)[6:8]+'-'+open(k+'.rsc', 'r').readlines()[-1].split()[1])
                # make corresponding VRT file, if it doesn't exist
                if not os.path.exists(k+'.vrt'):
                    tropo_rsc_dict={}
                    for line in open(k+'.rsc', 'r').readlines(): tropo_rsc_dict[line.split()[0]]=line.split()[1]
                    gacos_prod=np.fromfile(k, dtype='float32').reshape(int(tropo_rsc_dict['FILE_LENGTH']),int(tropo_rsc_dict['WIDTH']))
                    # Save as GDAL file, using proj from first unwrappedPhase file
                    renderVRT(k, gacos_prod, geotrans=(float(tropo_rsc_dict['X_FIRST']), float(tropo_rsc_dict['X_STEP']), 0.0, float(tropo_rsc_dict['Y_FIRST']), 0.0, float(tropo_rsc_dict['Y_STEP'])), drivername=outputFormat, gdal_fmt='float32', proj=gdal.Open(os.path.join(outDir,'unwrappedPhase',product_dict[2][0][0])).GetProjection(), nodata=0.)
                    gacos_prod = None
                    log.debug("GACOS product %s successfully converted to GDAL-readable raster", k)

    # If multiple GACOS directories, merge products.
    tropo_products = list(set([os.path.dirname(i) if (i.endswith('.ztd') or i.endswith('.tif')) else i for i in tropo_products]))
    if len(tropo_products)>1:
        tropo_products=os.path.join(outDir,'merged_GACOS')
        log.info('Stitching/storing GACOS products in %s.', tropo_products)
        # If specified merged directory doesn't exist, create it
        if not os.path.exists(os.path.join(outDir,'merged_GACOS')):
            os.mkdir(os.path.join(outDir,'merged_GACOS'))

        for i in tropo_date_dict:
            # only make rsc/vrt files if valid product
            if 'UTC' not in i and tropo_date_dict[i]!=[]:
                outname=os.path.join(outDir,'merged_GACOS',i+'.ztd.vrt')
                # building the VRT
                gdal.BuildVRT(outname, tropo_date_dict[i])
                geotrans=gdal.Open(outname).GetGeoTransform()
                # Create merge rsc file
                with open(outname[:-4]+'.rsc','w') as merged_rsc:
                    merged_rsc.write('WIDTH %s\n'%(gdal.Open(outname).ReadAsArray().shape[1])) ; merged_rsc.write('FILE_LENGTH %s\n'%(gdal.Open(outname).ReadAsArray().shape[0]))
                    merged_rsc.write('XMIN %s\n'%(0)) ; merged_rsc.write('XMAX %s\n'%(gdal.Open(outname).ReadAsArray().shape[1]))
                    merged_rsc.write('YMIN %s\n'%(0)) ; merged_rsc.write('YMAX %s\n'%(gdal.Open(outname).ReadAsArray().shape[0]))
                    merged_rsc.write('X_FIRST %f\n'%(geotrans[0])) ; merged_rsc.write('Y_FIRST %f\n'%(geotrans[3]))
                    merged_rsc.write('X_STEP %f\n'%(geotrans[1])) ; merged_rsc.write('Y_STEP %f\n'%(geotrans[-1]))
                    merged_rsc.write('X_UNIT %s\n'%('degres')) ; merged_rsc.write('Y_UNIT %s\n'%('degres'))
                    merged_rsc.write('Z_OFFSET %s\n'%(0)) ; merged_rsc.write('Z_SCALE %s\n'%(1))
                    merged_rsc.write('PROJECTION %s\n'%('LATLON')) ; merged_rsc.write('DATUM %s\n'%('WGS84'))
                    merged_rsc.write('TIME_OF_DAY %s\n'%(''.join(tropo_date_dict[i+"_UTC"])))
    else:
        tropo_products=tropo_products[0]

    # Estimate percentage of overlap with tropospheric product
    for i in glob.glob(os.path.join(tropo_products,'*.vrt')):
        # create shapefile
        geotrans=gdal.Open(i).GetGeoTransform()
        bbox=[geotrans[3]+(gdal.Open(i).ReadAsArray().shape[0]*geotrans[-1]),geotrans[3],geotrans[0],geotrans[0]+(gdal.Open(i).ReadAsArray().shape[1]*geotrans[1])]
        bbox=Polygon(np.column_stack((np.array([bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]),
                            np.array([bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]))))
        save_shapefile(i+'.json', bbox, 'GeoJSON')
        per_overlap=((user_bbox.intersection(open_shapefile(i+'.json', 0, 0)).area)/(user_bbox.area))*100
        if per_overlap!=100. and per_overlap!=0.:
            log.warning("Common track extent only has %d overlap with tropospheric product %d\n", per_overlap, i)
        if per_overlap==0.:
            raise Exception('No spatial overlap between tropospheric product %s and defined bounding box. Resolve conflict and relaunch', i)

    # Iterate through all IFGs and apply corrections
    missing_products = []
    for i in range(len(product_dict[0])):
        outname=os.path.join(workdir,product_dict[2][i][0])
        unwname=os.path.join(outDir,'unwrappedPhase',product_dict[2][i][0])
        tropo_reference=os.path.join(tropo_products,product_dict[2][i][0][:8]+'.ztd.vrt')
        tropo_secondary=os.path.join(tropo_products,product_dict[2][i][0][9:]+'.ztd.vrt')
        # if .ztd products don't exist, check if .tif exists
        if not os.path.exists(tropo_reference):
            tropo_reference=os.path.join(tropo_products,product_dict[2][i][0][:8]+'.ztd.tif.vrt')
        if not os.path.exists(tropo_secondary):
            tropo_secondary=os.path.join(tropo_products,product_dict[2][i][0][9:]+'.ztd.tif.vrt')
        # skip if corrected already generated and does not need to be updated
        if os.path.exists(outname):
            # get unwrappedPhase geotrans and productbounding box
            unw_prodcheck = gdal.Open(unwname)
            unw_geotrans = unw_prodcheck.GetGeoTransform()
            unw_prodcheck = np.isfinite(unw_prodcheck.ReadAsArray())
            tropo_prodcheck = gdal.Open(outname)
            output_geotrans = tropo_prodcheck.GetGeoTransform()
            tropo_prodcheck = np.isfinite(tropo_prodcheck.ReadAsArray())
            if unw_geotrans == output_geotrans and np.array_equal(unw_prodcheck, tropo_prodcheck):
                continue
            del unw_prodcheck, tropo_prodcheck
        if os.path.exists(tropo_reference) and os.path.exists(tropo_secondary):
            # Check if tropo products are temporally consistent with IFG
            for j in [tropo_reference, tropo_secondary]:
                # Get ARIA product times
                aria_rsc_dict={}
                aria_rsc_dict['azimuthZeroDopplerMidTime']=[datetime.strptime(os.path.basename(j)[:4]+'-'+os.path.basename(j)[4:6]+'-'+os.path.basename(j)[6:8]+'-'+m[11:], "%Y-%m-%d-%H:%M:%S.%f") for m in metadata_dict[0][0]]
                # Get tropo product UTC times
                tropo_rsc_dict={}
                tropo_rsc_dict['TIME_OF_DAY']=open(j[:-4]+'.rsc', 'r').readlines()[-1].split()[1].split('UTC')[:-1]
                # If stitched tropo product, must account for date change (if applicable)
                if '-' in tropo_rsc_dict['TIME_OF_DAY'][0]:
                    tropo_rsc_dict['TIME_OF_DAY'] = [datetime.strptime(m[:10]+'-'+m[11:].split('.')[0]+'-'+str(round(float('0.'+m[11:].split('.')[1])*60)), "%Y-%m-%d-%H-%M") for m in tropo_rsc_dict['TIME_OF_DAY']]
                else:
                    tropo_rsc_dict['TIME_OF_DAY'] = [datetime.strptime(os.path.basename(j)[:4]+'-'+os.path.basename(j)[4:6]+'-'+os.path.basename(j)[6:8]+'-'+tropo_rsc_dict['TIME_OF_DAY'][0].split('.')[0]+'-'+str(round(float('0.'+tropo_rsc_dict['TIME_OF_DAY'][0].split('.')[-1])*60)), "%Y-%m-%d-%H-%M")]

                # Check and report if tropospheric product falls outside of standard product range
                latest_start = max(aria_rsc_dict['azimuthZeroDopplerMidTime']+[min(tropo_rsc_dict['TIME_OF_DAY'])])
                earliest_end = min(aria_rsc_dict['azimuthZeroDopplerMidTime']+[max(tropo_rsc_dict['TIME_OF_DAY'])])
                delta = (earliest_end - latest_start).total_seconds() + 1
                if delta<0:
                    log.warning("tropospheric product was generated %f secs outside of acquisition interval for scene %s in IFG %s",
                                        abs(delta), os.path.basename(j)[:8], product_dict[2][i][0])

            # Open unwrappedPhase and mask nodata
            unwphase=gdal.Open(unwname)
            geotrans=unwphase.GetGeoTransform() ; proj=unwphase.GetProjection() ; unwnodata=unwphase.GetRasterBand(1).GetNoDataValue()
            unwphase=unwphase.ReadAsArray()
            arrshape=unwphase.shape # shape of output array to match ifgs
            unwphase=np.ma.masked_where(unwphase == unwnodata, unwphase)

            # Open corresponding tropo products and pass the difference
            tropo_product = gdal.Warp('', tropo_reference, format="MEM", cutlineDSName=prods_TOTbbox,
                          outputBounds=bounds, width=arrshape[1], height=arrshape[0], resampleAlg='lanczos',
                          dstNodata=0., multithread=True, options=['NUM_THREADS=%s'%(num_threads)]).ReadAsArray()

            tropo_product = np.ma.masked_where(tropo_product == 0., tropo_product)
            tropo_secondary = gdal.Warp('', tropo_secondary, format="MEM", cutlineDSName=prods_TOTbbox,
                                            outputBounds=bounds, width=arrshape[1], height=arrshape[0],
                                            resampleAlg='lanczos', dstNodata=0., multithread=True,
                                            options=['NUM_THREADS=%s'%(num_threads)]).ReadAsArray()

            tropo_secondary = np.ma.masked_where(tropo_secondary == 0., tropo_secondary)
            tropo_product   = np.subtract(tropo_secondary,tropo_product)

            # Convert troposphere to rad
            tropo_product=np.divide(tropo_product,float(metadata_dict[1][i][0])/(4*np.pi))
            # Account for lookAngle
            # if in TS mode, only 1 lookfile would be generated, so check for this
            if os.path.exists(os.path.join(outDir,'lookAngle',product_dict[2][i][0])):
                lookfile=gdal.Open(os.path.join(outDir,'lookAngle',product_dict[2][i][0])).ReadAsArray()
            else:
                lookfile=gdal.Open(os.path.join(outDir,'lookAngle',product_dict[2][0][0])).ReadAsArray()
            lookfile=np.sin(np.deg2rad(np.ma.masked_where(lookfile == 0., lookfile)))
            tropo_product=np.divide(tropo_product,lookfile)

            #Correct phase and save to file
            unwphase=np.subtract(unwphase,tropo_product)
            np.ma.set_fill_value(unwphase, unwnodata); np.ma.set_fill_value(tropo_product, 0.)
            renderVRT(outname+'_tropodiff', tropo_product.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=0.)
            renderVRT(outname, unwphase.filled(), geotrans=geotrans, drivername=outputFormat, gdal_fmt='float32', proj=proj, nodata=unwnodata)

            del unwphase, tropo_product, tropo_reference, tropo_secondary, lookfile

        else:
            log.warning("Must skip IFG %s, because the tropospheric products corresponding to the reference and/or secondary products are not found in the specified folder %s", product_dict[2][i][0], tropo_products)
            for j in [tropo_reference, tropo_secondary]:
                if not os.path.exists(j) and j not in missing_products:
                    missing_products.append(j)
    # Print list of dates missing tropospheric corrections
    if len(missing_products) > 0:
        missing_products = [os.path.basename(i)[:8] for i in missing_products]
        log.debug("Tropo products for the following dates are missing:")
        log.debug(missing_products)

def main(inps=None):
    '''
       Main workflow for extracting layers from ARIA products
    '''
    from ARIAtools.ARIAProduct import ARIA_standardproduct
    print ('**********************************************************************')
    print (version.release_description)
    print ('*** Extract Product Function ***')
    print ('**********************************************************************')

    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    if not inps.layers and not inps.tropo_products:
        log.info('No layers specified; only creating bounding box shapes')

    elif inps.tropo_products:
        log.info('Tropospheric corrections will be applied, making sure at least unwrappedPhase and lookAngle are extracted.')
        # If no input layers specified, initialize list
        if not inps.layers: inps.layers=[]
        # If valid argument for input layers passed, parse to list
        if isinstance(inps.layers,str): inps.layers=list(inps.layers.split(',')) ; inps.layers=[i.replace(' ','') for i in inps.layers]
        if 'lookAngle' not in inps.layers: inps.layers.append('lookAngle')
        if 'unwrappedPhase' not in inps.layers: inps.layers.append('unwrappedPhase')

    elif inps.layers.lower()=='all':
        log.info('All layers are to be extracted, pass all keys.')
        inps.layers=list(standardproduct_info.products[1][0].keys())
        # Must remove productBoundingBoxes & pair-names because they are not raster layers
        inps.layers=[i for i in inps.layers if i not in ['productBoundingBox','productBoundingBoxFrames','pair_name']]

    else:
        inps.layers=list(inps.layers.split(','))
        inps.layers=[i.replace(' ','') for i in inps.layers]

    # pass number of threads for gdal multiprocessing computation
    if inps.num_threads.lower()=='all':
        import multiprocessing
        log.info('User specified use of all %s threads for gdal multiprocessing', str(multiprocessing.cpu_count()))
        inps.num_threads='ALL_CPUS'
    log.info('Thread count specified for gdal multiprocessing = %s', inps.num_threads)

    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[0], standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr, arrshape, proj = merged_productbbox(standardproduct_info.products[0], standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion, num_threads=inps.num_threads, minimumOverlap=inps.minimumOverlap, verbose=inps.verbose)
    # Load or download mask (if specified).
    if inps.mask is not None:
        inps.mask = prep_mask([[item for sublist in [list(set(d['amplitude']))
                        for d in standardproduct_info.products[1] if 'amplitude' in d] for item in sublist],
                        [item for sublist in [list(set(d['pair_name'])) for d in
                        standardproduct_info.products[1] if 'pair_name' in d] for item in sublist]],
                        inps.mask, standardproduct_info.bbox_file, prods_TOTbbox, proj, amp_thresh=inps.amp_thresh, arrshape=arrshape,
                        workdir=inps.workdir, outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile,
                    standardproduct_info.bbox_file, prods_TOTbbox,
                    prods_TOTbbox_metadatalyr, proj, arrshape=arrshape,
                    workdir=inps.workdir, outputFormat=inps.outputFormat, num_threads=inps.num_threads)
    else:
        demfile, Latitude, Longitude = None, None, None

    # Extract user expected layers
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file,
            prods_TOTbbox, inps.layers, inps.rankedResampling, dem=demfile, lat=Latitude,
            lon=Longitude, mask=inps.mask, outDir=inps.workdir, outputFormat=inps.outputFormat,
            stitchMethodType='overlap', verbose=inps.verbose, num_threads=inps.num_threads, multilooking=inps.multilooking)

    # If necessary, resample DEM/mask AFTER they have been used to extract metadata layers and mask output layers, respectively
    if inps.multilooking is not None:
        # Import functions
        from ARIAtools.vrtmanager import resampleRaster
        bounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds
        # Resample mask
        if inps.mask is not None:
            resampleRaster(inps.mask.GetDescription(), inps.multilooking, bounds, prods_TOTbbox,
                        inps.rankedResampling, outputFormat=inps.outputFormat, num_threads=inps.num_threads)
        # Resample DEM
        if demfile is not None:
            resampleRaster(demfile.GetDescription(), inps.multilooking, bounds,
                            prods_TOTbbox, inps.rankedResampling, outputFormat=inps.outputFormat, num_threads=inps.num_threads)

    # Perform GACOS-based tropospheric corrections (if specified).
    if inps.tropo_products:
        tropo_correction(standardproduct_info.products, inps.tropo_products,
                         standardproduct_info.bbox_file, prods_TOTbbox, outDir=inps.workdir,
                         outputFormat=inps.outputFormat, verbose=inps.verbose, num_threads=inps.num_threads)
