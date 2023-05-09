#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Extract and organize specified layer(s).
If no layer is specified, extract product bounding box shapefile(s)
"""

import os
import numpy as np
from copy import deepcopy
import glob
from dem_stitcher.datasets import DATASETS
from dem_stitcher import get_dem_tile_paths
from osgeo import gdal, osr
import logging
import requests
import shutil
from ARIAtools.logger import logger

from ARIAtools.shapefile_util import open_shapefile, chunk_area
from ARIAtools.mask_util import prep_mask
from ARIAtools.unwrapStitching import product_stitch_overlap, \
                                      product_stitch_2stage
from ARIAtools.vrtmanager import renderVRT, resampleRaster, layerCheck, \
                                 get_basic_attrs
from ARIAtools.sequential_stitching import product_stitch_sequential, \
                                           product_stitch_sequential_metadata
import pyproj
from pyproj import CRS, Transformer
from rioxarray import open_rasterio
import rasterio as rio

gdal.UseExceptions()
#Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

log = logging.getLogger(__name__)



def createParser():
    """Extract specified product layers. The default will export all layers."""
    import argparse
    parser = argparse.ArgumentParser(description= \
            'Program to extract data and meta-data layers from ARIA standard GUNW products.'\
            ' Program will handle cropping/stitching when needed. '\
            'By default, the program will crop all IFGs to bounds determined by the common intersection and bbox (if specified)')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./',
        help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-gp', '--gacos_products', dest='gacos_products',
                        type=str, default=None, help='Path to director(ies) '
                        'or tar file(s) containing GACOS products.')
    parser.add_argument('-l', '--layers', dest='layers', default=None,
                        help='Specify layers to extract as a comma '
                        'deliminated list bounded by single quotes. '
                        'Allowed keys are: "unwrappedPhase", "coherence", '
                        '"amplitude", "bPerpendicular", "bParallel", '
                        '"incidenceAngle", "lookAngle", "azimuthAngle", '
                        '"ionosphere", "troposphereWet", '
                        '"troposphereHydrostatic", "troposphereTotal", '
                        '"solidEarthTide". '
                        'If "all" specified, then all layers are extracted. '
                        'If blank, will only extract bounding box.')
    parser.add_argument('-tm', '--tropo_models', dest='tropo_models',
                        type=str, default=None, help='Provide list of '
                        'weather models you wish to extract. Refer to '
                        'ARIA_TROPO_MODELS for list of supported models')
    parser.add_argument('-d', '--demfile', dest='demfile', type=str,
            default=None, help='DEM file. To download new DEM, specify "Download".')
    parser.add_argument('-p', '--projection', dest='projection', default='WGS84', type=str,
            help='projection for DEM. By default WGS84.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None,
        help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None,
        help="Path to mask file or 'Download'. File needs to be GDAL compatabile, contain spatial reference information, and have invalid/valid data represented by 0/1, respectively. If 'Download', will use GSHHS water mask. If 'NLCD', will mask classes 11, 12, 90, 95; see: https://www.mrlc.gov/national-land-cover-database-nlcd-201://www.mrlc.gov/national-land-cover-database-nlcd-2016")
    parser.add_argument('-at', '--amp_thresh', dest='amp_thresh', default=None, type=str,
        help='Amplitude threshold below which to mask. Specify "None" to not use amplitude mask. By default "None".')
    parser.add_argument('-nt', '--num_threads', dest='num_threads', default='2', type=str,
        help='Specify number of threads for multiprocessing operation in gdal. By default "2". Can also specify "All" to use all available threads.')
    parser.add_argument('-sm', '--stitchMethod', dest='stitchMethodType',  type=str, default='overlap', help="Method applied to stitch the unwrapped data. Allowed methods are: 'overlap', '2stage', and 'sequential'. 'overlap' - product overlap is minimized, '2stage' - minimization is done on connected components, 'sequential' - sequential minimization of all overlapping connected components.  Default is 'overlap'.")
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str, default='VRT',
        help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By default files are generated virtually except for "bPerpendicular", "bParallel", "incidenceAngle", "lookAngle","azimuthAngle", "unwrappedPhase" as these are require either DEM intersection or corrections to be applied')
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion',
        help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-ml', '--multilooking', dest='multilooking', type=int, default=None,
        help='Multilooking factor is an integer multiple of standard resolution. E.g. 2 = 90m*2 = 180m')
    parser.add_argument('-rr', '--rankedResampling', action='store_true', dest='rankedResampling',
        help="If turned on, IFGs resampled based off of the average of pixels in a given resampling window corresponding to the connected component mode (if multilooking specified). Program defaults to lanczos resampling algorithm through gdal (if multilooking specified).")
    parser.add_argument('-mo', '--minimumOverlap', dest='minimumOverlap', type=float, default=0.0081,
        help='Minimum km\u00b2 area of overlap of scenes wrt specified bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single pixel at standard 90m resolution"')
    parser.add_argument('--version', dest='version',  default=None,
                        help='Specify version as str, e.g. 2_0_4 or all prods; '
                        'default: all')
    parser.add_argument('--nc_version', dest='nc_version',  default='1b',
                        help='Specify netcdf version as str, '
                        'e.g. 1c or all prods;'
                        'default: 1b')
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose',
        help="Toggle verbose mode on.")

    return parser


def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)


class InterpCube(object):
    """Class to interpolate intersection of cube with DEM."""
    def __init__(self, inobj, hgtobj, latobj, lonobj):
        """Init with h5py dataset."""
        self.data = inobj[:]
        self.hgts = hgtobj[:]
        self.offset = None
        self.interp = []
        self.latobj = latobj[:]
        self.lonobj = lonobj[:]

        self.createInterp()

    def createInterp(self):
        """Create interpolators."""
        from scipy.interpolate import RectBivariateSpline
        self.offset = np.mean(self.data)
        for i in range(len(self.hgts)):
            self.interp.append(RectBivariateSpline(self.latobj, self.lonobj,
                                                    self.data[i]-self.offset))

    def __call__(self, line, pix, h):
        """Interpolate at a single point."""
        from scipy.interpolate import interp1d

        vals  = np.array( [x(line,pix)[0,0] for x in self.interp])
        est = interp1d(self.hgts, vals, kind='cubic')
        return est(h) + self.offset


class metadata_qualitycheck:
    """Metadata quality control function.
    Artifacts recognized based off of covariance of cross-profiles.
    Bug-fix varies based off of layer of interest.
    Verbose mode generates a series of quality control plots with
    these profiles.
    """

    def __init__(self, data_array, prod_key, outname, verbose=None):
        # Pass inputs
        self.data_array = data_array
        self.prod_key = prod_key
        self.outname = outname
        self.verbose = verbose
        self.data_array_band=data_array.GetRasterBand(1).ReadAsArray()
        #mask by nodata value
        self.data_array_band=np.ma.masked_where(self.data_array_band ==
                self.data_array.GetRasterBand(1).GetNoDataValue(),
                self.data_array_band)

        logger.setLevel(logging.DEBUG) if self.verbose else ''

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
            Xmask,Ymask = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1),
                                    np.arange(0, self.data_array_band.shape[0], 1))
            self.data_array_band, Xmask, Ymask = self.__truncateArray__(
                                            self.data_array_band, Xmask, Ymask)

        #append prefix for plot names
        prof_direc = profprefix + prof_direc

        #iterate through transpose of matrix if looking in azimuth
        arrT=''
        if 'azimuth' in prof_direc:
            arrT='.T'
        # Cycle between range and azimuth profiles
        rsquaredarr = []
        std_errarr  = []
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
                if j==range(0, len(mid_line.tolist()), chunk_size)[-2] and \
                                    len(mid_line.tolist()) % chunk_size != 0:
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
            if (min(rsquaredarr) < 0.9 and max(std_errarr) > 0.01) or \
                        (i[0]==(len(eval('self.data_array_band%s'%(arrT)))-1)):
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
                self.data_array_band=np.ma.masked_where(self.data_array_band == \
                    self.data_array.GetRasterBand(i).GetNoDataValue(),
                                        self.data_array_band)
                negs_percent=((self.data_array_band < 0).sum() \
                                /self.data_array_band.size)*100

                # Unique bug-fix for bPerp layers with sign-flips
                if (self.prod_key=='bPerpendicular' and min(rsquaredarr) < 0.8 \
                                    and max(std_errarr) > 0.1) \
                                and (negs_percent != 100 or negs_percent != 0):
                    #Circumvent Bperp sign-flip bug by comparing percentage of positive and negative values
                    self.data_array_band=abs(self.data_array_band)
                    if negs_percent>50:
                        self.data_array_band*=-1
                else:
                    # regular grid covering the domain of the data
                    X,Y = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1),
                                np.arange(0, self.data_array_band.shape[0], 1))
                    Xmask,Ymask = np.meshgrid(np.arange(0, self.data_array_band.shape[1], 1),
                                    np.arange(0, self.data_array_band.shape[0], 1))
                    # best-fit linear plane: for very large artifacts, must mask array for outliers to get best fit
                    if min(rsquaredarr) < 0.85 and max(std_errarr) > 0.0015:
                        maj_percent=((self.data_array_band < \
                                self.data_array_band.mean()).sum() \
                                / self.data_array_band.size)*100
                        #mask all values above mean
                        if maj_percent>50:
                            self.data_array_band = np.ma.masked_where(
                                self.data_array_band > self.data_array_band.mean(),
                                self.data_array_band)
                        #mask all values below mean
                        else:
                            self.data_array_band = np.ma.masked_where(
                                self.data_array_band < self.data_array_band.mean(),
                                self.data_array_band)
                    # Mask columns/rows which are entirely made up of 0s
                    if self.data_array_band.mask.size!=1 and \
                                    True in self.data_array_band.mask:
                        self.data_array_band, Xmask, Ymask = self.__truncateArray__(
                                    self.data_array_band, Xmask, Ymask)

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
                    self.data_array_band=np.ma.masked_where(
                        self.data_array_band == self.data_array.GetRasterBand(i
                                    ).GetNoDataValue(), self.data_array_band)

                    np.ma.set_fill_value(self.data_array_band,
                            self.data_array.GetRasterBand(i).GetNoDataValue())
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


def prep_dem(demfilename, bbox_file, prods_TOTbbox, prods_TOTbbox_metadatalyr,
                        proj, arrshape=None, workdir='./',
                        outputFormat='ENVI', num_threads='2', dem_name: str = 'glo_90'):
    """Function to load and export DEM, lat, lon arrays.
    If "Download" flag is specified, DEM will be downloaded on the fly.
    """
    # If specified DEM subdirectory exists, delete contents
    workdir      = os.path.join(workdir,'DEM')
    aria_dem     = os.path.join(workdir, f'{dem_name}.dem')
    os.makedirs(workdir, exist_ok=True)

    bounds       = open_shapefile(bbox_file, 0, 0).bounds # bounds of user bbox

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    outputFormat = 'ENVI' if outputFormat == 'VRT' else outputFormat

    if demfilename.lower()=='download':
        if dem_name not in DATASETS:
            raise ValueError(f'"dem_name" must be in {", ".join(DATASETS)}')
        demfilename = dl_dem(aria_dem, prods_TOTbbox_metadatalyr, num_threads, dem_name)

    else: # checks for user specified DEM, ensure it's georeferenced
        demfilename = os.path.abspath(demfilename)
        assert os.path.exists(demfilename), f'Cannot open DEM at: {demfilename}'
        ds_u = gdal.Open(demfilename, gdal.GA_ReadOnly)
        epsg = osr.SpatialReference(wkt=ds_u.GetProjection()).GetAttrValue('AUTHORITY',1)
        assert epsg is not None, f'No projection information in DEM: {demfilename}'
        del ds_u

    # write cropped DEM
    if demfilename == os.path.abspath(aria_dem):
        log.warning('The DEM you specified already exists in %s, '\
                'using the existing one...', os.path.dirname(aria_dem))
        ds_aria = gdal.Open(aria_dem, gdal.GA_ReadOnly)

    else:
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
                                         xRes=abs(gt[1]), yRes=abs(gt[-1]),
                     multithread=True, options=['NUM_THREADS=%s'%(num_threads)])
        ds_aria.SetProjection(proj); ds_aria.SetDescription(aria_dem)

    # Delete temporary dem-stitcher directory
    if os.path.exists(f'{dem_name}_tiles'):
        shutil.rmtree(f'{dem_name}_tiles')

    # Define lat/lon arrays for fullres layers
    gt, xs, ys  = ds_aria.GetGeoTransform(), ds_aria.RasterXSize, ds_aria.RasterYSize
    Latitude    = np.linspace(gt[3], gt[3]+(gt[5]*ys), ys)
    Latitude    = np.repeat(Latitude[:, np.newaxis], xs, axis=1)
    Longitude   = np.linspace(gt[0], gt[0]+(gt[1]*xs), xs)
    Longitude   = np.repeat(Longitude[:, np.newaxis], ys, axis=1).T

    return aria_dem, ds_aria, Latitude, Longitude


def dl_dem(path_dem, path_prod_union, num_threads, dem_name: str = 'glo_90'):
    """Download the DEM over product bbox union."""

    root  = os.path.splitext(path_dem)[0]
    prod_shapefile = open_shapefile(path_prod_union, 0, 0)
    extent = prod_shapefile.bounds

    dem_tile_paths = get_dem_tile_paths(bounds=extent,
                                        dem_name=dem_name,
                                        localize_tiles_to_gtiff=(False if dem_name == 'glo_30' else True),
                                        tile_dir=f'{dem_name}_tiles')

    vrt_path = f'{root}_uncropped.vrt'
    ds = gdal.BuildVRT(vrt_path, dem_tile_paths)
    ds = None

    return vrt_path


def merged_productbbox(metadata_dict, product_dict, workdir='./',
                       bbox_file=None, croptounion=False, num_threads='2',
                       minimumOverlap=0.0081, verbose=None):
    """Extract/merge productBoundingBox layers for each pair.
    Also update dict, report common track bbox
    (default is to take common intersection, but user may specify union),
    report common track union to accurately interpolate metadata fields,
    and expected shape for DEM.
    """
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
            raise Exception(f'User bound box {bbox_file} has an area of only '
                            f'{overlap_area}km\u00b2, below specified '
                            f'minimum threshold area '
                            f'{minimumOverlap}km\u00b2')

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
            save_shapefile(outname, prods_bbox, 'GeoJSON')
        scene["productBoundingBox"]=[outname]

    prods_TOTbbox=os.path.join(workdir, 'productBoundingBox.json')
    # Need to track bounds of max extent
    # to avoid metadata interpolation issues
    prods_TOTbbox_metadatalyr = os.path.join(workdir,
        'productBoundingBox_croptounion_formetadatalyr.json')
    sceneareas = [open_shapefile(i['productBoundingBox'][0], 0, 0).area \
                  for i in product_dict]
    save_shapefile(prods_TOTbbox_metadatalyr,
        open_shapefile(product_dict[sceneareas.index(max(
                sceneareas))]['productBoundingBox'][0], 0, 0), 'GeoJSON')
    # Initiate intersection file with bbox, if bbox specified
    if bbox_file is not None:
        save_shapefile(prods_TOTbbox, open_shapefile(bbox_file, 0, 0),
                       'GeoJSON')
    # Intiate intersection with largest scene, if bbox NOT specified
    else:
        save_shapefile(prods_TOTbbox, open_shapefile(
            product_dict[sceneareas.index(max(
                    sceneareas))]['productBoundingBox'][0], 0, 0), 'GeoJSON')
    rejected_scenes=[]
    for scene in product_dict:
        scene_obj = scene['productBoundingBox'][0]
        prods_bbox=open_shapefile(scene_obj, 0, 0)
        total_bbox=open_shapefile(prods_TOTbbox, 0, 0)
        total_bbox_metadatalyr=open_shapefile(prods_TOTbbox_metadatalyr, 0, 0)
        # Generate footprint for the union of all products
        if croptounion:
            # Get union
            total_bbox=total_bbox.union(prods_bbox)
            total_bbox_metadatalyr=total_bbox_metadatalyr.union(prods_bbox)
            # Save to file
            save_shapefile(prods_TOTbbox, total_bbox, 'GeoJSON')
            save_shapefile(prods_TOTbbox_metadatalyr,
                           total_bbox_metadatalyr, 'GeoJSON')
        # Generate footprint for the common intersection of all products
        else:
            # Now pass track intersection for cutline
            prods_bbox=prods_bbox.intersection(total_bbox)
            # Estimate percentage of overlap with bbox
            if prods_bbox.geom_type == 'MultiPolygon':
                log.debug(f'Rejected scene {scene_obj} is type MultiPolygon')
                rejected_scenes.append(product_dict.index(scene))
                os.remove(scene_obj)
                continue
            if prods_bbox.bounds==() or prods_bbox.is_empty:
                log.debug(f'Rejected scene {scene_obj} '
                          f'has no common overlap with bbox')
                rejected_scenes.append(product_dict.index(scene))
                os.remove(scene_obj)
            else:
                overlap_area=shapefile_area(prods_bbox)
                # Kick out scenes below specified overlap threshold
                if overlap_area < minimumOverlap:
                    log.debug(f'Rejected scene {scene_obj} has only '
                              f'{overlap_area}km\u00b2 overlap with bbox')
                    rejected_scenes.append(product_dict.index(scene))
                    os.remove(scene_obj)
                else:
                    save_shapefile(prods_TOTbbox, prods_bbox, 'GeoJSON')
                    # Need to track bounds of max extent
                    # to avoid metadata interpolation issues
                    total_bbox_metadatalyr = total_bbox_metadatalyr.union( \
                        open_shapefile(scene['productBoundingBox'][0], 
                        0, 0))
                    save_shapefile(prods_TOTbbox_metadatalyr,
                                   total_bbox_metadatalyr, 'GeoJSON')

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
        outputBounds = OG_bounds, xRes = arrres[0], yRes = arrres[1], targetAlignedPixels = True, \
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


def prep_metadatalayers(outname, metadata_arr, dem, key, layers, driver):
    """ Wrapper to prep metadata layer for extraction """

    if dem is None:
        raise Exception('No DEM input specified. '
                        'Cannot extract 3D imaging geometry '
                        'layers without DEM to intersect with.')

    ifg = os.path.basename(outname)
    out_dir = os.path.dirname(outname)
    ref_outname = deepcopy(outname)
    # ionosphere layer, heights do not exist to exit
    if metadata_arr[0].split('/')[-1] == 'ionosphere':
        gdal.BuildVRT(outname +'.vrt', metadata_arr)
        return [0], None, outname

    # capture model if tropo product
    model_name = None
    if 'tropo' in key:
        model_name = metadata_arr[0].split('/')[-3]
        out_dir = os.path.join(out_dir, model_name)
        outname = os.path.join(out_dir, ifg)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    # Get height values
    zdim = gdal.Open(metadata_arr[0]).GetMetadataItem( \
                              'NETCDF_DIM_EXTRA')[1:-1]
    hgt_field = f'NETCDF_DIM_{zdim}_VALUES'
    # Check if height layers are consistent
    if not os.path.exists(outname +'.vrt') and \
         not len(set([gdal.Open(i).GetMetadataItem(hgt_field) \
         for i in metadata_arr]))==1:
        raise Exception('Inconsistent heights for '
             'metadata layer(s) ', metadata_arr, \
             ' corresponding heights: ', \
             [gdal.Open(i).GetMetadataItem( \
             hgt_field) \
             for i in metadata_arr])

    if 'tropo' in key or key == 'solidEarthTide':
        # get ref and sec paths
        date_dir = os.path.join(out_dir, 'dates')
        if not os.path.exists(date_dir):
            os.mkdir(date_dir)
        ref_outname = os.path.join(date_dir, ifg.split('_')[0])
        sec_outname = os.path.join(date_dir, ifg.split('_')[1])
        ref_str = 'reference/' + key
        sec_str = 'secondary/' + key
        sec_metadata_arr = [i[:-len(ref_str)] + sec_str for i in metadata_arr]
        tup_outputs = [
            (ref_outname, metadata_arr),
            (sec_outname, sec_metadata_arr)
        ]

        # write ref and sec files
        for i in tup_outputs:
            # delete temporary files to circumvent potential inconsistent dims
            for j in glob.glob(i[0]+'*'): os.remove(j)
            # handle expected offset between frames only for SET
            if key == 'solidEarthTide':
                product_stitch_sequential_metadata(i[1],
                                                   output_meta=i[0],
                                                   output_format=driver,
                                                   verbose=True)
            else:
                gdal.BuildVRT(i[0]+'.vrt', i[1])

            # write height layers
            gdal.Open(i[0]+'.vrt').SetMetadataItem(hgt_field, \
                 gdal.Open(i[1][0]).GetMetadataItem(hgt_field))

        # compute differential
        generate_diff(ref_outname, sec_outname, outname, key, key, False,
                      hgt_field, driver)

        # write raster to file if it does not exist
        if key in layers:
            for i in [ref_outname, sec_outname]:
                if not os.path.exists(i):
                    # write to file
                    da = open_rasterio(i + '.vrt')
                    da.rio.to_raster(i, driver=driver)
                    da.close()
    else:
        if not os.path.exists(outname+'.vrt'):
            gdal.BuildVRT(outname+'.vrt', metadata_arr)
            # write height layers
            gdal.Open(outname+'.vrt').SetMetadataItem(hgt_field, \
                 gdal.Open(metadata_arr[0]).GetMetadataItem(hgt_field))

    return hgt_field, model_name, ref_outname


def generate_diff(ref_outname, sec_outname, outname, key, OG_key, tropo_total,
                  hgt_field, driver):
    """ Compute differential from reference and secondary scenes """

    # if specified workdir doesn't exist, create it
    output_dir = os.path.dirname(outname)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # open intermediate files
    with open_rasterio(sec_outname + '.vrt') as da_sec:
        arr_sec = da_sec.data
    with open_rasterio(ref_outname + '.vrt') as da_ref:
        arr_ref = da_ref.data

    # make the total arr
    # if computing total delay, must add dry and wet
    if tropo_total:
        arr_total = arr_sec + arr_ref
    else:
        arr_total = arr_sec - arr_ref
    da_total = da_sec.copy()
    da_total.data = arr_total
    
    # update attributes
    da_total.name = key
    og_da_attrs = da_total.attrs
    da_attrs = {}
    for key in og_da_attrs:
        new_key = key.replace(OG_key, key)
        new_val = og_da_attrs[key]
        if isinstance(new_val, str):
            new_val = new_val.replace(OG_key, key)
        da_attrs[new_key] = new_val
    da_total = da_total.assign_attrs(da_attrs)

    # write initial array to file
    da_total.rio.to_raster(outname, driver=driver)
    ds = gdal.BuildVRT(f'{outname}.vrt', outname)
    # fix if numpy array not set properly
    if not isinstance(da_attrs[hgt_field], np.ndarray):
        da_attrs[hgt_field] = np.array(da_attrs[hgt_field])
    da_attrs[hgt_field] = da_attrs[hgt_field].tolist()
    ds.SetMetadata(da_attrs)
    del ds

    return


def handle_epoch_layers(layers,
            product_dict,
            lyr_path,
            key,
            sec_key,
            ref_key,
            tropo_total,
            workdir,
            prog_bar,
            bounds,
            dem_bounds,
            prods_TOTbbox,
            dem,
            lat,
            lon,
            mask,
            outputFormat,
            verbose,
            multilooking,
            rankedResampling,
            num_threads):
    """Manage reference/secondary components for correction layers.
    Specifically record reference/secondary components within a `dates` subdir
    and deposit the differential fields in the level above.
    """

    # Depending on type, set sec/ref output dirs
    if key == 'troposphereTotal':
        layers.append(key)
        sec_workdir = os.path.join(os.path.dirname(workdir),
                                   sec_key)
        ref_workdir = os.path.join(os.path.dirname(workdir),
                                   ref_key)
    else:
        sec_workdir = deepcopy(workdir)
        ref_workdir = deepcopy(workdir)
    # If specified workdir doesn't exist, create it
    all_workdirs = [workdir, sec_workdir, ref_workdir]
    all_workdirs = list(set(all_workdirs))
    existing_outputs = []
    for i in all_workdirs:
        if not os.path.exists(i):
            os.mkdir(i)
        # for dedup, record previous outputs to avoid reprocessing
        existing_outputs.extend(glob.glob(os.path.join( \
                     i, '*/*[0-9].vrt')))
        existing_outputs.extend(glob.glob(os.path.join( \
                     i, '*/dates/*[0-9].vrt')))
        existing_outputs.extend(glob.glob(os.path.join( \
                     i, '*[0-9].vrt')))
        existing_outputs.extend(glob.glob(os.path.join( \
                     i, 'dates/*[0-9].vrt')))
    existing_outputs = list(set(existing_outputs))

    # set driver
    driver = 'ENVI' if outputFormat == 'VRT' else outputFormat

    # Iterate through all IFGs
    all_outputs = []
    for i in enumerate(product_dict[0]):
        ifg         = product_dict[1][i[0]][0]
        outname = os.path.abspath(os.path.join(workdir, ifg))

        # create temp files for ref/sec components
        ref_outname = os.path.abspath(os.path.join(ref_workdir, ifg))
        hgt_field, model_name, ref_outname = prep_metadatalayers(ref_outname,
                                      i[1], dem,
                                      ref_key, layers, driver)
        # Update progress bar
        prog_bar.update(i[0]+1,suffix=ifg)

        # record output directories
        if model_name is not None:
            all_outputs.append(os.path.join(workdir, model_name))
            all_outputs.append(os.path.join(ref_workdir, model_name))
            all_outputs.append(os.path.join(sec_workdir, model_name))
        else:
            all_outputs.append(workdir)
            all_outputs.append(ref_workdir)
            all_outputs.append(sec_workdir)

        # capture if tropo and separate distinct wet and hydro layers
        if 'tropo' in key:
            sec_outname = os.path.abspath(os.path.join(sec_workdir, ifg))
            wet_path = os.path.join(lyr_path, model_name,
                                    'reference', ref_key)
            dry_path = os.path.join(lyr_path, model_name,
                                    'reference', sec_key)
            sec_comp = [j.replace(wet_path, dry_path) \
                                for j in i[1]]
            hgt_field, model_name, sec_outname = prep_metadatalayers(
                                               sec_outname,
                                               sec_comp, dem,
                                               sec_key, layers, driver)
            # if specified, compute total delay
            if tropo_total:
                model_dir = os.path.abspath(os.path.join(workdir, model_name))
                outname = os.path.join(model_dir, ifg)
                # compute reference diff
                ref_diff = ref_outname
                sec_diff = sec_outname
                outname_diff = os.path.join(model_dir, 'dates',
                    os.path.basename(ref_diff))
                print('ref_diff, sec_diff, outname_diff', ref_diff, sec_diff, outname_diff)
                if not os.path.exists(outname_diff):
                    generate_diff(ref_diff, sec_diff, outname_diff, key,
                      sec_key, tropo_total, hgt_field, driver)
                # compute secondary diff
                ref_diff = os.path.join(os.path.dirname(ref_outname),
                                        ifg.split('_')[1])
                sec_diff = os.path.join(os.path.dirname(sec_outname),
                                        ifg.split('_')[1])
                outname_diff = os.path.join(model_dir, 'dates',
                    os.path.basename(ref_diff))
                if not os.path.exists(outname_diff):
                    generate_diff(ref_diff, sec_diff, outname_diff, key,
                      sec_key, tropo_total, hgt_field, driver)
                # compute total diff
                ref_diff = os.path.join(ref_workdir, model_name, ifg)
                sec_diff = os.path.join(sec_workdir, model_name, ifg)
                generate_diff(ref_diff, sec_diff, outname, key,
                      sec_key, tropo_total, hgt_field, driver)
        else:
            sec_outname = os.path.dirname(ref_outname)
            sec_outname = os.path.abspath(os.path.join(sec_outname,
                                      ifg.split('_')[0]))

    # delete temporary files if layers not requested
    prod_ver_list = i[1]
    for i in all_workdirs:
        key_name = os.path.basename(i)
        if os.path.exists(i):
            if key_name not in layers or len(os.listdir(i)) == 0:
                shutil.rmtree(i)

    # interpolate and intersect epochs for user requested layers
    all_outputs = list(set(all_outputs))
    for i in enumerate(all_outputs):
        if os.path.exists(i[1]):
            # Update progress bar
            prog_bar.update(i[0]+1,suffix=ifg)

            record_epochs = []
            record_epochs.extend(glob.glob(os.path.join( \
                     i[1], '*[0-9].vrt')))
            record_epochs.extend(glob.glob(os.path.join( \
                     i[1], 'dates/*[0-9].vrt')))
            # dedup check for interpolating only new files
            record_epochs = [j for j in record_epochs \
                                 if j not in existing_outputs]
            for j in enumerate(record_epochs):
                # Interpolate/intersect with DEM before cropping
                finalize_metadata(j[1][:-4], bounds, dem_bounds,
                      prods_TOTbbox, dem, lat, lon, hgt_field,
                      prod_ver_list, mask, outputFormat,
                      verbose=verbose)
                # If necessary, resample raster
                if multilooking is not None:
                    resampleRaster(j[1][:-4], multilooking, bounds,
                          prods_TOTbbox, rankedResampling,
                          outputFormat=outputFormat, num_threads=num_threads)

                # Track consistency of dimensions
                if j[0] == 0:
                    ref_wid, ref_hgt,_,_,_ = get_basic_attrs(j[1][:-4])
                else:
                    prod_wid, prod_hgt,_,_,_ = get_basic_attrs(j[1][:-4])
                    if (ref_wid != prod_wid) or (ref_hgt != prod_hgt):
                        raise Exception(f'Inconsistent product dims '
                              'between products {outname} and {prev_outname}:'
                              'check respective width ({ref_wid}, {ref_hgt})'
                              'and height ({prod_wid}, {prod_hgt})')
                prev_outname = j[1][:-4]
        
    prog_bar.close()

    return

def export_products(full_product_dict, bbox_file, prods_TOTbbox, layers,
                    rankedResampling=False, dem=None, lat=None, lon=None,
                    mask=None, outDir='./', outputFormat='VRT',
                    stitchMethodType='overlap', verbose=None, num_threads='2',
                    multilooking=None, tropo_total=False, model_names=[]):
    """Export layer and 2D meta-data layers (at the product resolution).
    The function finalize_metadata is called to derive the 2D metadata layer.
    Dem/lat/lon arrays must be passed for this process.
    The keys specify which layer to extract from the dictionary.
    All products are cropped by the bounds from the input bbox_file,
    and clipped to the track extent denoted by the input prods_TOTbbox.
    Optionally, a user may pass a mask-file.
    """
    ##Progress bar
    from ARIAtools import progBar

    if not layers and not tropo_total: return # only bbox

    # create dictionary of all inputs needed for correction lyr extraction
    lyr_input_dict = {
            'layers': layers,
            'prods_TOTbbox': prods_TOTbbox,
            'dem': dem,
            'lat': lat,
            'lon': lon,
            'mask': mask,
            'verbose': verbose,
            'multilooking': multilooking,
            'rankedResampling': rankedResampling,
            'num_threads': num_threads
    }

    # get bounds
    bounds = open_shapefile(bbox_file, 0, 0).bounds
    lyr_input_dict['bounds'] = bounds
    if dem is not None:
        dem_bounds = [dem.GetGeoTransform()[0],dem.GetGeoTransform()[3]+ \
        (dem.GetGeoTransform()[-1]*dem.RasterYSize),dem.GetGeoTransform()[0]+ \
        (dem.GetGeoTransform()[1]*dem.RasterXSize),dem.GetGeoTransform()[3]]
        lyr_input_dict['dem_bounds'] = dem_bounds

    # Mask specified, so file must be physically extracted,
    # cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT' and mask is not None:
        outputFormat='ENVI'
    lyr_input_dict['outputFormat'] = outputFormat

    # If specified, extract tropo layers
    tropo_lyrs = ['troposphereWet', 'troposphereHydrostatic']
    if tropo_total or list(set.intersection(*map(set, \
                        [layers, tropo_lyrs]))) != []:
        # set input keys
        lyr_prefix = '/science/grids/corrections/external/troposphere/'
        key = 'troposphereTotal'
        wet_key = 'troposphereWet'
        dry_key   = 'troposphereHydrostatic'
        workdir = os.path.join(outDir, key)
        lyr_input_dict['lyr_path'] = lyr_prefix
        lyr_input_dict['key'] = key
        lyr_input_dict['sec_key'] = dry_key
        lyr_input_dict['ref_key'] = wet_key
        lyr_input_dict['tropo_total'] = tropo_total
        lyr_input_dict['workdir'] = workdir

        # loop through valid models
        for i in model_names:
            model = wet_key + f'_{i}'
            tropo_lyrs.append(model)
            tropo_lyrs.append(dry_key + f'_{i}')
            product_dict = [
              [j[model] for j in full_product_dict if model in j.keys()], \
              [j["pair_name"] for j in full_product_dict if model in j.keys()]
            ]

            # set iterative keys
            prog_bar = progBar.progressBar(maxValue=len(product_dict[0]),
                                       prefix=f'Generating: {model} {key} - ')
            lyr_input_dict['prog_bar'] = prog_bar
            lyr_input_dict['product_dict'] = product_dict

            # extract layers
            handle_epoch_layers(**lyr_input_dict)

    # If specified, extract solid earth tides
    tropo_lyrs = list(set(tropo_lyrs))
    ext_corr_lyrs = tropo_lyrs + ['solidEarthTide', 'troposphereTotal']
    if list(set.intersection(*map(set, \
                        [layers, ['solidEarthTide']]))) != []:
        lyr_prefix = '/science/grids/corrections/external/tides/solidEarth/'
        key = 'solidEarthTide'
        ref_key = key
        sec_key = key
        product_dict = [[j[key] for j in full_product_dict], \
                      [j["pair_name"] for j in full_product_dict]]

        workdir = os.path.join(outDir, key)
        prog_bar = progBar.progressBar(maxValue=len(product_dict[0]),
                                       prefix='Generating: '+key+' - ')

        # set input keys
        lyr_input_dict['prog_bar'] = prog_bar
        lyr_input_dict['product_dict'] = product_dict
        lyr_input_dict['lyr_path'] = lyr_prefix
        lyr_input_dict['key'] = key
        lyr_input_dict['sec_key'] = sec_key
        lyr_input_dict['ref_key'] = ref_key
        lyr_input_dict['tropo_total'] = False
        lyr_input_dict['workdir'] = workdir

        # extract layers
        handle_epoch_layers(**lyr_input_dict)

    # Loop through other user expected layers
    layers = [i for i in layers if i not in ext_corr_lyrs]
    for key in layers:
        product_dict=[[j[key] for j in full_product_dict], [j["pair_name"] for j in full_product_dict]]
        workdir=os.path.join(outDir,key)

        prog_bar = progBar.progressBar(maxValue=len(product_dict[0]),prefix='Generating: '+key+' - ')

        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        # Iterate through all IFGs
        for i in enumerate(product_dict[0]):
            ifg = product_dict[1][i[0]][0]
            outname=os.path.abspath(os.path.join(workdir, ifg))
            ##Update progress bar
            prog_bar.update(i[0]+1,suffix=ifg)

            # Extract/crop metadata layers
            if any(":/science/grids/imagingGeometry" \
                 in s for s in i[1]) or \
                 any(":/science/grids/corrections" \
                 in s for s in i[1]):
                # make VRT pointing to metadata layers in standard product
                hgt_field, model_name, outname = prep_metadatalayers(outname,
                                                     i[1], dem, key, layers,
                                                     outputFormat)

                # Interpolate/intersect with DEM before cropping
                finalize_metadata(outname, bounds, dem_bounds,
                                  prods_TOTbbox, dem, lat, lon, hgt_field,
                                  i[1], mask, outputFormat,
                                  verbose=verbose)

            # Extract/crop full res layers, except for "unw" and "conn_comp" which requires advanced stitching
            elif key!='unwrappedPhase' and \
                 key!='connectedComponents':
                if outputFormat=='VRT' and mask is None:
                    # building the virtual vrt
                    gdal.BuildVRT(outname+ "_uncropped" +'.vrt', i[1])
                    # building the cropped vrt
                    gdal.Warp(outname+'.vrt', outname+"_uncropped"+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, options=['NUM_THREADS=%s'%(num_threads)]))
                else:
                    # building the VRT
                    gdal.BuildVRT(outname +'.vrt', i[1])
                    gdal.Warp(outname, outname+'.vrt', options=gdal.WarpOptions(format=outputFormat, cutlineDSName=prods_TOTbbox, outputBounds=bounds, options=['NUM_THREADS=%s'%(num_threads)]))

                    # Update VRT
                    gdal.Translate(outname+'.vrt', outname, options=gdal.TranslateOptions(format="VRT"))

                    # Apply mask (if specified).
                    if mask is not None:
                        update_file=gdal.Open(outname,gdal.GA_Update)
                        update_file.GetRasterBand(1).WriteArray(mask.ReadAsArray()*gdal.Open(outname+'.vrt').ReadAsArray())
                        del update_file

            # Extract/crop phs and conn_comp layers
            else:
                # get connected component input files
                conn_files = full_product_dict[i[0]]['connectedComponents']
                prod_bbox_files = \
                    full_product_dict[i[0]]['productBoundingBoxFrames']
                outFileConnComp = \
                    os.path.join(outDir, 'connectedComponents', ifg)

                # Check if phs phase and conn_comp files are already generated
                outFilePhs = os.path.join(outDir, 'unwrappedPhase', ifg)
                if not os.path.exists(outFilePhs) or not \
                     os.path.exists(outFileConnComp):
                    phs_files = full_product_dict[i[0]]['unwrappedPhase']
                    # calling the stitching methods
                    if stitchMethodType == 'overlap':
                        product_stitch_overlap(phs_files, conn_files,
                                              prod_bbox_files, bounds,
                                              prods_TOTbbox,
                                              outFileUnw=outFilePhs,
                                              outFileConnComp=outFileConnComp,
                                              mask=mask,
                                              outputFormat=outputFormat,
                                              verbose=verbose)
                        
                    elif stitchMethodType == '2stage':
                        product_stitch_2stage(phs_files,
                                              conn_files,
                                              bounds,
                                              prods_TOTbbox,
                                              outFileUnw=outFilePhs,
                                              outFileConnComp=outFileConnComp,
                                              mask=mask,
                                              outputFormat=outputFormat,
                                              verbose=verbose)

                    elif stitchMethodType == 'sequential':
                        product_stitch_sequential(phs_files,
                                                 conn_files,
                                                 bounds=bounds,
                                                 clip_json=prods_TOTbbox,
                                                 output_unw=outFilePhs,
                                                 output_conn=outFileConnComp,
                                                 mask_file=mask, # str filename
                                                 output_format=outputFormat,
                                                 range_correction=True,
                                                 save_fig=False,
                                                 overwrite=True,
                                                 verbose=verbose)

                    # If necessary, resample phs/conn_comp file
                    if multilooking is not None:
                        resampleRaster(outFilePhs, multilooking, bounds,
                                       prods_TOTbbox, rankedResampling,
                                       outputFormat=outputFormat,
                                       num_threads=num_threads)

            # If necessary, resample raster
            if multilooking is not None and \
                 key!='unwrappedPhase' and \
                 key!='connectedComponents':
                resampleRaster(outname, multilooking, bounds, prods_TOTbbox,
                               rankedResampling, outputFormat=outputFormat,
                               num_threads=num_threads)

            # Track consistency of dimensions
            if i[0] == 0:
                ref_wid, ref_height, _, _, _ = get_basic_attrs(outname +
                                                               '.vrt')
            else:
                prod_wid, prod_height, _, _, _ = get_basic_attrs(outname +
                                                                 '.vrt')
                if (ref_wid != prod_wid) or (ref_height != prod_height):
                     raise Exception(f'Inconsistent product dims between'
                         'products {outname} and {prev_outname}:'
                         'check respective width ({ref_wid}, {ref_height})'
                         'and height ({prod_wid}, {prod_height})')
            prev_outname = os.path.abspath(os.path.join(workdir, ifg))

        prog_bar.close()

        # check directory for quality control plots
        plots_subdir = os.path.abspath(os.path.join(outDir,
                                       'metadatalyr_plots'))
        # delete directory if empty
        if os.path.exists(plots_subdir):
            if len(os.listdir(plots_subdir)) == 0:
                shutil.rmtree(plots_subdir)

    return


def finalize_metadata(outname, bbox_bounds, dem_bounds, prods_TOTbbox, dem, \
                      lat, lon, hgt_field, prod_list, mask=None, \
                      outputFormat='ENVI', verbose=None, num_threads='2'):
    """Interpolate and extract 2D metadata layer.
    2D metadata layer is derived by interpolating and then intersecting
    3D layers with a DEM.
    Lat/lon arrays must also be passed for this process.
    """
    # import dependencies
    from scipy.interpolate import RegularGridInterpolator

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    if outputFormat=='VRT':
        outputFormat='ENVI'

    # get final shape
    arrshape=gdal.Open(dem.GetDescription()).ReadAsArray().shape
    # load layered metadata array
    tmp_name = outname+'.vrt'
    data_array = gdal.Warp('', tmp_name,
                     options=gdal.WarpOptions(format="MEM",
                         options=['NUM_THREADS=%s'%(num_threads)]))

    # get minimum version
    version_check = []
    for i in prod_list:
        v_num = i.split(':')[-2].split('/')[-1]
        v_num = v_num.split('.nc')[0][-5:]
        version_check.append(v_num)
    version_check = min(version_check)

    # metadata layer quality check, correction applied if necessary
    # only apply to geometry layers and prods derived from older ISCE versions
    geom_lyrs = ['bPerpendicular', 'bParallel', 'incidenceAngle',
                 'lookAngle', 'azimuthAngle']
    metadatalyr_name = outname.split('/')[-2]
    if metadatalyr_name in geom_lyrs and version_check < '2_0_4':
        # create directory for quality control plots
        plots_subdir = os.path.abspath(os.path.join(outname, '../..',
                                       'metadatalyr_plots', metadatalyr_name))
        if not os.path.exists(plots_subdir):
            os.makedirs(plots_subdir)
        data_array = metadata_qualitycheck( \
            data_array, \
            os.path.basename(os.path.dirname(outname)), \
            outname, \
            verbose).data_array
        # delete directory if empty
        if len(os.listdir(plots_subdir)) == 0:
            shutil.rmtree(plots_subdir)

    # only perform DEM intersection for rasters with valid height levels
    nohgt_lyrs = ['ionosphere']
    if metadatalyr_name not in nohgt_lyrs:
        tmp_name = outname+'_temp'
        # Define lat/lon/height arrays for metadata layers
        heightsMeta = np.array(gdal.Open(outname+'.vrt').GetMetadataItem( \
             hgt_field)[1:-1].split(','), dtype='float32')

        latitudeMeta = np.linspace(data_array.GetGeoTransform()[3],
                           data_array.GetGeoTransform()[3] + \
                           (data_array.GetGeoTransform()[5] * \
                           data_array.RasterYSize), data_array.RasterYSize)
        longitudeMeta = np.linspace(data_array.GetGeoTransform()[0],
                           data_array.GetGeoTransform()[0] + \
                           (data_array.GetGeoTransform()[1] * \
                           data_array.RasterXSize), data_array.RasterXSize)

        da_dem = open_rasterio(dem.GetDescription(), 
                     band_as_variable=True)['band_1']

        # interpolate the DEM to the GUNW lat/lon
        da_dem1 = da_dem.interp(x=lon[0, :], 
            y=lat[:, 0]).fillna(dem.GetRasterBand(1).GetNoDataValue())

        # hack to get an stack of coordinates for the interpolator
        # to interpolate in the right shape
        pnts = transformPoints(lat, lon, da_dem1.data, 'EPSG:4326', 'EPSG:4326')

        # set up the interpolator with the GUNW cube
        interper = RegularGridInterpolator((latitudeMeta, longitudeMeta,
                   heightsMeta), data_array.ReadAsArray().transpose(1, 2, 0),
                   fill_value=np.nan, bounds_error=False)

        # interpolate cube to DEM points
        out_interpolated = interper(pnts.transpose(2, 1, 0))

        # Save file
        renderVRT(tmp_name, out_interpolated,
                  geotrans=dem.GetGeoTransform(),
                  drivername=outputFormat,
                  gdal_fmt=data_array.ReadAsArray().dtype.name,
                  proj=dem.GetProjection(),
                  nodata=data_array.GetRasterBand(1).GetNoDataValue())
        del out_interpolated

    # Since metadata layer extends at least one grid node
    # outside of the expected track bounds,
    # it must be cut to conform with these bounds.
    # Crop to track extents
    gdal.Warp(outname,
              tmp_name,
              options=gdal.WarpOptions(format=outputFormat,
                  cutlineDSName=prods_TOTbbox,
                  outputBounds=bbox_bounds,
                  dstNodata=data_array.GetRasterBand(1).GetNoDataValue(),
                  width=arrshape[1],
                  height=arrshape[0],
                  options=['NUM_THREADS=%s'%(num_threads)+' -overwrite']))
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
        del update_file, out_interpolated

    del data_array


def gacos_correction(full_product_dict, gacos_products, bbox_file,
                     prods_TOTbbox, outDir='./', outputFormat='VRT',
                     verbose=None, num_threads='2'):
    """Perform tropospheric corrections.
    Must provide valid path to GACOS products.
    All products are cropped by the bounds from the input bbox_file,
    and clipped to the track extent denoted by the input prods_TOTbbox.
    """
    # Import functions
    import tarfile
    from datetime import datetime
    from ARIAtools.vrtmanager import rscGacos, tifGacos
    from ARIAtools.shapefile_util import save_shapefile
    from shapely.geometry import Polygon

    # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
    outputFormat = 'ENVI' if outputFormat == 'VRT' else outputFormat

    user_bbox    = open_shapefile(bbox_file, 0, 0)
    bounds       = user_bbox.bounds

    product_dict = [[j['unwrappedPhase'] for j in full_product_dict[1]],
                  [j['lookAngle'] for j in full_product_dict[1]],
                  [j["pair_name"] for j in full_product_dict[1]]]
    metadata_dict = [[j['azimuthZeroDopplerMidTime'] for j in \
                    full_product_dict[0]],
                   [j['wavelength'] for j in full_product_dict[0]]]
    workdir = os.path.join(outDir,'gacos_corrections')

    # If specified workdir doesn't exist, create it
    os.makedirs(workdir, exist_ok=True)

    # Get list of all dates for which standard products exist
    date_list = []

    for i in product_dict[2]:
        date_list.append(i[0][:8]); date_list.append(i[0][9:])
    date_list = list(set(date_list)) # trim to unique dates only

    ### Determine if file input is single file, a list, or wildcard.
    # If list of files
    if len([str(val) for val in gacos_products.split(',')]) > 1:
        gacos_products = [str(i) for i in gacos_products.split(',')]
        # Sort and parse tropo products
        tropo_sublist = []
        for i in gacos_products:
            # If wildcard
            if '*' in i:
                tropo_sublist.append( \
                     glob.glob(os.path.expanduser(os.path.expandvars(i))))
            else:
                tropo_sublist.append([i])
        # Finalize list of tropo products
        gacos_products = []
        for sublist in tropo_sublist:
            for item in sublist:
                gacos_products.append(os.path.abspath(item))

    # If single file or wildcard
    else:
        # If single file
        if os.path.isfile(gacos_products):
            gacos_products = [gacos_products]
        # If wildcard
        else:
            gacos_products = glob.glob(os.path.expanduser( \
                                       os.path.expandvars(gacos_products)))
        # Convert relative paths to absolute paths
        gacos_products = [os.path.abspath(i) for i in gacos_products]
    if len(gacos_products) == 0:
        raise Exception('No file match found')

    ###Extract tarfiles
    # Setup dictionary to track for products that are to be merged
    tropo_date_dict = {}
    for i in date_list:
        tropo_date_dict[i] = []
        tropo_date_dict[i+"_UTC"] = []
    for i in enumerate(gacos_products):
        if not os.path.isdir(i[1]) and not i[1].endswith('.ztd') \
            and not i[1].endswith('.tif'):
            untar_dir = os.path.join(os.path.abspath(os.path.join(i[1], \
                         os.pardir)),
                        os.path.basename(i[1]).split('.')[0] + '_extracted')
            if not tarfile.is_tarfile(i[1]):
                raise Exception('Cannot extract %s because it is not a ' \
                                'valid tarfile. Resolve this ' \
                                'and relaunch'%(i[1]))
            log.info('Extracting GACOS tarfile %s to %s.',
                os.path.basename(i[1]), untar_dir)
            tarfile.open(i[1]).extractall(path = untar_dir)
            gacos_products[i[0]] = untar_dir
        # Loop through each GACOS product file, differentiating between direct input list of GACOS products vs parent directory
        if i[1].endswith('.ztd') or i[1].endswith('.tif'):
            ztd_list = [i[1]]
        else:
            ztd_list = glob.glob(os.path.join(gacos_products[i[0]], '*.ztd')) \
                + glob.glob(os.path.join(gacos_products[i[0]], '*.tif'))
            # prioritize older .ztd over .tif duplicates if the former exists
            # older .ztd files contain UTC time info
            ztd_list = [i for i in ztd_list if not (i.endswith('.tif') and i.split('.tif')[0] in ztd_list)]
        for k in ztd_list:
            # Only check files corresponding to standard product dates
            ztd_basename = os.path.basename(k)
            ztd_basename = ztd_basename.split('.tif')[0].split('.ztd')[0]
            if ztd_basename in date_list:
                tropo_date_dict[ztd_basename].append(k)
                if os.path.exists(k+'.rsc'):
                    tropo_date_dict[ztd_basename+"_UTC"].append( \
                         ztd_basename[:4] + '-' \
                         + ztd_basename[4:6] \
                         + '-' + ztd_basename[6:8] + '-' \
                         + open(k+'.rsc', 'r').readlines()[-1].split()[1])
                else:
                    tropo_date_dict[ztd_basename+"_UTC"].append( \
                         ztd_basename[:4] + '-' \
                         + ztd_basename[4:6] \
                         + '-' + ztd_basename[6:8] + '-' \
                         + 'NoneUTC')
                # make corresponding VRT file, if it doesn't exist
                if not os.path.exists(k+'.vrt'):
                    # parse metadata from RSC if it exists
                    if os.path.exists(k+'.rsc'):
                        tropo_rsc_dict={}
                        for line in open(k+'.rsc', 'r').readlines():
                            tropo_rsc_dict[line.split()[0]]=line.split()[1]
                        gacos_prod = np.fromfile(k, dtype='float32').reshape( \
                             int(tropo_rsc_dict['FILE_LENGTH']), \
                             int(tropo_rsc_dict['WIDTH']))
                    else:
                        tropo_rsc_dict = tifGacos(k)
                        gacos_prod = gdal.Open(k).ReadAsArray()
                    # Save as GDAL file, using proj from first unwrappedPhase file
                    renderVRT(k, gacos_prod,
                        geotrans=(float(tropo_rsc_dict['X_FIRST']),
                        float(tropo_rsc_dict['X_STEP']),
                        0.0,
                        float(tropo_rsc_dict['Y_FIRST']),
                        0.0,
                        float(tropo_rsc_dict['Y_STEP'])),
                        drivername=outputFormat,
                        gdal_fmt='float32',
                        proj=gdal.Open(os.path.join(outDir,'unwrappedPhase', \
                         product_dict[2][0][0])).GetProjection(),
                        nodata=0.)
                    gacos_prod = None
                    log.debug('GACOS product %s successfully converted to ' \
                              'GDAL-readable raster', k)
                # make corresponding RSC file, if it doesn't exist
                if not os.path.exists(k+'.rsc'):
                    rscGacos(os.path.join(k+'.vrt'), os.path.join(k+'.rsc'),
                             tropo_date_dict)

    # If multiple GACOS directories, merge products.
    gacos_products = list(set([os.path.dirname(i) \
                          if (i.endswith('.ztd') or i.endswith('.tif')) \
                          else i for i in gacos_products]))
    if len(gacos_products)>1:
        gacos_products=os.path.join(outDir,'merged_GACOS')
        log.info('Stitching/storing GACOS products in %s.', gacos_products)
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
                # Create merged rsc file
                rscGacos(outname, merged_rsc, tropo_date_dict)
    else:
        gacos_products = gacos_products[0]

    # Estimate percentage of overlap with tropospheric product
    for i in glob.glob(os.path.join(gacos_products,'*.vrt')):
        # create shapefile
        geotrans = gdal.Open(i).GetGeoTransform()
        bbox = [geotrans[3] \
                 + (gdal.Open(i).ReadAsArray().shape[0]*geotrans[-1]),
                 geotrans[3],
                 geotrans[0],
                 geotrans[0] \
                 + (gdal.Open(i).ReadAsArray().shape[1]*geotrans[1])
        ]
        bbox = Polygon(np.column_stack( \
                       (np.array([bbox[2], bbox[3], bbox[3],bbox[2],bbox[2]]),
                        np.array([bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]))))
        save_shapefile(i+'.json', bbox, 'GeoJSON')
        per_overlap = ((user_bbox.intersection( \
                      open_shapefile(i+'.json', 0, 0)).area) \
                      / (user_bbox.area))*100
        if per_overlap != 100. and per_overlap != 0.:
            log.warning('Common track extent only has %d overlap with' \
                        'tropospheric product %s\n', per_overlap, i[0])
        if per_overlap == 0.:
            raise Exception('No spatial overlap between tropospheric ' \
                            'product %s and defined bounding box. ' \
                            'Resolve conflict and relaunch', i[1])

    # Iterate through all IFGs and apply corrections
    missing_products = []
    for i in range(len(product_dict[0])):
        ifg     = product_dict[2][i][0]
        outname = os.path.join(workdir, ifg)
        gacos_epochs_dir = os.path.join(workdir, 'dates')
        ref_outname = os.path.join(gacos_epochs_dir, f'{ifg[:8]}')
        sec_outname = os.path.join(gacos_epochs_dir, f'{ifg[9:]}')
        outname = os.path.join(workdir, ifg)
        unwname = os.path.join(outDir,'unwrappedPhase', ifg)
        if i == 0:
            meta = gdal.Info(unwname, format='json')
            geoT = meta['geoTransform']
            proj = meta['coordinateSystem']['wkt']
            arrshape = list(reversed(meta['size']))

        tropo_reference = os.path.join(gacos_products, f'{ifg[:8]}.ztd.vrt')
        tropo_secondary = os.path.join(gacos_products, f'{ifg[9:]}.ztd.vrt')

        # if .ztd products don't exist, check if .tif exists
        if not os.path.exists(tropo_reference):
            tropo_reference = os.path.join(gacos_products, f'{ifg[:8]}.ztd.tif.vrt')
        if not os.path.exists(tropo_secondary):
            tropo_secondary = os.path.join(gacos_products, f'{ifg[9:]}.ztd.tif.vrt')


        # skip if corrected already generated and does not need to be updated
        if os.path.exists(outname):
            # get unwrappedPhase geotrans and productbounding box
            unw_prodcheck = gdal.Open(unwname)
            unw_geotrans = unw_prodcheck.GetGeoTransform()
            unw_prodcheck = np.isfinite(unw_prodcheck.ReadAsArray())
            tropo_prodcheck = gdal.Open(outname)
            output_geotrans = tropo_prodcheck.GetGeoTransform()
            tropo_prodcheck = np.isfinite(tropo_prodcheck.ReadAsArray())
            if unw_geotrans == output_geotrans and np.array_equal( \
                               unw_prodcheck, tropo_prodcheck):
                continue
            del unw_prodcheck, tropo_prodcheck

        if os.path.exists(tropo_reference) and os.path.exists(tropo_secondary):
            # Check if tropo products are temporally consistent with IFG
            for j in [tropo_reference, tropo_secondary]:
                # Get ARIA product times
                aria_rsc_dict = {}
                aria_rsc_dict['azimuthZeroDopplerMidTime'] = \
                    [datetime.strptime(os.path.basename(j)[:4]
                    + '-' + os.path.basename(j)[4:6]
                    + '-' + os.path.basename(j)[6:8]
                    + '-' + m[11:], "%Y-%m-%d-%H:%M:%S.%f") \
                     for m in metadata_dict[0][0]]
                # Get tropo product UTC times
                tropo_rsc_dict = {}
                tropo_rsc_dict['TIME_OF_DAY'] = open(j[:-4]+'.rsc', \
                     'r').readlines()[-1].split()[1].split('UTC')[:-1]
                # If new TIF product, UTC times not available
                if 'None' in tropo_rsc_dict['TIME_OF_DAY'][0]:
                    tropo_rsc_dict['TIME_OF_DAY'] = [ \
                         max(aria_rsc_dict['azimuthZeroDopplerMidTime'])]
                # If stitched tropo product, must account for date change (if applicable)
                elif '-' in tropo_rsc_dict['TIME_OF_DAY'][0]:
                    tropo_rsc_dict['TIME_OF_DAY'] = [datetime.strptime(m[:10]
                        + '-' + m[11:].split('.')[0] + '-'
                        + str(round(float('0.' + m[11:].split('.')[1])*60)),
                        "%Y-%m-%d-%H-%M") \
                         for m in tropo_rsc_dict['TIME_OF_DAY']
                    ]
                else:
                    tropo_rsc_dict['TIME_OF_DAY'] = [datetime.strptime( \
                         os.path.basename(j)[:4] + '-'
                        + os.path.basename(j)[4:6] + '-'
                        + os.path.basename(j)[6:8] + '-'
                        + tropo_rsc_dict['TIME_OF_DAY'][0].split('.')[0]
                        + '-' + str(round(float('0.'
                        + tropo_rsc_dict['TIME_OF_DAY'][0].split('.')[-1]) \
                         *60)), "%Y-%m-%d-%H-%M")
                    ]

                # Check and report if tropospheric product falls outside of standard product range
                latest_start = max(aria_rsc_dict['azimuthZeroDopplerMidTime']
                                   + [min(tropo_rsc_dict['TIME_OF_DAY'])])
                earliest_end = min(aria_rsc_dict['azimuthZeroDopplerMidTime']
                                   + [max(tropo_rsc_dict['TIME_OF_DAY'])])
                delta = (earliest_end - latest_start).total_seconds() + 1
                if delta<0:
                    log.warning('tropospheric product was generated %f ' \
                                'secs outside of acquisition interval for ' \
                                'scene %s in IFG %s',
                                abs(delta), os.path.basename(j)[:8],
                                product_dict[2][i][0])

            # Open corresponding tropo products and pass the difference
            tropo_reference = gdal.Warp('', tropo_reference, format="MEM",
                                      outputBounds=bounds, width=arrshape[1],
                                      height=arrshape[0]).ReadAsArray()
            tropo_secondary = gdal.Warp('', tropo_secondary, format="MEM",
                                        outputBounds=bounds, width=arrshape[1],
                                        height=arrshape[0]).ReadAsArray()
            tropo_product  = np.subtract(tropo_secondary, tropo_reference)

            # Convert troposphere from m to rad
            scale         = float(metadata_dict[1][i][0]) / (4*np.pi)
            tropo_product /= scale

            # Account for incAngle
            # if in TS mode, only 1 incfile would be generated, so check for this
            path_inc = os.path.join(outDir, 'incidenceAngle', ifg)
            if os.path.exists(path_inc):
                da = rio.open(path_inc)
            else:
                da = rio.open(path_inc.replace(ifg, product_dict[2][0][0]))
            inc_arr = da.read().squeeze()
            inc_arr = np.where(np.isclose(inc_arr, da.nodata), np.nan, inc_arr)
            cos_inc = np.cos(np.deg2rad(inc_arr))

            tropo_product /= cos_inc

            # Save differential field to file
            tropo_product = np.where(np.isnan(tropo_product), 0., tropo_product)
            renderVRT(outname, tropo_product,
                      geotrans=geoT, drivername=outputFormat,
                      gdal_fmt='float32', proj=proj, nodata=0.)

            # check if reference and secondary scenes are written to file
            if not os.path.exists(ref_outname):
                tropo_reference /= scale
                tropo_reference /= cos_inc
                tropo_reference = np.where(np.isnan(tropo_reference), 0., tropo_reference)
                renderVRT(ref_outname, tropo_reference,
                      geotrans=geoT, drivername=outputFormat,
                      gdal_fmt='float32', proj=proj, nodata=0.)
            if not os.path.exists(sec_outname):
                tropo_secondary /= scale
                tropo_secondary /= cos_inc
                tropo_secondary = np.where(np.isnan(tropo_secondary), 0., tropo_secondary)
                renderVRT(sec_outname, tropo_secondary,
                      geotrans=geoT, drivername=outputFormat,
                      gdal_fmt='float32', proj=proj, nodata=0.)


            del tropo_product, tropo_reference, \
                tropo_secondary, da, inc_arr, cos_inc

            # Track consistency of dimensions
            if i[0] == 0:
                ref_wid, ref_height, _ = get_basic_attrs(outname)
            else:
                prod_wid, prod_height, _ = get_basic_attrs(outname)
                if (ref_wid != prod_wid) or (ref_height != prod_height):
                     raise Exception(f'Inconsistent product dims between'
                         'products {outname} and {prev_outname}:'
                         'check respective width ({ref_wid}, {ref_height})'
                         'and height ({prod_wid}, {prod_height})')
            prev_outname = os.path.join(workdir, ifg)

        else:
            log.warning('Must skip IFG %s, because the tropospheric ' \
                        'products corresponding to the reference and/or ' \
                        'secondary products are not found in the ' \
                        'specified folder %s',
                        ifg, gacos_products)
            for j in [tropo_reference, tropo_secondary]:
                if not os.path.exists(j) and j not in missing_products:
                    missing_products.append(j)
    # Print list of dates missing tropospheric corrections
    if len(missing_products) > 0:
        missing_products = [os.path.basename(i)[:8] for i in missing_products]
        log.debug("Tropo products for the following dates are missing:")
        log.debug(missing_products)


def main(inps=None):
    """Main workflow for extracting layers from ARIA products."""
    from ARIAtools.ARIAProduct import ARIA_standardproduct
    print ('*****************************************************************')
    print ('*** Extract Product Function ***')
    print ('*****************************************************************')

    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile,
                                                bbox=inps.bbox,
                                                workdir=inps.workdir,
                                                num_threads=inps.num_threads,
                                                url_version=inps.version,
                                                nc_version=inps.nc_version,
                                                verbose=inps.verbose)

    # Perform initial layer, product, and correction sanity checks
    inps.layers, inps.tropo_total, \
        model_names = layerCheck(standardproduct_info.products[1],
                                      inps.layers,
                                      inps.nc_version,
                                      inps.gacos_products,
                                      inps.tropo_models,
                                      extract_or_ts = 'extract')

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

    # Extract
    # aria_extract default parms
    export_dict = {
        'full_product_dict': standardproduct_info.products[1],
        'bbox_file': standardproduct_info.bbox_file,
        'prods_TOTbbox': prods_TOTbbox,
        'layers': inps.layers,
        'rankedResampling': inps.rankedResampling,
        'dem': demfile,
        'lat': Latitude,
        'lon': Longitude,
        'mask': inps.mask,
        'outDir': inps.workdir,
        'outputFormat': inps.outputFormat,
        'stitchMethodType': inps.stitchMethodType,
        'verbose': inps.verbose,
        'num_threads': inps.num_threads,
        'multilooking': inps.multilooking,
        'tropo_total': inps.tropo_total,
        'model_names': model_names
    }

    # Extract user expected layers
    export_products(**export_dict)

    # If necessary, resample DEM/mask AFTER they have been used to extract metadata layers and mask output layers, respectively
    if inps.multilooking is not None:
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
    if inps.gacos_products:
        gacos_correction(standardproduct_info.products, inps.gacos_products,
                         standardproduct_info.bbox_file, prods_TOTbbox, outDir=inps.workdir,
                         outputFormat=inps.outputFormat, verbose=inps.verbose, num_threads=inps.num_threads)


def transformPoints(lats: np.ndarray, lons: np.ndarray, hgts: np.ndarray, old_proj: CRS, new_proj: CRS) -> np.ndarray:
    '''
    Transform lat/lon/hgt data to an array of points in a new
    projection
    Args:
        lats: ndarray   - WGS-84 latitude (EPSG: 4326)
        lons: ndarray   - ditto for longitude
        hgts: ndarray   - Ellipsoidal height in meters
        old_proj: CRS   - the original projection of the points
        new_proj: CRS   - the new projection in which to return the points
    Returns:
        ndarray: the array of query points in the weather model coordinate system (YX)
    '''
    t = Transformer.from_crs(old_proj, new_proj)

    # Flags for flipping inputs or outputs
    if not isinstance(new_proj, pyproj.CRS):
        new_proj = CRS.from_epsg(new_proj.lstrip('EPSG:'))
    if not isinstance(old_proj, pyproj.CRS):
        old_proj = CRS.from_epsg(old_proj.lstrip('EPSG:'))

    in_flip = old_proj.axis_info[0].direction
    out_flip = new_proj.axis_info[0].direction

    if in_flip == 'east':
        res = t.transform(lons, lats, hgts)
    else:
        res = t.transform(lats, lons, hgts)

    if out_flip == 'east':
        return np.stack((res[1], res[0], res[2]), axis=-1).T
    else:
        return np.stack(res, axis=-1).T
