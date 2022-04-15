#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import pdb
import os
import time
import sys
import logging

import numpy as np
from osgeo import gdal,ogr
from osgeo.gdalconst import *
import tempfile
import shutil
from shapely.geometry import Point,Polygon,shape, LinearRing
from shapely.ops import nearest_points
import copy
import json
from joblib import Parallel, delayed, dump, load
import random
import glob
import collections

from ARIAtools.logger import logger
from ARIAtools.shapefile_util import open_shapefile, save_shapefile

log = logging.getLogger(__name__)

solverTypes = ['pulp', 'glpk', 'gurobi']
redarcsTypes = {'MCF':-1, 'REDARC0':0, 'REDARC1':1, 'REDARC2':2}
stitchMethodTypes = ['overlap', '2stage']

class Stitching:
    """This is the parent class of all stiching codes.

    It is whatever is shared between the different stitching method variants
    e.g. - setting of the main input arguments,
         - functions to verify GDAL compatibility,
         - function to write new connected component of merged product
               based on a mapping table
         - function to write unwrapped phase of merged product based on
               a mapping table
    """

    def __init__(self):
        """Setting the default arguments needed by the class.

        Parse the filenames and bbox as None as they need to be set
            by the user,
        which will be caught when running the child classes of the
        respective stitch method.
        """
        self.inpFile = None
        self.ccFile = None
        self.prodbboxFile = None
        self.prodbbox = None
        self.solver = 'pulp'
        self.redArcs = -1
        self.mask = None
        self.outFileUnw = './unwMerged'
        self.outFileConnComp = './connCompMerged'
        self.outputFormat = 'ENVI'

        # stitching methods, by default leverage product overlap method
        # other options would be to leverage connected component
        self.setStitchMethod("overlap")

    def setInpFile(self, input):
        """Set the input Filename for stitching/unwrapping."""
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(input, np.str):
            input = [input]
        self.inpFile = input

        # The number of files that needs to be merged/ unwrapped
        self.nfiles = np.shape(self.inpFile)[0]

    def setOutFile(self, output):
        """Set the output File name."""
        self.outFile = output

    def setConnCompFile(self, connCompFile):
        """Set the connected Component file."""
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(connCompFile, np.str):
            connCompFile = [connCompFile]
        self.ccFile = connCompFile

    def setProdBBoxFile(self, ProdBBoxFile):
        """Set the product bounding box file(s)."""
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(ProdBBoxFile, np.str):
            ProdBBoxFile = [ProdBBoxFile]
        self.prodbboxFile = ProdBBoxFile

    def setBBoxFile(self, BBoxFile):
        """Set bounds of bbox."""
        self.bbox_file = BBoxFile

    def setTotProdBBoxFile(self, prods_TOTbbox):
        """Set common track bbox file."""
        self.setTotProdBBoxFile = prods_TOTbbox

    def setStitchMethod(self,stitchMethodType):
        """Set stitch method to be used to handle parent class internals."""
        if stitchMethodType not in stitchMethodTypes:
            raise ValueError(stitchMethodType + ' must be in '
                + str(stitchMethodTypes))
        else:
            self.stitchMethodType = stitchMethodType

        if self.stitchMethodType == '2stage':
            self.description="Two-stage corrected/stiched Unwrapped Phase"
        elif self.stitchMethodType == 'overlap':
            self.description = "Overlap-based stiched Unwrapped Phase"

    def setRedArcs(self, redArcs):
        """Set the Redundant Arcs to use for LP unwrapping."""
        self.redArcs = redArcs

    def setSolver(self, solver):
        """Set the solver to use for unwrapping."""
        self.solver = solver

    def setMask(self, mask):
        """Set the mask file."""
        self.mask = mask

    def setOutputFormat(self,outputFormat):
        """Set the output format of the files to be generated."""
        # File must be physically extracted, cannot proceed with VRT format.
        # Defaulting to ENVI format.
        self.outputFormat = outputFormat
        if self.outputFormat=='VRT':
            self.outputFormat='ENVI'

    def setOutFileUnw(self,outFileUnw):
        """Set the output file name for the unw stitched file."""
        self.outFileUnw = outFileUnw

    def setOutFileConnComp(self,outFileConnComp):
        """Set the output file name for the conncomp stitched file."""
        self.outFileConnComp = outFileConnComp

    def setVerboseMode(self,verbose):
        """Set verbose output mode."""
        logger.setLevel(logging.DEBUG)

    def __verifyInputs__(self):
        """Verify if the unw and conncomp inputs are GDAL-compatible
        and that the provided shape files are well-formed.

        If not remove them from the list to be stitched.
        If a VRT exists and is GDAL-compatible update the file to be a VRT.
        """
        # Track a list of files to keep
        inpFile_keep = []
        ccFile_keep = []
        prodbboxFile_keep = []
        bbox_keep = []
        for k_file in range(self.nfiles):
            # unw and corresponding conncomponent file
            inFile = self.inpFile[k_file]
            ccFile = self.ccFile[k_file]
            # Shape file in case passed through
            if self.stitchMethodType == "overlap":
                prodbboxFile = self.prodbboxFile[k_file]

            # vrt files prefered as they contain proj and transf information
            # Convert inputs to vrt if exist, and check for gdal compatibility
            inFile_temp = gdalTest(inFile)
            ccFile_temp = gdalTest(ccFile)
            # check if it exist and well-formed and pass back shapefile of it
            if self.stitchMethodType == "overlap":
                bbox_temp = open_shapefile(
                                prodbboxFile,'productBoundingBox',
                                1)
                prodbboxFile_temp = prodbboxFile
            else:
                bbox_temp= "Pass"
                prodbboxFile_temp = "Pass"

            # if one of the two files fails,
            # then do not try and merge them later on as GDAL is leveraged
            if inFile_temp is None or ccFile_temp is None:
                log.info('Removing following pair combination '
                         '(not GDAL compatible)')
                log.info('UNW: %s', inFile)
                log.info('CONN: %s', ccFile)
            elif prodbboxFile_temp is None:
                log.info('Removing following pair combination '
                         '(malformed shapefile)')
                log.info('UNW: %s', inFile)
                log.info('CONN: %s', ccFile)
            else:
                inpFile_keep.append(inFile_temp)
                ccFile_keep.append(ccFile)
                # shape file in case passed through
                if self.stitchMethodType == "overlap":
                    bbox_keep.append(bbox_temp)
                    prodbboxFile_keep.append(prodbboxFile_temp)

        # Update the input files and only keep those that are GDAL-compatible
        self.inpFile=inpFile_keep
        self.ccFile=ccFile_keep

        # Shape file in case passed through
        if self.stitchMethodType == "overlap":
            self.prodbbox=bbox_keep
            self.prodbboxFile=prodbboxFile_keep

        # Update the number of files in case some got removed
        self.nfiles = np.shape(self.inpFile)[0]

        if self.nfiles == 0:
            log.info('No files left after GDAL compatibility check')
            sys.exit(0)

    def __createImages__(self):
        """This function will write the final merged unw and conncomp file.
        
        As intermediate step tiff files are generated
        with integer values which will represent the shift to be applied
        to connected componenet and the moduli shift to be applied to
        the unwrapped phase.
        """
        ## Will first make intermediate files in a temp folder.
        # For each product there will be 1 connComp and 3 related unw files
        # The connComp file will show the unique wrt to all merged files
        # For the unw related files, there will be an integer offset tif file,
        # a vrt file which scale this integer map by 2pi,
        # and a vrt which combines the orginal unw file with the scaled map.
        # The latter will be used for merging of the unwrapped phase.
        tempdir = tempfile.mkdtemp(prefix='IntermediateFiles_', dir='.')

        # Will try multi-core version,
        # but default to for loop in case of failure
        try:
            # Need to combine all inputs together as single argument tuple
            all_inputs = ()
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_" + str(counter)
                fileMappingDict['description'] = self.description
                # Parse inputs as a tuple
                inputs = (fileMappingDict)
                # Append all tuples in a single tuple
                all_inputs = all_inputs + (inputs,)
            # Compute the phase value using multi-thread functionality
            intermediateFiles = Parallel(n_jobs=-1, max_nbytes=1e6) \
                (delayed(createConnComp_Int)(ii) for ii in all_inputs)

        except:
            log.info('Multi-core version failed, will try single for loop')
            intermediateFiles = []
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_n" + str(counter)
                fileMappingDict['description'] = self.description
                # Parse inputs as a tuple
                inputs = (fileMappingDict)
                # Compute the phase value
                intermediateFiles_temp = createConnComp_Int(inputs)
                intermediateFiles.append(intermediateFiles_temp)

        # Combine all connComp and unw files that need to be blended
        connCompFiles = []
        unwFiles = []
        for intermediateFile in intermediateFiles:
            connCompFiles.append(intermediateFile[0])
            unwFiles.append(intermediateFile[1])

        # Check if the folder to which files are being generated exists
        outPathUnw = os.path.dirname(os.path.abspath(self.outFileUnw))
        outPathConnComp = os.path.dirname(
                              os.path.abspath(
                              self.outFileConnComp))
        if not os.path.isdir(outPathUnw):
            os.makedirs(outPathUnw)
        if not os.path.isdir(outPathConnComp):
            os.makedirs(outPathConnComp)

        ## Will now merge the unwrapped and connected component files
        # Remove existing unwrapped output file(s)
        for file in glob.glob(self.outFileUnw + "*"):
            os.remove(file)

        # Create new output file(s)
        vrtName = '{:s}.vrt'.format(self.outFileUnw)
        gdal.BuildVRT(vrtName, unwFiles,
            options=gdal.BuildVRTOptions(srcNodata=0))
        gdal.Warp(self.outFileUnw, vrtName,
            options=gdal.WarpOptions(format=self.outputFormat,
                                    cutlineDSName=self.setTotProdBBoxFile,
                                    outputBounds=self.bbox_file))
        # Update VRT
        gdal.Translate(vrtName, self.outFileUnw,
            options=gdal.TranslateOptions(format="VRT"))
        # Apply mask (if specified)
        if self.mask is not None:
            update_file = gdal.Open(self.outFileUnw,gdal.GA_Update)
            update_file = update_file.GetRasterBand(1).WriteArray(
                             self.mask.ReadAsArray() *
                             gdal.Open(self.outFileUnw+'.vrt').ReadAsArray())
            update_file = None

        # Remove existing connected component output file(s)
        for file in glob.glob(self.outFileConnComp + "*"):
            os.remove(file)
        gdal.BuildVRT(self.outFileConnComp+'.vrt',
                      connCompFiles,
                      options=gdal.BuildVRTOptions(srcNodata=-1))
        gdal.Warp(self.outFileConnComp,
                  self.outFileConnComp+'.vrt',
                  options=gdal.WarpOptions(
                      format=self.outputFormat,
                      cutlineDSName=self.setTotProdBBoxFile,
                      outputBounds=self.bbox_file))
        # Update VRT
        gdal.Translate(self.outFileConnComp+'.vrt',
                       self.outFileConnComp,
                       options=gdal.TranslateOptions(format="VRT"))
        # Apply mask (if specified).
        if self.mask is not None:
            update_file = gdal.Open(self.outFileConnComp, gdal.GA_Update)
            #mask value for conncomp must be set to nodata value -1
            ma_update_file = np.ma.masked_where(
                self.mask.ReadAsArray() == 0.,
                    gdal.Open(
                        self.outFileConnComp+'.vrt').ReadAsArray())
            np.ma.set_fill_value(
                ma_update_file,
                update_file.GetRasterBand(1).GetNoDataValue())
            update_file = update_file.GetRasterBand(1).WriteArray(
                ma_update_file.filled())
            update_file = None

        cmd = "gdal_translate -of png -scale -ot Byte -q " \
                + self.outFileUnw + ".vrt " \
                + self.outFileUnw + ".png"
        os.system(cmd)

        # Remove directory with intermediate files, which are no longer needed
        shutil.rmtree(tempdir)


class UnwrapOverlap(Stitching):
    """Stiching/unwrapping using product overlap minimization."""

    def __init__(self):
        """Inherit properties from the parent class.

        Parse the filenames and bbox as None as they need to be set by
        the user, which will be caught when running the class.
        """
        Stitching.__init__(self)

    def UnwrapOverlap(self):
        ## Set the method
        self.setStitchMethod("overlap")

        ## Check if required inputs are set
        if self.inpFile is None:
            log.error("Input unwrapped file(s) is (are) not set.")
            raise Exception
        if self.ccFile is None:
            log.error("Input Connected Components file(s) is (are) not set.")
            raise Exception
        if self.prodbboxFile is None:
            log.error("Input product Bounding box file(s) is (are) not set.")
            raise Exception

        ## Verify if all the inputs are well-formed/GDAL compatible
        # Update files to be vrt if they exist and remove files which
        # failed the gdal compatibility
        self.__verifyInputs__()

        ## Calculating the number of phase cycles needed to miminize
        # the residual between products
        self.__calculateCyclesOverlap__()

        ## Write out merged phase and connected component files
        self.__createImages__()

        return

    def __calculateCyclesOverlap__(self):
        """Function that will calculate the number of cycles each
        component needs to be shifted in order to minimize the two-pi
        modulo residual between a neighboring component.

        Outputs a fileMappingDict with a file number as a key. Within
        fileMappingDict with a integer phase shift value for each unique
        connected component.
        """
        # Only need to comptue the minimize the phase offset if the number
        # of files is larger than 2
        if self.nfiles > 1:
            # Initiate the residuals and design matrix
            residualcycles = np.zeros((self.nfiles-1, 1))
            residualrange = np.zeros((self.nfiles-1, 1))
            A = np.zeros((self.nfiles-1, self.nfiles))

            # The files are already sorted in the ARIAproduct class, will
            # make consecutive overlaps between these sorted products
            for counter in range(self.nfiles-1):
                # Getting the two neighboring frames
                bbox_frame1 = self.prodbbox[counter]
                bbox_frame2 = self.prodbbox[counter+1]

                # Determine the intersection between the two frames
                if not bbox_frame1.intersects(bbox_frame2):
                    log.error("Products do not overlap or were not "
                        "provided in a contigious sorted list.")
                    raise Exception
                polyOverlap = bbox_frame1.intersection(bbox_frame2)

                # Will save the geojson under a temp local filename
                # Do this just to get the file outname
                tmfile = tempfile.NamedTemporaryFile(mode='w+b',
                             suffix='.json',
                             prefix='Overlap_', dir='.')
                outname = tmfile.name

                # Remove it as GDAL polygonize function cannot overwrite files
                tmfile.close()
                tmfile = None

                # Save the temp geojson
                save_shapefile(outname, polyOverlap, 'GeoJSON')

                # Calculate the mean of the phase for each product in
                # the overlap region alone.
                # Will first attempt to mask out connected component 0,
                # and default to complete overlap if this fails.
                # Cropping the unwrapped phase and connected component
                # to the overlap region alone, inhereting the no-data.

                # Connected component
                _, connCompNoData1, _, _ = GDALread(
                        self.ccFile[counter], data_band=1, loadData=False)
                _, connCompNoData2, _, _ = GDALread(
                        self.ccFile[counter+1], data_band=1, loadData=False)

                connCompFile1 = gdal.Warp('', self.ccFile[counter],
                        options=gdal.WarpOptions(format="MEM",
                                cutlineDSName=outname,
                                dstAlpha=True,
                                outputBounds=polyOverlap.bounds,
                                dstNodata=connCompNoData1))
                connCompFile2 = gdal.Warp('', self.ccFile[counter+1],
                        options=gdal.WarpOptions(format="MEM",
                                cutlineDSName=outname,
                                dstAlpha=True,
                                outputBounds=polyOverlap.bounds,
                                dstNodata=connCompNoData2))

                # Reformat output bounds for GDAL translate
                ulx, lry, lrx, uly = polyOverlap.bounds
                projWin = (ulx, uly, lrx, lry)

                # Unwrapped phase
                _, unwNoData1, _, _ = GDALread(
                        self.inpFile[counter], data_band=1, loadData=False)
                _, unwNoData2, _, _ = GDALread(
                        self.inpFile[counter+1], data_band=1, loadData=False)

                unwFile1 = gdal.Translate('', self.inpFile[counter],
                        options=gdal.TranslateOptions(format="MEM",
                        projWin=projWin,
                        noData=unwNoData1))
                unwFile2 = gdal.Translate('', self.inpFile[counter+1],
                        options=gdal.TranslateOptions(format="MEM",
                        projWin=projWin,
                        noData=unwNoData2))

                # Find the component with the largest overlap
                connCompData1 = connCompFile1.GetRasterBand(1).ReadAsArray()
                connCompData1[((connCompData1==connCompNoData1)
                                | (connCompData1==0))] = np.nan
                connCompData2 = connCompFile2.GetRasterBand(1).ReadAsArray()
                connCompData2[((connCompData2==connCompNoData2)
                                | (connCompData2==0))] = np.nan
                connCompData2_temp = (connCompData2*100)
                connCompDiff = (connCompData2_temp.astype(np.int)
                        - connCompData1.astype(np.int))
                connCompDiff[(connCompDiff<0) | (connCompDiff>2000)] = 0
                temp_count = collections.Counter(connCompDiff.flatten())
                maxKey = 0
                maxCount = 0
                for key, keyCount in temp_count.items():
                    if key != 0:
                        if keyCount  > maxCount:
                            maxKey   = key
                            maxCount = keyCount
                print('Max key, count:', maxKey, maxCount)

                # If the max key count is 0, this means there is no good
                # overlap region between products.
                # In that scenario default to different stitching approach.
                if maxKey != 0 and maxCount > 75:
                    # Masking the unwrapped phase and only use the
                    # largest overlapping connected component
                    unwData1 = unwFile1.GetRasterBand(1).ReadAsArray()
                    unwData2 = unwFile2.GetRasterBand(1).ReadAsArray()

                    cutlineMask1 = connCompFile1.GetRasterBand(
                                       2).ReadAsArray()
                    cutlineMask2 = connCompFile2.GetRasterBand(
                                       2).ReadAsArray()

                    unwData1[((unwData1==unwNoData1)
                                | (connCompDiff!=maxKey)
                                | (cutlineMask1==0))] = np.nan
                    unwData2[((unwData2==unwNoData2)
                                | (connCompDiff!=maxKey)
                                | (cutlineMask2==0))] = np.nan

                    # Calculation of the range correction
                    unwData1_wrapped = (unwData1
                                    - np.round(unwData1/(2*np.pi))*(2*np.pi))
                    unwData2_wrapped = (unwData2
                                    - np.round(unwData2/(2*np.pi))*(2*np.pi))
                    arr = unwData1_wrapped - unwData2_wrapped

                    # Data is not fully decorrelated
                    arr = arr - np.round(arr/(2*np.pi))*2*np.pi
                    range_temp = np.angle(np.nanmean(np.exp(1j*arr)))

                    # Calculation of the number of 2PI cycles accounting
                    # for range correction
                    corrected_range = unwData1 - (unwData2+range_temp)
                    cycles_temp = np.round((np.nanmean(corrected_range))
                                      / (2*np.pi))
                    print(cycles_temp)

                else:
                    # Account for the case that no data was left, e.g.
                    # fully decorrelated
                    # In that scenario use all data and estimate from
                    # wrapped, histogram will be broader ...
                    unwData1 = unwFile1.GetRasterBand(1).ReadAsArray()
                    unwData2 = unwFile2.GetRasterBand(1).ReadAsArray()

                    cutlineMask1 = connCompFile1.GetRasterBand(
                                       2).ReadAsArray()
                    cutlineMask2 = connCompFile2.GetRasterBand(
                                       2).ReadAsArray()

                    unwData1[((unwData1==unwNoData1)
                                | (cutlineMask1==0))] = np.nan
                    unwData2[((unwData2==unwNoData2)
                                | (cutlineMask2==0))] = np.nan

                    # Calculation of the range correction
                    unwData1_wrapped = (unwData1
                                        - np.round(unwData1/(2*np.pi))
                                        * (2*np.pi))
                    unwData2_wrapped = (unwData2
                                        - np.round(unwData2/(2*np.pi))
                                        * (2*np.pi))
                    arr = unwData1_wrapped - unwData2_wrapped
                    arr = arr - np.round(arr/(2*np.pi))*2*np.pi
                    range_temp = np.angle(np.nanmean(np.exp(1j*arr)))

                    # Data is decorelated assume no 2-pi cycles
                    cycles_temp = 0

                # Closing the files
                unwFile1 = None
                unwFile2 = None
                connCompFile1 = None
                connCompFile2 = None

                # Remove the tempfile
                shutil.os.remove(outname)

                # Store the residual and populate the design matrix
                residualcycles[counter] = cycles_temp
                residualrange[counter] = range_temp
                A[counter,counter] = 1
                A[counter,counter+1] = -1

            # Invert the offsets with respect to the first product
            cycles = np.round(np.linalg.lstsq(A[:,1:],
                            residualcycles,rcond=None)[0])
            rangesoffset = np.linalg.lstsq(A[:,1:],
                            residualrange,rcond=None)[0]
            #pdb.set_trace()

            # Force first product to have 0 as offset
            cycles = -1*np.concatenate((np.zeros((1,1)), cycles), axis=0)
            rangesoffset = -1*np.concatenate((np.zeros((1,1)), rangesoffset),
                               axis=0)

        else:
            # Nothing to be done, i.e. no phase cycles to be added
            cycles = np.zeros((1,1))
            rangesoffset = np.zeros((1,1))

        # Build the mapping dictionary
        fileMappingDict = {}
        connCompOffset = 0

        for fileCounter in range(self.nfiles):
            # Get the number of connected components
            n_comp = 20

            # The original connected components
            connComp = np.array(range(0,n_comp+1))
            connComp = connComp.reshape((n_comp+1,1))

            # The cyclemapping of product using a single offset
            # (i.e. all connectedComponents have the same mapping)
            cycleMapping = np.ones((n_comp+1,1))*cycles[fileCounter,0]
            rangeOffsetMapping = np.ones((n_comp+1,1)) \
                                     * rangesoffset[fileCounter,0]

            # Generate the mapping of connectedComponents
            # Increment based on unique components of the merged product,
            # 0 is grouped for all products
            connCompMapping = connComp
            connCompMapping[1:] = connComp[1:]+connCompOffset

            # concatenate and generate a connComponent based mapping matrix
            # [original comp, new component, 2pi unw offset]
            connCompMapping = np.concatenate(
                (connComp, connCompMapping, cycleMapping, rangeOffsetMapping),
                axis=1)

            # Increment the count of total number of unique components
            connCompOffset = connCompOffset + n_comp

            # Populate mapping dictionary for each product
            fileMappingDict_temp = {}
            fileMappingDict_temp['connCompMapping'] = connCompMapping
            fileMappingDict_temp['connFile'] = self.ccFile[fileCounter]
            fileMappingDict_temp['unwFile'] = self.inpFile[fileCounter]

            # Store it in the general mapping dictionary
            fileMappingDict[fileCounter] = fileMappingDict_temp

        # Pass the fileMapping back into self
        self.fileMappingDict = fileMappingDict


class UnwrapComponents(Stitching):
    """Stiching/unwrapping using 2-Stage Phase Unwrapping."""

    def __init__(self):
        """
        Inheret properties from the parent class.
        Parse the filenames and bbox as None as they need to be set by
        the user, which will be caught when running the class.
        """
        Stitching.__init__(self)


    def unwrapComponents(self):
        ## Setting the method
        self.setStitchMethod("2stage")
        self.region=5

        if self.inpFile is None:
            log.error("Input unwrapped file(s) is (are) not set.")
            raise Exception

        if self.ccFile is None:
            log.error("Input Connected Components file(s) is (are) not set.")
            raise Exception

        if self.solver not in solverTypes:
            raise ValueError(self.treeType + ' must be in '
                + str(unwTreeTypes))

        if self.redArcs not in redarcsTypes.keys():
            raise ValueError(self.redArcs + ' must be in '
                + str(redarcsTypes))


        ## Verify if all the inputs are well-formed/GDAL compatible
        # Update files to be vrt if they exist and remove files
        # which failed the gdal compatibility
        self.__verifyInputs__()

        ## Loop over the number of valid files and populate
        ##     a combined dictionary with keys being polygons
        ##     around connected component edges in "polyTableDict" variable.
        #  Key is a unique polygon
        #      [1, max number of unique polygons of all products]
        #  Note: there can be components with mutiple polygons
        #      (e.g. inner/outer edges)
        #  Each key will have a corresponding dictionary with keys:
        #       "connFile" = Original connected component file
        #       "connCompID" = Original ID connected component id
        #           (wrt original connected component file)
        #       "unwFile" = Original unwrapped file
        self.__populatePolyTable__()

        ## Populate a table where each unique polygon is matches with all
        ##     other (not same component) polygons, and keep in the table
        ##     the points that represent the closest two points
        ##     between these polygons.
        # The table is populated with columns as:
        #       Source Unique CC ID, Source CC ID, Source CC FILE,
        #       Source UNW FILE, Source X-geocoordinate,
        #       Source Y-geocoordinate, Source geoTrans, Partner Unique CC,
        #       Dist to partner, Source unw phase, Source X-coordinate,
        #       Source Y-coordinate
        self.__populatePointTable__()

        ## Calculating the number of phase cycles needed to miminize
        ##     the residual between connected components.
        ##     Minimization is only done between the tablePoints
        ##     (i.e. X-coordinate, Y-coordinate, unw phase).
        self.__calculateCycles__()

        ## Write out merged phase and connected component files
        self.__createImages__()

        return

    def __populatePolyTable__(self):
        """
        Build a dictionary with connected component polygon information
        columns: key is a unique polygon
            [1, max number of unique polygons of all products]
        each key will have a corresponding dictionary with keys:
            "connCompFile" = Original connected component file
            "conmCompID" = Original ID connected component id
                (wrt original connected component file)
            "unwFile" = Original unwrapped file
        """

        # template of the dictionary for each polygon
        polygonDict_tmpl = {}
        polygonDict_tmpl['connFile'] = []
        polygonDict_tmpl['connCompID'] = []
        polygonDict_tmpl['connCompID_unique'] = []
        polygonDict_tmpl['unwFile'] = []
        polygonDict_tmpl['poly'] = []
        polygonDict_tmpl['geoTrans'] = []

        # template of the dictionary for each unique connected component
        componentDict_tmpl = {}
        componentDict_tmpl['connFile'] = []
        componentDict_tmpl['connCompMapping'] = []
        componentDict_tmpl['unwFile'] = []


        # populate table dictionary of all polygons
        tableDict = {}

        fileMappingDict ={}
        polygon_offset = 0
        connComp_offset = 0
        for k_file in range(self.nfiles):

            # polygonize each connected component file
            # outputs polygons for all connected components differnt than 0 and no-data value
            poly,ccomp,geoTrans,proj,sizeData = polygonizeConn(self.ccFile[k_file])

            # tracking for each file the mapping of connected component to a unique connected component defined below.
            componentDict = copy.deepcopy(componentDict_tmpl)
            componentDict['connFile'] = self.ccFile[k_file]
            componentDict['unwFile'] = self.inpFile[k_file]
            connCompMapping = [[0 , 0]]

            connComp_offset_increment = 0
            for poly_counter in range(len(poly)):
                # populate the dict with information specific to this polygon
                polygonDict =  copy.deepcopy(polygonDict_tmpl)
                polygonDict['connFile'] = self.ccFile[k_file]
                polygonDict['unwFile'] = self.inpFile[k_file]
                polygonDict['connCompID'] = ccomp[poly_counter]
                polygonDict['connCompID_unique'] = ccomp[poly_counter]+connComp_offset
                polygonDict['poly'] = poly[poly_counter]
                polygonDict['geoTrans'] = geoTrans
                polygonDict['sizeData'] = sizeData

                # Track the mapping from connCompID to connCompID_unique
                connCompMapping.append([ccomp[poly_counter] , ccomp[poly_counter]+connComp_offset])

                # tracking how much the next file needs to be offsetted in connComp count
                connComp_offset_increment = np.max([connComp_offset_increment,ccomp[poly_counter]])

                # store the polygonDict for each unigue polygon
                tableDict[poly_counter+polygon_offset]=polygonDict

            # storing the information for the file mapping
            # as there might be some repetition in the mapping (gdal polygonize might have >1 polygon for a given component)
            connCompMapping = np.array(connCompMapping)
            connCompMapping = np.unique(connCompMapping, axis=0)
            componentDict['connCompMapping'] = connCompMapping
            fileMappingDict[k_file]=componentDict

            # increment the polygon count for polygons to remain unique
            polygon_offset = polygon_offset+poly_counter+1
            # increment the connected component count for it to remain unique
            connComp_offset = connComp_offset + connComp_offset_increment

        # return the dictionary of all unique polygons
        self.polyTableDict = tableDict
        # return the dictionary of all unqiue connected component id
        self.fileMappingDict = fileMappingDict


    def __populatePointTable__(self):
        """
        Function that populates a table of points which maps
        two closest points between two connected components.
        Loops over all possible connected component pairings.
        The table is populated with columns as:
            Source Unique CC ID
            Source CC ID
            Source CC FILE
            Source UNW FILE
            Source X-geocoordinate
            Source Y-geocoordinate
            Source geoTrans
            Partner Unique CC
            Dist to partner
            Source unw phase
            Source X-coordinate
            Source Y-coordinate
        """

        # only proceed if there are polygons that needs mapping
        #     to a closest polygon
        if len(self.polyTableDict) <=1:
            raise Exception("Nothing to do as all are from "
                            "same connected component")

        # populate table  of all polygons
        polygon_pair_counter = 0
        # loop over each polygon and compute the distance to all polygons.
        #print('Calculating the phase values for polygon pairs')
        for poly_counter_1 in range(0,len(self.polyTableDict)):
            # polygon 1
            poly1 = self.polyTableDict[poly_counter_1]['poly']
            connCompID1 = self.polyTableDict[poly_counter_1]['connCompID']
            connCompID_unique1 = self.polyTableDict[
                                     poly_counter_1]['connCompID_unique']
            connFile1 = self.polyTableDict[poly_counter_1]['connFile']
            unwFile1 = self.polyTableDict[poly_counter_1]['unwFile']
            geoTrans1 = self.polyTableDict[poly_counter_1]['geoTrans']
            sizeData1 = self.polyTableDict[poly_counter_1]['sizeData']

            # loop over all poygons and keep track of min distance polygon
            #     not of own component
            for poly_counter_2 in range(poly_counter_1,len(self.polyTableDict)):
                # polygon 2
                poly2 = self.polyTableDict[poly_counter_2]['poly']
                connCompID2 = self.polyTableDict[poly_counter_2]['connCompID']
                connCompID_unique2 = self.polyTableDict[
                                         poly_counter_2]['connCompID_unique']
                connFile2 = self.polyTableDict[poly_counter_2]['connFile']
                unwFile2 = self.polyTableDict[poly_counter_2]['unwFile']
                geoTrans2 = self.polyTableDict[poly_counter_2]['geoTrans']
                sizeData2 = self.polyTableDict[poly_counter_2]['sizeData']

                # will skip the pairing of the two polygons under the
                #     following two conditions:
                # 1) when mapping two of the same polyons, i.e. i==j or
                #     a diagonal element of the matrix of combinations
                # 2) when the 2 polygons are from the same connected component
                #     and originate from the same connected comp file
                if poly_counter_1==poly_counter_2 and \
                       connCompID1==connCompID2 and \
                       connFile1==connFile2:
                    continue

                # finding the minMatch of closest points between
                #     two polygons lists
                # will try multi-core version and default to for
                #     loop in case of failure
                try:
                    log.info('Multi-core')
                    start = time.time()
                    points = []
                    pointDistance = []
                    pairings =  ()
                    for polygon1 in poly1:
                        for polygon2 in poly2:
                            # parse inputs as a tuple
                            pairing = (polygon1,polygon2,connFile1,connFile2)
                            # append all tuples in a single tuple
                            pairings = pairings + (pairing,)
                    temp = Parallel(n_jobs=-1,max_nbytes=1e6)(delayed(
                               minDistancePoints)(ii) for ii in pairings)
                    # unpacking the tuple back to a numpy array of
                    #     two variables
                    for line in temp:
                        pointDistance.append(line[0])
                        points.append(line[1])
                    stop = time.time()
                except:
                    log.info('For loop')
                    # for loop instead of multi-core
                    start = time.time()
                    points = []
                    pointDistance = []
                    for polygon1 in poly1:
                        for polygon2 in poly2:
                            # parse inputs as a tuple
                            pairing = (polygon1,polygon2,connFile1,connFile2)
                            pointDistance_temp, points_temp = \
                                minDistancePoints(pairing)
                            points.append(points_temp)
                            pointDistance.append(pointDistance_temp)
                    stop = time.time()
                log.info('DONE min distance points')
                # now find the point pair with the shortest distance and
                #     keep the correpsonding points
                min_ix = np.where(pointDistance ==np.min(pointDistance))[0]
                dist = pointDistance[min_ix[0]]
                points = points[min_ix[0]]

                # track the coordinates of the closest points
                #     in geocoordinates
                point1 = np.array(points[0])
                point2 = np.array(points[1])

                # The list of points are defined from one polygon point
                #     to another.
                # Therefore each polygon mapping will consist out of
                #     two Source points.
                # The columns for each point are populated as:
                    # Source Unique CC ID, Source CC ID, Source CC FILE,
                    # Source UNW FILE, Source X-coordinate,
                    # Source Y-coordinate, Source geoTrans,
                    # Partner Unique CC, Dist to partner
                # Starting from point 1
                point1_array=[]
                point1_array.append(connCompID_unique1)
                point1_array.append(connCompID1)
                point1_array.append(connFile1)
                point1_array.append(unwFile1)
                point1_array.append(point1[0])
                point1_array.append(point1[1])
                point1_array.append(str(geoTrans1))
                point1_array.append(connCompID_unique2)
                point1_array.append(dist)
                point1_array = np.array(point1_array)
                point1_array = point1_array.reshape((1,point1_array.shape[0]))

                # starting from point 2
                point2_array=[]
                point2_array.append(connCompID_unique2)
                point2_array.append(connCompID2)
                point2_array.append(connFile2)
                point2_array.append(unwFile2)
                point2_array.append(point2[0])
                point2_array.append(point2[1])
                point2_array.append(str(geoTrans2))
                point2_array.append(connCompID_unique1)
                point2_array.append(dist)
                point2_array = np.array(point2_array)
                point2_array = point2_array.reshape((1,point2_array.shape[0]))

                # append the list of points
                if 'tablePoints' in locals():
                    tablePoints = np.concatenate((tablePoints,
                                                  point1_array,
                                                  point2_array), axis=0)
                else:
                    tablePoints = np.concatenate((point1_array,
                                                  point2_array), axis=0)


        ## making sure the points are unique based on coordinate
        ##    and unique connected compoenent id
        _, unique_indices = np.unique(
                                                tablePoints[:, (0,4,5)],
                                                axis=0,return_index=True)
        self.tablePoints = tablePoints[unique_indices[:],:]

        # compute the phase values
        log.info('Calculating the phase values for %s points',
                     self.tablePoints.shape[0])
        # will try multi-core version and default to for loop
        #     in case of failure
        try:
            start = time.time()
            # need to combine all inputs together as single argument tuple
            all_inputs = ()
            for counter in range(self.tablePoints.shape[0]):
                connCompID1_temp = self.tablePoints[counter,1]
                connFile1_temp = self.tablePoints[counter,2]
                unwFile1_temp = self.tablePoints[counter,3]
                point1_temp = np.array([self.tablePoints[counter,4],
                                       self.tablePoints[counter,5]])
                geoTrans1_temp =self.tablePoints[counter,6]
                # parse inputs as a tuple
                inputs = (connCompID1_temp,
                          connFile1_temp,
                          unwFile1_temp,
                          point1_temp,
                          geoTrans1_temp,
                          self.region)
                # append all tuples in a single tuple
                all_inputs = all_inputs + (inputs,)
            # compute the phase value using multi-thread functionality
            unwPhase_value = Parallel(n_jobs=-1,max_nbytes=1e6)(
                                 delayed(point2unwPhase)(ii)
                                     for ii in all_inputs)
            # end timer
            end = time.time()
        except:
            log.info('Multi-core version failed, '
                     'will try single for loop')
            # will track processing time
            start = time.time()
            unwPhase_value = []
            for counter in range(self.tablePoints.shape[0]):
                connCompID1_temp = self.tablePoints[counter,1]
                connFile1_temp = self.tablePoints[counter,2]
                unwFile1_temp = self.tablePoints[counter,3]
                point1_temp = np.array([self.tablePoints[counter,4],
                                       self.tablePoints[counter,5]])
                geoTrans1_temp =self.tablePoints[counter,6]
                # parse inputs as a tuple
                inputs = (connCompID1_temp,
                          connFile1_temp,
                          unwFile1_temp,
                          point1_temp,
                          geoTrans1_temp,
                          self.region)
                # compute the phase value
                unwPhase_temp = point2unwPhase(inputs)
                unwPhase_value.append(unwPhase_temp)
            # end timer
            end = time.time()
        log.info("runtime phase value calculation: %s sec", end - start)

        # Appending the phase to the table stored in self with the points
        unwPhase_value = np.array(unwPhase_value)
        unwPhase_value = unwPhase_value.reshape((unwPhase_value.shape[0],1))
        self.tablePoints = np.concatenate((self.tablePoints, unwPhase_value),
                                           axis=1)

        ## The two-stager only works with integer coordinate system
        ## in its triangulation, and needs unique points for each node.
        ## Will apply two steps to ensure this
        # STEP 1) will convert the geo-coordinates WRT the most
           # southwest point and normalize WRT original grid sampling
        # STEP 2) for non-unique points will shift by a pixel
           # to get uniquesness back

        # step 1: generating raster coordinates
        geoTrans_global = list(geoTrans1)
        geoTrans_global[0] = np.min(self.tablePoints[:,4].astype(np.float64))
        geoTrans_global[3] = np.max(self.tablePoints[:,5].astype(np.float64))
        X_coordinate = ((self.tablePoints[:,4].astype(np.float64)
                         - geoTrans_global[0])
                         / geoTrans_global[1]).astype(np.int)
        Y_coordinate = ((self.tablePoints[:,5].astype(np.float64)
                         - geoTrans_global[3])
                         / geoTrans_global[5]).astype(np.int)
        X_coordinate = X_coordinate.reshape((X_coordinate.shape[0],1))
        Y_coordinate = Y_coordinate.reshape((Y_coordinate.shape[0],1))

        # step 2: making the nodes unqiue
        test_unique =  np.array([])
        test = np.array([1,1])
        while_count = 1
        while test_unique.shape!=test.shape:
            for count in range(X_coordinate.shape[0]):
                # find all points that have the same coordiante as this one
                point = np.array([X_coordinate[count][0],
                                 Y_coordinate[count][0]])
                points = np.concatenate((X_coordinate, Y_coordinate), axis=1)
                matches = np.where(np.all(point==points,axis=1))[0]
                if matches.shape[0]>1:
                    #print('Found ' + str(matches.shape[0]) + ' non-unique points')
                    X_coordinate[matches[1:],0] = \
                        X_coordinate[matches[1:],0] \
                        + random.sample(range(-3,3),matches.shape[0]-1)
                    Y_coordinate[matches[1:],0] = \
                        Y_coordinate[matches[1:],0] \
                        + random.sample(range(-3,3),matches.shape[0]-1)

            test = np.concatenate((X_coordinate, Y_coordinate), axis=1)
            test_unique, _ = np.unique(test,
                                                         axis=0,
                                                         return_index=True)
            log.info('Shifting coordinates by a pixel to make '
                     'all points unique: iteration {}'.format(while_count))
            while_count = while_count+1

        if test_unique.shape!=test.shape:
            log.error('Coordinates are not-unique, need to implement '
                      'a while loop, run again for another random sample...')
            import pdb
            pdb.set_trace()
        else:
            log.info('Coordinates are unique')

        # appending the X and Y-coordinates to the points table.
        self.tablePoints = np.concatenate((self.tablePoints, X_coordinate),
                                           axis=1)
        self.tablePoints = np.concatenate((self.tablePoints, Y_coordinate),
                                           axis=1)


    def __compToCycle__(self, cycle, compnum):
        """
        Generate a mapping from unique connected component to the
        number of phase cycles that needs to be applied
        """

        compMap = {}
        for n, c in zip(cycle, compnum):
            try:
                compN = compMap[c]
                # Test if same cycle
                if (compN == n):
                    continue
                else:
                    raise ValueError('Incorrect phaseunwrap output: '
                                     'Different cycles in same components')
            except:
                # Component cycle doesn't exist in the dictionary
                compMap[c] = n

        return compMap


    def __calculateCycles__(self):
        """
        Function that will calculate the number of cycles each component needs
        to be shifted in order to minimize the two-pi modulu residual
        between a neighboring component.
        Outputs a fileMappingDict with as key a file number.
        Within fileMappingDict with a integer phase shift value for each
        unique connected component.
        """

        # import dependecies, looks like the general import
        # does not percolate down to this level
        from ARIAtools.phaseMinimization import PhaseUnwrap

        # generating a list to be parsed in the two-stager,
        # not that table originally was a mixture of string and floats,
        # so all were temporaly mapped to a string for easy parsing.
        # Will need to undo that now when calling two stager
        x = list(self.tablePoints[:,10].astype(np.int)+1)
        y = list(self.tablePoints[:,11].astype(np.int)+1)
        phase = (self.tablePoints[:,9].astype(np.float32))
        compNum = list(self.tablePoints[:,7].astype(np.int))
        phaseunwrap = PhaseUnwrap(
                                  x=x,
                                  y=y,
                                  phase=phase,
                                  compNum=compNum,
                                  redArcs=redarcsTypes[self.redArcs])
        phaseunwrap.solve(self.solver)
        cycles = phaseunwrap.unwrapLP()

        # Map unique component to integer number of cycles
        compMap = self.__compToCycle__(cycles, compNum)

        # Map component number for each file to integer number of cycles
        for countFile in self.fileMappingDict:
            fileMappingDict = self.fileMappingDict[countFile]

            # generate an array mapping connected component of each file
            # to the integer phase cycles it needs to be shifted
            cycleMapping = []
            for connCompID_unique in fileMappingDict['connCompMapping'][:,1]:
                if connCompID_unique in compMap:
                    cycleMapping.append(compMap[connCompID_unique])
                else:
                    cycleMapping.append(0)
                    if connCompID_unique != 0 :
                        log.error('Was not expecting this')
                        pdb.set_trace()
            cycleMapping = np.array(cycleMapping)
            cycleMapping = cycleMapping.reshape((cycleMapping.shape[0],1))
            connCompMapping = np.concatenate((
                fileMappingDict['connCompMapping'],
                    cycleMapping), axis=1)

            # update the original dictionary and store it back in self
            fileMappingDict['connCompMapping'] = connCompMapping
            self.fileMappingDict[countFile]=fileMappingDict
        log.info('MAPPING complete.')



def minDistancePoints(pairing):
    """Finding closest two points between a list of two polygons."""
    poly1= pairing[0]
    poly2 = pairing[1]
    poly1_source = pairing[2]
    poly2_source = pairing[3]

    # Shapely nearest_point behaves well when two polygons do not overlap.
    # in case one polygon is contained within another,
    # it will return the point on the inner polygon
    # the latter is not good behavior in case of
    # nested components within components
    # will handle this manully by checking if the polygons have intersection.
    # However, note that in case of stiching multiple products,
    # there will be a case that polygons intersect in along the edges
    # in the overlap region.
    # We will track this by also including a check if the polygons originate
    # from the same connected component file.

    ## check of the polygon has an intersection, two scenarios exist
    # 1) if from the same connected component file, then one component
       # is within another, nearest points will return incorrect points
    # 2) if from different connected component files, then it is from
       # the overlap region, nearest points will not give an incorrect
       # estimate as there is overlap
    
    #pdb.set_trace()
    Poly2 = Polygon(poly2)
    Poly1 = Polygon(poly1)

    if Poly1.intersection(Poly2).area > 0 and  poly1_source == poly2_source:
    #if poly1.intersects(poly2) and poly1_source == poly2_source:
        # will manully find the nearest points by looping over each vertex
        # of polygon 2 and finding the minimum distance to polygon 1.

        point1_list = []
        point2_list = []
        dist_list = []
        for point_coordinate in poly2.coords:
            point2 = Point(point_coordinate)
            d = poly1.project(point2)
            point1 = poly1.interpolate(d)
            dist = point1.distance(point2)
            # track the paramters in a list
            dist_list.append(dist)
            point1_list.append(point1)
            point2_list.append(point2)

        # selecting only the points closest to each other
        min_ix = np.where(dist_list ==np.min(dist_list))[0]
        dist = dist_list[min_ix[0]]
        points = tuple([point1_list[min_ix[0]], point2_list[min_ix[0]]])
        log.info('intersect: min distance %s', dist)

    else:
        if Poly1.intersection(Poly2).area > 0:
            log.info("product overlap")
            pdb.set_trace()

        # leverage shapely distance function
        dist = poly1.distance(poly2)

        # track the coordinates of the closest points in geocoordinates
        points = nearest_points(poly1,poly2)
        point1 = np.array(points[0])
        point2 = np.array(points[1])

    #pdb.set_trace()

    return dist, points

def polygonizeConn(ccFile):
    """
    Polygonize a connected component image.
    Ignores polygons for no-data value and connected component 0
    Takes a GDAL compatible file as input, and returns shaply polygons
    and a list of corresponding connected components
    """

    # will save the geojson under a temp local filename
    tmfile = tempfile.NamedTemporaryFile(mode='w+b',
                                         suffix='.json',
                                         prefix='ConnPolygonize_',
                                         dir='.')
    outname = tmfile.name
    # will remove it as GDAL polygonize function cannot overwrite files
    tmfile.close()
    tmfile = None

    ## Run the GDAL polygonize functionality
    # The no-data value will be ignored by default
    drv = ogr.GetDriverByName("GeoJSON")
    dst_ds = drv.CreateDataSource(outname)
    dst_layer = dst_ds.CreateLayer("out", srs = None)
    dst_layer.CreateField(ogr.FieldDefn("DN", ogr.OFTInteger))
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("DN")
    src_ds = gdal.Open(ccFile)
    srcband = src_ds.GetRasterBand(1)
    gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
    del drv, dst_ds, dst_layer, dst_field, src_ds, srcband
    dt = json.load(open(outname))
    os.remove(outname)

    # create a list of polygons and connected components asociated with them
    poly = []
    ccomp = []
    # default structure is features/properties/geometry/coordinates
    for ft in dt['features']:
        # Connected component 0 is an actual data value,
        #     but not used for 2-stager and is removed here
        if ft['properties']['DN'] <= 0:
            continue
        ccomp.append(ft['properties']['DN'])
        polygon_complete = shape(ft['geometry'])
        # only take the outer shell of the polygon in consideration
        poly.append(LinearRing(polygon_complete.exterior.coords))

    # Generate a list of line objects, with one list per connected component
    poly_unique = []
    ccomp_unique = []
    for ccomp_it in np.unique(ccomp):
        ccomp_unique.append(ccomp_it)
        # find all the polygons for this connected component
        #     and add them as a list
        index = np.where(np.array(ccomp)==ccomp_it)[0]
        polylist = []
        for it in index:
            polylist.append(poly[it])
        poly_unique.append(polylist)

    # track the projection and geotransform
    data =  gdal.Open(ccFile, gdal.GA_ReadOnly)
    geoTrans = data.GetGeoTransform()
    proj = data.GetProjection()
    cols = data.RasterXSize
    rows = data.RasterYSize
    sizeData = [cols,rows]
    data = None

    return poly_unique,ccomp_unique,geoTrans,proj,sizeData

def GDALread(filename,data_band=1,loadData=True):
    '''
        Script to load GDAL data
        Returns data, no-data value, projection, transformation
    '''

    # open the GDAL file and get typical data information
    try:
        data =  gdal.Open(filename, gdal.GA_ReadOnly)
    except:
        raise Exception(filename + " is not a GDAL supported")

    # loading the requested band or by default all
    raster = data.GetRasterBand(data_band)
    if loadData:
        out_data = raster.ReadAsArray()
    else:
        out_data = None
    # parsing no-data
    try:
        NoData = raster.GetNoDataValue()
    except:
        log.warning('Could not find a no-data value...')
        NoData = None

    # getting the gdal transform and projection
    geoTrans = data.GetGeoTransform()
    proj = data.GetProjection()

    return out_data,NoData,geoTrans,proj


def createConnComp_Int(inputs):
    """Function to generate intermediate connected component files and
    unwrapped VRT files that have with interger 2pi pixel shift applied.

    Will parse inputs in a single argument as it allows for parallel
    processing.

    Return a list of files in a unqiue temp folder.
    """
    # Parse the inputs to variables
    saveDir = inputs['saveDir']
    saveNameID = inputs['saveNameID']
    connFile = inputs['connFile']
    unwFile = inputs['unwFile']
    connCompMapping = inputs['connCompMapping']

    ## Generating the intermediate files
    ## STEP 1: set-up the mapping functions
    # Load the connected component
    connData, connNoData, connGeoTrans, connProj = GDALread(connFile)

    ## Define the mapping tables
    # Set up the connected component unique ID mapping
    connIDMapping = connCompMapping[:,1]

    # Add the no-data to the mapping as well (such that we can handle a
    # no-data region)
    # I.e. num comp + 1 = Nodata connected component which gets mapped
    # to no-data value again
    NoDataMapping = len(connIDMapping)
    connIDMapping = np.append(connIDMapping, [connNoData])

    # Set up the connected component integer 2PI shift mapping
    intMapping = connCompMapping[:,2]

    # Add the no-data to the mapping as well (such that we can handle a
    # no-data region)
    # I.e. max comp ID + 1 = Nodata connected component which gets
    # mapped to 0 integer shift such that the no-data region remains
    # unaffected
    intMapping = np.append(intMapping, [0])

    # Update the connected component with the new no-data value used
    # in the mapping
    connData[connData==connNoData] = NoDataMapping

    ## STEP 2: apply the mapping functions
    # Interger 2PI scaling mapping for unw phase
    intShift = intMapping[connData.astype('int')]

    # Connected component mapping to unique ID
    connData = connIDMapping[connData.astype('int')]

    ## STEP 3: write out the datasets
    # Write out the unqiue ID connected component file
    connDataName = os.path.abspath(os.path.join(saveDir, 'connComp',
            saveNameID+'_connComp.tif'))
    write_ambiguity(connData, connDataName, connProj, connGeoTrans, connNoData)

    # Write out the integer map as tiff file
    intShiftName = os.path.abspath(os.path.join(saveDir, 'unw',
            saveNameID+'_intShift.tif'))
    write_ambiguity(intShift, intShiftName, connProj, connGeoTrans)

    # Write out the scalled vrt => 2PI * integer map
    length = intShift.shape[0]
    width = intShift.shape[1]
    scaleVRTName = os.path.abspath(
            os.path.join(saveDir, 'unw', saveNameID+'_scale.vrt'))
    build2PiScaleVRT(scaleVRTName, intShiftName, length=length, width=width)

    # Offset the VRT for the range offset correction
    unwRangeOffsetVRTName = os.path.abspath(
        os.path.join(saveDir, 'unw', saveNameID + '_rangeOffset.vrt'))
    buildScaleOffsetVRT(
                        unwRangeOffsetVRTName,
                        unwFile,connProj,
                        connGeoTrans,
                        File1_offset = connCompMapping[1,3],
                        length = length,width=width)

    # writing out the corrected unw phase vrt => phase + 2PI * integer map
    unwVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID
                                 + '.vrt'))
    buildSumVRT(unwVRTName,
                unwRangeOffsetVRTName,
                scaleVRTName,
                connProj,
                connGeoTrans,
                length,
                width,
                description= inputs['description'])

    return [connDataName, unwVRTName]

def write_ambiguity(data, outName,proj, geoTrans, noData=False):
    """Write out an integer mapping in the Int16/Byte data range of values."""
    # GDAL precision support in tif
    Byte = gdal.GDT_Byte
    Int16 = gdal.GDT_Int16

    # Check if the path to the file needs to be created
    dirname = os.path.dirname(outName)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Get the GEOTIFF driver
    driver = gdal.GetDriverByName('GTIFF')
    # Leverage the compression option to ensure small file size
    dst_options = ['COMPRESS=LZW']
    # Create the dataset
    ds = driver.Create(outName , data.shape[1], data.shape[0], 1, Int16,
        dst_options)
    # Set the proj and transformation
    ds.SetGeoTransform(geoTrans)
    ds.SetProjection(proj)
    # Populate the first band with data
    bnd = ds.GetRasterBand(1)
    bnd.WriteArray(data)
    # Set the no-data value
    if noData is not None:
        bnd.SetNoDataValue(noData)
    bnd.FlushCache()
    # Close the file
    ds = None


def build2PiScaleVRT(output, File, width=False, length=False):
    """Build a VRT file which scales a GDAL byte file with 2PI."""
    # DBTODO: The datatype should be loaded by default from the
    #     source raster to be applied.
    # should be ok for now, but could be an issue for large connected com

    # The VRT template with 2PI scaling functionality
    vrttmpl = '''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <VRTRasterBand dataType="Float32" band="1">
    <ComplexSource>
        <SourceFilename relativeToVRT="1">{File}</SourceFilename>
        <SourceBand>1</SourceBand>
        <SourceProperties RasterXSize="{width}" RasterYSize="{length}" DataType="Int16"/>
        <ScaleOffset>0</ScaleOffset>
        <ScaleRatio>6.28319</ScaleRatio>
    </ComplexSource>
    </VRTRasterBand>
</VRTDataset>'''

    # The inputs needed to build the VRT
    # Load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # Check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Write out the VRT file
    with open(output, 'w') as fid:
        fid.write(vrttmpl.format(width = width, length = length, File = File))


def buildScaleOffsetVRT(output,File1,proj,geoTrans,File1_offset=0,
                        File1_scale = 1, width=False,length=False,
                        description='Scalled and offsetted VRT'):
    """Build VRT file summing 2 files together using pixel functionality."""
    # the vrt template with sum pixel functionality
    vrttmpl = '''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <VRTRasterBand dataType="Float32" band="1">
    <Description>{description}</Description>
        <ComplexSource>
            <SourceFilename relativeToVRT="1">{File1}</SourceFilename>
            <SourceBand>1</SourceBand>
            <ScaleOffset>{File1_offset}</ScaleOffset>
            <ScaleRatio>{File1_scale}</ScaleRatio>
            <SrcRect xOff="0" yOff="0" xSize="{width}" ySize="{length}"/>
            <DstRect xOff="0" yOff="0" xSize="{width}" ySize="{length}"/>
            <NoData>0</NoData>
        </ComplexSource>
    </VRTRasterBand>
</VRTDataset>'''
    #<SourceProperties RasterXSize="{lengthwidth}" RasterYSize="{length}" DataType="Float32"/>

    # The inputs needed to build the VRT
    # Load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File1, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # Check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Write out the VRT file
    with open('{0}'.format(output) , 'w') as fid:
        fid.write( vrttmpl.format(width = width,
                                  length = length,
                                  File1 = File1,
                                  File1_offset = File1_offset,
                                  File1_scale = File1_scale,
                                  proj = proj,
                                  geoTrans=str(geoTrans)[1:-1],
                                  description = description))

def buildSumVRT(output,File1,File2,proj,geoTrans,length=False,
                width=False,description='Unwrapped Phase'):
    """Building a VRT file which sums two files together using pixel
    functionality.
    """

    # The VRT template with sum pixel functionality
    vrttmpl = '''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <SRS>{proj}</SRS>
    <GeoTransform>{geoTrans}</GeoTransform>
    <VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">
    <Description>{description}</Description>
        <PixelFunctionType>sum</PixelFunctionType>
        <SimpleSource>
            <SourceFilename>{File1}</SourceFilename>
            <NoData>0</NoData>
        </SimpleSource>
        <SimpleSource>
            <SourceFilename>{File2}</SourceFilename>
            <NoData>0</NoData>
        </SimpleSource>
    </VRTRasterBand>
</VRTDataset>'''

    # The inputs needed to build the vrt
    # Load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File1, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # Check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Write out the VRT file
    with open('{0}'.format(output) , 'w') as fid:
        fid.write( vrttmpl.format(width = width, length = length,
                                  File1 = File1, File2 = File2,
                                  proj = proj, geoTrans = str(geoTrans)[1:-1],
                                  description = description))

def point2unwPhase(inputs):
    """Function to calculate the unwrapped phase value of a point
       and a region around it from only a single connected component
    Will parse inputs in a single argument as it allows
    for parallel processing.
    """

    # parsing the inputs to variables
    connCompID = int(inputs[0])
    connFile = str(inputs[1])
    unwFile = str(inputs[2])
    point = inputs[3].astype('float64')
    geoTrans = tuple(map(np.float64, inputs[4][1:-1].split(',')))
    region = int(inputs[5])

    # the region needs to be odd
    if region % 2 == 0:
        region = region+1

    # Will first convert the point to the closed grid point
    pixel_x = int((point[0] - geoTrans[0]) / geoTrans[1])
    pixel_y = int((point[1] - geoTrans[3]) / geoTrans[5])

    # setting up the origin of the chunk to read
    xoff = pixel_x-(region-1)/2
    yoff = pixel_y-(region-1)/2
    # updating to a valid part of the raster
    if xoff<0:
        xoff=0
    if yoff<0:
        yoff=0

    # loading a chunk of the connected component data
    ds_connFile = gdal.Open(connFile, gdal.GA_ReadOnly)
    connData = ds_connFile.ReadAsArray(xoff = xoff,
                                       yoff = yoff,
                                       xsize = region,
                                       ysize = region)
    ds_connFile = None

    # loading a chunk of the unwrapped phase data
    ds_unwFile = gdal.Open(unwFile, gdal.GA_ReadOnly)
    unwData = ds_unwFile.ReadAsArray(xoff = xoff,
                                     yoff = yoff,
                                     xsize = region,
                                     ysize = region)
    ds_unwFile = None

    # taking the average phase for all the pixels of the connected component
    try:
        unwPhase = np.mean(unwData[connData==connCompID])
        np.where(connData==connCompID)
    except:
        log.error('Looks like the region is not large enough,'
                  'cannot find your connected component: {}.'
                  'Entering debug mode.'.format(connCompID))
        print(connData)
        pdb.set_trace()

    # return back the phase value
    return unwPhase


def gdalTest(file, verbose=False):
    """Checks if the file is GDAL compatible and return
    compatible VRT in case it exists
    """
    # GDAL by default does not error out, so catch exception as error.
    # note warnings are not captured, e.g. warning in case netcdf file exist,
    #     but variable is not contained.
    gdal.UseExceptions()

    # DBTODO: Beyond GDAL compatible we should check if the file is GEO-CODED!
    #         e.g. use of proj, if not remove the file from list.
    #             Do not rely on .geo string in name

    file_success = None
    # check first if the file is GDAL compatible
    try:
        ds = gdal.Open(file, gdal.GA_ReadOnly)
        ds = None
    except:
        log.debug('%s is not GDAL compatible', file)
        return file_success

    # ideally use vrt file
    filepart1, filepart2 = os.path.splitext(file)
    if filepart2 != '.vrt' and filepart2 != '.hdr':
        filetest = file + ".vrt"
    else:
        filetest = filepart1 + ".vrt"

    # try if the vrt file can be loaded with GDAL.
    # if it does not exist or fails, then return original file else the vrt
    try:
        # will first check if file exist, as for a netcdf file one will modify
        #     the variable and gdal will give a warning not and error
        if os.path.isfile(filetest):
            ds = gdal.Open(filetest, gdal.GA_ReadOnly)
            ds = None
            log.debug("%s is GDAL compatible", filetest)
            return filetest
        else:
            log.debug("%s is GDAL compatible", file)
            return file
    except:
        log.debug("%s is GDAL compatible", file)
        return file



def product_stitch_overlap(unw_files, conn_files, prod_bbox_files, bbox_file,
    prods_TOTbbox, outFileUnw = './unwMerged',
    outFileConnComp = './connCompMerged', outputFormat='ENVI',
    mask=None, verbose=False):
    """Stitching of products minimizing overlap between products"""
    # report method to user
    # print('STITCH Settings: Product overlap approach')
    # print('Solver: Minimize overlap')

    # Hand over to product overlap stitch code
    unw = UnwrapOverlap()
    unw.setInpFile(unw_files)
    unw.setConnCompFile(conn_files)
    unw.setOutFileConnComp(outFileConnComp)
    unw.setOutFileUnw(outFileUnw)
    unw.setProdBBoxFile(prod_bbox_files)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setOutFileUnw(outFileUnw)
    unw.setOutFileConnComp(outFileConnComp)
    unw.setVerboseMode(verbose)
    unw.UnwrapOverlap()

def product_stitch_2stage(unw_files, conn_files, bbox_file, prods_TOTbbox,
                          unwrapper_2stage_name = None, solver_2stage = None,
                          outFileUnw = './unwMerged',
                          outFileConnComp = './connCompMerged',
                          outputFormat='ENVI',mask=None, verbose=False):
    """Stitching of products using the two-stage unwrapper approach
    i.e. minimize the discontinuities between connected components
    """
    # The solver used in minimizing the stiching of products
    if unwrapper_2stage_name is None:
        unwrapper_2stage_name = 'REDARC0'

    if solver_2stage is None:
        # If unwrapper_2state_name is MCF then solver is ignored
        # and relaxIV MCF solver is used by default
        solver_2stage = 'pulp'

    # report method to user
    log.info('STITCH Settings: Connected component approach')
    log.info('Name: %s', unwrapper_2stage_name)
    log.info('Solver: %s', solver_2stage)

    # Hand over to 2Stage unwrap
    unw = UnwrapComponents()
    unw._legacy_flag = True
    unw.setInpFile(unw_files)
    unw.setConnCompFile(conn_files)
    unw.setOutFileConnComp(outFileConnComp)
    unw.setOutFileUnw(outFileUnw)
    unw.setSolver(solver_2stage)
    unw.setRedArcs(unwrapper_2stage_name)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setOutFileUnw(outFileUnw)
    unw.setOutFileConnComp(outFileConnComp)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setVerboseMode(verbose)
    unw.unwrapComponents()
