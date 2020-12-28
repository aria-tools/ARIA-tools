#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
from scipy.ndimage import binary_dilation,binary_erosion,generate_binary_structure


from ARIAtools.logger import logger
from ARIAtools.shapefile_util import open_shapefile, save_shapefile

log = logging.getLogger(__name__)

solverTypes = ['pulp', 'glpk', 'gurobi']
redarcsTypes = {'MCF':-1, 'REDARC0':0, 'REDARC1':1, 'REDARC2':2}
stitchMethodTypes = ['overlap','2stage','morph']

class Stitching:
    '''
    This is the parent class of all stiching codes.
    It is whatever is shared between the different variants of stitching methods
    e.g. - setting of the main input arguments,
         - functions to verify GDAL compatibility,
         - function to write new connected component of merged product based on a mapping table,
         - fucntion to write unwrapped phase of merged product based on a mapping table
    '''

    def __init__(self):
        '''
        Setting the default arguments needed by the class.
        Parse the filenames and bbox as None as they need to be set by the user,
        which will be caught when running the child classes of the respective stitch method
        '''
        self.inpFile = None
        self.ccFile = None
        self.prodbboxFile = None
        self.prodbbox = None
        self.solver ='pulp'
        self.redArcs =-1
        self.mask = None
        self.outFileUnw = './unwMerged'
        self.outFileConnComp = './connCompMerged'
        self.outputFormat='ENVI'

        # stitching methods, by default leverage product overlap method
        # other options would be to leverage connected component
        self.setStitchMethod("overlap")


    def setInpFile(self, input):
        """ Set the input Filename for stitching/unwrapping """
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(input, np.str):
            input = [input]
        self.inpFile = input

        # the number of files that needs to be merged/ unwrapped
        self.nfiles = np.shape(self.inpFile)[0]

    def setOutFile(self, output):
        """ Set the output File name """
        self.outFile = output

    def setConnCompFile(self, connCompFile):
        """ Set the connected Component file """
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(connCompFile, np.str):
            connCompFile = [connCompFile]
        self.ccFile = connCompFile

    def setProdBBoxFile(self, ProdBBoxFile):
        """ Set the product bounding box file(s) """
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(ProdBBoxFile, np.str):
            ProdBBoxFile = [ProdBBoxFile]
        self.prodbboxFile = ProdBBoxFile

    def setBBoxFile(self, BBoxFile):
        """ Set bounds of bbox """
        self.bbox_file = BBoxFile

    def setTotProdBBoxFile(self, prods_TOTbbox):
        """ Set common track bbox file"""
        self.setTotProdBBoxFile = prods_TOTbbox

    def setStitchMethod(self,stitchMethodType):
        """ Set the stitch method to be used to handle parant class internals """
        if stitchMethodType not in stitchMethodTypes:
            raise ValueError(stitchMethodType + ' must be in ' + str(stitchMethodTypes))
        else:
            self.stitchMethodType =stitchMethodType

        if self.stitchMethodType=='2stage':
            self.description="Two-stage corrected/stiched Unwrapped Phase"
        elif self.stitchMethodType=='overlap':
            self.description = "Overlap-based stiched Unwrapped Phase"
        elif self.stitchMethodType=='morph':
            self.description = "Morphology corrected/stiched Unwrapped Phase"
                
    def setRedArcs(self, redArcs):
        """ Set the Redundant Arcs to use for LP unwrapping """
        self.redArcs = redArcs

    def setSolver(self, solver):
        """ Set the solver to use for unwrapping """
        self.solver = solver

    def setMask(self, mask):
        """ Set the mask file """
        self.mask = mask

    def setOutputFormat(self,outputFormat):
        """ Set the output format of the files to be generated """
        # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
        self.outputFormat = outputFormat
        if self.outputFormat=='VRT':
            self.outputFormat='ENVI'

    def setOutFileUnw(self,outFileUnw):
        """ Set the output file name for the unwrapped stiched file to be generated"""
        self.outFileUnw = outFileUnw

    def setOutFileConnComp(self,outFileConnComp):
        """ Set the output file name for the connected component stiched file to be generated"""
        self.outFileConnComp = outFileConnComp

    def setVerboseMode(self,verbose):
        """ Set verbose output mode"""
        logger.setLevel(logging.DEBUG)

    def __verifyInputs__(self):
        '''
            Verify if the unwrapped and connected component inputs are gdal compatible.
            That the provided shape files are well-formed.
            If not remove them from the list to be stiched.
            If a vrt exist and gdalcompatible update the file to be a vrt
        '''
        # track a list of files to keep
        inpFile_keep = []
        ccFile_keep = []
        prodbboxFile_keep = []
        bbox_keep = []
        for k_file in range(self.nfiles):
            # unw and corresponding conncomponent file
            inFile = self.inpFile[k_file]
            ccFile = self.ccFile[k_file]
            # shape file in case passed through
            if self.stitchMethodType == "overlap":
                prodbboxFile = self.prodbboxFile[k_file]

            # vrt files are prefered as they contain proj and transf information
            # Convert to inputs to vrt if exist, and check for gdal compatibility
            inFile_temp = gdalTest(inFile)
            ccFile_temp = gdalTest(ccFile)
            # check if it exist and well-formed and pass back shapefile of it
            if self.stitchMethodType == "overlap":
                bbox_temp = open_shapefile(prodbboxFile,'productBoundingBox',1)
                prodbboxFile_temp = prodbboxFile
            else:
                bbox_temp= "Pass"
                prodbboxFile_temp = "Pass"

            # if one of the two files fails then do not try and merge them later on as GDAL is leveraged
            if inFile_temp is None or ccFile_temp is None:
                log.info('Removing following pair combination (not GDAL compatible)')
                log.info('UNW: %s', inFile)
                log.info('CONN: %s', ccFile)
            elif prodbboxFile_temp is None:
                log.info('Removing following pair combination (malformed shapefile)')
                log.info('UNW: %s', inFile)
                log.info('CONN: %s', ccFile)
            else:
                inpFile_keep.append(inFile_temp)
                ccFile_keep.append(ccFile)
                # shape file in case passed through
                if self.stitchMethodType == "overlap":
                    bbox_keep.append(bbox_temp)
                    prodbboxFile_keep.append(prodbboxFile_temp)

        # update the input files and only keep those being GDAL compatible
        self.inpFile=inpFile_keep
        self.ccFile=ccFile_keep
        # shape file in case passed through
        if self.stitchMethodType == "overlap":
            self.prodbbox=bbox_keep
            self.prodbboxFile=prodbboxFile_keep

        # update the number of file in case some got removed
        self.nfiles = np.shape(self.inpFile)[0]

        if self.nfiles==0:
            log.info('No files left after GDAL compatibility check')
            sys.exit(0)
    def __createImages__(self):
        '''
            This function will write the final merged unw and conencted component file. As intermediate step tiff   files are generated with integer values which will represent the shift to be applied to connected componenet and the moduli shift to be applied to the unwrapped phase.
        '''

        ## Will first make intermediate files in a temp folder.
        # For each product there will be a connComp file and 3 files related unw files.
        # The connected component file will show the unique wrt to all merged files
        # For the unwrapped related files, there will be an integer offset tif file, a vrt file which scale this integer map by 2pi, and a vrt which combines the orginal unw phase file with the scaled map. The latter will be used for merging of the unwrapped phase.
        tempdir = tempfile.mkdtemp(prefix='IntermediateFiles_',dir='.')

        # will try multi-core version and default to for loop in case of failure
        try:
            # need to combine all inputs together as single argument tuple
            all_inputs = ()
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_" + str(counter)
                fileMappingDict['description'] = self.description
                # parse inputs as a tuple
                inputs = (fileMappingDict)
                # append all tuples in a single tuple
                all_inputs = all_inputs + (inputs,)
            # compute the phase value using multi-thread functionality
            intermediateFiles = Parallel(n_jobs=-1,max_nbytes=1e6)(delayed(createConnComp_Int)(ii) for ii in all_inputs)

        except:
            log.info('Multi-core version failed, will try single for loop')
            intermediateFiles = []
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_n" + str(counter)
                fileMappingDict['description'] = self.description
                # parse inputs as a tuple
                inputs = (fileMappingDict)
                # compute the phase value
                intermediateFiles_temp = createConnComp_Int(inputs)
                intermediateFiles.append(intermediateFiles_temp)

        # combining all conComp and unw files that need to be blended
        conCompFiles = []
        unwFiles = []
        for intermediateFile in intermediateFiles:
            conCompFiles.append(intermediateFile[0])
            unwFiles.append(intermediateFile[1])

        # check if the folder exist to which files are being generated.
        outPathUnw = os.path.dirname(os.path.abspath(self.outFileUnw))
        outPathConnComp = os.path.dirname(os.path.abspath(self.outFileConnComp))
        if not os.path.isdir(outPathUnw):
            os.makedirs(outPathUnw)
        if not os.path.isdir(outPathConnComp):
            os.makedirs(outPathConnComp)

        ## Will now merge the unwrapped and connected component files
        # remove existing output file(s)
        for file in glob.glob(self.outFileUnw + "*"):
            os.remove(file)
        gdal.BuildVRT(self.outFileUnw+'.vrt', unwFiles, options=gdal.BuildVRTOptions(srcNodata=0))
        gdal.Warp(self.outFileUnw, self.outFileUnw+'.vrt', options=gdal.WarpOptions(format=self.outputFormat, cutlineDSName=self.setTotProdBBoxFile, outputBounds=self.bbox_file))
        # Update VRT
        gdal.Translate(self.outFileUnw+'.vrt', self.outFileUnw, options=gdal.TranslateOptions(format="VRT"))
        # Apply mask (if specified).
        if self.mask is not None:
            update_file=gdal.Open(self.outFileUnw,gdal.GA_Update)
            update_file=update_file.GetRasterBand(1).WriteArray(self.mask.ReadAsArray()*gdal.Open(self.outFileUnw+'.vrt').ReadAsArray())
            update_file=None

        # remove existing output file(s)
        for file in glob.glob(self.outFileConnComp + "*"):
            os.remove(file)
        gdal.BuildVRT(self.outFileConnComp+'.vrt', conCompFiles, options=gdal.BuildVRTOptions(srcNodata=-1))
        gdal.Warp(self.outFileConnComp, self.outFileConnComp+'.vrt', options=gdal.WarpOptions(format=self.outputFormat, cutlineDSName=self.setTotProdBBoxFile, outputBounds=self.bbox_file))
        # Update VRT
        gdal.Translate(self.outFileConnComp+'.vrt', self.outFileConnComp, options=gdal.TranslateOptions(format="VRT"))
        # Apply mask (if specified).
        if self.mask is not None:
            update_file=gdal.Open(self.outFileConnComp,gdal.GA_Update)
            #mask value for conncomp must be set to nodata value -1
            ma_update_file=np.ma.masked_where(self.mask.ReadAsArray() == 0., gdal.Open(self.outFileConnComp+'.vrt').ReadAsArray())
            np.ma.set_fill_value(ma_update_file, update_file.GetRasterBand(1).GetNoDataValue())
            update_file=update_file.GetRasterBand(1).WriteArray(ma_update_file.filled())
            update_file=None

        cmd = "gdal_translate -of png -scale -ot Byte -q " + self.outFileUnw + ".vrt " + self.outFileUnw + ".png"
        os.system(cmd)

        # Remove the directory with intermediate files as they are no longer needed
        shutil.rmtree(tempdir)


class UnwrapOverlap(Stitching):
    '''
        Stiching/unwrapping using product overlap minimization
    '''

    def __init__(self):
        '''
            Inheret properties from the parent class
            Parse the filenames and bbox as None as they need to be set by the user, which will be caught when running the class
        '''
        Stitching.__init__(self)

    def UnwrapOverlap(self):

        ## setting the method
        self.setStitchMethod("overlap")

        ## check if required inputs are set
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
        # Update files to be vrt if they exist and remove files which failed the gdal compatibility
        self.__verifyInputs__()

        ## Calculating the number of phase cycles needed to miminize the residual between products
        self.__calculateCyclesOverlap__()

        ## Write out merged phase and connected component files
        self.__createImages__()


        return

    def __calculateCyclesOverlap__(self):
        '''Function that will calculate the number of cycles each component needs to be shifted in order to minimize the two-pi modulu residual between a neighboring component. Outputs a fileMappingDict with as key a file number. Within fileMappingDict with a integer phase shift value for each unique connected component.
        '''

        # only need to comptue the minimize the phase offset if the number of files is larger than 2
        if self.nfiles>1:

            # initiate the residuals and design matrix
            residualcycles = np.zeros((self.nfiles-1,1))
            residualrange = np.zeros((self.nfiles-1,1))
            A = np.zeros((self.nfiles-1,self.nfiles))

            # the files are already sorted in the ARIAproduct class, will make consecutive overlaps between these sorted products
            for counter in range(self.nfiles-1):
                # getting the two neighboring frames
                bbox_frame1 = self.prodbbox[counter]
                bbox_frame2 = self.prodbbox[counter+1]

                # determining the intersection between the two frames
                if not bbox_frame1.intersects(bbox_frame2):
                    log.error("Products do not overlap or were not provided in a contigious sorted list.")
                    raise Exception
                polyOverlap = bbox_frame1.intersection(bbox_frame2)

                # will save the geojson under a temp local filename
                tmfile = tempfile.NamedTemporaryFile(mode='w+b',suffix='.json', prefix='Overlap_', dir='.')
                outname = tmfile.name
                # will remove it as GDAL polygonize function cannot overwrite files
                tmfile.close()
                tmfile = None
                # saving the temp geojson
                save_shapefile(outname, polyOverlap, 'GeoJSON')

                # calculate the mean of the phase for each product in the overlap region alone
                # will first attempt to mask out connected component 0, and default to complete overlap if this fails.
                # Cropping the unwrapped phase and connected component to the overlap region alone, inhereting the no-data.
                # connected component
                out_data,connCompNoData1,geoTrans,proj = GDALread(self.ccFile[counter],data_band=1,loadData=False)
                out_data,connCompNoData2,geoTrans,proj = GDALread(self.ccFile[counter+1],data_band=1,loadData=False)
                connCompFile1 = gdal.Warp('', self.ccFile[counter], options=gdal.WarpOptions(format="MEM", cutlineDSName=outname, outputBounds=polyOverlap.bounds, dstNodata=connCompNoData1))
                connCompFile2 = gdal.Warp('', self.ccFile[counter+1], options=gdal.WarpOptions(format="MEM", cutlineDSName=outname, outputBounds=polyOverlap.bounds, dstNodata=connCompNoData2))


                # unwrapped phase
                out_data,unwNoData1,geoTrans,proj = GDALread(self.inpFile[counter],data_band=1,loadData=False)
                out_data,unwNoData2,geoTrans,proj = GDALread(self.inpFile[counter+1],data_band=1,loadData=False)
                unwFile1 = gdal.Warp('', self.inpFile[counter], options=gdal.WarpOptions(format="MEM", cutlineDSName=outname, outputBounds=polyOverlap.bounds, dstNodata=unwNoData1))
                unwFile2 = gdal.Warp('', self.inpFile[counter+1], options=gdal.WarpOptions(format="MEM", cutlineDSName=outname, outputBounds=polyOverlap.bounds, dstNodata=unwNoData2))


                # finding the component with the largest overlap
                connCompData1 =connCompFile1.GetRasterBand(1).ReadAsArray()
                connCompData1[(connCompData1==connCompNoData1) | (connCompData1==0)]=np.nan
                connCompData2 =connCompFile2.GetRasterBand(1).ReadAsArray()
                connCompData2[(connCompData2==connCompNoData2) | (connCompData2==0)]=np.nan
                connCompData2_temp = (connCompData2*100)
                temp = connCompData2_temp.astype(np.int)-connCompData1.astype(np.int)
                temp[(temp<0) | (temp>2000)]=0
                temp_count = collections.Counter(temp.flatten())
                maxKey = 0
                maxCount = 0
                for key, keyCount in temp_count.items():
                    if key!=0:
                        if keyCount>maxCount:
                            maxKey =key
                            maxCount=keyCount

                # if the max key count is 0, this means there is no good overlap region between products.
                # In that scenario default to different stitching approach.
                if maxKey!=0 and maxCount>75:
                    # masking the unwrapped phase and only use the largest overlapping connected component
                    unwData1 = unwFile1.GetRasterBand(1).ReadAsArray()
                    unwData1[(unwData1==unwNoData1) | (temp!=maxKey)]=np.nan
                    unwData2 =unwFile2.GetRasterBand(1).ReadAsArray()
                    unwData2[(unwData2==unwNoData2) | (temp!=maxKey)]=np.nan

                    # Calculation of the range correction
                    unwData1_wrapped = unwData1-np.round(unwData1/(2*np.pi))*(2*np.pi)
                    unwData2_wrapped =unwData2-np.round(unwData2/(2*np.pi))*(2*np.pi)
                    arr =unwData1_wrapped-unwData2_wrapped

                    # data is not fully decorrelated
                    arr = arr - np.round(arr/(2*np.pi))*2*np.pi
                    range_temp =  np.angle(np.nanmean(np.exp(1j*arr)))

                    # calculation of the number of 2 pi cycles accounting for range correction
                    cycles_temp = np.round((np.nanmean(unwData1-(unwData2+range_temp)))/(2*np.pi))

                else:
                    # account for the case that no-data was left, e.g. fully decorrelated
                    # in that scenario use all data and estimate from wrapped, histogram will be broader...
                    unwData1 = unwFile1.GetRasterBand(1).ReadAsArray()
                    unwData1[(unwData1==unwNoData1)]
                    unwData2 =unwFile2.GetRasterBand(1).ReadAsArray()
                    unwData2[(unwData2==unwNoData2)]
                    # Calculation of the range correction
                    unwData1_wrapped = unwData1-np.round(unwData1/(2*np.pi))*(2*np.pi)
                    unwData2_wrapped =unwData2-np.round(unwData2/(2*np.pi))*(2*np.pi)
                    arr =unwData1_wrapped-unwData2_wrapped
                    arr = arr - np.round(arr/(2*np.pi))*2*np.pi
                    range_temp =  np.angle(np.nanmean(np.exp(1j*arr)))

                    # data is decorelated assume no 2-pi cycles
                    cycles_temp = 0


                # closing the files
                unwFile1 = None
                unwFile2 = None
                connCompFile1 = None
                connCompFile2 = None

                # remove the tempfile
                shutil.os.remove(outname)

                # store the residual and populate the design matrix
                residualcycles[counter]=cycles_temp
                residualrange[counter]=range_temp
                A[counter,counter]=1
                A[counter,counter+1]=-1


            # invert the offsets with respect to the first product
            cycles = np.round(np.linalg.lstsq(A[:,1:], residualcycles,rcond=None)[0])
            rangesoffset = np.linalg.lstsq(A[:,1:], residualrange,rcond=None)[0]
            #pdb.set_trace()

            # force first product to have 0 as offset
            cycles = -1*np.concatenate((np.zeros((1,1)), cycles), axis=0)
            rangesoffset = -1*np.concatenate((np.zeros((1,1)), rangesoffset), axis=0)

        else:
            # nothing to be done, i.e. no phase cycles to be added
            cycles = np.zeros((1,1))
            rangesoffset = np.zeros((1,1))

        # build the mapping dictionary
        fileMappingDict = {}
        connCompOffset =0

        for fileCounter in range(self.nfiles):

            # get the number of connected components
            n_comp = 20

            # The original connected components
            connComp = np.array(range(0,n_comp+1))
            connComp = connComp.reshape((n_comp+1,1))

            # The cyclemapping of product using a single offset (i.e. all connectedComponents have the same mapping)
            cycleMapping = np.ones((n_comp+1,1))*cycles[fileCounter,0]
            rangeOffsetMapping = np.ones((n_comp+1,1))*rangesoffset[fileCounter,0]

            # Generate the mapping of connectedComponents
            # Increment based on unique components of the merged product, 0 is grouped for all products
            connCompMapping = connComp
            connCompMapping[1:]=connComp[1:]+connCompOffset

            # concatenate and generate a connComponent based mapping matrix
            # [original comp, new component, 2pi unw offset]
            connCompMapping = np.concatenate((connComp, connCompMapping,cycleMapping,rangeOffsetMapping), axis=1)

            # Increment the count of total number of unique components
            connCompOffset = connCompOffset+n_comp

            # populate mapping dictionary for each product
            fileMappingDict_temp = {}
            fileMappingDict_temp['connCompMapping'] = connCompMapping
            fileMappingDict_temp['connFile'] =  self.ccFile[fileCounter]
            fileMappingDict_temp['unwFile'] = self.inpFile[fileCounter]
            # store it in the general mapping dictionary
            fileMappingDict[fileCounter]=fileMappingDict_temp

        # pass the fileMapping back into self
        self.fileMappingDict = fileMappingDict


class UnwrapMorph(Stitching):
    '''
        Stiching/unwrapping using product overlap minimization
    '''
    
    
    def __init__(self):
        '''
            Inheret properties from the parent class
            Parse the filenames and bbox as None as they need to be set by the user, which will be caught when running the class
            '''
        Stitching.__init__(self)
    

    
    def UnwrapMorph(self):
        
        ## setting the method
        self.setStitchMethod("morph")
    
        #processing parameters
        self.proc_params = {'ni':8,
                        'nb':10,
                        'minsize':10,
                        'debug':False,
                        'max_guard':4}

        ## Verify if all the inputs are well-formed/GDAL compatible
        # Update files to be vrt if they exist and remove files which failed the gdal compatibility
        self.__verifyInputs__()
    
        #import pdb
        #pdb.set_trace()
        
        ## run Morphology minization between connected components
        runw_ph,ncc = self.__fix_unwrap__()
        #runw_ph,ncc = self.fix_unwrap(runw_ph,ncc,ni=ni,nb=nb,minsize=minsize,debug=debug,max_guard=max_guard)



    def __fix_unwrap__(self):
        ni = self.proc_params['ni']
        nb = self.proc_params['nb']
        minsize = self.proc_params['minsize']
        debug = self.proc_params['debug']
        max_guard = self.proc_params['max_guard']
        

        
        # loading the data
        connData,connNoData,connGeoTrans,connProj = GDALread(self.ccFile[0])
        unwData,unwNoData,unwGeoTrans,unwProj = GDALread(self.inpFile[0])

        ncc=connData
        runw_ph=unwData
        
        max_cc    = np.ptp(ncc)
        numcc = np.unique(ncc[ncc > 0]).size
        guard = 0
        while True:
            cycles1,ord_overlaps1,nover1 = compute_cycles(runw_ph,ncc,nbest=nb,niter=ni,minsize=minsize,debug=debug)
            runw_ph = add_cycles(runw_ph,ncc,cycles1,ord_overlaps1)
            linked,singles = get_linked_cc(ord_overlaps1,max_cc)
            ncc = adjust_cc(ncc,linked)
            numcc_now = np.unique(ncc[ncc > 0]).size
            #if only one good component left break
            if numcc_now == 1:
                break
            #if no improvement break
            elif numcc_now == numcc:
                break
            #safe guard just in case
            elif guard == max_guard:
                break
            numcc = numcc_now
            guard += 1
        runw_ph = adjust_to_zero_cc(runw_ph,ncc,struc=None,niter=2,minsize=20)
       
        
        write_ambiguity(ncc, self.outFileConnComp,connProj, connGeoTrans,noData=connNoData)
        write_data(runw_ph, self.outFileUnw,unwProj, unwGeoTrans,noData=unwNoData)
        cmd = "gdal_translate -of png -scale -ot Byte -q " + self.outFileConnComp + " " + self.outFileConnComp + ".png"
        os.system(cmd)
        cmd = "gdal_translate -of png -scale -ot Byte -q " + self.outFileUnw + " " + self.outFileUnw + ".png"
        os.system(cmd)
        
        return runw_ph,ncc

    '''
    def get_all_overlaps(cc,scc,nbest=3,struct=None,niter=2,minsize=20):
        """
            Given the array of ccomp sorted by size get all possible overlaps using the first "nbest"
            inputs:
            cc = ccomp file.
            scc = array with all of sorted ccomp according to size.
            nbests = use the first largest nbest cc in scc to get the overlap with others cc.
            struct = morphological structure.
            niter = number iterations of morphological operations.
            minsize = minimum number of points in the overlaps. If less than minsize then no overlap.
            outputs:
            overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            complete = bool True if all the cc have overlap
            dones: the list of cc that have overlap
        """
        
        #keep trak of the cc that have not been adjusted. initially are all except for -1 and 0
        check_complete = np.unique(cc)
        check_complete = np.setdiff1d(check_complete,np.array([-1,0]),True)
        overlaps = {}
        #get all the pairs that overlap
        overlaps0,complete0,dones0 = get_overlaps(cc,scc[:nbest],overlaps,check_complete,struct,niter,minsize)
        #get the new list of cc minus -1,0
        check_complete = np.unique(cc)
        check_complete = np.setdiff1d(check_complete,np.array([-1,0]),True)
        #get the ones that still need to be done
        not_dones = np.setdiff1d(check_complete,dones0)
        #redo the overlap searching for the not_dones instead of nbest
        overlaps,complete,dones = get_overlaps(cc,not_dones,overlaps,not_dones,struct,niter,minsize)
        dones = np.unique(np.array(np.append(dones0,dones)))
        #check if all have been fixed. sometime is not possible if there the cc is isolated
        if len(dones) == len(check_complete):
            complete = True
        else:
            complete = False
        return overlaps,complete,dones

    def order_overlaps(pairs):
        """
            Create an OrderedDict with key the cc and value a list of cc that have overlap with the key.
            The order is based on the size of the values i.e. how many cc overlap with the key.
            inputs:
            pairs = list of cc pairs that overlap. The form is i_j where i,j are the cc.
            output:
            ddict = OrderedDict with key a cc and value an array of the cc overlapping with key.
        """
        
        ddict = collections.defaultdict(list)
        for p in pairs:
            keys = p.split('_')
            ddict[int(keys[0])].append(int(keys[1]))
        sizes = collections.defaultdict(int)
        for k,v in ddict.items():
            sizes[k] = len(v)
        ssizes = collections.OrderedDict(sorted(sizes.items(),key = lambda t:t[1],reverse=True))
        ddict = collections.OrderedDict([(k,ddict[k]) for k in ssizes.keys()])
        return ddict

    def get_resized_imgs(imgs,mask,buf=50):
        """
            Resize each image in images only where mask is non zero. Add a +- buf
            inputs:
            imgs = images to resize.
            mask = the mask to find the extremes.
            buf = add buf to the extremes.
            outputs:
            ret_imgs = the resized images
            min_x = the minimun x location in the original image
            min_y = the minimum y location in the original image
        """
        
        #find the extremes where mask is to reduce the computation
        locs = np.where(mask)
        min_y = max(locs[0].min() - buf,0)
        max_y = max(locs[0].max() + buf,mask.shape[0])
        min_x = max(locs[1].min() - buf,0)
        max_x = max(locs[1].max() + buf,mask.shape[1])
        ret_imgs = []
        for im in imgs:
            ret_imgs.append(im[min_y:max_y,min_x:max_x])
        return ret_imgs,min_x,min_y

    def get_overlaps(cc,scc,overlaps,check_complete,struct=None,niter=2,minsize=20):
        """
            Given the array of ccomp sorted by size get all possible overlaps using the first "nbests"
            inputs:
            cc = ccomp file.
            scc = array with subset of sorted ccomp according to size.
            overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            check_complete = list of the cc that have no overlap.
            struct = morphological structure.
            niter = number iterations of morphological operations.
            minsize = minimum number of points in the overlaps
            outputs:
            overlaps = the updated input 'overlaps'.
            complete = bool True if all the cc have overlap.
            dones = the list of cc that have overlap.
        """
        
        if struct is None:
            struct = generate_binary_structure(2, 1)

        dones = []
        complete = False
        for cc_now in scc:
            mask = cc == cc_now
            #resize cc and mask only around the extremes of mask
            [ncc,mask],min_x,min_y = get_resized_imgs([cc,mask],mask,buf=50)
            #get what's not cc_now
            nmask = np.logical_not(mask)
            #dilate mask
            dmask = binary_dilation(mask,struct,niter)
            #get intersection between dilate cc and not cc
            inters = np.logical_and(dmask,nmask)
            #select all the cc that are in the dilated part only
            sel_now = np.unique(ncc[inters])
            for i in sel_now:
                if i == -1:
                    continue
                if i == 0:
                    continue
                okey = str(i) + '_' + str(cc_now)
                #if already present skip
                if okey in overlaps:
                    continue
                mask1 = ncc == i
                [mask1,mmask,ddmask],min_xx,min_yy = get_resized_imgs([mask1,mask,dmask],mask1)
                dmask1 = binary_dilation(mask1,struct,niter)
                inters10 = np.logical_and(ddmask,mask1)
                inters01 = np.logical_and(dmask1,mmask)
                if np.count_nonzero(inters01) >= minsize and np.count_nonzero(inters10) >= minsize:
                    #if enough points overlap that save the ovelapping points.
                    #since we cropped the images need to put back the correct indices
                    ninters01 = np.zeros_like(cc,np.bool)
                    ninters01[min_y + min_yy:min_y + min_yy + inters01.shape[0],min_x + min_xx:min_x + min_xx + inters01.shape[1]] = inters01
                    ninters10 = np.zeros_like(cc,np.bool)
                    ninters10[min_y + min_yy:min_y + min_yy + inters10.shape[0],min_x + min_xx:min_x + min_xx + inters10.shape[1]] = inters10
                    overlaps[str(cc_now) + '_' + str(i)] = [ninters01,ninters10]
                    #mark the 2 cc as done
                    dones.append(cc_now)
                    dones.append(i)
                    if len(check_complete) == len(np.intersect1d(check_complete,np.unique(np.array(dones)))):
                        complete = True
                if complete == True:
                    break
            if complete == True:
                break
        return overlaps,complete,np.unique(np.array(dones))

    def get_overlaps_obs(cc,scc,overlaps,check_complete,struct=None,niter=2,minsize=20):
        """
            This method is deprecated (obsolete)
            Given the array of ccomp sorted by size get all possible overlaps using the first "nbests"
            inputs:
            cc = ccomp file.
            scc = array with subset of sorted ccomp according to size.
            overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            check_complete = list of the cc that have no overlap.
            struct = morphological structure.
            niter = number iterations of morphological operations.
            minsize = minimum number of points in the overlaps
            outputs:
            overlaps = the updated input 'overlaps'.
            complete = bool True if all the cc have overlap.
            dones = the list of cc that have overlap.
        """
        
        if struct is None:
            struct = generate_binary_structure(2, 1)
        
        dones = []
        complete = False
        for cc_now in scc:
            mask = cc == cc_now
            nmask = np.logical_not(mask)
            dmask = binary_dilation(mask,struct,niter)
            inters = np.logical_and(dmask,nmask)
            sel_now = np.unique(cc[inters])
            for i in sel_now:
                if i == -1:
                    continue
                if i == 0:
                    continue
                okey = str(i) + '_' + str(cc_now)
                #if already present skip
                if okey in overlaps:
                    continue
                mask1 = cc == i
                dmask1 = binary_dilation(mask1,struct,niter)
                inters10 = np.logical_and(dmask,mask1)
                inters01 = np.logical_and(dmask1,mask)
                if np.nonzero(inters01)[0].size >= minsize and np.nonzero(inters10)[0].size >= minsize:
                    overlaps[str(cc_now) + '_' + str(i)] = [inters01,inters10]
                    dones.append(cc_now)
                    dones.append(i)
                    if len(check_complete) == len(np.intersect1d(check_complete,np.unique(np.array(dones)))):
                        complete = True
                if complete == True:
                    break
            if complete == True:
                break
        return overlaps,complete,np.unique(np.array(dones))

    def swap_single(ord_overlaps,overlaps):
        """
            Swap the key and value (when there is only one value) of the ordered overlaps when the key is not part
            of any other pair of cc overlapping but the value is so that the pair is connected to others.
            inputs:
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            outputs:
            nord_overlaps = updated ord_overlaps.
            noverlaps = updates overlaps.
        """
        
        #when one cc is inside the other the one inside has more overlap because of the dilation, but it's not
        #connected to anything. So check if for single element values either the key or the value if part
        #of some other cc set. If the value is and not the key then swap since most likely the key is inside the value.
        #get_ipython().run_line_magic('debug', '')
        n_ord_over = {}
        nover = {}
        to_delete = []
        if len(ord_overlaps) > 1:
            for k,v in ord_overlaps.items():
                v = copy.deepcopy(v)
                if len(v) == 1:
                    for k1,v1 in ord_overlaps.items():
                        if k1 == k:
                            continue
                        if v[0] in v1 and k not in v1:#the second check is to see if both v and k are in v1
                            #if this is the case leave it alone otherwise remove it
                            n_ord_over[v[0]] = [k]
                            nover[str(v[0]) + '_' + str(k)] =  overlaps[str(k) + '_' + str(v[0])][::-1]#swap order
                            to_delete.append(str(k) + '_' + str(v[0]))
                        elif v[0] in v1 and k in v1:
                            continue
                        else:
                            n_ord_over[k] = v
                            nover[str(k) + '_' + str(v[0])] = overlaps[str(k) + '_' + str(v[0])]
                else:
                    n_ord_over[k] = v
                    for vi in v:
                        nover[str(k) + '_' + str(vi)] = overlaps[str(k) + '_' + str(vi)]
            
            #remove the one swapped
            to_delete = np.unique(np.array(to_delete))
            for k in to_delete:
                if k in nover:
                    del nover[k]
                    k0 = int(k.split('_')[0])
                    k1 = int(k.split('_')[1])
                    n_ord_over[k0].remove(k1)
                    if len(n_ord_over[k0]) == 0:
                        n_ord_over.pop(k0, None)
        else:#fringe case where there is only one left and the above logig doesn't work
            n_ord_over = ord_overlaps
            nover = overlaps

        return collections.OrderedDict(sorted(n_ord_over.items(),key = lambda t:t[0])),nover

    def remove_dups(ord_overlaps,nover):
        """
            Remove cc that are in multiple places. Keep the instance that has larger overlap.
            inputs:
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            nover = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            outputs:
            nord_overlaps = updated ord_overlaps.
            nover = updates nover.
        """
        
        keys = list(ord_overlaps.keys())
        for i in range(len(keys) - 1):
            k = keys[i]
            for j in range(i + 1,len(keys)):
                k1 = keys[j]
                common = np.intersect1d(ord_overlaps[k],ord_overlaps[k1])
                for c in common:
                    key = str(k) + '_' + str(c)
                    key1 = str(k1) + '_' + str(c)
                    if len(nover[key][0]) > len(nover[key1][0]):
                        ord_overlaps[k1].remove(c)
                        del nover[key1]
                    else:
                        ord_overlaps[k].remove(c)
                        del nover[key]
        keys = list(ord_overlaps.keys())
        for k in keys:
            v = ord_overlaps[k]
            if len(v) == 0:
                del ord_overlaps[k]
        return ord_overlaps,nover

    def reorder_overlaps(ord_overlaps):
        """
            Order the inputs so that when cycles are added to a cc no prior modified cc is affected.
            inputs:
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            outputs:
            nord_overlaps = the updated ord_overlaps.
        """
        
        #make sure that the ovelarps are sorted so that when one is modified no other connected to it
        #get modified later one. For instance 10:[11] and 13:[10,20]. In this order 11 is modified w.r.t. 10 but then
        #10 is modified w.r.t 13 and 11 now is off
        keys = list(ord_overlaps.keys())
        ret = []
        for i,k in enumerate(keys):
            swapped = False
            for k1 in keys[i+1:]:
                if k in ord_overlaps[k1]:
                    #NOTE:k cannot be present in any other list because we remove dups
                    #     so it's safe to swap and break
                    if k1 in ret:#already present from previous swap
                        ret.append(k)
                    else:
                        ret.extend([k1,k])
                    swapped = True
                    break
            if not swapped:#there was now connection with any other cc so just add
                ret.append(k)
        return collections.OrderedDict([(k,ord_overlaps[k]) for k in ret])

    def get_cycles(unw_ph,ord_overlaps,nover):
        """
            Get the number of cycles that need to be added to a cc based on the difference of the phase in the
            overlapping region. Only the values in the ord_overlap dict are updated while the key has 0 cycles (unless
            a key is a value for another key in the dict)
            inputs:
            unw_ph: original unwrapped phase.
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            nover = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
            outputs:
            cycles = dictionary with key a cc and value the number of cycles to be added.
            stats = a dictionary with key the overlapping pair 'i_j' and value an array with the cycles,
            stds of the phase in the two overlapping regions (normalized to 2pi), the number of elementes
            in the overlapping regions.
        """
        
        cycles = collections.defaultdict(int)
        stats = collections.defaultdict(list)
        first_time = True
        for k,v in ord_overlaps.items():
            if first_time:
                cycles[k] = 0
                first_time = False
            for k1 in v:
                key = str(k) + '_' + str(k1)
                inters,inters1 = nover[key]
                cyc = round((unw_ph[inters].mean() - unw_ph[inters1].mean())/(2*np.pi))
                cycles[k1] =  cyc
                stats[key] = [cyc,unw_ph[inters].std()/(2*np.pi),unw_ph[inters1].std()/(2*np.pi),
                              np.count_nonzero(inters),np.count_nonzero(inters1)]
        return cycles,stats

    def alter_cycles(unw_ph,cc,selcc,ncyc):
        """
            Alterate the unwrapped phase by adding an interger times of 2pi to a given cc.
            inputs:
            unw_ph = original unwrapped phase.
            cc = the connected component image.
            selcc = array with the cc that need to be modified.
            ncyc = array with the cycles to add to the cc in selcc.
            outputs:
            nunw_ph = the modified unwrapped phase.
        """
        
        nunw_ph = unw_ph.copy()
        for i,j in zip(selcc,ncyc):
            seli = cc == i
            nunw_ph[seli] += 2*np.pi*j
        return nunw_ph

    def get_images(dir_now,unw_file,cc_file):
        """
            Load the cc and unwrapped phase images into a numpy array.
            inputs:
            dir_now = directory containing the image files. Needs the xml metadata and the binary.
            outputs:
            cc: 2D int8 numpy array with the connected component image values.
            unw_ph: 2D float32 numpy array with the unwrapped phase image values.
        """
       
        data = gdal.Open(unw_file,gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        unw_ph = data_band.ReadAsArray()
        data = None
        # doing the connected component
        data = gdal.Open(cc_file,gdal.GA_ReadOnly)
        data_band = data.GetRasterBand(1)
        cc = data_band.ReadAsArray()
        data = None
        img = None
        return cc,unw_ph,img

    def compute_cycles(unw_ph,cc,nbest=3,niter=2,minsize=20,debug=True):
        """
            Compute the cycles to be added to the unpwrapped phase.
            inputs:
            unw_ph = original unwrapped phase.
            cc = the connected component image.
            nbests = use the first largest nbest cc in scc to get the overlap with others cc.
            niter = number iterations of morphological operations.
            minsize = minimum number of points in the overlaps.
            debug = bool if True prints intermediate steps results.
            output:
            cycles = dictionary with key a cc and value the number of cycles to be added.
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            nover = dictionary with key = i_j and values the two arrays which are the overlaps after
            dilating one cc (i or j) and not the other (j or i).
        """
        
        struct = generate_binary_structure(2, 1) # to identify neighbors
        cc_sizes = []
        
        
        # remove both the no-data and connected component 0 from the list under consideration
        cc_left = np.setdiff1d(np.unique(cc),[-1,0],True)
        for c in cc_left:
            cc_sizes.append(np.count_nonzero(cc == c))
        cc_sizes = np.array(cc_sizes,np.int32)
        scc_sizes = cc_left[np.argsort(cc_sizes)[::-1]]

        overlaps,complete,dones = get_all_overlaps(cc,scc_sizes,nbest=nbest,struct=struct,niter=niter,minsize=minsize)
        #get_all_overlaps(cc,scc_sizes,2)
        tord_overlaps = order_overlaps(overlaps.keys())
        if debug:
            print('order_overlaps',tord_overlaps)
        ord_overlaps,nover = swap_single(tord_overlaps,overlaps)
        if debug:
            print('swap_single',ord_overlaps,'\n',nover.keys())
        ord_overlaps,nover = remove_dups(ord_overlaps,nover)
        if debug:
            print('remove_dups',ord_overlaps,'\n',nover.keys())
        ord_overlaps = reorder_overlaps(ord_overlaps)
        if debug:
            print('reorder_overlaps',ord_overlaps)
        cycles,stats = get_cycles(unw_ph,ord_overlaps,nover)
        return cycles,ord_overlaps,nover

    def add_cycles(unw_ph,cc,cycles,ord_overlaps):
        """
            Add cycles to the unwrapped image to remove discontinuity.
            inputs:
            unw_ph = original unwrapped phase.
            cc = the connected component image.
            cycles = dictionary with key a cc and value the number of cycles to be added.
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            outputs:
            nunw_ph = the fixed unwrapped phase.
        """
        
        cp_cycles = copy.deepcopy(cycles)
        keys = list(ord_overlaps.keys())
        #go one key at the time except for the first since it's the zero shift
        for i in range(1,len(keys)):
            key = keys[i]
            #go back to all the previous main cc and see if it belongs to one of the overlaps. if so need
            #to adjust all the cycles for this cc
            for j in range(i):
                key1 = keys[j]
                if key in ord_overlaps[key1]:#cc overlaps with previous one.
                    for v in ord_overlaps[key]:
                        #print(key1,key,cp_cycles[key],v,cp_cycles[v])
                        cp_cycles[v] += cp_cycles[key]
        to_change = []
        max_cc    = np.ptp(cc)
        for i in range(max_cc):
            if i in cp_cycles:
                to_change.append(cp_cycles[i])
            else:
                to_change.append(0)
        to_change = np.array(to_change)
        return unw_ph + 2*np.pi*to_change[cc]

    def get_linked_cc(ord_overlaps,max_cc):
        """
            Find all the cc that linked to each other.
            inputs:
            ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
            max_cc = maximum number of connected components in the cc image.
            outputs:
            linked = list of all the cc linked to each other.
        """
        
        linked = []
        keys = list(ord_overlaps.keys())
        done = []
        for i,(k,v) in enumerate(ord_overlaps.items()):
            #if k and it's v are already linked than skip them
            if k in done:
                continue
            done.append(k)
            link_now = [k] + v
            #go through the remaining keys and if k1 was part of the value of another key
            #it means that they are linked. ex. k,v = (1, [15, 16, 17, 18, 13, 19, 12]) and k1,v1 = (13, [10, 20])
            #at the end they should all belong to the same cc.
            for k1 in keys[i+1:]:
                if k1 in link_now:
                    link_now += ord_overlaps[k1]
                    done.append(k1)
            linked.append(link_now)
        not_singles = np.array([i for j in linked for i in j])
        singles = np.setdiff1d(np.arange(1,max_cc),not_singles)
        singles.sort()
        return linked,singles

    def adjust_cc(cc,linked):
        """
            Adjust the connected component values of the linked cc.
            inputs:
            cc = the connected component image.
            linked = list of all the cc linked to each other.
            outputs:
            ncc = adjusted connected component image.
            """
        ncc = cc.copy()
        for l in linked:
            n_value = l[0]
            for l1 in l[1:]:
                sel = ncc == l1
                ncc[sel] = n_value
        return ncc

    def adjust_to_zero_cc(unw_ph,cc,struc=None,niter=2,minsize=20):
        """
            Adjust non zero cc to the zero one. Used to make the phase image look smooth.
            inputs:
            unw_ph = original unwrapped phase.
            cc = the connected component image.
            struct = morphological structure.
            niter = number iterations of morphological operations.
            minsize = minimum number of points in the overlaps.
            outputs:
            runw_ph = modified pahse image.
        """
        
        check_complete = np.unique(cc)[1:]
        overlaps = {}
        overlaps,complete0,dones0 = get_overlaps(cc,[0],overlaps,check_complete,struct=struc,niter=niter,minsize=minsize)
        ord_overlaps = {0:[int(i.split('_')[1]) for i in overlaps.keys()]}
        cycles,stats = get_cycles(unw_ph,ord_overlaps,overlaps)
        runw_ph = add_cycles(unw_ph,cc,cycles,ord_overlaps)
        return runw_ph
    '''



class UnwrapComponents(Stitching):
    '''
        Stiching/unwrapping using 2-Stage Phase Unwrapping
    '''

    def __init__(self):
        '''
            Inheret properties from the parent class
            Parse the filenames and bbox as None as they need to be set by the user, which will be caught when running the class
        '''
        Stitching.__init__(self)


    def UnwrapComponents(self):

        ## setting the method
        self.setStitchMethod("2stage")
        self.region=5

        if self.inpFile is None:
            log.error("Input unwrapped file(s) is (are) not set.")
            raise Exception

        if self.ccFile is None:
            log.error("Input Connected Components file(s) is (are) not set.")
            raise Exception

        if self.solver not in solverTypes:
            raise ValueError(self.treeType + ' must be in ' + str(unwTreeTypes))

        if self.redArcs not in redarcsTypes.keys():
            raise ValueError(self.redArcs + ' must be in ' + str(redarcsTypes))


        ## Verify if all the inputs are well-formed/GDAL compatible
        # Update files to be vrt if they exist and remove files which failed the gdal compatibility
        self.__verifyInputs__()

        ## Loop over the number of valid files and populate a combined dictionary with keys being polygons around connected component edges in "polyTableDict" variable.
        #  Key is a unique polygon [1, max number of unique polygons of all products]
        #  Note: there can be components with mutiple polygons (e.g. inner/outer edges)
        #  Each key will have a corresponding dictionary with keys:
        #       "connFile" = Original connected component file
        #       "connCompID" = Original ID connected component id (wrt original connected component file)
        #       "unwFile" = Original unwrapped file
        self.__populatePolyTable__()

        ## Populate a table where each unique polygon is matches with all other (not same component) polygons, and keep in the table the points that represent the closest two points between these polygons.
        # The table is populated with columns as:
        #       Source Unique CC ID, Source CC ID, Source CC FILE, Source UNW FILE, Source X-geocoordinate, Source Y-geocoordinate, Source geoTrans, Partner Unique CC, Dist to partner, Source unw phase, Source X-coordinate, Source Y-coordinate
        self.__populatePointTable__()

        ## Calculating the number of phase cycles needed to miminize the residual between connected components. Minimization is only done between the tablePoints (i.e. X-coordinate, Y-coordinate, unw phase).
        self.__calculateCycles__()

        ## Write out merged phase and connected component files
        self.__createImages__()


        return

    def __populatePolyTable__(self):
        '''
            Build a dictionary with connected component polygon information
            columns: key is a unique polygon [1, max number of unique polygons of all products]
            each key will have a corresponding dictionary with keys:
                "connCompFile" = Original connected component file
                "conmCompID" = Original ID connected component id (wrt original connected component file)
                "unwFile" = Original unwrapped file
        '''

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
        '''
            Function that populates a table of points which maps two closest points between two connected components. Loops over all possible connected component pairings.
            The table is populated with columns as: Source Unique CC ID, Source CC ID, Source CC FILE, Source UNW FILE, Source X-geocoordinate, Source Y-geocoordinate, Source geoTrans, Partner Unique CC, Dist to partner, Source unw phase, Source X-coordinate, Source Y-coordinate
        '''

        # only proceed if there are polygons that needs mapping to a closest polygon
        if len(self.polyTableDict) <=1:
            raise Exception("Nothing to do as all are from same connected component")

        # populate table  of all polygons
        polygon_pair_counter = 0
        # loop over each polygon and compute the distance to all polygons.
        #print('Calculating the phase values for polygon pairs')
        for poly_counter_1 in range(0,len(self.polyTableDict)):
            # polygon 1
            poly1 = self.polyTableDict[poly_counter_1]['poly']
            connCompID1 = self.polyTableDict[poly_counter_1]['connCompID']
            connCompID_unique1 = self.polyTableDict[poly_counter_1]['connCompID_unique']
            connFile1 = self.polyTableDict[poly_counter_1]['connFile']
            unwFile1 = self.polyTableDict[poly_counter_1]['unwFile']
            geoTrans1 = self.polyTableDict[poly_counter_1]['geoTrans']
            sizeData1 = self.polyTableDict[poly_counter_1]['sizeData']

            # loop over all poygons and keep track of min distance polygon not of own component
            for poly_counter_2 in range(poly_counter_1,len(self.polyTableDict)):
                # polygon 2
                poly2 = self.polyTableDict[poly_counter_2]['poly']
                connCompID2 = self.polyTableDict[poly_counter_2]['connCompID']
                connCompID_unique2 = self.polyTableDict[poly_counter_2]['connCompID_unique']
                connFile2 = self.polyTableDict[poly_counter_2]['connFile']
                unwFile2 = self.polyTableDict[poly_counter_2]['unwFile']
                geoTrans2 = self.polyTableDict[poly_counter_2]['geoTrans']
                sizeData2 = self.polyTableDict[poly_counter_2]['sizeData']

                # will skip the pairing of the two polygons under the following two conditions:
                # 1) when mapping two of the same polyons, i.e. i==j or a diagonal element of the matrix of combinations
                # 2) when the two polygons are from the same connected component and originate from the same connected comp file
                if poly_counter_1==poly_counter_2 and connCompID1==connCompID2 and connFile1==connFile2:
                    continue

                # finding the minMatch of closest points between two polygons lists
                # will try multi-core version and default to for loop in case of failure
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
                    temp = Parallel(n_jobs=-1,max_nbytes=1e6)(delayed(minDistancePoints)(ii) for ii in pairings)
                    # unpackign the tuple back to a numpy array of two variables
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
                            pointDistance_temp, points_temp = minDistancePoints(pairing)
                            points.append(points_temp)
                            pointDistance.append(pointDistance_temp)
                    stop = time.time()
                log.info('DONE min distance points')
                # now find the point pair with the shortest distance and keep the correpsonding points
                min_ix = np.where(pointDistance ==np.min(pointDistance))[0]
                dist = pointDistance[min_ix[0]]
                points = points[min_ix[0]]

                # track the coordinates of the closest points in geocoordinates
                point1 = np.array(points[0])
                point2 = np.array(points[1])

                # The list of points are defined from one polygon point to another.
                # Therefore each polygon mapping will consist out of two Source points.
                # The columns for each point are populated as:
                    # Source Unique CC ID, Source CC ID, Source CC FILE, Source UNW FILE, Source X-coordinate, Source Y-coordinate, Source geoTrans, Partner Unique CC, Dist to partner
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
                try:
                    tablePoints = np.concatenate((tablePoints,point1_array, point2_array), axis=0)
                except:
                    tablePoints = np.concatenate((point1_array, point2_array), axis=0)


        ## making sure the points are unique based on coordinate and unique connected compoenent id
        tablePoints_unique, unique_indices = np.unique(tablePoints[:, (0,4,5) ], axis=0,return_index=True)
        self.tablePoints =tablePoints[unique_indices[:],:]

        # compute the phase values
        log.info('Calculating the phase values for %s points', self.tablePoints.shape[0])
        # will try multi-core version and default to for loop in case of failure
        try:
            start = time.time()
            # need to combine all inputs together as single argument tuple
            all_inputs = ()
            for counter in range(self.tablePoints.shape[0]):
                connCompID1_temp = self.tablePoints[counter,1]
                connFile1_temp = self.tablePoints[counter,2]
                unwFile1_temp = self.tablePoints[counter,3]
                point1_temp = np.array([self.tablePoints[counter,4], self.tablePoints[counter,5]])
                geoTrans1_temp =self.tablePoints[counter,6]
                # parse inputs as a tuple
                inputs = (connCompID1_temp,connFile1_temp,unwFile1_temp,point1_temp,geoTrans1_temp,self.region)
                # append all tuples in a single tuple
                all_inputs = all_inputs + (inputs,)
            # compute the phase value using multi-thread functionality
            unwPhase_value = Parallel(n_jobs=-1,max_nbytes=1e6)(delayed(point2unwPhase)(ii) for ii in all_inputs)
            # end timer
            end = time.time()
        except:
            log.info('Multi-core version failed, will try single for loop')
            # will track processing time
            start = time.time()
            unwPhase_value = []
            for counter in range(self.tablePoints.shape[0]):
                connCompID1_temp = self.tablePoints[counter,1]
                connFile1_temp = self.tablePoints[counter,2]
                unwFile1_temp = self.tablePoints[counter,3]
                point1_temp = np.array([self.tablePoints[counter,4], self.tablePoints[counter,5]])
                geoTrans1_temp =self.tablePoints[counter,6]
                # parse inputs as a tuple
                inputs = (connCompID1_temp,connFile1_temp,unwFile1_temp,point1_temp,geoTrans1_temp,self.region)
                # compute the phase value
                unwPhase_temp = point2unwPhase(inputs)
                unwPhase_value.append(unwPhase_temp)
            # end timer
            end = time.time()
        log.info("runtime phase value calculation: %s sec", end - start)

        # Appending the phase to the table stored in self with the points
        unwPhase_value = np.array(unwPhase_value)
        unwPhase_value = unwPhase_value.reshape((unwPhase_value.shape[0],1))
        self.tablePoints = np.concatenate((self.tablePoints, unwPhase_value), axis=1)

        ## The two-stager only works with integer coordinate system in its triangulation, and needs unique points for each node. Will apply two steps to ensure this
        # STEP 1) will convert the geo-coordiantes with respect to the most southwest point and normalize with respect to original grid sampling
        # STEP 2) for non-unique points will shift by a pixel to get uniquesness back

        # step 1: generating raster coordinates
        geoTrans_global = list(geoTrans1)
        geoTrans_global[0] = np.min(self.tablePoints[:,4].astype(np.float64))
        geoTrans_global[3] = np.max(self.tablePoints[:,5].astype(np.float64))
        X_coordinate = ((self.tablePoints[:,4].astype(np.float64) - geoTrans_global[0])/geoTrans_global[1]).astype(np.int)
        Y_coordinate = ((self.tablePoints[:,5].astype(np.float64) - geoTrans_global[3])/geoTrans_global[5]).astype(np.int)
        X_coordinate = X_coordinate.reshape((X_coordinate.shape[0],1))
        Y_coordinate = Y_coordinate.reshape((Y_coordinate.shape[0],1))

        # step 2: making the nodes unqiue
        test_unique =  np.array([])
        test = np.array([1,1])
        while_count = 1
        while test_unique.shape!=test.shape:
            for count in range(X_coordinate.shape[0]):
                # find all the points that have the same coordiante as this one
                point = np.array([X_coordinate[count][0], Y_coordinate[count][0]])
                points = np.concatenate((X_coordinate, Y_coordinate), axis=1)
                matches = np.where(np.all(point==points,axis=1))[0]
                if matches.shape[0]>1:
                    #print('Found ' + str(matches.shape[0]) + ' non-unique points')
                    X_coordinate[matches[1:],0] = X_coordinate[matches[1:],0] + random.sample(range(-3,3),matches.shape[0]-1)
                    Y_coordinate[matches[1:],0] = Y_coordinate[matches[1:],0] + random.sample(range(-3,3),matches.shape[0]-1)

            test = np.concatenate((X_coordinate, Y_coordinate), axis=1)
            test_unique, test_unique_indices = np.unique(test, axis=0,return_index=True)
            log.info('Shifting coordinates by a pixel to make all points unique: iteration %s', while_count)
            while_count = while_count+1

        if test_unique.shape!=test.shape:
            log.error('Coordinates are not-unique, need to implement a while loop, run again for another random sample...')
            import pdb
            pdb.set_trace()
        else:
            log.info('Coordinates are unique')

        # appending the X and Y-coordinates to the points table.
        self.tablePoints = np.concatenate((self.tablePoints, X_coordinate), axis=1)
        self.tablePoints = np.concatenate((self.tablePoints, Y_coordinate), axis=1)


    def __compToCycle__(self, cycle, compnum):
        '''
            Generate a mapping from unique connected component to the number of phase cycles that needs to be applied
        '''

        compMap = {}
        for n, c in zip(cycle, compnum):
            try:
                compN = compMap[c]
                # Test if same cycle
                if (compN == n):
                    continue
                else:
                    raise ValueError("Incorrect phaseunwrap output: Different cycles in same components")
            except:
                # Component cycle doesn't exist in the dictionary
                compMap[c] = n

        return compMap


    def __calculateCycles__(self):
        '''
            Function that will calculate the number of cycles each component needs to be shifted in order to minimize the two-pi modulu residual between a neighboring component. Outputs a fileMappingDict with as key a file number. Within fileMappingDict with a integer phase shift value for each unique connected component.
        '''

        # import dependecies, looks liek the general import does not percolate down to this level.
        from ARIAtools.phaseMinimization import PhaseUnwrap

        # generating a list to be parsed in the two-stager, not that table originally was a mixture of string and floats, so all were temporaly mapped to a string for easy parsing. Will need to undo that now when calling two stager
        x = list(self.tablePoints[:,10].astype(np.int)+1)
        y = list(self.tablePoints[:,11].astype(np.int)+1)
        phase = (self.tablePoints[:,9].astype(np.float32))
        compNum = list(self.tablePoints[:,7].astype(np.int))
        phaseunwrap = PhaseUnwrap(x=x, y=y, phase=phase, compNum=compNum, redArcs=redarcsTypes[self.redArcs])
        phaseunwrap.solve(self.solver)
        cycles = phaseunwrap.unwrapLP()

        # Map unique component to integer number of cycles
        compMap = self.__compToCycle__(cycles, compNum)

        # Map component number for each file to integer number of cycles
        for countFile in self.fileMappingDict:
            fileMappingDict = self.fileMappingDict[countFile]

            # generate an array mapping connected component of each file to the integer phase cycles it needs to be shifted
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
            connCompMapping = np.concatenate((fileMappingDict['connCompMapping'], cycleMapping), axis=1)

            # update the original dictionary and store it back in self
            fileMappingDict['connCompMapping'] = connCompMapping
            self.fileMappingDict[countFile]=fileMappingDict
        log.info('MAPPING complete.')



def minDistancePoints(pairing):
    '''
        Finding closest two points between a list of two polygons
    '''

    poly1= pairing[0]
    poly2 = pairing[1]
    poly1_source = pairing[2]
    poly2_source = pairing[3]

    # Shapely nearest_point behaves well when two polygons do not overlap.
    # in case one polygon is contained within another, it will return the point on the inner polygon
    # the latter is not good behaviour in case of nested components within components
    # whill handle this manully by checking if the polygons have intersection.
    # However, note that in case of stiching multiple products, there will be a case that polygons intersect in along the edges in the overlap region.
    # We will track this by also including a check if the polygons originate form the same connected component file.

    ## check of the polygon has an intersection, two scenario's exist
    # 1) if from the same connected component file, then one component is within anoher, nearest points will return incorrect points
    # 2) if from differnt connected component files, then its form the overlap region, nearest poitns will not give an incorrect estimate as there is overlap

    #pdb.set_trace()
    Poly2 = Polygon(poly2)
    Poly1 = Polygon(poly1)

    if Poly1.intersection(Poly2).area > 0 and  poly1_source == poly2_source:
    #if poly1.intersects(poly2) and poly1_source == poly2_source:
        # will manully find the nearest points by looping over each vertex of polygon 2 and finding the minimum distance to polygon 1.

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
        Takes a GDAL compatible file as input, and returns shaply polygons and a list of corresponding connected components
    """

    # will save the geojson under a temp local filename
    tmfile = tempfile.NamedTemporaryFile(mode='w+b',suffix='.json', prefix='ConnPolygonize_', dir='.')
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
        # Connected component 0 is an actual data value, but not used for 2-stager and is removed here
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
        # find all the polygons for this connected component and add them as a list
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
        out_data=None
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
    '''
        Function to generate intermediate connected component files and unwrapped VRT files that have with interger 2pi pixel shift applied.
        Will parse inputs in a single argument as it allows for parallel processing.
        Return a list of files in a unqiue temp folder
    '''

    # parsing the inputs to variables
    saveDir = inputs['saveDir']
    saveNameID = inputs['saveNameID']
    connFile = inputs['connFile']
    unwFile = inputs['unwFile']
    connCompMapping = inputs['connCompMapping']

    ## Generating the intermediate files
    ## STEP 1: set-up the mapping functions
    # loading the connected component
    connData,connNoData,connGeoTrans,connProj = GDALread(connFile)

    ## Defining the mapping tables
    # setting up the connected component unique ID mapping
    connIDMapping = connCompMapping[:,1]
    # Will add the no-data to the mapping as well (such we can handle no-data region)
    # i.e. num comp  + 1 = Nodata connected component which gets mapped to no-data value again
    NoDataMapping = len(connIDMapping)
    connIDMapping = np.append(connIDMapping, [connNoData])

    # setting up the connected component integer 2PI shift mapping
    intMapping = connCompMapping[:,2]
    # Will add the no-data to the mapping as well (such we can handle no-data region)
    # i.e. max comp ID + 1 = Nodata connected component which gets mapped to 0 integer shift such no-data region remains unaffected
    intMapping = np.append(intMapping, [0])

    # update the connected component with the new no-data value used in the mapping
    connData[connData==connNoData]=NoDataMapping

    ## STEP 2: apply the mapping functions
    # interger 2PI scaling mapping for unw phase

    intShift = intMapping[connData.astype('int')]
    # connected component mapping to unique ID
    connData = connIDMapping[connData.astype('int')]

    ## STEP 3: writing out the datasets
    # writing out the unqiue ID connected component file
    connDataName = os.path.abspath(os.path.join(saveDir,'connComp', saveNameID + '_connComp.tif'))
    write_ambiguity(connData,connDataName,connProj,connGeoTrans,connNoData)

    # writing out the integer map as tiff file
    intShiftName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '_intShift.tif'))
    write_ambiguity(intShift,intShiftName,connProj,connGeoTrans)


    # writing out the scalled vrt => 2PI * integer map
    length = intShift.shape[0]
    width = intShift.shape[1]
    scaleVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '_scale.vrt'))
    build2PiScaleVRT(scaleVRTName,intShiftName,length=length,width=width)

    # Offseting the vrt for the range offset correctiom
    unwRangeOffsetVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '_rangeOffset.vrt'))
    buildScaleOffsetVRT(unwRangeOffsetVRTName,unwFile,connProj,connGeoTrans,File1_offset=connCompMapping[1,3],length=length,width=width)

    # writing out the corrected unw phase vrt => phase + 2PI * integer map
    unwVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '.vrt'))
    buildSumVRT(unwVRTName,unwRangeOffsetVRTName,scaleVRTName,connProj,connGeoTrans,length,width,description= inputs['description'])

    return [connDataName, unwVRTName]

def write_ambiguity(data, outName,proj, geoTrans,noData=False):
    '''
       Write out an integer mapping in the Int16/Byte data range of values
    '''

    # GDAL precision support in tif
    Byte = gdal.GDT_Byte
    Int16 = gdal.GDT_Int16

    # check if the path to the file needs to be created
    dirname = os.path.dirname(outName)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # Getting the GEOTIFF driver
    driver = gdal.GetDriverByName('GTIFF')
    # leverage the compression option to ensure small file size
    dst_options = ['COMPRESS=LZW']
    # create the dataset
    ds = driver.Create(outName , data.shape[1], data.shape[0], 1, Int16, dst_options)
    # setting the proj and transformation
    ds.SetGeoTransform(geoTrans)
    ds.SetProjection(proj)
    # populate the first band with data
    bnd = ds.GetRasterBand(1)
    bnd.WriteArray(data)
    # setting the no-data value
    if noData is not None:
        bnd.SetNoDataValue(noData)
    bnd.FlushCache()
    # close the file
    ds = None

def write_data(data, outName,proj, geoTrans,noData=False):
    '''
        Write out an phase dataset in the Float32 data range of values
    '''
    
    # GDAL precision support in envi
    Float32 = gdal.GDT_Float32
    
    # check if the path to the file needs to be created
    dirname = os.path.dirname(outName)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    
    # Getting the GEOTIFF driver
    driver = gdal.GetDriverByName('ENVI')
    # leverage the compression option to ensure small file size
    dst_options = ['COMPRESS=LZW']
    # create the dataset
    ds = driver.Create(outName , data.shape[1], data.shape[0], 1, Float32, dst_options)
    # setting the proj and transformation
    ds.SetGeoTransform(geoTrans)
    ds.SetProjection(proj)
    # populate the first band with data
    bnd = ds.GetRasterBand(1)
    bnd.WriteArray(data)
    # setting the no-data value
    if noData is not None:
        bnd.SetNoDataValue(noData)
    bnd.FlushCache()
    # close the file
    ds = None

def build2PiScaleVRT(output,File,width=False,length=False):
    '''
        Building a VRT file which scales a GDAL byte file with 2PI
    '''

    # DBTODO: The datatype should be loaded by default from the source raster to be applied.
    # should be ok for now, but could be an issue for large connected com

    # the vrt template with 2-pi scaling functionality
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

    # the inputs needed to build the vrt
    # load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # write out the VRT file
    with open(output, 'w') as fid:
        fid.write(vrttmpl.format(width = width, length = length, File = File))


def buildScaleOffsetVRT(output,File1,proj,geoTrans,File1_offset=0, File1_scale = 1, width=False,length=False,description='Scalled and offsetted VRT'):
    '''
        Building a VRT file which sums two files together using pixel functionality
        '''

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


    # the inputs needed to build the vrt
    # load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File1, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # write out the VRT file
    with open('{0}'.format(output) , 'w') as fid:
        fid.write( vrttmpl.format(width=width,length=length,File1=File1,File1_offset=File1_offset, File1_scale = File1_scale, proj=proj,geoTrans=str(geoTrans)[1:-1],description=description))

def buildSumVRT(output,File1,File2,proj,geoTrans,length=False, width=False,description='Unwrapped Phase'):
    '''
        Building a VRT file which sums two files together using pixel functionality
    '''

    # the vrt template with sum pixel functionality
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

    # the inputs needed to build the vrt
    # load the width and length from the GDAL file in case not specified
    if not width or not length:
        ds = gdal.Open(File1, gdal.GA_ReadOnly)
        width = ds.RasterXSize
        ysize = ds.RasterYSize
        ds = None

    # check if the path to the file needs to be created
    dirname = os.path.dirname(output)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    # write out the VRT file
    with open('{0}'.format(output) , 'w') as fid:
        fid.write( vrttmpl.format(width=width,length=length,File1=File1,File2=File2, proj=proj,geoTrans=str(geoTrans)[1:-1],description=description))

def point2unwPhase(inputs):
    '''
        Function to calculate the unwrapped phase value of a point and a region around it from only a single connected component
        Will parse inputs in a single argument as it allows for parallel processing.
    '''

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
    connData=ds_connFile.ReadAsArray(xoff=xoff,yoff=yoff,xsize=region, ysize=region)
    ds_connFile = None

    # loading a chunk of the unwrapped phase data
    ds_unwFile = gdal.Open(unwFile, gdal.GA_ReadOnly)
    unwData=ds_unwFile.ReadAsArray(xoff=xoff,yoff=yoff,xsize=region, ysize=region)
    ds_unwFile = None

    # taking the average phase for all the pixels of the connected component
    try:
        unwPhase = np.mean(unwData[connData==connCompID])
        np.where(connData==connCompID)
    except:
        log.error("Looks like the region is not large enough, cannot find your connected component: %s. Entering debug mode.", connCompID)
        print(connData)
        pdb.set_trace()

    # return back the phase value
    return unwPhase


def gdalTest(file, verbose=False):
    '''
        Checks if the file is GDAL compatible and return compatible VRT in case it exists
    '''
    # GDAL by default does not error out, so catch exception as error.
    # note warnings are not captured, e.g. warning in case netcdf file exist, but variable is not contained.
    gdal.UseExceptions()

    # DBTODO: Beyond GDAL compatible we should check if the file is GEO-CODED!
    #         e.g. use of proj, if not remove the file from list. Do not rely on .geo string in name

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
        # will first check if file exist, as for a netcdf file one will modify the variable and gdal will give a warning not and error
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



def product_stitch_overlap(unw_files, conn_files, prod_bbox_files, bbox_file, prods_TOTbbox, outFileUnw = './unwMerged', outFileConnComp = './connCompMerged', outputFormat='ENVI', mask=None, verbose=False):
    '''
        Stitching of products minimizing overlap betnween products
    '''

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

def product_stitch_morph(unw_files, conn_files, prod_bbox_files, bbox_file, prods_TOTbbox, outFileUnw = './unwMerged', outFileConnComp = './connCompMerged',outputFormat='ENVI',mask=None, verbose=False):
    '''
        Stitching of products using the two-stage unwrapper approach
        i.e. minimize the discontinuities between connected components
        '''
    
    import pdb
    
    # report method to user
    print('STITCH Settings: Connected component approach')
    print('Solver: MORPHOLOGY')

    
    # First, run the regular sticher, this will
    # (1) correct for range offset
    # (2) minimize the phase jump between adjacent product
    unw = UnwrapOverlap()
    unw.setInpFile(unw_files)
    unw.setConnCompFile(conn_files)
    unw.setOutFileConnComp(outFileConnComp + "_intermediate")
    unw.setOutFileUnw(outFileUnw + "_intermediate")
    unw.setProdBBoxFile(prod_bbox_files)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setVerboseMode(verbose)
    unw.UnwrapOverlap()


    # make some changes here
    # Second, run the Morph unwrap to minimizes phase jumps across connected component boundaries
    unw = UnwrapMorph()
    unw.setInpFile(os.path.abspath(outFileUnw + "_intermediate"))
    unw.setConnCompFile(os.path.abspath(outFileConnComp + "_intermediate"))
    unw.setOutFileConnComp(outFileConnComp)
    unw.setOutFileUnw(outFileUnw)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setVerboseMode(verbose)
    unw.UnwrapMorph()


def product_stitch_2stage(unw_files, conn_files, bbox_file, prods_TOTbbox, unwrapper_2stage_name = None, solver_2stage = None, outFileUnw = './unwMerged', outFileConnComp = './connCompMerged',outputFormat='ENVI',mask=None, verbose=False):
    '''
        Stitching of products using the two-stage unwrapper approach
        i.e. minimize the discontinuities between connected components
    '''

    # The solver used in minimizing the stiching of products
    if unwrapper_2stage_name is None:
        unwrapper_2stage_name = 'REDARC0'
    
    if solver_2stage is None:
        # If unwrapper_2state_name is MCF then solver is ignored
        # and relaxIV MCF solver is used by default
        solver_2stage = 'pulp'
    
    # report method to user
    print('STITCH Settings: Connected component approach')
    print('Name: %s'%unwrapper_2stage_name)
    print('Solver: %s'%solver_2stage)


    # First, run the regular sticher, this will
    # (1) correct for range offset
    # (2) minimize the phase jump between adjacent product
    unw = UnwrapOverlap()
    unw.setInpFile(unw_files)
    unw.setConnCompFile(conn_files)
    unw.setOutFileConnComp(outFileConnComp + "_intermediate")
    unw.setOutFileUnw(outFileUnw + "_intermediate")
    unw.setProdBBoxFile(prod_bbox_files)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setVerboseMode(verbose)
    unw.UnwrapOverlap()
    
    #pdb.set_trace()
    # Second, run the 2stager to minimizes phase jumps across connected component boundaries
    unw = UnwrapComponents()
    unw.setInpFile(os.path.abspath(outFileUnw + "_intermediate"))
    unw.setConnCompFile(os.path.abspath(outFileConnComp + "_intermediate"))
    unw.setOutFileConnComp(outFileConnComp)
    unw.setOutFileUnw(outFileUnw)
    unw.setSolver(solver_2stage)
    unw.setRedArcs(unwrapper_2stage_name)
    unw.setMask(mask)
    unw.setOutputFormat(outputFormat)
    unw.setBBoxFile(bbox_file)
    unw.setTotProdBBoxFile(prods_TOTbbox)
    unw.setVerboseMode(verbose)
    unw.UnwrapComponents()


def get_all_overlaps(cc,scc,nbest=3,struct=None,niter=2,minsize=20):
    """
        Given the array of ccomp sorted by size get all possible overlaps using the first "nbest"
        inputs:
        cc = ccomp file.
        scc = array with all of sorted ccomp according to size.
        nbests = use the first largest nbest cc in scc to get the overlap with others cc.
        struct = morphological structure.
        niter = number iterations of morphological operations.
        minsize = minimum number of points in the overlaps. If less than minsize then no overlap.
        outputs:
        overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        complete = bool True if all the cc have overlap
        dones: the list of cc that have overlap
        """
    
    #keep trak of the cc that have not been adjusted. initially are all except for -1 and 0
    check_complete = np.unique(cc)
    check_complete = np.setdiff1d(check_complete,np.array([-1,0]),True)
    overlaps = {}
    #get all the pairs that overlap
    overlaps0,complete0,dones0 = get_overlaps(cc,scc[:nbest],overlaps,check_complete,struct,niter,minsize)
    #get the new list of cc minus -1,0
    check_complete = np.unique(cc)
    check_complete = np.setdiff1d(check_complete,np.array([-1,0]),True)
    #get the ones that still need to be done
    not_dones = np.setdiff1d(check_complete,dones0)
    #redo the overlap searching for the not_dones instead of nbest
    overlaps,complete,dones = get_overlaps(cc,not_dones,overlaps,not_dones,struct,niter,minsize)
    dones = np.unique(np.array(np.append(dones0,dones)))
    #check if all have been fixed. sometime is not possible if there the cc is isolated
    if len(dones) == len(check_complete):
        complete = True
    else:
        complete = False
    return overlaps,complete,dones

def order_overlaps(pairs):
    """
        Create an OrderedDict with key the cc and value a list of cc that have overlap with the key.
        The order is based on the size of the values i.e. how many cc overlap with the key.
        inputs:
        pairs = list of cc pairs that overlap. The form is i_j where i,j are the cc.
        output:
        ddict = OrderedDict with key a cc and value an array of the cc overlapping with key.
        """
    
    ddict = collections.defaultdict(list)
    for p in pairs:
        keys = p.split('_')
        ddict[int(keys[0])].append(int(keys[1]))
    sizes = collections.defaultdict(int)
    for k,v in ddict.items():
        sizes[k] = len(v)
    ssizes = collections.OrderedDict(sorted(sizes.items(),key = lambda t:t[1],reverse=True))
    ddict = collections.OrderedDict([(k,ddict[k]) for k in ssizes.keys()])
    return ddict

def get_resized_imgs(imgs,mask,buf=50):
    """
        Resize each image in images only where mask is non zero. Add a +- buf
        inputs:
        imgs = images to resize.
        mask = the mask to find the extremes.
        buf = add buf to the extremes.
        outputs:
        ret_imgs = the resized images
        min_x = the minimun x location in the original image
        min_y = the minimum y location in the original image
        """
    
    #find the extremes where mask is to reduce the computation
    locs = np.where(mask)
    min_y = max(locs[0].min() - buf,0)
    max_y = max(locs[0].max() + buf,mask.shape[0])
    min_x = max(locs[1].min() - buf,0)
    max_x = max(locs[1].max() + buf,mask.shape[1])
    ret_imgs = []
    for im in imgs:
        ret_imgs.append(im[min_y:max_y,min_x:max_x])
    return ret_imgs,min_x,min_y

def get_overlaps(cc,scc,overlaps,check_complete,struct=None,niter=2,minsize=20):
    """
        Given the array of ccomp sorted by size get all possible overlaps using the first "nbests"
        inputs:
        cc = ccomp file.
        scc = array with subset of sorted ccomp according to size.
        overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        check_complete = list of the cc that have no overlap.
        struct = morphological structure.
        niter = number iterations of morphological operations.
        minsize = minimum number of points in the overlaps
        outputs:
        overlaps = the updated input 'overlaps'.
        complete = bool True if all the cc have overlap.
        dones = the list of cc that have overlap.
        """
    
    if struct is None:
        struct = generate_binary_structure(2, 1)
    
    dones = []
    complete = False
    for cc_now in scc:
        mask = cc == cc_now
        #resize cc and mask only around the extremes of mask
        [ncc,mask],min_x,min_y = get_resized_imgs([cc,mask],mask,buf=50)
        #get what's not cc_now
        nmask = np.logical_not(mask)
        #dilate mask
        dmask = binary_dilation(mask,struct,niter)
        #get intersection between dilate cc and not cc
        inters = np.logical_and(dmask,nmask)
        #select all the cc that are in the dilated part only
        sel_now = np.unique(ncc[inters])
        for i in sel_now:
            if i == -1:
                continue
            if i == 0:
                continue
            okey = str(i) + '_' + str(cc_now)
            #if already present skip
            if okey in overlaps:
                continue
            mask1 = ncc == i
            [mask1,mmask,ddmask],min_xx,min_yy = get_resized_imgs([mask1,mask,dmask],mask1)
            dmask1 = binary_dilation(mask1,struct,niter)
            inters10 = np.logical_and(ddmask,mask1)
            inters01 = np.logical_and(dmask1,mmask)
            if np.count_nonzero(inters01) >= minsize and np.count_nonzero(inters10) >= minsize:
                #if enough points overlap that save the ovelapping points.
                #since we cropped the images need to put back the correct indices
                ninters01 = np.zeros_like(cc,np.bool)
                ninters01[min_y + min_yy:min_y + min_yy + inters01.shape[0],min_x + min_xx:min_x + min_xx + inters01.shape[1]] = inters01
                ninters10 = np.zeros_like(cc,np.bool)
                ninters10[min_y + min_yy:min_y + min_yy + inters10.shape[0],min_x + min_xx:min_x + min_xx + inters10.shape[1]] = inters10
                overlaps[str(cc_now) + '_' + str(i)] = [ninters01,ninters10]
                #mark the 2 cc as done
                dones.append(cc_now)
                dones.append(i)
                if len(check_complete) == len(np.intersect1d(check_complete,np.unique(np.array(dones)))):
                    complete = True
            if complete == True:
                break
        if complete == True:
            break
    return overlaps,complete,np.unique(np.array(dones))

def get_overlaps_obs(cc,scc,overlaps,check_complete,struct=None,niter=2,minsize=20):
    """
        This method is deprecated (obsolete)
        Given the array of ccomp sorted by size get all possible overlaps using the first "nbests"
        inputs:
        cc = ccomp file.
        scc = array with subset of sorted ccomp according to size.
        overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        check_complete = list of the cc that have no overlap.
        struct = morphological structure.
        niter = number iterations of morphological operations.
        minsize = minimum number of points in the overlaps
        outputs:
        overlaps = the updated input 'overlaps'.
        complete = bool True if all the cc have overlap.
        dones = the list of cc that have overlap.
        """
    
    if struct is None:
        struct = generate_binary_structure(2, 1)
    
    dones = []
    complete = False
    for cc_now in scc:
        mask = cc == cc_now
        nmask = np.logical_not(mask)
        dmask = binary_dilation(mask,struct,niter)
        inters = np.logical_and(dmask,nmask)
        sel_now = np.unique(cc[inters])
        for i in sel_now:
            if i == -1:
                continue
            if i == 0:
                continue
            okey = str(i) + '_' + str(cc_now)
            #if already present skip
            if okey in overlaps:
                continue
            mask1 = cc == i
            dmask1 = binary_dilation(mask1,struct,niter)
            inters10 = np.logical_and(dmask,mask1)
            inters01 = np.logical_and(dmask1,mask)
            if np.nonzero(inters01)[0].size >= minsize and np.nonzero(inters10)[0].size >= minsize:
                overlaps[str(cc_now) + '_' + str(i)] = [inters01,inters10]
                dones.append(cc_now)
                dones.append(i)
                if len(check_complete) == len(np.intersect1d(check_complete,np.unique(np.array(dones)))):
                    complete = True
            if complete == True:
                break
        if complete == True:
            break
    return overlaps,complete,np.unique(np.array(dones))

def swap_single(ord_overlaps,overlaps):
    """
        Swap the key and value (when there is only one value) of the ordered overlaps when the key is not part
        of any other pair of cc overlapping but the value is so that the pair is connected to others.
        inputs:
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        overlaps = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        outputs:
        nord_overlaps = updated ord_overlaps.
        noverlaps = updates overlaps.
        """
    
    #when one cc is inside the other the one inside has more overlap because of the dilation, but it's not
    #connected to anything. So check if for single element values either the key or the value if part
    #of some other cc set. If the value is and not the key then swap since most likely the key is inside the value.
    #get_ipython().run_line_magic('debug', '')
    n_ord_over = {}
    nover = {}
    to_delete = []
    if len(ord_overlaps) > 1:
        for k,v in ord_overlaps.items():
            v = copy.deepcopy(v)
            if len(v) == 1:
                for k1,v1 in ord_overlaps.items():
                    if k1 == k:
                        continue
                    if v[0] in v1 and k not in v1:#the second check is to see if both v and k are in v1
                        #if this is the case leave it alone otherwise remove it
                        n_ord_over[v[0]] = [k]
                        nover[str(v[0]) + '_' + str(k)] =  overlaps[str(k) + '_' + str(v[0])][::-1]#swap order
                        to_delete.append(str(k) + '_' + str(v[0]))
                    elif v[0] in v1 and k in v1:
                        continue
                    else:
                        n_ord_over[k] = v
                        nover[str(k) + '_' + str(v[0])] = overlaps[str(k) + '_' + str(v[0])]
            else:
                n_ord_over[k] = v
                for vi in v:
                    nover[str(k) + '_' + str(vi)] = overlaps[str(k) + '_' + str(vi)]
    
        #remove the one swapped
        to_delete = np.unique(np.array(to_delete))
        for k in to_delete:
            if k in nover:
                del nover[k]
                k0 = int(k.split('_')[0])
                k1 = int(k.split('_')[1])
                n_ord_over[k0].remove(k1)
                if len(n_ord_over[k0]) == 0:
                    n_ord_over.pop(k0, None)
    else:#fringe case where there is only one left and the above logig doesn't work
        n_ord_over = ord_overlaps
        nover = overlaps
    
    return collections.OrderedDict(sorted(n_ord_over.items(),key = lambda t:t[0])),nover

def remove_dups(ord_overlaps,nover):
    """
        Remove cc that are in multiple places. Keep the instance that has larger overlap.
        inputs:
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        nover = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        outputs:
        nord_overlaps = updated ord_overlaps.
        nover = updates nover.
        """
    
    keys = list(ord_overlaps.keys())
    for i in range(len(keys) - 1):
        k = keys[i]
        for j in range(i + 1,len(keys)):
            k1 = keys[j]
            common = np.intersect1d(ord_overlaps[k],ord_overlaps[k1])
            for c in common:
                key = str(k) + '_' + str(c)
                key1 = str(k1) + '_' + str(c)
                if len(nover[key][0]) > len(nover[key1][0]):
                    ord_overlaps[k1].remove(c)
                    del nover[key1]
                else:
                    ord_overlaps[k].remove(c)
                    del nover[key]
    keys = list(ord_overlaps.keys())
    for k in keys:
        v = ord_overlaps[k]
        if len(v) == 0:
            del ord_overlaps[k]
    return ord_overlaps,nover

def reorder_overlaps(ord_overlaps):
    """
        Order the inputs so that when cycles are added to a cc no prior modified cc is affected.
        inputs:
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        outputs:
        nord_overlaps = the updated ord_overlaps.
        """
    
    #make sure that the ovelarps are sorted so that when one is modified no other connected to it
    #get modified later one. For instance 10:[11] and 13:[10,20]. In this order 11 is modified w.r.t. 10 but then
    #10 is modified w.r.t 13 and 11 now is off
    keys = list(ord_overlaps.keys())
    ret = []
    for i,k in enumerate(keys):
        swapped = False
        for k1 in keys[i+1:]:
            if k in ord_overlaps[k1]:
                #NOTE:k cannot be present in any other list because we remove dups
                #     so it's safe to swap and break
                if k1 in ret:#already present from previous swap
                    ret.append(k)
                else:
                    ret.extend([k1,k])
                swapped = True
                break
        if not swapped:#there was now connection with any other cc so just add
            ret.append(k)
    return collections.OrderedDict([(k,ord_overlaps[k]) for k in ret])

def get_cycles(unw_ph,ord_overlaps,nover):
    """
        Get the number of cycles that need to be added to a cc based on the difference of the phase in the
        overlapping region. Only the values in the ord_overlap dict are updated while the key has 0 cycles (unless
        a key is a value for another key in the dict)
        inputs:
        unw_ph: original unwrapped phase.
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        nover = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        outputs:
        cycles = dictionary with key a cc and value the number of cycles to be added.
        stats = a dictionary with key the overlapping pair 'i_j' and value an array with the cycles,
        stds of the phase in the two overlapping regions (normalized to 2pi), the number of elementes
        in the overlapping regions.
        """
    
    cycles = collections.defaultdict(int)
    stats = collections.defaultdict(list)
    first_time = True
    for k,v in ord_overlaps.items():
        if first_time:
            cycles[k] = 0
            first_time = False
        for k1 in v:
            key = str(k) + '_' + str(k1)
            inters,inters1 = nover[key]
            cyc = round((unw_ph[inters].mean() - unw_ph[inters1].mean())/(2*np.pi))
            cycles[k1] =  cyc
            stats[key] = [cyc,unw_ph[inters].std()/(2*np.pi),unw_ph[inters1].std()/(2*np.pi),
                          np.count_nonzero(inters),np.count_nonzero(inters1)]
    return cycles,stats

def alter_cycles(unw_ph,cc,selcc,ncyc):
    """
        Alterate the unwrapped phase by adding an interger times of 2pi to a given cc.
        inputs:
        unw_ph = original unwrapped phase.
        cc = the connected component image.
        selcc = array with the cc that need to be modified.
        ncyc = array with the cycles to add to the cc in selcc.
        outputs:
        nunw_ph = the modified unwrapped phase.
        """
    
    nunw_ph = unw_ph.copy()
    for i,j in zip(selcc,ncyc):
        seli = cc == i
        nunw_ph[seli] += 2*np.pi*j
    return nunw_ph

def get_images(dir_now,unw_file,cc_file):
    """
        Load the cc and unwrapped phase images into a numpy array.
        inputs:
        dir_now = directory containing the image files. Needs the xml metadata and the binary.
        outputs:
        cc: 2D int8 numpy array with the connected component image values.
        unw_ph: 2D float32 numpy array with the unwrapped phase image values.
        """
    
    data = gdal.Open(unw_file,gdal.GA_ReadOnly)
    data_band = data.GetRasterBand(1)
    unw_ph = data_band.ReadAsArray()
    data = None
    # doing the connected component
    data = gdal.Open(cc_file,gdal.GA_ReadOnly)
    data_band = data.GetRasterBand(1)
    cc = data_band.ReadAsArray()
    data = None
    img = None
    return cc,unw_ph,img

def compute_cycles(unw_ph,cc,nbest=3,niter=2,minsize=20,debug=True):
    """
        Compute the cycles to be added to the unpwrapped phase.
        inputs:
        unw_ph = original unwrapped phase.
        cc = the connected component image.
        nbests = use the first largest nbest cc in scc to get the overlap with others cc.
        niter = number iterations of morphological operations.
        minsize = minimum number of points in the overlaps.
        debug = bool if True prints intermediate steps results.
        output:
        cycles = dictionary with key a cc and value the number of cycles to be added.
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        nover = dictionary with key = i_j and values the two arrays which are the overlaps after
        dilating one cc (i or j) and not the other (j or i).
        """
    
    struct = generate_binary_structure(2, 1) # to identify neighbors
    cc_sizes = []
    
    
    # remove both the no-data and connected component 0 from the list under consideration
    cc_left = np.setdiff1d(np.unique(cc),[-1,0],True)
    for c in cc_left:
        cc_sizes.append(np.count_nonzero(cc == c))
    cc_sizes = np.array(cc_sizes,np.int32)
    scc_sizes = cc_left[np.argsort(cc_sizes)[::-1]]

    overlaps,complete,dones = get_all_overlaps(cc,scc_sizes,nbest=nbest,struct=struct,niter=niter,minsize=minsize)
    #get_all_overlaps(cc,scc_sizes,2)
    tord_overlaps = order_overlaps(overlaps.keys())
    if debug:
        print('order_overlaps',tord_overlaps)
    ord_overlaps,nover = swap_single(tord_overlaps,overlaps)
    if debug:
        print('swap_single',ord_overlaps,'\n',nover.keys())
    ord_overlaps,nover = remove_dups(ord_overlaps,nover)
    if debug:
        print('remove_dups',ord_overlaps,'\n',nover.keys())
    ord_overlaps = reorder_overlaps(ord_overlaps)
    if debug:
        print('reorder_overlaps',ord_overlaps)
    cycles,stats = get_cycles(unw_ph,ord_overlaps,nover)
    return cycles,ord_overlaps,nover

def add_cycles(unw_ph,cc,cycles,ord_overlaps):
    """
        Add cycles to the unwrapped image to remove discontinuity.
        inputs:
        unw_ph = original unwrapped phase.
        cc = the connected component image.
        cycles = dictionary with key a cc and value the number of cycles to be added.
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        outputs:
        nunw_ph = the fixed unwrapped phase.
        """
    
    cp_cycles = copy.deepcopy(cycles)
    keys = list(ord_overlaps.keys())
    #go one key at the time except for the first since it's the zero shift
    for i in range(1,len(keys)):
        key = keys[i]
        #go back to all the previous main cc and see if it belongs to one of the overlaps. if so need
        #to adjust all the cycles for this cc
        for j in range(i):
            key1 = keys[j]
            if key in ord_overlaps[key1]:#cc overlaps with previous one.
                for v in ord_overlaps[key]:
                    #print(key1,key,cp_cycles[key],v,cp_cycles[v])
                    cp_cycles[v] += cp_cycles[key]
    to_change = []
    max_cc    = np.ptp(cc)
    for i in range(max_cc):
        if i in cp_cycles:
            to_change.append(cp_cycles[i])
        else:
            to_change.append(0)
    to_change = np.array(to_change)
    return unw_ph + 2*np.pi*to_change[cc]

def get_linked_cc(ord_overlaps,max_cc):
    """
        Find all the cc that linked to each other.
        inputs:
        ord_overlaps = OrderedDict with key a cc and value an array of the cc overlapping with key.
        max_cc = maximum number of connected components in the cc image.
        outputs:
        linked = list of all the cc linked to each other.
        """
    
    linked = []
    keys = list(ord_overlaps.keys())
    done = []
    for i,(k,v) in enumerate(ord_overlaps.items()):
        #if k and it's v are already linked than skip them
        if k in done:
            continue
        done.append(k)
        link_now = [k] + v
        #go through the remaining keys and if k1 was part of the value of another key
        #it means that they are linked. ex. k,v = (1, [15, 16, 17, 18, 13, 19, 12]) and k1,v1 = (13, [10, 20])
        #at the end they should all belong to the same cc.
        for k1 in keys[i+1:]:
            if k1 in link_now:
                link_now += ord_overlaps[k1]
                done.append(k1)
        linked.append(link_now)
    not_singles = np.array([i for j in linked for i in j])
    singles = np.setdiff1d(np.arange(1,max_cc),not_singles)
    singles.sort()
    return linked,singles

def adjust_cc(cc,linked):
    """
        Adjust the connected component values of the linked cc.
        inputs:
        cc = the connected component image.
        linked = list of all the cc linked to each other.
        outputs:
        ncc = adjusted connected component image.
        """
    ncc = cc.copy()
    for l in linked:
        n_value = l[0]
        for l1 in l[1:]:
            sel = ncc == l1
            ncc[sel] = n_value
    return ncc

def adjust_to_zero_cc(unw_ph,cc,struc=None,niter=2,minsize=20):
    """
        Adjust non zero cc to the zero one. Used to make the phase image look smooth.
        inputs:
        unw_ph = original unwrapped phase.
        cc = the connected component image.
        struct = morphological structure.
        niter = number iterations of morphological operations.
        minsize = minimum number of points in the overlaps.
        outputs:
        runw_ph = modified pahse image.
        """
    
    check_complete = np.unique(cc)[1:]
    overlaps = {}
    overlaps,complete0,dones0 = get_overlaps(cc,[0],overlaps,check_complete,struct=struc,niter=niter,minsize=minsize)
    ord_overlaps = {0:[int(i.split('_')[1]) for i in overlaps.keys()]}
    cycles,stats = get_cycles(unw_ph,ord_overlaps,overlaps)
    runw_ph = add_cycles(unw_ph,cc,cycles,ord_overlaps)
    return runw_ph

