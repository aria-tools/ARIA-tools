#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os

import argparse
from phaseMinimization import Vertex, Edge, PhaseUnwrap
import time

import sys
import logging

import numpy as np
from osgeo import gdal
from gdalconst import *
import pdb
import tempfile
import shutil
import subprocess as sp
from shapely.geometry import Point,Polygon,shape, LinearRing
from shapely.ops import nearest_points
import copy
import json
from joblib import Parallel, delayed, dump, load
import random
import glob


solverTypes = ['pulp', 'glpk', 'gurobi']
redarcsTypes = {'MCF':-1, 'REDARC0':0, 'REDARC1':1, 'REDARC2':2}

class UnwrapComponents:
    ''' 
        2-Stage Phase Unwrapping
    '''
    
  
    def __init__(self):
        '''
            Setting the default arguments needed by the class.
            Fill parse the file-names as None as they need to be set by the user, which will be caught when running the class
        '''
        self.inpFile = None
        self.ccFile = None
        self.solver ='pulp'
        self.redArcs =-1
    
    def unwrapComponents(self):
        
        if self.inpFile is None:
            print("Error. Input unwrapped file(s) is (are) not set.")
            raise Exception

        if self.ccFile is None:
            print("Error. Input Connected Components file(s) is (are) not set.")
            raise Exception

        if self.solver not in solverTypes:
            raise ValueError(self.treeType + ' must be in ' + str(unwTreeTypes))

        if self.redArcs not in redarcsTypes.keys():
            raise ValueError(self.redArcs + ' must be in ' + str(redarcsTypes))

        '''
            if self.outFile is None:
            self.outFile = self.inpFile.split('.')[0] + '_2stage.xml' 
        '''
                
        self.region=5
        self.__run__()


        return

    def setRedArcs(self, redArcs):
        """ Set the Redundant Arcs to use for LP unwrapping """
        self.redArcs = redArcs

    def setSolver(self, solver):
        """ Set the solver to use for unwrapping """
        self.solver = solver

    def setInpFile(self, input):
        """ Set the input Filename for 2-stage stitching/unwrapping """
        # Convert a string (i.e. user gave single file) to a list
        if isinstance(input, np.str):
            input = [input]
        self.inpFile = input
        
        # the number of files that needs to be merged/two-staged unwrapped
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
            

    
    
    def __verifyInputs__(self):
        '''
            Verify if the unwrapped and connected component inputs are gdal compatible.
            If not remove them from the list to be stiched.
            If a vrt exist and gdalcompatible update the file to be a vrt
        '''
        # track a list of files to keep
        inpFile_keep = []
        ccFile_keep = []
        for k_file in range(self.nfiles):
            # unw and corresponding conncomponent file
            inFile = self.inpFile[k_file]
            ccFile = self.ccFile[k_file]
            
            # vrt files are prefered as they contain proj and transf information
            # Convert to inputs to vrt if exist, and check for gdal compatibility
            inFile_temp = gdalTest(inFile)
            ccFile_temp = gdalTest(ccFile)
            
            # if one of the two files fails then do not try and merge them later on as GDAL is leveraged
            if inFile_temp is None or ccFile_temp is None:
                print('Removing following pair combination (not GDAL compatible)')
                print('UNW: ' + inFile)
                print('CONN: ' + ccFile)
            else:
                inpFile_keep.append(inFile_temp)
                ccFile_keep.append(ccFile)
        
        # update the input files and only keep those being GDAL compatible
        self.inpFile=inpFile_keep
        self.ccFile=ccFile_keep
        
        # update the number of file in case some got removed
        self.nfiles = np.shape(self.inpFile)[0]
        
        if self.nfiles==0:
            print('No files left after GDAL compatibility check')
            sys.exit(0)


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
                    print('Multi-core')
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
                    print('For loop')
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
                print('DONE min distance points')
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
        

        pdb.set_trace()  


        ## making sure the points are unique based on coordinate and unique connected compoenent id
        tablePoints_unique, unique_indices = np.unique(tablePoints[:, (0,4,5) ], axis=0,return_index=True)
        self.tablePoints =tablePoints[unique_indices[:],:]
        
        # compute the phase values
        print('Calculating the phase values for ' + str(self.tablePoints.shape[0]) + ' points')
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
            print('Multi-core version failed, will try single for loop')
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
        print("runtime phase value calculation: " + str(end - start) + " sec")

        # Appending the phase to the table stored in self with the points
        unwPhase_value = np.array(unwPhase_value)
        unwPhase_value = unwPhase_value.reshape((unwPhase_value.shape[0],1))
        self.tablePoints = np.concatenate((self.tablePoints, unwPhase_value), axis=1)
            
        ## The two-stager only works with integer coordinate system in its triangulation, and needs unique points for each node. Will apply two steps to ensure this
        # STEP 1) will convert the geo-coordiantes with respect to the most southwest point and normalize with respect to original grid sampling
        # STEP 2) for non-unqiue points will shift by a pixel to get uniquesness back
        
        # step 1: generating raster coordiantes
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
            print('Shifting coordinates by a pixel to make all points unique: iteration ' + str(while_count))
            while_count = while_count+1
            
        if test_unique.shape!=test.shape:
            print('Coordinates are not-unique, need to implement a while loop, run again for another random sample...')
            pdb.set_trace()
        else:
            print('Coordinates are unique')
        
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
        pdb.set_trace()
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
                        print('Was not expecting this')
                        pdb.set_trace()
            cycleMapping = np.array(cycleMapping)
            cycleMapping = cycleMapping.reshape((cycleMapping.shape[0],1))
            connCompMapping = np.concatenate((fileMappingDict['connCompMapping'], cycleMapping), axis=1)

            # update the original dictionary and store it back in self
            fileMappingDict['connCompMapping'] = connCompMapping
            self.fileMappingDict[countFile]=fileMappingDict
        print('MAPPING complete:')
        pdb.set_trace()
    
    def __createImages__(self):
        '''
            This function will write the final merged unw and conencted component file. As intermediate step tiff files are generated with integer values which will represent the shift to be applied to connected componenet and the moduli shift to be applied to the unwrapped phase.
        '''
    
        ## Will first make intermediate files in a temp folder.
        # For each product there will be a connComp file and 3 files related unw files.
        # The connected component file will show the unique wrt to all merged files
        # For the unwrapped related files, there will be an integer offset tif file, a vrt file which scale this integer map by 2pi, and a vrt which combines the orginal unw phase file with the scaled map. The latter will be used for merging of the unwrapped phase.
        tempdir = tempfile.mkdtemp(prefix='IntermediateFiles_',dir='.')
        
        # will try multi-core version and default to for loop in case of failure
        try:
            fdsdfsfsFAILINGtheCODE
            start = time.time()
            # need to combine all inputs together as single argument tuple
            all_inputs = ()
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_" + str(counter)
                # parse inputs as a tuple
                inputs = (fileMappingDict)
                # append all tuples in a single tuple
                all_inputs = all_inputs + (inputs,)
            pdb.set_trace()
            # compute the phase value using multi-thread functionality
            intermediateFiles = Parallel(n_jobs=-1,max_nbytes=1e6)(delayed(createConnComp_Int)(ii) for ii in all_inputs)
            # end timer
            end = time.time()
    
        except:
            print('Multi-core version failed, will try single for loop')
            # will track processing time
            start = time.time()
            intermediateFiles = []
            for counter in  range(len(self.fileMappingDict)):
                fileMappingDict = self.fileMappingDict[counter]
                fileMappingDict['saveDir'] = tempdir
                fileMappingDict['saveNameID'] = "Product_n" + str(counter)
                # parse inputs as a tuple
                inputs = (fileMappingDict)
                # compute the phase value
                intermediateFiles_temp = createConnComp_Int(inputs)
                intermediateFiles.append(intermediateFiles_temp)
            # end timer
            end = time.time()
        print("runtime phase value calculation: " + str(end - start) + " sec")

        pdb.set_trace()

        # combining all conComp and unw files that need to be blended
        conCompFiles = []
        unwFiles = []
        for intermediateFile in intermediateFiles:
            conCompFiles.append(intermediateFile[0])
            unwFiles.append(intermediateFile[1])

        ## Will now merge the unwrapped and connected component files
        cmd = 'gdalwarp -co TILED=YES -srcnodata 0 -of envi ' + ' '.join(unwFiles) + ' unwPhase_merged.unwrapped'
        if os.path.isfile('unwPhase_merged.unwrapped'):
            for file in glob.glob("unwPhase_merged.unwrapped*"):
                os.remove(file)
        os.system(cmd)
        cmd = 'gdal2isce_xml.py -i unwPhase_merged.unwrapped'
        os.system(cmd)
            
        cmd = 'gdalwarp -co TILED=YES -srcnodata -1 -of envi ' + ' '.join(conCompFiles) + ' connComp_merged.conn'
        if os.path.isfile('connComp_merged.conn'):
            for file in glob.glob("connComp_merged.conn*"):
                os.remove(file)
        os.system(cmd)
        cmd = 'gdal2isce_xml.py -i connComp_merged.conn'
        os.system(cmd)
        
        # Remove the directory with intermediate files as they are no longer needed
        shutil.rmtree(tempdir)

        print('DONE')

    def __run__(self):
        '''
            The main driver.
        '''
        
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

        ## inform user of succesfull completion
        print('--DONE--')


        return


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
        print('intersect: min distance ' + str(dist) )

    else:  
        if Poly1.intersection(Poly2).area > 0:
            print("product overlap")
            pdb.set_trace()

        # leverage shapely distance function
        dist = poly1.distance(poly2)

        # track the coordinates of the closest points in geocoordinates
        points = nearest_points(poly1,poly2)
        point1 = np.array(points[0])
        point2 = np.array(points[1])

    #pdb.set_trace() 
   
    return dist, points

def product_stitch(unw_files, conn_files,unwrapped2StageFilename = None, unwrapper_2stage_name = None, solver_2stage = None):
      # The solver used in minimizing the stiching of products
      
      if unwrapper_2stage_name is None:
          unwrapper_2stage_name = 'REDARC0'
  
      if solver_2stage is None:
          # If unwrapper_2state_name is MCF then solver is ignored
          # and relaxIV MCF solver is used by default
          solver_2stage = 'pulp'
      '''
      if unwrapped2StageFilename is None:
          unwrapped2StageFilename = unwrappedIntFilename.replace('.unw','_2stage.unw')
      '''

      print('STICH/UNWRAP Settings:')
      print('Name: %s'%unwrapper_2stage_name)
      print('Solver: %s'%solver_2stage)
  
  
      out_file = 'test.unw'
  
      # Hand over to 2Stage unwrap
      unw = UnwrapComponents()
      unw._legacy_flag = True
      unw.setInpFile(unw_files)
      unw.setConnCompFile(conn_files)
      unw.setOutFile(out_file)
      unw.setSolver(solver_2stage)
      unw.setRedArcs(unwrapper_2stage_name)
      unw.unwrapComponents()

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
    sp.check_call(['gdal_polygonize.py',ccFile,'-f', 'GeoJSON', outname])
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

def GDALread(filename,data_band=None):
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
    if data_band is not None:
        raster = data.GetRasterBand(data_band)
        out_data = raster.ReadAsArray()
    else:
        # load the data
        out_data = data.ReadAsArray()

    # getting the gdal transform and projection
    geoTrans = data.GetGeoTransform()
    proj = data.GetProjection()

    # getting the no-data value
    try:
        NoData = data.GetNoDataValue()
        print(NoData)
    except:
        print('DBTODO Could not find a no-data value...')

        NoData = None
        NoData=-1

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
    width = intShift.shape[0]
    length = intShift.shape[1]
    scaleVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '_scale.vrt'))
    build2PiScaleVRT(scaleVRTName,intShiftName,length,width)

    # writing out the corrected unw phase vrt => phase + 2PI * integer map
    unwVRTName = os.path.abspath(os.path.join(saveDir,'unw',saveNameID + '.vrt'))
    buildSumVRT(unwVRTName,unwFile,intShiftName,connProj,connGeoTrans,length,width)
    
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

def buildSumVRT(output,File1,File2,proj,geoTrans,width=False,length=False):
    '''
        Building a VRT file which sums two files together using pixel functionality
    '''
    
    # the vrt template with sum pixel functionality
    vrttmpl = '''<VRTDataset rasterXSize="{width}" rasterYSize="{length}">
    <SRS>{proj}</SRS>
    <GeoTransform>{geoTrans}</GeoTransform>
    <VRTRasterBand dataType="Float32" band="1" subClass="VRTDerivedRasterBand">
        <Description>Two-stage corrected unwrapped phase</Description>
        <PixelFunctionType>sum</PixelFunctionType>
        <SimpleSource>
            <SourceFilename>{File1}</SourceFilename>
            <NoData>0</NoData>
        </SimpleSource>
        <SimpleSource>
            <SourceFilename >{File2}</SourceFilename>
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
        fid.write( vrttmpl.format(width=width,length=length,File1=File1,File2=File2, proj=proj,geoTrans=str(geoTrans)[1:-1]))


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

    #pdb.set_trace()  

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
        pdb.set_trace()
        np.where(connData==connCompID)
    except:
        print("Looks like the region is not large enough, cannot find your connected component: " + str(connCompID))
        print(connData)
        pdb.set_trace()

    # return back the phase value
    return unwPhase


def gdalTest(file):
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
        print(file + " is Not GDAL compatible")
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
            print(filetest + " is GDAL compatible")
            return filetest
        else:
            print(file + " is GDAL compatible")
            return file
    except:
        print(file + " is GDAL compatible")
        return file

#end class
if __name__ == "__main__":


    os.chdir('/Users/dbekaert/Documents/ARIA/Standard_product_tests/tops_stich_python_test')

    unw_files= ['NETCDF:"Iran_A072_20171117_20171111.nc":unwrappedPhase','NETCDF:"Iran_A072_20171117_20171111_part2.nc":unwrappedPhase']
    conn_files= ['NETCDF:"Iran_A072_20171117_20171111.nc":connectedComponents','NETCDF:"Iran_A072_20171117_20171111_part2.nc":connectedComponents']

    unw_files= ['unwPhase.unw.vrt','unwPhase2.unw.vrt']
    conn_files = ['connComp.vrt','connComp.vrt']

    product_stitch(unw_files,conn_files)

