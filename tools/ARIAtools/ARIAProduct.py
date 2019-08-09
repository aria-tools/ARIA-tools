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

from osgeo import gdal
gdal.UseExceptions()

# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Import functions
from ARIAtools.shapefile_util import open_shapefile,save_shapefile


class ARIA_standardproduct: #Input file(s) and bbox as either list or physical shape file.
    '''
        Class which loads ARIA standard products and splits them into spatiotemporally contigeous interferograms.
    '''

    # import dependencies
    import netCDF4
    import glob

    def __init__(self, filearg, bbox=None, workdir='./', verbose=False):
        # If user wants verbose mode
        self.verbose=verbose
        # Parse through file(s)/bbox input
        self.files = []
        self.products = []
        # Track bbox file
        self.bbox_file = None
        # Pair name for layer extraction
        self.pairname = None

        ### Determine if file input is single file, a list, or wildcard.
        # If list of files
        if len([str(val) for val in filearg.split(',')])>1:
            self.files=[str(i) for i in filearg.split(',')]
            # If wildcard
            self.files=[os.path.abspath(item) for sublist in [self.glob.glob(os.path.expanduser(os.path.expandvars(i))) if '*' in i else [i] for i in self.files] for item in sublist]
        # If single file or wildcard
        else:
            # If single file
            if os.path.isfile(filearg):
                self.files=[filearg]
            # If wildcard
            else:
                self.files=self.glob.glob(os.path.expanduser(os.path.expandvars(filearg)))
            # Convert relative paths to absolute paths
            self.files=[os.path.abspath(i) for i in self.files]
        if len(self.files)==0:
            raise Exception('No file match found')
        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        ### Check if bbox input is valid list or shapefile.
        if bbox is not None:
            # If list
            if isinstance ([str(val) for val in bbox.split()], list) and not os.path.isfile(bbox):
                from shapely.geometry import Polygon
                try:
                    bbox = [float(val) for val in bbox.split()]
                except:
                    raise Exception('Cannot understand the --bbox argument. String input is incorrect or path does not exist.')
                # Use shapely to make list
                self.bbox = Polygon(np.column_stack((np.array([bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]),
                            np.array([bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]])))) #Pass lons/lats to create polygon
                # Save polygon in shapefile
                save_shapefile(os.path.join(workdir,'user_bbox.shp'), self.bbox, 'GeoJSON')
                self.bbox_file=os.path.join(workdir,'user_bbox.shp')
                print("Shapefile %s created for input user bounds."%os.path.join(workdir,'user_bbox.shp'))
            # If shapefile
            elif os.path.isfile(bbox):
                self.bbox = open_shapefile(bbox, 0, 0)                       ##SS => We should track the projection of the shapefile. i.e. if user provides this in e.g. UTM etc.
                self.bbox_file = bbox
            else:
                raise Exception('bbox input neither valid list nor file')
        else:
            self.bbox=None

        ### Report dictionaries for all valid products
        self.__run__()


    def __readproduct__(self, file):
        '''
            Read product, determine expected layer names based off of version number, and populate corresponding product dictionary accordingly.
        '''

        ### Get standard product version from file
        # If netcdf with groups
        try:
            version=str(gdal.Open(file).GetMetadataItem('NC_GLOBAL#version'))
        except:
            print ('{} is not a supported file type... skipping'.format(file))
            return []

        # If netcdf with nogroups
        if version==str(None):
            version=str(gdal.Open(file).GetMetadataItem('version'))

        ### Get lists of radarmetadata/layer keys for this file version
        rmdkeys, sdskeys = self.__mappingVersion__(file, version)
        if self.bbox is not None:
            # Open standard product bbox
            file_bbox = open_shapefile('NETCDF:"' + file + '":'+sdskeys[0], 'productBoundingBox', 1)                    ##SS => We should track the projection of the shapefile. i.e. in case this changes in the product
            # Only generate dictionaries if there is spatial overlap with user bbox
            if file_bbox.intersects(self.bbox):
                product_dicts = [self.__mappingData__(file, rmdkeys, sdskeys)]
            else:
                product_dicts = []
        # If no bbox specified, just pass dictionaries
        else:
            product_dicts = [self.__mappingData__(file, rmdkeys, sdskeys)]

        return product_dicts


    def __mappingVersion__(self, file, version):
        '''
            Track the mapping of ARIA standard product versions.
            The order of the keys needs to be consistent with the keys in the mappingData function.
            E.g. a new expected radar-metadata key can be added as XXX to the end of the list "rmdkeys" below, and correspondingly to the end of the list "radarkeys" inside the mappingData function. Same protocol for new expected layer keys in the list "sdskeys" below, and correspondingly in "layerkeys" inside the mappingData function.
        '''

        # ARIA standard product version 1a and 1b have same mapping
        if version=='1a' or version=='1b':
            # Radarmetadata names for these versions
            rmdkeys=['missionID', 'wavelength', 'centerFrequency', 'productType',
            'ISCEversion', 'unwrapMethod', 'DEM', 'ESDthreshold', 'azimuthZeroDopplerStartTime', 'azimuthZeroDopplerEndTime',
            'azimuthTimeInterval', 'slantRangeSpacing', 'slantRangeEnd', 'slantRangeStart']

            # Layer names for these versions
            sdskeys=['productBoundingBox','unwrappedPhase','coherence',
            'connectedComponents','amplitude','perpendicularBaseline',
            'parallelBaseline','incidenceAngle','lookAngle','azimuthAngle','ionosphere']

            #Pass pair name
            read_file=self.netCDF4.Dataset(file, keepweakref=True).groups['science'].groups['radarMetaData'].groups['inputSLC']
            self.pairname=read_file.groups['reference']['L1InputGranules'][:][0][17:25] +'_'+ read_file.groups['secondary']['L1InputGranules'][:][0][17:25]
            del read_file

        return rmdkeys, sdskeys


    def __mappingData__(self, file, rmdkeys, sdskeys):
        '''
            Output and group together 2 dictionaries containing the “radarmetadata info” and “data layer keys+paths”, respectively
            The order of the dictionary keys below needs to be consistent with the keys in the __mappingVersion__ function of the ARIA_standardproduct class (see instructions on how to appropriately add new keys there).
        '''

        # Expected radarmetadata
        radarkeys=['missionID', 'wavelength', 'centerFrequency', 'productType',
        'ISCEversion', 'unwrapMethod', 'DEM', 'ESDthreshold', 'azimuthZeroDopplerStartTime', 'azimuthZeroDopplerEndTime',
        'azimuthTimeInterval', 'slantRangeSpacing', 'slantRangeEnd', 'slantRangeStart']

        # Expected layers
        layerkeys=['productBoundingBox','unwrappedPhase',
        'coherence','connectedComponents','amplitude','bPerpendicular',
        'bParallel','incidenceAngle','lookAngle',
        'azimuthAngle','ionosphere']

        # Parse radarmetadata
        rdrmetadata = self.netCDF4.Dataset(file, keepweakref=True, diskless=True).groups['science'].groups['radarMetaData']
        rdrmetakeys = list(rdrmetadata.variables.keys())
        rdrmetadata_dict={}

        # Parse layers
        sdsdict = gdal.Open(file).GetMetadata('SUBDATASETS')
        sdsdict = {k:v for k,v in sdsdict.items() if 'NAME' in k}
        datalyr_dict={}

        # Setup rdrmetadata_dict
        # for i, j in enumerate(rdrmetakeys):
        for i, j in enumerate(rdrmetakeys):
            try: #If layer expected
                rdrmetadata_dict[radarkeys[rmdkeys.index(j)]]=rdrmetadata[j][0]
            except: #If new, unaccounted layer not expected in rdrmetakeys
                print("WARNING: Radarmetadata key %s not expected in rmdkeys"%(j))
        rdrmetadata_dict['pair_name']=self.pairname

        # Setup datalyr_dict
        for i, j in enumerate(sdsdict.items()):
            #If layer expected
            try:
                datalyr_dict[layerkeys[sdskeys.index(j[1].split(':')[-1].split('/')[-1])]]=j[1]
            #If new, unaccounted layer not expected in layerkeys
            except:
                print("WARNING: Data layer key %s not expected in sdskeys"%(j[1]))
        datalyr_dict['pair_name']=self.pairname
        # 'productBoundingBox' will be updated to point to shapefile corresponding to final output raster, so record of indivdual frames preserved here
        datalyr_dict['productBoundingBoxFrames']=datalyr_dict['productBoundingBox']

        # remove temp variables
        del rdrmetadata, sdsdict

        return [rdrmetadata_dict, datalyr_dict]


    def __continuous_time__(self):
        '''
            Split the products into spatiotemporally continuous groups (i.e. by individual, continuous interferograms). Input must be already sorted by pair and start-time to fit the logic scheme below.
            Using their time-tags, this function determines whether or not successive products are in the same orbit. If in the same orbit, the program determines whether or not they overlap in time and are therefore spatially contiguous, and rejects/reports cases for which there is no temporal overlap and therefore a spatial gap.
        '''

        # import dependencies
        from datetime import datetime, timedelta
        import itertools

        sorted_products=[]
        track_rejected_pairs=[]

        # Check for (and remove) duplicate products
        for i, j in enumerate(self.products[:-1]):
            # If scenes share >90% spatial overlap AND same dates, they MUST be duplicates. Reject the latter.
            if (self.products[i+1][0]['pair_name'][9:]==j[0]['pair_name'][9:]) and (self.products[i+1][0]['pair_name'][:8]==j[0]['pair_name'][:8]) and (open_shapefile(self.products[i+1][1]['productBoundingBox'], 'productBoundingBox', 1).intersection(open_shapefile(j[1]['productBoundingBox'], 'productBoundingBox', 1)).area)/(open_shapefile(j[1]['productBoundingBox'], 'productBoundingBox', 1).area)>0.9:
                print("WARNING: Duplicate product captured. Rejecting scene %s"%(self.products[i+1][1]['unwrappedPhase'].split('"')[1]))
                # Overwrite latter scene with former
                self.products[i+1]=j
        # Delete duplicate products
        self.products=list(self.products for self.products,_ in itertools.groupby(self.products))

        # If only one pair in list, add it to list.
        if len(self.products)==1:
            sorted_products.extend([[dict(zip(self.products[0][0].keys(), [list(a) for a in zip(self.products[0][0].values())])), dict(zip(self.products[0][1].keys(), [list(a) for a in zip(self.products[0][1].values())]))]])

        # If multiple pairs in list, cycle through and evaluate temporal connectivity.
        for i, j in enumerate(self.products[:-1]):
            # Get this reference product's times
            scene_start=datetime.strptime(j[0]['azimuthZeroDopplerStartTime'], "%Y-%m-%dT%H:%M:%S.%fZ")
            scene_end=datetime.strptime(j[0]['azimuthZeroDopplerEndTime'], "%Y-%m-%dT%H:%M:%S.%fZ")
            master=datetime.strptime(j[0]['pair_name'][9:], "%Y%m%d")
            new_scene_start=datetime.strptime(self.products[i+1][0]['azimuthZeroDopplerStartTime'], "%Y-%m-%dT%H:%M:%S.%fZ")
            new_scene_end=datetime.strptime(self.products[i+1][0]['azimuthZeroDopplerEndTime'], "%Y-%m-%dT%H:%M:%S.%fZ")
            slave=datetime.strptime(self.products[i+1][0]['pair_name'][9:], "%Y%m%d")

            # Determine if next product in time is in same orbit AND overlaps AND corresponds to same pair
            # If it is within same orbit cycle, try to append scene. This accounts for day change.
            if abs(new_scene_end-scene_end)<=timedelta(minutes=100) and abs(slave-master)<=timedelta(days=1):
                # Don't export product if it is already tracked as a rejected pair
                if j[0]['pair_name'] in track_rejected_pairs or self.products[i+1][0]['pair_name'] in track_rejected_pairs:
                    track_rejected_pairs.extend((j[0]['pair_name'],self.products[i+1][0]['pair_name']))

                # Only pass scene if it temporally overlaps with reference scene
                elif ((scene_end <= new_scene_start) and (new_scene_end <= scene_start)) or ((scene_end >= new_scene_start) and (new_scene_end >= scene_start)):
                    # Check if dictionary for IFG corresponding to reference product already exists, and if it does then append values
                    try:
                        dict_ind=sorted_products.index(next(item for item in sorted_products if j[1]['productBoundingBox'] in item[1]['productBoundingBox']))
                        sorted_products[dict_ind]=[dict(zip(j[0].keys(), [[subitem for item in a for subitem in (item if isinstance(item, list) else [item])] for a in zip(sorted_products[dict_ind][0].values(), self.products[i+1][0].values())])), dict(zip(j[1].keys(), [[subitem for item in a for subitem in (item if isinstance(item, list) else [item])] for a in zip(sorted_products[dict_ind][1].values(), self.products[i+1][1].values())]))]
                    # Match IFG corresponding to reference product NOT found, so initialize dictionary for new IFG
                    except:
                        sorted_products.extend([[dict(zip(j[0].keys(), [list(a) for a in zip(j[0].values(), self.products[i+1][0].values())])), dict(zip(j[1].keys(), [list(a) for a in zip(j[1].values(), self.products[i+1][1].values())]))]])

                #Else if scene doesn't overlap, this means there is a gap. Reject date from product list, and keep track of all failed dates
                else:
                    print("Warning! Gap for interferogram %s"%(j[0]['pair_name']))
                    track_rejected_pairs.extend((j[0]['pair_name'], self.products[i+1][0]['pair_name']))

            # Products correspond to different dates, so pass both as separate IFGs.
            else:
                # Check if dictionary for IFG corresponding to reference product already exists. If not, then pass as new IFG.
                if [item for item in sorted_products if j[1]['productBoundingBox'] in item[1]['productBoundingBox']]==[] and j[0]['pair_name'] not in track_rejected_pairs:
                    sorted_products.extend([[dict(zip(j[0].keys(), [list(a) for a in zip(j[0].values())])), dict(zip(j[1].keys(), [list(a) for a in zip(j[1].values())]))]])
                # Check if dictionary for IFG corresponding to next product already exists. If not, then pass as new IFG.
                if [item for item in sorted_products if self.products[i+1][1]['productBoundingBox'] in item[1]['productBoundingBox']]==[] and self.products[i+1][0]['pair_name'] not in track_rejected_pairs:
                    sorted_products.extend([[dict(zip(self.products[i+1][0].keys(), [list(a) for a in zip(self.products[i+1][0].values())])), dict(zip(self.products[i+1][1].keys(), [list(a) for a in zip(self.products[i+1][1].values())]))]])

        # Remove duplicate dates
        track_rejected_pairs=list(set(track_rejected_pairs))
        sorted_products=[[item[0] for item in sorted_products if (item[0]['pair_name'][0] not in track_rejected_pairs)], [item[1] for item in sorted_products if (item[1]['pair_name'][0] not in track_rejected_pairs)]]

        ###Report dictionaries for all valid products
        if sorted_products==[[], []]: #Check if pairs were successfully selected
            raise Exception('No valid interferogram meet spatial criteria due to gaps and/or invalid input, nothing to export.')
        if len(track_rejected_pairs)>0:
            print("%d out of %d interferograms rejected since stitched interferogram would have gaps"%(len(track_rejected_pairs),len(sorted_products[1])+len(track_rejected_pairs)))
            # Provide report of which files were kept vs. which were not.
            if self.verbose:
                print("Specifically, the following interferograms were rejected:")
                print([item[1]['productBoundingBox'].split('"')[1] for item in sorted_products[1] if (item[1]['pair_name'][0] in track_rejected_pairs)])
        else:
            print("All (%d) interferograms are spatially continuous."%(len(sorted_products[1])))

        return sorted_products


    def __run__(self):
        # Only populate list of dictionaries if the file intersects with bbox
        for file in self.files:
            self.products += self.__readproduct__(file)
            self.products = sorted(self.products, key=lambda k: (k[0]['pair_name'], k[0]['azimuthZeroDopplerStartTime'])) #Sort by pair and start time.

        self.products=list(self.products)

        # Check if any pairs meet criteria
        if self.products==[]:
            raise Exception('No valid pairs meet spatial criteria, nothing to export.')
        if len(self.products)!=len(self.files):
            print("%d out of %d GUNW products rejected for not meeting user's bbox spatial criteria"%(len(self.files)-len(self.products),len(self.files)))
            # Provide report of which files were kept vs. which weren't
            if self.verbose:
                print("Specifically, the following GUNW products were rejected:")
                print([i for i in self.files if i not in [i[1]['productBoundingBox'].split('"')[1] for i in self.products]])
        else:
            print("All (%d) GUNW products meet spatial bbox criteria."%(len(self.files)))

        ### Split products in spatiotemporally continuous groups
        print("Group GUNW products into spatiotemporally continuous interferograms.")
        self.products = self.__continuous_time__()

        return self.products
