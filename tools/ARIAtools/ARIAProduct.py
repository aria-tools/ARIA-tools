#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import re
import numpy as np
from joblib import Parallel, delayed
import logging
from osgeo import gdal

from ARIAtools.url_manager import url_versions
from ARIAtools.shapefile_util import open_shapefile, save_shapefile
from ARIAtools.logger import logger
from ARIAtools.ARIA_global_variables import ARIA_TROPO_INTERNAL

gdal.UseExceptions()
gdal.PushErrorHandler('CPLQuietErrorHandler')
# gdal.SetConfigOption('CPL_VSIL_CURL_USE_HEAD', 'NO')

log = logging.getLogger(__name__)


# Unpacking class fuction of readproduct to become
# global fucntion that can be called in parallel
def unwrap_self_readproduct(arg):
    """

    Arg is the self argument and the filename is of the file to be read

    """
    return ARIA_standardproduct.__readproduct__(arg[0], arg[1])[0]


def package_dict(scene, new_scene, scene_ind, \
                 sorted_dict=None, dict_ind=None):
    """

    Strip and prep keys and values for dictionary of sorted, spatiotemporally contiguous products

    'scene' = specific reference product
    'new_scene' = specific secondary product (same as above if appending scenes)
    'scene_ind' = pass metadata (0) vs data layer (1) values for a given scene
    'sorted_dict' = dictionary of sorted products to modify
    'dict_ind' = index of sorted products dictionary being queried

    """

    dict_keys = scene[scene_ind].keys()
    # initialize new entry for new IFG and extend dict
    if not sorted_dict:
        if scene != new_scene:
            dict_vals = [list(a) for a in zip(scene[scene_ind].values(), \
                new_scene[scene_ind].values())]
        else:
            dict_vals = [list(a) for a in zip(scene[scene_ind].values())]
    # IFG corresponding to reference product already exists, append to dict
    if sorted_dict:
        dict_vals = [[subitem for item in a \
            for subitem in (item if isinstance(item, list) else [item])] \
                for a in zip(sorted_dict[dict_ind][scene_ind].values(), \
                    new_scene[scene_ind].values())]
    new_dict = dict(zip(dict_keys, dict_vals))

    return new_dict


# Input file(s) and bbox as either list or physical shape file.
class ARIA_standardproduct:
    """
    Load ARIA standard products

     and split them into
    spatiotemporally contiguous interferograms.
    """
    import glob

    def __init__(self, filearg, bbox=None, workdir='./', num_threads=1,
                 url_version='None', nc_version='None', verbose=False):
        """

        Parse products and input bounding box (if specified)

        """
        # If user wants verbose mode
        # Parse through file(s)/bbox input
        if verbose: logger.setLevel(logging.DEBUG)
        self.files = []
        self.products = []
        # Track bbox file
        self.bbox_file = None
        # Pair name for layer extraction
        self.pairname = None
        # enforced netcdf version
        self.nc_version = nc_version
        # pass number of threads for multiprocessing computation
        if num_threads == 'all':
            import multiprocessing
            self.num_threads = multiprocessing.cpu_count()
        else:
            self.num_threads = int(num_threads)

        # Determine if file input is single file, a list, or wildcard
        # If list of files
        if len([str(val) for val in filearg.split(',')]) > 1:
            self.files = [str(i) for i in filearg.split(',')]
            # If wildcard
            self.files = [os.path.abspath(item) for sublist in \
                [self.glob.glob(os.path.expanduser(os.path.expandvars(i))) \
                if '*' in i else [i] for i in self.files] for item in sublist]

        # If list of URLs provided
        elif os.path.basename(filearg).endswith('.txt'):
            with open(filearg, 'r') as fh:
                self.files = [f.rstrip('\n') for f in fh.readlines()]
        # If single file or wildcard
        else:
            # If single file
            if os.path.isfile(filearg):
                self.files = [filearg]
            # If wildcard
            else:
                self.files = self.glob.glob(os.path.expanduser( \
                    os.path.expandvars(filearg)))
            # Convert relative paths to absolute paths
            self.files = [os.path.abspath(i) for i in self.files]

        # capture and remove duplicate files (if applicable)
        self.files = url_versions(self.files, url_version,
                                 os.path.dirname(self.files[0]))

        # remove files that arent .nc; iterate over copy of list thats edited
        tmp_files = self.files.copy()
        for f in tmp_files:
            ext = os.path.splitext(f)[1].lower()
            if not ext == '.nc':
                self.files.remove(f)
                log.warning('%s is not a supported NetCDF... skipping', f)

        # If URLs, append with '/vsicurl/'
        self.files = [f'/vsicurl/{i}' if 'https://' in i else i for i in self.files]
        #check if virtual file reader is being captured as netcdf
        if any("https://" in i for i in self.files):
            # must configure gdal to load URLs
            gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE', 'cookies.txt')
            gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', 'cookies.txt')
            # gdal.SetConfigOption('CPL_VSIL_CURL_CHUNK_SIZE','10485760')
            gdal.SetConfigOption('VSI_CACHE', 'YES')

            fmt = gdal.Open( \
                [s for s in self.files if 'https://' in s][0] \
                ).GetDriver().GetDescription()
            if fmt != 'netCDF': raise Exception('System update required to '
                'read requested virtual products: '
                'Linux kernel >=4.3 and libnetcdf >=4.5')
        #check if local file reader is being captured as netcdf
        if any("https://" not in i for i in self.files):
            fmt = gdal.Open( \
                [s for s in self.files if 'https://' not in s][0] \
                ).GetDriver().GetDescription()
            if fmt != 'netCDF': raise Exception('System update required to '
                'read requested local products: '
                'Linux kernel >=4.3 and libnetcdf >=4.5')

        if len(self.files)==0:
            raise Exception('No file match found')
        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        ### Check if bbox input is valid list or shapefile.
        if bbox is not None:
            # If list
            if isinstance ([str(val) for val in bbox.split()], list) \
                and not os.path.isfile(bbox):
                from shapely.geometry import Polygon
                try:
                    bbox = [float(val) for val in bbox.split()]
                except:
                    raise Exception('Cannot understand the --bbox argument.' \
                        'String input is incorrect or path does not exist.')
                # Use shapely to make list
                # Pass lons/lats to create polygon
                self.bbox = Polygon(np.column_stack((np.array( \
                    [bbox[2],bbox[3],bbox[3],bbox[2],bbox[2]]), np.array(\
                    [bbox[0],bbox[0],bbox[1],bbox[1],bbox[0]]))))
                # Save polygon in shapefile
                save_shapefile(os.path.join(workdir,'user_bbox.json'),
                              self.bbox, 'GeoJSON')
                self.bbox_file=os.path.join(workdir,'user_bbox.json')
                log.info('Shapefile %s created for input user bounds',
                         os.path.join(workdir,'user_bbox.json'))
            # If shapefile
            elif os.path.isfile(bbox):
                self.bbox = open_shapefile(bbox, 0, 0)
                self.bbox_file = bbox
            else:
                raise Exception('bbox input neither valid list nor file')
        else:
            self.bbox=None

        ### Report dictionaries for all valid products
        self.__run__()


    def __readproduct__(self, fname):
        """

        Read products.

        Read product, determine expected layer names based off of version
        number, and populate corresponding product dictionary accordingly.

        """
        ### Get standard product version from file
        try:
            #version accessed differently between URL vs downloaded product
            version=str(gdal.Open(fname).GetMetadataItem('NC_GLOBAL#version'))
        except:
            log.warning('%s is not a supported file type... skipping', fname)
            return []

        # Enforce forward-compatibility of netcdf versions
        if self.nc_version == '1a': nc_version_check = ['1a', '1b', '1c']
        if self.nc_version == '1b': nc_version_check = ['1b', '1c']
        if self.nc_version == '1c': nc_version_check = ['1c']

        if version not in nc_version_check:
            log.warning('input nc_version = %s, file %s rejected because '
                        'it is a version %s product', \
                        self.nc_version, fname, version)
            return []

        ### Get lists of radarmetadata/layer keys for this file version
        fname = 'NETCDF:"' + fname
        rmdkeys, sdskeys = self.__mappingVersion__(fname, version)

        # Open standard product bbox
        if self.bbox is not None:
            file_bbox = open_shapefile(fname + '":'+sdskeys[0],
                                       'productBoundingBox', 1)
            # Only generate dictionaries if there is spatial overlap
            # with user bbox
            if file_bbox.intersects(self.bbox):
                product_dicts = [self.__mappingData__( \
                                 fname, rmdkeys, sdskeys, version)]
            else:
                product_dicts = []
        # If no bbox specified, just pass dictionaries
        else:
            product_dicts = [self.__mappingData__( \
                             fname, rmdkeys, sdskeys, version)]

        return product_dicts


    def __OGmappingVersion__(self, fname, version):
        """

        Track the mapping of ARIA standard product versions.

        The order of the keys needs to be consistent with the keys in the
        mappingData function.
        E.g. a new expected radar-metadata key can be added as XXX to the end
        of the list "rmdkeys" below, and correspondingly to the end of the list
        "radarkeys" inside the mappingData function. Same protocol for new
        expected layer keys in the list "sdskeys" below, and correspondingly in
        "layerkeys" inside the mappingData function.

        """
        import netCDF4
        # ARIA standard product version 1a and 1b have same mapping
        if version=='1a' or version=='1b':
            # Radarmetadata names for these versions
            rmdkeys = ['missionID', 'wavelength', 'centerFrequency',
                'productType', 'ISCEversion', 'unwrapMethod', 'DEM',
                'ESDthreshold', 'azimuthZeroDopplerStartTime',
                'azimuthZeroDopplerEndTime', 'azimuthTimeInterval',
                'slantRangeSpacing', 'slantRangeEnd', 'slantRangeStart']

            # Layer names for these versions
            sdskeys=['productBoundingBox', 'unwrappedPhase', 'coherence',
                'connectedComponents', 'amplitude', 'perpendicularBaseline',
                'parallelBaseline', 'incidenceAngle', 'lookAngle',
                'azimuthAngle', 'ionosphere']

            #Pass pair name
            read_file=netCDF4.Dataset(fname, keepweakref=True \
                ).groups['science'].groups['radarMetaData'].groups['inputSLC']
            self.pairname = \
               read_file.groups['reference']['L1InputGranules'][:][0][17:25] \
               + '_'+ \
               read_file.groups['secondary']['L1InputGranules'][:][0][17:25]
            del read_file

        return rmdkeys, sdskeys


    def __OGmappingData__(self, fname, rmdkeys, sdskeys):
        """

        Track the mapping of ARIA standard product versions.

        Output and group together 2 dictionaries containing the
        “radarmetadata info” and “data layer keys+paths”, respectively
        The order of the dictionary keys below needs to be consistent with
        the keys in the __mappingVersion__ function of the ARIA_standardproduct
        class (see instructions on how to appropriately add new keys there).

        """
        import netCDF4
        # Expected radarmetadata
        radarkeys=['missionID', 'wavelength', 'centerFrequency',
        'productType', 'ISCEversion', 'unwrapMethod', 'DEM',
        'ESDthreshold', 'azimuthZeroDopplerStartTime',
        'azimuthZeroDopplerEndTime',  'azimuthTimeInterval',
        'slantRangeSpacing', 'slantRangeEnd', 'slantRangeStart']

        # Expected layers
        layerkeys=['productBoundingBox','unwrappedPhase',
        'coherence','connectedComponents','amplitude','bPerpendicular',
        'bParallel','incidenceAngle','lookAngle',
        'azimuthAngle','ionosphere']

        # Parse radarmetadata
        rdrmetadata = netCDF4.Dataset(fname, keepweakref=True, \
            diskless=True).groups['science'].groups['radarMetaData']
        rdrmetakeys = list(rdrmetadata.variables.keys())
        rdrmetadata_dict={}

        # Parse layers
        sdsdict = gdal.Open(fname).GetMetadata('SUBDATASETS')
        sdsdict = {k:v for k,v in sdsdict.items() if 'NAME' in k}
        datalyr_dict={}

        # Setup rdrmetadata_dict
        for i in rdrmetakeys:
            try: #If layer expected
                rdrmetadata_dict[radarkeys[rmdkeys.index(i)]]=rdrmetadata[i][0]
            except: #If new, unaccounted layer not expected in rdrmetakeys
                log.warning("Radarmetadata key %s not expected in rmdkeys", i)
        rdrmetadata_dict['pair_name']=self.pairname

        # Setup datalyr_dict
        for i in sdsdict.items():
            #If layer expected
            try:
                datalyr_dict[layerkeys[sdskeys.index( \
                    i[1].split(':')[-1].split('/')[-1])]] = i[1]
            #If new, unaccounted layer not expected in layerkeys
            except:
                log.warning("Data layer key %s not expected in sdskeys", i[1])
        datalyr_dict['pair_name']=self.pairname
        # 'productBoundingBox' will be updated to point to shapefile
        # corresponding to final output raster, so record of
        # individual frames preserved here
        datalyr_dict['productBoundingBoxFrames'] = \
            datalyr_dict['productBoundingBox']

        # remove temp variables
        del rdrmetadata, sdsdict

        return [rdrmetadata_dict, datalyr_dict]


    def __mappingVersion__(self, fname, version):
        """

        Track the mapping of ARIA standard product versions.

        The order of the keys needs to be consistent with the keys in the
        mappingData function.
        E.g. a new expected radar-metadata key can be added as XXX to the end
        of the list "rmdkeys" below, and correspondingly to the end of the list
        "radarkeys" inside the mappingData function. Same protocol for new
        expected layer keys in the list "sdskeys" below, and correspondingly in
        "layerkeys" inside the mappingData function.

        """
        # ARIA standard prod v1a and 1b have same mapping
        # ARIA standard prod v1c differs with inclusion of ionosphere layer
        rdrmetadata_dict={}
        if version.lower() in ['1a', '1b', '1c']:
            #Pass pair name
            basename     = os.path.basename(fname)
            self.pairname=basename.split('-')[6]

            # Radarmetadata names for these versions
            rdrmetadata_dict['pair_name'] = self.pairname
            rdrmetadata_dict['azimuthZeroDopplerMidTime'] = self.pairname[:4] \
                + '-' + self.pairname[4:6] + '-' + self.pairname[6:8] \
                + 'T' + basename.split('-')[7][:2] + ':' \
                + basename.split('-')[7][2:4] + ':' \
                + basename.split('-')[7][4:] + '.0'

            #hardcoded keys for a given sensor
            if basename.split('-')[0] == 'S1':
                rdrmetadata_dict['missionID'] = 'Sentinel-1'
                rdrmetadata_dict['productType'] = 'UNW GEO IFG'
                rdrmetadata_dict['wavelength'] = 0.05546576
                rdrmetadata_dict['centerFrequency'] = 5.4050007e+09
                rdrmetadata_dict['slantRangeSpacing'] = 2.329562187194824
                rdrmetadata_dict['slantRangeStart'] = 798980.125
                rdrmetadata_dict['slantRangeEnd'] = 956307.125
                #hardcoded key meant to gauge temporal connectivity of scenes
                # (i.e. seconds between start and end)
                rdrmetadata_dict['sceneLength'] = 35
            elif basename.split('-')[0] == 'ALOS2':
                rdrmetadata_dict['missionID'] = 'ALOS-2'
                rdrmetadata_dict['productType'] = 'UNW GEO IFG'
                rdrmetadata_dict['wavelength'] = 0.229
                rdrmetadata_dict['centerFrequency'] = 1.2364997e+09
                rdrmetadata_dict['slantRangeSpacing'] = 8.582534
                rdrmetadata_dict['slantRangeStart'] = 695397.
                rdrmetadata_dict['slantRangeEnd'] = 913840.2
                #hardcoded key meant to gauge temporal connectivity of scenes
                # (i.e. seconds between start and end)
                rdrmetadata_dict['sceneLength'] = 52
            else:
                raise Exception('Sensor %s for file %s not supported.' \
                                %(basename.split('-')[0],fname))

            # Layer names for these versions
            sdskeys = [
                'productBoundingBox',
                '/science/grids/data/unwrappedPhase',
                '/science/grids/data/coherence',
                '/science/grids/data/connectedComponents',
                '/science/grids/data/amplitude',
                '/science/grids/imagingGeometry/perpendicularBaseline',
                '/science/grids/imagingGeometry/parallelBaseline',
                '/science/grids/imagingGeometry/incidenceAngle',
                '/science/grids/imagingGeometry/lookAngle',
                '/science/grids/imagingGeometry/azimuthAngle'
            ]
            if version.lower()=='1c':
                lyr_pref = '/science/grids/corrections'
                sdskeys_addlyrs = [
                    lyr_pref + '/derived/ionosphere/ionosphere',
                    lyr_pref + '/external/tides/solidEarth'
                        '/reference/solidEarthTide'
                ]
                # get weather model name(s)
                meta = gdal.Info(fname)
                model_name = []
                for i in meta.split():
                    if '/science/grids/corrections/external/troposphere/' in i:
                        model_name.append(i.split('/')[-3])
                model_name = list(set(model_name))
                for i in model_name:
                    sdskeys_addlyrs.append(lyr_pref +
                        f'/external/troposphere/{i}/reference/troposphereWet')
                    sdskeys_addlyrs.append(lyr_pref +
                        f'/external/troposphere/{i}/reference/'
                        'troposphereHydrostatic')
                # remove keys not found in product
                sdskeys_addlyrs = [i for i in sdskeys_addlyrs if i in meta]
                sdskeys.extend(sdskeys_addlyrs)

        return rdrmetadata_dict, sdskeys


    def __mappingData__(self, fname, rdrmetadata_dict, sdskeys, version):
        """

        Pass product record of metadata and layers

        Output and group together 2 dictionaries containing the
        “radarmetadata info” and “data layer keys+paths”, respectively
        The order of the dictionary keys below needs to be consistent with the
        keys in the __mappingVersion__ function of the ARIA_standardproduct
        class (see instructions on how to appropriately add new
        keys there).

        """
        # Expected layers
        layerkeys=['productBoundingBox','unwrappedPhase',
        'coherence','connectedComponents','amplitude','bPerpendicular',
        'bParallel','incidenceAngle','lookAngle','azimuthAngle']
        if version.lower()=='1c':
            # remove references to keys not found in product
            addkeys = ['ionosphere', 'solidEarthTide']
            keys_reject = [i for i in addkeys \
                                    if i not in ''.join(sdskeys)]
            addkeys = [i for i in addkeys if i not in keys_reject]

            # check for tropo layers for each model
            tropo_lyrs = ['troposphereWet', 'troposphereHydrostatic']
            for i in ARIA_TROPO_INTERNAL:
                if i in ''.join(sdskeys):
                    addkeys.append(f'{tropo_lyrs[0]}_' + i)
                    addkeys.append(f'{tropo_lyrs[1]}_' + i)
            keys_reject.extend([i for i in tropo_lyrs \
                                    if i not in ''.join(addkeys)])
            for i in keys_reject:
                log.warning(f'Expected data layer key {i} '
                    f'not found in {fname}')

            layerkeys.extend(addkeys)

        # Setup datalyr_dict
        datalyr_dict={}
        for i in enumerate(layerkeys):
            datalyr_dict[i[1]]=fname + '":'+sdskeys[i[0]]
        datalyr_dict['pair_name']=self.pairname
        # 'productBoundingBox' will be updated to point to shapefile
        # corresponding to final output raster, so record of
        # indivdual frames preserved here
        datalyr_dict['productBoundingBoxFrames'] = \
            datalyr_dict['productBoundingBox']

        return [rdrmetadata_dict, datalyr_dict]


    def __continuous_time__(self):
        """

        Split the products into spatiotemporally continuous groups.

        Split products by individual, continuous interferograms.
        Input must be already sorted by pair and start-time to fit
        the logic scheme below.
        Using their time-tags, this function determines whether or not
        successive products are in the same orbit.
        If in the same orbit, the program determines whether or not they
        overlap in time and are therefore spatially contiguous,
        and rejects/reports cases for which there is no temporal overlap
        and therefore a spatial gap.

        """
        # import dependencies
        from datetime import datetime, timedelta
        import itertools

        sorted_products = []
        track_rejected_pairs = []

        # Check for (and remove) duplicate products
        num_prods = len(self.products)
        num_dups = []
        for i in enumerate(self.products[:-1]):
            scene = i[1]
            new_scene = self.products[i[0]+1]
            # If scenes share >90% spatial overlap AND same dates
            # they MUST be duplicates. Reject the latter.
            same_time = (new_scene[0]['pair_name'][9:] == \
                                      scene[0]['pair_name'][9:]) and \
                                  (new_scene[0]['pair_name'][:8] == \
                                      scene[0]['pair_name'][:8])
            scene_shape = open_shapefile( \
                scene[1]['productBoundingBox'], 'productBoundingBox', 1)
            new_scene_shape = open_shapefile( \
                new_scene[1]['productBoundingBox'], 'productBoundingBox', 1)
            same_area = new_scene_shape.intersection(scene_shape).area \
                                      / scene_shape.area
            inv_same_area = scene_shape.intersection(new_scene_shape).area \
                                      / new_scene_shape.area
            if same_time and same_area > 0.9:
                # If applicable, overwrite smaller scene with larger one
                if same_area > inv_same_area:
                    i = (i[0]-1, new_scene)

                # Try to use newer version
                # else overwrite latter scene with former (argmax=0 when no max)
                vers   = []
                scenes = [scene, new_scene]
                for sc in scenes:
                    path_bbox = sc[1]['productBoundingBox']
                    ver_str   = re.search(r'(v\d+_\d+_\d+.*)\.', os.path.basename(path_bbox)).group(1)
                    ver_num   = float(ver_str[1:].replace('_', ''))
                    vers.append(ver_num)

                use_scene = scenes[np.argmax(vers)]
                scenes.remove(use_scene) # now scenes has the rejected scene
                self.products[i[0]+1] = use_scene
                num_dups.append(i[0])

                log.debug("Duplicate product captured. Rejecting scene %s",
                    os.path.basename(scenes[0][1]['unwrappedPhase'].split(':')[1]))

        # Delete duplicate products
        self.products=list(self.products for self.products,_ in \
                           itertools.groupby(self.products))
        if num_dups:
            log.warning("%d products rejected since they are duplicates",
                        len(num_dups))

        # If only one pair in list, add it to list.
        if len(self.products)==1:
            scene = self.products[0]
            dict_1 = package_dict(scene, scene, 0)
            dict_2 = package_dict(scene, scene, 1)
            new_dict = [dict_1, dict_2]
            sorted_products.extend([new_dict])

        # If multiple pairs in list
        # cycle through and evaluate temporal connectivity
        for i in enumerate(self.products[:-1]):
            scene = i[1]
            new_scene = self.products[i[0]+1]
            scene_t_ref = datetime.strptime(scene[0]['pair_name'][:8],
                "%Y%m%d")
            new_scene_t_ref = datetime.strptime(new_scene[0]['pair_name'][:8],
                "%Y%m%d")
            scene_t = datetime.strptime(scene[0]['pair_name'][9:],
                "%Y%m%d")
            new_scene_t = datetime.strptime(new_scene[0]['pair_name'][9:],
                "%Y%m%d")
            scene_area = open_shapefile( \
                                       scene[1]['productBoundingBox'],  \
                                       'productBoundingBox', 1)
            new_scene_area = open_shapefile( \
                                               new_scene[1]['productBoundingBox'], \
                                               'productBoundingBox', 1)

            # Only pass scene if it temporally (i.e. in same orbit)
            # and spatially overlaps with reference scene
            if scene_area.intersection(new_scene_area).area > 0. and \
                 abs(new_scene_t-scene_t) <= timedelta(days=1) and \
                 abs(new_scene_t_ref-scene_t_ref) <= timedelta(days=1):
                # Do not export prod if already tracked as a rejected pair
                if scene[0]['pair_name'] in track_rejected_pairs or \
                     new_scene[0]['pair_name'] in track_rejected_pairs:
                    track_rejected_pairs.extend((scene[0]['pair_name'], \
                        new_scene[0]['pair_name']))
                    continue
                # Check if IFG dict corresponding to ref prod already exists
                # and if it does then append values
                try:
                    dict_ind = sorted_products.index(next(\
                                     item for item in sorted_products \
                                     if scene[1]['productBoundingBox'] in \
                                     item[1]['productBoundingBox']))
                    dict_1 = package_dict(scene, new_scene, 0, \
                                     sorted_products, dict_ind)
                    dict_2 = package_dict(scene, new_scene, 1, \
                                     sorted_products, dict_ind)
                    sorted_products[dict_ind] = [dict_1, dict_2]
                # Match IFG corresponding to reference product NOT found
                # so initialize dictionary for new IFG
                except:
                    dict_1 = package_dict(scene, new_scene, 0)
                    dict_2 = package_dict(scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

            # If pairs are within the same day
            # but do not intersect this means there is a gap
            # Reject date from prod list, and keep track of all failed dates
            elif scene_area.intersection(new_scene_area).area == 0. and \
                 abs(new_scene_t-scene_t) <= timedelta(days=1) and \
                 abs(new_scene_t_ref-scene_t_ref) <= timedelta(days=1):
                track_rejected_pairs.extend((scene[0]['pair_name'], \
                    new_scene[0]['pair_name']))
                log.debug("Gap for interferogram %s \n", scene[0]['pair_name'])

            # If prods correspond to different orbits entirely
            else:
                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                if [item for item in sorted_products \
                        if scene[1]['productBoundingBox'] in \
                            item[1]['productBoundingBox']]==[] \
                        and scene[0]['pair_name'] not in \
                            track_rejected_pairs:
                    dict_1 = package_dict(scene, scene, 0)
                    dict_2 = package_dict(scene, scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])
                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                if [item for item in sorted_products \
                        if new_scene[1]['productBoundingBox'] in \
                            item[1]['productBoundingBox']]==[] \
                        and new_scene[0]['pair_name'] not in \
                            track_rejected_pairs:
                    dict_1 = package_dict(new_scene, new_scene, 0)
                    dict_2 = package_dict(new_scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

        # Remove duplicate dates
        track_rejected_pairs=list(set(track_rejected_pairs))
        if len(track_rejected_pairs)>0:
            log.warning('%d out of %d interferograms rejected since '
                        'stitched interferogram would have gaps', \
                        len(track_rejected_pairs), \
                        len([item[0] for item in sorted_products]))
            # Provide report of which files were kept vs. which were not
            log.debug('Specifically, gaps were found between the '
                        'following interferograms:')
            record_rejected_scenes = []
            for item in self.products:
                if item[0]['pair_name'] in track_rejected_pairs:
                    record_rejected_scenes.append( \
                         item[1]['productBoundingBox'].split('"')[1])
            record_rejected_scenes = list(set(record_rejected_scenes))
            record_rejected_scenes = [os.path.basename(i) \
                 for i in record_rejected_scenes]
            for i in record_rejected_scenes:
                log.debug(i)
        else:
            log.info('All (%d) interferograms are spatially continuous.', \
                     len(sorted_products))

        sorted_products=[[item[0] for item in sorted_products \
                              if (item[0]['pair_name'][0] \
                              not in track_rejected_pairs)], \
                         [item[1] for item in sorted_products \
                              if (item[1]['pair_name'][0] \
                              not in track_rejected_pairs)]]

        ###Report dictionaries for all valid products
        if sorted_products==[[], []]: #Check if pairs successfully selected
            raise Exception('No valid interferogram meet spatial criteria '
                            'due to gaps and/or invalid input, '
                            'nothing to export.')

        return sorted_products


    def __run__(self):
        # Only populate list of dictionaries if the file intersects with bbox
        # will try multi-core version (if multiple files are passed)
        # and default to for loop in case of failure

        if len(self.files)>1:
            if self.num_threads > 1:
                try:
                    log.info('Multi-core version')
                    # would probably be better not to write to the same self,
                    # and concatenate at completion.
                    self.products += Parallel(n_jobs= -1, max_nbytes=1e6)(
                                delayed(unwrap_self_readproduct)(i) for i in
                                        zip([self]*len(self.files), self.files))
                except Exception as E:
                    log.warning('Multi-core version failed with error: %s', E)
                    log.info('Will try single loop')

                    for f in self.files:
                        self.products += self.__readproduct__(f)
            else:
                for f in self.files:
                    self.products += self.__readproduct__(f)
        else:
            self.products += self.__readproduct__(self.files[0])

        # Sort by pair and start time.
        self.products = [i for i in self.products if i != []]
        self.products = sorted(self.products, key=lambda k:
                    (k[0]['pair_name'], k[0]['azimuthZeroDopplerMidTime']))
        self.products = list(self.products)

        # Exit if products from different sensors were mixed
        if not all(i[0]['missionID'] == 'Sentinel-1' for i in self.products) \
            and not all(i[0]['missionID'] == 'ALOS-2' for i in self.products):
            raise Exception('Specified input contains standard products from '
                'different sensors, please proceed with homogeneous products')

        # Check if any pairs meet criteria
        if self.products==[]:
            raise Exception('No valid pairs meet spatial criteria, nothing '
                'to export.')
        if len(self.products)!=len(self.files):
            log.warning('%d out of %d GUNW products rejected', \
                len(self.files)-len(self.products), len(self.files))
            # Provide report of which files were kept vs. which weren't
            log.debug('Specifically, the following GUNW products '
                'were rejected:')
            log.debug([os.path.basename(i) for i in self.files if i not in \
                [i[1]['productBoundingBox'].split('"')[1] \
                for i in self.products]])
        else:
            log.info('All (%d) GUNW products meet spatial bbox criteria.', \
                len(self.files))

        ### Split products in spatiotemporally continuous groups
        log.info('Group GUNW products into spatiotemporally '
                'continuous interferograms.')
        self.products = self.__continuous_time__()

        return self.products