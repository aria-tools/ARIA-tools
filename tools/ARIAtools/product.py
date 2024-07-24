# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import re
import glob
import logging
import datetime
import itertools
import osgeo

import h5py
import netCDF4
import numpy as np
import shapely.geometry
import shapely.ops
import shapely.wkt
from pyproj import Transformer

import ARIAtools.constants
import ARIAtools.util.url
import ARIAtools.util.shp

osgeo.gdal.UseExceptions()
osgeo.gdal.PushErrorHandler('CPLQuietErrorHandler')

ONE_DAY = datetime.timedelta(days=1)
LOGGER = logging.getLogger(__name__)

# Unpacking class fuction of readproduct to become
# global fucntion that can be called in parallel


def unwrap_self_readproduct(arg):
    """
    Arg is the self argument and the filename is of the file to be read
    """
    return arg[0].__readproduct__(arg[1])[0]


def package_dict(scene, new_scene, scene_ind,
                 sorted_dict=None, dict_ind=None):
    """
    Strip and prep keys and values for dictionary of sorted, spatiotemporally
    contiguous products

    'scene' = specific reference product
    'new_scene' = specific secondary product
        (same as above if appending scenes)
    'scene_ind' = pass metadata (0) vs data layer (1) values for a given scene
    'sorted_dict' = dictionary of sorted products to modify
    'dict_ind' = index of sorted products dictionary being queried
    """
    dict_keys_ref = scene[scene_ind].keys()
    if sorted_dict:
        dict_keys_ref = sorted_dict[dict_ind][scene_ind].keys()
    dict_keys_sec = new_scene[scene_ind].keys()
    # avoid dict mismatch by passing the longer list of keys
    # which corresponds to a newer product increment
    if len(dict_keys_sec) > len(dict_keys_ref):
        dict_keys = dict_keys_sec
    else:
        dict_keys = dict_keys_ref

    # initialize new entry for new IFG and extend dict
    if not sorted_dict:
        if scene != new_scene:
            # first pass empty list if key does not exist for a scene
            for i in dict_keys:
                if i not in dict_keys_ref:
                    scene[scene_ind][i] = []
                if i not in dict_keys_sec:
                    new_scene[scene_ind][i] = []
            # then merge
            dict_vals = [
                list(a) for a in zip(
                    scene[scene_ind].values(), new_scene[scene_ind].values())]
        else:
            dict_vals = [
                list(a) for a in zip(scene[scene_ind].values())]

    # IFG corresponding to reference product already exists, append to dict
    if sorted_dict:
        dict_vals = [[
            subitem for item in a for subitem in (
                item if isinstance(item, list) else [item])]
            for a in zip(
                sorted_dict[dict_ind][scene_ind].values(),
                new_scene[scene_ind].values())]

    new_dict = dict(zip(dict_keys, dict_vals))
    return new_dict


def remove_scenes(products):
    """
    For ARIA-S1 GUNW products, identify and reject legacy v2 products if
    in the presence of v3 products for a given IFG.

    'products' = list of input products
    """
    track_legacy_products = []
    sorted_products = []

    # only check ARIA-S1 GUNW products
    if os.path.basename(products[0][1]['unwrappedPhase'].split(
            '"')[1])[:2] == 'S1':
        for i in enumerate(products[:-1]):
            scene = i[1]
            new_scene = products[i[0] + 1]
            scene_t_ref = datetime.datetime.strptime(
                scene[0]['pair_name'][:8], "%Y%m%d")
            new_scene_t_ref = datetime.datetime.strptime(
                new_scene[0]['pair_name'][:8], "%Y%m%d")
            scene_t = datetime.datetime.strptime(
                scene[0]['pair_name'][9:], "%Y%m%d")
            new_scene_t = datetime.datetime.strptime(
                new_scene[0]['pair_name'][9:], "%Y%m%d")

            # check temporal overlap
            sec_within_day = abs(new_scene_t - scene_t) <= ONE_DAY
            ref_within_day = abs(new_scene_t_ref - scene_t_ref) <= ONE_DAY

            # Only pass scene if it temporally (i.e. in same orbit)
            # overlaps with reference scene
            if sec_within_day and ref_within_day:

                # Check if IFG dict corresponding to ref prod already exists
                # and if it does then append values
                dict_item = None
                for item in sorted_products:
                    scene_in_ifg = scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        dict_item = item
                        break
                if dict_item is not None:
                    dict_ind = sorted_products.index(dict_item)
                    dict_1 = package_dict(
                        scene, new_scene, 0, sorted_products, dict_ind)
                    dict_2 = package_dict(
                        scene, new_scene, 1, sorted_products, dict_ind)
                    sorted_products[dict_ind] = [dict_1, dict_2]
                # Match IFG corresponding to reference product NOT found
                # so initialize dictionary for new IFG
                else:
                    dict_1 = package_dict(scene, new_scene, 0)
                    dict_2 = package_dict(scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

            # If prods correspond to different orbits entirely
            else:

                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                track_existing_ifg = []
                for item in sorted_products:
                    scene_in_ifg = scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        track_existing_ifg.append(item)
                if track_existing_ifg == []:
                    dict_1 = package_dict(scene, scene, 0)
                    dict_2 = package_dict(scene, scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])
                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                track_existing_ifg = []
                for item in sorted_products:
                    scene_in_ifg = new_scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        track_existing_ifg.append(item)
                if track_existing_ifg == []:
                    dict_1 = package_dict(new_scene, new_scene, 0)
                    dict_2 = package_dict(new_scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

        sorted_products = [
            [item[0] for item in sorted_products],
            [item[1] for item in sorted_products]]

        # Go through each pair and check if there are any
        # legacy products to remove
        for i in enumerate(sorted_products[1]):
            # If any version 3 products exist for a given IFG
            # remove all version 2 products
            vers = []
            for unw_f in i[1]['unwrappedPhase']:
                basename_unw = os.path.basename(
                    unw_f.split('"')[1])
                ver_str = re.search(
                    r'(v\d+_\d+_\d+.*)\.',
                    basename_unw).group(1)
                ver_num = float(ver_str[1:].replace('_', ''))
                vers.append(ver_num)

            # determine if there is a mix of incompatible versions
            v2_prods = any(item < 300. for item in vers)
            v3_prods = any(item >= 300. for item in vers)

            if v2_prods and v3_prods:
                legacy_indices = [ind for ind, val in enumerate(vers)
                                  if val < 300.]
                track_legacy_products.extend([val for ind, val in enumerate(
                    i[1]['unwrappedPhase']) if ind in legacy_indices])

                # Iterate over each dictionary to remove these keys
                LOGGER.debug('The following v2 products were rejected '
                             'to ensure they are not mixed with v3 products:')
                for item in track_legacy_products:
                    LOGGER.debug(os.path.basename(item.split('"')[1]))

        # Remove legacy products in the presence of v3 products
        if track_legacy_products != []:
            filt_products = []
            for i in products:
                if i[1]['unwrappedPhase'] not in track_legacy_products:
                    filt_products.append(i)

            products = filt_products

    return products


# Input file(s) and bbox as either list or physical shape file.
class Product:
    """
    Load ARIA standard products and split them into spatiotemporally
    contiguous interferograms.
    """

    def __init__(self, filearg, bbox=None, workdir='./', num_threads=1,
                 url_version='None', nc_version='None', projection='4326',
                 verbose=False, tropo_models=None, layers=None):
        """
        Parse products and input bounding box (if specified)
        """
        # Parse through file(s)/bbox input
        self.files = []
        self.products = []

        # Track bbox file
        self.bbox_file = None

        # Pair name for layer extraction
        self.pairname = None

        # enforced netcdf version
        self.nc_version = nc_version

        # enforced projection for output rasters
        self.projection = projection

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
            self.files = [
                os.path.abspath(item) for sublist in
                [glob.glob(os.path.expanduser(os.path.expandvars(i)))
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
                self.files = glob.glob(os.path.expanduser(
                    os.path.expandvars(filearg)))

            # Convert relative paths to absolute paths
            self.files = [os.path.abspath(i) for i in self.files]

        # capture and remove duplicate files (if applicable)
        self.files = ARIAtools.util.url.url_versions(
            self.files, url_version, os.path.dirname(self.files[0]))

        # remove files that arent .nc; iterate over copy of list thats edited
        tmp_files = self.files.copy()
        for f in tmp_files:
            ext = os.path.splitext(f)[1].lower()
            if ext not in ['.nc', '.h5']:
                self.files.remove(f)
                LOGGER.warning('%s is not a supported NetCDF... skipping', f)

        # If URLs, append with '/vsicurl/'
        self.files = [
            f'/vsicurl/{i}' if 'https://' in i else i for i in self.files]

        # check if virtual file reader is being captured as netcdf
        if any("https://" in i for i in self.files):
            # must configure osgeo.gdal to load URLs
            osgeo.gdal.SetConfigOption('GDAL_HTTP_COOKIEFILE', 'cookies.txt')
            osgeo.gdal.SetConfigOption('GDAL_HTTP_COOKIEJAR', 'cookies.txt')
            osgeo.gdal.SetConfigOption('VSI_CACHE', 'YES')

            this_file = [s for s in self.files if 'https://' in s][0]
            fmt = osgeo.gdal.Open(this_file).GetDriver().GetDescription()

            if fmt != 'netCDF':
                raise Exception(
                    'System update required to read requested virtual '
                    'products: Linux kernel >=4.3 and libnetcdf >=4.5')

        # check if local file reader is being captured as netcdf
        check_for_urls = any("https://" not in i for i in self.files)
        check_for_h5 = any(".h5" in i for i in self.files)
        if check_for_urls and not check_for_h5:
            fmt = osgeo.gdal.Open(
                [s for s in self.files if 'https://' not in s][0]
            ).GetDriver().GetDescription()
            if fmt != 'netCDF':
                raise Exception('System update required to '
                                'read requested local products: '
                                'Linux kernel >=4.3 and libnetcdf >=4.5')

        if len(self.files) == 0:
            raise Exception('No file match found')

        # exit if user does not specify a valid tropo model name
        self.tropo_models = tropo_models
        fname = self.files[0]
        basename = os.path.basename(fname)
        if self.tropo_models is not None:
            if basename.startswith('S1_'):
                if isinstance(self.tropo_models, str):
                    if self.tropo_models.lower() == 'all':
                        self.tropo_models = (
                            ARIAtools.constants.ARIA_TROPO_INTERNAL
                        )
                    else:
                        self.tropo_models = list(
                            self.tropo_models.split(',')
                        )
                        self.tropo_models = [
                            i.replace(' ', '')
                            for i in self.tropo_models]
                for i in self.tropo_models:
                    if i not in ARIAtools.constants.ARIA_TROPO_INTERNAL:
                        error_msg = 'User-requested tropo model ' \
                                    '%s will not be generated ' \
                                    'as it is not one of the ' \
                                    'following valid models: %s' % \
                                    (i, ', '.join(
                                     ARIAtools.constants.ARIA_TROPO_INTERNAL)
                                    )
                        LOGGER.error(error_msg)
                        raise Exception(error_msg)

        # set variables to check for tropo extract logic
        self.tropo_extract = False
        if layers is not None:
            if 'troposphere' in layers:
                self.tropo_extract = True

        # If specified workdir doesn't exist, create it
        if not os.path.exists(workdir):
            os.mkdir(workdir)

        # find native projection, if specified
        if self.projection.lower() == 'native':
            if basename.startswith('NISAR_'):
                record_proj = []
                pol_dict = {}
                pol_dict['SV'] = 'VV'
                pol_dict['SH'] = 'HH'
                pol_dict['HHNA'] = 'HH'
                for i in self.files:
                    fname = 'NETCDF:"' + i
                    basename = os.path.basename(i)
                    file_pol = pol_dict[basename.split('_')[10]]
                    lyr_pref = '/science/LSAR/GUNW/grids/frequencyA'
                    lyr_pref += f'/unwrappedInterferogram/{file_pol}/'
                    proj_location = lyr_pref + 'projection'
                    with h5py.File(fname[8:], 'r') as hdf_gunw:
                        file_proj = int(hdf_gunw[proj_location][()])
                    record_proj.append(file_proj)
                self.projection = int(np.median(record_proj))
            else:
                self.projection = 4326

        else:
            self.projection = int(self.projection)

        # Check if bbox input is valid list or shapefile.
        if bbox is not None:
            # If list
            bbox_is_list = isinstance(
                [str(val) for val in bbox.split()], list)
            if bbox_is_list and not os.path.isfile(bbox):

                try:
                    bbox = [float(val) for val in bbox.split()]
                except ValueError:
                    raise Exception(
                        'Cannot understand the --bbox argument. String input '
                        'is incorrect or path does not exist.')

                # Use shapely to make list
                # Pass lons/lats to create polygon
                self.bbox = shapely.geometry.Polygon(np.column_stack((
                    np.array([bbox[2], bbox[3], bbox[3], bbox[2], bbox[2]]),
                    np.array([bbox[0], bbox[0], bbox[1], bbox[1], bbox[0]]))))

                # Save polygon in shapefile
                ARIAtools.util.shp.save_shp(
                    os.path.join(workdir, 'user_bbox.json'),
                    self.bbox, self.projection, drivername='GeoJSON')
                self.bbox_file = os.path.join(workdir, 'user_bbox.json')

                LOGGER.info(
                    'Shapefile %s created for input user bounds',
                    os.path.join(workdir, 'user_bbox.json'))

            # If shapefile
            elif os.path.isfile(bbox):
                self.bbox = ARIAtools.util.shp.open_shp(bbox)
                self.bbox_file = bbox

            else:
                raise Exception('bbox input neither valid list nor file')

        else:
            self.bbox = None

        # Report dictionaries for all valid products
        self.__run__()

    def __readproduct__(self, fname):
        """
        Read products.

        Read product, determine expected layer names based off of version
        number, and populate corresponding product dictionary accordingly.
        """
        # enforce NETCDF driver to access metadata
        basename = os.path.basename(fname)
        fname = 'NETCDF:"' + fname

        # Get standard product version from file
        # version accessed differently between URL vs downloaded product
        # vs NISAR and S1 GUNWs
        if basename.startswith('NISAR_'):
            version = basename.split('_')[-1][:-3]
            version = '.'.join(version)
            nc_version_check = [version]
            if not basename.endswith('_P_F_J_001.h5'):
                LOGGER.warning(
                    'input file %s rejected because it downloaded from asf '
                    'and not a supported version', fname)
                return []

        else:
            # version accessed differently between URL vs local product
            version = str(
                osgeo.gdal.Open(fname).GetMetadataItem('NC_GLOBAL#version'))
            if version == 'None':
                LOGGER.warning(
                    '%s is not a supported file type... skipping', fname)
                return []

            # Enforce forward-compatibility of netcdf versions
            if self.nc_version == '1a':
                nc_version_check = ['1a', '1b', '1c']

            if self.nc_version == '1b':
                nc_version_check = ['1b', '1c']

            if self.nc_version == '1c':
                nc_version_check = ['1c']

        if version not in nc_version_check:
            LOGGER.warning(
                'input nc_version = %s, file %s rejected because it is a '
                'version %s product', self.nc_version, fname, version)
            return []

        # Get lists of radarmetadata/layer keys for this file version
        # separate NISAR reader
        file_bbox_intersect = True
        if basename.split('_')[0] == 'NISAR':
            rmdkeys, sdskeys, file_bbox = self.__NISARmappingVersion__(
                fname, version)

        else:
            rmdkeys, sdskeys, file_bbox = self.__mappingVersion__(
                fname, version)

        # Open standard product bbox
        if self.bbox is not None:
            file_bbox_intersect = file_bbox.intersects(self.bbox)

            # Only generate dictionaries if there is spatial overlap
            # with user bbox
            if file_bbox_intersect is False:
                product_dicts = []

        # If no bbox specified, just pass dictionaries
        if file_bbox_intersect is not False:
            # separate NISAR dict convention
            if basename.split('_')[0] == 'NISAR':
                product_dicts = [
                    self.__NISARmappingData__(
                        fname, rmdkeys, sdskeys, version)
                ]

            else:
                product_dicts = [
                    self.__mappingData__(fname, rmdkeys, sdskeys, version)]

            # assign product bounding box object to dictionary
            product_dicts[0][1]['productBoundingBox'] = file_bbox

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
        # ARIA standard product version 1a and 1b have same mapping
        if version == '1a' or version == '1b':

            # Radarmetadata names for these versions
            rmdkeys = ['missionID', 'wavelength', 'centerFrequency',
                       'productType', 'ISCEversion', 'unwrapMethod', 'DEM',
                       'ESDthreshold', 'azimuthZeroDopplerStartTime',
                       'azimuthZeroDopplerEndTime', 'azimuthTimeInterval',
                       'slantRangeSpacing', 'slantRangeEnd', 'slantRangeStart']

            # Layer names for these versions
            sdskeys = [
                'productBoundingBox', 'unwrappedPhase', 'coherence',
                'connectedComponents', 'amplitude', 'perpendicularBaseline',
                'parallelBaseline', 'incidenceAngle', 'lookAngle',
                'azimuthAngle', 'ionosphere']

            # Pass pair name
            read_file = netCDF4.Dataset(
                fname, keepweakref=True).groups['science'].groups[
                'radarMetaData'].groups['inputSLC']
            self.pairname = (read_file.groups['reference'][
                'L1InputGranules'][:][0][17:25] + '_' +
                read_file.groups['secondary']['L1InputGranules'][:][0][17:25])
        return rmdkeys, sdskeys

    def __OGmappingData__(self, fname, rmdkeys, sdskeys):
        """
        Track the mapping of ARIA standard product versions.

        Output and group together 2 dictionaries containing the
        “radarmetadata info” and “data layer keys+paths”, respectively
        The order of the dictionary keys below needs to be consistent with
        the keys in the __mappingVersion__ function of the Product
        class (see instructions on how to appropriately add new keys there).
        """
        # Expected radarmetadata
        RADAR_KEYS = [
            'missionID', 'wavelength', 'centerFrequency', 'productType',
            'ISCEversion', 'unwrapMethod', 'DEM', 'ESDthreshold',
            'azimuthZeroDopplerStartTime', 'azimuthZeroDopplerEndTime',
            'azimuthTimeInterval', 'slantRangeSpacing', 'slantRangeEnd',
            'slantRangeStart']

        # Expected layers
        LAYER_KEYS = [
            'productBoundingBox', 'unwrappedPhase', 'coherence',
            'connectedComponents', 'amplitude', 'bPerpendicular', 'bParallel',
            'incidenceAngle', 'lookAngle', 'azimuthAngle', 'ionosphere']

        # Parse radarmetadata
        rdrmetadata = netCDF4.Dataset(
            fname, keepweakref=True, diskless=True).groups['science'].groups[
            'radarMetaData']
        rdrmetakeys = list(rdrmetadata.variables.keys())
        rdrmetadata_dict = {}

        # Parse layers
        sdsdict = osgeo.gdal.Open(fname).GetMetadata('SUBDATASETS')
        sdsdict = {k: v for k, v in sdsdict.items() if 'NAME' in k}
        datalyr_dict = {}

        # Setup rdrmetadata_dict
        for i in rdrmetakeys:
            # If layer expected
            try:
                rdrmetadata_dict[RADAR_KEYS[rmdkeys.index(
                    i)]] = rdrmetadata[i][0]

            # If new, unaccounted layer not expected in rdrmetakeys
            except BaseException:
                LOGGER.warning(
                    "Radarmetadata key %s not expected in rmdkeys", i)
        rdrmetadata_dict['pair_name'] = self.pairname

        # Setup datalyr_dict
        for i in sdsdict.items():
            # If layer expected
            try:
                datalyr_dict[LAYER_KEYS[sdskeys.index(
                    i[1].split(':')[-1].split('/')[-1])]] = i[1]
            # If new, unaccounted layer not expected in LAYER_KEYS
            except BaseException:
                LOGGER.warning(
                    "Data layer key %s not expected in sdskeys", i[1])

        datalyr_dict['pair_name'] = self.pairname

        # 'productBoundingBox' will be updated to point to shapefile
        # corresponding to final output raster, so record of
        # individual frames preserved here
        datalyr_dict[
            'productBoundingBoxFrames'] = datalyr_dict['productBoundingBox']
        return [rdrmetadata_dict, datalyr_dict]

    def __mappingVersion__(self, fname, version):
        """
        Track the mapping of ARIA standard product versions.

        The order of the keys needs to be consistent with the keys in the
        mappingData function.
        E.g. a new expected radar-metadata key can be added as XXX to the
        end of the list "rmdkeys" below, and correspondingly to the end of
        the list "radarkeys" inside the mappingData function. Same protocol
        for new expected layer keys in the list "sdskeys" below, and
        correspondingly in "layerkeys" inside the mappingData function.
        """
        # ARIA standard prod v1a and 1b have same mapping
        # ARIA standard prod v1c differs with inclusion of ionosphere layer
        rdrmetadata_dict = {}
        if version.lower() in ['1a', '1b', '1c']:

            # Pass pair name
            basename = os.path.basename(fname)
            self.pairname = basename.split('-')[6]

            # Radarmetadata names for these versions
            rdrmetadata_dict['pair_name'] = self.pairname
            rdrmetadata_dict['azimuthZeroDopplerMidTime'] = (
                self.pairname[:4] + '-' + self.pairname[4:6] + '-' +
                self.pairname[6:8] + 'T' + basename.split('-')[7][:2] + ':' +
                basename.split('-')[7][2:4] + ':' +
                basename.split('-')[7][4:] + '.0')

            # assign latitude to assist with sorting
            rdrmetadata_dict[
                'centerLatitude'] = basename.split('-')[8].split('_')[1]
            if rdrmetadata_dict['centerLatitude'][-1] == 'S':
                rdrmetadata_dict['centerLatitude'] = (
                    -1 * int(rdrmetadata_dict['centerLatitude'][:-1]))

            else:
                rdrmetadata_dict['centerLatitude'] = int(rdrmetadata_dict[
                    'centerLatitude'][:-1])

            # hardcoded keys for a given sensor
            rdrmetadata_dict['projection'] = self.projection
            if basename.startswith('S1'):
                rdrmetadata_dict['missionID'] = 'Sentinel-1'
                rdrmetadata_dict['productType'] = 'UNW GEO IFG'
                rdrmetadata_dict['wavelength'] = 0.05546576
                rdrmetadata_dict['centerFrequency'] = 5.4050007e+09
                rdrmetadata_dict['slantRangeSpacing'] = 2.329562187194824
                rdrmetadata_dict['slantRangeStart'] = 798980.125
                rdrmetadata_dict['slantRangeEnd'] = 956307.125
                # hardcoded key meant to gauge temporal connectivity of scenes
                # (i.e. seconds between start and end)
                rdrmetadata_dict['sceneLength'] = 35

            elif basename.startswith('ALOS2'):
                rdrmetadata_dict['missionID'] = 'ALOS-2'
                rdrmetadata_dict['productType'] = 'UNW GEO IFG'
                rdrmetadata_dict['wavelength'] = 0.229
                rdrmetadata_dict['centerFrequency'] = 1.2364997e+09
                rdrmetadata_dict['slantRangeSpacing'] = 8.582534
                rdrmetadata_dict['slantRangeStart'] = 695397.
                rdrmetadata_dict['slantRangeEnd'] = 913840.2
                # hardcoded key meant to gauge temporal connectivity of scenes
                # (i.e. seconds between start and end)
                rdrmetadata_dict['sceneLength'] = 52

            else:
                raise Exception('Sensor %s for file %s not supported.'
                                % (basename.split('-')[0], fname))

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
                '/science/grids/imagingGeometry/azimuthAngle']

            # get product bounding box
            file_bbox = ARIAtools.util.shp.open_shp(
                fname + '":' + sdskeys[0], 'productBoundingBox', 1)

            if version.lower() == '1c':
                lyr_pref = '/science/grids/corrections'
                sdskeys_addlyrs = [
                    lyr_pref + '/derived/ionosphere/ionosphere',
                    lyr_pref + '/external/tides/solidEarth'
                    '/reference/solidEarthTide']

                # get weather model name(s)
                meta = osgeo.gdal.Info(fname)
                model_name = []
                for i in meta.split():
                    if '/science/grids/corrections/external/troposphere/' in i:
                        model_name.append(i.split('/')[-3])

                # exit if user wishes to extract a tropo layer
                # but no valid tropo model name is specified by user
                if self.tropo_extract is True and self.tropo_models is None:
                    error_msg = 'User specifies extraction of tropo layer, ' \
                                'but no valid tropo model input specified ' \
                                'with the --tropo_models option'
                    LOGGER.error(error_msg)
                    raise Exception(error_msg)

                # exit if specifies a tropo model name
                # but does not explicitly specify to extract a tropo layer
                if (self.tropo_extract is False and
                        self.tropo_models is not None):
                    TROPO_OPTIONS = {'troposphereWet',
                                     'troposphereHydrostatic',
                                     'troposphereTotal'}
                    error_msg = 'User specifies tropo model with ' \
                                'the --tropo_models option, but does not ' \
                                'specify extraction of a tropo layer ' \
                                'with the --layers option, specifically ' \
                                'any of %s' % ', '.join(TROPO_OPTIONS)
                    LOGGER.error(error_msg)
                    raise Exception(error_msg)

                model_name = list(set(model_name))
                for i in model_name:
                    sdskeys_addlyrs.append(
                        lyr_pref +
                        f'/external/troposphere/{i}/reference/troposphereWet')
                    sdskeys_addlyrs.append(
                        lyr_pref + f'/external/troposphere/{i}/reference/'
                        'troposphereHydrostatic')

                # remove keys not found in product
                sdskeys_addlyrs = [i for i in sdskeys_addlyrs if i in meta]
                sdskeys.extend(sdskeys_addlyrs)
        return rdrmetadata_dict, sdskeys, file_bbox

    def __mappingData__(self, fname, rdrmetadata_dict, sdskeys, version):
        """
        Pass product record of metadata and layers

        Output and group together 2 dictionaries containing the
        “radarmetadata info” and “data layer keys+paths”, respectively
        The order of the dictionary keys below needs to be consistent with the
        keys in the __mappingVersion__ function of the Product
        class (see instructions on how to appropriately add new
        keys there).
        """
        # Expected layers
        LAYER_KEYS = [
            'productBoundingBox', 'unwrappedPhase', 'coherence',
            'connectedComponents', 'amplitude', 'bPerpendicular',
            'bParallel', 'incidenceAngle', 'lookAngle', 'azimuthAngle']

        these_layer_keys = LAYER_KEYS.copy()
        if version.lower() == '1c':

            # remove references to keys not found in product
            addkeys = ['ionosphere', 'solidEarthTide']
            keys_reject = [
                i for i in addkeys if i not in ''.join(sdskeys)]
            addkeys = [i for i in addkeys if i not in keys_reject]

            # check for tropo layers for each model
            tropo_lyrs = ['troposphereWet', 'troposphereHydrostatic']
            for i in ARIAtools.constants.ARIA_TROPO_INTERNAL:
                if i in ''.join(sdskeys):
                    addkeys.append(f'{tropo_lyrs[0]}_' + i)
                    addkeys.append(f'{tropo_lyrs[1]}_' + i)

            keys_reject.extend([
                i for i in tropo_lyrs if i not in ''.join(addkeys)])
            for i in keys_reject:
                LOGGER.warning(
                    'Expected data layer key %s not found in %s' % (i, fname))
            these_layer_keys.extend(addkeys)

        # Setup datalyr_dict
        datalyr_dict = {}
        datalyr_dict['pair_name'] = self.pairname
        # 'productBoundingBox' will be updated to point to shapefile
        # corresponding to final output raster, so record of
        # individual frames preserved here
        datalyr_dict[
            'productBoundingBoxFrames'] = fname + '":' + sdskeys[0]
        for i in enumerate(these_layer_keys):
            datalyr_dict[i[1]] = fname + '":' + sdskeys[i[0]]

        return [rdrmetadata_dict, datalyr_dict]

    def __NISARmappingVersion__(self, fname, version):
        """

        Track the mapping of NISAR ARIA standard product versions.

        The order of the keys needs to be consistent with the keys in the
        mappingData function.
        E.g. a new expected radar-metadata key can be added as XXX to the end
        of the list "rmdkeys" below, and correspondingly to the end of the list
        "radarkeys" inside the mappingData function. Same protocol for new
        expected layer keys in the list "sdskeys" below, and correspondingly in
        "layerkeys" inside the mappingData function.

        """
        # initiate variables
        rdrmetadata_dict = {}
        sdskeys = ['/science/LSAR/identification/boundingPolygon']
        # Pass pair name
        basename = os.path.basename(fname)
        self.pairname = basename.split('_')[11][:8] + '_'
        self.pairname += basename.split('_')[13][:8]

        # Get polarization
        pol_dict = {}
        pol_dict['SV'] = 'VV'
        pol_dict['SH'] = 'HH'
        pol_dict['HHNA'] = 'HH'
        file_pol = pol_dict[basename.split('_')[10]]

        # Radarmetadata names for these versions
        rdrmetadata_dict['pair_name'] = self.pairname
        # get mid azimuth time
        ref_doppler_time = datetime.datetime.strptime(basename.split('_')[11],
                                                      '%Y%m%dT%H%M%S')
        sec_doppler_time = datetime.datetime.strptime(basename.split('_')[12],
                                                      '%Y%m%dT%H%M%S')
        mid_dt = ref_doppler_time
        mid_dt += (sec_doppler_time - ref_doppler_time) / 2
        mid_datetime_str = mid_dt.strftime('%Y-%m-%dT%H:%M:%S')
        mid_datetime_str += '.0'
        rdrmetadata_dict['azimuthZeroDopplerMidTime'] = mid_datetime_str

        # assign latitude to assist with sorting
        # get product bounding box
        # and other variables
        lyr_pref = '/science/LSAR/GUNW/grids/frequencyA'
        lyr_pref += f'/unwrappedInterferogram/{file_pol}/'
        center_freq = '/science/LSAR/GUNW/grids/frequencyA/centerFrequency'
        with h5py.File(fname[8:], 'r') as hdf_gunw:
            # get bbox
            latlon_file_bbox = hdf_gunw[sdskeys[0]]
            latlon_file_bbox = latlon_file_bbox[()]
            latlon_file_bbox = shapely.wkt.loads(latlon_file_bbox)
            # get center frequency
            center_freq_var = float(hdf_gunw[center_freq][()])
            # get slant range info
            rdr_slant_range = hdf_gunw[
                '/science/LSAR/GUNW/metadata/' +
                'radarGrid/slantRange'][()].flatten()
            min_range = min(rdr_slant_range)
            max_range = max(rdr_slant_range)
            rdr_slant_range_spac = hdf_gunw['/science/LSAR/GUNW/grids/' +
                                            'frequencyA/' +
                                            'unwrappedInterferogram/' +
                                            'xCoordinateSpacing'][()]
        rdrmetadata_dict['centerFrequency'] = center_freq_var
        rdrmetadata_dict[
            'wavelength'] = 299792458 / rdrmetadata_dict['centerFrequency']
        rdrmetadata_dict['centerLatitude'] = int(latlon_file_bbox.centroid.y)
        rdrmetadata_dict['projection'] = self.projection
        pyproj_transformer = Transformer.from_crs(
            'EPSG:4326', f'EPSG:{self.projection}', always_xy=True)
        file_bbox = shapely.ops.transform(
            lambda x, y, z=None: pyproj_transformer.transform(x, y),
            latlon_file_bbox)

        # hardcoded keys
        rdrmetadata_dict['missionID'] = 'NISAR'
        rdrmetadata_dict['productType'] = 'UNW GEO IFG'
        rdrmetadata_dict['slantRangeSpacing'] = rdr_slant_range_spac
        rdrmetadata_dict['slantRangeStart'] = min_range
        rdrmetadata_dict['slantRangeEnd'] = max_range
        # hardcoded key meant to gauge temporal connectivity of scenes
        # (i.e. seconds between start and end)
        rdrmetadata_dict['sceneLength'] = 16

        # assigning full-res layer names
        sdskeys.extend([
            lyr_pref + 'unwrappedPhase',
            lyr_pref + 'coherenceMagnitude',
            lyr_pref + 'connectedComponents',
            lyr_pref + 'ionospherePhaseScreen',
            lyr_pref + 'ionospherePhaseScreenUncertainty'
        ])
        lyr_pref = '/science/LSAR/GUNW/metadata/radarGrid/'
        sdskeys.extend([
            lyr_pref + 'perpendicularBaseline',
            lyr_pref + 'parallelBaseline',
            lyr_pref + 'incidenceAngle',
            lyr_pref + 'losUnitVectorX',  # derive azimuthAngle from this
            lyr_pref + 'losUnitVectorY',  # derive azimuthAngle from this
            lyr_pref + 'elevationAngle'
        ])
        # track and add additional correction layers, if they exist
        sdskeys_addlyrs = [
            lyr_pref + 'slantRangeSolidEarthTidesPhase',
            lyr_pref + 'alongTrackSolidEarthTidesPhase',
            lyr_pref + 'hydrostaticTroposphericPhaseScreen',
            lyr_pref + 'wetTroposphericPhaseScreen'
        ]
        meta = osgeo.gdal.Info(fname)
        # remove keys not found in product
        sdskeys_addlyrs = [i for i in sdskeys_addlyrs if i in meta]
        sdskeys.extend(sdskeys_addlyrs)

        return rdrmetadata_dict, sdskeys, file_bbox

    def __NISARmappingData__(self, fname, rdrmetadata_dict, sdskeys, version):
        """

        Pass product record of metadata and layers

        Output and group together 2 dictionaries containing the
        “radarmetadata info” and “data layer keys+paths”, respectively
        The order of the dictionary keys below needs to be consistent with the
        keys in the __NISARmappingVersion__ function of the
        ARIA_standardproduct class (see instructions on how to appropriately
        add new keys there).

        """
        # Expected layers
        layerkeys = [
            'productBoundingBox', 'unwrappedPhase', 'coherence',
            'connectedComponents', 'ionospherePhaseScreen',
            'ionospherePhaseScreenUncertainty', 'bPerpendicular', 'bParallel',
            'incidenceAngle', 'losUnitVectorX', 'losUnitVectorY',
            'elevationAngle', 'slantRangeSolidEarthTidesPhase',
            'alongTrackSolidEarthTidesPhase',
            'hydrostaticTroposphericPhaseScreen', 'wetTroposphericPhaseScreen']

        # Setup datalyr_dict
        datalyr_dict = {}
        datalyr_dict['pair_name'] = self.pairname
        # 'productBoundingBox' will be updated to point to shapefile
        # corresponding to final output raster, so record of
        # individual frames preserved here
        datalyr_dict[
            'productBoundingBoxFrames'] = fname + '":' + sdskeys[0]
        for i in enumerate(layerkeys):
            datalyr_dict[i[1]] = fname + '":' + sdskeys[i[0]]

        # Rewrite tropo and iono keys
        datalyr_dict['ionosphere'] = datalyr_dict.pop(
            'ionospherePhaseScreen')
        datalyr_dict['troposphereHydrostatic'] = datalyr_dict.pop(
            'hydrostaticTroposphericPhaseScreen')
        datalyr_dict['troposphereWet'] = datalyr_dict.pop(
            'wetTroposphericPhaseScreen')

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
        sorted_products = []
        track_rejected_pairs = []

        # For ARIA-S1 GUNW, ensure version 2 and 3 products
        # are not mixed for each IFG
        self.products = remove_scenes(self.products)

        # Check for (and remove) duplicate products
        num_dups = []
        for i in enumerate(self.products[:-1]):
            scene = i[1]
            new_scene = self.products[i[0] + 1]

            # If scenes share >90% spatial overlap AND same dates
            # they MUST be duplicates. Reject the latter.
            same_time = (new_scene[0]['pair_name'][9:] ==
                         scene[0]['pair_name'][9:]) and \
                (new_scene[0]['pair_name'][:8] ==
                 scene[0]['pair_name'][:8])
            scene_shape = scene[1]['productBoundingBox']
            new_scene_shape = new_scene[1]['productBoundingBox']
            same_area = new_scene_shape.intersection(scene_shape).area \
                / scene_shape.area
            inv_same_area = scene_shape.intersection(new_scene_shape).area \
                / new_scene_shape.area

            if same_time and same_area > 0.9:
                # If applicable, overwrite smaller scene with larger one
                if same_area > inv_same_area:
                    i = (i[0] - 1, new_scene)

                # Try to use newer version
                # else overwrite latter scene with former (argmax=0 when no
                # max)
                vers = []
                scenes = [scene, new_scene]
                for sc in scenes:
                    unw_f = sc[1]['unwrappedPhase']
                    basename_unw = os.path.basename(
                        unw_f.split('"')[1])
                    # parse version from NISAR and S1 GUNWs appropriately
                    ext = os.path.splitext(basename_unw)[1].lower()
                    if ext == '.nc':
                        ver_str = re.search(
                            r'(v\d+_\d+_\d+.*)\.',
                            basename_unw).group(1)
                        ver_num = float(ver_str[1:].replace('_', ''))
                    elif ext == '.h5':
                        ver_str = basename_unw.split('_')[-1][:-3]
                        ver_num = float(ver_str)
                    else:
                        LOGGER.error(
                            "Unable to determine version of product: %s",
                            basename_unw)
                        raise Exception(
                            "Unable to determine version of product: %s" %
                            basename_unw)
                    vers.append(ver_num)

                use_scene = scenes[np.argmax(vers)]

                # now scenes has the rejected scene
                scenes.remove(use_scene)
                self.products[i[0] + 1] = use_scene
                num_dups.append(i[0])

                LOGGER.debug(
                    "Duplicate product captured. Rejecting scene %s",
                    os.path.basename(scenes[0][1]['unwrappedPhase'].split(
                        ':')[1]))

        # Delete duplicate products
        self.products = list(
            self.products for self.products, _ in
            itertools.groupby(self.products))

        if len(num_dups) > 0:
            LOGGER.warning(
                "%d products rejected since they are duplicates",
                len(num_dups))

        # If only one pair in list, add it to list.
        if len(self.products) == 1:
            scene = self.products[0]
            dict_1 = package_dict(scene, scene, 0)
            dict_2 = package_dict(scene, scene, 1)
            new_dict = [dict_1, dict_2]
            sorted_products.extend([new_dict])

        # If multiple pairs in list
        # cycle through and evaluate temporal connectivity
        for i in enumerate(self.products[:-1]):
            scene = i[1]
            new_scene = self.products[i[0] + 1]
            scene_t_ref = datetime.datetime.strptime(
                scene[0]['pair_name'][:8], "%Y%m%d")
            new_scene_t_ref = datetime.datetime.strptime(
                new_scene[0]['pair_name'][:8], "%Y%m%d")
            scene_t = datetime.datetime.strptime(
                scene[0]['pair_name'][9:], "%Y%m%d")
            new_scene_t = datetime.datetime.strptime(
                new_scene[0]['pair_name'][9:], "%Y%m%d")
            scene_area = scene[1]['productBoundingBox']
            new_scene_area = new_scene[1]['productBoundingBox']

            # check spatiotemporal overlap
            scene_intersects = scene_area.intersection(
                new_scene_area).area > 0.
            sec_within_day = abs(new_scene_t - scene_t) <= ONE_DAY
            ref_within_day = abs(new_scene_t_ref - scene_t_ref) <= ONE_DAY

            # check spatiotemporal overlap
            scene_intersects = scene_area.intersection(
                new_scene_area).area > 0.
            sec_within_day = abs(new_scene_t - scene_t) <= ONE_DAY
            ref_within_day = abs(new_scene_t_ref - scene_t_ref) <= ONE_DAY

            # Only pass scene if it temporally (i.e. in same orbit)
            # and spatially overlaps with reference scene
            if scene_intersects and sec_within_day and ref_within_day:

                # Do not export prod if already tracked as a rejected pair
                ref_rejected = scene[0][
                    'pair_name'] in track_rejected_pairs
                sec_rejected = new_scene[0][
                    'pair_name'] in track_rejected_pairs
                if ref_rejected or sec_rejected:
                    track_rejected_pairs.extend((
                        scene[0]['pair_name'], new_scene[0]['pair_name']))
                    continue

                # Check if IFG dict corresponding to ref prod already exists
                # and if it does then append values
                dict_item = None
                for item in sorted_products:
                    scene_in_ifg = scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        dict_item = item
                        break
                if dict_item is not None:
                    dict_ind = sorted_products.index(dict_item)
                    dict_1 = package_dict(
                        scene, new_scene, 0, sorted_products, dict_ind)
                    dict_2 = package_dict(
                        scene, new_scene, 1, sorted_products, dict_ind)
                    sorted_products[dict_ind] = [dict_1, dict_2]
                # Match IFG corresponding to reference product NOT found
                # so initialize dictionary for new IFG
                else:
                    dict_1 = package_dict(scene, new_scene, 0)
                    dict_2 = package_dict(scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

            # If pairs are within the same day
            # but do not intersect this means there is a gap
            # Reject date from prod list, and keep track of all failed dates
            elif (scene_area.intersection(new_scene_area).area == 0. and
                  abs(new_scene_t - scene_t) <= ONE_DAY and
                  abs(new_scene_t_ref - scene_t_ref) <= ONE_DAY):

                track_rejected_pairs.extend((
                    scene[0]['pair_name'], new_scene[0]['pair_name']))
                LOGGER.debug(
                    "Gap for interferogram %s", scene[0]['pair_name'])

            # If prods correspond to different orbits entirely
            else:

                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                track_existing_ifg = []
                for item in sorted_products:
                    scene_in_ifg = scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        track_existing_ifg.append(item)
                ref_not_rejected = scene[0][
                    'pair_name'] not in track_rejected_pairs
                if track_existing_ifg == [] and ref_not_rejected:
                    dict_1 = package_dict(scene, scene, 0)
                    dict_2 = package_dict(scene, scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])
                # Check if IFG dict corresponding to ref prod already exists
                # and if it does not then pass as new IFG
                track_existing_ifg = []
                for item in sorted_products:
                    scene_in_ifg = new_scene[1][
                        'unwrappedPhase'] in item[1]['unwrappedPhase']
                    if scene_in_ifg:
                        track_existing_ifg.append(item)
                scene_not_rejected = new_scene[0][
                    'pair_name'] not in track_rejected_pairs
                if track_existing_ifg == [] and scene_not_rejected:
                    dict_1 = package_dict(new_scene, new_scene, 0)
                    dict_2 = package_dict(new_scene, new_scene, 1)
                    new_dict = [dict_1, dict_2]
                    sorted_products.extend([new_dict])

        # Remove duplicate dates
        track_rejected_pairs = list(set(track_rejected_pairs))
        if len(track_rejected_pairs) > 0:
            LOGGER.warning(
                '%d out of %d interferograms rejected since stitched '
                'interferogram would have gaps' % (
                    len(track_rejected_pairs),
                    len([item[0] for item in sorted_products])))

            # Provide report of which files were kept vs. which were not
            LOGGER.debug('Specifically, gaps were found between the '
                         'following interferograms:')

            record_rejected_scenes = []
            for item in self.products:
                if item[0]['pair_name'] in track_rejected_pairs:
                    record_rejected_scenes.append(
                        item[1]['unwrappedPhase'].split('"')[1])

            for record_rejected_scene in list(set(record_rejected_scenes)):
                LOGGER.debug(os.path.basename(record_rejected_scene))

        else:
            LOGGER.info(
                'All (%d) interferograms are spatially continuous.',
                len(sorted_products))

        sorted_products = [
            [item[0] for item in sorted_products if (item[0]['pair_name'][0]
             not in track_rejected_pairs)],
            [item[1] for item in sorted_products if (item[1]['pair_name'][0]
             not in track_rejected_pairs)]]

        # Report dictionaries for all valid products
        # Check if pairs successfully selected
        if sorted_products == [[], []]:
            raise Exception(
                'No valid interferogram meet spatial criteria due to gaps '
                'and/or invalid input, nothing to export.')

        return sorted_products

    def __run__(self):
        # Only populate list of dictionaries if the file intersects with bbox
        for f in self.files:
            self.products += self.__readproduct__(f)

        # Sort by pair, start time, and latitude
        self.products = list(sorted(
            [i for i in self.products if i != []], key=lambda i: (
                i[0]['pair_name'], i[0]['centerLatitude'],
                i[0]['azimuthZeroDopplerMidTime'])))

        # determine if there is a mix of different sensors
        s1_prods = all(i[0]['missionID'] ==
                       'Sentinel-1' for i in self.products)
        alos2_prods = all(i[0]['missionID'] == 'ALOS-2' for i in self.products)
        nisar_prods = all(i[0]['missionID'] == 'NISAR' for i in self.products)

        # Exit if products from different sensors were mixed
        if not s1_prods and not alos2_prods and not nisar_prods:

            raise Exception(
                'Specified input contains standard products from different '
                'sensors, please proceed with homogeneous products')

        # Check if any pairs meet criteria
        if self.products == []:
            raise Exception('No valid pairs meet spatial criteria, nothing '
                            'to export.')

        if len(self.products) != len(self.files):
            LOGGER.warning(
                '%d out of %d GUNW products rejected',
                len(self.files) - len(self.products), len(self.files))

            # Provide report of which files were kept vs. which weren't
            LOGGER.debug(
                'Specifically, the following GUNW products were rejected:')
            for i in self.files:
                product_bboxes = [i[1]['unwrappedPhase'].split('"')[1]
                                  for i in self.products]
                if i not in product_bboxes:
                    LOGGER.debug(os.path.basename(i))

        else:
            LOGGER.info(
                'All (%d) GUNW products meet spatial bbox criteria.',
                len(self.files))

        # Split products in spatiotemporally continuous groups
        LOGGER.info(
            'Group GUNW products into spatiotemporally continuous '
            'interferograms.')
        self.products = self.__continuous_time__()
        return self.products
