#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Marin Govorcin

import os
import fiona
import shutil
import pickle
import warnings
import geopandas as gpd
import numpy as np
from pathlib import Path
from shapely.geometry import Polygon, box, mapping
from shapely import intersection_all
from .product.dataframe import get_gunws_df, get_df_d12stats, get_unioned_df, get_df_date12
from .product.utils import get_duplicates, get_union_extent, ds_get_extent, get_product_dict 
from .tsSetup_dask import exportUnwrappedPhase, exportCoherenceAmplitude, exportImagingGeometry

# ARIAtools
from ARIAtools.extractProduct import prep_dem
from ARIAtools.mask_util import prep_mask
from ARIAtools.tsSetup import generate_stack
from ARIAtools.ARIA_global_variables import ARIA_STACK_OUTFILES, ARIA_STACK_DEFAULTS

# NOTE: multiprocessing or joblib sometime hangs 
# TODO: https://stackoverflow.com/questions/70465276/multiprocessing-hanging-at-join

class ARIA_product():
    def __init__(self, work_dir: str, gunw_dir:str = None):
        # Initialize folder structure
        work_directory = Path(work_dir)
        work_directory.mkdir(parents=True, exist_ok=True)
        if gunw_dir is None:
            products_dir = work_directory / 'products'
            self.aria_dir = work_directory
        else:
            self.aria_dir = work_directory / gunw_dir
            gunw_dir += '/products'
            products_dir = work_directory / gunw_dir
        products_dir.mkdir(parents=True, exist_ok=True)

        # Initialize class keys
        self.work_dir = work_directory
        self.product_dir = products_dir
        self.aoi = None
        self.user_json = str(self.aria_dir / 'user_bbox.json')
        product_json_file = 'prods_TOTbbox_metadatalyr.json'
        self.product_json = str(self.aria_dir/ product_json_file)

        # Dataframe
        self.df = None # all 
        self.dataframe = None
        self.dataframe_filt = None
        self.dataframe_fin = None
        self.df_date12 = None
        self.df_duplicates = None
        self.df_rejected_disconnected = None
        self.df_rejected_aoi_coverage = None
        self.proj = None
        self.arres = None

        #Aria products
        self.products = []
        self.product_dict = []
        self.files = []
        self.dem = None
        self.dem_extent = None
        self.mask = None


    def load_gunws(self, n_jobs=10, 
                   duplicate_thresh=60,
                   verbose=True, 
                   overwrite=False):
        # duplicate threshold, flag duplicate if adjacent frames have coverage
        # more than duplicate treshold as area percentage
        df = get_gunws_df(self.product_dir, n_jobs, verbose, overwrite)
        df_gunw_dir = Path(df.PATH.iloc[0]).parent
        if self.product_dir != df_gunw_dir:
            print('GUNW directory is different than in the dataframe, reload!')
            print(f' GUNW dir: {self.product_dir} !=  Dataframe {df_gunw_dir}')
            df = get_gunws_df(self.product_dir, n_jobs, verbose, overwrite=True)

        if verbose is True:
            print('Get duplicates!')
            duplicates = get_duplicates(df, threshold=duplicate_thresh)
            print(f' Found {len(duplicates)} duplicates!')
        
        self.df = df
        self.df_duplicates = df[df.PATH.isin(duplicates)]
        self.dataframe = df[~df.PATH.isin(duplicates)]
        self.proj = df.PROJ.iloc[0]
        # NOTE; this should be better to set to fixed values
        #       same as Mintpy [-0.000833334, 0.000833334]
        #       for 90 meters
        #self.arres = [-round(np.mean(df.LAT_SPACING.values), 15), 
        #               round(np.mean(df.LON_SPACING.values), 15)]
        self.arres = [0.000833334, 0.000833334] 

    def find_aoi_intersection(self, south, north):
        # Creat aoi polygon
        if self.df_date12 is None:
            self.df_date12 = get_df_date12(self)
        extent = get_union_extent(self.df_date12)
        aoi = box(extent[0]-0.1, south, extent[2]+0.1, north)

        # Convert to GeoDataframe
        aoi_gdf = gpd.GeoDataFrame([1], geometry=[aoi],
                                   crs=self.df_date12.crs)
        
        # Find frame that intersect with aoi box
        intersection = self.dataframe.intersects(aoi_gdf.unary_union)
        self.dataframe_filt = self.dataframe[intersection]
        self.aoi = aoi
    
    def get_df_date12(self, df):
        df12 = get_df_date12(df)
        self.df_date12 = df12
        return df12 
    
    def filter_min_aoi_coverage(self, min_coverage_thresh=70, verbose=True):
        warnings.simplefilter("ignore")
        if self.aoi is None:
            raise ValueError('AOI does not exist, define it!')
        # Find union geometries per pair
        unioned_gdf = get_unioned_df(self.dataframe_filt)
        rejected = unioned_gdf[unioned_gdf['geometry'].geom_type == 'MultiPolygon'] 
        selected = unioned_gdf[unioned_gdf['geometry'].geom_type != 'MultiPolygon']
        flag = self.dataframe_filt.DATE1_DATE2.isin(rejected.DATE1_DATE2)
        self.df_rejected_disconnected = self.dataframe_filt[flag] 

        if verbose is True:
            print(f'Disconnected pairs along track')
            print(f'    Number of kept scenes: {selected.shape[0]}')
            print(f'    Number of rejected scenes: {rejected.shape[0]}')


        # Clip to aoi
        gdf_date12_clipped = selected.clip(self.aoi)
        # Get norm of all unioned pairs areas
        norm = gdf_date12_clipped.area / gdf_date12_clipped.area.max()
        min_coverage_thresh /= 100
        gdf_selected = gdf_date12_clipped[norm >= min_coverage_thresh]
        gdf_rejected = gdf_date12_clipped[norm < min_coverage_thresh]
        reject_filt = self.dataframe_filt.DATE1_DATE2.isin(gdf_rejected.DATE1_DATE2)
        select_filt = self.dataframe_filt.DATE1_DATE2.isin(gdf_selected.DATE1_DATE2) 
        self.df_rejected_aoi_coverage = self.dataframe_filt[reject_filt] 
        self.dataframe_fin = self.dataframe_filt[select_filt] 

        if verbose is True:
            print(f'\nPairs not fulfilling coverage requirement')
            print(f'    Number of kept pairs: {gdf_selected.shape[0]}')
            print(f'    Number of rejected pairs: {gdf_rejected.shape[0]}')
        return gdf_selected, gdf_rejected

    def save_aria_bbox(self, overwrite=False):
        if self.dataframe_fin is None:
            raise ValueError('Final dataframe does not exist!')
        unioned_gdf = get_unioned_df(self.dataframe_fin)

        if overwrite:
            print('Overwrite aoi json files!')
            user_json = Path(self.user_json)
            product_json = Path(self.product_json) 
            if user_json.exists(): user_json.unlink()
            if product_json.exists(): product_json.unlink()

        # Get the min common area
        bbox_shp = intersection_all(unioned_gdf.geometry)

        # Write a new Shapefile
        geojson_dict = dict(index=0, geometry='Polygon')
        with fiona.open(self.user_json, 'w', 'GeoJSON', geojson_dict) as c:
            c.write({
                'geometry': mapping(bbox_shp),
            })

        with fiona.open(self.product_json, 'w', 'GeoJSON', geojson_dict) as c:
            c.write({
                'geometry': mapping(unioned_gdf.unary_union),
            })

    def generate_product_dict(self):
        scenes = self.dataframe_fin.groupby(['DATE1_DATE2'], group_keys=True)
        scenes = scenes.apply(lambda x: x).index.levels[0]
        gdf_date12 = get_df_date12(self.dataframe_fin)
        # Generate product dict
        product_dict = [get_product_dict(gdf_date12, date) for date in scenes]

        # Split it for export
        product_dict1 = [p[0] for p in product_dict]
        product_dict2 = [p[1] for p in product_dict]

        self.products = [product_dict1, product_dict2]
        self.product_dict = product_dict2
        self.files = self.dataframe_fin.PATH.to_list()

    def prepare_dem(self, n_jobs=10, dem_option=None):
        dem_filename = Path(self.aria_dir) / 'DEM/glo_90.dem'

        if dem_option is None:
            if  dem_filename.exists():
                dem_option = str(dem_filename) #'DEM/' + dem_filename.name
            else:
                dem_option = 'download'

        # Overwrite
        #dem_option = 'download' # Uncomment if you want to skip Downloading

        # Download/Load DEM & Lat/Lon arrays, providing bbox,
        # expected DEM shape, and output dir as input.
        (dem, demfile,
        Latitude, Longitude) = prep_dem(dem_option, self.user_json, 
                                        self.user_json, self.product_json, 
                                        self.proj, arrres=self.arres,
                                        workdir=self.aria_dir, outputFormat='ISCE',
                                        num_threads=n_jobs)

        self.dem = demfile
        self.Latitude = Latitude 
        self.Longitude = Longitude
        self.dem_extent = ds_get_extent(demfile)

    def prepare_watermask(self, n_jobs=10, mask_option=None): 
        
        mask_filename = Path(self.aria_dir) / 'mask/watermask.msk'
        
        if mask_option is None:
            if  mask_filename.exists():
                mask_option = str(mask_filename) #'DEM/' + dem_filename.name
            else:
                mask_option = 'download'

        #Prepare mask
        amplitude_products = []
        for d in self.product_dict:
            if 'amplitude' in d:
                for item in list(set(d['amplitude'])):
                    amplitude_products.append(item)

        mask = prep_mask(amplitude_products,
                        mask_option,
                        self.user_json,
                        self.user_json,
                        self.proj, 
                        amp_thresh=None,
                        arrres=self.arres,
                        workdir=self.aria_dir,
                        outputFormat='ISCE',
                        num_threads=n_jobs)
        self.mask = mask.GetDescription()

    
    def export_layers(self, layer='unwrappedPhase', 
                      n_jobs=10, export_bperp=False,
                      mask_conn0 = True, verbose=True):

        SUPPORTED_LYRS = ['unwrappedPhase', 'coherence', 'incidenceAngle',
                          'azimuthAngle', 'bPerpendicular']

        if layer not in SUPPORTED_LYRS:
            msg = f'Selected layer: {layer} is not suported.'
            msg += f'Options: {SUPPORTED_LYRS}!'
            raise ValueError(msg)

        # Prepare dem and water mask
        if self.mask is None:
            prepare_watermask(self, n_jobs=n_jobs)
        if self.dem is None:
            prepare_dem(self, n_jobs=n_jobs)

        if layer == 'unwrappedPhase':
            exportUnwrappedPhase(self.product_dict, 
                                 self.user_json, self.user_json,
                                 self.arres, self.aria_dir,
                                 mask_zero_component=mask_conn0, 
                                 verbose=verbose,
                                 mask=self.mask, 
                                 n_jobs=n_jobs)

        elif layer == 'coherence':
            exportCoherenceAmplitude(self.product_dict, 
                                     self.user_json,
                                     self.user_json, 
                                     self.arres, self.aria_dir,
                                     mask=self.mask, 
                                     n_threads=1, n_jobs=n_jobs)
        
        elif layer == 'incidenceAngle' or layer == 'azimuthAngle':
            exportImagingGeometry(self.product_dict[:1], 
                                  self.user_json,
                                  self.user_json,    
                                  self.dem,
                                  self.Latitude,
                                  self.Longitude, 
                                  self.aria_dir,
                                  layer=layer, 
                                  mask=self.mask, 
                                  n_threads=1,
                                  n_jobs=1)

        elif layer == 'bPerpendicular':
            ## Take a lot of RAM memory per worker, 9GB per scene
            ## Dask reports leak - functions need restructuring
            if export_bperp:
                max_jobs = len(self.product_dict)

                # Workaround solution
                for n in range(0, max_jobs, n_jobs):
                    if n + n_jobs > len(self.product_dict):
                        print('Loop:', [n, max_jobs])
                        product_subset = self.product_dict[n:max_jobs]
                    else:
                        print('Loop:', [n, n + n_jobs])
                        product_subset = self.product_dict[n:n+n_jobs]

                    exportImagingGeometry(product_subset,
                                          self.user_json,
                                          self.user_json,    
                                          self.dem,
                                          self.Latitude,
                                          self.Longitude, 
                                          self.aria_dir,
                                          layer=layer, 
                                          mask=self.mask, 
                                          n_threads=1,
                                          n_jobs=n_jobs)

            else:
                # Store local file with bperp info
                cols = ['AVG_COHERENCE', 'BPERP', 'BTEMP']
                stack_stats = self.dataframe_fin.groupby(['DATE1_DATE2'])[cols].mean()
                stack_stats = stack_stats.reset_index()
                stack_stats.to_csv(str(self.aria_dir / 'stack_stats.csv'))

    def prepare_stack(self):

        self.check_exported_products('unwrappedPhase')
        self.check_exported_products('connectedComponents')
        ref_dlist = generate_stack(self, 'unwrappedPhase',
                                  'unwrapStack', 
                                  bperp_file='stack_stats.csv', 
                                  workdir=self.aria_dir)
        stack_dict = {
            'workdir': self.aria_dir,
            'ref_dlist': ref_dlist
        }

        # prepare additional stacks for other layers
        layers = ARIA_STACK_DEFAULTS
        layers = layers.copy()
        if 'unwrappedPhase' in layers: layers.remove('unwrappedPhase')
        remove_lyrs = []
        for i in layers:
            lyr_dir = os.path.join(str(self.aria_dir), i)
            if not os.path.exists(lyr_dir):
                if i in layers:
                    remove_lyrs.append(i)
        layers = [i for i in layers if i not in remove_lyrs] 

        # Create
        for layer in layers:
            print('')
            if layer in ARIA_STACK_OUTFILES.keys():
                self.check_exported_products(layer)
                generate_stack(self,
                            layer,
                            ARIA_STACK_OUTFILES[layer],
                                **stack_dict)

    def prep_mintpy(self, execute=False):
        # Create MIntpy directory
        mintpy_dir = self.work_dir / 'MINTPY'
        mintpy_dir.mkdir(parents=True, exist_ok=True)

        cmd = f"prep_aria.py -s {self.aria_dir}/stack -d {self.aria_dir}/DEM/glo_90.dem"
        cmd += f" -i '{self.aria_dir}/incidenceAngle/*.vrt' -a '{self.aria_dir}/azimuthAngle/*.vrt'"
        cmd += f" -w {self.aria_dir}/mask/watermask.msk"

        # Run prep_aria
        print(cmd)
        if execute:
            os.chdir(str(mintpy_dir))
            os.system(cmd)

    def save2pickle(self, fname):
        pickle = Path(fname)
        pickle.mkdir(parents=True, exist_ok=True)
        file = open(str(pickle), 'w')
        pickle.dump(self, file)

    def load_pickle(self, fname):
        pickle = Path(fname)
        file = open(str(pickle), 'r') 
        self = pickle.load(file)

    def check_exported_products(self, layer):
        # Ensure the correct number of layers:
        int_list = (Path(self.aria_dir) / layer)
        int_list = list(int_list.glob('[0-9]*[0-9].vrt'))
        product_list = [p['pair_name'][0] for p in self.products[0]]
        flag = np.array([ip.name.split('.')[0] in product_list for ip in int_list])

        # remove 
        for file in np.array(int_list)[~flag]:
            file.unlink()

    def clean_aria_directories(self):
        dir_remove_list = ['azimuthAngle', 'connectedComponents',
                           'incidenceAngle', 'coherence', 'DEM',
                           'mask', 'unwrappedPhase', 'stack']
        for rdir in dir_remove_list:
            if (Path(self.aria_dir) / rdir).exists():
                print(f'Removing {rdir}') 
                shutil.rmtree(Path(self.aria_dir) / rdir)

        for dfile in ['stack_stats.csv','user_bbox.json',
                     'prods_TOTbbox_metadatalyr.json']:
            if (Path(self.aria_dir) / dfile).exists():
                print(f'Removing {dfile}')
                (Path(self.aria_dir) / dfile).unlink()
