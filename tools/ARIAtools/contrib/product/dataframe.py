#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Marin Govorcin

import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
import logging

from pathlib import Path
from tqdm import tqdm
from osgeo import gdal, ogr
from shapely.wkt import loads
from datetime import datetime as dt

# Parallel processing
from multiprocessing import Pool

# Start logging
logging.basicConfig(filename='warning.log',
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

logging.info("Running ARIA-Product")

def get_dataframe(filename, get_stats=True):
    def _get_mean_stats(filename : str, nodata: float=None) -> np.float32: 
        # Read
        ds = gdal.Open(filename)
        ds_array = ds.ReadAsArray()
        if nodata:
            ds_array = np.ma.masked_equal(ds_array, nodata)
        ds = None
        return np.mean(ds_array)
        
    def _get_boundingBox(filename :str):
        # Supress warning messages reading netcdf
        gdal.UseExceptions()
        gdal.PushErrorHandler('CPLQuietErrorHandler')
        ds = ogr.Open(filename)
        feature = ds.GetLayer().GetFeature(1)
        wkt = feature.GetGeometryRef().ExportToWkt()
        ds, feature = None, None
        return wkt

    # Initialize output
    keyList = ['SENSOR', 'ORBIT', 'TRACK', 
               'DATE1_DATE2', 'AZ_DOP0_MIDTIME',
               'AVG_COHERENCE', 'BPERP', 
               'BTEMP', 'GEOMETRY', 'VERSION', 'PATH',
               'LAT_SPACING','LON_SPACING', 'LENGTH', 'WIDTH',
               'Y_ORIGIN', 'X_ORIGIN', 'PROJ', 'LAYERS']

    out_dict = dict(zip(keyList, [None]*len(keyList))) 
    
    # NETCDF filename
    if type(filename) == str:
        filename = Path(filename)
        
    nc_file = f'NETCDF:"{str(filename)}"' 

    # Get bounding box polygon
    try:
        boundBox_layer = nc_file + ':productBoundingBox'
        out_dict['GEOMETRY'] = _get_boundingBox(boundBox_layer) 

        # Get average spatial coherence
        coh_layer = nc_file + ':/science/grids/data/coherence'
        out_dict['AVG_COHERENCE'] = _get_mean_stats(coh_layer, nodata=0)

        # Get perpendicular baseline
        bperp_layer = nc_file
        bperp_layer += ':/science/grids/imagingGeometry/perpendicularBaseline'
        out_dict['BPERP'] = _get_mean_stats(bperp_layer)

        # Get GUNW version
        unw_layer = nc_file
        unw_layer += ':/science/grids/data/unwrappedPhase' 
        ds = gdal.Open(str(unw_layer), gdal.GA_ReadOnly)
        out_dict['VERSION'] = ds.GetMetadata()['NC_GLOBAL#version']
        out_dict['LAT_SPACING'] = ds.GetGeoTransform()[5]
        out_dict['LON_SPACING'] = ds.GetGeoTransform()[1]
        out_dict['X_ORIGIN'] = ds.GetGeoTransform()[0] 
        out_dict['Y_ORIGIN'] = ds.GetGeoTransform()[3] 
        out_dict['WIDTH'] = ds.RasterXSize
        out_dict['LENGTH'] = ds.RasterYSize
        out_dict['PROJ'] = ds.GetProjectionRef()
        ds = None
        
        # update dict with information
        name = filename.name
        out_dict['SENSOR'] = name.split('-')[0] 
        out_dict['ORBIT'] = name.split('-')[2]
        out_dict['TRACK'] = name.split('-')[4]
        out_dict['DATE1_DATE2'] = name.split('-')[6]
        midtime = dt.strptime(name.split('-')[7],'%H%M%S')
        out_dict['AZ_DOP0_MIDTIME'] = midtime.strftime('%H:%M:%S.0') 
        out_dict['PATH'] = str(filename)

        # Get all layers
        ds = gdal.Info(str(filename), format='json')
        out_dict['LAYERS'] = list(ds['metadata']['SUBDATASETS'].values())[::2]
        ds = None
        return out_dict
    except Exception:
        skipped_dir = (filename.parent / 'skipped').resolve()
        skipped_dir.mkdir(parents=True, exist_ok=True)
        logging.warning(f'Error reading {filename}')
        filename.rename(skipped_dir / filename.name)
        return None


def get_gunws_df(work_dir, n_jobs=1, verbose=False, overwrite=False):
    def _run_getdf(filenames, n_jobs=1):
        # initalize multiprocessing
        pool = Pool(processes=n_jobs)
        # Prepare jobs
        out = []
        print('Generate ARIA Dataframe:')
        for result in tqdm(pool.imap(func=get_dataframe, 
                                    iterable=filenames), 
                        total=len(filenames)):
            if result != None :
                out.append(result)
        pool.close()

        # Combine dataframe in one
        print(out)
        df = pd.DataFrame(out)
        # update Dataframe with temporal baseline
        df['DATE1'] = np.vstack(df.DATE1_DATE2.apply(get_data12))[:,0]
        df['DATE2'] = np.vstack(df.DATE1_DATE2.apply(get_data12))[:,1]
        df['BTEMP'] = (df['DATE1'] - df['DATE2']).dt.days

        return df

    # verbose printing
    vprint = lambda x: print(x) if verbose == True else None

    vprint(f'GUNW directory: {work_dir}')
    gunws = list(work_dir.glob('*S1*.nc'))
    vprint(f'Number of GUNW products: {len(gunws)}')

    # Get track number and pickle filename
    track = gunws[0].name.split('-')[4]
    pickle_file = str(work_dir / f"gunws_{track}.pkl")

    # Update mode
    # Delete if overwrite is on 
    if overwrite:  Path(pickle_file).unlink(missing_ok=True) 

    if Path(pickle_file).exists():
        vprint(f'  Pickle {Path(pickle_file).name} exists.')
        df = pd.read_pickle(pickle_file)
        if df.count()[0] != len(gunws):
            vprint(f'  Run update mode, pickle: {df.count()[0]} vs dir:{len(gunws)}')
            gunws_list = list(map(str, gunws))
            missing_gunws = ~np.array([g in df.PATH.values 
                                       for g in gunws_list])
            gunws = list(compress(gunws, missing_gunws.tolist()))

            # Run extracting
            df1 = _run_getdf(gunws, n_jobs=n_jobs)
            df = pd.concat([df, df1], ignore_index=True)
            df.to_pickle(pickle_file)
    
    else:
        # Run it all
        df = _run_getdf(gunws, n_jobs=n_jobs)
        # save to pickle
        df.to_pickle(pickle_file)


    # Convert to geoDataframe
    gdf = gpd.GeoDataFrame(df, crs = "EPSG:4326",
                           geometry=df.GEOMETRY.apply(loads))

    with warnings.catch_warnings():
        warnings.simplefilter(action='ignore', category=UserWarning)
        centroid = gdf["geometry"].centroid
        gdf['CENT_LAT'] = [ix.xy[1][0] for ix in centroid]

    return gdf 

def get_data12(date: str):
    d1, d2 = date.split('_')
    d1 = np.array(dt.strftime(dt.strptime(d1, '%Y%m%d'),'%Y-%m-%d'))
    d2 = np.array(dt.strftime(dt.strptime(d2, '%Y%m%d'),'%Y-%m-%d'))
    return d1.astype('datetime64[D]'), d2.astype('datetime64[D]')

def get_df_date12(df):
    column = ['DATE1_DATE2']
    grouped = df.groupby(column, group_keys=True)
    return grouped.apply(lambda x: x)


def get_df_d12stats(df, save=None):
    scenes2export = df.groupby(['DATE2','DATE1'])[['AVG_COHERENCE', 'BPERP', 'BTEMP']].mean().reset_index()
    scenes2export['BPERP'] *= -1 # Mintpy has reverse order of ref and sec scene
    seasons = {1:'DJF', 2: 'MAM', 3:'JJA', 4:'SON'}
    scenes2export['SEASON'] = ((scenes2export.DATE2.dt.month % 12 + 3) // 3).map(seasons)
    return scenes2export

def get_unioned_df(df):
    # Find union geometries per pair
    unioned_geometries = []
    for group_name, group_data in df.groupby('DATE1_DATE2', group_keys=True):
        unioned_geometry = group_data.unary_union
        unioned_geometries.append(unioned_geometry)

    unioned_gdf = gpd.GeoDataFrame(
        {'DATE1_DATE2': df.groupby('DATE1_DATE2', group_keys=True).groups.keys(),
        'geometry': unioned_geometries},
        geometry='geometry')
    return unioned_gdf