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
import warnings
from datetime import datetime as dt
from shapely.geometry import box

def get_duplicates(df, threshold=80):
    def _check_duplicates(df, date12, threshold=90):
        warnings.simplefilter("ignore")
        g12 = df[df.DATE1_DATE2 == date12]
        g12 = g12.sort_values('CENT_LAT', ascending=False)

        # Loop through frames
        remove_list = [] 
        for ix in np.arange(g12.PATH.count() - 1):
            intersect = g12.iloc[ix+1].geometry.intersection(g12.iloc[ix].geometry)
            overlap1 = (intersect.area / g12.iloc[ix+1].geometry.area) * 100
            overlap2 = (intersect.area / g12.iloc[ix].geometry.area) * 100
            
            if (overlap1 > threshold) or (overlap2 > threshold):
                if g12.iloc[ix+1].geometry.area > g12.iloc[ix].geometry.area:
                    remove_list.append(g12.PATH.iloc[ix])
                else:
                    remove_list.append(g12.PATH.iloc[ix+1])
        return remove_list

    # Loop through aquisition dates to check for duplicates
    dates = df.groupby('DATE1_DATE2', group_keys=True).apply(lambda x: x).index.levels[0]

    duplicates = [_check_duplicates(df, date12, threshold=threshold) for date12 in dates]
    duplicates12 = []
    for duplicate in duplicates:
        if duplicate != []:
            duplicates12.extend(duplicate)

    return duplicates12

def get_union_extent(df):
    gunw_aoi = box(*df.unary_union.bounds)
    return gunw_aoi.bounds

def ds_get_extent(gdal_ds):
    E = gdal_ds.GetGeoTransform()[0]
    W = gdal_ds.GetGeoTransform()[0] + gdal_ds.GetGeoTransform()[1] * gdal_ds.RasterXSize
    N = gdal_ds.GetGeoTransform()[3] 
    S = gdal_ds.GetGeoTransform()[3] + gdal_ds.GetGeoTransform()[5] * gdal_ds.RasterYSize
    return [E, W, S, N]

def get_product_dict(df, date12, model=None):
    unw_list, conn_list, coh_list, amp_list = [], [], [], []
    inc_list, azi_list, look_list, set_list = [], [], [], []
    set_list, iono_list, tropow_list, tropod_list = [], [], [], []
    parB_list, perpB_list = [], []

    # Get GUNW layers
    for product in df.loc[date12].LAYERS.values:
        for layer in product:
            unw_list.append(layer) if 'unwrappedPhase' in layer else None
            conn_list.append(layer) if 'connectedComponents' in layer else None
            coh_list.append(layer) if 'coherence' in layer else None
            amp_list.append(layer) if 'amplitude' in layer else None
            inc_list.append(layer) if 'incidenceAngle' in layer else None
            azi_list.append(layer) if 'azimuthAngle' in layer else None
            look_list.append(layer) if 'lookAngle' in layer else None
            parB_list.append(layer) if 'parallelBaseline' in layer else None
            perpB_list.append(layer) if 'perpendicularBaseline' in layer else None
            set_list.append(layer) if 'reference/solidEarthTide' in layer else None
            iono_list.append(layer) if 'ionosphere/ionosphere' in layer else None
            tropow_list.append(layer) if 'reference/troposphereWet' in layer else None
            tropod_list.append(layer) if 'reference/troposphereHydrostatic' in layer else None

    name_list = df.loc[date12].DATE1_DATE2.values.tolist()

    # Get az time
    az_times = []
    for i in df.loc[date12].iterrows():
        hms = i[1].AZ_DOP0_MIDTIME.split(':')
        t = pd.Timestamp(i[1].DATE1)
        t = t.replace(hour=int(hms[0]), minute=int(hms[1]), second=int(float(hms[2])))
        az_times.append(dt.strftime(t, '%Y-%m-%dT%H:%M:%S.%f'))

    # Hardcoded values
    n = len(unw_list)
    cfreq_list = [5405000700.0] * n
    srange_start = [798980.125] * n
    srange_end = [956307.125] * n
    srange_spacing = [2.329562187194824] * n
    wavelength = [0.05546576] * n
    sceneLength = [35] * n
    missionID = ['Sentinel-1'] * n


    # make dict
    gunw_layers_dict = {'unwrappedPhase': unw_list, 'coherence':coh_list, 
                 'connectedComponents': conn_list, 'incidenceAngle': inc_list,
                 'azimuthAngle': azi_list, 'pair_name':name_list}
    
    gunw_dict = {'pair_name':name_list,
                 'azimuthZeroDopplerMidTime':az_times, 'centerFrequency':cfreq_list,
                 'slantRangeSpacing':srange_spacing, 'slantRangeStart':srange_start,
                 'slantRangeEnd':srange_end, 'missionID':missionID, 
                 'wavelength':wavelength, 'sceneLength':sceneLength}

    return gunw_dict, gunw_layers_dict 