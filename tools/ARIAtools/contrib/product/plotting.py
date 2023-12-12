#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# ARIA product utils
from .dataframe import get_df_d12stats


def plot_network(df, coh_thresh = 0.7, min_coh=None, max_coh=None, rejected_df=None):
    '''
    df - grouped by ['DATE2','DATE1']
         than take mean of ['AVG_COHERENCE', 'BPERP', 'BTEMP'] 
    '''

    group_ix = ['DATE2','DATE1']
    group_ixs = ['AVG_COHERENCE', 'BPERP', 'BTEMP']

    # Get dataframe from plotting
    scenes2export = df.groupby(group_ix)[group_ixs].mean().reset_index()
    scenes2export['BPERP'] *= -1 # Mintpy has reverse order of ref and sec scene

    # Rejected dataframe
    if rejected_df is not None:
        rejected_gdf = rejected_df.groupby(group_ix)[group_ixs].mean().reset_index()
        rejected_gdf['BPERP'] *= -1 # Mintpy has reverse order of ref and sec scene

    if min_coh is None:
        min_coh = scenes2export.AVG_COHERENCE.min()
    if max_coh is None:
        max_coh = scenes2export.AVG_COHERENCE.max()
    
    if coh_thresh > max_coh:
        coh_thresh = np.mean([min_coh, max_coh])
    
    # Colormap
    cmap_lut = 2560
    vlist = [min_coh, coh_thresh, max_coh]
    n1_ratio = (vlist[1] - vlist[0]) / (vlist[2] - vlist[0])
    n1 = np.rint(cmap_lut * n1_ratio).astype('int')
    n2 = cmap_lut - n1
    colors1 = mpl.cm.coolwarm_r(np.linspace(0.0, 0.3, n1), alpha=0.8)
    colors2 = mpl.cm.coolwarm_r(np.linspace(0.6, 1.0, n2), alpha=0.8)
    colors = np.vstack((colors1, colors2))
    cmap = mpl.colors.LinearSegmentedColormap.from_list('network_cmap', 
                                                         colors, N=cmap_lut)

    ## Plot
    fig, ax = plt.subplots(figsize=(12,9))
    cax = make_axes_locatable(ax).append_axes("right", '3%', pad='3%')
    norm = mpl.colors.Normalize(vmin=min_coh, vmax=max_coh)
    cbar = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label('Average Spatial Coherence', fontsize=12)

    # Get unique dates
    dates = np.hstack([scenes2export[ix].unique() for ix in ['DATE1', 'DATE2']])
    dates = np.unique(dates)

    # Get perpendicular Baseline
    A = np.zeros((scenes2export.count()[0], len(dates)), np.float32)
    for index, row in scenes2export.reset_index().iterrows():
        A[index,:] = (np.int32(row['DATE1'] == dates) 
                    + np.int32(row['DATE2'] == dates)*-1)
    pbaseList = np.zeros(A.shape[1] + 1, dtype=np.float32)
    pbaseList= np.linalg.lstsq(A, scenes2export.BPERP.values, rcond=None)[0]
    print(f'Number of SAR scenes: {A.shape[1]}')
    print(f'Number of GUNWs:  {A.shape[0]}')

    # plot dates
    ax.plot(dates, pbaseList, 'ko', alpha=0.55, ms=6, mfc='black', zorder=2)
    ax.set_ylabel('Perpendicular baseline [m]')
    ax.set_xlabel('Time [year]')

    # Plot pairs
    for index, row in scenes2export.iterrows():
        x = np.array([row['DATE2'], row['DATE1']])
        ix1 = np.where(row['DATE2'] == dates)[0]
        ix2 = np.where(row['DATE1'] == dates)[0]
        val_norm = (row['AVG_COHERENCE'] - min_coh) / (max_coh - min_coh)
        y = np.array([pbaseList[ix1], pbaseList[ix2]])
        ax.plot(x, y, '-', lw=1, alpha=0.8, c=cmap(val_norm), zorder=1)

    if rejected_df is not None:
        # Get unique dates
        dates = np.hstack([rejected_gdf[ix].unique() for ix in ['DATE1', 'DATE2']])
        dates = np.unique(dates)

        A = np.zeros((rejected_gdf.count()[0], len(dates)), np.float32)
        for index, row in rejected_gdf.reset_index().iterrows():
            A[index,:] = (np.int32(row['DATE1'] == dates) 
                        + np.int32(row['DATE2'] == dates)*-1)
        pbaseList = np.zeros(A.shape[1] + 1, dtype=np.float32)
        pbaseList= np.linalg.lstsq(A, rejected_gdf.BPERP.values, rcond=None)[0]

        ax.plot(dates, pbaseList, 'ro', alpha=0.55, ms=6, mfc='black', zorder=2)

        for index, row in rejected_gdf.iterrows():
            x = np.array([row['DATE2'], row['DATE1']])
            ix1 = np.where(row['DATE2'] == dates)[0]
            ix2 = np.where(row['DATE1'] == dates)[0]
            val_norm = (row['AVG_COHERENCE'] - min_coh) / (max_coh - min_coh)
            # fix the rejected pbaseList
            y = np.array([pbaseList[ix1], pbaseList[ix2]])
            #ax.plot(x, y, '--', lw=0.5, alpha=0.8, c=cmap(val_norm), zorder=1)
            ax.plot(x, y, '--', lw=0.8, alpha=0.8, c= 'black', zorder=1)


def plot_pairing(df, color='blue', ax=None):
    if ax is None:
        fig, ax = plt.subplots(1)
    ax.plot(get_df_d12stats(df).DATE2, 
            get_df_d12stats(df).BTEMP,
            'o', c=color, ms=1, alpha=0.5)

    ax.set_title('GUNW pairs')
    ax.set_ylabel('Temporal baseline [days]')

def plot_gaps(df, min_dt=6, ax=None):
    dates = np.hstack([df[ix].unique() for ix in ['DATE1', 'DATE2']])
    dates = np.unique(dates)

    if ax is None:
        fig, ax = plt.subplots(1)
    ax.plot(dates[1:], np.diff(dates).astype('timedelta64[D]'), 'o-', ms=3, zorder=1)
    ax.plot(dates[1:], np.ones(dates[1:].shape) * min_dt, '-', ms=3, label=f'{min_dt} days delta_t', zorder=2)
    ax.set_ylabel('Temporal baseline [days]')
    ax.set_title('Vizualize gaps')
    ax.legend()

def hist_stats(df, season_boxplot=False):
    # Get states per aquisition date
    scenes2export = get_df_d12stats(df)

    #plot
    if season_boxplot is True:
        fig, ax = plt.subplots(1, figsize=(12,7))
        scenes2export.groupby(['SEASON']).boxplot(column='AVG_COHERENCE', ax=ax)
    else:
        fig, ax = plt.subplots(1,3, figsize=(12,4))
        scenes2export.hist('AVG_COHERENCE', ax=ax[0], color='orange')
        scenes2export.hist('BPERP', ax=ax[1], color='green')
        scenes2export.hist('BTEMP', ax=ax[2], color='blue')
        for a, txt in zip(ax, ['[0-1]', '[m]', '[days]']): a.set_xlabel(txt)