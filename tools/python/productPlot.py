#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import sys
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Import functions
from ARIAProduct import ARIA_standardproduct
from shapefile_util import open_shapefile

def createParser():
    '''
        Make any of the following specified plot(s): ⊥ baseline + histogram, coherence + histogram + average coherence raster, ⊥ baseline & coherence combo, and track extents. The default is to generate all of these.
    '''
    import argparse

    parser = argparse.ArgumentParser(description='Function to generate various quality control and baseline figures of the spatial-temporal network of products.')
    parser.add_argument('-f', '--file', dest='imgfile', type=str,
            required=True, help='ARIA file')
    parser.add_argument('-w', '--workdir', dest='workdir', default='./', help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('-b', '--bbox', dest='bbox', type=str, default=None, help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- Example : '19 20 -99.5 -98.5'")
    parser.add_argument('-m', '--mask', dest='mask', type=str, default=None, help="Mask file. This file needs to be GDAL compatible and contains projection information its meta-data.")
    parser.add_argument('-of', '--outputFormat', dest='outputFormat', type=str, default='ENVI', help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By default files are generated with ENVI format.')
    parser.add_argument('-croptounion', '--croptounion', action='store_true', dest='croptounion', help="If turned on, IFGs cropped to bounds based off of union and bbox (if specified). Program defaults to crop all IFGs to bounds based off of common intersection and bbox (if specified).")
    parser.add_argument('-plottracks', '--plottracks', action='store_true', dest='plottracks', help="Make plot of track latitude extents vs bounding bbox/common track extent.")
    parser.add_argument('-plotbperp', '--plotbperp', action='store_true', dest='plotbperp', help="Make a baseline plot, and a histogram of perpendicular baseline.")
    parser.add_argument('-plotbperpcoh', '--plotbperpcoh', action='store_true', dest='plotbperpcoh', help="Make a baseline plot that is color-coded based on average IFG coherence.")
    parser.add_argument('-plotcoh', '--plotcoh', action='store_true', dest='plotcoh', help="Make an average IFG coherence plot in time, and histogram of IFG average coherence.")
    parser.add_argument('-makeavgoh', '--makeavgoh', action='store_true', dest='makeavgoh', help="Generate a 2D raster of average IFG coherence.")
    parser.add_argument('-plotall', '--plotall', action='store_true', dest='plotall', help="Generate all above plots.")
    parser.add_argument('-verbose', '--verbose', action='store_true', dest='verbose', help="Toggle verbose mode on.")
    return parser

def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



class plot_class:
    '''
    Class to generate standard plots for ARIA products.
    '''

    # importing dependencies
    from datetime import datetime, date
    from dateutil.relativedelta import relativedelta
    import matplotlib as mpl
#    mpl.use("TkAgg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    import pandas as pd
    from pandas.plotting import register_matplotlib_converters
    import warnings
    register_matplotlib_converters()

    def __init__(self, product_dict, workdir='./', bbox_file=None, prods_TOTbbox=None, mask=None, outputFormat='ENVI', croptounion=False):
        # Pass inputs, and initialize list of pairs
        self.product_dict = product_dict
        self.bbox_file = bbox_file
        if self.bbox_file:
            self.bbox_file = open_shapefile(bbox_file, 0, 0).bounds
        self.workdir = workdir
        self.prods_TOTbbox = prods_TOTbbox
        self.mask = mask
        self.outputFormat = outputFormat
        self.croptounion = croptounion
        # File must be physically extracted, cannot proceed with VRT format. Defaulting to ENVI format.
        if  self.outputFormat=='VRT':
            self.outputFormat='ENVI'

        # create workdir if it doesn't exist
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)

        self.pairs=None

    def __date_list__(self):
        '''
            Make dictionary of time differences between successive epochs.
        '''

        # importing dependencies
        import time

        dateList = []
        tbase = []
        # Get list of epochs
        for di in self.pairs:
            dates = os.path.basename(di).split('_')
            if not dates[1] in dateList: dateList.append(dates[1])
            if not dates[0] in dateList: dateList.append(dates[0])

        dateList.sort()
        d1 = self.datetime(*time.strptime(dateList[0],"%Y%m%d")[0:5])
        for ni in range(len(dateList)):
            d2 = self.datetime(*time.strptime(dateList[ni],"%Y%m%d")[0:5])
            diff = d2-d1
            tbase.append(diff.days)
        dateDict = {}
        for i in range(len(dateList)): dateDict[dateList[i]] = tbase[i]

        return dateDict

    def __design_matrix__(self):
        '''
            Make the design matrix for the inversion.
        '''

        dateDict = self.__date_list__()
        numDates = len(dateDict)
        numIfgrams = len(self.pairs)

        #Initialize matrices for inversion
        A = np.zeros((numIfgrams,numDates))
        B = np.zeros(np.shape(A))
        L = np.zeros((numIfgrams,1))
        daysList = []
        for day in dateDict.values():
            daysList.append(day)

        tbase = np.array(list(dateDict.values()))
        t = np.zeros((numIfgrams,2))
        baseline_hist = []

        #Iterate through all IFGs
        for i, j in enumerate(self.product_dict[0]):
            date12 = self.product_dict[1][i][0]
            date = date12.split('_')
            ndxt2 = daysList.index(dateDict[date[0]])
            ndxt1 = daysList.index(dateDict[date[1]])
            A[i,ndxt1] = -1
            A[i,ndxt2] = 1
            B[i,ndxt1:ndxt2] = tbase[ndxt1+1:ndxt2+1]-tbase[ndxt1:ndxt2]
            t[i,:] = [dateDict[date[1]],dateDict[date[0]]]

            # Report mean
            pbaseline_nodata=gdal.Open(j[0])
            pbaseline_nodata=pbaseline_nodata.GetRasterBand(1).GetNoDataValue()
            pbaseline_val=gdal.BuildVRT('', j).ReadAsArray()
            pbaseline_val=np.ma.masked_where(pbaseline_val == pbaseline_nodata, pbaseline_val)
            pbaseline_val=pbaseline_val.mean()
            # Record baseline val for histogram
            baseline_hist.append(pbaseline_val)
            L[i] = float(pbaseline_val)
            if (np.isnan(L[i])):
                L[i] = 0.0
            del pbaseline_val, pbaseline_nodata

        A = A[:,1:]
        B = B[:,:-1]

        ind=~np.isnan(L)
        return A[ind[:,0],:],B[ind[:,0],:],L[ind],baseline_hist

    def plot_pbaselines(self):
        '''
            Make baseline plot + histogram of baselines.
        '''

        ax=self.plt.figure().add_subplot(111)
        self.pairs=[i[0] for i in self.product_dict[1]]
        dateDict = self.__date_list__()
        A,B,L,baseline_hist = self.__design_matrix__()

        # Perform inversion
        B1 = np.linalg.pinv(B)
        B1 = np.array(B1,np.float32)
        dS = np.dot(B1,L)
        dtbase = np.diff(list(dateDict.values()))
        dt = np.zeros((len(dtbase),1))
        zero = np.array([0.],np.float32)

        S = np.concatenate((zero,np.cumsum([dS*dtbase])))
        residual = L-np.dot(B,dS)

        RMSE = np.sqrt(np.sum(residual**2)/len(residual))
        if np.linalg.matrix_rank(B)!=len(list(dateDict.keys()))-1:
            print('Baseline plot warning!')
            print('Design matrix is rank deficient. Network is disconnected.')
            print('Using a fully connected network is recommended.')

        offset_dict={}
        # Plot dot for each date
        for i in range(len(list(dateDict.keys()))):
            offset_dict[list(dateDict.keys())[i]]=S[i]
            master=self.pd.to_datetime(list(dateDict.keys())[i][:8])
            ax.plot(master, S[i], 'k.', markeredgewidth = 3, markersize=15, linestyle='None', zorder=10)

        # Plot lines for each pair
        for i in self.pairs:
            slave=self.pd.to_datetime(i[:8])
            master=self.pd.to_datetime(i[9:])
            ax.plot([master, slave], [offset_dict[i[9:]], offset_dict[i[:8]]], 'b')

        # Make Baseline plot
        ax.set_ylabel(r'$\perp$'+' Baseline (m)',weight='bold')
        ax.set_xlabel('Time',weight='bold')
        xticks, labels = self._adaptive_xticks(list(set(dateDict.keys())))
        ax.set_xlim(min(xticks),max(xticks))
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
        ax.xaxis.set_major_formatter(self.mpl.dates.DateFormatter('%Y-%m'))
        for label in ax.get_xticklabels():
                 label.set_ha('center')
                 label.set_rotation(20.)

        # If plotting parameters are the default, must adjust X-axis
        if self.mpl.rcParams==self.mpl.rcParamsDefault:
            # X-axis widened by 1 inch.
            self.plt.gcf().set_size_inches([self.plt.gcf().get_size_inches()[0]+5,self.plt.gcf().get_size_inches()[1]+1])
        self.plt.savefig(os.path.join(self.workdir,'bperp_plot.eps'))
        self.plt.close()

        # Make Baseline histogram
        ax1=self.plt.figure().add_subplot(111)
        ax1.hist(baseline_hist)
        ax1.set_xlabel(r'$\perp$'+' Baseline (m)', weight='bold')
        ax1.set_ylabel('Number of Interferograms', weight='bold')
        ax1.yaxis.set_major_locator(self.MaxNLocator(integer=True)) #Force y-axis to only use ints
        self.plt.tight_layout()
        self.plt.savefig(os.path.join(self.workdir,'bperp_histogram.eps'))
        self.plt.close()

        return

    def plot_extents(self):
        '''
            Make plot of track extents vs bounding bbox/common track extent.
        '''

        ax=self.plt.figure().add_subplot(111)
        #Iterate through all IFGs
        S_extent=[]
        N_extent=[]
        for i, j in enumerate(self.product_dict[0]):
            prods_bbox=open_shapefile(j[0], 0, 0).bounds
            S_extent.append(prods_bbox[1])
            N_extent.append(prods_bbox[3])
            # Plot IFG extent bounds in latitude
            ax.plot([self.product_dict[1][i][0]]*2,list(prods_bbox[1::2]),'ko',markersize=10)
            # Plot IFG extent line connecting bounds in latitude
            ax.plot([self.product_dict[1][i][0]]*2,list(prods_bbox[1::2]), color='0.5', linestyle='--')

        # Plot bounds of common track extent
        if self.croptounion:
            S_extent=min(S_extent)
            N_extent=max(N_extent)
            if [self.bbox_file[1], self.bbox_file[3]]!=[S_extent, N_extent] and S_extent!=N_extent:
                ax.axhline(y=S_extent, color='r', linestyle=':', label="extent of union")
                ax.axhline(y=N_extent, color='r', linestyle=':')
        else:
            S_extent=max(S_extent)
            N_extent=min(N_extent)
            if [self.bbox_file[1], self.bbox_file[3]]!=[S_extent, N_extent] and S_extent!=N_extent:
                    ax.axhline(y=S_extent, color='r', linestyle=':', label="extent of intersection")
                    ax.axhline(y=N_extent, color='r', linestyle=':')

        # Plot bounds of final track extent all IFGs will be cropped to
        ax.axhline(y=self.bbox_file[1], color='b', linestyle='--', label="bounding box")
        ax.axhline(y=self.bbox_file[3], color='b', linestyle='--')

        # add legend
        self.plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

        # defining the axis labels
        ax.set_ylabel('Latitude',weight='bold')
        ax.set_xlabel('Interferograms', weight='bold')
        ax.set_title('Interferogram lat extents', weight='bold')
        self.plt.xticks(rotation=90)
        self.plt.tight_layout()

        # saving the figure
        self.plt.savefig(os.path.join(self.workdir,'lat_extents.eps'))
        self.plt.close()

        return

    def plot_coherence(self):
        '''
            Make coherence plot + histogram.
        '''

        fig, ax = self.plt.subplots()
        # ax=self.plt.figure().add_subplot(111)
        coh_hist = []

        # Iterate through all IFGs
        masters = []; slaves = []
        for i,j in enumerate(self.product_dict[0]):
            coh_file=gdal.Warp('', j, options=gdal.WarpOptions(format="MEM", cutlineDSName=self.prods_TOTbbox, outputBounds=self.bbox_file))
            coh_file_arr=np.ma.masked_where(coh_file.ReadAsArray() == coh_file.GetRasterBand(1).GetNoDataValue(), coh_file.ReadAsArray())

            # Apply mask (if specified).
            if self.mask is not None:
                coh_file_arr=np.ma.masked_where(self.mask == 0.0, coh_file_arr)

            # Report mean
            coh_val=coh_file_arr.mean()

            # Record average coherence val for histogram
            coh_hist.append(coh_val)
            slaves.append(self.pd.to_datetime(self.product_dict[1][i][0][:8]))
            masters.append(self.pd.to_datetime(self.product_dict[1][i][0][9:]))

        # Plot average coherence per IFG
        cols, mapper = self._create_colors_coh(coh_hist)
        ax.set_prop_cycle(color=cols)
        lines       = ax.plot([masters, slaves],
                          [coh_hist]*2)
        scatter     = ax.scatter(slaves, coh_hist, c='k', zorder=100)
        master      = ax.scatter(masters, coh_hist, c='k', zorder=100)

        cbar_ax     = fig.add_axes([0.91, 0.12, 0.02, 0.75])
        self.warnings.filterwarnings("ignore",category=UserWarning)
        cbar        = fig.colorbar(mapper, cbar_ax, spacing='proportional')
        # cbar.set_label(lbl, rotation=90, labelpad=15)

        ### Make average coherence plot
        ax.set_ylabel('Average Coherence',weight='bold')
        ax.set_xlabel('Time',weight='bold')
        ax.set_ylim(0, 1)
        self.pairs=[i[0] for i in self.product_dict[1]]; xticks=self.__date_list__()
        xticks, labels = self._adaptive_xticks(list(set(xticks.keys())))
        ax.set_xlim(min(xticks),max(xticks))
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
        ax.xaxis.set_major_formatter(self.mpl.dates.DateFormatter('%Y-%m'))
        for label in ax.get_xticklabels():
                 label.set_ha('center')
                 label.set_rotation(20.)

        # If plotting parameters are the default, must adjust X-axis
        if self.mpl.rcParams==self.mpl.rcParamsDefault:
            # X-axis widened by 1 inch.
            self.plt.gcf().set_size_inches([self.plt.gcf().get_size_inches()[0]+5,self.plt.gcf().get_size_inches()[1]+1])
        self.plt.savefig(os.path.join(self.workdir,'avgcoherence_plot.eps'))
        self.plt.close()

        ### Make average coherence histogram
        ax1=self.plt.figure().add_subplot(111)
        ax1.hist(coh_hist)
        ax1.set_xlabel('Average Coherence', weight='bold')
        ax1.set_ylabel('Number of Interferograms',weight='bold')
        ax1.yaxis.set_major_locator(self.MaxNLocator(integer=True)) #Force y-axis to only use ints
        self.plt.tight_layout()
        self.plt.savefig(os.path.join(self.workdir,'avgcoherence_histogram.eps'))
        self.plt.close()

        return

    def plot_avgcoherence(self):
        '''
            Generate average coherence raster.
        '''

        # Import functions
        from vrtmanager import renderVRT
        outname=os.path.join(self.workdir,'avgcoherence')

        # Iterate through all IFGs
        for i,j in enumerate(self.product_dict[0]):
            coh_file=gdal.Warp('', j, options=gdal.WarpOptions(format="MEM", cutlineDSName=self.prods_TOTbbox, outputBounds=self.bbox_file))
            coh_file_arr=np.ma.masked_where(coh_file.ReadAsArray() == coh_file.GetRasterBand(1).GetNoDataValue(), coh_file.ReadAsArray())

            # Apply mask (if specified).
            if self.mask is not None:
                coh_file_arr=np.ma.masked_where(self.mask == 0.0, coh_file_arr)

            # Iteratively update average coherence file
            # If looping through first coherence file, nothing to sum so just save to file
            if os.path.exists(outname):
                coh_file=gdal.Open(outname,gdal.GA_Update)
                coh_file=coh_file.GetRasterBand(1).WriteArray(coh_file_arr+coh_file.ReadAsArray())
            else:
                renderVRT(outname, coh_file_arr, geotrans=coh_file.GetGeoTransform(), drivername=self.outputFormat, gdal_fmt=coh_file_arr.dtype.name, proj=coh_file.GetProjection(), nodata=coh_file.GetRasterBand(1).GetNoDataValue())
            coh_file = coh_val = coh_file_arr = None

        # Take average of coherence sum
        coh_file=gdal.Open(outname,gdal.GA_Update)
        coh_file=coh_file.GetRasterBand(1).WriteArray(coh_file.ReadAsArray()/len(self.product_dict[0]))
        coh_file = None

        return

    def plotbperpcoh(self):
        '''
            Make pbaseline plot that is color-coded w.r.t. coherence.
        '''

        # importing dependencies
        from matplotlib import cm
        fig, ax = self.plt.subplots()
        self.pairs=[i[0] for i in self.product_dict[1]]
        dateDict = self.__date_list__()
        A,B,L,baseline_hist = self.__design_matrix__()

        B1 = np.linalg.pinv(B)
        B1 = np.array(B1,np.float32)
        dS = np.dot(B1,L)
        dtbase = np.diff(list(dateDict.values()))
        dt = np.zeros((len(dtbase),1))
        zero = np.array([0.],np.float32)
        S = np.concatenate((zero,np.cumsum([dS*dtbase])))
        residual = L-np.dot(B,dS)

        RMSE = np.sqrt(np.sum(residual**2)/len(residual))
        if np.linalg.matrix_rank(B)!=len(list(dateDict.keys()))-1:
            print('Baseline plot warning!')
            print('Design matrix is rank deficient. Network is disconnected.')
            print('Using a fully connected network is recommended.')

        offset_dict={}
        # Plot dot for each date
        for i in range(len(list(dateDict.keys()))):
            offset_dict[list(dateDict.keys())[i]]=S[i]
            master=self.pd.to_datetime(list(dateDict.keys())[i][:8])
            ax.plot(master, S[i], 'k.', markeredgewidth = 3, markersize=15, linestyle='None', zorder=10)
        slaves = []; masters = []; coh_vals = []; y1 = []; y2 = []
        for i,j in enumerate(self.pairs): #Plot lines for each pair
            slaves.append(self.pd.to_datetime(j[:8]))
            masters.append(self.pd.to_datetime(j[9:]))
            # Open coherence file
            coh_file=gdal.Warp('', self.product_dict[2][i], options=gdal.WarpOptions(format="MEM", cutlineDSName=self.prods_TOTbbox, outputBounds=self.bbox_file))
            coh_file_arr=np.ma.masked_where(coh_file.ReadAsArray() == coh_file.GetRasterBand(1).GetNoDataValue(), coh_file.ReadAsArray())

            # Apply mask (if specified).
            if self.mask is not None:
                coh_file_arr=np.ma.masked_where(self.mask == 0.0, coh_file_arr)

            # Report mean
            coh_val=coh_file_arr.mean()
            coh_vals.append(coh_val)
            y1.append(offset_dict[j[9:]])
            y2.append(offset_dict[j[:8]])
            coh_file = coh_file_arr = coh_val = None

        cols, mapper = self._create_colors_coh(coh_vals, 'autumn') # don't use hot
        ax.set_prop_cycle(color=cols)
        lines       = ax.plot([masters, slaves], [y1, y2])
        cbar_ax     = fig.add_axes([0.91, 0.12, 0.02, 0.75])
        self.warnings.filterwarnings("ignore",category=UserWarning)
        cbar        = fig.colorbar(mapper, cbar_ax)
        # cbar_ax.set_title('Coherence', loc='left')
        cbar_ax.set_ylabel('Coherence', rotation=-90, labelpad=17)
        ax.set_ylabel(r'$\perp$'+' Baseline (m)',weight='bold')
        ax.set_xlabel('Time',weight='bold')

        ### Make baseline plot
        xticks, labels = self._adaptive_xticks(list(set(dateDict.keys())))
        ax.set_xlim(min(xticks),max(xticks))
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels)
        ax.xaxis.set_major_formatter(self.mpl.dates.DateFormatter('%Y-%m'))
        for label in ax.get_xticklabels():
                 label.set_ha('center')
                 label.set_rotation(20.)

        # If plotting parameters are the default, must adjust X-axis
        if self.mpl.rcParams==self.mpl.rcParamsDefault:
            # X-axis widened by 1 inch.
            self.plt.gcf().set_size_inches([self.plt.gcf().get_size_inches()[0]+5,self.plt.gcf().get_size_inches()[1]+1])
        self.plt.savefig(os.path.join(self.workdir,'bperp_coh_plot.eps'))
        self.plt.close()

        return

    def _create_colors_coh(self, vals, cm='autumn'):
        """ create colors from a set of values between 0/1"""
        if not isinstance(vals, np.ndarray):
            vals = np.array(vals)
        norm        = self.mpl.colors.Normalize(vmin=0, vmax=1)# clip=False)
        cmap        = self.plt.cm.get_cmap(cm)
        mapper      = self.plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        mapper._A   = [] # necessary dummy array for cbar
        colors      =  [mapper.to_rgba(v) for v in vals]
        return colors, mapper

    def _adaptive_xticks(self, dates):
        """ Adjust the number of xticks based on the time interval """
        # dates = ['20150310', '20160410', '20160410', '20220511']
        dates = [self.datetime.strptime(i, '%Y%m%d') for i in dates]
        dates.sort()
        st   = self.datetime(min(dates).year, 1, 1)
        en   = max(dates) + self.relativedelta(months=1)
        elap = en - min(dates)

        if len(dates) == 2 or elap.days <= 365*2.5:
            st = min(dates).replace(day=1)
            labels = self.pd.date_range(st, en, freq='MS')
        elif elap.days > 365*2.5 and elap.days <= 365*5.5:
            labels = self.pd.date_range(st, en, freq='3MS')
        elif elap.days > 365*5.5 and elap.days <= 365*8.5:
            labels = self.pd.date_range(st, en, freq='6MS')
        else:
            labels = self.pd.date_range(st, en, freq='AS')

        xticks = [x.toordinal() for x in labels]
        return xticks, labels

if __name__ == '__main__':
    '''
        Main driver.
    '''
    inps = cmdLineParse()

    print("***Plotting Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    # If user requests to generate all plots.
    if inps.plotall:
        print('"-plotall"==True. All plots will be made.')
        inps.plottracks=True
        inps.plotbperp=True
        inps.plotcoh=True
        inps.plotbperpcoh=True
        inps.makeavgoh=True


    if inps.plottracks or inps.plotcoh or inps.makeavgoh or inps.plotbperpcoh:
        # Import functions
        from extractProduct import merged_productbbox
        # extract/merge productBoundingBox layers for each pair and update dict,
        # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
        standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)
        # Load mask (if specified).
        if inps.mask is not None:
            inps.mask=gdal.Warp('', inps.mask, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds, dstNodata=0))
            inps.mask.SetProjection(proj)
            # If no data value
            if inps.mask.GetRasterBand(1).GetNoDataValue():
                inps.mask=np.ma.masked_where(inps.mask.ReadAsArray() == inps.mask.GetRasterBand(1).GetNoDataValue(), inps.mask.ReadAsArray())
            else:
                inps.mask=inps.mask.ReadAsArray()


    # Make spatial extent plot
    if inps.plottracks:
        print("- Make plot of track latitude extents vs bounding bbox/common track extent.")
        make_plot=plot_class([[j['productBoundingBox'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'productBoundingBox'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, croptounion=inps.croptounion)
        make_plot.plot_extents()


    # Make pbaseline plot
    if inps.plotbperp:
        print("- Make baseline plot and histogram.")
        make_plot=plot_class([[j['bPerpendicular'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'bPerpendicular'))
        make_plot.plot_pbaselines()


    # Make average land coherence plot
    if inps.plotcoh:
        print("- Make average IFG coherence plot in time, and histogram of average IFG coherence.")
        make_plot=plot_class([[j['coherence'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'coherence'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask)
        make_plot.plot_coherence()


    # Generate average land coherence raster
    if inps.makeavgoh:
        print("- Generate 2D raster of average coherence.")
        make_plot=plot_class([[j['coherence'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'coherence'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask, outputFormat=inps.outputFormat)
        make_plot.plot_avgcoherence()


    # Make pbaseline/coherence combo plot
    if inps.plotbperpcoh:
        print("- Make baseline plot that is color-coded with respect to mean IFG coherence.")
        make_plot=plot_class([[j['bPerpendicular'] for j in standardproduct_info.products[1]],  [j["pair_name"] for j in standardproduct_info.products[1]], [j['coherence'] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'bPerpendicular'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask)
        make_plot.plotbperpcoh()
