#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For a given set of interferograms, compute the cumulative phase
#  misclosure.
# 
# Rob Zinke 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### IMPORT MODULES ---
import os
from glob import glob
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import mode
from datetime import datetime


### PARSER ---
Description='''Compute the cumulative misclosure of phase triplets based on a set of interferograms saved in the MintPy HDF5 data stucture.
This program assumes that all files are cover the same area in space, with the same spatial resolution. During triplet computation,
values at a reference point are removed from the interferograms prior to misclosure computation to account for abiguities in 
the unwrapped phase that might arise during pairwise computation.

The code works by reading in a list of triplets and their date pairs, and from those:

    1. Formulate a list of viable triplet combinations (IJ, JK, IK)
    2. Subtract the reference point from those interferograms
    3. Compute the disagreement based on IJ+JK-IK
       ... continue through the list of triplets.

Both the misclosure values (positive or negative) and the absolute misclosure values (always positive) are computed and stored in 3D arrays, where each slice represents a triplet. The "cumulative" misclosure and absolute misclosure are computed by summing the 3D arrays.

Once the (absolute) misclosure is calculated, the user can view the time history of any pixel by clicking on either map.
'''

Examples='''EXAMPLES
[from a MintPy "inputs" directory]
computeCumMisclosure.py ifgramStack.h5 --exclude-pairs 20190103_20190115 20160530_20170724 20160530_20170618 20160530_20170606 20160530_20160810 -refX 38 -refY 212

[from an ARIA-tools "unwrappedPhase" directory]
computeCumMisclosure.py '2017.vrt' -refX 40 -refY 220 --plot-pairs

[from an ISCE stack "merged/interferograms" directory]
computeCumMisclosure.py '*' --print-files
'''

def createParser():
    import argparse
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Input data
    parser.add_argument(dest='dataset', type=str, help='Name of folders, files, or HDF5 dataset')
    parser.add_argument('-t','--dataType', dest='dataType', type=str, default=None, help='(Recommended) Manually specify data type (ARIA, ISCE, MintPy, [None])')
    parser.add_argument('-s','--subDS', dest='subDS', type=str, default='unwrapPhase', help='Sub-data set [e.g., unwrapPhase_phaseClosure, default=unwrapPhase]')
    parser.add_argument('--exclude-pairs', dest='exclPairs', type=str, nargs='+', default=None, help='Pairs to exclude \"YYYYMMDD_YYYYMMDD\"')

    # Reference point
    parser.add_argument('-refX', dest='refX', default='auto', help='Reference X pixel')
    parser.add_argument('-refY', dest='refY', default='auto', help='Reference Y pixel')

    # Triplet formulation
    parser.add_argument('-l','--lags', dest='lags', type=int, default=1, help='Number of lags, e.g., 1 lags = [n1-n0, n2-n1, n2-n0]')
    parser.add_argument('--mintime', dest='minTime', type=str, default=None, help='Minimum time span of pairs in triplets (days)')
    parser.add_argument('--maxtime', dest='maxTime', type=str, default=None, help='Maximum time span of pairs in triplets (days)')

    # Vocalization
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true', help='Verbose mode')
    parser.add_argument('--print-files', dest='printFiles', action='store_true', help='Print list of detected files')
    parser.add_argument('--print-dates', dest='printDates', action='store_true', help='Print list of unique dates')
    parser.add_argument('--print-poss-triplets', dest='printPossTriplets', action='store_true', help='Print list of all possible triplets (not just those included in the data set)')
    parser.add_argument('--print-triplets', dest='printTriplets', action='store_true', help='Print list of triplets validated against existing date pairs')

    # Plots
    parser.add_argument('--pctmin', dest='pctMinClip', type=float, default=1, help='Minimum percent clip value for cumulative misclosure plot')
    parser.add_argument('--pctmax', dest='pctMaxClip', type=float, default=99, help='Maximum percent clip value for cumulative misclosure plot')

    parser.add_argument('--plot-inputs', dest='plotInputs', action='store_true', help='Plot input interferograms')

    # Misclosure map formatting
    parser.add_argument('--misclosure-limits', dest='miscLims', type=float, nargs=2, default=None, help='Cumulative misclosure plot color limits')
    parser.add_argument('--abs-misclosure-limit', dest='absMiscLim', type=float, default=None, help='Absolute cumulative misclosure plot color limits')

    return parser


def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### LOAD DATA ---
## Detect data type
def detectDataType(inps):
    # Search for specified data files
    inps.files=glob(inps.dataset)
    if inps.printFiles == True: print('Files detected:\n{}'.format(inps.files))

    # Summarize file detection
    inps.nFiles=len(inps.files)
    if inps.verbose == True: print('{} files detected'.format(inps.nFiles))

    # Determine data type (ARIA, ISCE, MintPy)
    if inps.dataType:
        # User-specified data type
        assert inps.dataType.lower() in ['aria','isce','mintpy'], \
            'Data type not found. Must be ARIA, ISCE, MintPy'
    else:
        # Auto-detect data type
        if inps.dataset[-3:]=='.h5': 
            inps.dataType='MintPy'
        elif inps.dataset[-4:]=='.vrt':
            inps.dataType='ARIA'
        else:
            # Check if files are directories
            try:
                os.listdir(inps.files[0])
                inps.dataType='ISCE'
            except:
                print('Data type unrecognized. Please specify explicitly'); exit()

    # Report if requested
    if inps.verbose == True: print('Data type: {}'.format(inps.dataType))


## Load data as 3D array
def loadData(inps):
    if inps.verbose == True: print('Loading data...')

    # Detect data type (ARIA, ISCE, MintPy)
    detectDataType(inps)

    # Load data based on filetype
    if inps.dataType.lower()=='aria':
        data=loadARIA(inps)
    elif inps.dataType.lower()=='isce':
        data=loadISCE(inps)
    elif inps.dataType.lower()=='mintpy':
        data=loadMintPy(inps)
    else:
        print('Could not identify data type - please specify explicityly.\nNo data loaded.'); exit()

    # Update variables to reflect map size, epochs, etc.
    data.updateVariables()

    return data


## Load ARIA data set
def loadARIA(inps):
    from osgeo import gdal
    data=dataSet(verbose=inps.verbose) # instatiate data set object

    # Exclude files
    if inps.exclPairs:
        inps.files=[file for file in inps.files if file.split('.')[0] not in inps.exclPairs]

    # Loop through to load each file and append to list
    for file in inps.files:
        # Load using gdal
        DS=gdal.Open(file,gdal.GA_ReadOnly)

        # Append IFG to list
        data.ifgs.append(DS.GetRasterBand(1).ReadAsArray())

        # Append pair name to list
        fname=os.path.basename(file)
        pairLabel=fname.split('.')[0] # remove extension
        pair=pairLabel.split('_')[::-1] # reverse for older-younger
        data.pairLabels.append(pairLabel)
        data.pairs.append(pairName)

    # Convert IFG list to array
    data.ifgs=np.array(data.ifgs)
    data.M,data.N=data.ifgs[0,:,:].shape # map dimensions
    if inps.verbose == True: print('Data array shape: {}'.format(data.ifgs.shape))

    return data


## Load ISCE data set
def loadISCE(inps):
    from osgeo import gdal
    data=dataSet(verbose=inps.verbose) # instatiate data set object

    # Exclude files
    if inps.exclPairs:
        inps.files=[file for file in inps.files if file not in inps.exclPairs]

    # Loop through to load each file and append to list
    for file in inps.files:
        # Load using gdal
        ifgName='{}/filt_fine.unw.vrt'.format(file) # format name
        DS=gdal.Open(ifgName,gdal.GA_ReadOnly) # open

        # Append IFG to list
        data.ifgs.append(DS.GetRasterBand(1).ReadAsArray())

        # Append pair name to list
        fname=os.path.basename(file)
        pairLabel=fname
        pairName=pairLabel.split('_')
        data.pairLabels.append(pairLabel)
        data.pairs.append(pairName)

    # Convert IFG list to array
    data.ifgs=np.array(data.ifgs)
    data.M,data.N=data.ifgs[0,:,:].shape # map dimensions
    if inps.verbose == True: print('Data array shape: {}'.format(data.ifgs.shape))

    return data


## Load MintPy data set
def loadMintPy(inps):
    import h5py
    data=dataSet(verbose=inps.verbose) # instatiate dataset object

    # Load HDF5 file
    with h5py.File(inps.files[0],'r') as DS:
        # Report keys and data sets if requested
        if inps.verbose == True: 
            print('Data sets: {}'.format(DS.keys()))
            print('Using subset: {}'.format(inps.subDS))

        # Load data cube
        data.ifgs=DS[inps.subDS][:]

        # Load dates
        pairs=DS['date'][:].astype(str)
        data.pairs=[[pair[0],pair[1]] for pair in pairs]
        data.pairLabels=['{}_{}'.format(pair[0],pair[1]) for pair in pairs]

        # Close data set
        DS.close()    

    # Exclude interferograms
    if inps.exclPairs:
        # Formatted list of pair names
        pairList=['{}_{}'.format(pair[0],pair[1]) for pair in data.pairs]

        # Rebuild data cube and date list based on non-excluded pairs
        validPairs=[]; validIFGs=[]
        for n,pair in enumerate(pairList):
            if pair not in inps.exclPairs:
                validPairs.append(data.pairs[n])
                validIFGs.append(data.ifgs[n,:,:])
        # Reassign valid pairs and ifgs to data object
        data.pairs=validPairs
        data.ifgs=np.array(validIFGs)

    # Format IFGs
    data.M,data.N=data.ifgs[0,:,:].shape # map dimensions
    if inps.verbose == True: print('Data array shape: {}'.format(data.ifgs.shape))

    return data



### DATA CLASS ---
## Class containing input data and dates, etc.
class dataSet:
    ## Initialize
    def __init__(self,verbose):
        # Initial values, to be changed later
        # Inputs
        self.verbose=verbose

        self.refY=0
        self.refX=0

        # Data
        self.ifgs=[]
        self.pairs=[]
        self.pairLabels=[]
        self.dates=[]


    ## Update variables from data set
    def updateVariables(self,printDates=False):
        # List of epochs (unique dates of acquisition)
        allDates=[] # all dates, including redundant
        [allDates.extend(pair) for pair in self.pairs]
        [self.dates.append(date) for date in allDates if date not in self.dates]
        self.dates.sort() # chronological order
        self.nDates=len(self.dates)

        if printDates == True: print('Unique dates:\n{}'.format(self.dates))
        if self.verbose == True: print('Nb unique dates: {}'.format(self.nDates))

        # Data dimensions
        self.K,self.M,self.N=self.ifgs.shape

        if self.verbose == True: 
            print('Nb ifgs: {}'.format(self.K))
            print('Map dimensions: {} x {}'.format(self.M,self.N))


    ## Create list of triplets
    # Create list of all possible triplets given date list
    def createTriplets(self,lags=1,minTime=None,maxTime=None,printTriplets=False):
        '''
            Provide a list of unique dates in format YYYYMMDD. This
             function will create a list of the (n1-n0, n2-n1, n2-n0)
             date combinations. 
            Lags is the minimum interval from one acquisition to the 
             next. For instance:
                lags=1 gives [n1-n0, n2-n1, n0-n2]
                lags=2 gives [n2-n0, n4-n2, n0-n4]
                lags=3 gives [n3-n0, n6-n3, n0-n6]
        '''
        if self.verbose == True: print('Creating list of all possible triplets')

        # Loop through dates to create possible triplet combinations
        triplets=[]
        for n in range(self.nDates-2*lags):
            dateI=self.dates[n] # first date in sequence
            dateJ=self.dates[n+lags] # second date in sequence
            dateK=self.dates[n+2*lags] # third date in sequence
            pairList=[[dateI,dateJ],[dateJ,dateK],[dateI,dateK]]
            triplets.append(pairList) # add to list

        # Check that pairs meet time requirements
        if minTime:
            # Convert pairs to intervals in days
            intervals=[]
            for triplet in triplets:
                intervalSet=[(datetime.strptime(pair[1],'%Y%m%d')-datetime.strptime(pair[0],'%Y%m%d')).days for pair in triplet]
                intervals.append(min(intervalSet))
            validTriplets=[triplet for ndx,triplet in enumerate(triplets) if intervals[ndx]>=int(minTime)]
            triplets=validTriplets

        if maxTime:
            # Convert pairs to intervals in days
            intervals=[]
            for triplet in triplets:
                intervalSet=[(datetime.strptime(pair[1],'%Y%m%d')-datetime.strptime(pair[0],'%Y%m%d')).days for pair in triplet]
                intervals.append(max(intervalSet))
            print(intervals)
            validTriplets=[triplet for ndx,triplet in enumerate(triplets) if intervals[ndx]<=int(maxTime)]
            triplets=validTriplets

        # Write to data object
        self.triplets=triplets
        self.nTriplets=len(self.triplets)

        # Report if requested
        if printTriplets == True:
            print('All possible triplets:')
            [print(triplet) for triplet in self.triplets]
        if self.verbose == True:
            print('{} possible triplets'.format(self.nTriplets))

        # Trim list to only "valid" triplets
        self.validateTriplets()

    # Limit list of triplets to only combinations for which data exist
    def validateTriplets(self):
        validTriplets=[]
        for triplet in self.triplets:
            count=0
            for pair in triplet:
                if pair in self.pairs: count+=1
            if count==3: validTriplets.append(triplet)
        self.triplets=validTriplets
        self.nTriplets=len(self.triplets)

        if self.verbose == True: print('{} valid triplets'.format(self.nTriplets))


    ## Compute misclosure
    def computeMisclosure(self):
        if self.verbose == True: print('Calculating misclosure...')

        # Empty placeholders
        self.miscStack=[]
        self.absMiscStack=[]

        # Compute misclosure map for each triplet
        for triplet in self.triplets:
            # Triplet date pairs
            IJdates=triplet[0]
            JKdates=triplet[1]
            IKdates=triplet[2]

            # Triplet ifg indices
            IJndx=self.pairs.index(IJdates)
            JKndx=self.pairs.index(JKdates)
            IKndx=self.pairs.index(IKdates)

            # Interferograms
            IJ=self.ifgs[IJndx,:,:]
            JK=self.ifgs[JKndx,:,:]
            IK=self.ifgs[IKndx,:,:]

            # Normalize to reference point
            IJ-=IJ[self.refY,self.refX]
            JK-=JK[self.refY,self.refX]
            IK-=IK[self.refY,self.refX]

            # Compute (abs)misclosure
            misclosure=IJ+JK-IK
            absMisclosure=np.abs(misclosure)

            # Append to stack
            self.miscStack.append(misclosure)
            self.absMiscStack.append(absMisclosure)

        # Convert lists to 3D arrays
        self.miscStack=np.array(self.miscStack)
        self.absMiscStack=np.array(self.absMiscStack)

        # Cumulative misclosure
        self.cumMisclosure=np.sum(self.miscStack,axis=0)
        self.cumAbsMisclosure=np.sum(self.absMiscStack,axis=0)

        # Use first datum from each triplet as reference date
        self.refDates=[str(triplet[0][0]) for triplet in self.triplets]
        self.refDatetimes=[datetime.strptime(str(triplet[0][0]),'%Y%m%d') for triplet in self.triplets]


    ## Plot misclosure
    # Detect background values
    def maskBackground(self,img):
        edges=np.concatenate([img[:,0].flatten(),
            img[0,:].flatten(),
            img[-1,:].flatten(),
            img[:,-1].flatten()])
        background=mode(edges)[0][0]

        img=np.ma.array(img,mask=(img==background))

        return img

    # Image percentiles
    def imgMinMax(self,img,percentiles):
        clipVals={}
        clipVals['min'],clipVals['max']=np.percentile(img.flatten(),percentiles)

        print('Clipping values to:\nMin {} % = {}\nMax {} % = {}'.format(\
                percentiles[0],clipVals['min'],percentiles[1],clipVals['max']))        

        return clipVals

    # Plotting function for cumulative misclosure
    def miscMapPlot(self,Fig,ax,img,title,cmap,minVal,maxVal):
        cax=ax.imshow(img,cmap=cmap,vmin=minVal,vmax=maxVal)
        ax.plot(self.refX,self.refY,'ks')
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(title)

        return cax

    # Plot timeseries
    def plotSeries(self,Fig,ax,series,name,timeAxis=False):
        # Plot data
        ax.plot(self.refDatetimes,series,'-k.')

        # Format axis
        ax.set_ylabel(name+'\nradians')
        if timeAxis == False:
            ax.set_xticklabels([])
        else:
            ax.set_xticks(self.refDatetimes)
            ax.set_xticklabels(self.refDates,rotation=80)

    # Plot cumulative misclosure
    def plotMisclosure(self,inps):
        ## Plot cumulative misclosure
        # Prep image by masking background and determining percentiles
        if self.verbose == True: print('Plotting cumulative misclosure')
        self.cumMiscFig=plt.figure(); self.cumMiscAx=self.cumMiscFig.add_subplot(111)
        self.cumMisclosure=self.maskBackground(self.cumMisclosure)
        self.cumMiscStats=self.imgMinMax(self.cumMisclosure,[inps.pctMinClip,inps.pctMaxClip])

        # Plot cumulative misclosure
        cax=self.miscMapPlot(self.cumMiscFig,self.cumMiscAx,self.cumMisclosure,title=None,
            cmap='plasma',minVal=self.cumMiscStats['min'],maxVal=self.cumMiscStats['max'])
        cbar=self.cumMiscFig.colorbar(cax,orientation='horizontal')
        cbar.set_label('radians')


        ## Plot cumulative absolute misclosure
        # Prep image by masking background and determining percentiles
        if self.verbose == True: print('Plotting cumulative absolute misclosure')
        self.cumAbsMiscFig=plt.figure(); self.cumAbsMiscAx=self.cumAbsMiscFig.add_subplot(111)
        self.cumAbsMisclosure=self.maskBackground(self.cumAbsMisclosure)
        self.cumAbsMiscStats=self.imgMinMax(self.cumAbsMisclosure,[0,inps.pctMaxClip])

        # Plot cumulative absolute misclosure
        cax=self.miscMapPlot(self.cumAbsMiscFig,self.cumAbsMiscAx,self.cumMisclosure,title=None,
            cmap='plasma',minVal=self.cumAbsMiscStats['min'],maxVal=self.cumAbsMiscStats['max'])
        cbar=self.cumAbsMiscFig.colorbar(cax,orientation='horizontal')
        cbar.set_label('radians')


        ## Plot timeseries points
        self.miscSeriesFig=plt.figure('Misclosure',figsize=(8,8))
        self.miscSeriesAx=self.miscSeriesFig.add_subplot(411)
        self.cumMiscSeriesAx=self.miscSeriesFig.add_subplot(412)
        self.absMiscSeriesAx=self.miscSeriesFig.add_subplot(413)
        self.cumAbsMiscSeriesAx=self.miscSeriesFig.add_subplot(414)


        # Link canvas to plots for interaction
        self.cumMiscFig.canvas.mpl_connect('button_press_event', self.misclosureAnalysis)
        self.cumAbsMiscFig.canvas.mpl_connect('button_press_event', self.misclosureAnalysis)



    ## Misclosure analysis
    def misclosureAnalysis(self,event):
        px=event.xdata; py=event.ydata
        px=int(round(px)); py=int(round(py))

        # Report position and cumulative values
        print('px {} py {}'.format(px,py)) # report position
        print('Cumulative misclosure: {}'.format(self.cumMisclosure[py,px]))
        print('Abs cumulative misclosure: {}'.format(self.cumAbsMisclosure[py,px]))

        # Plot query points on maps
        self.cumMiscAx.cla()
        self.miscMapPlot(self.cumMiscFig,self.cumMiscAx,self.cumMisclosure,title=None,
            cmap='plasma',minVal=self.cumMiscStats['min'],maxVal=self.cumMiscStats['max'])
        self.cumMiscAx.plot(px,py,color='k',marker='o',markerfacecolor='w',zorder=3)

        self.cumAbsMiscAx.cla()
        self.miscMapPlot(self.cumAbsMiscFig,self.cumAbsMiscAx,self.cumAbsMisclosure,title=None,
            cmap='plasma',minVal=self.cumAbsMiscStats['min'],maxVal=self.cumAbsMiscStats['max'])
        self.cumAbsMiscAx.plot(px,py,color='k',marker='o',markerfacecolor='w',zorder=3)

        # Plot misclosure over time
        print('Misclosure: {}'.format(self.miscStack[:,py,px]))
        self.miscSeriesAx.cla() # misclosure
        self.plotSeries(self.miscSeriesFig,self.miscSeriesAx,self.miscStack[:,py,px],'misclosure')
        self.cumMiscSeriesAx.cla() # cumulative misclosure
        self.plotSeries(self.miscSeriesFig,self.cumMiscSeriesAx,np.cumsum(self.miscStack[:,py,px]),'cum. miscl.')
        self.absMiscSeriesAx.cla() # absolute misclosure
        self.plotSeries(self.miscSeriesFig,self.absMiscSeriesAx,self.absMiscStack[:,py,px],'abs. miscl')
        self.cumAbsMiscSeriesAx.cla() # cumulative absolute misclosure
        self.plotSeries(self.miscSeriesFig,self.cumAbsMiscSeriesAx,np.cumsum(self.absMiscStack[:,py,px]),'cum. abs. miscl.',timeAxis=True)

        # Draw outcomes
        self.cumMiscFig.canvas.draw()
        self.cumAbsMiscFig.canvas.draw()
        self.miscSeriesFig.canvas.draw()



### ANCILLARY FUNCTIONS ---
## Imagette plotting
def plotImagettes(imgs,mRows,nCols,cmap='viridis',downsampleFactor=0,
    pctmin=0,pctmax=100,
    titleList=None,supTitle=None):

    # Number of imagettes
    nImgs=imgs.shape[0]

    # Number of imagettes per figure
    nbImagettes=mRows*nCols

    # Loop through image list
    x=1 # position variable
    for i in range(nImgs):
        # Generate new figure if needed
        if x%nbImagettes==1:
            Fig=plt.figure() # new figure
            x=1 # reset counter

        # Format image
        img=imgs[i,:,:] # current image from list
        ds=int(2**downsampleFactor) # downsample factor
        img=img[::ds,::ds] # downsample image
        img=np.ma.array(img,mask=(img==0)) # mask background

        imgMin,imgMax=np.percentile(img.compressed(),(pctmin,pctmax))

        # Plot image as subplot
        ax=Fig.add_subplot(mRows,nCols,x)
        ax.imshow(img,cmap=cmap,vmin=imgMin,vmax=imgMax)

        # Plot formatting
        ax.set_title(titleList[i])
        ax.set_xticks([]); ax.set_yticks([])
        Fig.suptitle(supTitle)

        x+=1 # update counter



### MAIN CALL ---
def main(inps=None):
    ## Gather arguments
    inps=cmdLineParse()

    ## Load data based on data type
    data=loadData(inps)

    # Plot data if requested
    if inps.plotInputs == True:
        plotImagettes(data.ifgs,4,5,cmap='jet',downsampleFactor=3,
            pctmin=1,pctmax=99,
            titleList=data.pairLabels,supTitle='Input data')


    ## Reference point
    if inps.refY == 'auto':
        data.refY=np.random.randint(0,data.M,1)
    else:
        data.refY=int(inps.refY)
    if inps.refX == 'auto':
        data.refX=np.random.randint(0,data.N,1)
    else:
        data.refX=int(inps.refX)


    ## Formulate valid triplets
    data.createTriplets(minTime=inps.minTime,maxTime=inps.maxTime,printTriplets=inps.printTriplets)


    ## Calculate misclosure
    data.computeMisclosure()


    ## Plot misclosure
    data.plotMisclosure(inps)


    plt.show()
