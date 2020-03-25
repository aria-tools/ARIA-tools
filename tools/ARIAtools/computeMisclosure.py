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
from datetime import datetime


### PARSER ---
Examples="""
"""

def createParser():
	import argparse
	parser = argparse.ArgumentParser(description='Compute the cumulative misclosure of phase triplets based on a \
		set of interferograms saved in the MintPy HDF5 data stucture.\nThis program assumes that all files are \
		coregistered and cover the same area.',
		epilog=Examples)

	# Input data
	parser.add_argument(dest='dataset', type=str, help='Name of folders, files, or HDF5 dataset')
	parser.add_argument('-t','--dataType', dest='dataType', type=str, default=None, help='(Recommended) Manually specify data type (ARIA, ISCE, MintPy, [None])')
	parser.add_argument('-s','--subDS', dest='subDS', type=str, default='unwrapPhase', help='Sub-data set [e.g., unwrapPhase_phaseClosure, default=unwrapPhase]')
	parser.add_argument('--exclude-pairs', dest='exclPairs', type=str, nargs='+', default=None, help='Pairs to exclude \"YYYYMMDD_YYYYMMDD\"')

	# Triplet formulation
	parser.add_argument('-l','--lags', dest='lags', type=int, default=1, help='Number of lags, e.g., 1 lags = [n1-n0, n2-n1, n2-n0]')
	parser.add_argument('--mintime','--min-time', dest='minTime', type=str, default=None, help='Minimum time span of pairs in triplets (days)')
	parser.add_argument('--maxtime','--max-time', dest='maxTime', type=str, default=None, help='Maximum time span of pairs in triplets (days)')

	# Reference point
	parser.add_argument('-refX', dest='refX', default='auto', help='Reference X pixel')
	parser.add_argument('-refY', dest='refY', default='auto', help='Reference Y pixel')

	# Vocalization
	parser.add_argument('-v','--verbose', dest='verbose', action='store_true', help='Verbose mode')
	parser.add_argument('--print-files', dest='printFiles', action='store_true', help='Print list of detected files')
	parser.add_argument('--print-dates', dest='printDates', action='store_true', help='Print list of unique dates')

	# Plots
	parser.add_argument('--plot-pairs', dest='plotPairs', action='store_true', help='Plot interferogram pairs in schematic form')
	parser.add_argument('--plot-inputs', dest='plotInputs', action='store_true', help='Plot input interferograms')
	parser.add_argument('--plot-triplets', dest='plotTriplets', action='store_true', help='Plot triplets and misclosure.')

	# Misclosure map formatting
	parser.add_argument('--misclosure-limits', dest='miscLims', type=float, nargs=2, default=None, help='Cumulative misclosure plot color limits')
	parser.add_argument('--abs-misclosure-limit', dest='absMiscLim', type=float, default=None, help='Absolute cumulative misclosure plot color limits')

	return parser


def cmdLineParse(iargs = None):
	parser = createParser()
	return parser.parse_args(args=iargs)



### ANCILLARY FUNCTIONS ---
## Find unique dates given a list of date pairs
def udatesFromPairs(datePairs,verbose=False):
	allDates=[] # empty list of all dates
	[allDates.extend(pair) for pair in datePairs]
	nAllDates=len(allDates)

	uniqueDates=[] # empty list for unique dates
	[uniqueDates.append(d) for d in allDates if d not in uniqueDates]
	nUniqueDates=len(uniqueDates)

	if verbose is True:
		print('Individual dates: {}'.format(nAllDates))
		print('Unique dates: {}'.format(nUniqueDates))

	return uniqueDates


## Image background
def imgBackground(I):
	# Use mode of background values
	from scipy.stats import mode
	edgeValues=np.concatenate([I[0,:],I[-1,:],I[:,0],I[:,-1]])
	background=mode(edgeValues).mode[0] # most common value
	return background


## Map stats
class mapStats:
	def __init__(self,I,pctmin=0,pctmax=100,verbose=False,hist=False): 
		# Guess at background value
		try:
			self.background=imgBackground(I)
		except:
			self.background=None
		# Check if masked array 
		try: 
			I=I.compressed() 
		except: 
			pass 
		# Convert to 1D array
		I=np.reshape(I,(1,-1)).squeeze(0) # 1D array 
		# Stats 
		self.min=np.min(I)	     # min 
		self.max=np.max(I)	     # max 
		self.mean=np.mean(I)	 # mean 
		self.median=np.median(I) # median 
		self.std=np.std(I)       # standard deviation 
		self.vmin,self.vmax=np.percentile(I,(pctmin,pctmax)) 
		# Print stats 
		if verbose is True: 
			print('Image stats:')
			print('\tmin: {}, max: {}'.format(self.min,self.max)) 
			print('\tmean: {}'.format(self.mean)) 
			print('\tmedian: {}'.format(self.median)) 
			print('\tvmin: {}, vmax: {}'.format(self.vmin,self.vmax)) 
			print('\tlikely background: {}'.format(self.background))
		# Histogram 
		if hist is not False: 
			if type(hist)==int: 
				nbins=hist 
			else: 
				nbins=50 
			# All values 
			H0,H0edges=np.histogram(I,bins=nbins) 
			H0cntrs=H0edges[:-1]+np.diff(H0edges)/2 
			# Clipped values 
			I=I[(I>=self.vmin) & (I<=self.vmax)] 
			H,Hedges=np.histogram(I,bins=nbins) 
			Hcntrs=Hedges[:-1]+np.diff(Hedges)/2 
			# Plot 
			plt.figure() 
			# Plot CDF 
			plt.subplot(2,1,1) 
			plt.axhline(pctmin/100,color=(0.5,0.5,0.5))
			plt.axhline(pctmax/100,color=(0.5,0.5,0.5)) 
			plt.plot(H0cntrs,np.cumsum(H0)/np.sum(H0),'k') 
			# Pad 
			H0cntrs=np.pad(H0cntrs,(1,1),'edge')
			H0=np.pad(H0,(1,1),'constant') 
			Hcntrs=np.pad(Hcntrs,(1,1),'edge')
			H=np.pad(H,(1,1),'constant') 
			# Plot PDF 
			plt.subplot(2,1,2)  
			plt.fill(H0cntrs,H0,color=(0.4,0.5,0.5),alpha=1,label='orig') 
			plt.bar(Hcntrs,H,color='r',alpha=0.5,label='new') 
			plt.legend()


## Imagette plotting
def imagettes(imgs,mRows,nCols,cmap='viridis',downsampleFactor=0,vmin=None,vmax=None,pctmin=None,pctmax=None,colorbarOrientation=None,background=None,
	extent=None,showExtent=False,titleList=None,supTitle=None):
	# If images are in a list, convert to 3D data cube
	if len(imgs)>1:
		imgs=np.array(imgs)

	# Number of imagettes
	nImgs=imgs.shape[0]

	# Number of imagettes per figure
	nbImagettes=mRows*nCols

	# Loop through image list
	x=1 # position variable
	for i in range(nImgs):
		# Generate new figure if needed
		if x%nbImagettes==1:
			F=plt.figure() # new figure
			x=1 # reset counter

		# Format image
		img=imgs[i,:,:] # current image from list

		ds=int(2**downsampleFactor) # downsample factor
		img=img[::ds,::ds] # downsample image

		if background is not None:
			if background=='auto':
				backgroundValue=imgBackground(img)
			else:
				backgroundValue=background
			img=np.ma.array(img,mask=(img==backgroundValue))

		if pctmin is not None:
			assert vmin is None and vmax is None, 'Specify either vmin/max or pctmin/max, not both'
			stats=mapStats(img,pctmin=pctmin,pctmax=pctmax); vminValue=stats.vmin; vmaxValue=stats.vmax 
		else:
			vminValue=vmin; vmaxValue=vmax

		# Plot image as subplot
		ax=F.add_subplot(mRows,nCols,x)
		cax=ax.imshow(img,extent=extent,
			cmap=cmap,vmin=vminValue,vmax=vmaxValue)

		# Plot formatting
		if titleList:
			if type(titleList[i]) in [float,np.float64]:
				titleList[i]=round(titleList[i],2)
			ax.set_title(titleList[i])
		else:
			ax.set_title(i)

		if supTitle:
			F.suptitle(supTitle)

		if showExtent is False:
			ax.set_xticks([]); ax.set_yticks([])

		if colorbarOrientation:
			F.colorbar(cax,orientation=colorbarOrientation)

		x+=1 # update counter


def plotDatePairs(pairs):
	# Provide pairs as a nested list of lists, e.g.,
	#  [[20190101,20181201]
	#   [20181201,20181101]]

	## Order dates earliest-latest based on ref (first) date
	pairs.sort(key=lambda d: d[0])

	## Convert to datetime objects
	int2date=lambda d: datetime.strptime(str(d),'%Y%m%d')

	timepairs=[[int2date(pair[0]),int2date(pair[1])] for pair in pairs]


	## Plot
	Fig=plt.figure()
	ax=Fig.add_subplot(111)
	n=0
	for timepair in timepairs:
		ax.plot(timepair,[n,n],'k-')
		n+=1

	## Formatting
	ax.set_yticks([])
	plt.xticks(rotation=90)
	ax.set_xlabel('Date')
	ax.set_title('Date pairs (n = {})'.format(len(pairs)))
	Fig.tight_layout()



### DATA CLASS ---
## Class containing input data and dates, etc.
class dataSet:
	def __init__(self):
		self.ifgs=[]
		self.pairs=[]
		self.dates=[]



### LOAD DATA ---
## Detect data type
def detectDataType(inps):
	# Search for specified data files
	inps.files=glob(inps.dataset)
	if inps.printFiles is True: print('Files detected:\n{}'.format(inps.files))

	# Summarize file detection
	inps.nFiles=len(inps.files)
	if inps.verbose is True: print('{} files detected'.format(inps.nFiles))

	# Determine data type (ARIA, ISCE, MintPy)
	if inps.dataType:
		assert inps.dataType.lower() not in ['aria','isce','mintpy'], \
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
				print('Data type unrecognized. Please specify explicitly')

	# Report if requested
	if inps.verbose is True: print('Data type: {}'.format(inps.dataType))


## Load data as 3D array
def loadData(inps):
	if inps.verbose is True: print('Loading data...')

	# Detect data type (ARIA, ISCE, MintPy)
	detectDataType(inps)

	# Load data based on filetype
	if inps.dataType=='ARIA':
		data=loadARIA(inps)
	elif inps.dataType=='ISCE':
		data=loadISCE(inps)
	elif inps.dataType=='MintPy':
		data=loadMintPy(inps)
	else:
		print('Could not identify data type - please specify explicityly.\nNo data loaded.'); exit()

	return data


## Load ARIA data set
def loadARIA(inps):
	from osgeo import gdal
	data=dataSet() # instatiate data set object

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
		fname=fname.split('.')[0] # remove extension
		pairName=fname.split('_')[::-1] # reverse for older-younger
		data.pairs.append(pairName)

	# Convert IFG list to array
	data.ifgs=np.array(data.ifgs)
	data.M,data.N=data.ifgs[0,:,:].shape # map dimensions
	if inps.verbose is True: print('Data array shape: {}'.format(data.ifgs.shape))

	# Unique dates from pairs
	data.dates=udatesFromPairs(data.pairs,verbose=inps.verbose); data.dates.sort()
	if inps.printDates is True: print('Dates: {}'.format(data.dates))

	return data


## Load ISCE data set
def loadISCE(inps):
	from osgeo import gdal
	data=dataSet() # instatiate data set object

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
		pairName=fname.split('_')
		data.pairs.append(pairName)

	# Convert IFG list to array
	data.ifgs=np.array(data.ifgs)
	data.M,data.N=data.ifgs[0,:,:].shape # map dimensions
	if inps.verbose is True: print('Data array shape: {}'.format(data.ifgs.shape))

	# Unique dates from pairs
	data.dates=udatesFromPairs(data.pairs,verbose=inps.verbose); data.dates.sort()
	if inps.printDates is True: print('Dates: {}'.format(data.dates))

	return data


## Load MintPy data set
def loadMintPy(inps):
	import h5py
	data=dataSet() # instatiate dataset object

	# Load HDF5 file
	with h5py.File(inps.files[0],'r') as DS:
		# Report keys and data sets if requested
		if inps.verbose is True: 
			print('Data sets: {}'.format(DS.keys()))
			print('Using subset: {}'.format(inps.subDS))

		# Load data cube
		data.ifgs=DS[inps.subDS][:]

		# Load dates
		pairs=DS['date'][:].astype(str)
		data.pairs=[[pair[0],pair[1]] for pair in pairs]

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
	if inps.verbose is True: print('Data array shape: {}'.format(data.ifgs.shape))

	# Unique dates from pairs
	data.dates=udatesFromPairs(data.pairs,verbose=inps.verbose); data.dates.sort()
	if inps.printDates is True: print('Dates: {}'.format(data.dates))

	return data



### Triplets
def createTriplets(dates,lags=1,minTime=None,maxTime=None,verbose=False):
	"""
		Provide a list of unique dates in format YYYYMMDD. This
		 function will create a list of the (n1-n0, n2-n1, n2-n0)
		 phase triplets. 
		It does not accept a list of pairs like the "formatTriplets"
		 function below.

		Lags is the minimum interval from one acquisition to the 
		 next. For instance:
			lags=1 gives [n1-n0, n2-n1, n0-n2]
			lags=2 gives [n2-n0, n4-n2, n0-n4]
			lags=3 gives [n3-n0, n6-n3, n0-n6]
	"""

	# Loop through dates to create valid triplet combinations
	nDates=len(dates)
	triplets=[]
	for n in range(nDates-2*lags):
		dateI=dates[n] # first date in sequence
		dateJ=dates[n+lags] # second date in sequence
		dateK=dates[n+2*lags] # third date in sequence
		pairList=[[dateI,dateJ],[dateJ,dateK],[dateI,dateK]]
		triplets.append(pairList) # add to list

	# Check that pairs meet time requirements
	if minTime:
		# Convert pairs to intervals in days
		intervals=[]
		for triplet in triplets:
			intervalSet=[daysBetween(pair[0],pair[1]) for pair in triplet]
			intervals.append(min(intervalSet))
		validTriplets=[triplet for ndx,triplet in enumerate(triplets) if intervals[ndx]>=int(minTime)]
		triplets=validTriplets

	if maxTime:
		# Convert pairs to intervals in days
		intervals=[]
		for triplet in triplets:
			intervalSet=[daysBetween(pair[0],pair[1]) for pair in triplet]
			intervals.append(max(intervalSet))
		print(intervals)
		validTriplets=[triplet for ndx,triplet in enumerate(triplets) if intervals[ndx]<=int(maxTime)]
		triplets=validTriplets

	# Print if requested
	if verbose is True:
		print('Triplets...')
		print('{} unique dates for triplet formulation'.format(nDates))
		print('Triplets:'); [print(triplet) for triplet in triplets]
		print('{} triplets created'.format(len(triplets)))

	return triplets



### CALCULATE MISCLOSURE ---
def calcMisclosure(inps,data):
	if inps.verbose is True: print('Calculating misclosure...')

	# Empty placeholders
	data.miscStack=[]
	data.absMiscStack=[]

	for triplet in inps.triplets:
		if inps.verbose is True: print(triplet)
		# Triplet date pairs
		IJdates=triplet[0]
		JKdates=triplet[1]
		IKdates=triplet[2]

		# Triplet ifg indices
		IJndx=data.pairs.index(IJdates)
		JKndx=data.pairs.index(JKdates)
		IKndx=data.pairs.index(IKdates)

		# Interferograms
		IJ=data.ifgs[IJndx,:,:]
		JK=data.ifgs[JKndx,:,:]
		IK=data.ifgs[IKndx,:,:]

		# Normalize to reference point
		IJ-=IJ[inps.refY,inps.refX]
		JK-=JK[inps.refY,inps.refX]
		IK-=IK[inps.refY,inps.refX]

		# Compute misclosure
		misclosure=IJ+JK-IK
		absMisclosure=np.abs(misclosure)

		# Plot if requested
		if inps.plotTriplets:
			tripletFig=plotTriplets(IJ,JK,IK,misclosure,inps.refX,inps.refY)
			tripletFig.suptitle(triplet)

		# Append to stack
		data.miscStack.append(misclosure)
		data.absMiscStack.append(absMisclosure)

	# Convert lists to 3D arrays
	data.miscStack=np.array(data.miscStack)
	data.absMiscStack=np.array(data.absMiscStack)

	# Cumulative misclosure
	data.cumMisclosure=np.sum(data.miscStack,axis=0)
	data.cumAbsMisclosure=np.sum(data.absMiscStack,axis=0)

	# Use first data as reference date
	data.refDates=[triplet[0][0] for triplet in inps.triplets]
	data.refDatetimes=[datetime.strptime(str(triplet[0][0]),"%Y%m%d") for triplet in inps.triplets]



### PLOTTING FUNCTIONS ---
## Plot triplets
def plotTriplets(IJ,JK,IK,misclosure,refX,refY):
	# Create figure and axes
	tripletFig=plt.figure()

	# Map stats
	IJstats=mapStats(IJ,pctmin=2,pctmax=98)
	JKstats=mapStats(JK,pctmin=2,pctmax=98)
	IKstats=mapStats(IK,pctmin=2,pctmax=98)
	misclosureStats=mapStats(misclosure,pctmin=2,pctmax=98)

	# Plot maps
	cbar_orient='horizontal'
	axIJ=tripletFig.add_subplot(221)
	caxIJ=axIJ.imshow(np.ma.array(IJ,mask=(IJ==IJstats.background)),
		cmap='jet',vmin=IJstats.vmin,vmax=IJstats.vmax)
	tripletFig.colorbar(caxIJ,orientation=cbar_orient)
	axIJ.plot(refX,refY,'ks')

	axJK=tripletFig.add_subplot(222)
	caxJK=axJK.imshow(np.ma.array(JK,mask=(JK==JKstats.background)),
		cmap='jet',vmin=IJstats.vmin,vmax=IJstats.vmax)
	tripletFig.colorbar(caxJK,orientation=cbar_orient)
	axJK.plot(refX,refY,'ks')

	axIK=tripletFig.add_subplot(223)
	caxIK=axIK.imshow(np.ma.array(IK,mask=(IK==IKstats.background)),
		cmap='jet',vmin=IJstats.vmin,vmax=IJstats.vmax)
	tripletFig.colorbar(caxIK,orientation=cbar_orient)
	axIK.plot(refX,refY,'ks')

	axMisclosure=tripletFig.add_subplot(224)
	caxMisclosure=axMisclosure.imshow(np.ma.array(misclosure,mask=(misclosure==misclosureStats.background)),
		cmap='jet',vmin=misclosureStats.vmin,vmax=misclosureStats.vmax)
	tripletFig.colorbar(caxMisclosure,orientation=cbar_orient)
	axMisclosure.plot(refX,refY,'ks')

	# Formatting
	axIJ.set_xticks([]); axIJ.set_yticks([]); axIJ.set_title('I_J')
	axJK.set_xticks([]); axJK.set_yticks([]); axJK.set_title('J_K')
	axIK.set_xticks([]); axIK.set_yticks([]); axIK.set_title('I_K')
	axMisclosure.set_xticks([]); axMisclosure.set_yticks([])
	axMisclosure.set_title('Misclosure')

	return tripletFig


## Plot cumulative misclosure
class misclosureMap:
	def __init__(self,Fig,ax,img,refX,refY,title=None,background=None,vmin=None,vmax=None):
		self.Fig=Fig; self.ax=ax
		self.img=img
		self.refX=refX; self.refY=refY
		self.title=title
		self.background=background
		self.vmin=vmin
		self.vmax=vmax
		self.cbar=None # format colorbar only once

		# Mask background if specified
		if background is not None:
			self.img=np.ma.array(img,mask=(img==background))

	def plotax(self):
		# Plot
		cax=self.ax.imshow(self.img,vmin=self.vmin,vmax=self.vmax,cmap='jet')
		self.ax.plot(self.refX,self.refY,'ks')
		self.ax.set_xticks([]); self.ax.set_yticks([])
		if self.title:
			self.ax.set_title(self.title)
		# Format colorbar once only
		if not self.cbar:
			self.cbar=self.Fig.colorbar(cax,orientation='horizontal')
			self.cbar.set_label('radians')


## Plot misclosure
def plotMisclosure(inps,data):
	## Map stats - clip to inner 96% unless specified
	cumMiscStats=mapStats(data.cumMisclosure,pctmin=2,pctmax=98)
	if inps.miscLims: # Specified limits
		cumMiscStats.vmin=inps.miscLims[0]; cumMiscStats.vmax=inps.miscLims[1]
	cumAbsMiscStats=mapStats(data.cumAbsMisclosure,pctmin=2,pctmax=98)
	if inps.absMiscLim: # Specified upper limit
		cumAbsMiscStats.vmin=0; cumAbsMiscStats.vmax=inps.absMiscLim


	## Plot cumulative misclosure
	title='Cumulative misclosure ({} triplets)'.format(inps.nTriplets)
	figCumMiscMap=plt.figure(); axCumMiscMap=figCumMiscMap.add_subplot(111)
	cumMiscMap=misclosureMap(figCumMiscMap,axCumMiscMap,data.cumMisclosure,inps.refX,inps.refY,title=title,
		background=cumMiscStats.background, vmin=cumMiscStats.vmin,vmax=cumMiscStats.vmax)
	cumMiscMap.plotax()


	## Plot cumulative absolute misclosure
	title='Cumulative absolute misclosure ({} triplets)'.format(inps.nTriplets)
	figCumAbsMiscMap=plt.figure(); axCumAbsMiscMap=figCumAbsMiscMap.add_subplot(111)
	cumAbsMiscMap=misclosureMap(figCumAbsMiscMap,axCumAbsMiscMap,data.cumAbsMisclosure,inps.refX,inps.refY,title=title,
		background=cumAbsMiscStats.background, vmin=cumAbsMiscStats.vmin,vmax=cumAbsMiscStats.vmax)
	cumAbsMiscMap.plotax()

	return cumMiscMap, cumAbsMiscMap



### MISCLOSURE ANALYSIS ---
def plotSeries(name,ax,data,series,timeAxis=False):
	ax.plot(data.refDatetimes,series,'-k.')

	ax.set_ylabel(name+'\nradians')
	if timeAxis is False:
		ax.set_xticklabels([])
	else:
		ax.set_xticks(data.refDatetimes)
		ax.set_xticklabels(data.refDates,rotation=80)
	

## Analyze misclosure stack
def analyzeStack(event):
	print('Stack analysis')

	# Location
	px=event.xdata; py=event.ydata
	px=int(round(px)); py=int(round(py))

	# Report position and cumulative values
	print('px {} py {}'.format(px,py)) # report position
	print('Cumulative misclosure: {}'.format(data.cumMisclosure[py,px]))
	print('Abs cumulative misclosure: {}'.format(data.cumAbsMisclosure[py,px]))

	# Plot query points on maps
	cumMiscMap.ax.cla(); cumMiscMap.plotax() # clear and replot map
	cumMiscMap.ax.plot(px,py,color='k',marker='o',markerfacecolor='w',zorder=3)

	cumAbsMiscMap.ax.cla(); cumAbsMiscMap.plotax()
	cumAbsMiscMap.ax.plot(px,py,color='k',marker='o',markerfacecolor='w',zorder=3)

	# Timeseries
	# Plot misclosure over time
	print('Misclosure: {}'.format(data.miscStack[:,py,px,]))
	miscSeriesAx.cla() # misclosure
	plotSeries('misclosure',miscSeriesAx,data,data.miscStack[:,py,px])
	cumMiscSeriesAx.cla() # cumulative misclosure
	plotSeries('cum. miscl.',cumMiscSeriesAx,data,np.cumsum(data.miscStack[:,py,px]))
	absMiscSeriesAx.cla() # absolute misclosure
	plotSeries('(abs. miscl.)',absMiscSeriesAx,data,data.absMiscStack[:,py,px])
	cumAbsMiscSeriesAx.cla() # cumulative absolute misclosure
	plotSeries('(cum. abs. miscl.)',cumAbsMiscSeriesAx,data,np.cumsum(data.absMiscStack[:,py,px]),
		timeAxis=True)

	# Draw outcomes
	cumMiscMap.Fig.canvas.draw()
	cumAbsMiscMap.Fig.canvas.draw()
	miscSeriesFig.canvas.draw()



### MAIN CALL ---
def main(inps=None):
	## Gather arguments
	inps=cmdLineParse()

	## Load data based on data type
	data=loadData(inps)
	if inps.plotPairs is True:
		from viewingFunctions import plotDatePairs
		plotDatePairs(data.pairs)

	# Plot data if requested
	if inps.plotInputs is True:
		imagettes(data.ifgs,4,5,cmap='jet',downsampleFactor=3,pctmin=2,pctmax=98,background='auto',	titleList=data.pairs,supTitle=None)


	## Formulate valid triplets
	# List of triplets based on all dates
	inps.triplets=createTriplets(data.dates,
		lags=inps.lags,minTime=inps.minTime,maxTime=inps.maxTime,
		verbose=inps.verbose)

	# Validate triplets against list of pairs
	validTriplets=[]
	for triplet in inps.triplets:
		count=0
		for pair in triplet:
			if pair in data.pairs: count+=1
		if count==3: validTriplets.append(triplet)
	inps.triplets=validTriplets
	inps.nTriplets=len(inps.triplets)

	if inps.verbose is True:
		print('{} triplets validated'.format(inps.nTriplets))


	## Detect reference point
	if inps.refY is 'auto' or inps.refX is 'auto':
		inps.refY=np.random.randint(0,data.M,1)[0]
		inps.refX=np.random.randint(0,data.N,1)[0]
	else:
		inps.refX=int(inps.refX); inps.refY=int(inps.refY)

	if inps.verbose is True: print('Reference points: Y {}; X {}'.format(inps.refY,inps.refX))
	

	## Calculate misclosure
	calcMisclosure(inps,data)

	# Plot misclosure
	cumMiscMap,cumAbsMiscMap=plotMisclosure(inps,data)


	## Misclosure analysis
	# Spawn misclosure figure
	miscSeriesFig=plt.figure('Misclosure',figsize=(8,8))
	miscSeriesAx=miscSeriesFig.add_subplot(411)
	cumMiscSeriesAx=miscSeriesFig.add_subplot(412)
	absMiscSeriesAx=miscSeriesFig.add_subplot(413)
	cumAbsMiscSeriesAx=miscSeriesFig.add_subplot(414)

	cumMiscMap.Fig.canvas.mpl_connect('button_press_event',analyzeStack)
	cumAbsMiscMap.Fig.canvas.mpl_connect('button_press_event',analyzeStack)


	plt.show()
