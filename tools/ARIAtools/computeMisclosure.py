#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For a given set of interferograms, compute the cumulative phase
#  misclosure.
#
# Rob Zinke 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### IMPORT MODULES ---
import argparse
import os
from glob import glob
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import register_matplotlib_converters
from osgeo import gdal
from scipy.stats import mode

import logging
from ARIAtools.logger import logger

register_matplotlib_converters()
log = logging.getLogger(__name__)


### PARSER ---
Description='''Compute the cumulative misclosure of phase triplets based on a set of interferograms
saved in the stack/unwrapStack.vrt data set. During triplet computation, values at a reference point
are removed from the interferograms prior to misclosure computation to account for abiguities in the
unwrapped phase that might arise during pairwise computation.

The code works by generating a list of triplets from the available pairs.
From the list of triplets, the routine will:

    1. Formulate a list of viable triplet combinations (IJ, JK, IK)
    2. Subtract the reference point from those interferograms
    3. Compute the disagreement based on IJ+JK-IK
       ... continue through the list of triplets.

Both the misclosure values (positive or negative) and the absolute misclosure values
(always positive) are computed and stored in 3D arrays, where each slice represents a triplet. The
"cumulative" misclosure and absolute misclosure are computed by summing the 3D arrays.

Once the (absolute) misclosure is calculated, the user can view the time history of any pixel by
clicking on the maps.

Thumbnail images of the misclosure associated with any given triplet are saved in the MisclosureFigs
folder. Additionally, georeferenced tiffs of the cumulative misclosure and absolute cumulative
misclosure maps are saved.
'''

Examples='''EXAMPLES

# Using the unwrapStack.vrt to call all interferograms, automatically find a reference point
ariaMisclosure.py -f stack/unwrapStack.vrt

# Provide a predefined reference point
ariaMisclosure.py -f stack/unwrapStack.vrt -refLon 89.358 -refLat 32.621

# Limit triplet selection by time interval (12 days to 48 days)
ariaMisclosure.py -f stack/unwrapStack.vrt --mintime 12 --maxtime 48

'''

def createParser():
    parser = argparse.ArgumentParser(description=Description,
        formatter_class=argparse.RawTextHelpFormatter, epilog=Examples)

    # Input data
    parser.add_argument('-f', '--file', dest='unwFile', type=str, required=True,
        help='ARIA files. Specify the stack/unwrapStack.vrt file, or a wildcard operator in the unwrappedPhase folder (see EXAMPLES)')
    parser.add_argument('-w', '--workdir', dest='workdir', type=str, default='./',
        help='Specify directory to deposit all outputs. Default is local directory where script is launched.')
    parser.add_argument('--coherence', dest='cohFile', type=str, default=None,
        help='Coherence stack for use in automatic reference point selection.')

    parser.add_argument('--startdate', dest='startDate', type=str, default=None,
        help='Start date for data series')
    parser.add_argument('--enddate', dest='endDate', type=str, default=None,
        help='End date for data series')
    parser.add_argument('--exclude-pairs', dest='excludePairs', type=str, default=None,
        help='List of pairs to exclude, e.g., \'20160116_20160101 20171031_20161030\'. This can also be provided in as a text file with one line per date pair.')
    parser.add_argument('--plot-pairs', dest='plotPairs', action='store_true',
        help='Plot the timespans of date pairs')

    # Triplet formulation
    parser.add_argument('--mintime', dest='minTime', type=int, default=None,
        help='Minimum time span of pairs in triplets (days)')
    parser.add_argument('--maxtime', dest='maxTime', type=int, default=None,
        help='Maximum time span of pairs in triplets (days)')
    parser.add_argument('--print-triplets', dest='printTriplets', action='store_true',
        help='Print list of existing triplets (i.e., those included in the data set).')
    parser.add_argument('--plot-triplets', dest='plotTriplets', action='store_true',
        help='Plot existing triplets')

    # Reference point
    parser.add_argument('-refX', dest='refX', type=int, default=None, help='Reference X pixel')
    parser.add_argument('-refY', dest='refY', type=int, default=None, help='Reference Y pixel')
    parser.add_argument('-refLon', dest='refLon', type=float, default=None,
        help='Reference longitude')
    parser.add_argument('-refLat', dest='refLat', type=float, default=None,
        help='Reference latitude')

    # Query point
    parser.add_argument('--queryX', dest='queryX', type=int, default=None, help='Query point X pixel')
    parser.add_argument('--queryY', dest='queryY', type=int, default=None, help='Query point Y pixel')
    parser.add_argument('--queryLon', dest='queryLon', type=float, default=None, help='Query point longitude')
    parser.add_argument('--queryLat', dest='queryLat', type=float, default=None, help='Query point latitude')

    # Vocalization
    parser.add_argument('-v','--verbose', dest='verbose', action='store_true', help='Verbose mode')

    # Misclosure map formatting
    parser.add_argument('--plot-time-intervals', dest='plotTimeIntervals', action='store_true',
        help='Plot triplet intervals in misclosure analysis figure.')

    return parser


def cmdLineParse(iargs = None):
    parser = createParser()
    return parser.parse_args(args=iargs)



### STACK OBJECT ---
class MisclosureStack:
    '''Class for loading and storing stack data for phase triplet
    misclosure analysis.
    '''
    ## Load data
    def __init__(self,
            unwFile,
            cohFile=None,
            workdir='./',
            startDate=None, endDate=None,
            excludePairs=None,
            verbose=False):
        '''Initialize object. Store essential info for posterity.
        loadStackData() - Load data from unwrapStack.vrt using gdal
        formatDates() - Determine the IFG pairs and list of unique dates
        from the data set
        formatExcludePairs() - Load and format pairs to exclude, if provided
        '''
        if verbose == True:
            print('Initializing stack for misclosure analysis')

        # Verbosity
        self.verbose = verbose
        if self.verbose:
            logger.setLevel(logging.DEBUG)

        # Files and directories
        self.unwFile = os.path.abspath(unwFile)
        self.workdir = os.path.abspath(workdir)

        # Check if output directory exists
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)

        # Read stack data and retrieve list of dates
        self.__loadUnwStack__()

        # Load coherence file if specified
        if cohFile is not None:
            self.__loadCohStack__(cohFile)

        # Format dates
        self.__formatDates__(startDate, endDate)

        # Format pairs to exclude, if provided
        self.__formatExcludePairs__(excludePairs)

    def __loadUnwStack__(self):
        '''Load data from unwrapStack.vrt file.'''
        # Open data set
        self.unwStack = gdal.Open(self.unwFile, gdal.GA_ReadOnly)

        # Format extent
        self.__formatGeoInfo__()

        # Report if requested
        log.debug('%s bands detected', self.unwStack.RasterCount)

    def __loadCohStack__(self, cohFile):
        '''Load data from cohStack.vrt file.'''
        cohFile = os.path.abspath(cohFile)

        # Load GDAL data set
        cohStack = gdal.Open(cohFile, gdal.GA_ReadOnly)

        # Image bands
        imgs = np.array([cohStack.GetRasterBand(i+1).ReadAsArray() for i in \
                range(cohStack.RasterCount)])

        # Mean coherence
        self.meanCoh = np.mean(imgs, axis=0)

        log.debug('Coherence stack loaded and averaged')


    ## Geographic information
    def __formatGeoInfo__(self):
        '''Parse the spatial metadata associated with the GDAL data set.'''
        # Image sizes
        self.N = self.unwStack.RasterXSize
        self.M = self.unwStack.RasterYSize

        # Get projection for later
        self.proj = self.unwStack.GetProjection()

        # Get geographic transform
        self.tnsf = self.unwStack.GetGeoTransform()

        # Parse geotransform
        left, xstep, xskew, top, yskew, ystep = self.tnsf

        # Re-format geotransform as matrix
        self.tnsfMatrix = np.array([[xstep, yskew],
                                    [xskew, ystep]])

        # Origin coordinates as vector
        self.tnsfOrigin = np.array([[left, top]]).T

        # Plot extent
        right = left + xstep*self.N
        bottom = top + ystep*self.M
        self.extent = (left, right, bottom, top)

    def xy2lola(self, x, y):
        '''Convert X/Y to lon/lat.'''
        # Reshape points as vector
        p = np.array([[x, y]]).T

        # Calculate geographic position
        lon, lat = self.tnsfMatrix.dot(p) + self.tnsfOrigin

        # Flatten arrays
        lon = lon.flatten()
        lat = lat.flatten()

        return lon, lat

    def lola2xy(self, lon, lat):
        '''Convert lon/lat coordinates to XY.'''
        # Reshape points as vector
        L = np.array([[lon, lat]]).T

        # Calculate point coordinates
        px, py = np.linalg.inv(self.tnsfMatrix).dot(L - self.tnsfOrigin)

        # Convert pixel coordinates to integers
        px = int(px.flatten())
        py = int(py.flatten())

        return px, py


    ## Date and date pair formatting
    def __formatDates__(self, startDate, endDate):
        '''Retrieve list of date pairs and unique dates (epochs).
        The "pairs" attribute is a formatted list of interferogram date
        pairs **in the order in which they are stored in the unwrapStack.vrt**
        file. This list should not be modified.

        Constrain the list of epochs available for triplet determination
        using the "startDate" and "endDate" provided in the __init__
        function.
        '''
        log.debug('Formatting dates')

        # PairNames - list of pairs composing the data set in the order
        # they are written
        pairNames = [os.path.basename(fname) for fname in \
                self.unwStack.GetFileList()]

        # Remove extra file name
        pairNames.remove('unwrapStack.vrt')

        # Remove extensions
        pairNames = [pairName.strip('.vrt') for pairName in pairNames]

        # Convert pair name strings to datetime objects
        self.datePairs = self.__pairNames2datePairs__(pairNames)

        # Get unique dates from date pairs
        self.dates = self.__uniqueDatesFromPairs__(self.datePairs)

        # Format start and end dates
        self.__applyStartEndDates__(startDate, endDate)

        # Number of dates within start-end range
        self.nDates = len(self.dates)

        log.debug('%s unique dates detected', self.nDates)

    def __pairNames2datePairs__(self, pairNames):
        '''Convert list of pairs in format ['master_slave','master_slave',...]
        to dates in format [[master, slave], [master, slave], ...]
        '''
        datePairs = [self.__strPair2datePair__(pairName) for pairName \
                in pairNames]

        return datePairs

    def __strPair2datePair__(self,pair):
        '''Convert pair in format 'master_slave' to date in format
        [master, slave]
        '''
        pair = pair.split('_')
        masterDate = datetime.strptime(pair[0], '%Y%m%d')
        secondaryDate = datetime.strptime(pair[1], '%Y%m%d')
        datePair = [masterDate, secondaryDate]
        return datePair

    def __datePair2strPair__(self, datePair):
        '''Convert pair in format [master, slave] to date in format
        [master_slave]
        '''
        masterStr = datePair[0].strftime('%Y%m%d')
        secondaryStr = datePair[1].strftime('%Y%m%d')
        strPair = '{:s}_{:s}'.format(masterStr, secondaryStr)
        return strPair

    def __uniqueDatesFromPairs__(self, datePairs):
        '''Get a list of unique datetimes representing epochs of
        acqusition.
        '''
        # List of all dates
        dates = []
        [dates.extend(pair) for pair in datePairs]

        # Filter for unique dates
        dates = list(set(dates))

        # Sort oldest-youngest
        dates.sort()

        return dates

    def __applyStartEndDates__(self, startDate, endDate):
        '''Format the start and end dates as datetime objects.
        Trim list of dates to start and end date limits.'''

        # Start date
        if startDate:
            # Format start date as datetime
            startDate = datetime.strptime(startDate, '%Y%m%d')
        else:
            # Use first date
            startDate = self.dates[0]

        # End date
        if endDate:
            # Format end date as datetime
            endDate = datetime.strptime(endDate, '%Y%m%d')
        else:
            # Use last date
            endDate = self.dates[-1]

        log.debug('Start date: %s; end date %s', startDate, endDate)

        # Crop list of dates to start and end
        self.dates = [date for date in self.dates if date >= startDate]
        self.dates = [date for date in self.dates if date <= endDate]

        # Create array of dates for plotting
        self.__createDateAxis__(startDate, endDate)

    def __createDateAxis__(self, startDate, endDate):
        '''Create an array of dates for plotting misclosure values.'''
        # Plot x-ticks
        self.dateTicks = pd.date_range(
                            startDate - timedelta(days=30),
                            endDate + timedelta(days=30),
                            freq='MS')

        # Plot x-tick labels
        self.dateTickLabels = [date.strftime('%Y-%m') for date in
                                self.dateTicks]

    def __formatExcludePairs__(self, excludePairs):
        '''Check that exclude dates are in one of two formats:
        1. a string containing the pairs in YOUNGER_OLDER format,
           space-separated
        2. a .txt file with lines of the same formatting

        Formatting should match "pair" formatting: [[master,slave]]
        '''
        if excludePairs is not None:
            # Determine whether list or text file
            if excludePairs.endswith('.txt'):
                # Treat as text file with list
                with open(excludePairs, 'r') as exclFile:
                    excludePairs = exclFile.readlines()
                    excludePairs = [pair.strip('\n') for pair in excludePairs]
                    self.excludePairs = self.__pairNames2datePairs__(excludePairs)
            else:
                # Treat as list - split at spaces
                excludePairs = excludePairs.split()
                self.excludePairs = self.__pairNames2datePairs__(excludePairs)
        else:
            # Include as empty list
            self.excludePairs = []

    def plotPairs(self):
        '''Plot the timespans of interferogram pairs.'''
        # Copy the list of pairs and sort them in time order
        datePairs = self.datePairs[:]  # copy to separate object
        datePairs.sort(key=lambda s: s[1])  # sort by secondary date

        # Spawn figure and axis
        pairFig, pairAx = plt.subplots()

        # Loop through date pairs and plot them in time
        for n, pair in enumerate(datePairs):
            # Color based on whether included or excluded
            if pair not in self.excludePairs:
                color = 'k'
                label = 'valid pair'
                linestyle = '-'
            else:
                color = (0.6, 0.6, 0.6)
                label = 'excluded pair'
                linestyle = '--'

            # Convert to datetime format
            pairAx.plot([pair[0],pair[1]], [n,n],
                color=color, label=label, linestyle=linestyle)

        # Format x-axis
        pairAx.set_xticks(self.dateTicks)
        pairAx.set_xticklabels(self.dateTickLabels, rotation=90)

        # Legend
        handles, labels=pairAx.get_legend_handles_labels()
        uniqueLabels = dict(zip(labels,handles))
        pairAx.legend(uniqueLabels.values(), uniqueLabels.keys(),
            bbox_to_anchor=(0.005, 0.99), loc='upper left', borderaxespad=0.)

        # Other formatting
        pairAx.set_yticks([])
        pairAx.set_title('IFG pairs')
        pairFig.tight_layout()


    ## Create triplet list
    def createTriplets(self, minTime=None, maxTime=None, printTriplets=False):
        '''Create a list of triplets given the date list and user-specified
        parameters.

        First generate a list of all possible triplets based on the
        available dates. Then validate that list across the list of
        existing pairs.

        The stack object retains an ordered list of dates in both
        YYYYMMDD format and datetime format, and pairList, based on the
        __format_dates__ function.
        '''
        log.debug('Creating list of triplets')

        # Create a list of all possible triplets
        self.__createAllTriplets__()

        # Remove triplets with pairs in "exclude pairs" list
        self.__checkExcludedTriplets__()

        # Check that pairs meet time requirements
        self.__checkTripletsMinTime__(minTime)

        self.__checkTripletsMaxTime__(maxTime)

        # Check against list of existing interferograms
        self.__checkTripletsExist__()

        # Finished sorting
        self.nTriplets = len(self.triplets)

        log.debug('%s valid triplets identified', self.nTriplets)

        # Retrieve triplet reference dates
        self.__retreiveTripletReferenceDates__()

        # Save triplets to text file
        self.__saveTriplets__()

        # Print triplets to screen if requested
        if printTriplets == True:
            self.__printTriplets__()

    def __createAllTriplets__(self):
        '''Loop through dates to create all possible triplet combinations.'''
        log.debug('Listing all possible triplets')

        # Create empty list
        self.triplets = []

        # Loop through first date in ordered list (date 1)
        for i in range(self.nDates-2):
            # Loop through second dates (start from first date after date 1)
            for j in range(i+1, self.nDates-1):
                # Loop through third dates (start from first date after date 2)
                for k in range(j+1, self.nDates):
                    dateI = self.dates[i]  # first date in sequence
                    dateJ = self.dates[j]  # second date in sequence
                    dateK = self.dates[k]  # third date in sequence
                    self.triplets.append([
                            [dateJ, dateI],
                            [dateK, dateJ],
                            [dateK, dateI]
                            ])

    def __checkExcludedTriplets__(self):
        '''Check triplet list against excluded pairs list. Remove the
        triplet if any of the pairs is listed in "exclude pairs".
        '''
        log.debug('Checking triplets against excluded pairs')

        # Empty list of non-excluded triplets
        validTriplets = []

        # Loop through all possible triplets
        for triplet in self.triplets:
            # Invalid date pairs in triplet
            invalidTriplets = 0  # reset counter
            for datePair in triplet:
                # Check against list of pairs to exclude
                if datePair in self.excludePairs:
                    invalidTriplets += 1

            # If no pairs are excluded, append to the valid triplets list
            if invalidTriplets == 0:
                validTriplets.append(triplet)

        # Update triplets list
        self.triplets = validTriplets

    # Check triplets against minTime
    def __checkTripletsMinTime__(self, minTime):
        '''Check that all pairs in a triplet are longer in duration
        than the minimum time interval specified.
        '''
        log.debug('Checking triplet pairs are longer than %s days', minTime)

        if minTime:
            validTriplets = []
            for triplet in self.triplets:
                # Determine intervals between dates in days
                intervals = [(datePair[0]-datePair[1]).days for \
                            datePair in triplet]

                # Check against minimum allowable time interval
                if min(intervals) >= minTime:
                    validTriplets.append(triplet)

            # Update triplets list
            self.triplets = validTriplets

    # Check triplets against maxTime
    def __checkTripletsMaxTime__(self, maxTime):
        '''Check that all pairs in a triplet are shorter in duration than
        the maximum time interval specified.
        '''
        log.debug('Checking triplet pairs are shorter than %s days', maxTime)

        if maxTime:
            validTriplets = []
            for triplet in self.triplets:
                # Determine intervals between dates in days
                intervals = [(datePair[0]-datePair[1]).days for \
                            datePair in triplet]

                # Check against maximum allowable time interval
                if max(intervals) <= maxTime:
                    validTriplets.append(triplet)

            # Update triplets list
            self.triplets = validTriplets

    # Check triplets exist
    def __checkTripletsExist__(self):
        '''Check list of all possible triplets against the list of pairs
        that actually exist.
        '''
        log.debug('Checking triplets provided in data set.')

        existingTriplets = []
        for triplet in self.triplets:
            # Reset count of existing pairs
            existing = 0

            # Check that each pair of the triplet has a corresponding
            # interferogram
            for datePair in triplet:
                if datePair in self.datePairs:
                    # Update if IFG exists
                    existing += 1
            if existing == 3:
                existingTriplets.append(triplet)

        # Update triplet list
        self.triplets = existingTriplets

    def __retreiveTripletReferenceDates__(self):
        '''Create a list of reference dates with each triplet.'''
        log.debug('Retreiving list of reference dates from valid triplet')

        # Reference dates
        self.tripletRefDates = []

        # Loop through triplets
        for triplet in self.triplets:
            # All dates in triplet
            tripletDates = []
            [tripletDates.extend(datePair) for datePair in triplet]

            # Triplet unique dates
            tripletDates = list(set(tripletDates))

            # Sort earliest to latest
            tripletDates.sort(key=lambda date: date.strftime("%Y%m%d"))

            # Append to list
            self.tripletRefDates.append(tripletDates)

    def __saveTriplets__(self):
        '''Save the list of valid triplets to a text file.'''
        with open(os.path.join(self.workdir, 'ValidTriplets.txt'), 'w') \
                as tripletFile:
            # Loop through valid triplets
            for triplet in self.triplets:
                strPair = [self.__datePair2strPair__(pair) for pair in triplet]
                tripletFile.write('{}\n'.format(strPair))

    def __printTriplets__(self):
        '''Print the list of valid triplets.'''
        log.info('Existing triplets:')

        # Loop through triplets
        for triplet in self.triplets:
            log.info([self.__datePair2strPair__(pair) for pair in triplet])

        # Final statistic
        log.info('%s existing triplets found based on search criteria',
            self.nTriplets)

    def plotTriplets(self):
        '''Plot triplets.'''
        # Setup figure
        tripletFig, tripletAx = plt.subplots()

        # Plot triplets
        for i in range(self.nTriplets):
            tripletAx.plot(self.triplets[i], [i,i,i], 'k', marker='o')

        # Format x-axis
        tripletAx.set_xticks(self.dateTicks)
        tripletAx.set_xticklabels(self.dateTickLabels, rotation=90)

        # Other formatting
        tripletAx.set_yticks([])
        tripletAx.set_title('Triplets in data set')
        tripletFig.tight_layout()


    ## Compute misclosure
    def computeMisclosure(self, refXY=None, refLoLa=None):
        '''Compute the misclosure of the phase triplets.

        A common reference point is required because the ifgs are not
        coregistered.
        '''
        log.debug('Computing misclosure')

        # Create background value mask
        self.__createMask__()

        # Determine reference point
        self.__referencePoint__(refXY, refLoLa)

        # Misclosure placeholders
        self.netMscStack = []
        self.absMscStack = []

        # Compute phase triplets
        log.debug('Calculating triplet misclosure')

        for triplet in self.triplets:
            # Triplet date pairs
            JIdates = triplet[0]
            KJdates = triplet[1]
            KIdates = triplet[2]

            # Triplet indices - add 1 because raster bands start at 1
            JIndx = self.datePairs.index(JIdates)+1
            KJndx = self.datePairs.index(KJdates)+1
            KIndx = self.datePairs.index(KIdates)+1

            # Interferograms
            JI = self.unwStack.GetRasterBand(JIndx).ReadAsArray()
            KJ = self.unwStack.GetRasterBand(KJndx).ReadAsArray()
            KI = self.unwStack.GetRasterBand(KIndx).ReadAsArray()

            # Normalize to reference point
            JI -= JI[self.refY, self.refX]
            KJ -= KJ[self.refY, self.refX]
            KI -= KI[self.refY, self.refX]

            # Compute (abs)misclosure
            netMisclosure = JI + KJ - KI
            absMisclosure = np.abs(netMisclosure)

            # Append to stack
            self.netMscStack.append(netMisclosure)
            self.absMscStack.append(absMisclosure)

        # Convert lists to 3D arrays
        self.netMscStack = np.array(self.netMscStack)
        self.absMscStack = np.array(self.absMscStack)

        # Cumulative misclosure
        self.cumNetMisclosure = np.sum(self.netMscStack, axis=0)
        self.cumAbsMisclosure = np.sum(self.absMscStack, axis=0)

        # Apply mask
        self.cumNetMisclosure[self.mask==0] = 0
        self.cumAbsMisclosure[self.mask==0] = 0

    def __createMask__(self):
        '''Create a mask based on the nodata value.'''
        log.debug('Creating no data mask')

        # Retrieve first image from stack
        img = self.unwStack.GetRasterBand(1).ReadAsArray()

        # Mask no data values
        self.mask = np.ones((self.M, self.N))
        self.mask[img==0] = 0

    def __referencePoint__(self, refXY, refLoLa):
        '''Determine the reference point in XY coordinates. The reference
        point can be automatically or manually selected by the user and
        is subtracted from each interferogram.

        The point can be given in pixels or lon/lat coordinates. If
        given in Lat/Lon, determine the location in XY, and vice-versa.
        '''
        log.debug('Determining reference point...')

        if refLoLa.count(None) == 0:
            # Record reference lon/lat
            self.refLon = refLoLa[0]
            self.refLat = refLoLa[1]

            # Determine the x/y coordinates from the given lon/lat
            x, y = self.lola2xy(self.refLon, self.refLat)
            self.refX = int(x)
            self.refY = int(y)

            log.debug('Reference point given as: X %s / Y %s; Lon %s / Lat %s',
                                 self.refX, self.refY, self.refLon, self.refLat)

        elif refXY.count(None) == 0:
            # Record reference x/y
            self.refX = refXY[0]
            self.refY = refXY[1]

            # Determine the lon/lat coordinates from the given x/y
            self.refLon, self.refLat = self.xy2lola(self.refX, self.refY)
            log.debug('Reference point given as: X %s / Y %s; Lon %.4f / Lat %.4f',
                                 self.refX, self.refY, self.refLon, self.refLat)

        else:
            # Use a random reference point
            self.__autoReferencePoint__()

    def __autoReferencePoint__(self):
        '''Use the coherence stack to automatically determine a suitable
        reference point.
        '''
        # Try to determine the reference point using coherence map
        if self.__autoReferenceCoherence__() == False:
            # Otherwise, resort to picking a random point
            self.__autoReferenceRandom__()

    def __autoReferenceCoherence__(self):
        '''Attempt to find cohStack.vrt file based on unwStack directory.
        If found, determine a random high-coherence point.
        If file not found, return False.
        '''
        log.debug('Attempting to find high-coherence reference point')

        # Check if coherence data have already been loaded
        if not hasattr(self, 'meanCoh'):
            # Check unwStack folder for cohStack
            dirName = os.path.dirname(self.unwFile)

            try:
                # Standard coherence file name
                cohFile = os.path.join(dirName, 'cohStack.vrt')

                # Load average coherence map
                self.__loadCohStack__(cohFile)

            except:
                log.debug('Could not automatically find coherence stack')
                return False

            # Randomly sample for high-coherence values
            nTries = 10000

            while nTries > 0:
                # Pick random points
                x = np.random.randint(0, self.N)
                y = np.random.randint(0, self.M)

                # Check coherence at those points
                if self.meanCoh[y,x] >= 0.7:
                    # Point is high-coherence, stop trying
                    break

                # Decrement counter
                nTries -= 1

            # Assign reference x/y
            self.refX, self.refY = x, y

            # Convert to lon/lat
            self.refLon, self.refLat = self.xy2lola(self.refX, self.refY)

            log.debug('Reference point chosen based on coherence as: X %s / Y %s; Lon %.4f / Lat %.4f.',
                                     self.refX, self.refY, self.refLon, self.refLat)

            return True

    def __autoReferenceRandom__(self):
        '''Choose random, non-masked point for reference.'''
        log.debug('Choosing random reference point')

        # Number of tries
        nTries = 100

        while nTries > 0:
            # Pick random points
            x = np.random.randint(0, self.N)
            y = np.random.randint(0, self.M)

            # Check if that point is masked
            if self.mask[y,x] == 1:
                # Point is not masked, stop trying
                break

            # Decrement counter
            nTries -= 1

        # Assign reference x/y
        self.refX, self.refY = x, y

        # Convert to lon/lat
        self.refLon, self.refLat = self.xy2lola(self.refX, self.refY)

        log.debug('Reference point chosen randomly as: X %s / Y %s; Lon %.4f / Lat %.4f.',
                                 self.refX, self.refY, self.refLon, self.refLat)


    ## Plot misclosure
    def plotMisclosure(self, queryXY=None, queryLoLa=None,
            plotTimeIntervals=False):
        '''Map-view plots of cumulative misclosure.'''
        log.debug('Begin misclosure analysis')

        # Parameters
        self.plotTimeIntervals = plotTimeIntervals

        # Set up interactive plots
        self.netMscFig, self.netMscAx = plt.subplots()
        self.absMscFig, self.absMscAx = plt.subplots()

        # Mask arrays to ignore no data values
        self.cumNetMisclosure = np.ma.array(self.cumNetMisclosure,
            mask=(self.cumNetMisclosure==0))
        self.cumAbsMisclosure = np.ma.array(self.cumAbsMisclosure,
            mask=(self.cumAbsMisclosure==0))

        # Auto-detect values for clipping color scale
        self.cumNetMscClips = self.__imgClipValues__(self.cumNetMisclosure)
        self.cumAbsMscClips = self.__imgClipValues__(self.cumAbsMisclosure)

        # Plot maps
        cax = self.__plotCumNetMisclosure__()
        cbar = self.netMscFig.colorbar(cax, orientation='vertical')
        cbar.set_label('cum. misclosure (radians)')

        cax = self.__plotCumAbsMisclosure__()
        cbar = self.absMscFig.colorbar(cax,orientation='vertical')
        cbar.set_label('cum. abs. misclosure (radians)')

        # Plot timeseries points
        self.mscSeriesFig = plt.figure('Misclosure', figsize=(8,8))
        self.netMscSeriesAx = self.mscSeriesFig.add_subplot(411)
        self.cumNetMscSeriesAx = self.mscSeriesFig.add_subplot(412)
        self.absMscSeriesAx = self.mscSeriesFig.add_subplot(413)
        self.cumAbsMscSeriesAx = self.mscSeriesFig.add_subplot(414)

        # Pre-specified query points
        if (queryLoLa.count(None) == 0) or (queryXY.count(None) == 0):
            self.__misclosureQuery__(queryXY, queryLoLa)

        # Link canvas to plots for interaction
        self.netMscFig.canvas.mpl_connect('button_press_event',
            self.__samplePixel__)
        self.absMscFig.canvas.mpl_connect('button_press_event',
            self.__samplePixel__)

    def __imgClipValues__(self, img):
        '''Find values at which to clip the images (min/max) based on
        histogram percentiles.
        '''
        # Determine clip values based on percentiles
        clipValues={}
        clipValues['min'], clipValues['max'] = \
            np.percentile(img.compressed(), (2, 98))

        return clipValues

    def __plotCumNetMisclosure__(self):
        '''Plot cumulative misclosure map.'''
        # Plot map
        cax = self.netMscAx.imshow(self.cumNetMisclosure, cmap='plasma',
            vmin=self.cumNetMscClips['min'], vmax=self.cumNetMscClips['max'],
            zorder=1)

        # Plot reference point
        self.netMscAx.plot(self.refX, self.refY, 'ks', zorder=2)

        # Format axis
        self.netMscAx.set_title('Cumulative misclosure')

        return cax

    def __plotCumAbsMisclosure__(self):
        '''Plot cumulative absolute misclosure map.'''
        # Plot map
        cax = self.absMscAx.imshow(self.cumAbsMisclosure, cmap='plasma',
            vmin=self.cumAbsMscClips['min'], vmax=self.cumAbsMscClips['max'],
            zorder=1)

        # Plot reference point
        self.absMscAx.plot(self.refX, self.refY, 'ks', zorder=2)

        # Format axis
        self.absMscAx.set_title('Cumulative absolute misclosure')

        return cax

    def __plotSeries__(self, ax, data, title):
        '''Plot misclosure timeseries.'''
        # Plot data
        if self.plotTimeIntervals == False:
            ax.plot([tripletRefDate[1] for tripletRefDate in self.tripletRefDates],
                data, '-k.')
        else:
            for n in range(self.nTriplets):
                ax.plot([self.tripletRefDates[n][0],
                        self.tripletRefDates[n][2]],
                        [data[n], data[n]],'k')
                ax.plot(self.tripletRefDates[n][1], data[n], 'ko')

        # Formatting
        ax.set_xticks(self.dateTicks)
        ax.set_xticklabels([])
        ax.set_ylabel(title)


    ## Misclosure analysis
    def __misclosureQuery__(self, queryXY=None, queryLoLa=None):
        '''Show the time history of each pixel based on pre-specified
        selection.
        '''
        log.debug('Pre-specified query point...')

        # Convert bewteen lon/lat and image coordinates
        if queryLoLa.count(None) == 0:
            # Determine the XY coordinates from the given lon/lat
            qLon, qLat = queryLoLa
            qx, qy = self.lola2xy(queryLoLa[0], queryLoLa[1])

        elif queryXY.count(None) == 0:
            # Use the provided XY coordinates
            qx, qy = queryXY
            qLon, qLat = self.xy2lola(queryXY[0], queryXY[1])

        log.debug('Query point: X %s / Y %s; Lon %.4f / Lat %.4f',
            qx, qy, qLon, qLat)

        # Plot query points on map
        self.netMscAx.plot(qx, qy,
            color='k', marker='o', markerfacecolor='w', zorder=3)
        self.absMscAx.plot(qx, qy,
            color='k', marker='o', markerfacecolor='w', zorder=3)

        # Plot misclosure over time
        self.__plotSeries__(self.netMscSeriesAx,
                            self.netMscStack[:,qy,qx],
                            'misclosure')
        self.__plotSeries__(self.cumNetMscSeriesAx,
                            np.cumsum(self.netMscStack[:,qy,qx]),
                            'cum. miscl.')
        self.__plotSeries__(self.absMscSeriesAx,
                            self.absMscStack[:,qy,qx],
                            'abs. miscl')
        self.__plotSeries__(self.cumAbsMscSeriesAx,
                            np.cumsum(self.absMscStack[:,qy,qx]),
                            'cum. abs. miscl.')

        # Format x-axis
        self.cumAbsMscSeriesAx.set_xticks(self.dateTicks)
        self.cumAbsMscSeriesAx.set_xticklabels(self.dateTickLabels,
            rotation=90)

    def __samplePixel__(self, event):
        '''Show the time history of each pixel based on interactive map.'''
        log.debug('Sampling point')

        # Retrieve values from map
        px = event.xdata
        py = event.ydata
        px = int(round(px))
        py = int(round(py))

        # Convert pixels to lat/lon
        lon, lat = self.xy2lola(px, py)

        # Report position and cumulative values
        log.info('px %s py %s', px, py)
        log.info('lon %s lat %s', lon, lat)
        log.info('Cumulative misclosure: %s', self.cumNetMisclosure[py,px])
        log.info('Abs cumulative misclosure: %s', self.cumAbsMisclosure[py,px])

        # Plot query points on maps
        self.netMscAx.cla()
        self.__plotCumNetMisclosure__()
        self.netMscAx.plot(px, py,
            color='k', marker='o', markerfacecolor='w', zorder=3)

        self.absMscAx.cla()
        self.__plotCumAbsMisclosure__()
        self.absMscAx.plot(px, py,
            color='k', marker='o', markerfacecolor='w', zorder=3)

        # Plot misclosure over time
        log.info('Misclosure: %s', self.netMscStack[:,py,px])
        self.netMscSeriesAx.cla()  # misclosure
        self.__plotSeries__(self.netMscSeriesAx,
                            self.netMscStack[:,py,px],
                            'misclosure')

        self.cumNetMscSeriesAx.cla()  # cumulative misclosure
        self.__plotSeries__(self.cumNetMscSeriesAx,
                            np.cumsum(self.netMscStack[:,py,px]),
                            'cum. miscl.')

        self.absMscSeriesAx.cla()  # absolute misclosure
        self.__plotSeries__(self.absMscSeriesAx,
                            self.absMscStack[:,py,px],
                            'abs. miscl')

        self.cumAbsMscSeriesAx.cla()  # cumulative absolute misclosure
        self.__plotSeries__(self.cumAbsMscSeriesAx,
                            np.cumsum(self.absMscStack[:,py,px]),
                            'cum. abs. miscl.')

        # Format x-axis
        self.cumAbsMscSeriesAx.set_xticks(self.dateTicks)
        self.cumAbsMscSeriesAx.set_xticklabels(self.dateTickLabels,
            rotation=90)

        # Draw outcomes
        self.netMscFig.canvas.draw()
        self.absMscFig.canvas.draw()
        self.mscSeriesFig.canvas.draw()


    ## Plot triplet misclosure maps
    def plotTripletMaps(self):
        '''Plot the misclosure measurements for each triplet to figures.'''
        log.debug('Saving incremental misclosure maps to image files')

        # Parameters
        subplotDims = (2, 2)
        maxSubplots = subplotDims[0]*subplotDims[1]

        # Output directory/subdirectory
        self.figdir = os.path.join(self.workdir, 'MisclosureFigs')
        try:
            os.mkdir(self.figdir)
        except:
            pass

        figNb = 0  # start figure counter
        plotNb = 1  # start subplot counter
        for i in range(self.nTriplets):
            # Plot number and subplot position
            if plotNb % maxSubplots == 1:
                # Spawn new figure
                figNb += 1  # update figure counter
                plotNb = 1  # reset subplot counter
                Fig = plt.figure(figsize=(8,6))
            ax = Fig.add_subplot(subplotDims[0], subplotDims[1], plotNb)

            # Format misclosure map
            mscMap = np.ma.array(self.netMscStack[i,:,:], mask=(self.mask==0))
            mscMapClips = self.__imgClipValues__(mscMap)

            # Plot misclosure
            cax = ax.imshow(mscMap, cmap='plasma',
                vmin=mscMapClips['min'], vmax=mscMapClips['max'])

            # Format axis
            ax.set_xticks([])
            ax.set_yticks([])
            Fig.colorbar(cax, orientation='horizontal')
            dates = [date.strftime('%Y%m%d') for date in self.tripletRefDates[i]]
            ax.set_title('_'.join(dates))

            if plotNb == maxSubplots:
                # Save old figure
                Fig.suptitle('Phase misclosure (radians)')
                figname = 'MisclosureValues_fig{:d}.png'.format(figNb)
                figpath = os.path.join(self.figdir, figname)
                Fig.savefig(figpath, dpi=300)

            plotNb += 1

            # Save final figure
            Fig.suptitle('Phase misclosure (radians)')
            figname = 'MisclosureValues_fig{:d}.png'.format(figNb)
            figpath = os.path.join(self.figdir, figname)
            Fig.savefig(figpath, dpi=300)


    ## Save cumulative misclosure plots to geotiffs
    def saveCumMisclosure(self):
        '''Save cumulative (/absolute) misclosure plots to georeferenced
        tiff files. Use metadata from unwrapStack.vrt file.
        '''
        log.debug('Saving misclosure maps to geotiffs')

        # Fix background values
        self.cumNetMisclosure[self.mask==0] = 0
        self.cumAbsMisclosure[self.mask==0] = 0

        # Save cumulative misclosure
        cumNetMscSavename = os.path.join(self.figdir, 'CumulativeMisclosure.tif')
        self.__saveGeoTiff__(cumNetMscSavename, self.cumNetMisclosure)

        # Save cumulative absolute misclosure
        cumAbsMscSavename = os.path.join(self.figdir, 'CumulativeAbsoluteMisclosure.tif')
        self.__saveGeoTiff__(cumAbsMscSavename, self.cumAbsMisclosure)

    def __saveGeoTiff__(self, savename, img):
        '''Template for saving geotiffs.'''
        driver = gdal.GetDriverByName('GTiff')
        DSout = driver.Create(savename,
                self.N, self.M,
                1, gdal.GDT_Float32)
        DSout.GetRasterBand(1).WriteArray(img)
        DSout.GetRasterBand(1).SetNoDataValue(0)
        DSout.SetProjection(self.unwStack.GetProjection())
        DSout.SetGeoTransform(self.unwStack.GetGeoTransform())
        DSout.FlushCache()



### MAIN CALL ---
def main(inps=None):
    ## Gather arguments
    inps = cmdLineParse()


    ## Load data based on data type
    dataStack = MisclosureStack(
            unwFile=inps.unwFile,
            cohFile=inps.cohFile,
            workdir=inps.workdir,
            startDate=inps.startDate, endDate=inps.endDate,
            excludePairs=inps.excludePairs,
            verbose=inps.verbose)

    # Plot pairs if requested
    if inps.plotPairs == True:
        dataStack.plotPairs()


    ## Create list of triplets
    dataStack.createTriplets(minTime=inps.minTime, maxTime=inps.maxTime,
        printTriplets=inps.printTriplets)

    # Plot triplets if requested
    if inps.plotTriplets == True:
        dataStack.plotTriplets()


    ## Compute misclosure
    dataStack.computeMisclosure(refXY=[inps.refX, inps.refY],
        refLoLa=[inps.refLon, inps.refLat])

    # Plot and analyze data
    dataStack.plotMisclosure(
        queryXY=[inps.queryX, inps.queryY],
        queryLoLa=[inps.queryLon, inps.queryLat],
        plotTimeIntervals=inps.plotTimeIntervals)
    plt.show()

    # Save misclosure map for each triplet to figures
    dataStack.plotTripletMaps()

    # Save misclosure maps to geotiffs
    dataStack.saveCumMisclosure()
