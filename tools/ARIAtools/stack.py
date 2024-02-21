# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For a given set of interferograms, compute the cumulative phase
#  misclosure.
#
# Rob Zinke 2020
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# IMPORT MODULES ---
import os
import logging
import datetime
import scipy.stats

import osgeo
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)

# STACK OBJECT ---
class Stack:
    '''
    Class for loading and storing stack data.
    '''
    # Load data
    def __init__(
            self, imgfile, workdir='./', startDate='20140615', endDate=None,
            excludePairs=None, verbose=False):
        '''
        Initialize object. Store essential info for posterity.
        loadStackData() - Load data from unwrapStack.vrt using osgeo.gdal.
        formatDates() - Determine the IFG pairs and list of unique dates
        from the data set.

        formatExcludePairs() - Load and format pairs to exclude, if provided.
        '''
        # Files and directories
        self.imgfile = os.path.abspath(imgfile)
        self.basename = os.path.basename(self.imgfile)
        self.imgdir = os.path.dirname(self.imgfile)
        self.workdir = os.path.abspath(workdir)
        self.verbose = verbose

        # Check if output directory exists
        if not os.path.exists(self.workdir):
            os.mkdir(self.workdir)

        # Dates and pairs
        self.startDate = datetime.datetime.strptime(startDate, '%Y%m%d')
        if not endDate:
            self.endDate = datetime.datetime.now()

        else:
            self.endDate = datetime.datetime.strptime(endDate, '%Y%m%d')

        self.excludePairs = excludePairs

        # Read stack data and retrieve list of dates
        self.__loadStackData__()

        # Format dates
        self.__formatDates__()

        # Format pairs to exclude, if provided
        self.__formatExcludePairs__()

    # Load data from unwrapStack.vrt
    def __loadStackData__(self):
        '''
        Load data from unwrapStack.vrt file.
        '''
        # Open dataset
        self.IFGs = osgeo.gdal.Open(self.imgfile)

        # Format extent
        N = self.IFGs.RasterXSize
        M = self.IFGs.RasterYSize
        tnsf = self.IFGs.GetGeoTransform()
        left = tnsf[0]
        xstep = tnsf[1]
        right = left + N * xstep
        top = tnsf[3]
        ystep = tnsf[5]
        bottom = top + M * ystep
        self.extent = (left, right, bottom, top)

        # Report if requested
        LOGGER.debug('%s bands detected', self.IFGs.RasterCount)

    # Convert date pair to string
    def __datePair2strPair__(self, datePair):
        '''
        Convert pair in format [master, slave] to date in format 'master_slave'
        '''
        masterStr = datePair[0].strftime('%Y%m%d')
        secondaryStr = datePair[1].strftime('%Y%m%d')
        strPair = '{}_{}'.format(masterStr, secondaryStr)
        return strPair

    # Convert pair string to date list
    def __strPair2datePair__(self, pair):
        '''
        Convert pair in format 'master_slave' to date in format [master, slave]
        '''
        pair = pair.split('_')
        masterDate = datetime.datetime.strptime(pair[0], '%Y%m%d')
        secondaryDate = datetime.datetime.strptime(pair[1], '%Y%m%d')
        datePair = [masterDate, secondaryDate]
        return datePair

    # Convert pair list to date lists
    def __pairList2dateList__(self, pairList):
        '''
        Convert list of pairs in format ['master_slave','master_slave',...]
        to dates in format [[master, slave], [master, slave], ...]
        '''
        pairDates = []
        for pair in pairList:
            datePair = self.__strPair2datePair__(pair)
            pairDates.append(datePair)
        return pairDates

    # Date pairs and unique dates
    def __formatDates__(self):
        '''
        Retrieve list of date pairs and unique dates (epochs).
        The "pairs" attribute is a formatted list of interferogram date
        pairs **in the order in which they are stored in the
        unwrapStack.vrt** file. This list should not be modified.
        Constrain the list of epochs available for triplet determination
        using the "startDate" and "endDate" provided in the __init__ function.
        '''
        # Pairs - list of pairs composing the data set, in the order they are
        # written
        pairs = [os.path.basename(fname) for fname in self.IFGs.GetFileList()]
        pairs = [pair.split('.')[0] for pair in pairs]  # remove extensions
        pairs.remove('unwrapStack')  # remove extra file name
        self.pairs = self.__pairList2dateList__(pairs)

        # Get unique dates from date pairs
        self.epochs = []
        [self.epochs.extend(pair) for pair in self.pairs]
        self.epochs = list(set(self.epochs))  # unique dates only
        self.epochs.sort()  # sort oldest-youngest

        self.nEpochs = len(self.epochs)
        LOGGER.debug('%s unique dates detected', self.nEpochs)

        # Limit dates available for triplet formulation by start and end date
        self.tripletEpochs = [
            epoch for epoch in self.epochs if epoch >= self.startDate]
        self.tripletEpochs = [
            epoch for epoch in self.epochs if epoch <= self.endDate]
        self.nTripletEpochs = len(self.tripletEpochs)

        # Update start and end dates
        self.startDate = self.tripletEpochs[0]
        self.endDate = self.tripletEpochs[-1]

        # Ticks for plots
        self.dateTicks = pd.date_range(
            self.startDate - timedelta(days=30),
            self.endDate + timedelta(days=30), freq='MS')
        self.dateTickLabels = [
            date.strftime('%Y-%m') for date in self.dateTicks]

    # Format dates in list to exclude
    def __formatExcludePairs__(self):
        '''
        Check that exclude dates are in one of two formats:
        1. a string containing the pairs in YOUNGER_OLDER format, space-separated
        2. a .txt file with lines of the same formatting
        Formatting should match "pair" formatting: [[master,slave]]
        '''
        if self.excludePairs is not None:

            # Determine whether list or text file
            if self.excludePairs[-4:] == '.txt':

                # Treat as text file with list
                with open(self.excludePairs, 'r') as exclFile:
                    excludePairs = exclFile.readlines()
                    excludePairs = [pair.strip('\n') for pair in excludePairs]
                    self.excludePairs = self.__pairList2dateList__(
                        excludePairs)

            else:
                # Treat as list - split at spaces
                excludePairs = self.excludePairs.split(' ')
                self.excludePairs = self.__pairList2dateList__(excludePairs)

        else:
            # Include as empty list
            self.excludePairs = []

    # Plot pairs
    def plotPairs(self):
        '''
        Plot the timespans of interferogram pairs.
        '''
        # Copy the list of pairs and sort them in time order
        pairs = self.pairs[:]  # copy to separate object
        pairs.sort(key=lambda s: s[1])  # sort by secondary date

        # Plot pairs in time
        pairFig = plt.figure()
        pairAx = pairFig.add_subplot(111)
        for n, pair in enumerate(pairs):

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
            pairAx.plot([pair[0], pair[1]], [n, n], color=color,
                        label=label, linestyle=linestyle)

        # Format x-axis
        pairAx.set_xticks(self.dateTicks)
        pairAx.set_xticklabels(self.dateTickLabels, rotation=90)

        # Legend
        handles, labels = pairAx.get_legend_handles_labels()
        uniqueLabels = dict(zip(labels, handles))
        pairAx.legend(uniqueLabels.values(), uniqueLabels.keys(),
                      bbox_to_anchor=(0.005, 0.99), loc='upper left', borderaxespad=0.)

        # Other formatting
        pairAx.set_yticks([])
        pairAx.set_title('IFG pairs')
        pairFig.tight_layout()

    # Create triplet list
    def createTriplets(self, minTime=None, maxTime=None, printTriplets=False):
        '''
        Create a list of triplets given the date list and user-specified
        parameters. First generate a list of all possible triplets based on
        the available dates. Then validate that list across the list of
        existing pairs.

        The stack object retains an ordered list of dates in both YYYYMMDD
        format and datetime format, and pairList, based on the
        __format_dates__ function.
        '''
        LOGGER.debug('Creating list of all possible triplets')

        # Loop through dates to create all possible triplet combinations
        self.triplets = []
        for i in range(self.nTripletEpochs - 2):
            for j in range(i + 1, self.nTripletEpochs - 1):
                for k in range(j + 1, self.nTripletEpochs):
                    epochI = self.tripletEpochs[i]  # first date in sequence
                    epochJ = self.tripletEpochs[j]  # second date in sequence
                    epochK = self.tripletEpochs[k]  # third date in sequence
                    self.triplets.append(
                        [[epochJ, epochI], [epochK, epochJ], [epochK, epochI]])

        # Remove triplets with pairs in "exclude pairs" list
        self.__checkExcludedTriplets__()

        # Check that pairs meet time requirements
        self.__checkTripletsMinTime__(minTime)

        self.__checkTripletsMaxTime__(maxTime)

        # Check against list of existing interferograms
        self.__checkTripletsExist__()

        # Finished sorting
        self.nTriplets = len(self.triplets)

        # Print to text file
        with open(os.path.join(self.workdir, 'ValidTriplets.txt'), 'w') as tripletFile:
            for triplet in self.triplets:
                strPair = [self.__datePair2strPair__(pair) for pair in triplet]
                tripletFile.write('{}\n'.format(strPair))
            tripletFile.close()

        # Report if requested
        if printTriplets == True:

            # Print to screen
            LOGGER.info('Existing triplets:')
            for triplet in self.triplets:
                LOGGER.info([
                    self.__datePair2strPair__(pair) for pair in triplet])

        if self.verbose == True:
            LOGGER.info(
                '%s existing triplets found based on search criteria',
                self.nTriplets)

        # Reference dates
        self.tripletDates = [[triplet[0][1], triplet[1]
                              [1], triplet[2][0]] for triplet in self.triplets]

    # Check triplets against excluded pairs
    def __checkExcludedTriplets__(self):
        '''
        Check triplet list against excluded pairs list. Remove the triplet
        if any of the pairs is listed in "exclude pairs".
        '''
        validTriplets = []
        for triplet in self.triplets:

            # If no pairs are excluded, append to the valid triplets list
            invalidTriplets = 0  # reset counter
            for pair in triplet:

                if pair in self.excludePairs:
                    invalidTriplets += 1

            if invalidTriplets == 0:
                validTriplets.append(triplet)

        self.triplets = validTriplets  # update triplets list

    # Check triplets against minTime
    def __checkTripletsMinTime__(self, minTime):
        '''
        Check that all pairs in a triplet are longer in duration than the
        minimum time interval specified.
        '''
        if minTime:
            validTriplets = []
            for triplet in self.triplets:

                # Determine intervals between dates in days
                intervals = [(pair[0] - pair[1]).days for pair in triplet]
                if min(intervals) >= minTime:
                    validTriplets.append(triplet)
            self.triplets = validTriplets  # update triplets list

    # Check triplets against maxTime
    def __checkTripletsMaxTime__(self, maxTime):
        '''
        Check that all pairs in a triplet are shorter in duration than
        the maximum time interval specified.
        '''
        if maxTime:
            validTriplets = []
            for triplet in self.triplets:

                # Determine intervals between dates in days
                intervals = [(pair[0] - pair[1]).days for pair in triplet]
                if max(intervals) <= maxTime:
                    validTriplets.append(triplet)

            self.triplets = validTriplets  # update triplets list

    # Check triplets exist
    def __checkTripletsExist__(self):
        '''
        Check list of all possible triplets against the list of pairs
        that actually exist.
        '''
        existingTriplets = []
        for triplet in self.triplets:

            existing = 0  # reset count of existing pairs
            # Check that each pair of the triplet has a corresponding
            # interferogram
            for tripletPair in triplet:
                if tripletPair in self.pairs:
                    existing += 1  # update if ifg exists

            if existing == 3:
                existingTriplets.append(triplet)

        self.triplets = existingTriplets  # update triplet list

    # Plot triplets

    def plotTriplets(self):
        '''
        Plot triplets.
        '''
        # Setup figure
        tripletFig = plt.figure()
        tripletAx = tripletFig.add_subplot(111)

        # Plot triplets
        for i in range(self.nTriplets):
            tripletAx.plot(self.triplets[i], [i, i, i], 'k', marker='o')

        # Format x-axis
        tripletAx.set_xticks(self.dateTicks)
        tripletAx.set_xticklabels(self.dateTickLabels, rotation=90)

        # Other formatting
        tripletAx.set_yticks([])
        tripletAx.set_title('Triplets in data set')
        tripletFig.tight_layout()

    # Compute misclosure
    # Geo to map coordinates
    def LoLa2XY(self, lon, lat):
        '''
        Convert lon/lat coordinates to XY.
        '''
        tnsf = self.IFGs.GetGeoTransform()
        x = (lon - tnsf[0]) / tnsf[1]
        y = (lat - tnsf[3]) / tnsf[5]
        return x, y

    # Map to geo coordinates
    def XY2LoLa(self, x, y):
        '''
        Convert X/Y to lon/lat.
        '''
        tnsf = self.IFGs.GetGeoTransform()
        lon = tnsf[0] + tnsf[1] * x
        lat = tnsf[3] + tnsf[5] * y
        return lon, lat

    # Reference point formatting
    def __referencePoint__(self, refXY, refLoLa):
        '''
        Determine the reference point in XY coordinates. The reference point can be
        automatically or manually selected by the user and is subtracted
        from each interferogram.
        The point can be given in pixels or lon/lat coordinates. If given in Lat/Lon, determine
        the location in XY.
        '''
        LOGGER.debug('Determining reference point...')

        if refLoLa.count(None) == 0:
            # Determine the XY coordinates from the given lon/lat
            self.refLon = refLoLa[0]
            self.refLat = refLoLa[1]
            x, y = self.LoLa2XY(refLoLa[0], refLoLa[1])
            self.refX = int(x)
            self.refY = int(y)
            LOGGER.debug(
                'Reference point given as: X %s / Y %s; Lon %s / Lat %s',
                self.refX, self.refY, self.refLon, self.refLat)

        elif refXY.count(None) == 0:
            # Use the provided XY coordinates
            self.refX = refXY[0]
            self.refY = refXY[1]
            self.refLon, self.refLat = self.XY2LoLa(refXY[0], refXY[1])
            LOGGER.debug(
                'Reference point given as: X %s / Y %s; Lon %.4f / Lat %.4f',
                self.refX, self.refY, self.refLon, self.refLat)

        else:
            # Use a random reference point
            self.__autoReferencePoint__()

    # Random reference point
    def __autoReferencePoint__(self):
        '''
        Use the coherence stack to automatically determine a suitable reference point.
        '''
        # Load coherence data from cohStack.vrt
        cohfile = os.path.join(self.imgdir, 'cohStack.vrt')
        cohDS = osgeo.gdal.Open(cohfile)
        cohMap = np.zeros((cohDS.RasterYSize, cohDS.RasterXSize))
        coh_min = 0.7

        for n in range(1, cohDS.RasterCount + 1):
            cohMap += cohDS.GetRasterBand(n).ReadAsArray()
        aveCoherence = cohMap / cohDS.RasterCount
        cohMask = (aveCoherence >= coh_min)

        # Start with initial guess for reference point
        self.refX = np.random.randint(cohDS.RasterXSize)
        self.refY = np.random.randint(cohDS.RasterYSize)

        # Loop until suitable reference point is found
        n = 0
        while cohMask[self.refY, self.refX] == False:

            # Reselect reference points
            self.refX = np.random.randint(cohDS.RasterXSize)
            self.refY = np.random.randint(cohDS.RasterYSize)
            n += 1  # update counter

            # Break loop after 10000 iterations
            if n == 10000:
                msg = f'No reference point with coherence >= {coh_min} found'
                LOGGER.error(msg)
                raise Exception(msg)

        # Convert to lon/lat
        self.refLon, self.refLat = self.XY2LoLa(self.refX, self.refY)

        LOGGER.debug(
            'Reference point chosen randomly as: X %s / Y %s; Lon %.4f / '
            'Lat %.4f.', self.refX, self.refY, self.refLon, self.refLat)

    # Compute misclosure
    def computeMisclosure(self, refXY=None, refLoLa=None):
        '''
        Compute the misclosure of the phase triplets.
        A common reference point is required because the ifgs are not
        coregistered.
        '''
        # Determine reference point
        self.__referencePoint__(refXY, refLoLa)

        # Misclosure placeholders
        self.netMscStack = []
        self.absMscStack = []

        # Compute phase triplets
        LOGGER.debug('Calculating misclosure')

        for triplet in self.triplets:
            # Triplet date pairs
            JIdates = triplet[0]
            KJdates = triplet[1]
            KIdates = triplet[2]

            # Triplet indices - add 1 because raster bands start at 1
            JIndx = self.pairs.index(JIdates) + 1
            KJndx = self.pairs.index(KJdates) + 1
            KIndx = self.pairs.index(KIdates) + 1

            # Interferograms
            JI = self.IFGs.GetRasterBand(JIndx).ReadAsArray()
            KJ = self.IFGs.GetRasterBand(KJndx).ReadAsArray()
            KI = self.IFGs.GetRasterBand(KIndx).ReadAsArray()

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

    # Plot and analyze misclosure
    # Plotting miscellaneous functions
    def __backgroundDetect__(self, img):
        '''
        Detect the background value of an image.
        '''
        edges = np.concatenate([
            img[:, 0].flatten(), img[0, :].flatten(), img[-1, :].flatten(),
            img[:, -1].flatten()])
        backgroundValue = scipy.stats.mode(edges)[0][0]
        return backgroundValue

    def __imgClipValues__(self, img, percentiles):
        '''
        Find values at which to clip the images (min/max) based on
        histogram percentiles.
        '''
        clipValues = {}
        clipValues['min'], clipValues['max'] = np.percentile(
            img.flatten(), percentiles)
        return clipValues

    def __plotCumNetMisclosure__(self):
        '''
        Plot cumulative misclosure.
        '''
        cax = self.netMscAx.imshow(
            self.cumNetMsc, cmap='plasma', vmin=self.cumNetMscClips['min'],
            vmax=self.cumNetMscClips['max'], zorder=1)
        self.netMscAx.plot(self.refX, self.refY, 'ks', zorder=2)
        self.netMscAx.set_xticks([])
        self.netMscAx.set_yticks([])
        self.netMscAx.set_title('Cumulative misclosure')
        return cax

    def __plotCumAbsMisclosure__(self):
        '''
        Plot cumulative absolute misclosure.
        '''
        cax = self.absMscAx.imshow(
            self.cumAbsMsc, cmap='plasma', vmin=self.cumAbsMscClips['min'],
            vmax=self.cumAbsMscClips['max'], zorder=1)
        self.absMscAx.plot(self.refX, self.refY, 'ks', zorder=2)
        self.absMscAx.set_xticks([])
        self.absMscAx.set_yticks([])
        self.absMscAx.set_title('Cumulative absolute misclosure')
        return cax

    def __plotSeries__(self, ax, data, title):
        '''
        Plot misclosure timeseries.
        '''
        # Plot data
        if self.plotTimeIntervals == False:
            ax.plot([tripletDate[1]
                    for tripletDate in self.tripletDates], data, '-k.')
        else:
            for n in range(self.nTriplets):
                ax.plot([self.tripletDates[n][0], self.tripletDates[n][2]],
                        [data[n], data[n]], 'k')
                ax.plot(self.tripletDates[n][1], data[n], 'ko')

        # Formatting
        ax.set_xticks(self.dateTicks)
        ax.set_xticklabels([])
        ax.set_ylabel(title)

    # Plot misclosure
    def plotCumMisclosure(self, queryXY=None, queryLoLa=None,
                          pctmin=1, pctmax=99, plotTimeIntervals=False):
        '''
        Map-view plot of cumulative misclosure.
        '''
        LOGGER.debug('Begin misclosure analysis')

        # Parameters
        self.plotTimeIntervals = plotTimeIntervals

        # Set up interactive plots
        self.netMscFig = plt.figure()
        self.netMscAx = self.netMscFig.add_subplot(111)

        self.absMscFig = plt.figure()
        self.absMscAx = self.absMscFig.add_subplot(111)

        # Auto-detect background and clip values
        self.cumNetMscBackground = self.__backgroundDetect__(
            self.cumNetMisclosure)
        self.cumNetMsc = np.ma.array(
            self.cumNetMisclosure,
            mask=(self.cumNetMisclosure == self.cumNetMscBackground))
        self.cumNetMscClips = self.__imgClipValues__(
            self.cumNetMsc, percentiles=[pctmin, pctmax])

        self.cumAbsMscBackground = self.__backgroundDetect__(
            self.cumAbsMisclosure)
        self.cumAbsMsc = np.ma.array(
            self.cumAbsMisclosure,
            mask=(self.cumAbsMisclosure == self.cumAbsMscBackground))
        self.cumAbsMscClips = self.__imgClipValues__(
            self.cumAbsMsc, percentiles=[pctmin, pctmax])

        # Plot maps
        cax = self.__plotCumNetMisclosure__()
        cbar = self.netMscFig.colorbar(cax, orientation='vertical')
        cbar.set_label('cum. misclosure (radians)')

        cax = self.__plotCumAbsMisclosure__()
        cbar = self.absMscFig.colorbar(cax, orientation='vertical')
        cbar.set_label('cum. abs. misclosure (radians)')

        # Plot timeseries points
        self.mscSeriesFig = plt.figure('Misclosure', figsize=(8, 8))
        self.netMscSeriesAx = self.mscSeriesFig.add_subplot(411)
        self.cumNetMscSeriesAx = self.mscSeriesFig.add_subplot(412)
        self.absMscSeriesAx = self.mscSeriesFig.add_subplot(413)
        self.cumAbsMscSeriesAx = self.mscSeriesFig.add_subplot(414)

        # Pre-specified query points
        self.__misclosureQuery__(queryXY, queryLoLa)

        # Link canvas to plots for interaction
        self.netMscFig.canvas.mpl_connect(
            'button_press_event', self.__misclosureAnalysis__)
        self.absMscFig.canvas.mpl_connect(
            'button_press_event', self.__misclosureAnalysis__)

    # Misclosure analysis
    def __misclosureAnalysis__(self, event):
        '''
        Show the time history of each pixel based on interactive map.
        '''
        px = event.xdata
        py = event.ydata
        px = int(round(px))
        py = int(round(py))

        # Report position and cumulative values
        LOGGER.info('px %s py %s', px, py)  # report position
        LOGGER.info('Cumulative misclosure: %s', self.cumNetMisclosure[py, px])
        LOGGER.info(
            'Abs cumulative misclosure: %s', self.cumAbsMisclosure[py, px])

        # Plot query points on maps
        self.netMscAx.cla()
        self.__plotCumNetMisclosure__()
        self.netMscAx.plot(
            px, py, color='k', marker='o', markerfacecolor='w', zorder=3)

        self.absMscAx.cla()
        self.__plotCumAbsMisclosure__()
        self.absMscAx.plot(
            px, py, color='k', marker='o', markerfacecolor='w', zorder=3)

        # Plot misclosure over time
        LOGGER.info('Misclosure: %s', self.netMscStack[:, py, px])
        self.netMscSeriesAx.cla()  # misclosure
        self.__plotSeries__(
            self.netMscSeriesAx, self.netMscStack[:, py, px], 'misclosure')

        self.cumNetMscSeriesAx.cla()  # cumulative misclosure
        self.__plotSeries__(self.cumNetMscSeriesAx, np.cumsum(
            self.netMscStack[:, py, px]), 'cum. miscl.')

        self.absMscSeriesAx.cla()  # absolute misclosure
        self.__plotSeries__(
            self.absMscSeriesAx, self.absMscStack[:, py, px], 'abs. miscl')

        self.cumAbsMscSeriesAx.cla()  # cumulative absolute misclosure
        self.__plotSeries__(
            self.cumAbsMscSeriesAx, np.cumsum(self.absMscStack[:, py, px]),
            'cum. abs. miscl.')

        # Format x-axis
        dates = pd.date_range(
            self.startDate - timedelta(days=30),
            self.endDate + timedelta(days=30), freq='MS')
        dateLabels = [date.strftime('%Y-%m') for date in dates]
        self.cumAbsMscSeriesAx.set_xticklabels(
            self.dateTickLabels, rotation=90)

        # Draw outcomes
        self.netMscFig.canvas.draw()
        self.absMscFig.canvas.draw()
        self.mscSeriesFig.canvas.draw()

    # Misclosure query
    def __misclosureQuery__(self, queryXY=None, queryLoLa=None):
        '''
        Show the time history of each pixel based on pre-specified selection.
        '''
        LOGGER.debug('Pre-specified query point...')

        # Convert bewteen lon/lat and image coordinates
        if queryLoLa.count(None) == 0:
            # Determine the XY coordinates from the given lon/lat
            qLon, qLat = queryLoLa
            qx, qy = self.LoLa2XY(queryLoLa[0], queryLoLa[1])
            qx = int(qx)
            qy = int(qy)  # convert to integer values

        elif queryXY.count(None) == 0:
            # Use the provided XY coordinates
            qx, qy = queryXY
            qLon, qLat = self.XY2LoLa(queryXY[0], queryXY[1])

        LOGGER.debug(
            'Query point: X %s / Y %s; Lon %.4f / Lat %.4f', qx, qy, qLon, qLat)

        # Plot query points on map
        self.netMscAx.plot(
            qx, qy, color='k', marker='o', markerfacecolor='w', zorder=3)
        self.absMscAx.plot(
            qx, qy, color='k', marker='o', markerfacecolor='w', zorder=3)

        # Plot misclosure over time
        self.__plotSeries__(
            self.netMscSeriesAx, self.netMscStack[:, qy, qx], 'misclosure')
        self.__plotSeries__(
            self.cumNetMscSeriesAx,
            np.cumsum(self.netMscStack[:, qy, qx]), 'cum. miscl.')
        self.__plotSeries__(
            self.absMscSeriesAx, self.absMscStack[:, qy, qx], 'abs. miscl')
        self.__plotSeries__(
            self.cumAbsMscSeriesAx, np.cumsum(self.absMscStack[:, qy, qx]),
            'cum. abs. miscl.')

        # Format x-axis
        dates = pd.date_range(
            self.startDate - timedelta(days=30),
            self.endDate + timedelta(days=30), freq='MS')
        dateLabels = [date.strftime('%Y-%m') for date in dates]
        self.cumAbsMscSeriesAx.set_xticklabels(
            self.dateTickLabels, rotation=90)

    # Plot triplet misclosure maps
    def plotTripletMaps(self, pctmin=1, pctmax=99):
        '''
        Plot the misclosure measurements for each triplet to figures.
        '''
        LOGGER.debug('Saving incremental misclosure maps to image files')

        # Parameters
        subplotDims = (2, 2)
        maxSubplots = subplotDims[0] * subplotDims[1]

        # Output directory/subdirectory
        self.figdir = os.path.join(self.workdir, 'MisclosureFigs')
        try:
            os.mkdir(self.figdir)

        except BaseException:
            pass

        figNb = 0  # start figure counter
        plotNb = 1  # start subplot counter
        for i in range(self.nTriplets):

            # Plot number and subplot position
            if plotNb % maxSubplots == 1:

                # Spawn new figure
                figNb += 1  # update figure counter
                plotNb = 1  # reset subplot counter
                Fig = plt.figure(figsize=(8, 6))
            ax = Fig.add_subplot(subplotDims[0], subplotDims[1], plotNb)

            # Format misclosure map
            mscMapBackground = self.__backgroundDetect__(
                self.netMscStack[i, :, :])
            mscMap = np.ma.array(self.netMscStack[i, :, :], mask=(
                self.netMscStack[i, :, :] == mscMapBackground))
            mscMapClips = self.__imgClipValues__(mscMap, [pctmin, pctmax])

            # Plot misclosure
            cax = ax.imshow(
                mscMap, cmap='plasma', vmin=mscMapClips['min'],
                vmax=mscMapClips['max'])

            # Format axis
            ax.set_xticks([])
            ax.set_yticks([])
            Fig.colorbar(cax, orientation='horizontal')
            dates = [date.strftime('%Y%m%d') for date in self.tripletDates[i]]
            ax.set_title('_'.join(dates))

            if plotNb == maxSubplots:
                # Save old figure
                Fig.suptitle('Phase misclosure (radians)')
                figname = 'MisclosureValues_fig{}.png'.format(figNb)
                figpath = os.path.join(self.figdir, figname)
                Fig.savefig(figpath, dpi=300)

            plotNb += 1

            # Save final figure
            Fig.suptitle('Phase misclosure (radians)')
            figname = 'MisclosureValues_fig{}.png'.format(figNb)
            figpath = os.path.join(self.figdir, figname)
            Fig.savefig(figpath, dpi=300)

    # Save cumulative misclosure plots to geotiffs

    def saveCumMisclosure(self):
        '''
        Save cumulative (/absolute) misclosure plots to georeferenced tiff
        files. Use metadata from unwrapStack.vrt file.
        '''
        LOGGER.debug('Saving misclosure maps to geotiffs')

        # Fix background values
        self.cumNetMisclosure[
            self.cumNetMisclosure == self.cumNetMscBackground] = 0
        self.cumAbsMisclosure[
            self.cumAbsMisclosure == self.cumAbsMscBackground] = 0

        # Save cumulative misclosure
        cumNetMscSavename = os.path.join(
            self.figdir, 'CumulativeMisclosure.tif')
        self.__saveGeoTiff__(cumNetMscSavename, self.cumNetMisclosure)

        # Save cumulative absolute misclosure
        cumAbsMscSavename = os.path.join(
            self.figdir, 'CumulativeAbsoluteMisclosure.tif')
        self.__saveGeoTiff__(cumAbsMscSavename, self.cumAbsMisclosure)

    def __saveGeoTiff__(self, savename, img):
        '''
        Template for saving geotiffs.
        '''
        driver = osgeo.gdal.GetDriverByName('GTiff')
        DSout = driver.Create(
            savename, self.IFGs.RasterXSize, self.IFGs.RasterYSize, 1,
            osgeo.gdal.GDT_Float32)
        DSout.GetRasterBand(1).WriteArray(img)
        DSout.GetRasterBand(1).SetNoDataValue(0)
        DSout.SetProjection(self.IFGs.GetProjection())
        DSout.SetGeoTransform(self.IFGs.GetGeoTransform())
        DSout.FlushCache()
