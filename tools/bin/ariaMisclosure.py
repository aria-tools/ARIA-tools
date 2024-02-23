#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Compute and analyze the phase triplet misclosure for a series of
# interferograms.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import functions
import argparse
import logging

import ARIAtools.computeMisclosure
import ARIAtools.util.logger

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from pkg_resources import get_distribution

LOGGER = logging.getLogger('ariaMisclosure.py')

DESCRIPTION = '''
Compute the cumulative misclosure of phase triplets based on a set of
interferograms saved in the stack/unwrapStack.vrt data set. During triplet
computation, values at a reference point are removed from the interferograms
prior to misclosure computation to account for abiguities in the
unwrapped phase that might arise during pairwise computation.

The code works by generating a list of triplets from the available pairs.
From the list of triplets, the routine will:

    1. Formulate a list of viable triplet combinations (IJ, JK, IK)
    2. Subtract the reference point from those interferograms
    3. Compute the disagreement based on IJ+JK-IK
       ... continue through the list of triplets.

Both the misclosure values (positive or negative) and the absolute misclosure
values (always positive) are computed and stored in 3D arrays, where each
slice represents a triplet. The "cumulative" misclosure and absolut
misclosure are computed by summing the 3D arrays.

Once the (absolute) misclosure is calculated, the user can view the time
history of any pixel by clicking on the maps.

Thumbnail images of the misclosure associated with any given triplet are
saved in the MisclosureFigs folder. Additionally, georeferenced tiffs of the cumulative misclosure and absolute cumulative
misclosure maps are saved.
'''

EXAMPLES = '''EXAMPLES

# Using the unwrapStack.vrt to call all interferograms, automatically find a
reference point
ariaMisclosure.py -f stack/unwrapStack.vrt

# Provide a predefined reference point
ariaMisclosure.py -f stack/unwrapStack.vrt -refLon 89.358 -refLat 32.621

# Limit triplet selection by time interval (12 days to 48 days)
ariaMisclosure.py -f stack/unwrapStack.vrt --mintime 12 --maxtime 48
'''

def create_parser():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter,
        epilog=EXAMPLES)

    # Input data
    parser.add_argument(
        '-f', '--file', dest='imgfile', type=str, required=True,
        help='ARIA files. Specify the stack/unwrapStack.vrt file, or a '
             'wildcard operator in the unwrappedPhase folder (see EXAMPLES)')
    parser.add_argument(
        '-w', '--workdir', dest='workdir', type=str, default='./',
        help='Specify directory to deposit all outputs. Default is local '
             'directory where script is launched.')
    parser.add_argument(
        '--startdate', dest='startDate', type=str, default='20140615',
        help='Start date for data series')
    parser.add_argument(
        '--enddate', dest='endDate', type=str, default=None,
        help='End date for data series')
    parser.add_argument(
        '--exclude-pairs', dest='excludePairs', type=str, default=None,
        help='List of pairs to exclude, e.g., "20160116_20160101 '
             '20171031_20161030". This can also be provided in as a text '
             'file with one line per date pair.')
    parser.add_argument(
        '--plot-pairs', dest='plotPairs', action='store_true',
        help='Plot the timespans of date pairs')

    # Triplet formulation
    parser.add_argument(
        '--mintime', dest='minTime', type=int, default=None,
        help='Minimum time span of pairs in triplets (days)')
    parser.add_argument(
        '--maxtime', dest='maxTime', type=int, default=None,
        help='Maximum time span of pairs in triplets (days)')
    parser.add_argument(
        '--print-triplets', dest='printTriplets', action='store_true',
        help='Print list of existing triplets (i.e., those included in the '
             'data set).')
    parser.add_argument(
        '--plot-triplets', dest='plotTriplets', action='store_true',
        help='Plot existing triplets')

    # Reference point
    parser.add_argument(
        '-refX', dest='refX', type=int, default=None, help='Reference X pixel')
    parser.add_argument(
        '-refY', dest='refY', type=int, default=None, help='Reference Y pixel')
    parser.add_argument(
        '-refLon', dest='refLon', type=float, default=None,
        help='Reference longitude')
    parser.add_argument(
        '-refLat', dest='refLat', type=float, default=None,
        help='Reference latitude')

    # Query point
    parser.add_argument(
        '--queryX', dest='queryX', type=int, default=None,
        help='Query point X pixel')
    parser.add_argument(
        '--queryY', dest='queryY', type=int, default=None,
        help='Query point Y pixel')
    parser.add_argument(
        '--queryLon', dest='queryLon', type=float, default=None,
        help='Query point longitude')
    parser.add_argument(
        '--queryLat', dest='queryLat', type=float, default=None,
        help='Query point latitude')

    # Vocalization
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='Verbose mode')

    # Misclosure map formatting
    parser.add_argument(
        '--pctmin', dest='pctMinClip', type=float, default=1,
        help='Minimum percent clip value for cumulative misclosure plot')
    parser.add_argument(
        '--pctmax', dest='pctMaxClip', type=float, default=99,
        help='Maximum percent clip value for cumulative misclosure plot')
    parser.add_argument(
        '--plot-time-intervals', dest='plotTimeIntervals', action='store_true',
        help='Plot triplet intervals in misclosure analysis figure.')
    return parser


def main(inps=None):
    try:
        print('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except BaseException:
        pass

    # Gather arguments
    parser = createParser()
    args = parser.parse_args()

    if args.verbose:
        ARIAtools.util.logger.logger.setLevel(logging.DEBUG)

    # Load data based on data type
    dataStack = ARIATools.stack.Stack(
        imgfile=args.imgfile, workdir=args.workdir, startDate=args.startDate,
        endDate=args.endDate, excludePairs=args.excludePairs,
        verbose=args.verbose)

    # Plot pairs if requested
    if args.plotPairs == True:
        dataStack.plotPairs()

    # Create list of triplets
    dataStack.createTriplets(
        minTime=args.minTime, maxTime=args.maxTime,
        printTriplets=args.printTriplets)

    # Plot triplets if requested
    if args.plotTriplets == True:
        dataStack.plotTriplets()

    # Compute misclosure
    dataStack.computeMisclosure(
        refXY=[args.refX, args.refY], refLoLa=[args.refLon, args.refLat])

    # Plot and analyze data
    dataStack.plotCumMisclosure(
        queryXY=[args.queryX, args.queryY],
        queryLoLa=[args.queryLon, args.queryLat], pctmin=args.pctMinClip,
        pctmax=args.pctMaxClip,
        plotTimeIntervals=args.plotTimeIntervals)
    plt.show()

    # Save misclosure map for each triplet to figures
    dataStack.plotTripletMaps(pctmin=args.pctMinClip, pctmax=args.pctMaxClip)

    # Save misclosure maps to geotiffs
    dataStack.saveCumMisclosure()


if __name__ == '__main__':
    main()
