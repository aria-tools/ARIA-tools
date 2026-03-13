#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright (c) 2026, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Command line tool to assist with ordering ARIA Sentinel-1 GUNW products
from ASF. This tool helps users identify which ARIA frames intersect their
area of interest before placing an order, and generate interferogram network
configurations for ordering.

Note: This tool supports Sentinel-1 GUNW products only, not NISAR.
"""

import argparse
import csv
import datetime
import json
import logging
import os

import matplotlib.dates as mdates
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import requests
from shapely.geometry import shape, Polygon

try:
    import asf_search
    HAS_ASF_SEARCH = True
except ImportError:
    HAS_ASF_SEARCH = False

try:
    from asf_enumeration import aria_s1_gunw
    HAS_ASF_ENUMERATION = True
except ImportError:
    HAS_ASF_ENUMERATION = False

try:
    import hyp3_sdk
    HAS_HYP3_SDK = True
except ImportError:
    HAS_HYP3_SDK = False

try:
    import ARIAtools.util.log
    _LOG_FORMAT = ARIAtools.util.log.FORMAT
except ImportError:
    _LOG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

LOGGER = logging.getLogger('ariaOrderASF.py')

# URL for ARIA frames geojson from ASF enumeration repository
ARIA_FRAMES_URL = (
    'https://raw.githubusercontent.com/ASFHyP3/asf-enumeration/develop/'
    'src/asf_enumeration/frame_maps/aria_frames.geojson'
)

# Network generation types
NETWORK_TYPES = ['sequential', 'seasonal', 'annual']

# HyP3 credit costs
CREDITS_PER_PAIR = 60
ASF_MONTHLY_QUOTA = 8000


def create_parser():
    """
    Create argument parser for ariaOrderASF tool.

    Three mutually exclusive modes:
        --getframes : Discover ARIA frames covering an area of interest.
        --getpairs  : Generate interferogram pair lists for a frame.
        --orderpairs: Submit pairs to HyP3 for processing.

    Note: Only Sentinel-1 GUNW products are supported.
    """
    parser = argparse.ArgumentParser(
        description='Command line interface to identify ARIA Sentinel-1 GUNW '
                    'frames, generate interferogram pair lists, and order '
                    'products from ASF via HyP3.\n\n'
                    'NOTE: This tool supports Sentinel-1 GUNW only, '
                    'not NISAR.',
        epilog='Examples of use:\n\n'
               '\t # Step 1: Find frames intersecting a bounding box\n'
               '\t ariaOrderASF.py --getframes -b "36.0 37.0 -118.0 '
               '-117.0"\n\n'
               '\t # Narrow by track (returns both asc/desc for that '
               'track)\n'
               '\t ariaOrderASF.py --getframes -b "36.0 37.0 -118.0 '
               '-117.0" -t 064\n\n'
               '\t # Step 2a: Generate sequential nearest-3-neighbor '
               'pairs\n'
               '\t #   (identifies which pairs already exist at ASF)\n'
               '\t ariaOrderASF.py --getpairs --frame 25050 '
               '--num-neighbors 3 \\\n'
               '\t     -s 20230101 -e 20230601\n\n'
               '\t # Step 2b: Generate annual pairs over multiple '
               'years\n'
               '\t ariaOrderASF.py --getpairs --frame 22440 '
               '--network annual \\\n'
               '\t     -s 20200101 -e 20240101\n\n'
               '\t # Step 2c: Build both sequential + annual for the '
               'same frame\n'
               '\t #   Run each network type separately into its own '
               'work directory:\n'
               '\t ariaOrderASF.py --getpairs --frame 25050 '
               '--num-neighbors 1 \\\n'
               '\t     -s 20220101 -e 20230601 -w ./seq\n'
               '\t ariaOrderASF.py --getpairs --frame 25050 '
               '--network annual \\\n'
               '\t     -s 20220101 -e 20230601 -w ./annual\n\n'
               '\t # Seasonal pairs (same season across years)\n'
               '\t ariaOrderASF.py --getpairs --frame 25050 '
               '--network seasonal \\\n'
               '\t     --seasonal-window 30\n\n'
               '\t # Step 3a: Preview order (dry run — shows what would be\n'
               '\t #   submitted and the credit cost, but spends nothing)\n'
               '\t ariaOrderASF.py --orderpairs --frame 25050 '
               '--pairs-file ./aria_pairs_frame25050.csv --dry-run\n\n'
               '\t # Step 3b: Submit only the first N pairs (safety cap).\n'
               '\t #   Useful for testing with a single pair before\n'
               '\t #   committing to a larger batch.\n'
               '\t ariaOrderASF.py --orderpairs --frame 25050 '
               '--pairs-file ./aria_pairs_frame25050.csv --limit 1\n\n'
               '\t # Step 3c: Submit all new pairs (prompts for\n'
               '\t #   confirmation with credit cost before submitting).\n'
               '\t #   The program will stop if the cost exceeds your\n'
               '\t #   remaining credits and suggest a --limit value.\n'
               '\t ariaOrderASF.py --orderpairs --frame 25050 '
               '--pairs-file ./aria_pairs_frame25050.csv\n\n'
               '\t # Step 4a: Check job status by the job name shown\n'
               '\t #   after submission\n'
               '\t ariaOrderASF.py --statusjobs '
               '--status-name ARIA_frame25050_20260311_120000\n\n'
               '\t # Step 4b: Check status of specific job IDs\n'
               '\t ariaOrderASF.py --statusjobs '
               '--job-ids "abc123-def456,ghi789-jkl012"\n',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # =========================================================================
    # Mode selection (mutually exclusive)
    # =========================================================================
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        '--getframes', action='store_true',
        help='Mode 1: Discover ARIA frames covering the bounding box.')
    mode_group.add_argument(
        '--getpairs', action='store_true',
        help='Mode 2: Generate interferogram pair list for a frame.')
    mode_group.add_argument(
        '--orderpairs', action='store_true',
        help='Mode 3: Submit pairs to HyP3 for processing.')
    mode_group.add_argument(
        '--statusjobs', action='store_true',
        help='Mode 4: Check status of previously submitted HyP3 jobs.')

    # =========================================================================
    # Frame discovery arguments (--getframes)
    # =========================================================================
    frame_group = parser.add_argument_group(
        'Frame Discovery (--getframes)',
        'Identify which ARIA frames intersect your area of interest.')
    frame_group.add_argument(
        '-b', '--bbox', default=None, type=str,
        help='Lat/Lon Bounding box SNWE (e.g., "36.0 37.0 -118.0 -117.0"), '
             'or GDAL-readable file containing POLYGON geometry. '
             'Required for --getframes.')
    frame_group.add_argument(
        '-t', '--track', default=None, type=str,
        help='Optional filter by track number; single number (e.g., 064) '
             'or comma-separated list (e.g., 064,137). '
             'Returns both ascending and descending unless -d is given.')
    frame_group.add_argument(
        '-d', '--direction', dest='flightdir', default=None, type=str,
        choices=['ascending', 'descending', 'a', 'd', 'A', 'D',
                 'ASCENDING', 'DESCENDING'],
        help='Optional filter by flight direction: ascending (a) or '
             'descending (d).')

    # =========================================================================
    # Pair generation arguments (--getpairs)
    # =========================================================================
    pair_group = parser.add_argument_group(
        'Pair Generation (--getpairs)',
        'Generate interferogram pairs for a single ARIA frame.')
    pair_group.add_argument(
        '--frame', default=None, type=int,
        help='ARIA frame ID (e.g., 7423). Required for --getpairs and '
             '--orderpairs.')
    pair_group.add_argument(
        '--network', default='sequential', type=str.lower,
        choices=NETWORK_TYPES,
        help='Type of interferogram network: "sequential" (nearest N '
             'neighbors), "seasonal" (same season across years), "annual" '
             '(nearest neighbor ~1 year apart). Default: sequential.')
    pair_group.add_argument(
        '--num-neighbors', dest='num_neighbors', default=1, type=int,
        help='Number of nearest temporal neighbors for sequential/annual '
             'network. Default: 1.')
    pair_group.add_argument(
        '--seasonal-window', dest='seasonal_window', default=30, type=int,
        help='Window in days (+/-) around same day-of-year for seasonal '
             'network. Default: 30.')
    pair_group.add_argument(
        '-s', '--start', default='20140101', type=str,
        help='Start date YYYYMMDD for acquisition search. '
             'Default: 20140101 (Sentinel-1 launch).')
    pair_group.add_argument(
        '-e', '--end', default=None, type=str,
        help='End date YYYYMMDD for acquisition search. Default: today.')

    # =========================================================================
    # Order arguments (--orderpairs)
    # =========================================================================
    order_group = parser.add_argument_group(
        'Ordering (--orderpairs)',
        'Submit interferogram pairs to HyP3 for processing.')
    order_group.add_argument(
        '--pairs-file', dest='pairs_file', default=None, type=str,
        help='Path to pairs CSV file (output of --getpairs). '
             'Required for --orderpairs. Users may edit this file to '
             'remove unwanted pairs before ordering.')
    order_group.add_argument(
        '--job-name', dest='job_name', default=None, type=str,
        help='HyP3 job name for tracking. If not provided, an '
             'auto-generated name will be used: '
             'ARIA_frame{ID}_{YYYYMMDD_HHMMSS}.')
    order_group.add_argument(
        '--dry-run', dest='dry_run', action='store_true',
        help='Show what would be submitted without actually submitting. '
             'Use this to preview the order before committing credits.')
    order_group.add_argument(
        '--limit', default=None, type=int,
        help='Maximum number of pairs to submit in this invocation. '
             'Acts as a safety cap to avoid accidentally submitting '
             'large batches. Pairs are submitted in CSV row order.')

    # =========================================================================
    # Status arguments (--statusjobs)
    # =========================================================================
    status_group = parser.add_argument_group(
        'Job Status (--statusjobs)',
        'Check the status of previously submitted HyP3 jobs.')
    status_group.add_argument(
        '--status-name', dest='status_name', default=None, type=str,
        help='HyP3 job name to query (the name used during --orderpairs). '
             'Returns all jobs matching this name.')
    status_group.add_argument(
        '--job-ids', dest='job_ids', default=None, type=str,
        help='Comma-separated list of HyP3 job IDs to query '
             '(e.g., "abc123,def456").')

    # =========================================================================
    # Common arguments
    # =========================================================================
    common_group = parser.add_argument_group('Common Options')
    common_group.add_argument(
        '-w', '--workdir', dest='wd', default='./', type=str,
        help='Output directory for files and plots. Default: ./')
    common_group.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable verbose output.')
    common_group.add_argument(
        '--log-level', default='info',
        choices=['debug', 'info', 'warning', 'error'],
        help='Logger log level. Default: info.')

    return parser


# =============================================================================
# Utility helpers
# =============================================================================

# Property name fallback maps for the ARIA frames GeoJSON
# The official ASF enumeration GeoJSON uses: id, path, dir
_FRAME_ID_KEYS = ('frame_id', 'frameID', 'id', 'frame')
_TRACK_KEYS = ('path', 'track', 'track_number', 'relativeOrbit',
               'relative_orbit')
_DIRECTION_KEYS = ('dir', 'direction', 'flightDirection',
                   'flight_direction', 'orbit_direction')


def get_frame_property(feature, keys, default='N/A'):
    """
    Extract a property from a frame feature, trying multiple key names.

    Parameters
    ----------
    feature : dict
        GeoJSON feature.
    keys : tuple
        Property names to try in order.
    default : object
        Value to return if no key matches.

    Returns
    -------
    object
        The first matching property value, or *default*.
    """
    props = feature.get('properties', {})
    for key in keys:
        val = props.get(key)
        if val is not None:
            return val
    return default


def make_bbox(inp_bbox):
    """
    Create a shapely Polygon from a SNWE bounding box string or shapefile.

    Parameters
    ----------
    inp_bbox : str
        Either a space-separated string of "S N W E" coordinates,
        or a path to a GDAL-readable file with polygon geometry.

    Returns
    -------
    shapely.geometry.Polygon
        Polygon representing the bounding box or loaded geometry.
    """
    if inp_bbox is None:
        return None

    if os.path.exists(os.path.abspath(inp_bbox)):
        from ARIAtools.util.shp import open_shp
        ring = open_shp(inp_bbox, 0, 0).exterior
        poly = Polygon(ring)
    else:
        try:
            S, N, W, E = [float(i) for i in inp_bbox.split()]

            # Adjust for degrees easting / northing (0 - 360 / 0:180)
            if W > 180:
                W -= 360
                LOGGER.info('Adjusted W from 0-360 to -180-180 range')

            if E > 180:
                E -= 360
                LOGGER.info('Adjusted E from 0-360 to -180-180 range')

            if N > 90:
                N -= 90
                S -= 90
                LOGGER.info('Adjusted N/S from 0-180 to -90-90 range')

            # Validate bounds
            if S >= N:
                raise ValueError(
                    f'South ({S}) must be less than North ({N})')
            if W >= E:
                raise ValueError(
                    f'West ({W}) must be less than East ({E})')

            poly = Polygon([(W, N), (W, S), (E, S), (E, N)])

        except ValueError as e:
            raise ValueError(
                f'Cannot parse --bbox argument: {e}. '
                'Input format should be "S N W E" (space-separated) '
                'or a valid shapefile path.'
            ) from e
        except Exception as e:
            raise Exception(
                'Cannot understand the --bbox argument. Input string was '
                'entered incorrectly or path does not exist.'
            ) from e

    return poly


def _ariaframe_to_feature(af):
    """Convert an ``asf_enumeration.aria_s1_gunw.AriaFrame`` to a
    GeoJSON-like feature dict consumable by the display helpers."""
    from shapely.geometry import mapping
    return {
        'type': 'Feature',
        'geometry': mapping(af.polygon),
        'properties': {
            'id': af.id,
            'path': af.path,
            'dir': af.flight_direction,
        },
    }


def fetch_aria_frames(url=ARIA_FRAMES_URL):
    """
    Fetch ARIA frames geojson from the ASF enumeration repository.

    Parameters
    ----------
    url : str, optional
        URL to the ARIA frames geojson file.

    Returns
    -------
    dict
        GeoJSON FeatureCollection containing all ARIA frames.
    """
    LOGGER.info('Fetching ARIA frames from ASF enumeration repository...')
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        geojson_data = response.json()
        LOGGER.info(
            'Successfully fetched %d ARIA frames.',
            len(geojson_data.get('features', []))
        )
        return geojson_data
    except requests.exceptions.RequestException as e:
        raise RuntimeError(
            f'Failed to fetch ARIA frames from {url}: {e}'
        ) from e


def filter_frames_by_bbox(frames_geojson, bbox_poly):
    """
    Filter ARIA frames by intersection with a bounding box polygon.

    Parameters
    ----------
    frames_geojson : dict
        GeoJSON FeatureCollection of ARIA frames.
    bbox_poly : shapely.geometry.Polygon
        Bounding box polygon to filter by.

    Returns
    -------
    list
        List of frame features that intersect the bounding box.
    """
    matching_frames = []

    for feature in frames_geojson.get('features', []):
        try:
            frame_geom = shape(feature['geometry'])
            if frame_geom.intersects(bbox_poly):
                matching_frames.append(feature)
        except Exception as e:
            LOGGER.warning(
                'Could not process frame geometry: %s', e
            )
            continue

    return matching_frames


def filter_frames_by_ids(frames_geojson, frame_ids):
    """
    Filter ARIA frames by specific frame IDs.

    Parameters
    ----------
    frames_geojson : dict
        GeoJSON FeatureCollection of ARIA frames.
    frame_ids : list
        List of frame IDs to select.

    Returns
    -------
    list
        List of frame features matching the specified IDs.
    """
    matching_frames = []
    frame_ids_set = set(frame_ids)

    for feature in frames_geojson.get('features', []):
        frame_id = get_frame_property(feature, _FRAME_ID_KEYS, default=None)
        if frame_id is not None and int(frame_id) in frame_ids_set:
            matching_frames.append(feature)

    return matching_frames


def filter_frames_by_track(frames, tracks):
    """
    Filter frames by track number(s).

    Parameters
    ----------
    frames : list
        List of frame features to filter.
    tracks : list
        List of track numbers to keep.

    Returns
    -------
    list
        Filtered list of frames matching the specified tracks.
    """
    if not tracks:
        return frames

    tracks_set = set(tracks)
    filtered = []

    for feature in frames:
        track = get_frame_property(feature, _TRACK_KEYS, default=None)
        if track is not None and int(track) in tracks_set:
            filtered.append(feature)

    return filtered


def filter_frames_by_direction(frames, direction):
    """
    Filter frames by flight direction (ascending/descending).

    Parameters
    ----------
    frames : list
        List of frame features to filter.
    direction : str
        Flight direction: 'ascending' or 'descending'.

    Returns
    -------
    list
        Filtered list of frames matching the specified direction.
    """
    if not direction:
        return frames

    # Normalize direction input
    direction_normalized = (
        'ascending' if direction.lower().startswith('a') else 'descending'
    )

    filtered = []

    for feature in frames:
        frame_dir = get_frame_property(feature, _DIRECTION_KEYS, default=None)
        if frame_dir and frame_dir.lower() == direction_normalized:
            filtered.append(feature)

    return filtered


# =============================================================================
# Network Generation Functions
# =============================================================================

def get_acquisitions_for_frame(frame_id, start_date, end_date):
    """
    Get all Sentinel-1 acquisition dates for an ARIA frame within a
    date range using ``asf_enumeration.aria_s1_gunw``.

    Parameters
    ----------
    frame_id : int
        ARIA frame ID.
    start_date : datetime.date
        Start date for acquisition search.
    end_date : datetime.date
        End date for acquisition search.

    Returns
    -------
    list
        Sorted list of ``datetime.date`` objects.
    """
    if not HAS_ASF_ENUMERATION:
        raise ImportError(
            'asf_enumeration is required for acquisition queries.\n'
            'Install with:  pip install "asf_search[asf-enumeration]"')

    LOGGER.info('Querying acquisitions for frame %s via '
                'asf_enumeration ...', frame_id)
    acquisitions = aria_s1_gunw.get_acquisitions(int(frame_id))
    dates = sorted(
        acq.date for acq in acquisitions
        if start_date <= acq.date <= end_date
    )
    LOGGER.info('Found %d acquisitions in %s – %s.',
                len(dates), start_date, end_date)
    return dates


def generate_sequential_pairs(dates, num_neighbors=3):
    """
    Generate sequential (nearest N neighbor) interferogram pairs.

    Parameters
    ----------
    dates : list
        List of acquisition dates (sorted).
    num_neighbors : int
        Number of nearest temporal neighbors to connect.

    Returns
    -------
    list
        List of tuples (reference_date, secondary_date) for each pair.
    """
    pairs = []
    dates = sorted(dates)

    for i, ref_date in enumerate(dates):
        # Connect to the next N acquisitions
        for j in range(1, num_neighbors + 1):
            if i + j < len(dates):
                sec_date = dates[i + j]
                pairs.append((ref_date, sec_date))

    LOGGER.info('Generated %d sequential pairs (num_neighbors=%d).',
                len(pairs), num_neighbors)
    return pairs


def generate_seasonal_pairs(dates, seasonal_window=30):
    """
    Generate seasonal interferogram pairs (same season across years).

    Connects acquisitions from the same day-of-year (+/- window) across
    different years.

    Parameters
    ----------
    dates : list
        List of acquisition dates (sorted).
    seasonal_window : int
        Window in days around same day-of-year to consider as "same season".

    Returns
    -------
    list
        List of tuples (reference_date, secondary_date) for each pair.
    """
    pairs = []
    dates = sorted(dates)

    for i, ref_date in enumerate(dates):
        ref_doy = ref_date.timetuple().tm_yday

        for j in range(i + 1, len(dates)):
            sec_date = dates[j]
            sec_doy = sec_date.timetuple().tm_yday

            # Check if within seasonal window (handle year boundary)
            doy_diff = min(
                abs(ref_doy - sec_doy),
                365 - abs(ref_doy - sec_doy)
            )

            # Only pair if in different years and within seasonal window
            if ref_date.year != sec_date.year and doy_diff <= seasonal_window:
                pairs.append((ref_date, sec_date))

    LOGGER.info('Generated %d seasonal pairs (window=%d days).',
                len(pairs), seasonal_window)
    return pairs


def generate_annual_pairs(dates, num_neighbors=1):
    """
    Generate annual nearest neighbor interferogram pairs.

    Connects each acquisition to nearest acquisition(s) approximately one
    year later.

    Parameters
    ----------
    dates : list
        List of acquisition dates (sorted).
    num_neighbors : int
        Number of nearest annual neighbors to connect.

    Returns
    -------
    list
        List of tuples (reference_date, secondary_date) for each pair.
    """
    pairs = []
    dates = sorted(dates)

    for ref_date in dates:
        # Target date approximately one year later
        target_date = ref_date + datetime.timedelta(days=365)

        # Find nearest dates to the target
        candidates = []
        for sec_date in dates:
            if sec_date > ref_date:
                days_from_target = abs((sec_date - target_date).days)
                # Consider candidates within 60 days of annual target
                if days_from_target <= 60:
                    candidates.append((days_from_target, sec_date))

        # Sort by proximity to annual target and take num_neighbors
        candidates.sort(key=lambda x: x[0])
        for _, sec_date in candidates[:num_neighbors]:
            pairs.append((ref_date, sec_date))

    LOGGER.info('Generated %d annual pairs (num_neighbors=%d).',
                len(pairs), num_neighbors)
    return pairs


def filter_pairs_by_temporal_baseline(pairs, max_baseline_days):
    """
    Filter interferogram pairs by maximum temporal baseline.

    Parameters
    ----------
    pairs : list
        List of (reference_date, secondary_date) tuples.
    max_baseline_days : int
        Maximum allowed temporal baseline in days.

    Returns
    -------
    list
        Filtered list of pairs.
    """
    if max_baseline_days is None:
        return pairs

    filtered = []
    for ref_date, sec_date in pairs:
        baseline = abs((sec_date - ref_date).days)
        if baseline <= max_baseline_days:
            filtered.append((ref_date, sec_date))

    LOGGER.info('Filtered to %d pairs with temporal baseline <= %d days.',
                len(filtered), max_baseline_days)
    return filtered


def check_existing_products(pairs, frame_id, max_results=5000):
    """
    Check which interferogram pairs already exist at ASF.

    Uses a single bulk ``asf_search.search`` query to fetch all existing
    GUNW products for the frame within the relevant date range, then
    matches pairs locally.  This is dramatically faster than calling
    ``product_exists`` per pair (1 HTTP request vs N).

    The CMR backend paginates at 250 results per HTTP request.  To avoid
    unbounded queries (and potential timeouts on very large archives),
    a ``max_results`` safety cap is applied.  If the cap is reached, we
    log a warning and fall back to narrower chunked queries (250 products
    at a time — the CMR default page size) for the unmatched pairs.

    Parameters
    ----------
    pairs : list
        List of ``(reference_date, secondary_date)`` tuples.
        Our convention: ``(earlier_date, later_date)``.
    frame_id : int
        ARIA frame ID.
    max_results : int, optional
        Upper bound on the number of products returned by the bulk query.
        Prevents runaway pagination for very large archives.  Default 5000
        (well above the ~1000 products expected for a single frame's full
        history).  ``asf_search`` pages at 250 results / request, so
        5000 corresponds to at most 20 HTTP round-trips.

    Returns
    -------
    dict
        ``{'existing': [...], 'new': [...]}``
    """
    existing = []
    new = []

    if not HAS_ASF_ENUMERATION:
        LOGGER.warning('asf_enumeration not available — cannot check '
                       'which products already exist; marking all as new.')
        return {'existing': [], 'new': list(pairs)}

    if not pairs:
        return {'existing': [], 'new': []}

    # ------------------------------------------------------------------
    # Bulk query: fetch existing GUNW products for this frame in one
    # request (paginated internally by asf_search at 250/page).
    # The ASF `start`/`end` parameters filter on the *reference*
    # (= later) date, so we derive the range from the later date in
    # each pair.  A 1-day buffer mirrors what
    # asf_enumeration.aria_s1_gunw.get_product uses internally.
    #
    # max_results caps the total to avoid unbounded pagination.  CMR
    # itself has no hard result-count limit, but network timeouts
    # (30 s / page) and memory become a concern for very large sets.
    # ------------------------------------------------------------------
    later_dates = [max(a, b) for a, b in pairs]
    date_buffer = datetime.timedelta(days=1)
    query_start = min(later_dates) - date_buffer
    query_end = max(later_dates) + date_buffer

    LOGGER.info('Checking ASF archive for existing products '
                '(%d pairs, bulk query %s → %s, max_results=%d) ...',
                len(pairs), query_start, query_end, max_results)

    truncated = False
    try:
        results = asf_search.search(
            dataset=asf_search.constants.DATASET.ARIA_S1_GUNW,
            frame=int(frame_id),
            start=query_start,
            end=query_end,
            maxResults=max_results,
        )
        if len(results) >= max_results:
            truncated = True
            LOGGER.warning(
                'Bulk query hit the %d-result safety cap — results may '
                'be incomplete.  Unmatched pairs will be checked '
                'individually.', max_results)
        LOGGER.info('Bulk query returned %d existing products.', len(results))
    except Exception as exc:
        LOGGER.warning('Bulk ASF query failed (%s); '
                       'marking all pairs as new.', exc)
        return {'existing': [], 'new': list(pairs)}

    # Build a set of (reference_date, secondary_date) from scene names.
    # Scene name format: ...-YYYYMMDD_YYYYMMDD-...
    # The first date is the reference (later), the second is secondary.
    existing_set = set()
    for product in results:
        scene = product.properties.get('sceneName', '')
        try:
            date_part = scene.split('-')[6]  # e.g. '20230327_20230315'
            ref_str, sec_str = date_part.split('_')
            ref_dt = datetime.datetime.strptime(ref_str, '%Y%m%d').date()
            sec_dt = datetime.datetime.strptime(sec_str, '%Y%m%d').date()
            existing_set.add((ref_dt, sec_dt))
        except (IndexError, ValueError) as exc:
            LOGGER.debug('Could not parse dates from scene name %r: %s',
                         scene, exc)

    LOGGER.debug('Parsed %d unique (ref, sec) pairs from ASF results.',
                 len(existing_set))

    # Match each of our pairs against the bulk result set.
    # Our pairs are (earlier, later); GUNW convention is
    # (reference=later, secondary=earlier).
    unmatched = []
    for earlier, later in pairs:
        if (later, earlier) in existing_set:
            existing.append((earlier, later))
        else:
            unmatched.append((earlier, later))

    # ------------------------------------------------------------------
    # Fallback: if the bulk query was truncated (hit max_results cap),
    # the unmatched set may contain pairs whose products were not in the
    # truncated result.  Re-query in chunks of CMR_FALLBACK_PAGE (250,
    # the CMR default page size) scoped to only the unmatched reference
    # dates, rather than falling back to slow per-pair look-ups.
    # ------------------------------------------------------------------
    CMR_FALLBACK_PAGE = 250
    if truncated and unmatched:
        LOGGER.info(
            'Bulk results were truncated — re-querying %d unmatched '
            'pair(s) in chunks of %d ...', len(unmatched), CMR_FALLBACK_PAGE)

        # Group unmatched pairs by their reference (later) date so we
        # can issue narrower date-range queries.
        unmatched_later = sorted({max(a, b) for a, b in unmatched})
        for chunk_start in range(0, len(unmatched_later), CMR_FALLBACK_PAGE):
            chunk_dates = unmatched_later[
                chunk_start:chunk_start + CMR_FALLBACK_PAGE]
            q_start = min(chunk_dates) - date_buffer
            q_end = max(chunk_dates) + date_buffer
            try:
                extra = asf_search.search(
                    dataset=asf_search.constants.DATASET.ARIA_S1_GUNW,
                    frame=int(frame_id),
                    start=q_start,
                    end=q_end,
                )
                for product in extra:
                    scene = product.properties.get('sceneName', '')
                    try:
                        dp = scene.split('-')[6]
                        rs, ss = dp.split('_')
                        existing_set.add((
                            datetime.datetime.strptime(rs, '%Y%m%d').date(),
                            datetime.datetime.strptime(ss, '%Y%m%d').date(),
                        ))
                    except (IndexError, ValueError):
                        pass
            except Exception as exc:
                LOGGER.debug('Fallback chunk query failed: %s', exc)

        # Re-evaluate unmatched pairs with the augmented existing_set.
        still_new = []
        for earlier, later in unmatched:
            if (later, earlier) in existing_set:
                existing.append((earlier, later))
            else:
                still_new.append((earlier, later))
        new = still_new
    else:
        new = unmatched

    LOGGER.info('Existing: %d | New: %d', len(existing), len(new))
    return {'existing': existing, 'new': new}


def get_stack_for_frame(frame_id, start_date, end_date):
    """
    Retrieve SLC acquisition dates and perpendicular baselines for a frame.

    Uses ``asf_enumeration.aria_s1_gunw.get_acquisitions`` for dates and
    optionally ``asf_search`` baseline stack for *bperp* values.

    Parameters
    ----------
    frame_id : int
        ARIA frame ID.
    start_date : datetime.date
        Start date (inclusive).
    end_date : datetime.date
        End date (inclusive).

    Returns
    -------
    dict
        ``{datetime.date: float}``  —  acquisition date → *bperp* (m).
        *bperp* is 0.0 when baseline info is unavailable.
    """
    # Get acquisition dates via asf_enumeration (primary)
    dates = get_acquisitions_for_frame(frame_id, start_date, end_date)
    dates_bperp = {d: 0.0 for d in dates}

    # Enrich with perpendicular baselines from asf_search baseline stack
    if HAS_ASF_SEARCH:
        try:
            from asf_search import ARIAS1GUNWProduct
            LOGGER.info('Computing perpendicular baselines for frame %s '
                        'via asf_search ...', frame_id)
            groups = ARIAS1GUNWProduct.get_aria_groups_for_frame(
                int(frame_id))
            if groups:
                stack = list(groups)
                # Use a product near the middle of the stack as reference
                reference = stack[len(stack) // 2]
                enriched, _warnings = (
                    asf_search.baseline_search.get_baseline_from_stack(
                        reference, stack))
                for prod in enriched:
                    props = prod.properties
                    acq_str = (props.get('startTime') or
                               props.get('sceneDateString', ''))[:10]
                    if not acq_str:
                        continue
                    acq_date = datetime.datetime.strptime(
                        acq_str, '%Y-%m-%d').date()
                    if acq_date in dates_bperp:
                        bperp = props.get('perpendicularBaseline', 0) or 0
                        dates_bperp[acq_date] = float(bperp)
                bperp_count = sum(1 for v in dates_bperp.values() if v != 0)
                LOGGER.info('Enriched %d / %d dates with bperp.',
                            bperp_count, len(dates_bperp))
        except Exception as exc:
            LOGGER.warning('Baseline enrichment failed: %s — '
                           'bperp will be 0.', exc)

    return dates_bperp


def save_pairs_csv(pairs_dict, frame_id, dates_bperp, output_dir='./'):
    """
    Write interferogram pairs to a CSV file.

    Columns: ``reference_date, secondary_date, temporal_baseline_days,
    perpendicular_baseline_m, exists_at_asf``

    Parameters
    ----------
    pairs_dict : dict
        ``{'existing': [...], 'new': [...]}`` pair lists.
    frame_id : int
        ARIA frame ID (used in the filename).
    dates_bperp : dict
        ``{datetime.date: float}`` perpendicular baseline map.
    output_dir : str
        Output directory.

    Returns
    -------
    str
        Path to the written CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir,
                               f'aria_pairs_frame{frame_id}.csv')

    with open(output_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['reference_date', 'secondary_date',
                         'temporal_baseline_days',
                         'perpendicular_baseline_m',
                         'exists_at_asf'])
        for exists_flag, pairs in [('True', pairs_dict['existing']),
                                   ('False', pairs_dict['new'])]:
            for ref, sec in pairs:
                temp_bl = (sec - ref).days
                bperp_ref = dates_bperp.get(ref, 0)
                bperp_sec = dates_bperp.get(sec, 0)
                bperp_pair = round(bperp_sec - bperp_ref, 2)
                writer.writerow([ref.strftime('%Y-%m-%d'),
                                 sec.strftime('%Y-%m-%d'),
                                 temp_bl, bperp_pair, exists_flag])

    LOGGER.info('Saved pairs CSV to: %s', output_path)
    return output_path


def plot_baseline(dates_bperp, pairs_dict, frame_id, output_dir='./',
                  network_type='sequential'):
    """
    Create a perpendicular-baseline vs. time plot for the pair network.

    Acquisitions are shown as markers; existing pairs as light-grey lines
    and new pairs as coloured lines.

    Parameters
    ----------
    dates_bperp : dict
        ``{datetime.date: float}`` of acquisition date → bperp (m).
    pairs_dict : dict
        ``{'existing': [...], 'new': [...]}`` pair lists.
    frame_id : int
        ARIA frame ID.
    output_dir : str
        Output directory.
    network_type : str
        Label for the title.
    """
    existing = pairs_dict['existing']
    new = pairs_dict['new']

    fig, ax = plt.subplots(figsize=(14, 6))

    # ----- pair lines ----- #
    for ref, sec in existing:
        ax.plot([ref, sec],
                [dates_bperp.get(ref, 0), dates_bperp.get(sec, 0)],
                color='lightgrey', linewidth=0.8, zorder=1)

    for ref, sec in new:
        ax.plot([ref, sec],
                [dates_bperp.get(ref, 0), dates_bperp.get(sec, 0)],
                color='steelblue', linewidth=1.0, alpha=0.7, zorder=2)

    # ----- acquisition markers ----- #
    acq_dates = sorted(dates_bperp.keys())
    acq_bperp = [dates_bperp[d] for d in acq_dates]
    ax.scatter(acq_dates, acq_bperp, c='black', s=20, zorder=5)

    # ----- legend ----- #
    handles = [
        mlines.Line2D([], [], color='black', marker='o', linestyle='None',
                      markersize=5,
                      label=f'Acquisitions ({len(acq_dates)})'),
        mlines.Line2D([], [], color='lightgrey', linewidth=2,
                      label=f'Existing ({len(existing)})'),
        mlines.Line2D([], [], color='steelblue', linewidth=2,
                      label=f'New ({len(new)})'),
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=10)

    # ----- axes formatting ----- #
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    
    # Prevent duplicate month labels on short timeframes by forcing a month locator
    time_span_days = (acq_dates[-1] - acq_dates[0]).days
    if time_span_days < 180:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
    else:
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        
    fig.autofmt_xdate(rotation=45, ha='right')

    ax.set_xlabel('Date', fontsize=12)
    ax.set_ylabel('Perpendicular Baseline (m)', fontsize=12)
    ax.set_title(f'ARIA S1 GUNW — Frame {frame_id}  '
                 f'({network_type.title()} network)',
                 fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.3)

    # ----- save ----- #
    os.makedirs(output_dir, exist_ok=True)
    base_stem = os.path.join(output_dir,
                             f'aria_baseline_frame{frame_id}')
    # Find a unique stem (avoid overwriting previous runs)
    stem = base_stem
    counter = 1
    while (os.path.exists(f'{stem}.png')
           or os.path.exists(f'{stem}.eps')):
        stem = f'{base_stem}_{counter}'
        counter += 1

    for ext in ('.png', '.eps'):
        out = f'{stem}{ext}'
        fig.savefig(out, dpi=150, bbox_inches='tight')
        LOGGER.info('Saved baseline plot to: %s', out)
    plt.close(fig)


def print_pairs_summary(pairs_dict, dates, frame_id):
    """
    Print a console summary of the pair network including credit costs.

    Parameters
    ----------
    pairs_dict : dict
        ``{'existing': [...], 'new': [...]}`` pair lists.
    dates : list
        Sorted acquisition dates.
    frame_id : int
        ARIA frame ID.
    """
    existing = pairs_dict['existing']
    new = pairs_dict['new']
    total = len(existing) + len(new)
    credit_cost = len(new) * CREDITS_PER_PAIR

    print('\n' + '=' * 70)
    print('INTERFEROGRAM PAIR SUMMARY')
    print('=' * 70)
    print(f'Frame ID          : {frame_id}')
    print(f'Date range        : {min(dates)} to {max(dates)}')
    print(f'Total acquisitions: {len(dates)}')
    print(f'Total pairs       : {total}')
    print(f'  Existing at ASF : {len(existing)}')
    print(f'  New (to order)  : {len(new)}')
    print('-' * 70)

    if new:
        baselines = [(sec - ref).days for ref, sec in new]
        print(f'New pair statistics:')
        print(f'  Min temporal baseline : {min(baselines)} days')
        print(f'  Max temporal baseline : {max(baselines)} days')
        print(f'  Mean temporal baseline: {np.mean(baselines):.1f} days')
        print()

    print(f'Estimated credit cost   : {credit_cost:,}  '
          f'({len(new)} pairs x {CREDITS_PER_PAIR} credits)')
    print(f'ASF monthly quota       : {ASF_MONTHLY_QUOTA:,} credits')
    if credit_cost > ASF_MONTHLY_QUOTA:
        pct_quota = credit_cost / ASF_MONTHLY_QUOTA * 100
        print(f'  ** WARNING: cost exceeds monthly quota by '
              f'{credit_cost - ASF_MONTHLY_QUOTA:,} credits '
              f'({pct_quota:.0f}% of quota) **')
    elif credit_cost > 0:
        pct_quota = credit_cost / ASF_MONTHLY_QUOTA * 100
        print(f'  Uses {pct_quota:.1f}% of monthly quota')

    # Report user's current remaining credits if hyp3_sdk is available
    if HAS_HYP3_SDK:
        try:
            hyp3 = hyp3_sdk.HyP3()
            credits_remaining = hyp3.check_credits()
            pct_remaining = (credits_remaining / ASF_MONTHLY_QUOTA * 100
                             if ASF_MONTHLY_QUOTA > 0 else 0)
            print(f'Your remaining credits  : {credits_remaining:,}  '
                  f'({pct_remaining:.1f}% of monthly quota)')
            if credit_cost > 0 and credits_remaining > 0:
                pct_of_remaining = credit_cost / credits_remaining * 100
                print(f'  This order would use {pct_of_remaining:.1f}% '
                      f'of your remaining credits')
            if credit_cost > credits_remaining:
                print(f'  ** WARNING: order cost ({credit_cost:,}) exceeds '
                      f'your remaining credits ({credits_remaining:,}). **')
                print(f'  Consider using --limit with --orderpairs to '
                      f'submit in smaller batches.')
        except Exception:
            LOGGER.debug('Could not query HyP3 credits (not logged in?).')
    print('=' * 70 + '\n')


# =============================================================================
# Frame Information Functions
# =============================================================================

def print_frame_info(frames):
    """
    Print information about matching frames to stdout.

    Parameters
    ----------
    frames : list
        List of frame features to display.
    """
    if not frames:
        LOGGER.info('No frames found matching the criteria.')
        return

    print('\n' + '=' * 70)
    print('MATCHING ARIA SENTINEL-1 GUNW FRAMES')
    print('=' * 70)
    print(f'Found {len(frames)} frames matching criteria:\n')

    # Print header
    print(f'{"Frame ID":<12} {"Track":<8} {"Direction":<12} '
          f'{"Bounds (W, S, E, N)"}')
    print('-' * 70)

    for feature in frames:
        geom = shape(feature['geometry'])
        frame_id = get_frame_property(feature, _FRAME_ID_KEYS)
        track = get_frame_property(feature, _TRACK_KEYS)
        direction = get_frame_property(feature, _DIRECTION_KEYS)

        bounds = geom.bounds  # (W, S, E, N)
        bounds_str = (f'({bounds[0]:.2f}, {bounds[1]:.2f}, '
                      f'{bounds[2]:.2f}, {bounds[3]:.2f})')

        print(f'{str(frame_id):<12} {str(track):<8} {str(direction):<12} '
              f'{bounds_str}')

    print('-' * 70)
    print(f'Total: {len(frames)} frames')
    print('=' * 70 + '\n')


def save_frames_csv(frames, output_dir='./'):
    """
    Save frame information to a CSV file.

    Parameters
    ----------
    frames : list
        List of frame features.
    output_dir : str
        Directory to write the CSV file.

    Returns
    -------
    str
        Path to the written CSV file.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'aria_frames.csv')

    with open(output_path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['frame_id', 'track', 'direction',
                         'bounds_W', 'bounds_S', 'bounds_E', 'bounds_N'])
        for feature in frames:
            geom = shape(feature['geometry'])
            frame_id = get_frame_property(feature, _FRAME_ID_KEYS)
            track = get_frame_property(feature, _TRACK_KEYS)
            direction = get_frame_property(feature, _DIRECTION_KEYS)
            bounds = geom.bounds
            writer.writerow([frame_id, track, direction,
                             f'{bounds[0]:.4f}', f'{bounds[1]:.4f}',
                             f'{bounds[2]:.4f}', f'{bounds[3]:.4f}'])

    LOGGER.info('Saved frame list to: %s', output_path)
    return output_path


# Natural Earth GeoJSON URLs (110m resolution — lightweight)
_NATURAL_EARTH_URLS = {
    'coastline': (
        'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
        'master/geojson/ne_110m_coastline.geojson'),
    'land': (
        'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
        'master/geojson/ne_110m_land.geojson'),
    'borders': (
        'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/'
        'master/geojson/ne_110m_admin_0_boundary_lines_land.geojson'),
}


def _get_naturalearth(name, cache_dir='.'):
    """
    Return a list of shapely geometries for a Natural Earth dataset.

    Downloads the GeoJSON on first call and caches it in ``cache_dir``.
    Returns an empty list if the download fails (no network, timeout,
    write-permission issues, etc.).

    Parameters
    ----------
    name : str
        One of ``'coastline'``, ``'land'``, ``'borders'``.
    cache_dir : str
        Directory in which to cache downloaded GeoJSON files.
        Defaults to the current working directory.

    Returns
    -------
    list[shapely.geometry.base.BaseGeometry]
    """
    if name not in _NATURAL_EARTH_URLS:
        LOGGER.warning('Unknown Natural Earth dataset: %s', name)
        return []

    cache_path = os.path.join(cache_dir, f'ne_110m_{name}.geojson')

    # Try reading from cache first
    if os.path.isfile(cache_path):
        try:
            with open(cache_path) as fh:
                data = json.load(fh)
            return [shape(f['geometry']) for f in data['features']]
        except Exception:
            LOGGER.debug('Cached file corrupt, re-downloading: %s',
                         cache_path)

    # Download
    url = _NATURAL_EARTH_URLS[name]
    try:
        LOGGER.debug('Downloading Natural Earth %s from %s', name, url)
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        data = resp.json()
    except Exception as exc:
        LOGGER.warning(
            'Could not download Natural Earth %s data: %s. '
            'Map will be plotted without coastlines/borders.', name, exc)
        return []

    # Cache for next time
    try:
        os.makedirs(cache_dir, exist_ok=True)
        with open(cache_path, 'w') as fh:
            json.dump(data, fh)
        LOGGER.debug('Cached Natural Earth %s to %s', name, cache_path)
    except Exception as exc:
        LOGGER.debug('Could not cache Natural Earth data: %s', exc)

    return [shape(f['geometry']) for f in data['features']]


def plot_frames(frames, bbox_poly=None, output_dir='./'):
    """
    Create a map showing ARIA frame outlines with a unique colour per
    frame and a legend placed outside the axes.

    Parameters
    ----------
    frames : list
        List of GeoJSON frame features to plot.
    bbox_poly : shapely.geometry.Polygon, optional
        User-supplied bounding box (drawn as an orange rectangle).
    output_dir : str
        Directory in which to save the plot.
    """
    if not frames:
        LOGGER.warning('No frames to plot.')
        return

    # ----- compute extent from frames + bbox ----- #
    all_bounds = []
    if bbox_poly:
        all_bounds.append(bbox_poly.bounds)
    for feat in frames:
        all_bounds.append(shape(feat['geometry']).bounds)

    min_x = min(b[0] for b in all_bounds)
    min_y = min(b[1] for b in all_bounds)
    max_x = max(b[2] for b in all_bounds)
    max_y = max(b[3] for b in all_bounds)

    x_buf = (max_x - min_x) * 0.10 or 1.0
    y_buf = (max_y - min_y) * 0.10 or 1.0
    extent = [min_x - x_buf, max_x + x_buf,
              min_y - y_buf, max_y + y_buf]

    # ----- create figure ----- #
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, linestyle='--', alpha=0.5)

    # ----- Natural Earth background (coastlines, land, borders) ----- #
    land_geoms = _get_naturalearth('land', cache_dir=output_dir)
    for geom in land_geoms:
        polys = geom.geoms if geom.geom_type == 'MultiPolygon' else [geom]
        for poly in polys:
            x, y = poly.exterior.xy
            ax.fill(x, y, facecolor='lightgray', alpha=0.4, edgecolor='none')

    coast_geoms = _get_naturalearth('coastline', cache_dir=output_dir)
    for geom in coast_geoms:
        lines = (geom.geoms if geom.geom_type == 'MultiLineString'
                 else [geom])
        for line in lines:
            x, y = line.xy
            ax.plot(x, y, color='black', linewidth=0.6, alpha=0.7)

    border_geoms = _get_naturalearth('borders', cache_dir=output_dir)
    for geom in border_geoms:
        lines = (geom.geoms if geom.geom_type == 'MultiLineString'
                 else [geom])
        for line in lines:
            x, y = line.xy
            ax.plot(x, y, color='gray', linewidth=0.4, linestyle='--',
                    alpha=0.6)

    # ----- assign a unique colour to each frame ----- #
    n_frames = len(frames)
    cmap = plt.colormaps.get_cmap('tab20' if n_frames <= 20 else 'hsv')
    colours = [cmap(i / max(n_frames - 1, 1)) for i in range(n_frames)]

    legend_handles = []

    for idx, feature in enumerate(frames):
        geom = shape(feature['geometry'])
        frame_id = get_frame_property(feature, _FRAME_ID_KEYS)
        colour = colours[idx]

        # Draw outline(s) – no fill
        polys = (geom.geoms if geom.geom_type == 'MultiPolygon'
                 else [geom])
        for poly in polys:
            x, y = poly.exterior.xy
            ax.plot(x, y, color=colour, linewidth=2)

        # Label at centroid
        centroid = geom.centroid
        ax.text(centroid.x, centroid.y, str(frame_id),
                fontsize=8, ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.2',
                          facecolor='white', alpha=0.7))

        # Legend entry
        legend_handles.append(
            mlines.Line2D([], [], color=colour, linewidth=2,
                          label=str(frame_id)))

    # ----- bounding box rectangle ----- #
    if bbox_poly:
        bx, by = bbox_poly.exterior.xy
        ax.plot(bx, by, color='orange', linewidth=3, linestyle='-')
        legend_handles.insert(
            0, mlines.Line2D([], [], color='orange', linewidth=3,
                             label='Bounding Box'))

    # ----- legend outside the axes ----- #
    ax.legend(handles=legend_handles, loc='upper left',
              bbox_to_anchor=(1.02, 1.0), fontsize=9, framealpha=0.9,
              borderaxespad=0)

    # ----- title & labels ----- #
    ax.set_title(f'ARIA Sentinel-1 GUNW Frames '
                 f'({n_frames} matching)',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # ----- save ----- #
    os.makedirs(output_dir, exist_ok=True)
    base_stem = os.path.join(output_dir, 'aria_frames_map')
    stem = base_stem
    counter = 1
    while (os.path.exists(f'{stem}.png')
           or os.path.exists(f'{stem}.eps')):
        stem = f'{base_stem}_{counter}'
        counter += 1

    for ext in ('.png', '.eps'):
        out = f'{stem}{ext}'
        fig.savefig(out, dpi=150, bbox_inches='tight')
        LOGGER.info('Saved frame map to: %s', out)
    plt.close(fig)


# =============================================================================
# Ordering Functions
# =============================================================================

def _read_pairs_csv(pairs_file):
    """
    Read the pairs CSV and return rows where ``exists_at_asf == False``.

    Parameters
    ----------
    pairs_file : str
        Path to the pairs CSV file (output of ``--getpairs``).

    Returns
    -------
    list[tuple[str, str]]
        List of ``(reference_date, secondary_date)`` string pairs to order.
        Dates are in ``YYYY-MM-DD`` format with reference = earlier date.
    """
    if not os.path.isfile(pairs_file):
        raise FileNotFoundError(f'Pairs file not found: {pairs_file}')

    pairs_to_order = []
    with open(pairs_file, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get('exists_at_asf', '').strip().lower() == 'false':
                pairs_to_order.append(
                    (row['reference_date'], row['secondary_date']))
    return pairs_to_order


def _save_jobs_csv(jobs_info, frame_id, job_name, output_dir='./'):
    """
    Save submitted job details to a CSV for tracking.

    Parameters
    ----------
    jobs_info : list[dict]
        Each dict has keys: ``job_id``, ``reference_date``,
        ``secondary_date``, ``status_code``.
    frame_id : int
        ARIA frame ID.
    job_name : str
        HyP3 job name.
    output_dir : str
        Output directory.

    Returns
    -------
    str
        Path to the written CSV.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(
        output_dir, f'aria_jobs_frame{frame_id}.csv')

    # Append if file exists, otherwise create with header
    file_exists = os.path.isfile(output_path)
    with open(output_path, 'a', newline='') as fh:
        writer = csv.writer(fh)
        if not file_exists:
            writer.writerow(['job_id', 'job_name', 'frame_id',
                             'reference_date', 'secondary_date',
                             'status_code', 'submitted_at'])
        now_str = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
        for info in jobs_info:
            writer.writerow([
                info['job_id'], job_name, frame_id,
                info['reference_date'], info['secondary_date'],
                info['status_code'], now_str])

    LOGGER.info('Saved job tracking CSV to: %s', output_path)
    return output_path


def order_pairs(frame_id, pairs_file, job_name=None, dry_run=False,
                limit=None, output_dir='./'):
    """
    Submit interferogram pairs to HyP3 for ARIA S1 GUNW processing.

    Reads the pairs CSV produced by ``--getpairs``, filters to rows where
    ``exists_at_asf == False``, and submits them via the HyP3 SDK.

    Safety features:
    - ``dry_run=True``: preview what would be submitted without spending
      credits.
    - ``limit``: cap the number of pairs submitted per invocation.
    - Interactive confirmation prompt before actual submission.

    Parameters
    ----------
    frame_id : int
        ARIA frame ID.
    pairs_file : str
        Path to the pairs CSV file (output of ``--getpairs``).
    job_name : str, optional
        HyP3 job name.  If *None*, auto-generates
        ``ARIA_frame{ID}_{YYYYMMDD_HHMMSS}``.
    dry_run : bool
        If True, show what would be submitted but do not submit.
    limit : int or None
        Maximum number of pairs to submit.  None means no limit.
    output_dir : str
        Directory for the job-tracking CSV output.

    Returns
    -------
    object or None
        The HyP3 batch object, or *None* if nothing was submitted.
    """
    if not HAS_HYP3_SDK:
        raise ImportError(
            'hyp3_sdk is required for --orderpairs.  Install it with:\n'
            '  pip install hyp3_sdk')

    # --- read pairs CSV --- #
    pairs_to_order = _read_pairs_csv(pairs_file)

    if not pairs_to_order:
        print('All pairs already exist at ASF — nothing to order.')
        return None

    # --- enforce date order: (earlier, later) --- #
    # The CSV may have been hand-edited.  Normalise so that ref_str is
    # always the earlier date and sec_str the later date.  The swap to
    # match the ARIA GUNW convention (reference = later, secondary =
    # earlier) happens at submission time.
    normalised = []
    for ref_str, sec_str in pairs_to_order:
        d1 = datetime.datetime.strptime(ref_str, '%Y-%m-%d').date()
        d2 = datetime.datetime.strptime(sec_str, '%Y-%m-%d').date()
        if d1 > d2:
            LOGGER.debug('Swapping misordered pair: %s / %s', ref_str, sec_str)
            ref_str, sec_str = sec_str, ref_str
        normalised.append((ref_str, sec_str))
    pairs_to_order = normalised

    # --- apply limit --- #
    total_new = len(pairs_to_order)
    if limit is not None and limit > 0:
        pairs_to_order = pairs_to_order[:limit]

    # --- auto job name --- #
    if not job_name:
        now = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d_%H%M%S')
        job_name = f'ARIA_frame{frame_id}_{now}'

    # --- authenticate --- #
    hyp3 = hyp3_sdk.HyP3()  # uses ~/.netrc by default

    credits_available = hyp3.check_credits()
    credit_cost = len(pairs_to_order) * CREDITS_PER_PAIR

    # --- summary --- #
    print(f'\n{"=" * 60}')
    print(f'  ARIA S1 GUNW Order Summary — Frame {frame_id}')
    print(f'{"=" * 60}')
    print(f'  Pairs file     : {pairs_file}')
    print(f'  Total new pairs: {total_new}')
    if limit is not None:
        print(f'  Limit applied  : {limit}')
    print(f'  Submitting     : {len(pairs_to_order)} pair(s)')
    print(f'  Job name       : {job_name}')
    print(f'  Credit cost    : {credit_cost:,}  '
          f'({len(pairs_to_order)} x {CREDITS_PER_PAIR})')
    print(f'  Credits avail. : {credits_available:,}')
    remaining = credits_available - credit_cost
    print(f'  After submit   : {remaining:,} credits remaining')
    print(f'{"=" * 60}')

    # --- pair listing --- #
    print(f'\n  {"#":<4} {"Reference (earlier)":<22} {"Secondary (later)":<22}')
    print(f'  {"-" * 4} {"-" * 22} {"-" * 22}')
    for i, (ref_str, sec_str) in enumerate(pairs_to_order, 1):
        print(f'  {i:<4} {ref_str:<22} {sec_str:<22}')

    # --- dry run exit --- #
    if dry_run:
        print(f'\n  ** DRY RUN — no jobs submitted, no credits spent. **')
        print(f'  Remove --dry-run to submit these {len(pairs_to_order)} '
              f'job(s).\n')
        return None

    # --- credit check --- #
    if credit_cost > credits_available:
        raise SystemExit(
            f'\n  ** ORDER BLOCKED: insufficient credits. **\n'
            f'  Credits needed   : {credit_cost:,}\n'
            f'  Credits available: {credits_available:,}\n'
            f'  Shortfall        : {credit_cost - credits_available:,}\n\n'
            f'  Reduce the number of pairs with --limit.  For example:\n'
            f'    ariaOrderASF.py --orderpairs --frame {frame_id} '
            f'--pairs-file {pairs_file} --limit '
            f'{int(credits_available // CREDITS_PER_PAIR)}\n')

    # --- interactive confirmation --- #
    try:
        answer = input(
            f'\n  Submit {len(pairs_to_order)} job(s) for '
            f'{credit_cost:,} credits? [y/N] ').strip().lower()
    except EOFError:
        answer = ''
    if answer not in ('y', 'yes'):
        print('  Submission cancelled.')
        return None

    # --- prepare and submit --- #
    prepared = []
    for ref_str, sec_str in pairs_to_order:
        # CSV stores pairs as (earlier, later).  HyP3 follows the ARIA
        # GUNW convention where reference = later date, secondary =
        # earlier date, so we swap the order.
        prepared.append(
            hyp3.prepare_aria_s1_gunw_job(
                reference_date=sec_str,
                secondary_date=ref_str,
                frame_id=int(frame_id),
                name=job_name))

    batch = hyp3.submit_prepared_jobs(prepared)

    # --- extract job details --- #
    jobs_info = []
    for job, (ref_str, sec_str) in zip(batch, pairs_to_order):
        job_dict = job.to_dict()
        jobs_info.append({
            'job_id': job_dict.get('job_id', 'unknown'),
            'reference_date': ref_str,
            'secondary_date': sec_str,
            'status_code': job_dict.get('status_code', 'PENDING'),
        })

    # --- print results --- #
    print(f'\n  Submitted {len(batch)} job(s).')
    print(f'  Job name: {job_name}')
    print(f'\n  {"Job ID":<40} {"Ref Date":<14} {"Sec Date":<14} '
          f'{"Status"}')
    print(f'  {"-" * 40} {"-" * 14} {"-" * 14} {"-" * 10}')
    for info in jobs_info:
        print(f'  {info["job_id"]:<40} {info["reference_date"]:<14} '
              f'{info["secondary_date"]:<14} {info["status_code"]}')

    # --- save job-tracking CSV --- #
    _save_jobs_csv(jobs_info, frame_id, job_name, output_dir=output_dir)

    print(f'\n  Monitor at: '
          f'https://search.asf.alaska.edu/#/?searchType=On%20Demand')
    print(f'  Or run:  ariaOrderASF.py --statusjobs '
          f'--status-name {job_name}\n')

    return batch


# =============================================================================
# Job Status Functions
# =============================================================================

def status_jobs(status_name=None, job_ids=None):
    """
    Query and display the status of HyP3 jobs.

    At least one of ``status_name`` or ``job_ids`` must be provided.
    Jobs can be looked up by their shared job name (used during submission)
    or by individual job IDs.

    Parameters
    ----------
    status_name : str or None
        HyP3 job name to query (returns all matching jobs).
    job_ids : str or None
        Comma-separated string of HyP3 job IDs to look up individually.
    """
    if not HAS_HYP3_SDK:
        raise ImportError(
            'hyp3_sdk is required for --statusjobs.  Install it with:\n'
            '  pip install hyp3_sdk')

    if not status_name and not job_ids:
        raise ValueError(
            '--statusjobs requires at least one of '
            '--status-name or --job-ids.')

    hyp3 = hyp3_sdk.HyP3()

    jobs = []

    # --- query by job name --- #
    if status_name:
        LOGGER.info('Querying jobs with name: %s', status_name)
        batch = hyp3.find_jobs(name=status_name)
        jobs.extend(batch)
        LOGGER.info('Found %d job(s) for name "%s".', len(batch), status_name)

    # --- query by individual job IDs --- #
    if job_ids:
        id_list = [jid.strip() for jid in job_ids.split(',') if jid.strip()]
        LOGGER.info('Querying %d individual job ID(s).', len(id_list))
        for jid in id_list:
            try:
                job = hyp3.get_job_by_id(jid)
                # Avoid duplicates if also found by name
                existing_ids = {j.to_dict().get('job_id') for j in jobs}
                if jid not in existing_ids:
                    jobs.append(job)
            except Exception as exc:
                LOGGER.warning('Could not retrieve job %s: %s', jid, exc)

    if not jobs:
        print('\n  No jobs found.\n')
        return

    # --- build status summary --- #
    status_counts = {}
    job_details = []
    for job in jobs:
        d = job.to_dict()
        status = d.get('status_code', 'UNKNOWN')
        status_counts[status] = status_counts.get(status, 0) + 1

        # Extract job parameters (reference/secondary dates, frame_id)
        params = d.get('job_parameters', {})
        job_details.append({
            'job_id': d.get('job_id', 'unknown'),
            'name': d.get('name', ''),
            'status_code': status,
            'reference_date': params.get('reference_date', ''),
            'secondary_date': params.get('secondary_date', ''),
            'frame_id': params.get('frame_id', ''),
            'credit_cost': d.get('credit_cost', ''),
            'request_time': d.get('request_time', ''),
            'expiration_time': d.get('expiration_time', ''),
        })

    # --- display --- #
    print(f'\n{"=" * 72}')
    print(f'  HyP3 Job Status — {len(jobs)} job(s)')
    print(f'{"=" * 72}')

    # Status counts
    for status, count in sorted(status_counts.items()):
        print(f'  {status:<12}: {count}')

    # Job details table
    print(f'\n  {"Job ID":<40} {"Frame":<8} {"Ref Date":<14} '
          f'{"Sec Date":<14} {"Status":<12} {"Credits"}')
    print(f'  {"-" * 40} {"-" * 8} {"-" * 14} {"-" * 14} '
          f'{"-" * 12} {"-" * 8}')
    for info in job_details:
        print(f'  {info["job_id"]:<40} {str(info["frame_id"]):<8} '
              f'{info["reference_date"]:<14} '
              f'{info["secondary_date"]:<14} '
              f'{info["status_code"]:<12} {info["credit_cost"]}')

    # Total credit cost
    total_credits = sum(
        info['credit_cost'] for info in job_details
        if isinstance(info['credit_cost'], (int, float)))
    if total_credits > 0:
        print(f'\n  Total credits: {total_credits:,}')

    print(f'\n  Monitor at: '
          f'https://search.asf.alaska.edu/#/?searchType=On%20Demand\n')


def main():
    """Main entry point for ariaOrderASF."""
    parser = create_parser()
    args = parser.parse_args()

    # ----- logging ----- #
    log_level = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
    }[args.log_level]
    logging.basicConfig(level=log_level, format=_LOG_FORMAT)
    if args.verbose:
        LOGGER.setLevel(logging.DEBUG)

    LOGGER.info('=' * 70)
    LOGGER.info('ARIA Sentinel-1 GUNW Frame Query / Order Tool')
    LOGGER.info('NOTE: This tool supports Sentinel-1 only, not NISAR.')
    LOGGER.info('=' * 70)

    output_dir = os.path.abspath(args.wd)
    os.makedirs(output_dir, exist_ok=True)

    # =====================================================================
    # Mode 1: --getframes
    # =====================================================================
    if args.getframes:
        if args.bbox is None:
            parser.error('--getframes requires -b / --bbox.')

        bbox_poly = make_bbox(args.bbox)
        LOGGER.info(
            'Bounding box: S=%.4f  N=%.4f  W=%.4f  E=%.4f',
            bbox_poly.bounds[1], bbox_poly.bounds[3],
            bbox_poly.bounds[0], bbox_poly.bounds[2])

        # --- resolve track / direction args for upstream query --- #
        path_arg = None
        if args.track:
            tracks = [int(t.strip()) for t in args.track.split(',')]
            if len(tracks) == 1:
                path_arg = tracks[0]

        dir_arg = None
        if args.flightdir:
            dir_arg = ('ASCENDING'
                       if args.flightdir.lower().startswith('a')
                       else 'DESCENDING')

        # --- primary: asf_enumeration.aria_s1_gunw.get_frames --- #
        if HAS_ASF_ENUMERATION:
            LOGGER.info('Querying frames via asf_enumeration ...')
            af_list = aria_s1_gunw.get_frames(
                bbox_poly,
                path=path_arg,
                flight_direction=dir_arg)
            # Apply multi-track filter client-side when >1 track given
            if args.track and path_arg is None:
                af_list = [af for af in af_list
                           if af.path in tracks]
            matching = [_ariaframe_to_feature(af) for af in af_list]
        else:
            # fallback: HTTP fetch of GeoJSON
            LOGGER.info('asf_enumeration not available — fetching '
                        'frames GeoJSON from GitHub ...')
            frames_geojson = fetch_aria_frames()
            matching = filter_frames_by_bbox(frames_geojson, bbox_poly)
            if args.track:
                tracks = [int(t.strip()) for t in args.track.split(',')]
                matching = filter_frames_by_track(matching, tracks)
            if args.flightdir:
                direction = ('ascending'
                             if args.flightdir.lower().startswith('a')
                             else 'descending')
                matching = filter_frames_by_direction(matching, direction)

        LOGGER.info('Found %d frames.', len(matching))

        # Outputs: console table, CSV, map
        print_frame_info(matching)
        save_frames_csv(matching, output_dir=output_dir)
        if matching:
            plot_frames(matching, bbox_poly=bbox_poly,
                        output_dir=output_dir)
        else:
            LOGGER.warning('No frames to plot.')

    # =====================================================================
    # Mode 2: --getpairs
    # =====================================================================
    elif args.getpairs:
        if args.frame is None:
            parser.error('--getpairs requires --frame.')

        frame_id = args.frame

        # Date range
        start_date = datetime.datetime.strptime(args.start, '%Y%m%d').date()
        end_date = (datetime.datetime.strptime(args.end, '%Y%m%d').date()
                    if args.end else datetime.date.today())
        LOGGER.info('Frame %s | %s → %s', frame_id, start_date, end_date)

        # Acquisition dates + perpendicular baselines
        dates_bperp = get_stack_for_frame(frame_id, start_date, end_date)
        dates = sorted(dates_bperp.keys())

        if len(dates) < 2:
            raise SystemExit(
                f'Only {len(dates)} acquisition(s) found for frame '
                f'{frame_id} — at least 2 are needed.')

        LOGGER.info('Found %d acquisitions.', len(dates))

        # Generate pairs
        network_type = args.network

        if network_type == 'sequential':
            pairs = generate_sequential_pairs(
                dates, num_neighbors=args.num_neighbors)
        elif network_type == 'seasonal':
            pairs = generate_seasonal_pairs(
                dates, seasonal_window=args.seasonal_window)
        elif network_type == 'annual':
            pairs = generate_annual_pairs(
                dates, num_neighbors=args.num_neighbors)
        else:
            raise ValueError(f'Unknown network type: {network_type}')

        if not pairs:
            raise SystemExit('No pairs generated with current settings.')

        # Check existing products at ASF
        pairs_dict = check_existing_products(pairs, frame_id)

        # Outputs: CSV, baseline plot, console summary
        save_pairs_csv(pairs_dict, frame_id, dates_bperp,
                       output_dir=output_dir)
        plot_baseline(dates_bperp, pairs_dict, frame_id,
                      output_dir=output_dir, network_type=network_type)
        print_pairs_summary(pairs_dict, dates, frame_id)

    # =====================================================================
    # Mode 3: --orderpairs
    # =====================================================================
    elif args.orderpairs:
        if args.frame is None:
            parser.error('--orderpairs requires --frame.')
        if args.pairs_file is None:
            parser.error('--orderpairs requires --pairs-file.')

        order_pairs(args.frame, args.pairs_file,
                    job_name=args.job_name,
                    dry_run=args.dry_run,
                    limit=args.limit,
                    output_dir=output_dir)

    # =====================================================================
    # Mode 4: --statusjobs
    # =====================================================================
    elif args.statusjobs:
        if not args.status_name and not args.job_ids:
            parser.error(
                '--statusjobs requires --status-name and/or --job-ids.')

        status_jobs(status_name=args.status_name,
                    job_ids=args.job_ids)

    LOGGER.info('Done.')


if __name__ == '__main__':
    main()
