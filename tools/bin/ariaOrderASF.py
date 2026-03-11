#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
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
import datetime
import json
import logging
import os
import warnings
from collections import defaultdict
from itertools import combinations

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.dates as mdates
import numpy as np
import requests
import shapely
import shapely.geometry
from shapely.geometry import shape, Polygon, MultiPolygon

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

import ARIAtools.util.log
from ARIAtools.util.shp import open_shp

LOGGER = logging.getLogger('ariaOrderASF.py')

# URL for ARIA frames geojson from ASF enumeration repository
ARIA_FRAMES_URL = (
    'https://raw.githubusercontent.com/ASFHyP3/asf-enumeration/develop/'
    'src/asf_enumeration/frame_maps/aria_frames.geojson'
)

# Network generation types
NETWORK_TYPES = ['sequential', 'seasonal', 'annual']


def create_parser():
    """
    Create argument parser for ariaOrderASF tool.

    This tool assists users in ordering ARIA Sentinel-1 GUNW products by
    identifying frames that intersect their area of interest.

    Note: Only Sentinel-1 GUNW products are supported.
    """
    parser = argparse.ArgumentParser(
        description='Command line interface to identify ARIA Sentinel-1 GUNW '
                    'frames for ordering from ASF.\n\n'
                    'NOTE: This tool supports Sentinel-1 GUNW only, '
                    'not NISAR.',
        epilog='Examples of use:\n\n'
               '\t # Find frames intersecting a bounding box\n'
               '\t ariaOrderASF.py --bbox "36.0 37.0 -118.0 -117.0"\n\n'
               '\t # Find ascending frames for a specific track\n'
               '\t ariaOrderASF.py --bbox "36.0 37.0 -118.0 -117.0" '
               '-t 064 -d ascending\n\n'
               '\t # Generate sequential network with 3 nearest neighbors\n'
               '\t ariaOrderASF.py --frames 7423 --network sequential '
               '--num-neighbors 3 --start 20200101 --end 20210101\n\n'
               '\t # Generate seasonal interferograms (same season pairs)\n'
               '\t ariaOrderASF.py --frames 7423 --network seasonal '
               '--seasonal-window 30 --start 20190101 --end 20220101\n\n'
               '\t # Generate annual nearest neighbor interferograms\n'
               '\t ariaOrderASF.py --frames 7423 --network annual '
               '--start 20200101 --end 20230101\n',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Frame selection arguments
    frame_group = parser.add_argument_group('Frame Selection')
    frame_group.add_argument(
        '-b', '--bbox', default=None, type=str,
        help='Lat/Lon Bounding box SNWE (e.g., "36.0 37.0 -118.0 -117.0"), '
             'or GDAL-readable file containing POLYGON geometry.')
    frame_group.add_argument(
        '-f', '--frames', default=None, type=str,
        help='Comma-separated list of specific ARIA frame IDs to query '
             '(e.g., "7423,7424,7425").')
    frame_group.add_argument(
        '-t', '--track', default=None, type=str,
        help='Filter by track number; single number (e.g., 064) or '
             'comma-separated list (e.g., 064,137).')
    frame_group.add_argument(
        '-d', '--direction', dest='flightdir', default=None, type=str,
        choices=['ascending', 'descending', 'a', 'd', 'A', 'D',
                 'ASCENDING', 'DESCENDING'],
        help='Filter by flight direction: ascending (a) or descending (d).')

    # Network generation arguments
    network_group = parser.add_argument_group('Network Generation')
    network_group.add_argument(
        '--network', default=None, type=str.lower,
        choices=NETWORK_TYPES,
        help='Type of interferogram network to generate: '
             '"sequential" (nearest N neighbors), '
             '"seasonal" (same season across years), '
             '"annual" (annual nearest neighbor pairs).')
    network_group.add_argument(
        '--num-neighbors', dest='num_neighbors', default=3, type=int,
        help='Number of nearest temporal neighbors for sequential network. '
             'Default is 3.')
    network_group.add_argument(
        '--seasonal-window', dest='seasonal_window', default=30, type=int,
        help='Window in days (+/-) around same day-of-year for seasonal '
             'network. Default is 30 days.')
    network_group.add_argument(
        '--temporal-baseline-max', dest='temp_baseline_max', default=None,
        type=int,
        help='Maximum temporal baseline in days. If specified, pairs '
             'exceeding this will be excluded.')
    network_group.add_argument(
        '-s', '--start', default='20140101', type=str,
        help='Start date as YYYYMMDD for acquisition search. '
             'Default is 20140101 (Sentinel-1 launch).')
    network_group.add_argument(
        '-e', '--end', default=None, type=str,
        help='End date as YYYYMMDD for acquisition search. '
             'Default is today.')

    # Output arguments
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument(
        '-w', '--workdir', dest='wd', default='./', type=str,
        help='Directory to save output plots. Default is current directory.')
    output_group.add_argument(
        '-o', '--output', default='plot', type=str.lower,
        choices=['plot', 'frames', 'both', 'network'],
        help='Output type: "plot" generates frame map, "frames" prints '
             'frame info, "both" does both, "network" generates network '
             'plot (requires --network). Default is "plot".')
    output_group.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable verbose output.')
    output_group.add_argument(
        '--log-level', default='info',
        choices=['debug', 'info', 'warning', 'error'],
        help='Logger log level. Default is "info".')

    return parser


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
        props = feature.get('properties', {})
        # Try common property names for frame ID
        frame_id = props.get('frame_id') or props.get('frameID') or \
            props.get('id') or props.get('frame')

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
        props = feature.get('properties', {})
        # Try common property names for track/relative orbit
        track = props.get('track') or props.get('track_number') or \
            props.get('relativeOrbit') or props.get('relative_orbit')

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
        props = feature.get('properties', {})
        # Try common property names for direction
        frame_dir = props.get('direction') or props.get('flightDirection') or \
            props.get('flight_direction') or props.get('orbit_direction')

        if frame_dir and frame_dir.lower() == direction_normalized:
            filtered.append(feature)

    return filtered


# =============================================================================
# Network Generation Functions
# =============================================================================

def get_acquisitions_for_frame(frame_id, start_date, end_date, track=None,
                                direction=None):
    """
    Get all Sentinel-1 acquisitions for an ARIA frame within a date range.

    Uses asf_enumeration if available, otherwise falls back to asf_search.

    Parameters
    ----------
    frame_id : int
        ARIA frame ID.
    start_date : datetime.date
        Start date for acquisition search.
    end_date : datetime.date
        End date for acquisition search.
    track : int, optional
        Track number for fallback search.
    direction : str, optional
        Flight direction for fallback search.

    Returns
    -------
    list
        List of acquisition dates (datetime.date objects), sorted.
    """
    acquisition_dates = []

    # Try using asf_enumeration first (preferred)
    if HAS_ASF_ENUMERATION:
        try:
            LOGGER.info('Using asf_enumeration to get acquisitions for '
                        'frame %s...', frame_id)
            acquisitions = aria_s1_gunw.get_acquisitions(int(frame_id))
            for acq in acquisitions:
                if start_date <= acq.date <= end_date:
                    acquisition_dates.append(acq.date)
            LOGGER.info('Found %d acquisitions via asf_enumeration.',
                        len(acquisition_dates))
            return sorted(acquisition_dates)
        except Exception as e:
            LOGGER.warning('asf_enumeration failed: %s. '
                           'Falling back to asf_search.', e)

    # Fallback to asf_search
    if HAS_ASF_SEARCH:
        LOGGER.info('Using asf_search to get acquisitions...')
        try:
            # Search for SLC products over the frame's time period
            # We use a general search since we need acquisition dates
            flight_direction = None
            if direction:
                flight_direction = (
                    'ascending' if direction.lower().startswith('a')
                    else 'descending'
                )

            results = asf_search.geo_search(
                platform=asf_search.PLATFORM.SENTINEL1,
                processingLevel=asf_search.PRODUCT_TYPE.SLC,
                relativeOrbit=track,
                flightDirection=flight_direction,
                start=start_date - datetime.timedelta(days=1),
                end=end_date + datetime.timedelta(days=1),
            )

            for result in results:
                props = result.properties
                # Get acquisition date
                start_time = props.get('startTime', '')
                if start_time:
                    acq_date = datetime.datetime.strptime(
                        start_time[:10], '%Y-%m-%d').date()
                    if start_date <= acq_date <= end_date:
                        acquisition_dates.append(acq_date)

            # Remove duplicates and sort
            acquisition_dates = sorted(set(acquisition_dates))
            LOGGER.info('Found %d unique acquisition dates via asf_search.',
                        len(acquisition_dates))
            return acquisition_dates

        except Exception as e:
            LOGGER.error('asf_search failed: %s', e)

    # If no search library available, generate approximate dates
    # Sentinel-1 repeat cycle is ~12 days
    LOGGER.warning('No ASF search library available. Generating approximate '
                   'acquisition dates based on 12-day repeat cycle.')
    current = start_date
    while current <= end_date:
        acquisition_dates.append(current)
        current += datetime.timedelta(days=12)

    return acquisition_dates


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


def check_existing_products(pairs, frame_id):
    """
    Check which interferogram pairs already exist in the ASF archive.

    Parameters
    ----------
    pairs : list
        List of (reference_date, secondary_date) tuples.
    frame_id : int
        ARIA frame ID.

    Returns
    -------
    dict
        Dictionary with keys 'existing' and 'new', each containing list of
        pairs.
    """
    existing = []
    new = []

    if HAS_ASF_ENUMERATION:
        LOGGER.info('Checking ASF archive for existing products...')
        for i, (ref_date, sec_date) in enumerate(pairs):
            try:
                exists = aria_s1_gunw.product_exists(ref_date, sec_date,
                                                     int(frame_id))
                if exists:
                    existing.append((ref_date, sec_date))
                else:
                    new.append((ref_date, sec_date))
            except Exception as e:
                LOGGER.debug('Could not check product %s_%s: %s',
                             ref_date, sec_date, e)
                new.append((ref_date, sec_date))

            # Progress update
            if (i + 1) % 50 == 0:
                LOGGER.info('Checked %d/%d pairs...', i + 1, len(pairs))
    else:
        # If asf_enumeration not available, try asf_search
        if HAS_ASF_SEARCH:
            LOGGER.info('Checking ASF archive via asf_search...')
            for ref_date, sec_date in pairs:
                try:
                    # Format interferogram pair name pattern
                    ifg_str = f'{ref_date:%Y%m%d}_{sec_date:%Y%m%d}'
                    results = asf_search.granule_search(
                        [f'*{ifg_str}*'],
                        asf_search.ASFSearchOptions(
                            dataset=asf_search.constants.ARIA_S1_GUNW
                        )
                    )
                    if results:
                        existing.append((ref_date, sec_date))
                    else:
                        new.append((ref_date, sec_date))
                except Exception:
                    new.append((ref_date, sec_date))
        else:
            LOGGER.warning('Cannot check existing products: '
                           'asf_enumeration or asf_search required.')
            new = pairs.copy()

    LOGGER.info('Found %d existing products, %d new products to create.',
                len(existing), len(new))
    return {'existing': existing, 'new': new}


def print_network_summary(pairs_dict, dates, frame_id):
    """
    Print summary of the interferogram network.

    Parameters
    ----------
    pairs_dict : dict
        Dictionary with 'existing' and 'new' pair lists.
    dates : list
        List of acquisition dates.
    frame_id : int
        ARIA frame ID.
    """
    existing = pairs_dict['existing']
    new = pairs_dict['new']
    total = len(existing) + len(new)

    print('\n' + '=' * 70)
    print('INTERFEROGRAM NETWORK SUMMARY')
    print('=' * 70)
    print(f'Frame ID: {frame_id}')
    print(f'Date Range: {min(dates)} to {max(dates)}')
    print(f'Total Acquisitions: {len(dates)}')
    print(f'Total Pairs: {total}')
    print(f'  - Existing in ASF Archive: {len(existing)}')
    print(f'  - New (to be created): {len(new)}')
    print('-' * 70)

    if new:
        # Calculate temporal baseline statistics for new pairs
        baselines = [(sec - ref).days for ref, sec in new]
        print(f'\nNew Pairs Statistics:')
        print(f'  Min Temporal Baseline: {min(baselines)} days')
        print(f'  Max Temporal Baseline: {max(baselines)} days')
        print(f'  Mean Temporal Baseline: {np.mean(baselines):.1f} days')

    print('=' * 70 + '\n')


def plot_network(dates, pairs_dict, frame_id, output_dir='./',
                 network_type='sequential'):
    """
    Create a network plot showing acquisitions and interferogram pairs.

    Parameters
    ----------
    dates : list
        List of acquisition dates.
    pairs_dict : dict
        Dictionary with 'existing' and 'new' pair lists.
    frame_id : int
        ARIA frame ID.
    output_dir : str
        Directory to save the plot.
    network_type : str
        Type of network for title.
    """
    existing = pairs_dict['existing']
    new = pairs_dict['new']

    # Convert dates to matplotlib format
    dates_sorted = sorted(dates)
    date_to_idx = {d: i for i, d in enumerate(dates_sorted)}

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10),
                                    gridspec_kw={'height_ratios': [3, 1]})

    # =========================================================================
    # Top plot: Network diagram (date vs temporal baseline)
    # =========================================================================
    # Plot acquisitions as dots
    for date in dates_sorted:
        ax1.plot(date, 0, 'ko', markersize=6, zorder=10)

    # Plot pairs as arcs
    for ref, sec in existing:
        baseline = (sec - ref).days
        mid_date = ref + datetime.timedelta(days=baseline / 2)
        # Draw arc
        ax1.annotate('', xy=(sec, 0), xytext=(ref, 0),
                     arrowprops=dict(arrowstyle='-', color='green',
                                     connectionstyle=f'arc3,rad=0.3',
                                     alpha=0.6, linewidth=1.5))

    for ref, sec in new:
        baseline = (sec - ref).days
        mid_date = ref + datetime.timedelta(days=baseline / 2)
        # Draw arc
        ax1.annotate('', xy=(sec, 0), xytext=(ref, 0),
                     arrowprops=dict(arrowstyle='-', color='red',
                                     connectionstyle=f'arc3,rad=0.3',
                                     alpha=0.6, linewidth=1.5))

    ax1.set_xlim(min(dates_sorted) - datetime.timedelta(days=30),
                 max(dates_sorted) + datetime.timedelta(days=30))

    # Adjust y-axis to show arcs
    ax1.set_ylim(-0.5, 1.5)
    ax1.set_yticks([])

    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax1.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Legend
    existing_line = mlines.Line2D([], [], color='green', linewidth=2,
                                   label=f'Existing ({len(existing)})')
    new_line = mlines.Line2D([], [], color='red', linewidth=2,
                              label=f'New ({len(new)})')
    acq_marker = mlines.Line2D([], [], color='black', marker='o',
                                linestyle='None', markersize=6,
                                label=f'Acquisitions ({len(dates)})')
    ax1.legend(handles=[acq_marker, existing_line, new_line],
               loc='upper right', fontsize=10)

    ax1.set_title(f'ARIA S1 GUNW Network - Frame {frame_id}\n'
                  f'Network Type: {network_type.title()}',
                  fontsize=14, fontweight='bold')
    ax1.set_xlabel('Date', fontsize=12)

    # =========================================================================
    # Bottom plot: Temporal baseline histogram
    # =========================================================================
    all_pairs = existing + new
    baselines_existing = [(sec - ref).days for ref, sec in existing]
    baselines_new = [(sec - ref).days for ref, sec in new]

    bins = np.arange(0, max([max(baselines_existing or [0]),
                             max(baselines_new or [0])]) + 30, 12)

    if baselines_existing:
        ax2.hist(baselines_existing, bins=bins, alpha=0.7, color='green',
                 label=f'Existing ({len(existing)})', edgecolor='darkgreen')
    if baselines_new:
        ax2.hist(baselines_new, bins=bins, alpha=0.7, color='red',
                 label=f'New ({len(new)})', edgecolor='darkred')

    ax2.set_xlabel('Temporal Baseline (days)', fontsize=12)
    ax2.set_ylabel('Number of Pairs', fontsize=12)
    ax2.legend(loc='upper right', fontsize=10)
    ax2.set_title('Temporal Baseline Distribution', fontsize=12)

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir,
                               f'aria_network_frame{frame_id}.png')

    # Avoid overwriting
    counter = 1
    base_path = output_path
    while os.path.exists(output_path):
        output_path = base_path.replace('.png', f'_{counter}.png')
        counter += 1

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    LOGGER.info('Saved network plot to: %s', output_path)

    # Also save as PDF
    pdf_path = output_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    LOGGER.info('Saved network plot to: %s', pdf_path)

    plt.close()


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

    LOGGER.info('=' * 70)
    LOGGER.info('MATCHING ARIA SENTINEL-1 GUNW FRAMES')
    LOGGER.info('=' * 70)
    LOGGER.info('Found %d frames matching criteria:\n', len(frames))

    # Print header
    print(f'{"Frame ID":<12} {"Track":<8} {"Direction":<12} '
          f'{"Bounds (W, S, E, N)"}')
    print('-' * 70)

    for feature in frames:
        props = feature.get('properties', {})
        geom = shape(feature['geometry'])

        # Extract properties with fallbacks
        frame_id = props.get('frame_id') or props.get('frameID') or \
            props.get('id') or props.get('frame') or 'N/A'
        track = props.get('track') or props.get('track_number') or \
            props.get('relativeOrbit') or props.get('relative_orbit') or 'N/A'
        direction = props.get('direction') or \
            props.get('flightDirection') or \
            props.get('flight_direction') or \
            props.get('orbit_direction') or 'N/A'

        # Get bounds
        bounds = geom.bounds  # (minx, miny, maxx, maxy) = (W, S, E, N)
        bounds_str = f'({bounds[0]:.2f}, {bounds[1]:.2f}, ' \
                     f'{bounds[2]:.2f}, {bounds[3]:.2f})'

        print(f'{str(frame_id):<12} {str(track):<8} {str(direction):<12} '
              f'{bounds_str}')

    print('-' * 70)
    print(f'Total: {len(frames)} frames')


def plot_frames(frames, bbox_poly=None, output_dir='./', verbose=False):
    """
    Create a plot showing the bounding box and matching ARIA frames
    with coastlines and country borders.

    Parameters
    ----------
    frames : list
        List of frame features to plot.
    bbox_poly : shapely.geometry.Polygon, optional
        Bounding box polygon to display.
    output_dir : str
        Directory to save the plot.
    verbose : bool
        Whether to show additional information.
    """
    # Try to use cartopy for coastlines and borders
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        use_cartopy = True
        LOGGER.debug('Using cartopy for coastlines and borders.')
    except ImportError:
        use_cartopy = False
        LOGGER.warning(
            'cartopy not installed. Plot will not include coastlines '
            'or country borders. Install cartopy for better visualization: '
            'conda install -c conda-forge cartopy'
        )

    # Determine plot extent from frames and bbox
    all_bounds = []
    if bbox_poly:
        all_bounds.append(bbox_poly.bounds)

    for feature in frames:
        geom = shape(feature['geometry'])
        all_bounds.append(geom.bounds)

    if not all_bounds:
        LOGGER.warning('No geometries to plot.')
        return

    # Calculate extent with buffer
    min_x = min(b[0] for b in all_bounds)
    min_y = min(b[1] for b in all_bounds)
    max_x = max(b[2] for b in all_bounds)
    max_y = max(b[3] for b in all_bounds)

    # Add 10% buffer
    x_buffer = (max_x - min_x) * 0.1 or 1.0
    y_buffer = (max_y - min_y) * 0.1 or 1.0
    extent = [min_x - x_buffer, max_x + x_buffer,
              min_y - y_buffer, max_y + y_buffer]

    # Create figure
    if use_cartopy:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
        ax.set_extent(extent, crs=ccrs.PlateCarree())

        # Add map features
        ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.5)
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.3)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.8, edgecolor='black')
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray',
                       linestyle='--')
        ax.add_feature(cfeature.LAKES, facecolor='lightblue', alpha=0.5)

        # Add gridlines
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray',
                          alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
    else:
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='--', alpha=0.5)

    # Color frames by direction
    asc_color = 'blue'
    desc_color = 'red'
    other_color = 'green'

    # Plot frames
    asc_patches = []
    desc_patches = []
    other_patches = []

    for feature in frames:
        props = feature.get('properties', {})
        geom = shape(feature['geometry'])

        direction = props.get('direction') or \
            props.get('flightDirection') or \
            props.get('flight_direction') or \
            props.get('orbit_direction') or ''

        if direction.lower() == 'ascending':
            color = asc_color
            patches_list = asc_patches
        elif direction.lower() == 'descending':
            color = desc_color
            patches_list = desc_patches
        else:
            color = other_color
            patches_list = other_patches

        # Plot the frame geometry
        if geom.geom_type == 'Polygon':
            x, y = geom.exterior.xy
            if use_cartopy:
                patch = ax.fill(x, y, alpha=0.3, facecolor=color,
                                edgecolor=color, linewidth=1.5,
                                transform=ccrs.PlateCarree())[0]
            else:
                patch = ax.fill(x, y, alpha=0.3, facecolor=color,
                                edgecolor=color, linewidth=1.5)[0]
            patches_list.append(patch)

            # Add frame ID label at centroid
            frame_id = props.get('frame_id') or props.get('frameID') or \
                props.get('id') or props.get('frame')
            if frame_id:
                centroid = geom.centroid
                if use_cartopy:
                    ax.text(centroid.x, centroid.y, str(frame_id),
                            fontsize=8, ha='center', va='center',
                            transform=ccrs.PlateCarree(),
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='white', alpha=0.7))
                else:
                    ax.text(centroid.x, centroid.y, str(frame_id),
                            fontsize=8, ha='center', va='center',
                            bbox=dict(boxstyle='round,pad=0.2',
                                      facecolor='white', alpha=0.7))

        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                x, y = poly.exterior.xy
                if use_cartopy:
                    patch = ax.fill(x, y, alpha=0.3, facecolor=color,
                                    edgecolor=color, linewidth=1.5,
                                    transform=ccrs.PlateCarree())[0]
                else:
                    patch = ax.fill(x, y, alpha=0.3, facecolor=color,
                                    edgecolor=color, linewidth=1.5)[0]
                patches_list.append(patch)

    # Plot bounding box
    bbox_patch = None
    if bbox_poly:
        x, y = bbox_poly.exterior.xy
        if use_cartopy:
            bbox_patch = ax.fill(x, y, alpha=0.1, facecolor='yellow',
                                 edgecolor='orange', linewidth=3,
                                 linestyle='-',
                                 transform=ccrs.PlateCarree())[0]
        else:
            bbox_patch = ax.fill(x, y, alpha=0.1, facecolor='yellow',
                                 edgecolor='orange', linewidth=3,
                                 linestyle='-')[0]

    # Create legend
    legend_elements = []
    if bbox_patch:
        legend_elements.append(
            mpatches.Patch(facecolor='yellow', edgecolor='orange',
                           alpha=0.3, linewidth=2, label='Bounding Box'))
    if asc_patches:
        legend_elements.append(
            mpatches.Patch(facecolor=asc_color, alpha=0.3, edgecolor=asc_color,
                           label=f'Ascending ({len(asc_patches)})'))
    if desc_patches:
        legend_elements.append(
            mpatches.Patch(facecolor=desc_color, alpha=0.3,
                           edgecolor=desc_color,
                           label=f'Descending ({len(desc_patches)})'))
    if other_patches:
        legend_elements.append(
            mpatches.Patch(facecolor=other_color, alpha=0.3,
                           edgecolor=other_color,
                           label=f'Other ({len(other_patches)})'))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right',
                  fontsize=10, framealpha=0.9)

    # Set title and labels
    ax.set_title(f'ARIA Sentinel-1 GUNW Frames\n'
                 f'({len(frames)} frames matching criteria)',
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Longitude', fontsize=12)
    ax.set_ylabel('Latitude', fontsize=12)

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'aria_frames_map.png')

    # Avoid overwriting existing files
    counter = 1
    base_path = output_path
    while os.path.exists(output_path):
        output_path = base_path.replace('.png', f'_{counter}.png')
        counter += 1

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    LOGGER.info('Saved plot to: %s', output_path)

    # Also save as PDF
    pdf_path = output_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, bbox_inches='tight')
    LOGGER.info('Saved plot to: %s', pdf_path)

    plt.close()


def main():
    """Main entry point for ariaOrderASF."""
    parser = create_parser()
    args = parser.parse_args()

    # Set up logging
    log_level = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR
    }[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    if args.verbose:
        LOGGER.setLevel(logging.DEBUG)

    # Validate inputs
    if args.bbox is None and args.frames is None:
        raise ValueError(
            'Must specify either --bbox or --frames. '
            'Use --help for usage information.'
        )

    # Validate network requirements
    if args.output == 'network' and args.network is None:
        raise ValueError(
            '--output network requires --network to be specified. '
            'Choose from: sequential, seasonal, annual.'
        )

    LOGGER.info('=' * 70)
    LOGGER.info('ARIA Sentinel-1 GUNW Frame Query / Order Tool')
    LOGGER.info('NOTE: This tool supports Sentinel-1 only, not NISAR.')
    LOGGER.info('=' * 70)

    # Fetch ARIA frames from ASF
    frames_geojson = fetch_aria_frames()

    # Parse track filter
    tracks = None
    if args.track:
        tracks = [int(t.strip()) for t in args.track.split(',')]
        LOGGER.info('Filtering by track(s): %s', tracks)

    # Parse direction filter
    direction = None
    if args.flightdir:
        direction = (
            'ascending' if args.flightdir.lower().startswith('a')
            else 'descending'
        )
        LOGGER.info('Filtering by direction: %s', direction)

    # Filter frames based on input
    bbox_poly = None
    matching_frames = []

    if args.bbox:
        bbox_poly = make_bbox(args.bbox)
        LOGGER.info(
            'Bounding box: %s',
            f'S={bbox_poly.bounds[1]:.4f}, N={bbox_poly.bounds[3]:.4f}, '
            f'W={bbox_poly.bounds[0]:.4f}, E={bbox_poly.bounds[2]:.4f}'
        )
        matching_frames = filter_frames_by_bbox(frames_geojson, bbox_poly)

    if args.frames:
        frame_ids = [int(f.strip()) for f in args.frames.split(',')]
        LOGGER.info('Querying specific frame IDs: %s', frame_ids)
        frames_by_id = filter_frames_by_ids(frames_geojson, frame_ids)
        # Merge with bbox results if both specified
        if matching_frames:
            existing_ids = {
                f.get('properties', {}).get('frame_id') or
                f.get('properties', {}).get('id')
                for f in matching_frames
            }
            for f in frames_by_id:
                fid = f.get('properties', {}).get('frame_id') or \
                      f.get('properties', {}).get('id')
                if fid not in existing_ids:
                    matching_frames.append(f)
        else:
            matching_frames = frames_by_id

    # Apply additional filters
    matching_frames = filter_frames_by_track(matching_frames, tracks)
    matching_frames = filter_frames_by_direction(matching_frames, direction)

    LOGGER.info('Found %d frames matching all criteria.', len(matching_frames))

    # Output directory
    output_dir = os.path.abspath(args.wd)
    os.makedirs(output_dir, exist_ok=True)

    # =========================================================================
    # Standard output modes: frames and/or plot
    # =========================================================================
    if args.output in ['frames', 'both']:
        print_frame_info(matching_frames)

    if args.output in ['plot', 'both']:
        if matching_frames:
            plot_frames(
                matching_frames,
                bbox_poly=bbox_poly,
                output_dir=output_dir,
                verbose=args.verbose
            )
        else:
            LOGGER.warning(
                'No frames to plot. Try adjusting your search criteria.'
            )

    # =========================================================================
    # Network generation mode
    # =========================================================================
    if args.network or args.output == 'network':
        if not matching_frames:
            raise ValueError(
                'No frames found for network generation. '
                'Adjust your search criteria.'
            )

        # Parse dates
        start_date = datetime.datetime.strptime(args.start, '%Y%m%d').date()
        if args.end:
            end_date = datetime.datetime.strptime(args.end, '%Y%m%d').date()
        else:
            end_date = datetime.date.today()

        LOGGER.info('Network generation date range: %s to %s',
                    start_date, end_date)

        # Process each frame
        for frame_feature in matching_frames:
            props = frame_feature.get('properties', {})
            frame_id = props.get('frame_id') or props.get('frameID') or \
                props.get('id') or props.get('frame')
            frame_track = props.get('track') or props.get('track_number') or \
                props.get('relativeOrbit')
            frame_dir = props.get('direction') or \
                props.get('flightDirection') or \
                props.get('flight_direction')

            LOGGER.info('\n' + '=' * 70)
            LOGGER.info('Processing Frame: %s', frame_id)
            LOGGER.info('=' * 70)

            # Get acquisitions for this frame
            dates = get_acquisitions_for_frame(
                frame_id,
                start_date,
                end_date,
                track=frame_track,
                direction=frame_dir
            )

            if len(dates) < 2:
                LOGGER.warning(
                    'Insufficient acquisitions for frame %s. '
                    'At least 2 dates required.', frame_id
                )
                continue

            LOGGER.info('Found %d acquisitions for frame %s.',
                        len(dates), frame_id)

            # Generate pairs based on network type
            network_type = args.network or 'sequential'

            if network_type == 'sequential':
                pairs = generate_sequential_pairs(
                    dates, num_neighbors=args.num_neighbors
                )
            elif network_type == 'seasonal':
                pairs = generate_seasonal_pairs(
                    dates, seasonal_window=args.seasonal_window
                )
            elif network_type == 'annual':
                pairs = generate_annual_pairs(
                    dates, num_neighbors=args.num_neighbors
                )

            # Apply temporal baseline filter if specified
            if args.temp_baseline_max:
                pairs = filter_pairs_by_temporal_baseline(
                    pairs, args.temp_baseline_max
                )

            if not pairs:
                LOGGER.warning(
                    'No pairs generated for frame %s with current settings.',
                    frame_id
                )
                continue

            # Check which products already exist
            pairs_dict = check_existing_products(pairs, frame_id)

            # Print summary
            print_network_summary(pairs_dict, dates, frame_id)

            # Plot network
            plot_network(
                dates,
                pairs_dict,
                frame_id,
                output_dir=output_dir,
                network_type=network_type
            )

    LOGGER.info('Done.')


if __name__ == '__main__':
    main()
