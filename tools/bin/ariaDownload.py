#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga, David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
import concurrent.futures
import datetime
import getpass
import logging
import math
import os
import re
import time
import warnings

import tqdm
import shapely
import asf_search
from requests.exceptions import RequestException

import ARIAtools.util.log
from ARIAtools.util.shp import open_shp
from ARIAtools.util.url import url_versions
import ARIAtools.util.s3

LOGGER = logging.getLogger('ariaDownload.py')

def createParser():
    """
    Download ARIA products using asf_search

    see: https://github.com/asfadmin/Discovery-asf_search
    """
    parser = argparse.ArgumentParser(
        description='Command line interface to download Sentinel-1/NISAR '
                    'GUNW products from the ASF DAAC. \nDownloading them '
                    'requires a NASA Earthdata URS user login',
        epilog='Examples of use:\n\n'
                '\t # Count Sentinel-1 products available for track 004\n'
                '\t ariaDownload.py --track 004 --output count\n\n'
                '\t # Download Sentinel-1 products within specified '
                'bounding box\n'
                '\t ariaDownload.py --bbox "36.75 37.225 -76.655 '
                '-75.928"\n\n'
                '\t # Count Sentinel-1 products for tracks 004 & 077 '
                'since Jan 2019\n'
                '\t ariaDownload.py --mission S1 -t 004,077 '
                '--start 20190101 -o count\n\n'
                '\t # Count all available NISAR products\n'
                '\t ariaDownload.py --mission NISAR -o count\n\n'
                '\t # Download globally available descending NISAR '
                'products\n'
                '\t ariaDownload.py --mission NISAR -d d '
                '-b "-90 90 -180 180"\n\n'
                '\t # Count all available NISAR products for track 172\n'
                '\t ariaDownload.py --mission NISAR -t 172 -o count\n\n'
                '\t # Download specific NISAR interferogram '
                'and query the globe for it\n'
                '\t ariaDownload.py --mission NISAR '
                '-b "-90 90 -180 180" -i 20251122_20251204\n',
        
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '-o', '--output', default='Download', type=str.title,
        choices=('Download', 'Count', 'Url'),
        help='Output type. Default="Download". Use "Url" for ingestion to '
             'aria*.py')
    parser.add_argument(
        '-t', '--track', default=None, type=str,
        help='track to download; single number or '
             'comma separated')
    parser.add_argument(
        '-b', '--bbox', default="-90 90 -180 180", type=str,
        help='Lat/Lon Bounding SNWE, or GDAL-readable file containing '
             'POLYGON geometry. Default is set to global scale')
    parser.add_argument(
        '-w', '--workdir', dest='wd', default='./products', type=str,
        help='Specify directory to deposit all outputs. Default is "products" '
             'in local directory where script is launched.')
    parser.add_argument(
        '-s', '--start', default='20100101', type=str,
        help='Start date as YYYYMMDD; If none provided, starts at beginning '
             'of 2010.')
    parser.add_argument(
        '-e', '--end', default='21000101', type=str,
        help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument(
        '-u', '--user', default=None, type=str,
        help='NASA Earthdata URS user login.')
    parser.add_argument(
        '-p', '--pass', dest='passw', default=None, type=str,
        help='NASA Earthdata URS user password.')
    parser.add_argument(
        '--mission', default='S1', type=str.upper, choices=('S1', 'NISAR'),
        help='Sentinel-1 (S1) or NISAR. Default is S1')
    parser.add_argument(
        '-l', '--daysless', dest='dayslt', default=math.inf, type=int,
        help='Take pairs with a temporal baseline -- days less than this '
             'value.')
    parser.add_argument(
        '-m', '--daysmore', dest='daysgt', default=0, type=int,
        help='Take pairs with a temporal baseline -- days greater than this '
             'value. Example, annual pairs: ariaDownload.py -t 004 '
             '--daysmore 364.')
    parser.add_argument(
        '-nt', '--num_threads', default='1', type=str,
        help='Specify number of threads for multiprocessing download. By '
             'default "1". Can also specify "All" to use all available '
             'threads.')
    parser.add_argument(
        '-i', '--ifg', default=None, type=str,
        help='Retrieve one interferogram by its start/end date, specified as '
             'YYYYMMDD_YYYYMMDD (order independent).')
    parser.add_argument(
        '-d', '--direction', dest='flightdir', default=None, type=str,
        help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument(
        '--version', default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods. All products '
             'are downloaded by default. If version is specified, only '
             'products which match that version are downloaded. '
             'Not supported for NISAR currently.')
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Print products to be downloaded to stdout')
    parser.add_argument(
        '--log-level', 
        choices=['debug', 'info', 'warning', 'error'], 
        default='info', 
        help='Logger log level. Default: info.'
    )
    return parser


def make_bbox(inp_bbox):
    """Make a WKT from SNWE or a shapefile"""
    if inp_bbox is None:
        return None

    if os.path.exists(os.path.abspath(inp_bbox)):
        ring = open_shp(inp_bbox, 0, 0).exterior
        poly = shapely.geometry.Polygon(ring)

    else:
        try:
            S, N, W, E = [float(i) for i in inp_bbox.split()]

            # adjust for degrees easting / northing (0 - 360 / 0:180)
            if W > 180:
                W -= 360
                LOGGER.info('AdjustedW')

            if E > 180:
                E -= 360
                LOGGER.info('AdjustedE')

            if N > 90:
                N -= 90
                S -= 90
                LOGGER.info('Adjusted N/S')

            # set poly object
            poly = shapely.geometry.Polygon([(W, N), (W, S), (E, S), (E, N)])

        except BaseException:
            raise Exception(
                'Cannot understand the --bbox argument. Input string was '
                'entered incorrectly or path does not exist.')

    return poly


def _get_s3_data_url(scene):
    """Extract the S3 data URL for the main product file from a scene.

    ASF search results include ``s3Urls`` with multiple files (browse,
    metadata, QA, etc.).  This returns the S3 URL matching the
    product's primary data file (same filename as the HTTPS URL).

    Returns
    -------
    str or None
        ``s3://bucket/path/file`` or ``None`` if not available.
    """
    props = scene.geojson()['properties']
    filename = props.get('fileName', '')
    s3_urls = props.get('s3Urls', [])
    for s3_url in s3_urls:
        if s3_url.endswith(filename):
            return s3_url
    return None


def get_url_ifg(scenes):
    """Get url, ifg of fetched ASF scene"""
    urls, ifgs = [], []
    for scene in scenes:
        s = scene.geojson()['properties']
        urls.append(s['url'])
        # NISAR files are formatted differently
        if s['fileID'].startswith('NISAR_'):
            f = s['fileID'].split('_')
            pairname = f[11][:8] + '_'
            pairname += f[13][:8]
            ifgs.append(pairname)
        else:
            f = s['fileID'].split('-')
            pairname = f[6]
            ifgs.append(pairname)

    # determine if NISAR GUNW
    is_nisar_file = False
    if urls != []:
        if '/NISAR_' in urls[0]:
            is_nisar_file = True

    return urls, ifgs, is_nisar_file


def fmt_dst(args):
    """Format the save name"""
    ext = '.kmz' if args.output == 'Kml' else '.txt'

    if args.track is not None:
        fn_track = f'track{args.track}'.replace(',', '-')
    else:
        fn_track = ''

    if args.bbox is not None:
        WSEN = make_bbox(args.bbox).bounds
        WSEN_fmt = []
        for i, coord in enumerate(WSEN):
            if i < 2:
                WSEN_fmt.append(math.floor(float(coord)))
            else:
                WSEN_fmt.append(math.ceil(float(coord)))
        fn_bbox = f'_bbox{WSEN_fmt[0]}W{WSEN_fmt[1]}S{WSEN_fmt[2]}E{WSEN_fmt[3]}N'
    else:
        fn_bbox = ''

    dst = os.path.join(args.wd, f'{fn_track}{fn_bbox}_0{ext}'.lstrip('_'))
    count = 1  # don't overwrite if already exists
    while os.path.exists(dst):
        basen = f'{re.split(str(count-1)+ext, os.path.basename(dst))[0]}' \
                f'{count}{ext}'
        dst = os.path.join(os.path.dirname(dst), basen)
        count += 1
    return dst


class Downloader:
    """Product Downloading Class."""
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.args.output = self.args.output.title()
        self.args.wd = os.path.abspath(self.args.wd)
        os.makedirs(self.args.wd, exist_ok=True)
        LOGGER.setLevel(logging.DEBUG if self.args.verbose else logging.INFO)

    def __call__(self):
        scenes = self.query_asf()
        urls, ifgs, is_nisar_file = get_url_ifg(scenes)

        # Subset everything by version
        if is_nisar_file and self.args.version is not None:
            raise Exception(
                'Version support not included for NISAR, remove the critera'
            )
        else:
            urls = url_versions(urls, self.args.version, self.args.wd)
        scenes = [scene for scene, url in zip(scenes, urls) if url in urls]
        ifgs = [ifg for ifg, url in zip(ifgs, urls) if url in urls]

        # Filter scenes based on date and elapsed time criteria
        scenes, urls, ifgs = self.filter_scenes(
            scenes,
            urls,
            ifgs,
            is_nisar_file
        )
        

        if self.args.output == "Count":
            LOGGER.info("Found -- %d -- products", len(scenes))
        elif self.args.output == "Url":
            self.write_urls(urls, scenes)
        elif self.args.output == "Download":
            self.download_scenes(scenes)
            
        if self.args.verbose:
            for scene in scenes:
                LOGGER.info(scene.geojson()["properties"]["sceneName"])

    def query_asf(self):
        """Query ASF for scenes."""
        bbox = make_bbox(self.args.bbox)
        bbox_wkt = bbox.wkt if bbox else None

        flight_direction = None
        if self.args.flightdir:
            flight_direction = (
                "ascending"
                if self.args.flightdir.lower().startswith("a")
                else "descending"
            )

        tracks = (
            [int(track) for track in self.args.track.split(",")]
            if self.args.track
            else None
        )

        start = self.args.start - datetime.timedelta(days=1)
        end = self.args.end + datetime.timedelta(days=1)

        if self.args.mission.upper() == "S1":
            return asf_search.geo_search(
                collections=["C2859376221-ASF", "C1261881077-ASF"],
                dataset=asf_search.constants.ARIA_S1_GUNW,
                processingLevel=asf_search.constants.GUNW_STD,
                relativeOrbit=tracks,
                flightDirection=flight_direction,
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
            )
        elif self.args.mission.upper() == "NISAR":
            # Authenticate so the private ephemeral archive
            # collection (C4052499921-ASF) is visible in CMR.
            session = self._get_asf_session()
            opts = asf_search.ASFSearchOptions(
                collections=[
                    "C2850261892-ASF",   # public NISAR GUNW
                    "C4052499921-ASF",   # private ephemeral archive
                ],
                dataset=asf_search.constants.NISAR,
                processingLevel=asf_search.constants.GUNW,
                relativeOrbit=tracks,
                flightDirection=flight_direction,
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
                session=session,
            )
            return asf_search.geo_search(opts=opts)

    def _get_asf_session(self):
        """Return an authenticated ASFSession.

        Uses explicit user/pass args when provided, otherwise falls
        back to ~/.netrc credentials for urs.earthdata.nasa.gov.
        Returns an unauthenticated session with a warning when no
        credentials are available.
        """
        session = asf_search.ASFSession()
        if self.args.user:
            session.auth_with_creds(
                self.args.user,
                self.args.passw or getpass.getpass("NASA Earthdata password: ")
            )
        else:
            try:
                import netrc as _netrc
                nrc = _netrc.netrc()
                auth = nrc.authenticators('urs.earthdata.nasa.gov')
                if auth:
                    session.auth_with_creds(auth[0], auth[2])
                else:
                    LOGGER.warning(
                        'No urs.earthdata.nasa.gov entry in ~/.netrc. '
                        'Private collections (e.g. NISAR ephemeral '
                        'archive) will not be visible.')
            except FileNotFoundError:
                LOGGER.warning(
                    '~/.netrc not found. Private collections (e.g. '
                    'NISAR ephemeral archive) will not be visible.')
        return session

    def filter_scenes(self, scenes, urls, ifgs, is_nisar_file):
        filtered_scenes, filtered_urls, filtered_ifgs = [], [], []

        for scene, url, ifg in zip(scenes, urls, ifgs):
            eni, sti = self.parse_dates(ifg, is_nisar_file)
            if self.args.ifg:
                if self.match_single_ifg(sti, eni):
                    filtered_scenes.append(scene)
                    filtered_urls.append(url)
                    filtered_ifgs.append(ifg)
            elif self.match_date_criteria(sti, eni):
                filtered_scenes.append(scene)
                filtered_urls.append(url)
                filtered_ifgs.append(ifg)
        return filtered_scenes, filtered_urls, filtered_ifgs

    def parse_dates(self, ifg, is_nisar_file):
        if is_nisar_file:
            sti, eni = [datetime.datetime.strptime(d, "%Y%m%d")
                        for d in ifg.split("_")]
        else:
            eni, sti = [datetime.datetime.strptime(d, "%Y%m%d")
                        for d in ifg.split("_")]
        return eni, sti

    def match_single_ifg(self, sti, eni):
        dates = [
            datetime.datetime.strptime(i, "%Y%m%d").date()
            for i in self.args.ifg.split("_")
        ]
        st1, en1 = sorted(dates)
        return st1 == sti.date() and en1 == eni.date()

    def match_date_criteria(self, sti, eni):
        sten_chk = sti >= self.args.start and eni <= self.args.end
        elap = (eni - sti).days
        elap_chk = self.args.daysgt <= elap <= self.args.dayslt
        return sten_chk and elap_chk

    def write_urls(self, urls, scenes):
        dst = fmt_dst(self.args)
        with open(dst, "w") as fh:
            for url, scene in zip(urls, scenes):
                s3_url = _get_s3_data_url(scene) or ''
                print(f'{url},{s3_url}', file=fh)
        LOGGER.info("Wrote -- %d -- product urls to: %s", len(urls), dst)

    def download_scenes(self, scenes):
        scenes = asf_search.ASFSearchResults(scenes)
        nt = int(self.args.num_threads)
        LOGGER.info("Downloading %d products...", len(scenes))

        # Check if we can use S3 direct download (on AWS)
        use_s3 = ARIAtools.util.s3.is_on_aws()
        s3_client = None
        if use_s3:
            # Determine credential endpoint from first S3 URL
            first_s3 = next(
                (_get_s3_data_url(s) for s in scenes
                 if _get_s3_data_url(s)), None)
            endpoint_key = ARIAtools.util.s3._endpoint_key_for_s3uri(
                first_s3) if first_s3 else 'default'
            try:
                s3_client = ARIAtools.util.s3.get_s3_client(
                    endpoint_key,
                    max_pool_connections=nt * 10)
                LOGGER.info('Using S3 direct download (endpoint: %s)',
                            endpoint_key)
            except Exception as exc:
                LOGGER.warning('S3 client setup failed, falling back '
                               'to HTTPS: %s', exc)
                use_s3 = False

        session = asf_search.ASFSession()
        if self.args.user:
            session.auth_with_creds(self.args.user, self.args.passw)

        def download_file(scene, max_retries=3, retry_delay=5):
            url = scene.properties['url']
            local_filename = url.split("/")[-1]
            filepath = os.path.join(self.args.wd, local_filename)

            # Skip download if file already exists
            if os.path.exists(filepath):
                pbar.update(1)
                LOGGER.info("Product already in directory: %s", filepath)
                return filepath

            # Try S3 download first if available
            s3_url = _get_s3_data_url(scene) if use_s3 else None
            if s3_client and s3_url:
                attempt = 0
                while attempt < max_retries:
                    attempt += 1
                    try:
                        bucket, key = \
                            ARIAtools.util.s3.parse_s3_uri(s3_url)
                        s3_client.download_file(
                            bucket, key, filepath)
                        LOGGER.debug('S3 download: %s', filepath)
                        return filepath
                    except Exception as exc:
                        LOGGER.warning(
                            'S3 download attempt %d failed: %s',
                            attempt, exc)
                        if os.path.exists(filepath):
                            os.remove(filepath)
                        if attempt < max_retries:
                            time.sleep(retry_delay)
                LOGGER.warning('S3 download failed, falling back '
                               'to HTTPS for %s', local_filename)

            # HTTPS download (default or fallback)
            attempt = 0
            while attempt < max_retries:
                attempt += 1
                try:
                    response = session.get(url, stream=True)
                    response.raise_for_status()

                    with open(filepath, "wb") as f:
                        for chunk in response.iter_content(
                                chunk_size=8192):
                            if chunk:
                                f.write(chunk)

                    # Verify download size
                    expected_size = int(
                        response.headers.get("Content-Length", 0)
                    )
                    file_size = os.path.getsize(filepath)
                    if expected_size > 0 and file_size < expected_size:
                        LOGGER.warning(
                            "Incomplete download detected "
                            "(%d/%d bytes). Retrying...",
                            file_size, expected_size
                        )
                        os.remove(filepath)
                        time.sleep(retry_delay)
                        continue

                except RequestException as e:
                    LOGGER.error("Error downloading %s: %s",
                                 url, e)

            return filepath

        # Create a progress bar
        pbar = tqdm.tqdm(total=len(scenes), unit="file",
                         desc="Downloading")
        try:
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=nt) as executor:
                future_to_scene = {
                    executor.submit(download_file, scene): scene
                    for scene in scenes
                }
                for future in concurrent.futures.as_completed(
                        future_to_scene):
                    scene = future_to_scene[future]
                    try:
                        filepath = future.result()
                        LOGGER.debug("Downloaded: %s", filepath)
                        pbar.update(1)
                    except Exception as exc:
                        LOGGER.error(
                            "%s generated an exception: %s",
                            scene.properties['url'], exc)
        finally:
            pbar.close()

        LOGGER.info(
            "Download complete. Wrote -- %d -- products to: %s",
            len(scenes),
            self.args.wd,
        )


def main():
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    print('*****************************************************************')
    LOGGER.info('*** Download Function ***')
    print('*****************************************************************')

    # format dates
    args.start = datetime.datetime.strptime(args.start, '%Y%m%d')
    args.end = datetime.datetime.strptime(args.end, '%Y%m%d')

    if not args.track and not args.bbox:
        raise Exception('Must specify either a bbox or track')
    Downloader(args)()


if __name__ == '__main__':
    main()
