#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import re
import math
import datetime
import argparse
import logging
import shapely
import warnings
from getpass import getpass
import asf_search

import ARIAtools.util.log
from ARIAtools.util.shp import open_shp
from ARIAtools.util.url import url_versions
from pkg_resources import get_distribution

LOGGER = logging.getLogger('ariaDownload.py')

def createParser():
    """
    Download ARIA products using asf_search

    see: https://github.com/asfadmin/Discovery-asf_search
    """
    parser = argparse.ArgumentParser(
        description='Command line interface to download Sentinel-1/NISAR GUNW '
                    'products from the ASF DAAC. \nDownloading them requires a '
                     'NASA Earthdata URS user login and requires users to add '
                     '"GRFN Door (PROD)" and "ASF Datapool Products" to their '
                     'URS approved applications. Access to NISAR products '
                     'requires an Earthdata Bearer token from: '
                     'https://urs.earthdata.nasa.gov/documentation/for_users/user_token',
        epilog='Examples of use:\n'
                '\t ariaDownload.py --track 004 --output count\n'
                '\t ariaDownload.py --bbox "36.75 37.225 -76.655 -75.928"\n'
                '\t ariaDownload.py -t 004,077 --start 20190101 -o count'
                '\t ariaDownload.py --mission NISAR -o count',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        '-o', '--output', default='Download', type=str.title,
        choices=('Download', 'Count', 'Url'),
        help='Output type. Default="Download". Use "Url" for ingestion to '
             'aria*.py')
    parser.add_argument(
        '-t', '--track', default=None, type=str,
        help='track to download; single number (including leading zeros) or '
             'comma separated')
    parser.add_argument(
        '-b', '--bbox', default=None, type=str,
        help='Lat/Lon Bounding SNWE, or GDAL-readable file containing POLYGON '
             'geometry.')
    parser.add_argument(
        '-w', '--workdir', dest='wd', default='./products', type=str,
        help='Specify directory to deposit all outputs. Default is "products" '
             'in local directory where script is launched.')
    parser.add_argument(
        '-s', '--start', default='20140101', type=str,
        help='Start date as YYYYMMDD; If none provided, starts at beginning '
             'of Sentinel record (2014).')
    parser.add_argument(
        '-e', '--end', default='21000101', type=str,
        help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument(
        '-u', '--user', default=None, type=str,
        help='NASA Earthdata URS user login. Users must add "GRFN Door (PROD)" '
             'and "ASF Datapool Products" to their URS approved applications.')
    parser.add_argument(
        '-p', '--pass', dest='passw', default=None, type=str,
        help='NASA Earthdata URS user password. Users must add "GRFN Door '
             '(PROD)" and "ASF Datapool Products" to their URS approved '
             'applications.')

    parser.add_argument(
        '--mission', default='S1', type=str.upper, choices=('S1', 'NISAR'),
        help='Sentinel-1 (S1) or NISAR. Default is S1')

    parser.add_argument(
        '-l', '--daysless', dest='dayslt', default=math.inf, type=int,
        help='Take pairs with a temporal baseline -- days less than this value.')
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
             'YYYYMMDD_YYYYMMDD (order independent)')
    parser.add_argument(
        '-d', '--direction', dest='flightdir', default=None, type=str,
        help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument(
        '--version', default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods. All products '
             'are downloaded by default. If version is specified, only '
             'products which match that version are downloaded.')
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Print products to be downloaded to stdout')
    parser.add_argument(
        '--log-level', default='warning', help='Logger log level')
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

        except BaseException:
            raise Exception('Cannot understand the --bbox argument. '
                            'Input string was entered incorrectly or path does not '
                            'exist.')

    return shapely.geometry.Polygon([(W, N), (W, S), (E, S), (E, N)])


def get_url_ifg(scenes):
    """Get url, ifg of fetched ASF scene"""
    urls, ifgs = [], []
    for scene in scenes:
        s = scene.geojson()['properties']
        f = s['fileID'].split('-')
        ifgs.append(f[6])
        urls.append(s['url'])

    return urls, ifgs


def fmt_dst(inps):
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


class Downloader(object):
    def __init__(self, args):
        self.args = args
        self.args.output = self.args.output.title()
        self.args.wd = os.path.abspath(self.args.wd)
        os.makedirs(self.args.wd, exist_ok=True)
        if self.args.verbose:
            LOGGER.setLevel('DEBUG')
        else:
            LOGGER.setLevel('INFO')


    def __call__(self):
        scenes = self.query_asf()

        urls, ifgs = get_url_ifg(scenes)

        # subset everything by version
        urls_new = url_versions(urls, self.args.version, self.args.wd)
        idx = [urls.index(url) for url in urls_new]
        scenes = [scenes[i] for i in idx]
        urls = [urls[i] for i in idx]
        ifgs = [ifgs[i] for i in idx]

        # loop over ifgs to check for subset options
        idx = []
        for i, ifg in enumerate(ifgs):
            eni, sti = [datetime.datetime.strptime(
                d, '%Y%m%d') for d in ifg.split('_')]

            # optionally match a single ifg
            if self.args.ifg is not None:
                dates1 = [datetime.datetime.strptime(i, '%Y%m%d').date() for
                          i in self.args.ifg.split('_')]
                st1, en1 = sorted(dates1)
                if st1 == sti.date() and en1 == eni.date():
                    idx.append(i)
                else:
                    continue

            # optionally match other conditions (st/end date, elapsed time)
            else:
                sten_chk = sti >= self.args.start and eni <= self.args.end
                elap = (eni - sti).days
                elap_chk = elap >= self.args.daysgt and elap <= self.args.dayslt
                if sten_chk and elap_chk:
                    idx.append(i)
                else:
                    continue

        scenes = [scenes[i] for i in idx]
        urls = [urls[i] for i in idx]
        ifgs = [ifgs[i] for i in idx]

        if self.args.output == 'Count':
            LOGGER.info('Found -- %d -- products', len(scenes))

        # elif self.args.output == 'Kml':
        #     dst    = fmt_dst(inps)
        #     LOGGER.error('Kml option is not yet supported. '\
        #                    'Revert to an older version of ARIAtools')

        elif self.args.output == 'Url':
            dst = fmt_dst(inps)
            with open(dst, 'w') as fh:
                for url in urls:
                    print(url, sep='\n', file=fh)
            LOGGER.info(f'Wrote -- {len(urls)} -- product urls to: {dst}')

        elif self.args.output == 'Download':
            # turn the list back into an ASF object
            scenes = asf_search.ASFSearchResults(scenes)
            nt = int(self.args.num_threads)  # so legacy works
            # allow a user to specify username / password
            LOGGER.info(f'Downloading {len(scenes)} products...')
            if self.args.user is not None:
                session = asf_search.ASFSession()
                session.auth_with_creds(self.args.user, self.args.passw)

                # Suppress warnings on files already downloaded
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    scenes.download(self.args.wd, processes=nt, session=session)

            else:
                # Suppress warnings on files already downloaded
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    scenes.download(self.args.wd, processes=nt)

            LOGGER.info(
                f'Download complete. Wrote -- {len(scenes)} -- products to: '
                '{self.args.wd}')

        if self.args.verbose:
            for scene in scenes:
                LOGGER.info(scene.geojson()['properties']['sceneName'])

        return


    def query_asf(self):
        """Get the scenes from ASF"""
        bbox = make_bbox(self.args.bbox)
        bbox = bbox.wkt if bbox is not None else None
        if self.args.flightdir is not None:
            if self.args.flightdir.lower()[0] == 'a':
                flight_direction = 'ascending'
            elif self.args.flightdir == 'd':
                self.args.flightdir.lower()[0] = 'descending'
        else:
            flight_direction = None

        if self.args.track is not None:
            tracks = self.args.track.split(',')
            tracks = [int(track) for track in tracks]
        else:
            tracks = self.args.track


        ## buffer a bit for loose subsetting and speedup
        start = self.args.start + datetime.timedelta(days=-90)
        end   = self.args.end + relativedelta(days=90)

        if self.args.mission.upper() == 'S1':
            dct_kw = dict(dataset='ARIA S1 GUNW',
                            processingLevel=asf_search.constants.GUNW_STD,
                            relativeOrbit=tracks,
                            flightDirection=flight_direction,
                            intersectsWith=bbox,
                            start=start,
                            end=end)
            scenes = asf_search.geo_search(**dct_kw)

        elif self.args.mission.upper() == 'NISAR':
            session = asf_search.ASFSession()
            session.auth_with_token(getpass('EDL Token:'))
            LOGGER.info('Token accepted.')
            ## populate once available on GUNWs are available
            # TODO:
            search_opts = asf_search.ASFSearchOptions(
                            shortName='NISAR_L2_GUNW_BETA_V1',
                            session=session
                            # processingLevel=asf_search.constants.GUNW_STD,
                            # relativeOrbit=tracks,
                            # flightDirection=flight_direction,
                            # intersectsWith=bbox,
                            # start=self.args.start,
                            # end=self.args.end)
                        )
            scenes = asf_search.search(opts=search_opts, maxResults=250)
            if len(scenes)>0:
                LOGGER.info('Found NISAR GUNW Betas.')
            else:
                LOGGER.warning('No NISAR GUNW Betas found.')
            raise Exception('Exiting, NISAR GUNWs not futher supported.')

        return scenes


def main():
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    # format dates
    args.start = datetime.datetime.strptime(args.start, '%Y%m%d')
    args.end = datetime.datetime.strptime(args.end, '%Y%m%d')

    if not args.track and not args.bbox:
        raise Exception('Must specify either a bbox or track')
    Downloader(args)()


if __name__ == '__main__':
    main()
