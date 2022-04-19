#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os, os.path as op
import argparse
import math
import re
from datetime import datetime
import logging
import asf_search as asf
from ARIAtools.logger import logger
from ARIAtools.url_manager import url_versions
from pkg_resources import get_distribution

log = logging.getLogger('ARIAtools')


def createParser():
    """Download ARIA products using asf_search

    see: https://github.com/asfadmin/Discovery-asf_search
    """
    parser = argparse.ArgumentParser(description=
        'Command line interface to download GUNW products from the ASF DAAC. '
        'GUNW products are hosted at the NASA ASF DAAC.\nDownloading them '
        'requires a NASA Earthdata URS user login and requires users to add '
        '"GRFN Door (PROD)" and "ASF Datapool Products" to their URS '
        'approved applications.',
        epilog='Examples of use:\n\t ariaDownload.py --track 004 '
                '--output count'
                '\n\t ariaDownload.py --bbox "36.75 37.225 -76.655 -75.928"'
                '\n\t ariaDownload.py -t 004,077 --start 20190101 -o count',
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', \
                        type=str,
        help='Output type, default is "Download". "Download", "Count", and "Url"'
             '"Kmz" are currently supported. Use "Url" for ingestion to '
             'aria*.py')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str,
        help='track to download; single number (including leading zeros) or '
             'comma separated')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str,
        help='Lat/Lon Bounding SNWE, or GDAL-readable file containing '
             'POLYGON geometry.')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', \
                        type=str,
        help='Specify directory to deposit all outputs. Default is '
             '"products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default='20140101', type=str,
        help='Start date as YYYYMMDD; If none provided, starts at beginning '
             'of Sentinel record (2014).')
    parser.add_argument('-e', '--end', dest='end', default='21000101', type=str,
        help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument('-u', '--user', dest='user', default=None, type=str,
        help='NASA Earthdata URS user login. Users must add "GRFN Door '
             '(PROD)" and "ASF Datapool Products" to their URS approved '
             'applications.')
    parser.add_argument('-p', '--pass', dest='passw', default=None, type=str,
        help='NASA Earthdata URS user password. Users must add "GRFN Door '
             '(PROD)" and "ASF Datapool Products" to their URS approved '
             'applications.')
    parser.add_argument('-l', '--daysless', dest='dayslt', default=math.inf, \
                        type=int,
        help='Take pairs with a temporal baseline -- days less than this '
             'value.')
    parser.add_argument('-m', '--daysmore', dest='daysgt', default=0, \
                        type=int,
        help='Take pairs with a temporal baseline -- days greater than this '
             'value. Example, annual pairs: ariaDownload.py -t 004 '
             '--daysmore 364.')
    parser.add_argument('-nt', '--num_threads', dest='num_threads', \
                        default='1', type=str,
        help='Specify number of threads for multiprocessing '
             'download. By default "1". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument('-i', '--ifg', dest='ifg', default=None, type=str,
        help='Retrieve one interferogram by its start/end date, specified as '
             'YYYYMMDD_YYYYMMDD (order independent)')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, \
                        type=str,
        help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('--version', dest='version',  default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods; default: '
             'newest')
    parser.add_argument('-v', '--verbose', dest='v', action='store_true',
        help='Print products to be downloaded to stdout')
    return parser


def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(2)

    inps = parser.parse_args(args=iargs)

    ## format dates
    inps.start = datetime.strptime(inps.start, '%Y%m%d')
    inps.end   = datetime.strptime(inps.end, '%Y%m%d')

    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    if not inps.output.lower() in ['count', 'kmz', 'kml', 'url', 'download']:
        raise Exception ('Incorrect output keyword. Choose "count", "kmz", '
                         '"url", or "download"')

    inps.output = 'Kml' if inps.output.lower() == 'kmz' or \
            inps.output.lower() == 'kml' else inps.output.title()
    return inps


def make_bbox(inp_bbox):
    """Make a WKT from SNWE or a shapefile"""
    if inp_bbox is None:
        return None
    from shapely.geometry import Polygon
    if op.exists(op.abspath(inps.bbox)):
        from ARIAtools.shapefile_util import open_shapefile
        ring = open_shapefile(inps.bbox, 0, 0).exterior
        poly = Polygon(ring)
    else:
        try:
            S, N, W, E = [float(i) for i in inps.bbox.split()]
            ## adjust for degrees easting / northing (0 - 360 / 0:180)
            if W > 180: W -= 360; print('AdjustedW')
            if E > 180: E -= 360; print('AdjustedE')
            if N > 90: N-=90; S-=90; print('Adjusted N/S')
        except:
            raise Exception('Cannot understand the --bbox argument. '
            'Input string was entered incorrectly or path does not '
            'exist.')

        poly = Polygon([(W, N), (W,S), (E,S), (E, N)])
    return poly


def get_url_ifg(scenes):
    """Get url, ifg of fetched ASF scene"""
    urls, ifgs =  [], []
    for scene in scenes:
        s   = scene.geojson()['properties']
        f   = s['fileID'].split('-')
        ifgs.append(f[6])
        urls.append(s['url'])

    return urls, ifgs


def fmt_dst(inps):
    """Format the save name"""
    ext  = '.kmz' if inps.output == 'Kml' else '.txt'

    if inps.track is not None:
        fn_track = f'track{inps.track}'.replace(',', '-')
    else:
        fn_track = ''

    if inps.bbox is not None:
        WSEN = make_bbox(inps.bbox).bounds
        WSEN_fmt = []
        for i, coord in enumerate(WSEN):
            if i < 2:
                WSEN_fmt.append(math.floor(float(coord)))
            else:
                WSEN_fmt.append(math.ceil(float(coord)))
        fn_bbox = f'_bbox{WSEN_fmt[0]}W{WSEN_fmt[1]}S{WSEN_fmt[2]}E{WSEN_fmt[3]}N'
    else:
        fn_bbox = ''

    dst   = op.join(inps.wd, f'{fn_track}{fn_bbox}_0{ext}'.lstrip('_'))
    count = 1 # don't overwrite if already exists
    while op.exists(dst):
        basen  = f'{re.split(str(count-1)+ext, op.basename(dst))[0]}' \
                 f'{count}{ext}'
        dst    = op.join(op.dirname(dst), basen)
        count += 1
    return dst


class Downloader(object):
    def __init__(self, inps):
        self.inps            = inps
        self.inps.output     = self.inps.output.title()
        self.inps.wd         = op.abspath(self.inps.wd)
        os.makedirs(self.inps.wd, exist_ok=True)


    def __call__(self):
        scenes     = self.query_asf()
        urls, ifgs = get_url_ifg(scenes)

        ## subset everything by version
        urls_new = url_versions(urls, self.inps.version, self.inps.wd)
        idx      = [urls.index(url) for url in urls_new]
        scenes   = [scenes[i] for i in idx]
        urls     = [urls[i] for i in idx]
        ifgs     = [ifgs[i] for i in idx]

        ## loop over ifgs to check for subset options
        idx = []
        for i, ifg in enumerate(ifgs):
            eni, sti = [datetime.strptime(d, '%Y%m%d') for d in ifg.split('_')]
            ## optionally match a single ifg
            if self.inps.ifg is not None:
                dates1    = [datetime.strptime(i, '%Y%m%d').date() for
                                    i in self.inps.ifg.split('_')]
                st1, en1 = sorted(dates1)
                if st1 == sti.date() and en1 == eni.date():
                    idx.append(i)
                else:
                    continue

            ## optionally match other conditions (st/end date, elapsed time)
            else:
                sten_chk = sti >= self.inps.start and  eni <= self.inps.end
                elap     = (eni-sti).days
                elap_chk =  elap >= self.inps.daysgt and elap <= self.inps.dayslt
                if sten_chk and elap_chk:
                    idx.append(i)
                else:
                    continue

        scenes   = [scenes[i] for i in idx]
        urls     = [urls[i] for i in idx]
        ifgs     = [ifgs[i] for i in idx]

        if self.inps.output == 'Count':
            log.info('\nFound -- %d -- products', len(scenes))

        elif self.inps.output == 'Kml':
            dst    = self._fmt_dst()
            self.log.error('Kml option is not yet supported. '\
                           'Revert to an older version of ARIAtools')

        elif self.inps.output == 'Url':
            dst  = fmt_dst(inps)
            with open(dst, 'w') as fh: [print(url, sep='\n', file=fh) \
                                        for url in urls]
            log.info(f'Wrote -- {len(urls)} -- product urls to: {dst}')

        elif self.inps.output == 'Download':
            ## turn the list back into an ASF object
            scenes = asf.ASFSearchResults(scenes)
            nt     = int(self.inps.num_threads) # so legacy works
            ## allow a user to specify username / password
            if self.inps.user is not None:
                session = asf.ASFSession()
                session.auth_with_creds(self.inps.user, self.inps.passw)
                scenes.download(self.inps.wd, processes=nt, session=session)

            else:
                scenes.download(self.inps.wd, processes=nt)
            log.info(f'Wrote -- {len(scenes)} -- products to: {self.inps.wd}')

        return


    def query_asf(self):
        """Get the scenes from ASF"""
        bbox   = make_bbox(self.inps.bbox)

        if self.inps.track is not None:
            tracks = self.inps.track.split(',')
            tracks = [int(track) for track in tracks]
        else:
            tracks = self.inps.track
        dct_kw = dict(platform=asf.constants.SENTINEL1,
                    processingLevel=asf.constants.GUNW_STD,
                    relativeOrbit=tracks,
                    lookDirection=self.inps.flightdir,
                    intersectsWith=bbox)
        scenes = asf.geo_search(**dct_kw)


        return scenes


if __name__ == '__main__':
    try:
        print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except:
        pass

    print ('this is a test')
    inps = cmdLineParse()
    Downloader(inps)()
