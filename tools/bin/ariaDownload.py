#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga, Marin Govorcin
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
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
import earthaccess
import numpy as np
import pandas as pd
from shapely.geometry import Polygon
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
    parser.add_argument('-o', '--output', default='Download', type=str.title,
                 choices=('Download', 'Count', 'Url'), help='Output type. '\
                     'Default="Download". Use "Url" for ingestion to aria*.py')
    parser.add_argument('-t', '--track', default=None, type=str,
        help='track to download; single number (including leading zeros) or '
             'comma separated')
    parser.add_argument('-b', '--bbox',  default=None, type=str,
        help='Lat/Lon Bounding SNWE, or GDAL-readable file containing '
             'POLYGON geometry.')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', \
                        type=str,
        help='Specify directory to deposit all outputs. Default is '
             '"products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', default='20140101', type=str,
        help='Start date as YYYYMMDD; If none provided, starts at beginning '
             'of Sentinel record (2014).')
    parser.add_argument('-e', '--end', default='21000101', type=str,
        help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument('-u', '--user', default=None, type=str,
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
    parser.add_argument('-nt', '--num_threads', default='1', type=str,
        help='Specify number of threads for multiprocessing '
             'download. By default "1". Can also specify "All" to use all '
             'available threads.')
    parser.add_argument('-i', '--ifg', default=None, type=str,
        help='Retrieve one interferogram by its start/end date, specified as '
             'YYYYMMDD_YYYYMMDD (order independent)')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, \
                        type=str,
        help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('-api', '--api',  dest='api', default='asf_search',
                choices=('asf_search', 'earthaccess'), help='Use different api to'\
                    'download data either from ASF Vertex or NASA Earthdata')
    parser.add_argument('--version',  default=None,
        help='Specify version as str, e.g. 2_0_4 or all prods.'
             'All products are downloaded by default. If version is specified, '
             'only products which match that version are downloaded.')
    parser.add_argument('-v', '--verbose', action='store_true',
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

# Earthdata query and download
# MG added, April 04, 2024,
def query_earthaccess(bbox,
                      flight_direction : str = None, tracks_id :list = None,
                      start_date : str = None, end_date : str = None,
                      min_btemp : int = None, max_btemp :int = None,
                      gunw_version : str = None):
    '''
    Parameters
        bbox : tuple                bounding box (West, South, North East)
        flight_direction : str      flight orbit direction, ascending or descending 
        track : list                list of path number to filter the results
        
        Filtering options:
        start_date : str            start date in format (YYYY-MM-DD)
        end_date : str              end date in format (YYYY-MM-DD)
        min_btemp : int             min temporal baseline in days
        max_btemp : int             max temporal baseline in days
        gunw_version : str          gunw version, e.g. 2_0_4          
  
    Return
        earthdata  : object         list of earthdata granules
        earthdata_df : Dataframe    dataframe of selected granules containing filenames
    '''

    def _get_scenes_id(earthaccess_results) -> list:
        return list(map(lambda x: x.items().mapping['meta']['native-id'],
                    earthaccess_results))
    
    def _get_scenes_bounds(earthaccess_results) -> list:
        ''' Get products bounds '''
        def _get_polygon(geom_dict):
            coords = geom_dict['HorizontalSpatialDomain']['Geometry']['GPolygons'][0]['Boundary']['Points']
            return Polygon([(point['Longitude'], point['Latitude']) for point in coords])

        geom =  list(map(lambda x: x.items().mapping['umm']['SpatialExtent'],
                        earthaccess_results))
        bounds = list(map(_get_polygon, geom))
        return bounds
    
    def _get_scenes_attributes(earthaccess_results) -> list:
        ''' Get additional attributes '''
        def _get_attribute(attr_dict):
            attr = [{key['Name']:key['Values'][0]} for key in attr_dict.items().mapping['umm']['AdditionalAttributes']]
            attr =  {k: v for d in attr for k, v in d.items()}
            return attr

        return list(map(lambda x: _get_attribute(x), earthaccess_results))
    
    def _get_scenes_url(earthaccess_results) -> list:
        ''' Get additional urls '''
        url = list(map(lambda x: x.items().mapping['umm']['RelatedUrls'][0]['URL'],
                       earthaccess_results))
        url_png = list(map(lambda x: x.items().mapping['umm']['RelatedUrls'][1]['URL'],
                           earthaccess_results))

        return url, url_png
    
    if bbox is None:
        # NOTE: this can add some time to query
        bbox = (-180,-90, 180,90)

    # Get the scenes from Earthdata
    print('Searching ARIA-GUNW on NASA Earthdata:')
    results = earthaccess.search_data(
                short_name='ARIA_S1_GUNW',
                bounding_box=bbox,
                temporal=(start_date, end_date),
                count=-1
            )
    
    # Get results scenes id
    scenes_id = _get_scenes_id(results)
    scenes_df = pd.DataFrame(list(map(lambda x: x.split('-'), scenes_id)),
                            columns=['sensor', 'datasetName', 'flight_direction',
                                    'look_direction', 'track_number', 'mode', 
                                    'reference_secondary', 'UTC_time', 'upper_left_latlon',
                                    'orbits', 'system_tag', 'version'])

    # Change track number dtype from str to int
    scenes_df = scenes_df.astype({'track_number':int})

    # Add bounds, attributes, urls
    scenes_bounds = _get_scenes_bounds(results)
    scenes_attributes = _get_scenes_attributes(results)
    scenes_df['geometry'] = scenes_bounds
    scenes_df = scenes_df.join(pd.DataFrame(scenes_attributes))
    url, url_png =_get_scenes_url(results)
    scenes_df['URL'] = url
    scenes_df['PNG_URL'] = url_png

    # split reference and secondary date and covert to datetime
    # needed for btemp filtering
    str2datetime = lambda x: datetime.strptime(x, '%Y%m%d')
    date12 = np.vstack(scenes_df.reference_secondary.apply(lambda x: x.split('_')))
    scenes_df['start'] = list(map(str2datetime, date12[:, 0]))
    scenes_df['end'] = list(map(str2datetime, date12[:, 1]))
    scenes_df['btemp'] = (scenes_df.start - scenes_df.end).dt.days

    # filter based on flight direction
    if flight_direction is not None: 
        if flight_direction[0].upper() in ['A', 'D']:
            print(f'Filter based on flight direction: {flight_direction}')
            flight_direction = flight_direction[0].upper()
            scenes_df = scenes_df[scenes_df.flight_direction == flight_direction]

    # filter based on track number
    if tracks_id is not None:
        print(f'Filter based on path id: {tracks_id}')
        for track in tracks_id:
            ix = scenes_df[scenes_df.track_number != int(track)].index
            scenes_df = scenes_df.drop(index = ix)

    # filter based on min temporal baseline
    if min_btemp is not None:
        print(f'Filter based on minimum temporal baseline {min_btemp}')
        scenes_df = scenes_df[scenes_df.btemp > min_btemp]
    # filter based on max temporal baseline
    if max_btemp is not None:
        print(f'Filter based on maximum temporal baseline {max_btemp}')
        scenes_df = scenes_df[scenes_df.btemp < max_btemp]
        
    # NOTE; GUNW version should be converted to number to allow 
    # more simple filtering, e.g all v2 would be < 3

    if gunw_version is not None:
        print(f'Filter based on GUNW version: {gunw_version}')
        scenes_df = scenes_df[scenes_df.version == gunw_version]    

    # Get final results
    results = [results[i] for i in scenes_df.index]
    print(f'Final granules number: {len(results)}')
    
    return results, scenes_df

def download_earthdata(earthdata_results:list, output_dir: str = '.', threads: int = 1):
    # Login 
    auth = earthaccess.login()

    if auth.authenticated is False:
        auth = earthaccess.login(strategy="interactive", persist=True)

    # Streaming data
    # https://earthaccess.readthedocs.io/en/latest/tutorials/emit-earthaccess/
    try:
        downloaded_files = earthaccess.download(earthdata_results[0],
                                                    local_path=output_dir,
                                                    threads=1)
    except:
        raise ValueError('Downloading failed, check earthdata credentials, e.g. .netrc')

    # if passed the check, start downloading the data
    downloaded_files = earthaccess.download(earthdata_results,
                                            local_path=output_dir,
                                            threads=threads) 

    return downloaded_files


class Downloader(object):
    def __init__(self, inps):
        self.inps            = inps
        self.inps.output     = self.inps.output.title()
        self.inps.wd         = op.abspath(self.inps.wd)
        os.makedirs(self.inps.wd, exist_ok=True)
        if self.inps.verbose:
            log.setLevel('DEBUG')
        else:
            log.setLevel('INFO')


    def __call__(self):
        if self.inps.api == 'asf_search':
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
                log.info('Found -- %d -- products', len(scenes))

            # elif self.inps.output == 'Kml':
            #     dst    = fmt_dst(inps)
            #     log.error('Kml option is not yet supported. '\
            #                    'Revert to an older version of ARIAtools')

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
                log.info (f'Downloading {len(scenes)} products...')
                if self.inps.user is not None:
                    session = asf.ASFSession()
                    session.auth_with_creds(self.inps.user, self.inps.passw)
                    scenes.download(self.inps.wd, processes=nt, session=session)

                else:
                    scenes.download(self.inps.wd, processes=nt)
                log.info(f'Download complete. Wrote -- {len(scenes)} -- products to: {self.inps.wd}')

            if inps.verbose:
                for scene in scenes:
                    print(scene.geojson()['properties']['sceneName'])

        elif self.inps.api == 'earthaccess':
            # Get bbox
            bbox = make_bbox(self.inps.bbox)
            bbox = bbox.bounds if bbox is not None else None

            # Start and end data
            start_date = self.inps.start.strftime('%Y-%m-%d')
            end_date = self.inps.end.strftime('%Y-%m-%d')

            # Tracks, make sure it is a list of int
            if len(self.inps.track.split(',')) > 0 :
                tracks = [int(track) for track in self.inps.track.split(',')] 
            else:
                tracks = [self.inps.track]
   
            #  fix min and max temporal baseline 
            min_temp = self.inps.daysgt if self.inps.daysgt != 0 else None
            max_btemp = self.inps.dayslt if self.inps.dayslt != math.inf else None

            scenes, query_df = query_earthaccess(bbox,
                                                 flight_direction=self.inps.flightdir,
                                                 tracks_id=tracks,
                                                 start_date=start_date,
                                                 end_date=end_date,
                                                 min_btemp=min_temp,
                                                 max_btemp=max_btemp,
                                                 gunw_version=self.inps.version)
 
            if self.inps.output == 'Count':
                log.info('Found -- %d -- products', len(scenes))
            elif self.inps.output == 'Url':
                dst  = fmt_dst(inps)
                with open(dst, 'w') as fh: [print(url, sep='\n', file=fh) \
                                            for url in query_df.URL]
                log.info(f'Wrote -- {len(urls)} -- product urls to: {dst}')
            elif self.inps.output == 'Download':
                log.info (f'Downloading {len(scenes)} products...')
                downloaded = download_earthdata(scenes, output_dir=self.inps.wd, threads=int(self.inps.num_threads))
                log.info(f'Download complete. Wrote -- {len(scenes)} -- products to: {self.inps.wd}')


        return


    def query_asf(self):
        """Get the scenes from ASF"""
        bbox = make_bbox(self.inps.bbox)
        bbox = bbox.wkt if bbox is not None else None
        if self.inps.flightdir is not None:
            if self.inps.flightdir.lower()[0] == 'a':
                flight_direction = 'ascending'
            elif self.inps.flightdir == 'd':
                self.inps.flightdir.lower()[0]  = 'descending'
        else:
            flight_direction = None

        if self.inps.track is not None:
            tracks = self.inps.track.split(',')
            tracks = [int(track) for track in tracks]
        else:
            tracks = self.inps.track

        dct_kw = dict(platform=asf.constants.SENTINEL1,
                      processingLevel=asf.constants.GUNW_STD,
                      relativeOrbit=tracks,
                      flightDirection=flight_direction,
                      intersectsWith=bbox)
        scenes = asf.geo_search(**dct_kw)


        return scenes


if __name__ == '__main__':
    try:
        print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except:
        pass

    inps = cmdLineParse()
    Downloader(inps)()
