#! /usr/bin/env python3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import os.path as op
import shutil
import math
import re
import json
import requests
import argparse
from datetime import datetime


def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Command line interface to download GUNW products from the ASF DAAC. GUNW products are hosted at the NASA ASF DAAC.\nDownloading them requires a NASA Earthdata URS user login and requires users to add "ARIA Product Search" to their URS approved applications.',
                                     epilog='Examples of use:\n\t ariaDownload.py --track 004 --output count\n\t ariaDownload.py --bbox "36.75 37.225 -76.655 -75.928"\n\t ariaDownload.py -t 004,077 --start 20190101 -o count',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', type=str, help='Output type, default is "Download". "Download", "Count", "Url" and "Kmz" are currently supported. Use "Url" for ingestion to aria*.py')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str, help='track to download; single number (including leading zeros) or comma separated')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str, help='Lat/Lon Bounding SNWE, or GDAL-readable file containing POLYGON geometry.')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', type=str, help='Specify directory to deposit all outputs. Default is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=None, type=str, help='Start date as YYYYMMDD; If none provided, starts at beginning of Sentinel record (2014).')
    parser.add_argument('-e', '--end', dest='end', default=None, type=str, help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument('-l', '--daysless', dest='dayslt', default=None, type=int, help='Take pairs with a temporal baseline -- days less than this value.')
    parser.add_argument('-m', '--daysmore', dest='daysgt', default=None, type=int, help='Take pairs with a temporal baseline -- days greater than this value. Example, annual pairs: ariaDownload.py -t 004 --daysmore 364.')
    parser.add_argument('-i', '--ifg', dest='ifg', default=None, type=str, help='Retrieve one interferogram by its start/end date, specified as YYYYMMDD_YYYYMMDD (order independent)')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, type=str, help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('-v', '--verbose', dest='v', action='store_true', help='Print products to be downloaded to stdout')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)

    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    if not inps.output.lower() in ['count', 'kml', 'kmz', 'url', 'download']:
        raise Exception ('Incorrect output keyword. Choose "count", "kmz", "url", or "download"')

    inps.output = 'Kml' if inps.output.lower() == 'kmz' else inps.output.title()
    return inps

class Downloader(object):
    def __init__(self, inps):
        self.inps        = inps
        self.inps.output = self.inps.output.title()
        self.inps.wd     = op.abspath(self.inps.wd)
        self.url_base    = 'https://api.daac.asf.alaska.edu/services/search/param?'

    def __call__(self):
        url              = self.form_url()
        dict_prod, urls  = self.parse_json(url)
        script           = requests.post('{}&output={}'.format(self.url_base,
                                        self.inps.output), data=dict_prod).text
        if self.inps.output == 'Count':
            print ('\nFound -- {} -- products'.format(len(urls)))

        elif self.inps.output == 'Kml':
            if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)
            dst = self._fmt_dst()
            print (script, file=open(dst, 'w'))
            print ('Wrote .KMZ to:\n\t {}'.format(dst))

        elif self.inps.output == 'Url':
            if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)
            dst = self._fmt_dst()
            with open(dst, 'w') as fh: [print(url, sep='\n', file=fh) for url in urls]
            print ('Wrote -- {} -- product urls to: {}'.format(len(urls), dst))

        elif self.inps.output == 'Download':
            if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)
            os.chdir(self.inps.wd)
            os.sys.argv = []
            fileName = os.path.abspath(os.path.join(self.inps.wd,'ASFDataDload.py'))
            with open(fileName, 'w') as f:
                f.write(script)

            os.sys.path.append(os.path.abspath(self.inps.wd))
            import ASFDataDload as AD
            downloader = AD.bulk_downloader()
            downloader.download_files()
            downloader.print_summary()
            # Delete temporary files
            shutil.rmtree(os.path.abspath(os.path.join(self.inps.wd,'__pycache__')))
            os.remove(fileName)

        return urls

    def form_url(self):
        url = '{}asfplatform=Sentinel-1%20Interferogram%20(BETA)&processingLevel=GUNW_STD&output=JSON'.format(self.url_base)
        if self.inps.track:
            url += '&relativeOrbit={}'.format(self.inps.track)
        if self.inps.bbox:
            url += '&bbox=' + self._get_bbox()
        if self.inps.flightdir:
            url += '&flightDirection={}'.format(self.inps.flightdir.upper())

        url = url.replace(' ', '+')
        print (url)
        return url

    def parse_json(self, url):
        j        = json.loads(requests.get(url).text)[0]
        if len(j) == 0:
            raise Exception('No products found with given url; check inputs for errors.')

        prod_ids = []; dl_urls = []
        i=0
        for prod in j:
            if 'layer' in prod['downloadUrl']: continue
            i+=1
            FileId = prod['product_file_id']
            match = re.search(r'\d{8}_\d{8}', FileId).group()
            dates = [datetime.strptime(i, '%Y%m%d').date() for i in match.split('_')]
            st, end = sorted(dates)

            if self.inps.ifg:
                dates1 = [datetime.strptime(i, '%Y%m%d').date() for i in self.inps.ifg.split('_')]
                st1, end1 = sorted(dates1)
                if st1 != st or end1 != end: continue

            if (self.inps.start or self.inps.end):
                if self.inps.start and not (st >= datetime.strptime(self.inps.start, '%Y%m%d').date()): continue
                if self.inps.end and not (end <= datetime.strptime(self.inps.end, '%Y%m%d').date()): continue

            if (self.inps.daysgt or self.inps.dayslt):
                elap = (end - st).days
                if self.inps.daysgt and not (elap > self.inps.daysgt): continue
                if self.inps.dayslt and not (elap < self.inps.dayslt): continue

            prod_ids.append(FileId); dl_urls.append(prod['downloadUrl'])

            if self.inps.v: print ('Found: {}'.format(FileId))

        if len(prod_ids) == 0:
            raise Exception('No products found that satisfy requested conditions.')

        data = {'product_list': ','.join(prod_ids)}
        return data, dl_urls

    ## utility functions
    def _get_bbox(self):
        if op.exists(op.abspath(self.inps.bbox)):
            from ARIAtools.shapefile_util import open_shapefile
            bounds = open_shapefile(self.inps.bbox, 0, 0).bounds
            W, S, E, N = [str(i) for i in bounds]
        else:
            try:
                S, N, W, E = self.inps.bbox.split()
            except:
                raise Exception('Cannot understand the --bbox argument. Input string was entered incorrectly or path does not exist.')
        return ','.join([W,S,E,N])

    def _fmt_dst(self):
        dst = op.join(self.inps.wd, 'download_products')
        ext      = '.kmz' if self.inps.output == 'Kmz' else '.txt'
        if self.inps.track:
            dst = '{}_{}track'.format(dst, self.inps.track).replace(',', '-')

        if self.inps.bbox:
            WSEN     = self._get_bbox().split(',')
            WSEN_fmt = []
            for i, coord in enumerate(WSEN):
                if i < 2:
                    WSEN_fmt.append(math.floor(float(coord)))
                else:
                    WSEN_fmt.append(math.ceil(float(coord)))
            dst = '{}_{}W{}S{}E{}Nbbox'.format(dst, str(WSEN_fmt[0]), str(WSEN_fmt[1]), str(WSEN_fmt[2]), str(WSEN_fmt[3]))
        dst  += '_0{}'.format(ext)
        count = 1 # don't overwrite if already exists
        while op.exists(dst):
            basen  = '{}{}{}'.format(re.split(str(count-1)+ext, op.basename(dst))[0], count, ext)
            dst    = op.join(op.dirname(dst), basen)
            count += 1
        return dst

if __name__ == '__main__':
    inps = cmdLineParse()
    Downloader(inps)()
