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
import math
import re
import json
import requests
import argparse
import importlib
import signal
from datetime import datetime, timedelta

def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Program to download all GUNW products for a single track.\nDetails on parameters at: https://www.asf.alaska.edu/get-data/learn-by-doing/',
                                     epilog='Examples of use:\n\t productAPI.py --track 004 --output count\n\t productAPI.py --bbox "36.75 37.225 -76.655 -75.928"\n\t productAPI.py -t 004 --start "1 year ago" --end now',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', type=str, help='Output type, default is "Download". "Download", "Count", and "Kml" are currently supported.')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str, help='track to download')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str, help='Lat/Lon Bounding SNWE')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', type=str, help='Specify directory to deposit all outputs. Default is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=None, type=str, help='Start date; If none provided, starts at beginning of Sentinel record (2015). See https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    parser.add_argument('-e', '--end', dest='end', default=None, type=str, help='End date; If none provided, ends today. See https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    parser.add_argument('-g', '--daysgt', dest='daysgt', default=None, type=int, help='Take pairs with a temporal baseline -- days greater than this value. Example, annual pairs: productAPI.py -t 004 --daysgt 364.')
    parser.add_argument('-l', '--dayslt', dest='dayslt', default=None, type=int, help='Take pairs with a temporal baseline -- days less than this value.')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, type=str, help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('-v', '--verbose', dest='v', default=True, type=bool, help='Print products to be downloaded to stdout')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(1)

    inps = parser.parse_args(args=iargs)

    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    i1 = []
    for i in [inps.start, inps.end]:
        if not i or i == 'now' or i == 'today' or i.endswith('ago'):
            i1.append(i)
        elif len(i) == 8:
            i1.append(datetime.strptime(i, '%Y%m%d').strftime('%Y-%m-%d'))
        elif len(i) == 10:
            i1.append(i)
        else:
            raise Exception('Date format not understood; see https://www.asf.alaska.edu/get-data/learn-by-doing/ for options')

    inps.start = i1[0]; inps.end = i1[1]
    return inps

class Downloader(object):
    def __init__(self, inps):
        self.inps        = inps
        self.inps.output = self.inps.output.title()
        self.inps.wd     = op.abspath(self.inps.wd)
        self.url_base    = 'https://api.daac.asf.alaska.edu/services/search/param?'
        if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)

    def __call__(self):
        url       = self.form_url()
        dict_prod = self.parse_json(url)

        script = requests.post('{}&output={}'.format(self.url_base,
                                        self.inps.output), data=dict_prod).text
        if self.inps.output == 'Count':
            print ('\nFound -- {} -- products'.format(script.strip()))

        elif self.inps.output == 'Kml':
            dst = self._fmt_dst()
            print (script, file=open(dst, 'w'))
            print ('Wrote .KMZ to:\n\t {}'.format(dst))

        elif self.inps.output == 'Download':
            os.chdir(self.inps.wd)
            os.sys.argv = []
            exec(script, globals())
        return


    def form_url(self):
        url = '{}asfplatform=Sentinel-1%20Interferogram%20(BETA)&output=JSON'.format(self.url_base)
        if self.inps.track:
            url += '&relativeOrbit={}'.format(self.inps.track)
        if self.inps.bbox:
            url += '&bbox=' + self._get_bbox()
        if self.inps.start:
            url += '&start={}'.format(self.inps.start)
        if self.inps.end:
            url += '&end={}'.format(self.inps.end)
        if self.inps.flightdir:
            url += '&flightDirection={}'.format(self.inps.flightdir.upper())

        url = url.replace(' ', '+')
        print (url)
        return url

    def parse_json(self, url):
        j        = json.loads(requests.get(url).text)[0]
        if len(j) == 0:
            raise Exception('No products found with given url')

        prod_ids = []
        for prod in j:
            id = prod['product_file_id']
            if 'layer' in prod['downloadUrl']: continue
            elif (self.inps.daysgt or self.inps.dayslt):
                match = re.search(r'\d{8}_\d{8}', id).group()
                end, st = [datetime.strptime(i, '%Y%m%d').date() for i in match.split('_')]
                elap = (end - st).days
                if self.inps.daysgt and not (elap > self.inps.daysgt): continue
                if self.inps.dayslt and not (elap < self.inps.dayslt): continue

            prod_ids.append(id)
            if self.inps.v: print ('Will process: {}'.format(id))

        data = {'product_list': ','.join(prod_ids)}
        return data

    ## utility functions
    def _get_bbox(self):
        if op.exists(op.abspath(self.inps.bbox)):
            from shapefile_util import open_shapefile
            bounds = open_shapefile(self.inps.bbox, 0, 0).bounds
            W, S, E, N = [str(i) for i in bounds]
        else:
            try:
                S, N, W, E = self.inps.bbox.split()
            except:
                raise Exception('Cannot understand the --bbox argument. Input string was entered incorrectly or path does not exist.')
        return ','.join([W,S,E,N])

    def _fmt_dst(self):
        dst_base = op.join(self.inps.wd, 'download_products')
        if self.inps.track:
            dst = '{}_{}'.format(dst_base, self.inps.track)
        elif self.inps.bbox:
            WSEN     = self._get_bbox().split(',')
            WSEN_fmt = []
            for i, coord in enumerate(WSEN):
                if i < 2:
                    WSEN_fmt.append(math.floor(float(coord)))
                else:
                    WSEN_fmt.append(math.ceil(float(coord)))
            dst = '{}_{}W{}S{}E{}N.kmz'.format(dst_base, *WSEN_fmt)

        dst += '.kmz'
        return dst

if __name__ == '__main__':
    inps = cmdLineParse()

    Downloader(inps)()
