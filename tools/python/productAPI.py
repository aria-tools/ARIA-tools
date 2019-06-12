#! /usr/bin/env python3

import os
import os.path as op
import math
import requests
import argparse
import importlib
import signal
from datetime import datetime

def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Program to download all GUNW products for a single track.\nDetails on parameters at: https://www.asf.alaska.edu/get-data/learn-by-doing/',
                                     epilog='Examples of use:\n\t productAPI.py --track 004 --output count\n\t productAPI.py --bbox "36.75 37.225 -76.655 -75.928"\n\t productAPI.py -t 004 --start "1 year ago" --now now',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', type=str, help='Output type, default is "Download". Only "Download" and "Count" are currently supported.')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str, help='track to download')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str, help='Lat/Lon Bounding SNWE')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', type=str, help='Specify directory to deposit all outputs. Default is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=None, type=str, help='Start date; If none provided, starts at beginning of Sentinel record (2015). See https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    parser.add_argument('-e', '--end', dest='end', default=None, type=str, help='End date; If none provided, ends today. See https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, type=str, help='Flight direction, options: ascending, a, descending, d')
    return parser

class Downloader(object):
    def __init__(self, inps):
        self.inps        = inps
        self.inps.output = self.inps.output.title()
        self.inps.wd     = op.abspath(self.inps.wd)

        if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)
        os.sys.path.append(self.inps.wd);

    def __call__(self):
        url         = self.form_url()
        path_script = self.write_dl_script(url)
        module = importlib.import_module(op.splitext(op.basename(path_script))[0])

        # setup and run the newly created download script
        os.chdir(self.inps.wd)
        signal.signal(signal.SIGINT, self.signal_handler)
        downloader = module.bulk_downloader()
        downloader.download_files()
        downloader.print_summary()

    def form_url(self):
        url_base  = 'https://api.daac.asf.alaska.edu/services/search/param?asfplatform=Sentinel-1%20Interferogram%20(BETA)'
        url = '{}&output={}'.format(url_base, self.inps.output)
        if self.inps.track:
            url += '&relativeOrbit={}'.format(self.inps.track)
        if self.inps.bbox:
            url += '&bbox=' + self.get_bbox()
        if self.inps.start:
            url += '&start={}'.format(self.inps.start)
        if self.inps.end:
            url += '&end={}'.format(self.inps.end)
        if self.inps.flightdir:
            url += '&flightDirection={}'.format(self.inps.flightdir.upper())

        url = url.replace(' ', '+')
        print (url)
        return url

    def write_dl_script(self, url):
        """Remove non gunw products from a download script """
        script = requests.get(url).text.split('\n')
        dst    = self._fmt_dst()
        if self.inps.output == 'Count':
            print ('\nFound -- {} -- products, INCLUDING the on-demand layers'.format(script[0]))
            quit()

        elif self.inps.output == 'Kml':
            print ('\n'.join(script), file=open(dst, 'w'))
            print ('Wrote .KMZ to:\n\t {}'.format(dst))
            quit()

        elif self.inps.output == 'Download':
            script_gunw = []
            for i, line in enumerate(script):
                if 'services' in line and 'product' in line:
                    if not script[i+1]:
                        script_gunw.append(']')
                        continue
                else:
                    script_gunw.append(line)
            script_exec = '\n'.join(script_gunw)
            print (script_exec, file=open(dst, 'w'))


        else:
            raise Exception('{} output type not currently supported'.format(self.inps.output))
        return dst

    # from bulk download script
    def signal_handler(sig, frame):
        """ A routine that handles trapped signals """
        global abort
        sys.stderr.output("\n > Caught Signal. Exiting!\n")
        abort = True # necessary to cause the program to stop
        raise SystemExit  # this will only abort the thread that the ctrl+c was caught in
        return

    def get_bbox(self):
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
            WSEN     = self.get_bbox().split(',')
            WSEN_fmt = []
            for i, coord in enumerate(WSEN):
                if i < 2:
                    WSEN_fmt.append(math.floor(float(coord)))
                else:
                    WSEN_fmt.append(math.ceil(float(coord)))
            dst = '{}_{}W{}S{}E{}N'.format(dst_base, *WSEN_fmt)

        if self.inps.output == 'Kml':
            dst += '.kmz'
        elif self.inps.output == 'Download':
            dst += '.py'
        return dst


if __name__ == '__main__':
    inps  = createParser().parse_args()
    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    Downloader(inps)()
