#! /usr/bin/env python3

import os
import os.path as op
import requests
import argparse
import importlib
import signal
from datetime import datetime


def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Program to download all GUNW products for a single track. Details on parameter at: https://www.asf.alaska.edu/get-data/learn-by-doing/')
    parser.add_argument('-t', '--track', dest='track', type=str, required=True, help='track to download')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', help='Specify directory to deposit all outputs. Default is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=datetime.now().strftime('%Y%m%d'), help='Start date; see https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    parser.add_argument('-e', '--end', dest='end', default='20150101', help='End date; see https://www.asf.alaska.edu/get-data/learn-by-doing/ for input options')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    return parser.parse_args(args=iargs)

class Downloader(object):
    def __init__(self, inps):
        self.inps  = inps
        self.inps.wd = op.abspath(self.inps.wd)
        if not op.exists(self.inps.wd): os.mkdir(self.inps.wd)
        os.sys.path.append(self.inps.wd); os.chdir(self.inps.wd)


    def __call__(self):
        url         = self.form_url()
        path_script = self.write_dl_script(url)
        module = importlib.import_module(op.splitext(op.basename(path_script))[0])

        # setup and run the newly created download script
        signal.signal(signal.SIGINT, self.signal_handler)

        downloader = module.bulk_downloader()
        downloader.download_files()
        downloader.print_summary()

    def form_url(self):
        url_base  = 'https://api.daac.asf.alaska.edu/services/search/param?asfplatform=Sentinel-1%20Interferogram%20(BETA)&output=Download'
        url = '{}&relativeOrbit={}'.format(url_base, self.inps.track)
        # start and end not currently supported
        # url += '&start={}'.format(self.inps.start) + '&end={}'.format(self.inps.end)
        print (url)
        return url

    def write_dl_script(self, url):
        """Remove non gunw products from a download script """
        script = requests.get(url).text.split('\n')
        script_gunw = []
        for i, line in enumerate(script):
            if 'services' in line and 'product' in line:
                if not script[i+1]:
                    script_gunw.append(']')
            else:
                script_gunw.append(line)
        script_exec = '\n'.join(script_gunw)

        path_script = op.join(self.inps.wd, 'download_products_{}.py'.format(self.inps.track))
        print (script_exec, file=open(path_script, 'w'))

        return path_script

    # from bulk download script
    def signal_handler(sig, frame):
        """ A routine that handles trapped signals """
        global abort
        sys.stderr.output("\n > Caught Signal. Exiting!\n")
        abort = True # necessary to cause the program to stop
        raise SystemExit  # this will only abort the thread that the ctrl+c was caught in

if __name__ == '__main__':
    ## example call to get all track 4 GUNW products and store in ./track_004/products
    # dl_wrapper.py -t 004 -w ./track_004
    track = '004'
    inps  = cmdLineParse()
    Downloader(inps)()
