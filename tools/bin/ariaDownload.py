#! /usr/bin/env python3

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os, os.path as op
import argparse
import shutil, math
import re, json, requests
import numpy as np
from datetime import datetime
import logging
from ARIAtools.logger import logger

log = logging.getLogger('ARIAtools')

def createParser():
    """ Download a bulk download script and execute it """
    parser = argparse.ArgumentParser(description='Command line interface to download GUNW products from the ASF DAAC. GUNW products are hosted at the NASA ASF DAAC.\nDownloading them requires a NASA Earthdata URS user login and requires users to add "GRFN Door (PROD)" and "ASF Datapool Products" to their URS approved applications.',
                                     epilog='Examples of use:\n\t ariaDownload.py --track 004 --output count\n\t ariaDownload.py --bbox "36.75 37.225 -76.655 -75.928"\n\t ariaDownload.py -t 004,077 --start 20190101 -o count',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', dest='output', default='Download', type=str, help='Output type, default is "Download". "Download", "Count", "Url" and "Kmz" are currently supported. Use "Url" for ingestion to aria*.py')
    parser.add_argument('-t', '--track', dest='track', default=None, type=str, help='track to download; single number (including leading zeros) or comma separated')
    parser.add_argument('-b', '--bbox', dest='bbox',  default=None, type=str, help='Lat/Lon Bounding SNWE, or GDAL-readable file containing POLYGON geometry.')
    parser.add_argument('-w', '--workdir', dest='wd', default='./products', type=str, help='Specify directory to deposit all outputs. Default is "products" in local directory where script is launched.')
    parser.add_argument('-s', '--start', dest='start', default=None, type=str, help='Start date as YYYYMMDD; If none provided, starts at beginning of Sentinel record (2014).')
    parser.add_argument('-e', '--end', dest='end', default=None, type=str, help='End date as YYYYMMDD. If none provided, ends today.')
    parser.add_argument('-u', '--user', dest='user', default=None, type=str, help='NASA Earthdata URS user login. Users must add "GRFN Door (PROD)" and "ASF Datapool Products" to their URS approved applications.')
    parser.add_argument('-p', '--pass', dest='passw', default=None, type=str, help='NASA Earthdata URS user password. Users must add "GRFN Door (PROD)" and "ASF Datapool Products" to their URS approved applications.')
    parser.add_argument('-l', '--daysless', dest='dayslt', default=None, type=int, help='Take pairs with a temporal baseline -- days less than this value.')
    parser.add_argument('-m', '--daysmore', dest='daysgt', default=None, type=int, help='Take pairs with a temporal baseline -- days greater than this value. Example, annual pairs: ariaDownload.py -t 004 --daysmore 364.')
    parser.add_argument('-nt', '--num_threads', dest='num_threads', default='1', type=str, help='Specify number of threads for multiprocessing download. By default "1". Can also specify "All" to use all available threads.')
    parser.add_argument('-i', '--ifg', dest='ifg', default=None, type=str, help='Retrieve one interferogram by its start/end date, specified as YYYYMMDD_YYYYMMDD (order independent)')
    parser.add_argument('-d', '--direction', dest='flightdir', default=None, type=str, help='Flight direction, options: ascending, a, descending, d')
    parser.add_argument('-v', '--verbose', dest='v', action='store_true', help='Print products to be downloaded to stdout')
    return parser

def cmdLineParse(iargs=None):
    parser = createParser()
    if len(os.sys.argv) < 2:
        parser.print_help()
        os.sys.exit(2)

    inps = parser.parse_args(args=iargs)

    if not inps.track and not inps.bbox:
        raise Exception('Must specify either a bbox or track')

    if not inps.output.lower() in ['count', 'kmz', 'kml', 'url', 'download']:
        raise Exception ('Incorrect output keyword. Choose "count", "kmz", "url", or "download"')

    inps.output = 'Kml' if inps.output.lower() == 'kmz' or inps.output.lower() == 'kml' else inps.output.title()
    return inps

class Downloader(object):
    def __init__(self, inps):
        self.inps        = inps
        self.inps.output = self.inps.output.title()
        self.inps.wd     = op.abspath(self.inps.wd)
        self.url_base    = 'https://api.daac.asf.alaska.edu/services/search/param?'
        if self.inps.v: log.setLevel('DEBUG')
        return

    def __call__(self):
        url              = self.form_url()
        dct_prod, urls  = self.parse_json(url)
        script = requests.post(f'{self.url_base}&output={self.inps.output}', data=dct_prod).text


        if self.inps.output == 'Count':
            log.info('\nFound -- %d -- products', len(urls))

        elif self.inps.output == 'Kml':
            os.makedirs(self.inps.wd, exist_ok=True)
            dst    = self._fmt_dst()
            script = requests.post(f'{self.url_base}&output={self.inps.output}', data=dct_prod).text
            print(script, file=open(dst, 'w'))
            log.info(f'Wrote .KMZ to:\n\t %s', dst)

        elif self.inps.output == 'Url':
            os.makedirs(self.inps.wd, exist_ok=True)
            dst = self._fmt_dst()
            with open(dst, 'w') as fh: [print(url, sep='\n', file=fh) for url in urls]
            log.info(f'Wrote -- {len(urls)} -- product urls to: {dst}')

        elif self.inps.output == 'Download':
            os.makedirs(self.inps.wd, exist_ok=True)
            os.chdir(self.inps.wd)
            os.sys.path.append(os.path.abspath(self.inps.wd))
            ## make a cookie from a .netrc that ASFDataDload will pick up
            self.make_nc_cookie()
            inps.url_base = self.url_base
            prod_dl(self.inps, dct_prod)

            # Delete temporary files
            shutil.rmtree(op.abspath(op.join(self.inps.wd,'__pycache__')))

        return urls


    def form_url(self):
        url = f'{self.url_base}asfplatform=Sentinel-1%20Interferogram%20(BETA)&processingLevel=GUNW_STD&output=JSON'
        if self.inps.track:
            url += f'&relativeOrbit={self.inps.track}'
        if self.inps.bbox:
            url += '&bbox=' + self._get_bbox()
        if self.inps.flightdir:
            url += f'&flightDirection={self.inps.flightdir.upper()}'

        url = url.replace(' ', '+')
        log.info(url)
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

            log.debug('Found: %s', FileId)

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
        ext      = '.kmz' if self.inps.output == 'Kml' else '.txt'
        if self.inps.track:
            dst = f'{dst}_{self.inps.track}track'.replace(',', '-')

        if self.inps.bbox:
            WSEN     = self._get_bbox().split(',')
            WSEN_fmt = []
            for i, coord in enumerate(WSEN):
                if i < 2:
                    WSEN_fmt.append(math.floor(float(coord)))
                else:
                    WSEN_fmt.append(math.ceil(float(coord)))
            dst = f'{dst}_{WSEN_fmt[0]}W{WSEN_fmt[1]}S{WSEN_fmt[2]}E{WSEN_fmt[3]}Nbbox'
        dst  += f'_0{ext}'
        count = 1 # don't overwrite if already exists
        while op.exists(dst):
            basen  = f'{re.split(str(count-1)+ext, op.basename(dst))[0]}{count}{ext}'
            dst    = op.join(op.dirname(dst), basen)
            count += 1
        return dst

    def make_nc_cookie(self):
       import netrc, base64
       from urllib.request import build_opener, Request
       from urllib.request import HTTPHandler, HTTPSHandler, HTTPCookieProcessor
       from urllib.error import HTTPError, URLError
       from http.cookiejar import MozillaCookieJar

       cookie_jar_path = op.join(op.expanduser('~'), '.bulk_download_cookiejar.txt')
       asf_urs4        = {'url': 'https://urs.earthdata.nasa.gov/oauth/authorize',
                         'client': 'BO_n7nTIlMljdvU6kRRB3g',
                         'redir': 'https://auth.asf.alaska.edu/login'}
       if self.inps.user is not None and self.inps.passw is not None:
           user  = self.inps.user
           passw = self.inps.passw
       elif op.exists(op.join(op.expanduser('~'), '.netrc')):
           log.info('Attempting to obtaining user/pass from .netrc')
           try:
               user, _,  passw = netrc.netrc().authenticators('urs.earthdata.nasa.gov')
           except:
               log.warning('Could not obtain credentials from existing .netrc')
               return None
       else:
           # resort to ASF credential checks (will prompt for input if can't find the cookiejar)
           return None

       # Build URS4 Cookie request
       auth_cookie_url = f"{asf_urs4['url']}?client_id={asf_urs4['client']}&redirect_uri={asf_urs4['redir']}&response_type=code&state="
       user_pass = base64.b64encode(bytes(f'{user}:{passw}', 'utf-8')).decode('utf-8')

       # Authenticate against URS, grab all the cookies
       cookie_jar = MozillaCookieJar()
       opener     = build_opener(HTTPCookieProcessor(cookie_jar), HTTPHandler(), HTTPSHandler(**{}))
       request    = Request(auth_cookie_url, headers={'Authorization': f'Basic {user_pass}'})

       # Watch out cookie rejection!
       try:
          response = opener.open(request)
       except HTTPError as e:
          if "WWW-Authenticate" in e.headers and "Please enter your Earthdata Login credentials" in e.headers["WWW-Authenticate"]:
             log.info(' > Username and Password combo was not successful. Please try again.')
             return False
          else:
             # If an error happens here, the user most likely has not confirmed EULA.
             log.info('\nThere was an error obtaining a download cookie')
             log.info('Most likely you lack permission to download data from the ASF Datapool.')
             log.info('\n\nNew users: you must first log into Vertex and accept the EULA. In addition, your Study Area must be set at Earthdata https://urs.earthdata.nasa.gov')
             os.sys.exit(1)

       except URLError:
          log.info('\nThere was a problem communicating with URS, unable to obtain cookie')
          log.info('Try cookie generation later.')
          os.sys.exit(1)

       # Did we get a cookie?
       if check_cookie_is_logged_in(cookie_jar):
          #COOKIE SUCCESS!
          cookie_jar.save(cookie_jar_path)
          return True

       # if we aren't successful generating the cookie, nothing will work. Stop here!
       log.warning('Could not generate new cookie! Cannot proceed. Please check credentials and/or .netrc.')
       log.info(f'Response was {response.getcode()}')
       log.info('\n\nNew users: you must first log into Vertex and accept the EULA. In addition, your Study Area must be set at Earthdata https://urs.earthdata.nasa.gov')
       os.sys.exit(1)

def prod_dl(inps, dct_prod):
    """ Perform downloading using ASF bulk dl; parallel processing supported """
    import multiprocessing
    max_threads = multiprocessing.cpu_count()

    if inps.num_threads == 'all':
        nt = max_threads
    else:
        nt = int(inps.num_threads)
    nt     = nt if nt < max_threads else max_threads
    log.info('Using  %s threads for parallel downloads', nt)

    files          = dct_prod['product_list'].split(',')
    chunks1        = np.array_split(files, np.ceil(len(files)/200)) # split by 200s
    check          = []
    for chunk in chunks1:
        chunks     = np.array_split(chunk, nt)
        lst_dcts   = [vars(inps).copy() for i in range(nt)]  # Namespace->dictionary, repeat it
        for i, chunk1 in enumerate(chunks):
            lst_dcts[i]['files'] = chunk1 # put split up files to the objects for threads
            lst_dcts[i]['ext']   = i      # for unique bulk downloader file name

        with multiprocessing.Pool(nt) as pool:
            try:
                info = pool.map(_dl_helper, lst_dcts)
            except Exception as E:
                print ('ASF bulk downloader error:', E)
                print ('Likely a bad handshake with the ASF DAAC. Try rerunning')
                os.sys.exit(1)

        check.extend(chunkc for chunkc in chunks) # in case products missed in split

    rewrite_summary(info)

    check = [item for sublist in check for item in sublist]
    for f in files:
        if not f in check:
            log.critical('File splitting missed: %s; subset your ARIA in command in time to try again', f)

def _dl_helper(inp_dct):
    """ Helper function for parallel processing """
    prod_dct = {'product_list': ','.join(inp_dct['files'])}
    log.debug ('# of files: %d', len(inp_dct['files']))
    script   = requests.post(f'{inp_dct["url_base"]}&output=Download', data=prod_dct).text

    fname    = f'ASFDataDload{inp_dct["ext"]}.py'
    fpath    = os.path.abspath(os.path.join(inp_dct['wd'], fname))
    with open(fpath, 'w') as f: f.write(script)
    AD = __import__(os.path.splitext(fname)[0])
    os.sys.argv = [] # gets around spurious messages
    dler  = AD.bulk_downloader()
    dler.download_files()

    return dler.success, dler.failed, dler.skipped, dler.total_time, dler.total_bytes

def check_cookie_is_logged_in(cj):
    """Make sure successfully logged into URS; try to get cookie"""
    for cookie in cj:
        if cookie.name == 'urs_user_already_logged':
            return True
    return False

def rewrite_summary(infos):
    succeed, failed, skipped, tot_time, tot_bytes = [], [], [], 0, 0
    for info in infos:
        if info[0]: succeed.extend(info[0])
        if info[1]: failed.extend(info[1])
        if info[2]: skipped.extend(info[2])
        tot_time += float(info[3])
        tot_bytes += float(info[4])


    log.info ('Successes: %d files, %s bytes', len(succeed), tot_bytes)
    [print ('\n\t-', sf['file'], f'{sf["size"]/1024**2:.2f}') for sf in succeed]

    if skipped: log.info('Skipped: %d files', len(skipped))
    [print ('\n\t-', sf) for sf in skipped]

    if failed: log.info('Failures: %d files', len(failed))
    [print ('\n\t-', sf) for sf in failed]

    if succeed: log.info('Average Rate: %.2f MB/sec', (tot_bytes/1024.0**2)/tot_time)

    if failed:
        log.critical('We recommend rerunning the same ariaDownload command to address the %d failures', len(failed))
    else:
        log.info ('All files have been downloaded successfully')

if __name__ == '__main__':
    inps = cmdLineParse()
    Downloader(inps)()
