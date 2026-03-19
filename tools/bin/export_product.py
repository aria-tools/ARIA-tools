#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import sys
import argparse
import logging
import json
import shutil

# --- INJECT THIS BLOCK BEFORE ANY ARIA OR GDAL IMPORTS ---
# 1. Inherit Earthdata cookies from parent to avoid unauthenticated 403s
reference_cookie = os.environ.get('GDAL_HTTP_COOKIEFILE', '/tmp/cookies.txt')
secondary_cookie = f"/tmp/cookies_worker_{os.getpid()}.txt"

if os.path.exists(reference_cookie):
    try:
        shutil.copy(reference_cookie, secondary_cookie)
    except Exception:
        pass

os.environ['GDAL_HTTP_COOKIEFILE'] = secondary_cookie
os.environ['GDAL_HTTP_COOKIEJAR'] = secondary_cookie

# 2. Aggressive GDAL HTTP Retry settings (Catches AWS 503 SlowDown / 429 Rate Limits)
os.environ['GDAL_HTTP_MAX_RETRY'] = '10'
os.environ['GDAL_HTTP_RETRY_DELAY'] = '3'
os.environ['GDAL_HTTP_MULTIPLEX'] = 'YES'

# 3. Disable HDF5 locking and remote directory scanning
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['VSI_CACHE'] = 'YES'
# ---------------------------------------------------------

import ARIAtools.extractProduct


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'export_products_json', help='Export products json input file')
    args = parser.parse_args()

    with open(args.export_products_json) as ifp:
        export_product_args = json.load(ifp)

    ii, ilayer, outname, prod_arr = \
        ARIAtools.extractProduct.export_product_worker(*export_product_args)

    outputs = {
        'ii': ii, 'ilayer': ilayer, 'outname': outname, 'prod_arr': prod_arr}

    base_dir = os.path.abspath(os.path.dirname(args.export_products_json))
    outfile = os.path.join(base_dir, 'outputs_%d_%d.json' % (
        ii, ilayer))

    with open(outfile, 'w') as ofp:
        json.dump(outputs, ofp)

    # Flush standard streams and force a hard exit for the gnu_parallel worker.
    # This completely bypasses Python's noisy GDAL C-binding garbage collection phase.
    try:
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(0)
    except Exception:
        pass


if __name__ == "__main__":
    main()
