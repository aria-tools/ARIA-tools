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
# 1. Unique cookies per process to prevent Auth collisions
cookie_path = f"/tmp/cookies_worker_{os.getpid()}.txt"
os.environ['GDAL_HTTP_COOKIEFILE'] = cookie_path
os.environ['GDAL_HTTP_COOKIEJAR'] = cookie_path

# 2. Disable HDF5 locking and enable heavy caching
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['VSI_CACHE'] = 'YES'
os.environ['VSI_CACHE_SIZE'] = '536870912' # 512MB cache to prevent thrashing

# 3. CRITICAL: Disable HTTP Multiplexing and Range Merging
# This prevents GDAL from merging the losX and losY byte requests, 
# which corrupts the HDF5 C-library's strict byte parser.
os.environ['GDAL_HTTP_MULTIPLEX'] = 'NO'
os.environ['GDAL_HTTP_MERGE_CONSECUTIVE_RANGES'] = 'NO'
os.environ['CPL_VSIL_CURL_USE_HEAD'] = 'NO'
os.environ['GDAL_MAX_DATASET_POOL_SIZE'] = '0'

# 4. Aggressive Retries for AWS Rate Limits
os.environ['GDAL_HTTP_MAX_RETRY'] = '10'
os.environ['GDAL_HTTP_RETRY_DELAY'] = '3'
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
