#!/usr/bin/env python3
"""
Test script for ARIA-tools virtual (vsicurl) data access.

Exercises the virtual access workflow end-to-end:
  1. Query ASF for Sentinel-1 GUNW products (Track 161, Feb 2026)
  2. Write product URLs (no download)
  3. Test metadata extraction and caching
  4. Test virtual bounding box extraction with a small spatial crop
  5. Verify that subsetting via vsicurl returns valid data

Prerequisites
-------------
- NASA Earthdata credentials in ~/.netrc
- ARIA-tools installed (pip install -e .)
- GDAL >= 3.8 with NetCDF support, libnetcdf >= 4.5

Usage
-----
    python tests/test_virtual_access.py
    python tests/test_virtual_access.py -v          # verbose / debug
    python tests/test_virtual_access.py --quick      # skip extract, metadata only
"""
import argparse
import datetime
import json
import logging
import os
import shutil
import sys
import tempfile
import time

# ---------------------------------------------------------------------------
# Ensure the tools/ tree is on the path so we can import ARIAtools even when
# the package is not installed via pip.
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'tools'))

import glob
import subprocess

import osgeo.gdal
import numpy as np

import ARIAtools.util.meta_cache as meta_cache

LOGGER = logging.getLogger('test_virtual_access')

# ── Test configuration ────────────────────────────────────────────────────
TRACK = 161
# BBox around Belgium / Netherlands area – Track 161, frame ~25050
# WKT polygon provided for ASF search
BBOX_WKT = 'POLYGON((2.7291 50.7556,5.3028 50.7556,5.3028 51.4959,2.7291 51.4959,2.7291 50.7556))'
# SNWE for Product class bbox cropping (derived from the polygon above)
BBOX_SNWE = '50.7556 51.4959 2.7291 5.3028'
START_DATE = '20260201'
END_DATE = '20260228'


# ── Helpers ───────────────────────────────────────────────────────────────
def query_urls(workdir):
    """Use ariaDownload machinery to query ASF and return product URLs."""
    import asf_search

    start = datetime.datetime.strptime(START_DATE, '%Y%m%d')
    end = datetime.datetime.strptime(END_DATE, '%Y%m%d')

    LOGGER.info('Querying ASF – Track %d, %s to %s, bbox=%s',
                TRACK, START_DATE, END_DATE, BBOX_WKT)

    scenes = asf_search.geo_search(
        collections=['C2859376221-ASF', 'C1261881077-ASF'],
        dataset=asf_search.constants.ARIA_S1_GUNW,
        processingLevel=asf_search.constants.GUNW_STD,
        relativeOrbit=[TRACK],
        intersectsWith=BBOX_WKT,
        start=start - datetime.timedelta(days=1),
        end=end + datetime.timedelta(days=1),
    )

    urls = [s.properties['url'] for s in scenes]
    LOGGER.info('ASF returned %d products', len(urls))

    # Log first few URLs for debugging
    for u in urls[:3]:
        LOGGER.debug('  %s', u)

    if not urls:
        LOGGER.warning(
            'No products found for Track %d in %s–%s. '
            'Try widening the date range or checking Earthdata credentials.',
            TRACK, START_DATE, END_DATE)
    return urls


def write_url_file(urls, workdir):
    """Write URLs to a text file, return path."""
    url_file = os.path.join(workdir, f'track{TRACK}_urls.txt')
    with open(url_file, 'w') as fh:
        for u in urls:
            fh.write(u + '\n')
    LOGGER.info('Wrote %d URLs to %s', len(urls), url_file)
    return url_file


# ── Test functions ────────────────────────────────────────────────────────
class Result:
    def __init__(self, name):
        self.name = name
        self.passed = False
        self.message = ''
        self.elapsed = 0.0

    def __repr__(self):
        status = 'PASS' if self.passed else 'FAIL'
        return f'[{status}] {self.name} ({self.elapsed:.1f}s) {self.message}'


def test_query(workdir):
    """Test 1: Query ASF for product URLs."""
    r = Result('ASF query')
    t0 = time.time()
    try:
        urls = query_urls(workdir)
        r.elapsed = time.time() - t0
        if len(urls) > 0:
            r.passed = True
            r.message = f'{len(urls)} products found'
        else:
            r.message = 'No products found'
    except Exception as exc:
        r.elapsed = time.time() - t0
        r.message = str(exc)
    return r, urls if r.passed else []


def test_gdal_config():
    """Test 2: Verify GDAL virtual access configuration."""
    from ARIAtools.product import _configure_gdal_virtual_access
    r = Result('GDAL config')
    t0 = time.time()
    try:
        _configure_gdal_virtual_access()
        checks = {
            'VSI_CACHE': 'YES',
            'VSI_CACHE_SIZE': '67108864',
            'CPL_VSIL_CURL_CHUNK_SIZE': '524288',
            'GDAL_HTTP_MAX_RETRY': '3',
            'GDAL_HTTP_MERGE_CONSECUTIVE_RANGES': 'YES',
            'GDAL_HTTP_COOKIEFILE': '/tmp/cookies.txt',
            'GDAL_HTTP_COOKIEJAR': '/tmp/cookies.txt',
        }
        failed = []
        for key, expected in checks.items():
            actual = osgeo.gdal.GetConfigOption(key)
            if actual != expected:
                failed.append(f'{key}: expected={expected}, got={actual}')
        r.elapsed = time.time() - t0
        if not failed:
            r.passed = True
            r.message = 'All GDAL config options set correctly'
        else:
            r.message = '; '.join(failed)
    except Exception as exc:
        r.elapsed = time.time() - t0
        r.message = str(exc)
    return r


def test_metadata_cache(urls, workdir):
    """Test 3: Extract metadata via cache, verify cache hit on second call."""
    r = Result('Metadata cache')
    t0 = time.time()
    if not urls:
        r.message = 'Skipped – no URLs'
        r.elapsed = time.time() - t0
        return r

    url_file = write_url_file(urls, workdir)
    cache_file = meta_cache._cache_path(url_file)
    cache_data = {}

    try:
        # First URL – should be cache miss
        test_url = f'/vsicurl/{urls[0]}'

        # Configure GDAL for virtual access
        from ARIAtools.product import _configure_gdal_virtual_access
        _configure_gdal_virtual_access()

        t_miss = time.time()
        meta = meta_cache.get_or_extract(test_url, cache_data)
        dt_miss = time.time() - t_miss

        if meta is None:
            r.message = 'Failed to extract metadata from first URL'
            r.elapsed = time.time() - t0
            return r

        # Validate extracted metadata
        version = meta.get('version')
        driver = meta.get('driver')
        subdatasets = meta.get('subdatasets', [])
        LOGGER.info('  Version: %s, Driver: %s, Subdatasets: %d',
                     version, driver, len(subdatasets))

        # Save then reload cache
        meta_cache.save_cache(cache_file, cache_data)
        cache_data_reloaded = meta_cache.load_cache(cache_file)

        # Second call should be a cache hit (instant)
        t_hit = time.time()
        meta2 = meta_cache.get_or_extract(test_url, cache_data_reloaded)
        dt_hit = time.time() - t_hit

        LOGGER.info('  Cache miss: %.2fs, cache hit: %.4fs', dt_miss, dt_hit)

        if meta2 is not None and meta2.get('version') == version:
            r.passed = True
            r.message = (
                f'version={version}, {len(subdatasets)} subdatasets, '
                f'miss={dt_miss:.1f}s, hit={dt_hit:.4f}s')
        else:
            r.message = 'Cache reload returned inconsistent data'

    except Exception as exc:
        r.message = str(exc)

    r.elapsed = time.time() - t0
    return r


def test_virtual_open_single(urls):
    """Test 4: Open a single product via vsicurl and read version."""
    r = Result('Virtual open (single)')
    t0 = time.time()
    if not urls:
        r.message = 'Skipped – no URLs'
        r.elapsed = time.time() - t0
        return r

    try:
        from ARIAtools.product import _configure_gdal_virtual_access
        _configure_gdal_virtual_access()

        vsi_path = f'/vsicurl/{urls[0]}'
        nc_path = f'NETCDF:"{vsi_path}'
        ds = osgeo.gdal.Open(nc_path)
        if ds is None:
            r.message = 'gdal.Open returned None'
            r.elapsed = time.time() - t0
            return r

        version = ds.GetMetadataItem('NC_GLOBAL#version')
        driver = ds.GetDriver().GetDescription()
        n_sub = len([k for k in ds.GetMetadata('SUBDATASETS')
                     if 'NAME' in k])
        ds = None

        r.passed = True
        r.message = f'driver={driver}, version={version}, subdatasets={n_sub}'
    except Exception as exc:
        r.message = str(exc)

    r.elapsed = time.time() - t0
    return r


def test_virtual_read_layer(urls):
    """Test 5: Read a small subset of coherence via vsicurl."""
    r = Result('Virtual layer read (coherence)')
    t0 = time.time()
    if not urls:
        r.message = 'Skipped – no URLs'
        r.elapsed = time.time() - t0
        return r

    try:
        from ARIAtools.product import _configure_gdal_virtual_access
        _configure_gdal_virtual_access()

        vsi_path = f'/vsicurl/{urls[0]}'
        nc_path = f'NETCDF:"{vsi_path}'

        ds = osgeo.gdal.Open(nc_path)
        sds = ds.GetMetadata('SUBDATASETS')
        ds = None

        # Find the coherence subdataset
        coh_sds = None
        for k, v in sorted(sds.items()):
            if 'NAME' in k and 'coherence' in v:
                coh_sds = v
                break

        if coh_sds is None:
            r.message = 'No coherence subdataset found'
            r.elapsed = time.time() - t0
            return r

        LOGGER.info('  Opening coherence: %s', coh_sds)
        coh_ds = osgeo.gdal.Open(coh_sds)
        if coh_ds is None:
            r.message = 'Failed to open coherence subdataset'
            r.elapsed = time.time() - t0
            return r

        xsize = coh_ds.RasterXSize
        ysize = coh_ds.RasterYSize
        LOGGER.info('  Coherence raster size: %d x %d', xsize, ysize)

        # Read a small 100x100 window from the center
        x_off = max(0, xsize // 2 - 50)
        y_off = max(0, ysize // 2 - 50)
        win_x = min(100, xsize - x_off)
        win_y = min(100, ysize - y_off)

        band = coh_ds.GetRasterBand(1)
        arr = band.ReadAsArray(x_off, y_off, win_x, win_y)
        coh_ds = None

        valid_pixels = np.count_nonzero(~np.isnan(arr.astype(float)))
        r.passed = arr is not None and arr.size > 0
        r.message = (
            f'Read {win_x}x{win_y} window, '
            f'min={np.nanmin(arr):.3f}, max={np.nanmax(arr):.3f}, '
            f'valid_pixels={valid_pixels}/{arr.size}')
    except Exception as exc:
        r.message = str(exc)

    r.elapsed = time.time() - t0
    return r


def test_virtual_bbox_subset(urls, workdir):
    """Test 6: Use Product class with bbox to verify spatial subsetting."""
    r = Result('Product class with bbox crop')
    t0 = time.time()
    if not urls:
        r.message = 'Skipped – no URLs'
        r.elapsed = time.time() - t0
        return r

    # Write URL file
    url_file = write_url_file(urls, workdir)

    try:
        from ARIAtools.product import Product

        LOGGER.info('  Initializing Product with bbox=%s', BBOX_SNWE)
        p = Product(
            filearg=url_file,
            bbox=BBOX_SNWE,
            workdir=workdir,
            num_threads=1,
            url_version=None,
            nc_version='1c',
        )
        n_products = len(p.products[0]) if len(p.products) == 2 else 0
        n_files = len(p.files)

        LOGGER.info('  Files: %d, IFG groups: %d', n_files, n_products)

        if n_products > 0:
            # After __continuous_time__, products is [list_of_rmd, list_of_lyr]
            rmd_list = p.products[0]
            lyr_list = p.products[1]
            sample_rmd = rmd_list[0]
            expected_keys = {'pair_name', 'missionID', 'wavelength'}
            found_keys = expected_keys.intersection(sample_rmd.keys())
            r.passed = len(found_keys) == len(expected_keys)
            r.message = (
                f'{len(rmd_list)} ifg groups from {n_files} files, '
                f'metadata keys present: {sorted(found_keys)}')
        else:
            r.message = f'{n_files} files loaded but 0 products matched bbox'

        # Check that cache file was created
        cache_file = meta_cache._cache_path(url_file)
        if cache_file and os.path.isfile(cache_file):
            cache_size = os.path.getsize(cache_file)
            r.message += f', cache={cache_size}B'

    except Exception as exc:
        import traceback
        r.message = f'{exc}\n{traceback.format_exc()}'

    r.elapsed = time.time() - t0
    return r


def test_cache_speedup(urls, workdir):
    """Test 7: Measure speedup from metadata cache on second Product init."""
    r = Result('Cache speedup')
    t0 = time.time()
    if not urls:
        r.message = 'Skipped – no URLs'
        r.elapsed = time.time() - t0
        return r

    url_file = os.path.join(workdir, f'track{TRACK}_urls.txt')
    if not os.path.isfile(url_file):
        url_file = write_url_file(urls, workdir)

    try:
        from ARIAtools.product import Product

        # First run (cache cold or partially warm from test 6)
        # Delete cache to ensure cold start
        cache_file = meta_cache._cache_path(url_file)
        if cache_file and os.path.isfile(cache_file):
            os.remove(cache_file)

        t1 = time.time()
        p1 = Product(
            filearg=url_file, bbox=BBOX_SNWE,
            workdir=workdir, num_threads=1, url_version=None,
            nc_version='1c')
        dt_cold = time.time() - t1

        # Second run (cache warm)
        t2 = time.time()
        p2 = Product(
            filearg=url_file, bbox=BBOX_SNWE,
            workdir=workdir, num_threads=1, url_version=None,
            nc_version='1c')
        dt_warm = time.time() - t2

        speedup = dt_cold / dt_warm if dt_warm > 0 else float('inf')
        r.passed = True
        r.message = (
            f'cold={dt_cold:.1f}s, warm={dt_warm:.1f}s, '
            f'speedup={speedup:.1f}x')

    except Exception as exc:
        import traceback
        r.message = f'{exc}\n{traceback.format_exc()}'

    r.elapsed = time.time() - t0
    return r


# ── Download helper ───────────────────────────────────────────────────────
def download_products(urls, outdir):
    """Download GUNW products from ASF using Earthdata credentials.

    Uses asf_search.ASFSession which handles Earthdata OAuth
    transparently (reads ~/.netrc).
    """
    import asf_search
    session = asf_search.ASFSession()
    paths = []
    for url in urls:
        fname = url.split('/')[-1]
        outpath = os.path.join(outdir, fname)
        if os.path.exists(outpath):
            LOGGER.info('  Already downloaded: %s', fname)
            paths.append(outpath)
            continue
        LOGGER.info('  Downloading %s ...', fname)
        resp = session.get(url, stream=True)
        resp.raise_for_status()
        with open(outpath, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=65536):
                f.write(chunk)
        size_mb = os.path.getsize(outpath) / 1e6
        LOGGER.info('  Saved %s (%.1f MB)', fname, size_mb)
        paths.append(outpath)
    return paths


def test_virtual_vs_local_extract(urls, workdir):
    """Test 8: Download products and compare virtual vs local ariaExtract.

    Downloads 2 products, runs ariaExtract -l coherence on both
    local files and virtual URLs, then compares output rasters
    pixel-by-pixel to verify they are identical.
    """
    r = Result('Virtual vs local extract comparison')
    t0 = time.time()

    if len(urls) < 2:
        r.message = 'Skipped – need at least 2 URLs'
        r.elapsed = time.time() - t0
        return r

    # Use exactly 2 products for the comparison
    test_urls = urls[:2]

    local_dl_dir = os.path.join(workdir, 'downloaded_products')
    local_extract_dir = os.path.join(workdir, 'extract_local')
    virtual_extract_dir = os.path.join(workdir, 'extract_virtual')
    os.makedirs(local_dl_dir, exist_ok=True)
    os.makedirs(local_extract_dir, exist_ok=True)
    os.makedirs(virtual_extract_dir, exist_ok=True)

    try:
        # Step 1 – Download products
        LOGGER.info('Step 1: Downloading %d products ...', len(test_urls))
        local_files = download_products(test_urls, local_dl_dir)
        if len(local_files) != len(test_urls):
            r.message = (
                f'Download incomplete: {len(local_files)}/{len(test_urls)}')
            r.elapsed = time.time() - t0
            return r

        # Step 2 – ariaExtract on LOCAL files
        LOGGER.info('Step 2: Running ariaExtract on local files ...')
        local_glob = os.path.join(local_dl_dir, '*.nc')
        cmd_local = (
            f'ariaExtract.py '
            f'-f "{local_glob}" '
            f'-l coherence '
            f'-b "{BBOX_SNWE}" '
            f'-w {local_extract_dir} '
            f'-of ENVI '
            f'--nc_version 1c '
            f'--log-level info'
        )
        LOGGER.info('  cmd: %s', cmd_local)
        cp_local = subprocess.run(
            cmd_local, shell=True,
            capture_output=True, text=True, timeout=600)
        if cp_local.returncode != 0:
            r.message = (
                f'Local extract failed (rc={cp_local.returncode}): '
                f'{cp_local.stderr[-500:]}')
            r.elapsed = time.time() - t0
            return r

        # Step 3 – ariaExtract on VIRTUAL URLs
        LOGGER.info('Step 3: Running ariaExtract on virtual URLs ...')
        url_file = os.path.join(workdir, 'compare_urls.txt')
        with open(url_file, 'w') as fh:
            for u in test_urls:
                fh.write(u + '\n')

        cmd_virtual = (
            f'ariaExtract.py '
            f'-f {url_file} '
            f'-l coherence '
            f'-b "{BBOX_SNWE}" '
            f'-w {virtual_extract_dir} '
            f'-of ENVI '
            f'--nc_version 1c '
            f'--log-level info'
        )
        LOGGER.info('  cmd: %s', cmd_virtual)
        cp_virtual = subprocess.run(
            cmd_virtual, shell=True,
            capture_output=True, text=True, timeout=600)
        if cp_virtual.returncode != 0:
            r.message = (
                f'Virtual extract failed (rc={cp_virtual.returncode}): '
                f'{cp_virtual.stderr[-500:]}')
            r.elapsed = time.time() - t0
            return r

        # Step 4 – Compare outputs
        LOGGER.info('Step 4: Comparing extracted rasters ...')
        local_coh_dir = os.path.join(local_extract_dir, 'coherence')
        virtual_coh_dir = os.path.join(virtual_extract_dir, 'coherence')

        # Find output VRT files (ENVI also produces .vrt alongside binary)
        local_vrts = sorted(glob.glob(
            os.path.join(local_coh_dir, '*.vrt')))
        virtual_vrts = sorted(glob.glob(
            os.path.join(virtual_coh_dir, '*.vrt')))

        if not local_vrts:
            r.message = (
                f'No local coherence outputs in {local_coh_dir}')
            r.elapsed = time.time() - t0
            return r
        if not virtual_vrts:
            r.message = (
                f'No virtual coherence outputs in {virtual_coh_dir}')
            r.elapsed = time.time() - t0
            return r

        # Match files by basename
        local_by_name = {os.path.basename(f): f for f in local_vrts}
        virtual_by_name = {os.path.basename(f): f for f in virtual_vrts}
        common = sorted(
            set(local_by_name.keys()) & set(virtual_by_name.keys()))

        if not common:
            r.message = (
                f'No matching outputs: '
                f'local={sorted(local_by_name.keys())}, '
                f'virtual={sorted(virtual_by_name.keys())}')
            r.elapsed = time.time() - t0
            return r

        LOGGER.info('  Comparing %d coherence pairs: %s', len(common), common)

        all_match = True
        match_details = []
        for name in common:
            ds_local = osgeo.gdal.Open(local_by_name[name])
            ds_virtual = osgeo.gdal.Open(virtual_by_name[name])

            if ds_local is None or ds_virtual is None:
                match_details.append(f'{name}: GDAL open failed')
                all_match = False
                continue

            # Check dimensions
            if (ds_local.RasterXSize != ds_virtual.RasterXSize or
                    ds_local.RasterYSize != ds_virtual.RasterYSize):
                match_details.append(
                    f'{name}: size mismatch '
                    f'({ds_local.RasterXSize}x{ds_local.RasterYSize} vs '
                    f'{ds_virtual.RasterXSize}x{ds_virtual.RasterYSize})')
                all_match = False
                ds_local = ds_virtual = None
                continue

            # Check geotransform
            gt_local = ds_local.GetGeoTransform()
            gt_virtual = ds_virtual.GetGeoTransform()
            if gt_local != gt_virtual:
                match_details.append(
                    f'{name}: geotransform mismatch')
                all_match = False

            # Compare raster data band by band
            for band_i in range(1, ds_local.RasterCount + 1):
                arr_l = ds_local.GetRasterBand(band_i).ReadAsArray()
                arr_v = ds_virtual.GetRasterBand(band_i).ReadAsArray()

                # Handle NaN: treat as equal
                nan_l = np.isnan(arr_l)
                nan_v = np.isnan(arr_v)
                if not np.array_equal(nan_l, nan_v):
                    match_details.append(
                        f'{name} band{band_i}: NaN pattern differs')
                    all_match = False
                    continue

                # Compare non-NaN values
                valid = ~nan_l
                if valid.any():
                    max_diff = np.max(
                        np.abs(arr_l[valid] - arr_v[valid]))
                    if max_diff > 0:
                        match_details.append(
                            f'{name} band{band_i}: max_diff={max_diff:.2e}')
                        all_match = False
                    else:
                        match_details.append(
                            f'{name} band{band_i}: identical '
                            f'({np.count_nonzero(valid)} pixels)')
                else:
                    match_details.append(
                        f'{name} band{band_i}: all NaN')

            ds_local = ds_virtual = None

        r.passed = all_match
        r.message = '; '.join(match_details)

    except subprocess.TimeoutExpired:
        r.message = 'ariaExtract subprocess timed out (600s)'
    except Exception as exc:
        import traceback
        r.message = f'{exc}\n{traceback.format_exc()}'

    r.elapsed = time.time() - t0
    return r


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='Test ARIA-tools virtual data access')
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Enable debug logging')
    parser.add_argument(
        '--quick', action='store_true',
        help='Quick mode: metadata tests only, skip layer reads')
    parser.add_argument(
        '--max-products', type=int, default=4,
        help='Max products to use for testing (default: 4)')
    parser.add_argument(
        '--compare', action='store_true',
        help='Run download + virtual-vs-local comparison test (slow)')
    parser.add_argument(
        '--workdir', default=None,
        help='Working directory (default: temp dir, auto-cleaned)')
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S')

    # Working directory
    cleanup = args.workdir is None
    workdir = args.workdir or tempfile.mkdtemp(prefix='aria_vtest_')
    os.makedirs(workdir, exist_ok=True)
    LOGGER.info('Working directory: %s', workdir)

    results = []
    urls = []

    try:
        # Test 1: ASF query
        r_query, urls = test_query(workdir)
        results.append(r_query)

        # Limit products for faster testing
        if urls and args.max_products:
            urls = urls[:args.max_products]
            LOGGER.info('Limited to %d products for testing', len(urls))

        # Test 2: GDAL config
        results.append(test_gdal_config())

        if urls:
            # Test 3: Metadata cache
            results.append(test_metadata_cache(urls, workdir))

            # Test 4: Virtual open
            results.append(test_virtual_open_single(urls))

            if not args.quick:
                # Test 5: Virtual layer read
                results.append(test_virtual_read_layer(urls))

                # Test 6: Product class with bbox crop
                results.append(test_virtual_bbox_subset(urls, workdir))

                # Test 7: Cache speedup measurement
                results.append(test_cache_speedup(urls, workdir))

                # Test 8: Virtual vs local comparison (download-heavy)
                if args.compare:
                    results.append(
                        test_virtual_vs_local_extract(urls, workdir))
        else:
            LOGGER.warning(
                'No products found – skipping remote access tests. '
                'Check ~/.netrc or widen date range.')

    finally:
        # Summary
        print('\n' + '=' * 70)
        print('VIRTUAL ACCESS TEST RESULTS')
        print('=' * 70)
        n_pass = 0
        n_fail = 0
        for r in results:
            print(r)
            if r.passed:
                n_pass += 1
            else:
                n_fail += 1
        print('-' * 70)
        print(f'Total: {n_pass} passed, {n_fail} failed '
              f'out of {len(results)} tests')
        print('=' * 70)

        if cleanup:
            LOGGER.info('Cleaning up %s', workdir)
            shutil.rmtree(workdir, ignore_errors=True)

    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
