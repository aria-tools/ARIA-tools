#!/usr/bin/env python3
"""
Benchmark: S3 direct access vs HTTPS (vsicurl) virtual processing.

Compares wall-clock times for the key ARIA-tools virtual processing steps
using ``/vsis3/`` (S3 direct, AWS only) vs ``/vsicurl/`` (HTTPS, anywhere).

Benchmark categories
--------------------
1. **GDAL Open**         – open a remote NetCDF/HDF5 with ``gdal.Open()``
2. **Metadata Extract**  – ``gdal.Info()`` metadata extraction (cache miss)
3. **Layer Read**        – read a coherence raster subset via GDAL
4. **Product Init**      – ``Product()`` full initialisation with bbox crop
5. **h5py Scalar Read**  – read an HDF5 scalar field via HTTP byte-range

Each operation is performed ``--reps`` times (default 3) and the median
is reported.  On AWS, both S3 and HTTPS paths are tested side-by-side;
off-AWS only HTTPS is available.

Prerequisites
-------------
- NASA Earthdata credentials in ``~/.netrc``
- ARIA-tools installed (``pip install -e .``)
- GDAL >= 3.8 with NetCDF support
- ``asf_search`` for product discovery
- ``boto3`` (for S3 path, optional)

Usage
-----
    # On AWS – full comparison:
    python tests/benchmark_s3_vs_https.py

    # Off-AWS – HTTPS-only baselines:
    python tests/benchmark_s3_vs_https.py

    # Quick run (1 rep, 2 products):
    python tests/benchmark_s3_vs_https.py --reps 1 --max-products 2

    # Specify products / custom dataset:
    python tests/benchmark_s3_vs_https.py --track 64 \\
        --start 20260101 --end 20260131

    # Use a pre-existing URL file instead of querying ASF:
    python tests/benchmark_s3_vs_https.py --url-file /path/to/urls.txt

    # Include download benchmark:
    python tests/benchmark_s3_vs_https.py --download

    # Save results to JSON:
    python tests/benchmark_s3_vs_https.py -o benchmark_results.json
"""
import argparse
import csv
import datetime
import json
import logging
import os
import shutil
import statistics
import sys
import tempfile
import time

# ---------------------------------------------------------------------------
# Ensure the tools/ tree is importable
# ---------------------------------------------------------------------------
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, 'tools'))

# ---------------------------------------------------------------------------
# Heavy imports are deferred to main() so --help works without deps
# ---------------------------------------------------------------------------


LOGGER = logging.getLogger('benchmark')

# ── Defaults ──────────────────────────────────────────────────────────────
DEFAULT_TRACK = 161
DEFAULT_BBOX_WKT = (
    'POLYGON((2.7291 50.7556,5.3028 50.7556,'
    '5.3028 51.4959,2.7291 51.4959,2.7291 50.7556))'
)
DEFAULT_BBOX_SNWE = '50.7556 51.4959 2.7291 5.3028'
DEFAULT_START = '20260201'
DEFAULT_END = '20260228'


# ══════════════════════════════════════════════════════════════════════════
#  Timing helpers
# ══════════════════════════════════════════════════════════════════════════
class Timer:
    """Reusable context-manager timer."""
    def __init__(self):
        self.elapsed = 0.0

    def __enter__(self):
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *_):
        self.elapsed = time.perf_counter() - self._t0


class BenchmarkResult:
    """Holds timings for a single benchmark across multiple repetitions."""
    def __init__(self, name, access_method):
        self.name = name
        self.access_method = access_method  # 'S3' or 'HTTPS'
        self.timings = []                   # list of floats (seconds)
        self.detail = ''                    # extra context

    @property
    def median(self):
        return statistics.median(self.timings) if self.timings else 0.0

    @property
    def mean(self):
        return statistics.mean(self.timings) if self.timings else 0.0

    @property
    def stdev(self):
        return statistics.stdev(self.timings) if len(self.timings) > 1 else 0.0

    @property
    def best(self):
        return min(self.timings) if self.timings else 0.0

    @property
    def worst(self):
        return max(self.timings) if self.timings else 0.0

    def to_dict(self):
        return {
            'name': self.name,
            'access_method': self.access_method,
            'median_s': round(self.median, 4),
            'mean_s': round(self.mean, 4),
            'stdev_s': round(self.stdev, 4),
            'best_s': round(self.best, 4),
            'worst_s': round(self.worst, 4),
            'reps': len(self.timings),
            'detail': self.detail,
        }


# ══════════════════════════════════════════════════════════════════════════
#  Product discovery
# ══════════════════════════════════════════════════════════════════════════
def discover_products(track, bbox_wkt, start_date, end_date, max_products):
    """Query ASF for GUNW products and return (https_urls, s3_urls)."""
    import asf_search

    start = datetime.datetime.strptime(start_date, '%Y%m%d')
    end = datetime.datetime.strptime(end_date, '%Y%m%d')

    LOGGER.info('Querying ASF – Track %d, %s to %s', track,
                start_date, end_date)
    scenes = asf_search.geo_search(
        collections=['C2859376221-ASF', 'C1261881077-ASF'],
        dataset=asf_search.constants.ARIA_S1_GUNW,
        processingLevel=asf_search.constants.GUNW_STD,
        relativeOrbit=[track],
        intersectsWith=bbox_wkt,
        start=start - datetime.timedelta(days=1),
        end=end + datetime.timedelta(days=1),
    )

    if not scenes:
        LOGGER.error('No products found. Check query parameters or '
                     'Earthdata credentials in ~/.netrc.')
        return [], []

    if max_products and len(scenes) > max_products:
        scenes = scenes[:max_products]

    https_urls = [s.properties['url'] for s in scenes]
    s3_urls = []
    for s in scenes:
        addl = s.properties.get('additionalUrls', [])
        s3 = next((u for u in addl if u.startswith('s3://')), None)
        s3_urls.append(s3)

    n_s3 = sum(1 for u in s3_urls if u)
    LOGGER.info('Found %d products (%d with S3 URLs)', len(https_urls), n_s3)
    return https_urls, s3_urls


def load_url_file(path):
    """Parse an existing URL file (1 or 2 columns) into lists."""
    https_urls, s3_urls = [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            h, s = s3_util.parse_url_line(line)
            https_urls.append(h)
            s3_urls.append(s)
    LOGGER.info('Loaded %d URLs from %s', len(https_urls), path)
    return https_urls, s3_urls


# ══════════════════════════════════════════════════════════════════════════
#  GDAL configuration helpers
# ══════════════════════════════════════════════════════════════════════════
def _configure_https():
    """Configure GDAL for HTTPS /vsicurl/ access."""
    from ARIAtools.product import _configure_gdal_virtual_access
    _configure_gdal_virtual_access()


def _configure_s3(s3_urls):
    """Configure GDAL for S3 direct access. Returns endpoint_key."""
    _configure_https()
    first_s3 = next((u for u in s3_urls if u), None)
    if first_s3 is None:
        return None
    endpoint_key = s3_util._endpoint_key_for_s3uri(first_s3)
    s3_util.configure_gdal_s3(endpoint_key)
    return endpoint_key


# ══════════════════════════════════════════════════════════════════════════
#  Individual benchmarks
# ══════════════════════════════════════════════════════════════════════════
def bench_gdal_open(url, access_method, reps):
    """Benchmark: gdal.Open() on a single product."""
    r = BenchmarkResult('gdal.Open', access_method)
    for i in range(reps):
        t = Timer()
        with t:
            ds = osgeo.gdal.Open(f'NETCDF:"{url}')
        r.timings.append(t.elapsed)
        if ds:
            n_sub = len([k for k in ds.GetMetadata('SUBDATASETS')
                         if 'NAME' in k])
            r.detail = f'subdatasets={n_sub}'
            ds = None
        else:
            r.detail = 'gdal.Open returned None'
            break
    return r


def bench_metadata_extract(url, workdir, access_method, reps):
    """Benchmark: full metadata extraction (cache miss each time)."""
    r = BenchmarkResult('metadata extract (gdal.Info)', access_method)
    for i in range(reps):
        cache_data = {}  # fresh cache each rep → forces cache miss
        t = Timer()
        with t:
            meta = meta_cache.get_or_extract(url, cache_data)
        r.timings.append(t.elapsed)
        if meta:
            r.detail = (
                f"version={meta.get('version')}, "
                f"subdatasets={len(meta.get('subdatasets', []))}")
        else:
            r.detail = 'metadata extraction returned None'
            break
    return r


def bench_layer_read(url, access_method, reps):
    """Benchmark: open coherence subdataset and read a 256x256 window."""
    r = BenchmarkResult('layer read (coherence 256x256)', access_method)

    # Discover coherence subdataset path
    ds = osgeo.gdal.Open(f'NETCDF:"{url}')
    if ds is None:
        r.detail = 'failed to open product'
        return r
    sds = ds.GetMetadata('SUBDATASETS')
    ds = None
    coh_path = None
    for k, v in sorted(sds.items()):
        if 'NAME' in k and 'coherence' in v:
            coh_path = v
            break
    if coh_path is None:
        r.detail = 'no coherence subdataset'
        return r

    for i in range(reps):
        t = Timer()
        with t:
            coh_ds = osgeo.gdal.Open(coh_path)
            if coh_ds:
                xsz = coh_ds.RasterXSize
                ysz = coh_ds.RasterYSize
                x_off = max(0, xsz // 2 - 128)
                y_off = max(0, ysz // 2 - 128)
                wx = min(256, xsz - x_off)
                wy = min(256, ysz - y_off)
                arr = coh_ds.GetRasterBand(1).ReadAsArray(
                    x_off, y_off, wx, wy)
                coh_ds = None
        r.timings.append(t.elapsed)
        if arr is not None:
            r.detail = (f'{wx}x{wy} window, '
                        f'range=[{np.nanmin(arr):.3f}, {np.nanmax(arr):.3f}]')

    return r


def bench_product_init(url_file, bbox, workdir, access_method, reps):
    """Benchmark: Product() initialisation (metadata, pairing, bbox)."""
    r = BenchmarkResult('Product() init', access_method)

    from ARIAtools.product import Product

    for i in range(reps):
        # Delete cache so each rep starts cold (measures full I/O)
        cache_file = meta_cache._cache_path(url_file)
        if cache_file and os.path.isfile(cache_file):
            os.remove(cache_file)
        t = Timer()
        with t:
            p = Product(
                filearg=url_file,
                bbox=bbox,
                workdir=workdir,
                num_threads=1,
                url_version=None,
                nc_version='1c',
            )
        r.timings.append(t.elapsed)
        n = len(p.products[0]) if len(p.products) == 2 else 0
        r.detail = f'{len(p.files)} files, {n} ifg groups'

    return r


def bench_product_init_cached(url_file, bbox, workdir, access_method, reps):
    """Benchmark: Product() init with warm cache."""
    r = BenchmarkResult('Product() init (cached)', access_method)

    from ARIAtools.product import Product

    # Warm the cache with a single init
    Product(
        filearg=url_file, bbox=bbox, workdir=workdir,
        num_threads=1, url_version=None, nc_version='1c')

    for i in range(reps):
        t = Timer()
        with t:
            p = Product(
                filearg=url_file, bbox=bbox, workdir=workdir,
                num_threads=1, url_version=None, nc_version='1c')
        r.timings.append(t.elapsed)
        n = len(p.products[0]) if len(p.products) == 2 else 0
        r.detail = f'{len(p.files)} files, {n} ifg groups (cached)'

    return r


def bench_h5py_open(https_url, access_method, reps):
    """Benchmark: open a remote HDF5/NetCDF via h5py HTTP byte-range.

    Uses the meta_cache ``open_gunw_h5`` helper which wraps h5py with
    an HTTP byte-range file-like object.  This measures the raw
    connection + HDF5 superblock read overhead.
    """
    r = BenchmarkResult('h5py remote open', access_method)
    vsi_url = f'/vsicurl/{https_url}'

    for i in range(reps):
        t = Timer()
        try:
            with t:
                h5f = meta_cache.open_gunw_h5(vsi_url)
                # Read the root group keys to force superblock parse
                keys = list(h5f.keys())
                h5f.close()
            r.timings.append(t.elapsed)
            r.detail = f'root keys: {keys}'
        except Exception as exc:
            r.detail = f'error: {exc}'
            LOGGER.warning('h5py open failed: %s', exc)
            break

    return r


def bench_download_single(https_url, s3_url, workdir, access_method,
                          endpoint_key, reps):
    """Benchmark: download a single product via S3 or HTTPS."""
    r = BenchmarkResult('single file download', access_method)
    fname = https_url.split('/')[-1]

    for i in range(reps):
        filepath = os.path.join(workdir, f'dl_bench_{i}_{fname}')
        t = Timer()
        try:
            with t:
                if access_method == 'S3' and s3_url:
                    client = s3_util.get_s3_client(
                        endpoint_key, max_pool_connections=10)
                    bucket, key = s3_util.parse_s3_uri(s3_url)
                    client.download_file(bucket, key, filepath)
                else:
                    import asf_search
                    session = asf_search.ASFSession()
                    resp = session.get(https_url, stream=True)
                    resp.raise_for_status()
                    with open(filepath, 'wb') as f:
                        for chunk in resp.iter_content(chunk_size=65536):
                            f.write(chunk)
            r.timings.append(t.elapsed)
            size_mb = os.path.getsize(filepath) / 1e6
            r.detail = f'{size_mb:.1f} MB'
        except Exception as exc:
            LOGGER.warning('Download bench failed (%s): %s',
                           access_method, exc)
            r.detail = f'error: {exc}'
        finally:
            if os.path.exists(filepath):
                os.remove(filepath)

    return r


# ══════════════════════════════════════════════════════════════════════════
#  Report formatting
# ══════════════════════════════════════════════════════════════════════════
def print_report(results, on_aws):
    """Pretty-print benchmark results as a table, with speedup column."""
    # Group by benchmark name
    by_name = {}
    for r in results:
        by_name.setdefault(r.name, {})[r.access_method] = r

    print()
    print('=' * 90)
    print(f'  ARIA-tools Virtual Access Benchmark')
    print(f'  On AWS: {on_aws}')
    print('=' * 90)

    hdr = (f'  {"Benchmark":<35s} {"Method":<7s} '
           f'{"Median":>8s} {"Best":>8s} {"Worst":>8s} '
           f'{"Speedup":>8s}  Detail')
    print(hdr)
    print('-' * 90)

    for name in dict.fromkeys(r.name for r in results):
        group = by_name[name]
        s3_r = group.get('S3')
        https_r = group.get('HTTPS')

        for method, r in sorted(group.items()):
            speedup = ''
            if s3_r and https_r and method == 'S3' and https_r.median > 0:
                sp = https_r.median / s3_r.median
                speedup = f'{sp:.1f}x'
            elif s3_r and https_r and method == 'HTTPS' and s3_r:
                speedup = '(baseline)'

            print(f'  {r.name:<35s} {r.access_method:<7s} '
                  f'{r.median:>7.3f}s {r.best:>7.3f}s {r.worst:>7.3f}s '
                  f'{speedup:>8s}  {r.detail}')
        print()

    print('=' * 90)


def save_results(results, filepath, on_aws, meta):
    """Save results to JSON for later analysis."""
    data = {
        'timestamp': datetime.datetime.utcnow().isoformat() + 'Z',
        'on_aws': on_aws,
        'metadata': meta,
        'benchmarks': [r.to_dict() for r in results],
    }
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)
    LOGGER.info('Results saved to %s', filepath)


def save_csv(results, filepath, on_aws):
    """Save results to CSV for easy import into spreadsheets."""
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'benchmark', 'method', 'median_s', 'mean_s', 'stdev_s',
            'best_s', 'worst_s', 'reps', 'on_aws', 'detail',
        ])
        for r in results:
            writer.writerow([
                r.name, r.access_method,
                f'{r.median:.4f}', f'{r.mean:.4f}', f'{r.stdev:.4f}',
                f'{r.best:.4f}', f'{r.worst:.4f}', len(r.timings),
                on_aws, r.detail,
            ])
    LOGGER.info('CSV saved to %s', filepath)


# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════
def main():
    parser = argparse.ArgumentParser(
        description='Benchmark S3 vs HTTPS virtual access for ARIA-tools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage\n-----')[0],
    )
    # Product discovery
    parser.add_argument('--url-file', default=None,
                        help='Pre-existing URL file (skip ASF query)')
    parser.add_argument('--track', type=int, default=DEFAULT_TRACK,
                        help=f'S1 track number (default: {DEFAULT_TRACK})')
    parser.add_argument('--bbox', default=DEFAULT_BBOX_WKT,
                        help='WKT polygon for ASF search')
    parser.add_argument('--bbox-snwe', default=DEFAULT_BBOX_SNWE,
                        help='SNWE bbox for Product class')
    parser.add_argument('--start', default=DEFAULT_START,
                        help=f'Start date YYYYMMDD (default: {DEFAULT_START})')
    parser.add_argument('--end', default=DEFAULT_END,
                        help=f'End date YYYYMMDD (default: {DEFAULT_END})')
    parser.add_argument('--max-products', type=int, default=4,
                        help='Max products to benchmark (default: 4)')

    # Benchmark control
    parser.add_argument('--reps', type=int, default=3,
                        help='Repetitions per benchmark (default: 3)')
    parser.add_argument('--download', action='store_true',
                        help='Include single-file download benchmark')
    parser.add_argument('--skip-product-init', action='store_true',
                        help='Skip Product() init benchmarks (slowest)')
    parser.add_argument('--force-s3', action='store_true',
                        help='Run S3 benchmarks even if is_on_aws() is False '
                             '(for testing with explicit AWS creds)')

    # Output
    parser.add_argument('-o', '--output', default=None,
                        help='Save results to JSON file')
    parser.add_argument('--csv', default=None,
                        help='Save results to CSV file')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Debug logging')
    parser.add_argument('--workdir', default=None,
                        help='Working directory (default: auto temp dir)')

    args = parser.parse_args()

    # Deferred heavy imports — after argparse so --help works without deps
    global np, osgeo, s3_util, meta_cache
    import numpy as np
    import osgeo.gdal
    import ARIAtools.util.s3 as s3_util
    import ARIAtools.util.meta_cache as meta_cache

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(name)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S')

    # ── Product discovery ─────────────────────────────────────────────
    if args.url_file:
        https_urls, s3_urls = load_url_file(args.url_file)
    else:
        https_urls, s3_urls = discover_products(
            args.track, args.bbox, args.start, args.end, args.max_products)

    if not https_urls:
        LOGGER.error('No products to benchmark. Exiting.')
        return 1

    # ── Environment detection ─────────────────────────────────────────
    on_aws = s3_util.is_on_aws()
    has_s3 = any(u is not None for u in s3_urls)
    run_s3 = (on_aws or args.force_s3) and has_s3

    LOGGER.info('Environment: on_aws=%s, has_s3_urls=%s, run_s3=%s',
                on_aws, has_s3, run_s3)

    # ── Working directory ─────────────────────────────────────────────
    cleanup = args.workdir is None
    workdir = args.workdir or tempfile.mkdtemp(prefix='aria_bench_')
    os.makedirs(workdir, exist_ok=True)
    LOGGER.info('Working directory: %s', workdir)

    # Prepare directories
    https_workdir = os.path.join(workdir, 'https')
    s3_workdir = os.path.join(workdir, 's3')
    os.makedirs(https_workdir, exist_ok=True)
    if run_s3:
        os.makedirs(s3_workdir, exist_ok=True)

    results = []
    endpoint_key = None

    try:
        # ── Configure GDAL ───────────────────────────────────────────
        _configure_https()

        # Build HTTPS paths
        test_https_url = https_urls[0]
        https_vsi = f'/vsicurl/{test_https_url}'

        # Build S3 paths if available
        s3_vsi = None
        if run_s3:
            endpoint_key = _configure_s3(s3_urls)
            test_s3_url = next((u for u in s3_urls if u), None)
            if test_s3_url:
                s3_vsi = s3_util.s3uri_to_vsis3(test_s3_url)

        # Write URL files for Product() init
        https_url_file = os.path.join(https_workdir, 'urls.txt')
        with open(https_url_file, 'w') as f:
            for h in https_urls:
                f.write(h + '\n')

        s3_url_file = None
        if run_s3:
            s3_url_file = os.path.join(s3_workdir, 'urls.txt')
            with open(s3_url_file, 'w') as f:
                for h, s in zip(https_urls, s3_urls):
                    line = f'{h},{s}' if s else h
                    f.write(line + '\n')

        # ══════════════════════════════════════════════════════════════
        #  Run benchmarks
        # ══════════════════════════════════════════════════════════════
        LOGGER.info('Starting benchmarks (%d reps each) ...', args.reps)

        # --- 1. gdal.Open ---
        LOGGER.info('[1/%d] gdal.Open ...', 5 + int(args.download))
        results.append(bench_gdal_open(
            https_vsi, 'HTTPS', args.reps))
        if run_s3 and s3_vsi:
            results.append(bench_gdal_open(
                s3_vsi, 'S3', args.reps))

        # --- 2. Metadata extract ---
        LOGGER.info('[2/%d] Metadata extract ...', 5 + int(args.download))
        results.append(bench_metadata_extract(
            https_vsi, https_workdir, 'HTTPS', args.reps))
        if run_s3 and s3_vsi:
            results.append(bench_metadata_extract(
                s3_vsi, s3_workdir, 'S3', args.reps))

        # --- 3. Layer read ---
        LOGGER.info('[3/%d] Layer read ...', 5 + int(args.download))
        results.append(bench_layer_read(
            https_vsi, 'HTTPS', args.reps))
        if run_s3 and s3_vsi:
            results.append(bench_layer_read(
                s3_vsi, 'S3', args.reps))

        # --- 4. Product() init ---
        if not args.skip_product_init:
            LOGGER.info('[4/%d] Product() init (cold) ...',
                        5 + int(args.download))
            results.append(bench_product_init(
                https_url_file, args.bbox_snwe, https_workdir,
                'HTTPS', args.reps))
            if run_s3 and s3_url_file:
                results.append(bench_product_init(
                    s3_url_file, args.bbox_snwe, s3_workdir,
                    'S3', args.reps))

            LOGGER.info('[5/%d] Product() init (cached) ...',
                        5 + int(args.download))
            results.append(bench_product_init_cached(
                https_url_file, args.bbox_snwe, https_workdir,
                'HTTPS', args.reps))
            if run_s3 and s3_url_file:
                results.append(bench_product_init_cached(
                    s3_url_file, args.bbox_snwe, s3_workdir,
                    'S3', args.reps))

        # --- 5. h5py remote open ---
        LOGGER.info('[%d/%d] h5py remote open ...',
                    (4 if args.skip_product_init else 6),
                    5 + int(args.download))
        results.append(bench_h5py_open(
            test_https_url, 'HTTPS', args.reps))
        # h5py always uses HTTPS internally (even on S3), so no S3 bench

        # --- 6. Download (optional) ---
        if args.download:
            LOGGER.info('[%d/%d] Single-file download ...',
                        (5 if args.skip_product_init else 7),
                        5 + int(args.download))
            dl_dir = os.path.join(workdir, 'dl_bench')
            os.makedirs(dl_dir, exist_ok=True)
            results.append(bench_download_single(
                test_https_url, s3_urls[0], dl_dir, 'HTTPS',
                endpoint_key, 1))  # download only 1 rep — it's big
            if run_s3 and s3_urls[0]:
                results.append(bench_download_single(
                    test_https_url, s3_urls[0], dl_dir, 'S3',
                    endpoint_key, 1))

        # ══════════════════════════════════════════════════════════════
        #  Results
        # ══════════════════════════════════════════════════════════════
        print_report(results, on_aws)

        # Metadata for output file
        run_meta = {
            'track': args.track,
            'bbox_wkt': args.bbox,
            'bbox_snwe': args.bbox_snwe,
            'start_date': args.start,
            'end_date': args.end,
            'n_products': len(https_urls),
            'reps': args.reps,
        }

        if args.output:
            save_results(results, args.output, on_aws, run_meta)
        if args.csv:
            save_csv(results, args.csv, on_aws)

    finally:
        if cleanup:
            LOGGER.info('Cleaning up %s', workdir)
            shutil.rmtree(workdir, ignore_errors=True)

    return 0


if __name__ == '__main__':
    sys.exit(main())
