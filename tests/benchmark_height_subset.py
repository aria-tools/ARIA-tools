#!/usr/bin/env python
"""End-to-end benchmark for the 3D cube height-subsetting optimisation.

Runs ariaExtract.py twice for each requested layer — once with the height
subsetting enabled (default) and once with it disabled — and compares
wall-clock time, output accuracy, and reports DEM/access details.

Usage:
  # Benchmark a single layer with a local GUNW + DEM
  python tests/benchmark_height_subset.py \\
      --file "products/*.nc" \\
      --dem DEM/glo_90.dem \\
      --layers incidenceAngle

  # Multiple layers, custom bbox
  python tests/benchmark_height_subset.py \\
      --file "products/*.nc" \\
      --dem DEM/glo_90.dem \\
      --layers incidenceAngle,lookAngle,azimuthAngle \\
      --bbox '34.5 35.0 -118.5 -117.5'

  # Test with S3 virtual access (vsicurl/vsis3)
  python tests/benchmark_height_subset.py \\
      --file urls.txt \\
      --dem DEM/glo_90.dem \\
      --layers incidenceAngle
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
import osgeo.gdal


def _detect_access_mode(file_arg):
    """Determine whether the GUNW input is local, vsicurl, or vsis3."""
    # If it's a text file, peek inside for URLs
    expanded = glob.glob(file_arg) if '*' in file_arg else [file_arg]
    for f in expanded:
        if f.startswith('/vsis3/') or f.startswith('s3://'):
            return 'vsis3'
        if f.startswith('/vsicurl/') or f.startswith('http'):
            return 'vsicurl'
        if f.endswith('.txt') and os.path.isfile(f):
            with open(f) as fh:
                first = fh.readline().strip()
            if first.startswith('s3://') or first.startswith('/vsis3/'):
                return 'vsis3'
            if first.startswith('http'):
                return 'vsicurl'
    return 'local'


def _get_dem_range(demfile):
    """Return (min, max) elevation from the DEM, skipping nodata."""
    ds = osgeo.gdal.Open(demfile)
    if ds is None:
        return None, None
    from ARIAtools.util.interp import _compute_dem_range
    try:
        vmin, vmax = _compute_dem_range(ds)
    except ValueError:
        vmin, vmax = None, None
    ds = None
    return vmin, vmax


def _get_height_levels(file_arg):
    """Read height levels from the first GUNW product found."""
    try:
        import h5py
        expanded = glob.glob(file_arg) if '*' in file_arg else [file_arg]
        # If text file with URLs, skip (can't h5py over network trivially)
        for f in expanded:
            if f.endswith('.nc') and os.path.isfile(f):
                with h5py.File(f, 'r') as h5:
                    base = '/science/grids/imagingGeometry'
                    if base + '/heightsMeta' in h5:
                        return h5[base + '/heightsMeta'][:].tolist()
    except Exception:
        pass
    return None


def _compare_outputs(dir_a, dir_b, layer):
    """Compare raster outputs between two workdirs for a given layer."""
    # Find output files
    pattern_a = os.path.join(dir_a, layer, '**', '*')
    pattern_b = os.path.join(dir_b, layer, '**', '*')
    files_a = sorted(
        f for f in glob.glob(pattern_a, recursive=True)
        if os.path.isfile(f) and not f.endswith('.vrt')
        and not f.endswith('.xml') and not f.endswith('.hdr'))
    files_b = sorted(
        f for f in glob.glob(pattern_b, recursive=True)
        if os.path.isfile(f) and not f.endswith('.vrt')
        and not f.endswith('.xml') and not f.endswith('.hdr'))

    if not files_a or not files_b:
        return {'n_files': 0, 'max_diff': float('nan'),
                'mean_diff': float('nan')}

    max_diff = 0.0
    diffs = []
    for fa, fb in zip(files_a, files_b):
        ds_a = osgeo.gdal.Open(fa)
        ds_b = osgeo.gdal.Open(fb)
        if ds_a is None or ds_b is None:
            continue
        arr_a = ds_a.ReadAsArray().astype('float64')
        arr_b = ds_b.ReadAsArray().astype('float64')
        ds_a = ds_b = None

        if arr_a.shape != arr_b.shape:
            diffs.append(float('nan'))
            continue

        mask = np.isfinite(arr_a) & np.isfinite(arr_b)
        if mask.any():
            d = np.max(np.abs(arr_a[mask] - arr_b[mask]))
            max_diff = max(max_diff, d)
            diffs.append(d)

    return {
        'n_files': len(files_a),
        'max_diff': max_diff,
        'mean_diff': float(np.nanmean(diffs)) if diffs else float('nan'),
    }


def _run_extract(file_arg, dem, bbox, layer, workdir, disable_subset,
                 force_https=False, extra_args=None):
    """Run ariaExtract.py and return wall-clock time."""
    cmd = [
        sys.executable, '-m', 'tools.bin.ariaExtract',
        '-f', file_arg,
        '-d', dem,
        '-l', layer,
        '-w', workdir,
        '-of', 'ENVI',
        '--log-level', 'warning',
    ]
    if bbox:
        cmd.extend(['-b', bbox])
    if extra_args:
        cmd.extend(extra_args)

    env = os.environ.copy()
    if disable_subset:
        env['ARIA_DISABLE_HEIGHT_SUBSET'] = '1'
    else:
        env.pop('ARIA_DISABLE_HEIGHT_SUBSET', None)
    if force_https:
        env['ARIA_FORCE_HTTPS'] = '1'
    else:
        env.pop('ARIA_FORCE_HTTPS', None)

    t0 = time.perf_counter()
    result = subprocess.run(
        cmd, env=env, capture_output=True, text=True,
        cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    elapsed = time.perf_counter() - t0

    return elapsed, result.returncode, result.stdout, result.stderr


def main():
    parser = argparse.ArgumentParser(
        description='End-to-end benchmark for height-subsetting optimisation.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    parser.add_argument(
        '-f', '--file', required=True,
        help='GUNW products (glob or txt file, same as ariaExtract.py -f)')
    parser.add_argument(
        '-d', '--dem', required=True,
        help='DEM file (same as ariaExtract.py -d)')
    parser.add_argument(
        '-l', '--layers', default='incidenceAngle',
        help='Comma-separated layers to benchmark (default: incidenceAngle)')
    parser.add_argument(
        '-b', '--bbox', default=None,
        help="Bounding box 'S N W E' (default: use product extent)")
    parser.add_argument(
        '--keep-outputs', action='store_true',
        help='Keep the temporary output directories for inspection')
    parser.add_argument(
        '--compare-access', action='store_true',
        help='When using URLs, run both S3 and vsicurl (HTTPS) access '
             'modes to compare performance. Only meaningful on AWS with '
             'S3-capable URL files.')
    args = parser.parse_args()

    layers = [l.strip() for l in args.layers.split(',')]

    # Gather info
    access_mode = _detect_access_mode(args.file)
    dem_min, dem_max = _get_dem_range(args.dem)
    height_levels = _get_height_levels(args.file)

    print(f"\n{'=' * 70}")
    print("  Height-Subset End-to-End Benchmark")
    print(f"{'=' * 70}")
    print(f"  GUNW input:    {args.file}")
    print(f"  Access mode:   {access_mode}")
    print(f"  DEM file:      {args.dem}")
    if dem_min is not None:
        print(f"  DEM range:     {dem_min:.1f} to {dem_max:.1f} m")
    if height_levels is not None:
        print(f"  Height levels: {height_levels} ({len(height_levels)} bands)")

        from ARIAtools.util.interp import _get_height_subset_indices
        if dem_min is not None:
            idx = _get_height_subset_indices(
                np.array(height_levels), dem_min, dem_max, pad=0)
            print(f"  Subset levels: {[height_levels[i] for i in idx]} "
                  f"({len(idx)} bands)")
            print(f"  Band reduction: {len(height_levels)} → {len(idx)} "
                  f"({100*(1-len(idx)/len(height_levels)):.0f}% saving)")

    if args.bbox:
        print(f"  Bounding box:  {args.bbox}")
    print(f"  Layers:        {', '.join(layers)}")
    if args.compare_access and access_mode != 'local':
        print(f"  Access compare: S3 vs vsicurl (HTTPS)")
    print(f"{'=' * 70}")

    # Define run configurations
    # Each is (label, disable_subset, force_https)
    if args.compare_access and access_mode != 'local':
        configs = [
            ('S3 full cube',        True,  False),
            ('S3 + subset',         False, False),
            ('vsicurl full cube',   True,  True),
            ('vsicurl + subset',    False, True),
        ]
    else:
        configs = [
            ('full cube',           True,  False),
            ('with subset',         False, False),
        ]

    results = []

    for layer in layers:
        print(f"\n--- Benchmarking layer: {layer} ---")

        run_results = {}
        dirs = {}

        for label, disable_sub, force_https in configs:
            tag = label.replace(' ', '_')
            d = tempfile.mkdtemp(prefix=f'bench_{tag}_{layer}_')
            dirs[label] = d

            mode_str = 'vsicurl' if force_https else access_mode
            sub_str = 'OFF' if disable_sub else 'ON'
            print(f"  [{mode_str}, subset={sub_str}] Running...")

            t, rc, _, stderr = _run_extract(
                args.file, args.dem, args.bbox, layer, d,
                disable_subset=disable_sub,
                force_https=force_https)

            if rc != 0:
                print(f"    ERROR: ariaExtract failed. stderr:")
                print(stderr[:500])
                run_results[label] = {'time': t, 'error': True}
                continue

            print(f"    Time: {t:.2f}s")
            run_results[label] = {'time': t, 'error': False}

        # Compare outputs: all runs vs the first successful one
        ref_label = next(
            (l for l, cfg in zip(
                [c[0] for c in configs],
                configs) if not run_results.get(l, {}).get('error')),
            None)

        for label in run_results:
            r = run_results[label]
            if r.get('error'):
                continue
            if label == ref_label:
                r['max_diff'] = 0.0
            else:
                comp = _compare_outputs(
                    dirs[ref_label], dirs[label], layer)
                r['max_diff'] = comp['max_diff']
                r['n_files'] = comp['n_files']

        # Build result entries
        for label, disable_sub, force_https in configs:
            r = run_results.get(label, {})
            mode_str = 'vsicurl' if force_https else (
                's3' if access_mode == 'vsis3' else access_mode)
            results.append({
                'layer': layer,
                'mode': mode_str,
                'subset': not disable_sub,
                'label': f"{layer} [{mode_str}{'+ subset' if not disable_sub else ''}]",
                'time': r.get('time', 0),
                'max_diff': r.get('max_diff', float('nan')),
                'error': r.get('error', True),
            })

        # Cleanup
        if not args.keep_outputs:
            for d in dirs.values():
                shutil.rmtree(d, ignore_errors=True)
        else:
            for label, d in dirs.items():
                print(f"  {label}: {d}")

    # Print summary
    print(f"\n{'=' * 78}")
    print("  SUMMARY")
    print(f"{'=' * 78}")
    if dem_min is not None:
        print(f"  DEM range: {dem_min:.1f} to {dem_max:.1f} m")
    if height_levels is not None:
        print(f"  Heights:   {len(height_levels)} bands: {height_levels}")
    print()
    print(f"  {'Layer':<22} {'Access':<10} {'Subset':<8} "
          f"{'Time (s)':>10} {'MaxDiff':>12} {'Status':>8}")
    print(f"  {'-'*22} {'-'*10} {'-'*8} {'-'*10} {'-'*12} {'-'*8}")
    for r in results:
        sub_str = 'ON' if r['subset'] else 'OFF'
        if r.get('error'):
            print(f"  {r['layer']:<22} {r['mode']:<10} {sub_str:<8} "
                  f"{r['time']:10.2f} {'':>12} {'FAILED':>8}")
        else:
            print(f"  {r['layer']:<22} {r['mode']:<10} {sub_str:<8} "
                  f"{r['time']:10.2f} {r['max_diff']:12.2e} {'OK':>8}")

    # Compute and show speedups
    print()
    layer_set = list(dict.fromkeys(r['layer'] for r in results))
    for layer in layer_set:
        lr = [r for r in results if r['layer'] == layer and not r.get('error')]
        if len(lr) < 2:
            continue
        ref = lr[0]  # first (full cube, default access) is baseline
        print(f"  Speedups for {layer} (vs {ref['mode']} full cube "
              f"@ {ref['time']:.2f}s):")
        for r in lr[1:]:
            speedup = ref['time'] / r['time'] if r['time'] > 0 else float('inf')
            sub_str = '+ subset' if r['subset'] else 'full'
            print(f"    {r['mode']:>10} {sub_str:<10} → "
                  f"{speedup:.2f}x  ({r['time']:.2f}s)")

    print(f"\n{'=' * 78}")


if __name__ == '__main__':
    main()
