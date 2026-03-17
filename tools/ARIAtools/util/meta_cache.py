# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright (c) 2026, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Metadata cache for ARIA GUNW products.

Caches per-product metadata extracted from remote (vsicurl) or local
NetCDF files so that repeated invocations avoid re-reading product
headers over HTTP. The cache is stored as a JSON sidecar file.
"""

import hashlib
import json
import logging
import os
import time

import osgeo.gdal

LOGGER = logging.getLogger(__name__)

# Sentinel used when a key is genuinely absent from a product
_MISSING = '__missing__'


def _cache_path(url_file):
    """Derive the cache JSON path from the URL list file directory.

    Uses a fixed name in the same directory as the URL file so that
    different URL files (e.g. when expanding a time series) share a
    single cache of per-product metadata.
    """
    if url_file:
        cache_dir = os.path.dirname(os.path.abspath(url_file))
        return os.path.join(cache_dir, 'aria_meta_cache.json')
    return None


def _file_key(fname):
    """Produce a stable key for a product file (URL or local path).

    For remote URLs we strip the /vsicurl/ prefix so the key is the
    canonical URL.  For local files we use the absolute path.
    """
    key = fname.replace('/vsicurl/', '')
    return key


def load_cache(cache_file):
    """Load existing cache from disk. Returns dict keyed by file URL/path."""
    if cache_file and os.path.isfile(cache_file):
        try:
            with open(cache_file) as fh:
                data = json.load(fh)
            LOGGER.debug('Loaded metadata cache with %d entries from %s',
                         len(data), cache_file)
            return data
        except (json.JSONDecodeError, OSError) as exc:
            LOGGER.warning('Could not load metadata cache %s: %s',
                           cache_file, exc)
    return {}


def save_cache(cache_file, cache_data):
    """Persist the cache dict to disk."""
    if cache_file is None:
        return
    try:
        with open(cache_file, 'w') as fh:
            json.dump(cache_data, fh, indent=1)
        LOGGER.debug('Saved metadata cache (%d entries) to %s',
                     len(cache_data), cache_file)
    except OSError as exc:
        LOGGER.warning('Could not write metadata cache %s: %s',
                       cache_file, exc)


def extract_metadata_gdal(fname):
    """Extract lightweight metadata from a GUNW product using only GDAL.

    Uses a single ``gdal.Info(..., options=['-json'])`` call which
    fetches the minimum number of HTTP range requests for remote files.

    Parameters
    ----------
    fname : str
        Path or ``/vsicurl/...`` path to the NetCDF product.

    Returns
    -------
    dict
        Dictionary with keys: version, driver, subdatasets, nc_global,
        gdal_info_ts (timestamp of extraction).
    """
    netcdf_fname = f'NETCDF:"{fname}'

    # gdal.Info with -json returns everything we need in one call
    info_str = osgeo.gdal.Info(netcdf_fname, options=['-json'])
    if info_str is None:
        LOGGER.warning('gdal.Info returned None for %s', fname)
        return None

    info = json.loads(info_str) if isinstance(info_str, str) else info_str

    # Extract version from NC_GLOBAL metadata
    metadata = info.get('metadata', {})
    nc_global = metadata.get('', {})
    version = nc_global.get('NC_GLOBAL#version', None)

    # Collect subdataset names
    subdatasets_raw = metadata.get('SUBDATASETS', {})
    subdatasets = [
        v for k, v in sorted(subdatasets_raw.items()) if 'NAME' in k]

    # Identify available troposphere models from subdataset paths
    tropo_models = set()
    for sd in subdatasets:
        parts = sd.split('/')
        if 'troposphere' in parts:
            idx = parts.index('troposphere')
            if idx + 1 < len(parts):
                tropo_models.add(parts[idx + 1])

    return {
        'version': version,
        'driver': info.get('driverShortName', None),
        'subdatasets': subdatasets,
        'tropo_models': sorted(tropo_models),
        'nc_global': nc_global,
        'gdal_info_ts': time.time(),
    }


def get_or_extract(fname, cache_data):
    """Return cached metadata for *fname*, extracting if not cached.

    Parameters
    ----------
    fname : str
        Product path (may include ``/vsicurl/`` prefix).
    cache_data : dict
        Mutable cache dict; will be updated in-place on cache miss.

    Returns
    -------
    dict or None
        Metadata dict for the product, or None on failure.
    """
    key = _file_key(fname)
    if key in cache_data:
        LOGGER.debug('Cache hit: %s', os.path.basename(key))
        return cache_data[key]

    LOGGER.debug('Cache miss – extracting metadata: %s',
                os.path.basename(key))
    meta = extract_metadata_gdal(fname)
    if meta is not None:
        cache_data[key] = meta
    return meta
