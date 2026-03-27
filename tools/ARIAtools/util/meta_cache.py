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

Also provides h5py-based access for scalar/string HDF5 metadata
that GDAL cannot read (e.g. boundingPolygon, centerFrequency).
"""

import hashlib
import http.cookiejar
import io
import json
import logging
import os
import time

import h5py
import numpy as np
import osgeo.gdal
import requests

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

    For remote URLs we strip the /vsicurl/ or /vsis3/ prefix so the
    key is the canonical URL.  For local files we use the absolute path.
    """
    import re
    key = re.sub(r'^/vsi(curl|s3)/', '', fname)
    # Strip the bucket prefix for S3 paths to recover the original URL
    # e.g. asf-cumulus-prod-nisar-gunw/path → path
    # The key should be the URL path or local path, not bucket-specific
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
    if fname.endswith('.nc'):
        subdatasets_raw = metadata.get('SUBDATASETS', {})
        subdatasets = [
            v for k, v in sorted(subdatasets_raw.items()) if 'NAME' in k]
    if fname.endswith('.h5'):
        subdatasets_raw = metadata.get('Subdatasets', {})
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


# ---------- h5py helpers for scalar/string HDF5 metadata ------------------


class _HTTPRangeFile(io.RawIOBase):
    """Seekable read-only file backed by HTTP byte-range requests.

    Handles the Earthdata OAuth redirect chain by resolving the final
    signed URL up front, then issues byte-range GET requests for each
    read.  h5py uses this to fetch only the HDF5 chunks it needs.
    """

    def __init__(self, url, session=None):
        super().__init__()
        self._session = session or requests.Session()
        self._pos = 0

        # Resolve the redirect chain to get the final signed URL +
        # determine file size.
        resp = self._session.get(
            url, headers={'Range': 'bytes=0-0'},
            allow_redirects=True, stream=True, timeout=60)
        resp.close()

        # Store the final (signed) URL so subsequent requests skip
        # the redirect chain entirely.
        self._url = resp.url

        # Parse total file size from Content-Range: bytes 0-0/<total>
        cr = resp.headers.get('Content-Range', '')
        if '/' in cr:
            self._size = int(cr.split('/')[-1])
        else:
            head = self._session.head(self._url, timeout=60)
            self._size = int(head.headers.get('Content-Length', 0))

        LOGGER.debug('Resolved URL (size=%d bytes)', self._size)

    def readable(self):
        return True

    def writable(self):
        return False

    def seekable(self):
        return True

    def tell(self):
        return self._pos

    def seek(self, offset, whence=io.SEEK_SET):
        if whence == io.SEEK_SET:
            self._pos = offset
        elif whence == io.SEEK_CUR:
            self._pos += offset
        elif whence == io.SEEK_END:
            self._pos = self._size + offset
        return self._pos

    def read(self, size=-1):
        if self._pos >= self._size:
            return b''
        if size < 0:
            size = self._size - self._pos
        end = min(self._pos + size - 1, self._size - 1)
        headers = {'Range': f'bytes={self._pos}-{end}'}
        resp = self._session.get(
            self._url, headers=headers, timeout=120)
        resp.raise_for_status()
        data = resp.content
        self._pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        n = len(data)
        b[:n] = data
        return n

    @property
    def size(self):
        return self._size


def _get_earthdata_session():
    """Create a ``requests.Session`` with Earthdata cookie auth."""
    cookie_file = osgeo.gdal.GetConfigOption(
        'GDAL_HTTP_COOKIEFILE', '/tmp/cookies.txt')

    session = requests.Session()

    if os.path.isfile(cookie_file):
        jar = http.cookiejar.MozillaCookieJar(cookie_file)
        try:
            jar.load(ignore_discard=True, ignore_expires=True)
            session.cookies.update(jar)
            LOGGER.debug('Loaded cookies from %s (%d cookies)',
                         cookie_file, len(jar))
        except Exception as exc:
            LOGGER.debug('Could not load cookies from %s: %s',
                         cookie_file, exc)

    return session


def open_gunw_h5(url_or_path):
    """Open a GUNW HDF5 file (local or remote) as an h5py.File.

    For remote URLs (https://...), resolves the Earthdata OAuth
    redirect chain, then uses HTTP byte-range requests so that h5py
    fetches only the HDF5 chunks it needs.

    Parameters
    ----------
    url_or_path : str
        HTTPS URL or local file path.  VSICURL/VSIS3 prefixes
        are stripped automatically.  For ``/vsis3/`` paths, the
        underlying HTTPS URL is resolved via the reverse mapping
        built by ``ARIAtools.util.s3.maybe_use_s3()``.

    Returns
    -------
    h5py.File
    """
    path = url_or_path.replace('/vsicurl/', '')

    # For /vsis3/ paths, resolve back to HTTPS for h5py access
    if path.startswith('/vsis3/'):
        from ARIAtools.util.s3 import vsis3_to_https
        https_url = vsis3_to_https(path)
        if https_url:
            path = https_url
        else:
            LOGGER.warning(
                'Cannot resolve /vsis3/ path to HTTPS: %s', path)

    if path.startswith('https://') or path.startswith('http://'):
        LOGGER.debug('Opening remote HDF5: %s', path)
        session = _get_earthdata_session()
        fh = _HTTPRangeFile(path, session=session)
        return h5py.File(fh, 'r')
    else:
        LOGGER.debug('Opening local HDF5: %s', path)
        return h5py.File(path, 'r')


# ---------- h5py scalar/string metadata caching --------------------------


def _extract_h5_fields(fname, h5_fields):
    """Read scalar/string metadata from an HDF5 file via h5py.

    Parameters
    ----------
    fname : str
        Product path or ``/vsicurl/...`` URL.
    h5_fields : dict
        Mapping of ``{field_name: hdf5_path}``.  Each path is read
        from the HDF5 file and returned as a Python scalar or string.

    Returns
    -------
    dict
        Keys matching *h5_fields* with their scalar values.
    """
    result = {}
    with open_gunw_h5(fname) as h5f:
        for field_name, h5_path in h5_fields.items():
            if h5_path in h5f:
                val = h5f[h5_path][()]
                # Decode bytes → str for string datasets
                if isinstance(val, bytes):
                    val = val.decode('utf-8')
                elif isinstance(val, np.generic):
                    val = val.item()
                result[field_name] = val
            else:
                LOGGER.warning('h5py field %s not found at %s',
                               field_name, h5_path)
    return result


def get_h5_field(fname, field_name, h5_fields, cache_data):
    """Return a single h5py metadata field, using cache.

    On the first call for a given product, all fields in *h5_fields*
    are read in one ``h5py.File`` open and cached together.
    Subsequent calls return from cache without any I/O.

    Parameters
    ----------
    fname : str
        Product path (may include ``/vsicurl/`` prefix).
    field_name : str
        One of the keys in *h5_fields*.
    h5_fields : dict
        Mapping of ``{field_name: hdf5_path}`` — the caller defines
        which datasets to read.
    cache_data : dict
        Mutable cache dict; updated in-place on cache miss.

    Returns
    -------
    value
        The scalar/string value for the requested field.
    """

    key = _file_key(fname)
    entry = cache_data.get(key, {})

    # handle subdatasets behavior differently
    # full filename path is expected in this workflow
    if field_name == 'subdatasets':
        
        # Check if this specific field is already cached
        if field_name in entry:
            LOGGER.debug(
                'Cache hit: %s [%s]', os.path.basename(key), field_name
            )
            return entry

        # Cache miss - extract using GDAL Info
        LOGGER.debug(
            'Cache miss - reading: %s [%s]', os.path.basename(key), field_name
        )
        
        # GDAL Info requires NETCDF prefix to properly read HDF5 subdatasets
        gdal_fname = f'NETCDF:"{fname}'
        meta = osgeo.gdal.Info(gdal_fname)
        
        # Filter the requested fields against the GDAL metadata
        sdskeys_addlyrs = [k for k in h5_fields if k in meta]
        
        # Merge into the existing cache entry safely
        if key not in cache_data:
            cache_data[key] = {}
        cache_data[key][field_name] = sdskeys_addlyrs

        return cache_data[key]

    # Check if this specific h5py field is already cached
    h5_cache_key = f'h5_{field_name}'
    if h5_cache_key in entry:
        LOGGER.debug('Cache hit: %s [%s]',
                     os.path.basename(key), field_name)
        return entry[h5_cache_key]

    # Cache miss — extract all requested h5py fields at once
    LOGGER.debug('Cache miss – reading: %s',
                 os.path.basename(key))
    h5_meta = _extract_h5_fields(fname, h5_fields)

    # Merge into the existing cache entry
    if key not in cache_data:
        cache_data[key] = {}
    for k, v in h5_meta.items():
        cache_data[key][f'h5_{k}'] = v

    if field_name not in h5_meta:
        raise RuntimeError(
            f'h5py field {field_name!r} not found in {fname}')
    return h5_meta[field_name]
