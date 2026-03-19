# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert
# Copyright (c) 2026
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
AWS S3 direct-access helpers for ARIA GUNW virtual processing.

When running on an AWS EC2 instance (any region), GDAL can access
NASA Earthdata products via ``/vsis3/`` instead of ``/vsicurl/``.
This bypasses the Earthdata OAuth redirect chain and uses the AWS
internal network, which is significantly faster.

The URL list file written by ``ariaDownload.py`` contains two
comma-separated columns per line::

    https://host/path/file.h5,s3://bucket/path/file.h5

On startup, ARIA-tools checks if it is running on AWS via the EC2
metadata service.  If yes, the S3 column is used (with ``/vsis3/``);
otherwise the HTTPS column is used (with ``/vsicurl/``).

The module provides:
- AWS instance detection via the EC2 metadata service (IMDSv2)
- Temporary S3 credential retrieval from the ASF DAAC endpoint
- GDAL configuration for S3 access
- Parsing of the 2-column URL file format
"""

import logging
import os
import time

import osgeo.gdal
import requests

LOGGER = logging.getLogger(__name__)

# ASF DAAC S3 credential endpoint — returns temporary AWS credentials
# for in-cloud access to NASA Earthdata products hosted on S3.
_ASF_S3_CREDS_URL = 'https://cumulus.asf.alaska.edu/s3credentials'

# Cached credentials (module-level) — refreshed when expired
_s3_creds = None
_s3_creds_expiry = 0

# Reverse mapping: /vsis3/ path → HTTPS URL, populated by maybe_use_s3.
# Used by open_gunw_h5 to resolve /vsis3/ paths back to HTTPS for h5py.
_vsis3_to_https = {}


def is_on_aws():
    """Detect whether the current environment is an AWS EC2 instance.

    Uses IMDSv2 (Instance Metadata Service v2) with a 1.5-second
    timeout.  Returns ``False`` on any failure, so this is safe to
    call from any environment.

    Returns
    -------
    bool
    """
    try:
        # IMDSv2 requires a PUT to get a session token first
        token_resp = requests.put(
            'http://169.254.169.254/latest/api/token',
            headers={'X-aws-ec2-metadata-token-ttl-seconds': '60'},
            timeout=1.5)
        token_resp.raise_for_status()
        token = token_resp.text

        # Verify we can read instance metadata
        meta_resp = requests.get(
            'http://169.254.169.254/latest/meta-data/instance-id',
            headers={'X-aws-ec2-metadata-token': token},
            timeout=1.5)
        meta_resp.raise_for_status()
        LOGGER.info('Running on AWS EC2 (instance %s)',
                    meta_resp.text)
        return True
    except Exception:
        LOGGER.debug('Not running on AWS EC2 (metadata service '
                     'unreachable)')
        return False


def _fetch_s3_credentials():
    """Fetch temporary S3 credentials from the ASF DAAC endpoint.

    Requires Earthdata Login credentials in ``~/.netrc``.  The
    returned credentials are valid for ~1 hour.

    Returns
    -------
    dict
        Keys: ``accessKeyId``, ``secretAccessKey``, ``sessionToken``,
        ``expiration``.

    Raises
    ------
    RuntimeError
        If credential retrieval fails.
    """
    LOGGER.info('Fetching temporary S3 credentials from ASF')
    try:
        resp = requests.get(_ASF_S3_CREDS_URL, timeout=30)
        resp.raise_for_status()
        creds = resp.json()
    except Exception as exc:
        raise RuntimeError(
            f'Failed to fetch S3 credentials from {_ASF_S3_CREDS_URL}: '
            f'{exc}') from exc

    required = ('accessKeyId', 'secretAccessKey', 'sessionToken')
    for key in required:
        if key not in creds:
            raise RuntimeError(
                f'S3 credential response missing key {key!r}')

    LOGGER.debug('S3 credentials obtained (expires: %s)',
                 creds.get('expiration', 'unknown'))
    return creds


def get_s3_credentials():
    """Return cached S3 credentials, refreshing if expired.

    Returns
    -------
    dict
        Keys: ``accessKeyId``, ``secretAccessKey``, ``sessionToken``.
    """
    global _s3_creds, _s3_creds_expiry

    # Refresh 5 minutes before expiry to avoid mid-operation failures
    if _s3_creds is None or time.time() > (_s3_creds_expiry - 300):
        _s3_creds = _fetch_s3_credentials()
        # Default to 1 hour if no expiration provided
        _s3_creds_expiry = time.time() + 3600
    return _s3_creds


def configure_gdal_s3():
    """Set GDAL config options for S3 access using ASF credentials.

    Should be called once before opening ``/vsis3/`` paths with GDAL.
    Automatically fetches/refreshes temporary credentials.
    """
    creds = get_s3_credentials()

    _set = osgeo.gdal.SetConfigOption
    _set('AWS_ACCESS_KEY_ID', creds['accessKeyId'])
    _set('AWS_SECRET_ACCESS_KEY', creds['secretAccessKey'])
    _set('AWS_SESSION_TOKEN', creds['sessionToken'])
    _set('AWS_REGION', 'us-west-2')
    _set('AWS_NO_SIGN_REQUEST', 'NO')

    LOGGER.info('GDAL configured for S3 direct access (region: us-west-2)')


def s3uri_to_vsis3(s3_uri):
    """Convert an ``s3://bucket/path`` URI to a GDAL ``/vsis3/`` path.

    Parameters
    ----------
    s3_uri : str
        S3 URI, e.g. ``s3://bucket-name/path/to/file.h5``

    Returns
    -------
    str
        ``/vsis3/bucket-name/path/to/file.h5``
    """
    return '/vsis3/' + s3_uri[len('s3://'):]


def parse_url_line(line):
    """Parse a line from the URL list file.

    Supports both the 2-column format (``https_url,s3_url``) and the
    legacy single-column format (``https_url`` only).

    Parameters
    ----------
    line : str
        A single line from the URL file (stripped of newline).

    Returns
    -------
    tuple of (str, str or None)
        ``(https_url, s3_url)`` where ``s3_url`` is ``None`` if not
        present or empty.
    """
    parts = line.split(',', 1)
    https_url = parts[0].strip()
    s3_url = parts[1].strip() if len(parts) > 1 and parts[1].strip() else None
    return https_url, s3_url


def maybe_use_s3(https_urls, s3_urls):
    """If on AWS and S3 URLs are available, use ``/vsis3/`` access.

    This is the main entry point.  Call it once with the URL lists
    parsed from the URL file before processing begins.

    Parameters
    ----------
    https_urls : list of str
        Product HTTPS URLs.
    s3_urls : list of (str or None)
        Corresponding S3 URIs (``s3://...``) from the URL file.
        ``None`` entries indicate no S3 URL was available.

    Returns
    -------
    list of str
        GDAL-ready paths: ``/vsis3/...`` if on AWS with S3 URLs
        available, otherwise the original HTTPS URLs (caller wraps
        with ``/vsicurl/``).
    bool
        ``True`` if S3 access is being used.
    """
    # Check if any S3 URLs are available
    has_s3 = any(u is not None for u in s3_urls)
    if not has_s3:
        LOGGER.info('No S3 URLs in URL file — using HTTPS access')
        return https_urls, False

    if not is_on_aws():
        LOGGER.info('Not on AWS — using HTTPS (vsicurl) access')
        return https_urls, False

    # We're on AWS with S3 URLs — set up credentials and convert
    configure_gdal_s3()

    converted = []
    for https_url, s3_url in zip(https_urls, s3_urls):
        if s3_url is not None:
            vsis3_path = s3uri_to_vsis3(s3_url)
            converted.append(vsis3_path)
            _vsis3_to_https[vsis3_path] = https_url
        else:
            # Fallback for products without S3 URLs
            converted.append(https_url)
    LOGGER.info('Using S3 direct access for %d/%d products',
                sum(1 for u in s3_urls if u is not None),
                len(https_urls))
    return converted, True


def vsis3_to_https(vsis3_path):
    """Look up the HTTPS URL for a ``/vsis3/`` path.

    Uses the reverse mapping built by ``maybe_use_s3()``.

    Returns
    -------
    str or None
    """
    return _vsis3_to_https.get(vsis3_path)


def get_s3_client():
    """Create a boto3 S3 client using ASF temporary credentials.

    Credentials are automatically fetched/refreshed.

    Returns
    -------
    boto3.client
    """
    import boto3
    creds = get_s3_credentials()
    return boto3.client(
        's3',
        aws_access_key_id=creds['accessKeyId'],
        aws_secret_access_key=creds['secretAccessKey'],
        aws_session_token=creds['sessionToken'],
        region_name='us-west-2')


def parse_s3_uri(s3_uri):
    """Split an ``s3://bucket/key`` URI into bucket and key.

    Parameters
    ----------
    s3_uri : str
        e.g. ``s3://sds-n-cumulus-prod-nisar-products/path/file.h5``

    Returns
    -------
    tuple of (str, str)
        ``(bucket, key)``
    """
    without_scheme = s3_uri[len('s3://'):]
    bucket, key = without_scheme.split('/', 1)
    return bucket, key
