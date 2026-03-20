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

# S3 credential endpoints — each NASA DAAC Cumulus deployment has its
# own endpoint that exchanges Earthdata Login tokens for temporary
# AWS credentials scoped to that DAAC's S3 bucket(s).
_S3_CREDS_ENDPOINTS = {
    'default': 'https://cumulus.asf.alaska.edu/s3credentials',
    'nisar':   'https://nisar.asf.earthdatacloud.nasa.gov/s3credentials',
}

# Map S3 bucket names to credential endpoint keys
_BUCKET_TO_ENDPOINT = {
    'sds-n-cumulus-prod-nisar-products': 'nisar',
}

# Cached credentials (module-level) — keyed by endpoint name
_s3_creds_cache = {}       # endpoint_key → creds dict
_s3_creds_expiry_cache = {}  # endpoint_key → expiry timestamp

# Reverse mapping: /vsis3/ path → HTTPS URL, populated by maybe_use_s3.
# Used by open_gunw_h5 to resolve /vsis3/ paths back to HTTPS for h5py.
_vsis3_to_https = {}


def is_on_aws():
    """Detect whether the current environment is running on AWS.

    Detection strategy (in order):

    1. **IMDSv2** — EC2 Instance Metadata Service v2.  Works on
       standard EC2 instances.
    2. **boto3 STS** — ``GetCallerIdentity``.  Works in containers,
       JupyterHub on EKS/ECS, SageMaker, Lambda, and any environment
       where IAM credentials are available (instance role, task role,
       env vars, or config files).

    Returns ``False`` on any failure, so this is safe to call from
    any environment.

    Returns
    -------
    bool
    """
    # --- Attempt 1: IMDSv2 (fast, no extra dependency) ---
    try:
        token_resp = requests.put(
            'http://169.254.169.254/latest/api/token',
            headers={'X-aws-ec2-metadata-token-ttl-seconds': '60'},
            timeout=1.5)
        token_resp.raise_for_status()
        token = token_resp.text

        meta_resp = requests.get(
            'http://169.254.169.254/latest/meta-data/instance-id',
            headers={'X-aws-ec2-metadata-token': token},
            timeout=1.5)
        meta_resp.raise_for_status()
        LOGGER.info('Running on AWS EC2 (instance %s)',
                    meta_resp.text)
        return True
    except Exception:
        LOGGER.debug('IMDSv2 unavailable — trying boto3 fallback')

    # --- Attempt 2: boto3 STS GetCallerIdentity ---
    try:
        import boto3
        import botocore.exceptions
        sts = boto3.client('sts', region_name='us-west-2')
        identity = sts.get_caller_identity()
        LOGGER.info(
            'Running on AWS (STS identity: %s)', identity.get('Arn'))
        return True
    except ImportError:
        LOGGER.debug('boto3 not installed — cannot use STS fallback')
    except Exception:
        LOGGER.debug('boto3 STS call failed — not on AWS or no '
                     'credentials available')

    LOGGER.debug('Not running on AWS (all detection methods failed)')
    return False


def _endpoint_key_for_bucket(bucket):
    """Return the credential endpoint key for a given S3 bucket."""
    return _BUCKET_TO_ENDPOINT.get(bucket, 'default')


def _endpoint_key_for_s3uri(s3_uri):
    """Return the credential endpoint key for an ``s3://`` URI."""
    bucket = s3_uri[len('s3://'):].split('/', 1)[0]
    return _endpoint_key_for_bucket(bucket)


def _fetch_s3_credentials(endpoint_key='default'):
    """Fetch temporary S3 credentials from a DAAC endpoint.

    Requires Earthdata Login credentials in ``~/.netrc``.  The
    returned credentials are valid for ~1 hour.

    Parameters
    ----------
    endpoint_key : str
        Key into ``_S3_CREDS_ENDPOINTS`` (e.g. ``'default'``,
        ``'nisar'``).

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
    url = _S3_CREDS_ENDPOINTS[endpoint_key]
    LOGGER.info('Fetching temporary S3 credentials from %s', url)
    try:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        creds = resp.json()
    except Exception as exc:
        raise RuntimeError(
            f'Failed to fetch S3 credentials from {url}: '
            f'{exc}') from exc

    required = ('accessKeyId', 'secretAccessKey', 'sessionToken')
    for key in required:
        if key not in creds:
            raise RuntimeError(
                f'S3 credential response missing key {key!r}')

    LOGGER.debug('S3 credentials obtained from %s (expires: %s)',
                 endpoint_key, creds.get('expiration', 'unknown'))
    return creds


def get_s3_credentials(endpoint_key='default'):
    """Return cached S3 credentials, refreshing if expired.

    Parameters
    ----------
    endpoint_key : str
        Key into ``_S3_CREDS_ENDPOINTS`` (e.g. ``'default'``,
        ``'nisar'``).

    Returns
    -------
    dict
        Keys: ``accessKeyId``, ``secretAccessKey``, ``sessionToken``.
    """
    expiry = _s3_creds_expiry_cache.get(endpoint_key, 0)

    # Refresh 5 minutes before expiry to avoid mid-operation failures
    if endpoint_key not in _s3_creds_cache or \
            time.time() > (expiry - 300):
        _s3_creds_cache[endpoint_key] = \
            _fetch_s3_credentials(endpoint_key)
        # Default to 1 hour if no expiration provided
        _s3_creds_expiry_cache[endpoint_key] = time.time() + 3600
    return _s3_creds_cache[endpoint_key]


def configure_gdal_s3(endpoint_key='default'):
    """Set GDAL config options for S3 access.

    Should be called once before opening ``/vsis3/`` paths with GDAL.
    Automatically fetches/refreshes temporary credentials from the
    appropriate DAAC endpoint.

    Credentials are also exported as environment variables so that
    child processes (GNU parallel workers, Dask process workers)
    inherit them automatically.

    Parameters
    ----------
    endpoint_key : str
        Key into ``_S3_CREDS_ENDPOINTS`` (e.g. ``'default'``,
        ``'nisar'``).
    """
    creds = get_s3_credentials(endpoint_key)

    _set = osgeo.gdal.SetConfigOption
    _set('AWS_ACCESS_KEY_ID', creds['accessKeyId'])
    _set('AWS_SECRET_ACCESS_KEY', creds['secretAccessKey'])
    _set('AWS_SESSION_TOKEN', creds['sessionToken'])
    _set('AWS_REGION', 'us-west-2')
    _set('AWS_NO_SIGN_REQUEST', 'NO')

    # Export as env vars so child processes inherit the credentials
    os.environ['AWS_ACCESS_KEY_ID'] = creds['accessKeyId']
    os.environ['AWS_SECRET_ACCESS_KEY'] = creds['secretAccessKey']
    os.environ['AWS_SESSION_TOKEN'] = creds['sessionToken']
    os.environ['AWS_DEFAULT_REGION'] = 'us-west-2'

    LOGGER.info('GDAL configured for S3 direct access '
                '(endpoint: %s, region: us-west-2)', endpoint_key)


def restore_gdal_s3_from_env():
    """Restore GDAL S3 config from environment variables.

    Called by worker processes (GNU parallel, Dask) that inherit
    AWS credentials via environment variables but need the GDAL
    config options set in their own process.
    """
    key_id = os.environ.get('AWS_ACCESS_KEY_ID')
    if not key_id:
        return

    _set = osgeo.gdal.SetConfigOption
    _set('AWS_ACCESS_KEY_ID', key_id)
    _set('AWS_SECRET_ACCESS_KEY', os.environ['AWS_SECRET_ACCESS_KEY'])
    _set('AWS_SESSION_TOKEN', os.environ['AWS_SESSION_TOKEN'])
    _set('AWS_REGION', os.environ.get('AWS_DEFAULT_REGION', 'us-west-2'))
    _set('AWS_NO_SIGN_REQUEST', 'NO')
    LOGGER.debug('GDAL S3 config restored from environment variables')


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

    # Determine credential endpoint from the first available S3 URL
    first_s3 = next(u for u in s3_urls if u is not None)
    endpoint_key = _endpoint_key_for_s3uri(first_s3)

    # We're on AWS with S3 URLs — set up credentials and convert
    configure_gdal_s3(endpoint_key)

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


def fixup_vrt_s3_paths(directory):
    """Replace ``/vsis3/`` paths with ``/vsicurl/`` paths in all VRT files.

    VRT files are persistent output artifacts consumed by downstream
    tools (e.g. MintPy) that don't have ARIA-tools' S3 credential
    setup.  This function rewrites ``/vsis3/bucket/key`` references
    to ``/vsicurl/https_url`` using the reverse mapping built by
    ``maybe_use_s3()``, making the VRTs portable.

    Should be called after all extraction/export is complete and
    before downstream tools read the VRTs.

    Parameters
    ----------
    directory : str
        Root directory to scan for ``.vrt`` files (recursive).
    """
    if not _vsis3_to_https:
        return  # No S3 paths were used — nothing to fix

    import glob as _glob
    vrt_files = _glob.glob(os.path.join(directory, '**', '*.vrt'),
                           recursive=True)
    n_fixed = 0
    for vrt_path in vrt_files:
        try:
            with open(vrt_path, 'r') as f:
                content = f.read()
        except (OSError, UnicodeDecodeError):
            continue

        if '/vsis3/' not in content:
            continue

        new_content = content
        for vsis3_path, https_url in _vsis3_to_https.items():
            new_content = new_content.replace(
                vsis3_path, f'/vsicurl/{https_url}')

        if new_content != content:
            with open(vrt_path, 'w') as f:
                f.write(new_content)
            n_fixed += 1

    if n_fixed:
        LOGGER.info('Fixed %d VRT file(s): replaced /vsis3/ with '
                     '/vsicurl/ for downstream compatibility', n_fixed)


def get_s3_client(endpoint_key='default', max_pool_connections=10):
    """Create a boto3 S3 client using temporary DAAC credentials.

    Parameters
    ----------
    endpoint_key : str
        Key into ``_S3_CREDS_ENDPOINTS`` (e.g. ``'default'``,
        ``'nisar'``).
    max_pool_connections : int
        Maximum number of connections in the urllib3 pool.
        Set this to match the number of concurrent download
        threads to avoid connection pool overflow warnings.

    Returns
    -------
    boto3.client
    """
    import boto3
    from botocore.config import Config
    creds = get_s3_credentials(endpoint_key)
    config = Config(
        max_pool_connections=max_pool_connections,
    )
    return boto3.client(
        's3',
        aws_access_key_id=creds['accessKeyId'],
        aws_secret_access_key=creds['secretAccessKey'],
        aws_session_token=creds['sessionToken'],
        region_name='us-west-2',
        config=config)


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
