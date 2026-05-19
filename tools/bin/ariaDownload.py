#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Brett A. Buzzanga, David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import argparse
import concurrent.futures
import datetime
import getpass
import hashlib
import json
import logging
import math
import os
import re
import time

import ARIAtools.util.log
import ARIAtools.util.s3
import asf_search
import osgeo.gdal
import shapely
import tqdm
from ARIAtools.util.shp import open_shp
from ARIAtools.util.url import url_versions
from requests.exceptions import RequestException

LOGGER = logging.getLogger("ariaDownload.py")
VALIDATION_METADATA_SUFFIX = ".aria-download.json"


def createParser():
    """
    Download ARIA products using asf_search

    see: https://github.com/asfadmin/Discovery-asf_search
    """
    parser = argparse.ArgumentParser(
        description="Command line interface to download Sentinel-1/NISAR "
        "GUNW products from the ASF DAAC. \nDownloading them "
        "requires a NASA Earthdata URS user login",
        epilog="Examples of use:\n\n"
        "\t # Count Sentinel-1 products available for track 004\n"
        "\t ariaDownload.py --track 004 --output count\n\n"
        "\t # Download Sentinel-1 products within specified "
        "bounding box\n"
        '\t ariaDownload.py --bbox "36.75 37.225 -76.655 '
        '-75.928"\n\n'
        "\t # Count Sentinel-1 products for tracks 004 & 077 "
        "since Jan 2019\n"
        "\t ariaDownload.py --mission S1 -t 004,077 "
        "--start 20190101 -o count\n\n"
        "\t # Count all available NISAR products\n"
        "\t ariaDownload.py --mission NISAR -o count\n\n"
        "\t # Download globally available descending NISAR "
        "products\n"
        "\t ariaDownload.py --mission NISAR -d d "
        '-b "-90 90 -180 180"\n\n'
        "\t # Count all available NISAR products for track 172\n"
        "\t ariaDownload.py --mission NISAR -t 172 -o count\n\n"
        "\t # Download specific NISAR interferogram "
        "and query the globe for it\n"
        "\t ariaDownload.py --mission NISAR "
        '-b "-90 90 -180 180" -i 20251122_20251204\n',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-o",
        "--output",
        default="Download",
        type=str.title,
        choices=("Download", "Count", "Url"),
        help='Output type. Default="Download". Use "Url" for ingestion to ' "aria*.py",
    )
    parser.add_argument(
        "-t",
        "--track",
        default=None,
        type=str,
        help="track to download; single number or " "comma separated",
    )
    parser.add_argument(
        "-b",
        "--bbox",
        default="-90 90 -180 180",
        type=str,
        help='Lat/Lon Bounding SNWE (e.g., "36.75 37.225 -76.655 -75.928"), '
        'WKT POLYGON string (e.g., "POLYGON((lon lat, ...))"), '
        "or GDAL-readable file containing POLYGON geometry. "
        "Default is set to global scale",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        dest="wd",
        default="./products",
        type=str,
        help='Specify directory to deposit all outputs. Default is "products" '
        "in local directory where script is launched.",
    )
    parser.add_argument(
        "-s",
        "--start",
        default="20100101",
        type=str,
        help="Start date as YYYYMMDD; If none provided, starts at beginning "
        "of 2010.",
    )
    parser.add_argument(
        "-e",
        "--end",
        default="21000101",
        type=str,
        help="End date as YYYYMMDD. If none provided, ends today.",
    )
    parser.add_argument(
        "-u", "--user", default=None, type=str, help="NASA Earthdata URS user login."
    )
    parser.add_argument(
        "-p",
        "--pass",
        dest="passw",
        default=None,
        type=str,
        help="NASA Earthdata URS user password.",
    )
    parser.add_argument(
        "--mission",
        default="S1",
        type=str.upper,
        choices=("S1", "NISAR"),
        help="Sentinel-1 (S1) or NISAR. Default is S1",
    )
    parser.add_argument(
        "-l",
        "--daysless",
        dest="dayslt",
        default=math.inf,
        type=int,
        help="Take pairs with a temporal baseline -- days less than this " "value.",
    )
    parser.add_argument(
        "-m",
        "--daysmore",
        dest="daysgt",
        default=0,
        type=int,
        help="Take pairs with a temporal baseline -- days greater than this "
        "value. Example, annual pairs: ariaDownload.py -t 004 "
        "--daysmore 364.",
    )
    parser.add_argument(
        "-nt",
        "--num_threads",
        default="1",
        type=str,
        help="Specify number of threads for multiprocessing download. By "
        'default "1". Can also specify "All" to use all available '
        "threads.",
    )
    parser.add_argument(
        "-i",
        "--ifg",
        default=None,
        type=str,
        help="Retrieve one interferogram by its start/end date, specified as "
        "YYYYMMDD_YYYYMMDD (order independent).",
    )
    parser.add_argument(
        "-d",
        "--direction",
        dest="flightdir",
        default=None,
        type=str,
        help="Flight direction, options: ascending, a, descending, d",
    )
    parser.add_argument(
        "--version",
        default=None,
        help="Specify version as str, e.g. 2_0_4 or all prods. All products "
        "are downloaded by default. If version is specified, only "
        "products which match that version are downloaded. "
        "Not supported for NISAR currently.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print products to be downloaded to stdout",
    )
    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning", "error"],
        default="info",
        help="Logger log level. Default: info.",
    )
    return parser


def make_bbox(inp_bbox):
    """Make a WKT from SNWE, WKT string, or a shapefile"""
    if inp_bbox is None:
        return None

    if os.path.exists(os.path.abspath(inp_bbox)):
        ring = open_shp(inp_bbox, 0, 0).exterior
        poly = shapely.geometry.Polygon(ring)

    else:
        # Try WKT string first (e.g., "POLYGON((...))")
        if inp_bbox.strip().upper().startswith("POLYGON"):
            try:
                poly = shapely.wkt.loads(inp_bbox)
                return poly
            except Exception:
                pass  # Fall through to SNWE parsing

        # Parse as SNWE string
        try:
            S, N, W, E = (float(i) for i in inp_bbox.split())

            # adjust for degrees easting / northing (0 - 360 / 0:180)
            if W > 180:
                W -= 360
                LOGGER.info("AdjustedW")

            if E > 180:
                E -= 360
                LOGGER.info("AdjustedE")

            if N > 90:
                N -= 90
                S -= 90
                LOGGER.info("Adjusted N/S")

            # set poly object
            poly = shapely.geometry.Polygon([(W, N), (W, S), (E, S), (E, N)])

        except BaseException:
            raise Exception(
                "Cannot understand the --bbox argument. Input string was "
                "entered incorrectly or path does not exist."
            )

    return poly


def _get_s3_data_url(scene):
    """Extract the S3 data URL for the main product file from a scene.

    ASF search results include ``s3Urls`` with multiple files (browse,
    metadata, QA, etc.).  This returns the S3 URL matching the
    product's primary data file (same filename as the HTTPS URL).

    Returns
    -------
    str or None
        ``s3://bucket/path/file`` or ``None`` if not available.
    """
    props = scene.geojson()["properties"]
    filename = props.get("fileName", "")
    s3_urls = props.get("s3Urls", [])
    for s3_url in s3_urls:
        if s3_url.endswith(filename):
            return s3_url
    return None


def _get_scene_checksum(scene):
    """
    Extract MD5 checksum from ASF scene metadata.

    Parameters
    ----------
    scene : asf_search.ASFProduct
        Scene object from ASF search

    Returns
    -------
    tuple
        (checksum_value, checksum_type) or (None, None) if unavailable
    """
    props = scene.geojson()["properties"]

    # ASF provides MD5 checksums
    md5sum = props.get("md5sum") or props.get("md5")
    if md5sum:
        return md5sum, "md5"

    # Fallback to other checksum types if available
    sha256sum = props.get("sha256sum") or props.get("sha256")
    if sha256sum:
        return sha256sum, "sha256"

    return None, None


def _compute_file_checksum(filepath, checksum_type="md5", block_size=8192):
    """
    Compute checksum of a file.

    Parameters
    ----------
    filepath : str
        Path to file
    checksum_type : str
        Type of checksum: 'md5' or 'sha256'
    block_size : int
        Size of blocks to read (default 8KB)

    Returns
    -------
    str
        Hex digest of checksum, or None if error
    """
    if checksum_type.lower() == "md5":
        hasher = hashlib.md5()
    elif checksum_type.lower() == "sha256":
        hasher = hashlib.sha256()
    else:
        LOGGER.warning("Unsupported checksum type: %s", checksum_type)
        return None

    try:
        with open(filepath, "rb") as f:
            while True:
                data = f.read(block_size)
                if not data:
                    break
                hasher.update(data)
        return hasher.hexdigest()
    except Exception as e:
        LOGGER.warning("Error computing checksum for %s: %s", filepath, e)
        return None


def _validation_metadata_path(filepath):
    """Return the sidecar path used to cache trusted local validation metadata."""
    return f"{filepath}{VALIDATION_METADATA_SUFFIX}"


def _load_validation_metadata(filepath):
    """Load cached local validation metadata for a file if present."""
    meta_path = _validation_metadata_path(filepath)
    if not os.path.exists(meta_path):
        return None

    try:
        with open(meta_path) as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        LOGGER.debug("Could not read validation metadata %s: %s", meta_path, exc)
        return None


def _persist_validation_metadata(filepath, checksum_type=None, checksum=None):
    """Persist local validation metadata after a successful download or check."""
    checksum_type = checksum_type or "md5"
    checksum = checksum or _compute_file_checksum(filepath, checksum_type)
    if checksum is None:
        LOGGER.debug("Skipping validation metadata write for %s", filepath)
        return

    payload = {
        "size": os.path.getsize(filepath),
        "checksum_type": checksum_type,
        "checksum": checksum,
    }
    meta_path = _validation_metadata_path(filepath)
    try:
        with open(meta_path, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
    except OSError as exc:
        LOGGER.debug("Could not write validation metadata %s: %s", meta_path, exc)


def _remote_file_info_from_response(response):
    """Extract total size and range support from an HTTP response."""
    expected_size = 0
    supports_range = False

    content_range = response.headers.get("Content-Range", "")
    if "/" in content_range:
        try:
            expected_size = int(content_range.rsplit("/", 1)[-1])
            supports_range = True
        except ValueError:
            expected_size = 0

    if expected_size == 0:
        try:
            expected_size = int(response.headers.get("Content-Length", 0))
        except (TypeError, ValueError):
            expected_size = 0

    accept_ranges = response.headers.get("Accept-Ranges", "").lower()
    if accept_ranges == "bytes":
        supports_range = True

    return expected_size, supports_range


def _probe_remote_file_info(url, session):
    """Return trusted remote file size metadata when available."""
    try:
        head_response = session.head(url, allow_redirects=True, timeout=10)
        head_response.raise_for_status()
        expected_size, supports_range = _remote_file_info_from_response(head_response)
        if expected_size > 0:
            return expected_size, supports_range
    except Exception as exc:
        LOGGER.debug("HEAD metadata probe failed for %s: %s", url, exc)

    try:
        probe_response = session.get(
            url,
            headers={"Range": "bytes=0-0"},
            allow_redirects=True,
            stream=True,
            timeout=10,
        )
        probe_response.raise_for_status()
        expected_size, supports_range = _remote_file_info_from_response(probe_response)
        probe_response.close()
        if expected_size > 0:
            return expected_size, supports_range
    except Exception as exc:
        LOGGER.debug("Range metadata probe failed for %s: %s", url, exc)

    return 0, False


class DownloadProgressBar:
    """
    Progress bar for individual file downloads with resume support.

    Wraps tqdm to show bytes downloaded, speed, and ETA.
    """

    def __init__(self, filename, total_size, initial_size=0, position=1):
        """
        Initialize progress bar.

        Parameters
        ----------
        filename : str
            Name of file being downloaded
        total_size : int
            Total expected file size in bytes
        initial_size : int
            Bytes already downloaded (for resume)
        position : int
            Vertical position for this progress bar (for parallel downloads)
        """
        self.filename = filename
        self.total_size = total_size
        self.initial_size = initial_size
        self.position = position

        # Create tqdm progress bar
        self.pbar = tqdm.tqdm(
            total=total_size,
            initial=initial_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=f"Downloading {filename[:40]}",
            miniters=1,
            position=position,  # Unique position for each parallel download
            leave=False,  # Remove when done
        )

        if initial_size > 0:
            resume_pct = (initial_size / total_size) * 100
            self.pbar.set_postfix_str(f"Resumed from {resume_pct:.1f}%")

    def update(self, bytes_downloaded):
        """Update progress bar with new bytes downloaded."""
        self.pbar.update(bytes_downloaded)

    def close(self):
        """Close progress bar."""
        self.pbar.close()

    def set_postfix(self, text):
        """Set postfix text (e.g., 'Validating...')."""
        self.pbar.set_postfix_str(text)


def validate_and_get_resume_info(
    filepath, url, session, expected_checksum=None, checksum_type=None
):
    """
    Validate existing file and determine if we can resume download.

    Parameters
    ----------
    filepath : str
        Local file path to validate
    url : str
        Remote URL to compare against
    session : requests.Session
        Session object for HTTP requests

    Returns
    -------
    tuple
        (is_complete, resume_from_byte, expected_size)
        - is_complete: True if file is valid and complete
        - resume_from_byte: byte position to resume from (0 = start fresh)
        - expected_size: expected total file size
    """
    if not os.path.exists(filepath):
        expected_size, _supports_range = _probe_remote_file_info(url, session)
        return False, 0, expected_size

    # Check file size is non-zero
    try:
        file_size = os.path.getsize(filepath)
        if file_size == 0:
            LOGGER.warning("Existing file %s has zero size, will restart", filepath)
            return False, 0, 0
    except OSError as e:
        LOGGER.warning("Cannot access file %s: %s, will restart", filepath, e)
        return False, 0, 0

    expected_size, supports_range = _probe_remote_file_info(url, session)

    if expected_size > 0:
        if file_size < expected_size:
            if supports_range:
                LOGGER.info(
                    "Partial file found (%d/%d bytes = %.1f%%), will resume download",
                    file_size,
                    expected_size,
                    (file_size / expected_size) * 100,
                )
                return False, file_size, expected_size
            LOGGER.warning(
                "Partial file found but server doesn't support resume, "
                "will restart from beginning"
            )
            return False, 0, expected_size

        if file_size == expected_size:
            # File appears complete, validate checksum if available
            if expected_checksum and checksum_type:
                LOGGER.info(
                    "Validating %s checksum for %s...",
                    checksum_type.upper(),
                    os.path.basename(filepath),
                )
                computed_checksum = _compute_file_checksum(filepath, checksum_type)

                if computed_checksum is None:
                    LOGGER.warning("Could not compute checksum, skipping validation")
                elif computed_checksum.lower() != expected_checksum.lower():
                    LOGGER.error(
                        "Checksum mismatch for %s!\n"
                        "  Expected: %s\n"
                        "  Got:      %s\n"
                        "File may be corrupted, will re-download",
                        filepath,
                        expected_checksum,
                        computed_checksum,
                    )
                    return False, 0, expected_size
                else:
                    LOGGER.info(
                        "Checksum validation passed (%s)", checksum_type.upper()
                    )

            # Verify with GDAL
            try:
                ds = osgeo.gdal.Open(filepath, osgeo.gdal.GA_ReadOnly)
                if ds is None:
                    LOGGER.warning(
                        "Complete file %s cannot be opened by GDAL, will restart",
                        filepath,
                    )
                    return False, 0, expected_size
                ds = None
                LOGGER.debug(
                    "Existing file %s validated (%d bytes)", filepath, file_size
                )
                return True, file_size, expected_size
            except Exception as e:
                LOGGER.warning(
                    "Error validating file %s with GDAL: %s, will restart",
                    filepath,
                    e,
                )
                return False, 0, expected_size

        LOGGER.warning(
            "Existing file %s is larger than expected (%d > %d bytes), " "will restart",
            filepath,
            file_size,
            expected_size,
        )
        return False, 0, expected_size

    if expected_checksum and checksum_type:
        LOGGER.info(
            "Remote size unavailable for %s; validating with %s checksum instead",
            os.path.basename(filepath),
            checksum_type.upper(),
        )
        computed_checksum = _compute_file_checksum(filepath, checksum_type)
        if computed_checksum is None:
            return False, 0, 0
        if computed_checksum.lower() != expected_checksum.lower():
            LOGGER.warning(
                "Checksum mismatch for %s without remote size metadata, will restart",
                filepath,
            )
            return False, 0, 0

        try:
            ds = osgeo.gdal.Open(filepath, osgeo.gdal.GA_ReadOnly)
            if ds is None:
                return False, 0, 0
            ds = None
            _persist_validation_metadata(
                filepath,
                checksum_type=checksum_type,
                checksum=computed_checksum,
            )
            LOGGER.debug("Existing file %s validated by checksum and GDAL", filepath)
            return True, file_size, file_size
        except Exception as exc:
            LOGGER.warning(
                "Checksum passed but GDAL validation failed for %s: %s",
                filepath,
                exc,
            )
            return False, 0, 0

    local_validation = _load_validation_metadata(filepath)
    if local_validation:
        stored_size = local_validation.get("size")
        stored_checksum_type = local_validation.get("checksum_type")
        stored_checksum = local_validation.get("checksum")

        if stored_size == file_size and stored_checksum_type and stored_checksum:
            LOGGER.info(
                "Remote metadata unavailable for %s; using cached local %s validation",
                os.path.basename(filepath),
                stored_checksum_type.upper(),
            )
            computed_checksum = _compute_file_checksum(filepath, stored_checksum_type)
            if (
                computed_checksum
                and computed_checksum.lower() == stored_checksum.lower()
            ):
                try:
                    ds = osgeo.gdal.Open(filepath, osgeo.gdal.GA_ReadOnly)
                    if ds is None:
                        return False, 0, 0
                    ds = None
                    return True, file_size, file_size
                except Exception as exc:
                    LOGGER.warning(
                        "Cached local validation passed but GDAL failed for %s: %s",
                        filepath,
                        exc,
                    )
                    return False, 0, 0

    LOGGER.warning(
        "Could not verify remote size or checksum for existing file %s; "
        "will restart to avoid trusting a possibly incomplete file",
        filepath,
    )
    return False, 0, 0


def download_file_resumable(
    filepath,
    url,
    session,
    expected_checksum=None,
    checksum_type=None,
    max_retries=3,
    retry_delay=5,
    show_progress=True,
    progress_position=1,
):
    """
    Download a file with resume capability using HTTP Range requests.

    Parameters
    ----------
    filepath : str
        Local path to save the file
    url : str
        URL to download from
    session : requests.Session
        Session object for HTTP requests
    expected_checksum : str, optional
        Expected checksum value (hex string)
    checksum_type : str, optional
        Type of checksum: 'md5' or 'sha256'
    max_retries : int
        Maximum number of retry attempts
    retry_delay : int
        Seconds to wait between retries
    show_progress : bool
        Whether to show per-file progress bar
    progress_position : int
        Vertical position for progress bar (for parallel downloads)

    Returns
    -------
    bool
        True if download successful, False otherwise
    """
    attempt = 0
    while attempt < max_retries:
        attempt += 1

        # Check if we can resume
        is_complete, resume_from, expected_size = validate_and_get_resume_info(
            filepath, url, session, expected_checksum, checksum_type
        )

        if is_complete:
            LOGGER.debug("File already complete: %s", filepath)
            return True

        # If resume_from is 0 and file exists, remove it to start fresh
        if resume_from == 0 and os.path.exists(filepath):
            try:
                os.remove(filepath)
            except OSError as e:
                LOGGER.error("Could not remove file %s: %s", filepath, e)
                return False

        progress = None
        try:
            # Set up headers for resume
            headers = {}
            if resume_from > 0:
                headers["Range"] = f"bytes={resume_from}-"
                LOGGER.info("Resuming download from byte %d", resume_from)

            response = session.get(url, stream=True, headers=headers, timeout=30)

            # Check if server supports range request
            if resume_from > 0 and response.status_code != 206:
                LOGGER.warning(
                    "Server returned %d instead of 206 for range request, "
                    "restarting download",
                    response.status_code,
                )
                if os.path.exists(filepath):
                    os.remove(filepath)
                resume_from = 0
                # Retry without range header
                response = session.get(url, stream=True, timeout=30)

            response.raise_for_status()

            # Open file in append mode if resuming, write mode if starting fresh
            mode = "ab" if resume_from > 0 else "wb"

            # Get content length from this response
            content_length = int(response.headers.get("Content-Length", 0))
            total_size = (
                resume_from + content_length if resume_from > 0 else content_length
            )

            progress_total = expected_size or total_size

            # Create progress bar for this file
            if show_progress and progress_total > 0:
                filename = os.path.basename(filepath)
                progress = DownloadProgressBar(
                    filename, progress_total, resume_from, progress_position
                )

            # Download with real-time progress
            bytes_downloaded = resume_from
            with open(filepath, mode) as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        chunk_size = len(chunk)
                        bytes_downloaded += chunk_size

                        # Update progress bar
                        if progress:
                            progress.update(chunk_size)

            # Close progress bar during validation
            if progress:
                progress.set_postfix("Validating size...")

            # Verify final size
            actual_size = os.path.getsize(filepath)

            if expected_size > 0 and actual_size < expected_size:
                if progress:
                    progress.close()
                LOGGER.warning(
                    "Download incomplete (%d/%d bytes). Attempt %d/%d",
                    actual_size,
                    expected_size,
                    attempt,
                    max_retries,
                )
                if attempt < max_retries:
                    time.sleep(retry_delay)
                    continue  # Will resume in next iteration
                else:
                    return False

            # Checksum validation (if provided)
            if expected_checksum and checksum_type:
                if progress:
                    progress.set_postfix(f"Validating {checksum_type.upper()}...")

                computed = _compute_file_checksum(filepath, checksum_type)
                if computed is None or computed.lower() != expected_checksum.lower():
                    if progress:
                        progress.close()
                    LOGGER.warning(
                        "Checksum validation failed. Attempt %d/%d",
                        attempt,
                        max_retries,
                    )
                    if attempt < max_retries:
                        os.remove(filepath)
                        time.sleep(retry_delay)
                        continue
                    else:
                        return False
                LOGGER.debug("Checksum validation passed")

            # Final GDAL validation
            if progress:
                progress.set_postfix("Validating with GDAL...")

            try:
                ds = osgeo.gdal.Open(filepath, osgeo.gdal.GA_ReadOnly)
                if ds is None:
                    if progress:
                        progress.close()
                    LOGGER.warning(
                        "Downloaded file failed GDAL validation. Attempt %d/%d",
                        attempt,
                        max_retries,
                    )
                    if attempt < max_retries:
                        # Don't delete - we can resume on next attempt
                        time.sleep(retry_delay)
                        continue
                    else:
                        return False
                ds = None

                _persist_validation_metadata(
                    filepath,
                    checksum_type=checksum_type,
                    checksum=computed if expected_checksum and checksum_type else None,
                )
                if progress:
                    progress.set_postfix("Complete ✓")
                    progress.close()
                LOGGER.info("Download complete and validated: %s", filepath)
                return True

            except Exception as e:
                if progress:
                    progress.close()
                LOGGER.warning(
                    "Error validating downloaded file: %s. Attempt %d/%d",
                    e,
                    attempt,
                    max_retries,
                )
                if attempt < max_retries:
                    time.sleep(retry_delay)
                    continue
                else:
                    return False

        except RequestException as e:
            if progress:
                progress.close()
            LOGGER.error(
                "Error downloading %s: %s. Attempt %d/%d", url, e, attempt, max_retries
            )
            if attempt < max_retries:
                time.sleep(retry_delay)
                continue
            else:
                return False

    return False


def download_file_s3_resumable(
    filepath,
    s3_client,
    bucket,
    key,
    url,
    session,
    expected_checksum=None,
    checksum_type=None,
    max_retries=3,
    retry_delay=5,
    show_progress=True,
    progress_position=1,
):
    """
    Download from S3 with resume capability.

    S3 supports range requests through boto3's get_object.

    Parameters
    ----------
    filepath : str
        Local path to save the file
    s3_client : boto3.client
        S3 client object
    bucket : str
        S3 bucket name
    key : str
        S3 object key
    url : str
        HTTP URL for size validation
    session : requests.Session
        Session object for HTTP HEAD requests
    max_retries : int
        Maximum number of retry attempts
    retry_delay : int
        Seconds to wait between retries

    Returns
    -------
    bool
        True if download successful, False otherwise
    """
    attempt = 0
    while attempt < max_retries:
        attempt += 1

        # Check if we can resume
        is_complete, resume_from, expected_size = validate_and_get_resume_info(
            filepath, url, session, expected_checksum, checksum_type
        )

        if is_complete:
            LOGGER.debug("File already complete: %s", filepath)
            return True

        progress = None
        try:
            progress_total = expected_size
            head_response = None

            if resume_from > 0:
                # S3 resumable download using get_object with Range
                LOGGER.info("Resuming S3 download from byte %d", resume_from)

                # Get object metadata to determine total size
                head_response = s3_client.head_object(Bucket=bucket, Key=key)
                total_size = head_response["ContentLength"]
                progress_total = progress_total or total_size

                if show_progress and progress_total > 0:
                    filename = os.path.basename(filepath)
                    progress = DownloadProgressBar(
                        filename, progress_total, resume_from, progress_position
                    )

                # Download remaining bytes
                response = s3_client.get_object(
                    Bucket=bucket, Key=key, Range=f"bytes={resume_from}-"
                )

                # Append to existing file
                with open(filepath, "ab") as f:
                    for chunk in response["Body"].iter_chunks(chunk_size=8192):
                        f.write(chunk)
                        if progress:
                            progress.update(len(chunk))

            else:
                # Fresh S3 download with progress
                if os.path.exists(filepath):
                    os.remove(filepath)

                if progress_total <= 0:
                    head_response = s3_client.head_object(Bucket=bucket, Key=key)
                    progress_total = head_response["ContentLength"]

                if show_progress and progress_total > 0:
                    filename = os.path.basename(filepath)
                    progress = DownloadProgressBar(
                        filename, progress_total, resume_from, progress_position
                    )

                if progress:
                    # Use callback for progress
                    def progress_callback(bytes_amount):
                        progress.update(bytes_amount)

                    s3_client.download_file(
                        bucket, key, filepath, Callback=progress_callback
                    )
                else:
                    s3_client.download_file(bucket, key, filepath)

            # Update progress for validation
            if progress:
                progress.set_postfix("Validating size...")

            # Verify download
            actual_size = os.path.getsize(filepath)
            if expected_size > 0 and actual_size < expected_size:
                if progress:
                    progress.close()
                LOGGER.warning(
                    "S3 download incomplete (%d/%d bytes). Attempt %d/%d",
                    actual_size,
                    expected_size,
                    attempt,
                    max_retries,
                )
                if attempt < max_retries:
                    time.sleep(retry_delay)
                    continue
                else:
                    return False

            # Checksum validation
            if expected_checksum and checksum_type:
                if progress:
                    progress.set_postfix(f"Validating {checksum_type.upper()}...")

                computed = _compute_file_checksum(filepath, checksum_type)
                if computed is None or computed.lower() != expected_checksum.lower():
                    if progress:
                        progress.close()
                    LOGGER.warning(
                        "S3 checksum validation failed. Attempt %d/%d",
                        attempt,
                        max_retries,
                    )
                    if attempt < max_retries:
                        os.remove(filepath)
                        time.sleep(retry_delay)
                        continue
                    else:
                        return False

            # GDAL validation
            if progress:
                progress.set_postfix("Validating with GDAL...")

            try:
                ds = osgeo.gdal.Open(filepath, osgeo.gdal.GA_ReadOnly)
                if ds is None:
                    if progress:
                        progress.close()
                    LOGGER.warning("S3 download failed GDAL validation")
                    if attempt < max_retries:
                        time.sleep(retry_delay)
                        continue
                    return False
                ds = None

                _persist_validation_metadata(
                    filepath,
                    checksum_type=checksum_type,
                    checksum=computed if expected_checksum and checksum_type else None,
                )
                if progress:
                    progress.set_postfix("Complete ✓")
                    progress.close()

                LOGGER.info("S3 download complete and validated: %s", filepath)
                return True
            except Exception as e:
                if progress:
                    progress.close()
                LOGGER.warning("Error validating S3 download: %s", e)
                if attempt < max_retries:
                    time.sleep(retry_delay)
                    continue
                return False

        except Exception as exc:
            if progress:
                progress.close()
            LOGGER.warning("S3 download attempt %d failed: %s", attempt, exc)
            if attempt < max_retries:
                time.sleep(retry_delay)
                continue
            return False

    return False


def get_url_ifg(scenes):
    """Get url, ifg of fetched ASF scene"""
    urls, ifgs = [], []
    for scene in scenes:
        s = scene.geojson()["properties"]
        urls.append(s["url"])
        # NISAR files are formatted differently
        if s["fileID"].startswith("NISAR_"):
            f = s["fileID"].split("_")
            pairname = f[11][:8] + "_"
            pairname += f[13][:8]
            ifgs.append(pairname)
        else:
            f = s["fileID"].split("-")
            pairname = f[6]
            ifgs.append(pairname)

    # determine if NISAR GUNW
    is_nisar_file = False
    if urls != []:
        if "/NISAR_" in urls[0]:
            is_nisar_file = True

    return urls, ifgs, is_nisar_file


def fmt_dst(args):
    """Format the save name"""
    ext = ".kmz" if args.output == "Kml" else ".txt"

    if args.track is not None:
        fn_track = f"track{args.track}".replace(",", "-")
    else:
        fn_track = ""

    if args.bbox is not None:
        WSEN = make_bbox(args.bbox).bounds
        WSEN_fmt = []
        for i, coord in enumerate(WSEN):
            if i < 2:
                WSEN_fmt.append(math.floor(float(coord)))
            else:
                WSEN_fmt.append(math.ceil(float(coord)))
        fn_bbox = f"_bbox{WSEN_fmt[0]}W{WSEN_fmt[1]}S{WSEN_fmt[2]}E{WSEN_fmt[3]}N"
    else:
        fn_bbox = ""

    dst = os.path.join(args.wd, f"{fn_track}{fn_bbox}_0{ext}".lstrip("_"))
    count = 1  # don't overwrite if already exists
    while os.path.exists(dst):
        basen = (
            f"{re.split(str(count-1)+ext, os.path.basename(dst))[0]}" f"{count}{ext}"
        )
        dst = os.path.join(os.path.dirname(dst), basen)
        count += 1
    return dst


class Downloader:
    """Product Downloading Class."""

    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.args.output = self.args.output.title()
        self.args.wd = os.path.abspath(self.args.wd)
        LOGGER.setLevel(logging.DEBUG if self.args.verbose else logging.INFO)

    def __call__(self):
        scenes = self.query_asf()
        urls, ifgs, is_nisar_file = get_url_ifg(scenes)

        # Subset everything by version
        if is_nisar_file and self.args.version is not None:
            raise Exception(
                "Version support not included for NISAR, remove the critera"
            )
        else:
            urls = url_versions(urls, self.args.version, self.args.wd)
        scenes = [scene for scene, url in zip(scenes, urls) if url in urls]
        ifgs = [ifg for ifg, url in zip(ifgs, urls) if url in urls]

        # Filter scenes based on date and elapsed time criteria
        scenes, urls, ifgs = self.filter_scenes(scenes, urls, ifgs, is_nisar_file)

        if self.args.output == "Count":
            LOGGER.info("Found -- %d -- products", len(scenes))
        elif self.args.output == "Url":
            self.write_urls(urls, scenes)
        elif self.args.output == "Download":
            self.download_scenes(scenes)

        if self.args.verbose:
            for scene in scenes:
                LOGGER.info(scene.geojson()["properties"]["sceneName"])

    def query_asf(self):
        """Query ASF for scenes."""
        bbox = make_bbox(self.args.bbox)
        bbox_wkt = bbox.wkt if bbox else None

        flight_direction = None
        if self.args.flightdir:
            flight_direction = (
                "ascending"
                if self.args.flightdir.lower().startswith("a")
                else "descending"
            )

        tracks = (
            [int(track) for track in self.args.track.split(",")]
            if self.args.track
            else None
        )

        start = self.args.start - datetime.timedelta(days=1)
        end = self.args.end + datetime.timedelta(days=1)

        if self.args.mission.upper() == "S1":
            return asf_search.geo_search(
                collections=["C2859376221-ASF", "C1261881077-ASF"],
                dataset=asf_search.constants.ARIA_S1_GUNW,
                processingLevel=asf_search.constants.GUNW_STD,
                relativeOrbit=tracks,
                flightDirection=flight_direction,
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
            )
        elif self.args.mission.upper() == "NISAR":
            # Authenticate so the private ephemeral archive
            # collection (C4052499921-ASF) is visible in CMR.
            session = self._get_asf_session()
            opts = asf_search.ASFSearchOptions(
                collections=[
                    "C2850261892-ASF",  # public NISAR GUNW
                    "C4052499921-ASF",  # private ephemeral archive
                ],
                dataset=asf_search.constants.NISAR,
                processingLevel=asf_search.constants.GUNW,
                relativeOrbit=tracks,
                flightDirection=flight_direction,
                intersectsWith=bbox_wkt,
                start=start,
                end=end,
                session=session,
            )
            return asf_search.geo_search(opts=opts)

    def _get_asf_session(self):
        """Return an authenticated ASFSession.

        Uses explicit user/pass args when provided, otherwise falls
        back to ~/.netrc credentials for urs.earthdata.nasa.gov.
        Returns an unauthenticated session with a warning when no
        credentials are available.
        """
        session = asf_search.ASFSession()
        if self.args.user:
            session.auth_with_creds(
                self.args.user,
                self.args.passw or getpass.getpass("NASA Earthdata password: "),
            )
        else:
            try:
                import netrc as _netrc

                nrc = _netrc.netrc()
                auth = nrc.authenticators("urs.earthdata.nasa.gov")
                if auth:
                    session.auth_with_creds(auth[0], auth[2])
                else:
                    LOGGER.warning(
                        "No urs.earthdata.nasa.gov entry in ~/.netrc. "
                        "Private collections (e.g. NISAR ephemeral "
                        "archive) will not be visible."
                    )
            except FileNotFoundError:
                LOGGER.warning(
                    "~/.netrc not found. Private collections (e.g. "
                    "NISAR ephemeral archive) will not be visible."
                )
        return session

    def filter_scenes(self, scenes, urls, ifgs, is_nisar_file):
        filtered_scenes, filtered_urls, filtered_ifgs = [], [], []

        for scene, url, ifg in zip(scenes, urls, ifgs):
            eni, sti = self.parse_dates(ifg, is_nisar_file)
            if self.args.ifg:
                if self.match_single_ifg(sti, eni):
                    filtered_scenes.append(scene)
                    filtered_urls.append(url)
                    filtered_ifgs.append(ifg)
            elif self.match_date_criteria(sti, eni):
                filtered_scenes.append(scene)
                filtered_urls.append(url)
                filtered_ifgs.append(ifg)
        return filtered_scenes, filtered_urls, filtered_ifgs

    def parse_dates(self, ifg, is_nisar_file):
        if is_nisar_file:
            sti, eni = (datetime.datetime.strptime(d, "%Y%m%d") for d in ifg.split("_"))
        else:
            eni, sti = (datetime.datetime.strptime(d, "%Y%m%d") for d in ifg.split("_"))
        return eni, sti

    def match_single_ifg(self, sti, eni):
        dates = [
            datetime.datetime.strptime(i, "%Y%m%d").date()
            for i in self.args.ifg.split("_")
        ]
        st1, en1 = sorted(dates)
        return st1 == sti.date() and en1 == eni.date()

    def match_date_criteria(self, sti, eni):
        sten_chk = sti >= self.args.start and eni <= self.args.end
        elap = (eni - sti).days
        elap_chk = self.args.daysgt <= elap <= self.args.dayslt
        return sten_chk and elap_chk

    def write_urls(self, urls, scenes):
        os.makedirs(self.args.wd, exist_ok=True)
        dst = fmt_dst(self.args)
        with open(dst, "w") as fh:
            for url, scene in zip(urls, scenes):
                s3_url = _get_s3_data_url(scene) or ""
                print(f"{url},{s3_url}", file=fh)
        LOGGER.info("Wrote -- %d -- product urls to: %s", len(urls), dst)

    def download_scenes(self, scenes):
        os.makedirs(self.args.wd, exist_ok=True)
        scenes = asf_search.ASFSearchResults(scenes)
        nt = int(self.args.num_threads)
        LOGGER.info("Downloading %d products...", len(scenes))

        # Check if we can use S3 direct download (on AWS)
        use_s3 = ARIAtools.util.s3.is_on_aws()
        s3_client = None
        if use_s3:
            # Determine credential endpoint from first S3 URL
            first_s3 = next(
                (_get_s3_data_url(s) for s in scenes if _get_s3_data_url(s)), None
            )
            endpoint_key = (
                ARIAtools.util.s3._endpoint_key_for_s3uri(first_s3)
                if first_s3
                else "default"
            )
            try:
                s3_client = ARIAtools.util.s3.get_s3_client(
                    endpoint_key, max_pool_connections=nt * 10
                )
                LOGGER.info("Using S3 direct download (endpoint: %s)", endpoint_key)
            except Exception as exc:
                LOGGER.warning(
                    "S3 client setup failed, falling back " "to HTTPS: %s", exc
                )
                use_s3 = False

        session = asf_search.ASFSession()
        if self.args.user:
            session.auth_with_creds(self.args.user, self.args.passw)

        # Track progress bar positions for parallel downloads
        import threading

        position_lock = threading.Lock()
        position_counter = [1]  # Start at 1 (position 0 is overall progress)
        active_positions = {}  # Maps thread_id -> position

        def get_progress_position():
            """Get a unique position for this thread's progress bar."""
            thread_id = threading.get_ident()
            with position_lock:
                if thread_id in active_positions:
                    return active_positions[thread_id]
                pos = position_counter[0]
                position_counter[0] += 1
                active_positions[thread_id] = pos
                return pos

        def release_progress_position():
            """Release this thread's progress bar position."""
            thread_id = threading.get_ident()
            with position_lock:
                if thread_id in active_positions:
                    del active_positions[thread_id]

        def download_file(scene, max_retries=3, retry_delay=5):
            """Download a single scene with resume support and validation."""
            url = scene.properties["url"]
            local_filename = url.split("/")[-1]
            filepath = os.path.join(self.args.wd, local_filename)

            # Get unique progress bar position for this thread
            progress_pos = get_progress_position()

            try:
                # Extract checksum from scene metadata
                expected_checksum, checksum_type = _get_scene_checksum(scene)
                if expected_checksum:
                    LOGGER.debug(
                        "Expected %s checksum: %s",
                        checksum_type.upper(),
                        expected_checksum,
                    )

                # Check if file is already complete and valid
                is_complete, resume_from, expected_size = validate_and_get_resume_info(
                    filepath, url, session, expected_checksum, checksum_type
                )

                if is_complete:
                    LOGGER.info("Valid product already in directory: %s", filepath)
                    return filepath

                # Try S3 download first if available
                s3_url = _get_s3_data_url(scene) if use_s3 else None
                if s3_client and s3_url:
                    try:
                        bucket, key = ARIAtools.util.s3.parse_s3_uri(s3_url)
                        success = download_file_s3_resumable(
                            filepath,
                            s3_client,
                            bucket,
                            key,
                            url,
                            session,
                            expected_checksum,
                            checksum_type,
                            max_retries,
                            retry_delay,
                            show_progress=True,
                            progress_position=progress_pos,
                        )
                        if success:
                            return filepath
                        else:
                            LOGGER.warning(
                                "S3 download failed, falling back to HTTPS for %s",
                                local_filename,
                            )
                    except Exception as exc:
                        LOGGER.warning(
                            "S3 download error: %s, falling back to HTTPS", exc
                        )

                # HTTPS download (default or fallback) with resume support
                success = download_file_resumable(
                    filepath,
                    url,
                    session,
                    expected_checksum,
                    checksum_type,
                    max_retries,
                    retry_delay,
                    show_progress=True,
                    progress_position=progress_pos,
                )

                return filepath if success else None

            finally:
                # Release progress bar position
                release_progress_position()

        # Create overall progress bar
        success_count = 0
        failure_count = 0
        pbar = tqdm.tqdm(
            total=len(scenes),
            unit="file",
            desc="Overall Progress",
            position=0,
            leave=True,
        )
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=nt) as executor:
                future_to_scene = {
                    executor.submit(download_file, scene): scene for scene in scenes
                }
                for future in concurrent.futures.as_completed(future_to_scene):
                    scene = future_to_scene[future]
                    try:
                        filepath = future.result()
                        if filepath is not None:
                            success_count += 1
                            LOGGER.debug("Downloaded: %s", filepath)
                        else:
                            failure_count += 1
                            LOGGER.error(
                                "Failed to download: %s",
                                scene.properties["url"],
                            )
                        pbar.update(1)
                    except Exception as exc:
                        failure_count += 1
                        LOGGER.error(
                            "%s generated an exception: %s",
                            scene.properties["url"],
                            exc,
                        )
                        pbar.update(1)
        finally:
            pbar.close()

        LOGGER.info(
            "Download complete. Wrote -- %d/%d -- products to: %s",
            success_count,
            len(scenes),
            self.args.wd,
        )
        if failure_count:
            LOGGER.warning("Download finished with %d failed products.", failure_count)


def main():
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)

    print("*****************************************************************")
    LOGGER.info("*** Download Function ***")
    print("*****************************************************************")

    # format dates
    args.start = datetime.datetime.strptime(args.start, "%Y%m%d")
    args.end = datetime.datetime.strptime(args.end, "%Y%m%d")

    if not args.track and not args.bbox:
        raise Exception("Must specify either a bbox or track")
    Downloader(args)()


if __name__ == "__main__":
    main()
