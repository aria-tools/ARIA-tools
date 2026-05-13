# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import logging
import os

LOGGER = logging.getLogger(__name__)


def _parse_version_token(version_token):
    """Return a comparable tuple for version strings like ``v2.0.0``."""

    normalized = version_token.lstrip("v").rstrip(".nc").replace("_", ".")
    return tuple(int(part) for part in normalized.split("."))


# Grab older version products if specified.
def url_versions(urls, user_version, wd):
    """For duplicate products (other than version number)

    Uses the the latest if user_version is None else use specified ver.
    """
    if user_version is not None and str(user_version).lower() not in ("all", "none"):
        user_version = user_version.lstrip("v")
        if not user_version[0].isdigit():
            raise Exception(
                f"Input version {user_version} not in format X* "
                "as expected (e.g., 3, 3_0_0)"
            )
        LOGGER.debug(f"Only using products version: {user_version}")
        urls_final = [url for url in urls if f"-v{user_version}" in url]
        if not urls_final:
            raise Exception(f"No products with user specified version: {urls_final}")
    else:
        urls_final = urls
    return urls_final


# Currently does not work as expected, as lat/lon coords of newer versions
# do not match older. Will only support specific versions or 'all'
def url_versions_full(urls, user_version, wd):
    """For duplicate products (other than version number)
    Uses the the latest if user_version is None else use specified ver.

    Optimized O(n) implementation using dict-based grouping instead of
    nested loops (was O(n²) before optimization).
    """
    if isinstance(user_version, str) and user_version.lower() == "all":
        return urls

    # Group URLs by base name in single pass - O(n)
    url_groups = {}
    for url in urls:
        url_base = "-".join(url.split("-")[:-1])
        if url_base not in url_groups:
            url_groups[url_base] = []
        url_groups[url_base].append(url)

    # Process each group to select version - O(n) total
    urls_final = []
    for url_base, duplicates in url_groups.items():
        if len(duplicates) == 1:
            urls_final.append(duplicates[0])
        else:
            parsed_duplicates = [
                (_parse_version_token(dupe.split("-")[-1]), dupe) for dupe in duplicates
            ]

            if user_version is None:
                _, selected_url = max(parsed_duplicates, key=lambda item: item[0])
            else:
                requested_version = _parse_version_token(str(user_version))
                matches = [
                    dupe
                    for version, dupe in parsed_duplicates
                    if version == requested_version
                ]
                if not matches:
                    raise Exception(
                        f"No products with user specified version: {user_version}"
                    )
                selected_url = matches[0]

            urls_final.append(selected_url)

            # move duplicates to a different folder
            dupe_folder = os.path.join(wd, "duplicated_products")
            os.makedirs(dupe_folder, exist_ok=True)
            for dupe in duplicates:
                dupe_path = os.path.join(dupe_folder, os.path.basename(dupe))
                wd_path = os.path.join(wd, os.path.basename(dupe))
                if os.path.basename(dupe) != os.path.basename(
                    urls_final[-1]
                ) and os.path.exists(wd_path):
                    os.rename(wd_path, dupe_path)
    return urls_final
