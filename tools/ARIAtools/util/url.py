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
            versions = []
            for dupe in duplicates:
                ver_str = dupe.split("-")[-1]
                versions.append(float(ver_str.lstrip("v").rstrip(".nc")))

            # default, use latest version
            if user_version is None:
                version = str(max(versions))
                version = f"v{version[0]}_{version[1]}_{version[2]}.nc"
            else:
                version = f"v{user_version}.nc"

            urls_final.append(f"{url_base}-{version}")

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
