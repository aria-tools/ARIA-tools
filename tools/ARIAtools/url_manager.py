#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os, os.path as op

###Capture and segregate duplicate products (if applicable)
def url_versions(urls, user_version, wd):
    """For duplicate products (other than version number)
    Uses the the latest if user_version is None else use specified ver."""
    if isinstance(user_version, str) and user_version.lower() == 'all':
        return urls

    url_bases = []
    for url in urls:
        url_base = '-'.join(url.split('-')[:-1])
        if not url_base in url_bases:
            url_bases.append(url_base)

    urls_final = []
    for url_base in url_bases:
        # get the possible matches
        duplicates = [url for url in urls if url_base in url]
        if len(duplicates) == 1:
            urls_final.append(duplicates[0])
        else:
            versions = []
            for dupe in duplicates:
                ver_str = dupe.split('-')[-1]
                versions.append(float(ver_str.lstrip('v').rstrip('.nc')))

            ## default, use latest version
            if user_version is None:
                version = str(max(versions))
                version = f'v{version[0]}_{version[1]}_{version[2]}.nc'
            else:
                version=f'v{user_version}.nc'

            urls_final.append(f'{url_base}-{version}')

            ## move duplicates to a different folder
            dupe_folder = op.join(wd, 'duplicated_products')
            os.makedirs(dupe_folder, exist_ok=True)
            for dupe in duplicates:
                dupe_path = os.path.join(dupe_folder, os.path.basename(dupe))
                wd_path = os.path.join(wd, os.path.basename(dupe))
                if os.path.basename(dupe) != \
                     os.path.basename(urls_final[-1]) and \
                     os.path.exists(wd_path):
                    os.rename(wd_path, dupe_path)
    return urls_final
