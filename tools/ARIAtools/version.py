#!/usr/bin/env python3
# grab version / date of the latest commit

import os
import subprocess


###########################################################################

## have to figure out where the ariaTools directory with setup.py probably
version = 'v1.1.1'
date    = '2021-01-01'
def get_release_info(version, date):
    """Grab version and date of the latest commit from a git repository"""
    # go to the repository directory
    breakpoint()
    dir_orig = os.getcwd()
    os.chdir(os.path.dirname(os.path.dirname(__file__)))

    # grab git info into string
    try:
        cmd = "git describe --tags"
        version = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        version = version.decode('utf-8').strip()

        #if there are new commits after the latest release
        if '-' in version:
            version, num_commit = version.split('-')[:2]
            version += '-{}'.format(num_commit)

        cmd = "git log -1 --date=short --format=%cd"
        date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
        date = date.decode('utf-8').strip()

    except:
        pass

    # go back to the original directory
    os.chdir(dir_orig)
    return version, date


###########################################################################

release_version, release_date = get_release_info(version, date)

website = 'https://github.com/aria-tools/ARIA-tools'

description = 'Advandced Rapid Imaging and Analysis (ARIA) Tools for '\
                                'manipulating ARIA standard InSAR products.'

release_description = f'ARIAtools release version {release_version} '\
                      f'release date {release_date}'
