#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert & Piyush Agram
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from distutils.core import setup, Extension
import os
import subprocess
from setuptools import find_packages
from datetime import datetime

# setup version; adapted from MintPy
cmd = "git describe --tags"
try:
    version = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
    version = version.decode('utf-8').strip()
except BaseException:
    # for circle ci
    version = 'v0'


# if there are new commits after the latest release
if '-' in version:
    version, num_commit = version.split('-')[:2]
    version0 = version
    version += f'-{num_commit}'
else:
    version0 = version

cmd = f"git log -1 --date=short --format=%cd {version0}"
try:
    date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
    date = date.decode('utf-8').strip()
except BaseException:
    # for circle ci
    date = datetime.strftime(datetime.today(), '%Y-%m-%d')

print(f'ARIA-tools {version} {date}')
print('Installing ARIA-tools')

setup(name='ARIAtools',
      version=version0,
      description='Installing ARIA-tools',
      packages=find_packages(where='tools'),
      package_dir={'': 'tools'},
      scripts=['tools/bin/ariaPlot.py', 'tools/bin/ariaDownload.py',
               'tools/bin/ariaExtract.py', 'tools/bin/ariaTSsetup.py',
               'tools/bin/ariaAOIassist.py', 'tools/bin/ariaMisclosure.py',
               'tools/bin/export_product.py'])
