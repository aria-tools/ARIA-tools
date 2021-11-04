#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: David Bekaert & Piyush Agram
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

from distutils.core import setup, Extension
import os
import subprocess
from datetime import datetime


# Path where relax third party is download to
# Note Min-Cost-Flow-Class has its own liscence agreement, users to ensure proper use of third party license agreements.
relaxPath='tools/thirdParty/Min-Cost-Flow-Class'

## setup version; adapted from MintPy
cmd = "git describe --tags"
try:
    version  = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
    version  = version.decode('utf-8').strip()
except:
    # for circle ci
    version  = 'v0'


# if there are new commits after the latest release
if '-' in version:
    version, num_commit = version.split('-')[:2]
    version0 = version
    version += f'-{num_commit}'
else:
    version0 = version

cmd  = f"git log -1 --date=short --format=%cd {version0}"
try:
    date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
    date = date.decode('utf-8').strip()
except:
    # for circle ci
    date = datetime.strftime(datetime.today(), '%Y-%m-%d')

print (f'ARIA-tools {version} {date}')

# If third party package is found compile as well
if os.path.isdir(relaxPath):
    print('Installing ARIA-tools with support for RelaxIV')
    module1 = Extension('ARIAtools.demo', sources = ['tools/bindings/relaxIVdriver.cpp',
                                                     'tools/bindings/unwcompmodule.cpp',
                                                     os.path.join(relaxPath,'RelaxIV/RelaxIV.C')], include_dirs=['tools/include' ,os.path.join(relaxPath,'MCFClass'),os.path.join(relaxPath,'OPTUtils'),os.path.join(relaxPath,'RelaxIV')])

    setup (name = 'ARIAtools',
           version = version0,
           description = 'This is the ARIA tools package with RelaxIV support',
           ext_modules = [module1],
           packages=['ARIAtools'],
           package_dir={'': 'tools'},
           scripts=['tools/bin/ariaPlot.py','tools/bin/ariaDownload.py','tools/bin/ariaExtract.py','tools/bin/ariaTSsetup.py','tools/bin/ariaAOIassist.py','tools/bin/ariaMisclosure.py'])
else:
    # Third party package RelaxIV not found
    print('Installing ARIA-tools without support for RelaxIV')

    setup (name = 'ARIAtools',
           version = version0,
           description = 'This is the ARIA tools package without RelaxIV support',
           packages=['ARIAtools'],
           package_dir={'': 'tools'},
           scripts=['tools/bin/ariaPlot.py','tools/bin/ariaDownload.py','tools/bin/ariaExtract.py','tools/bin/ariaTSsetup.py','tools/bin/ariaAOIassist.py','tools/bin/ariaMisclosure.py'])
