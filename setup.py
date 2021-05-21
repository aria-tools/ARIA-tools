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


# Path where relax third party is download to
# Note Min-Cost-Flow-Class has its own liscence agreement, users to ensure proper use of third party license agreements.
root_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tools')
bind_path  = os.path.join(root_path, 'bindings')
bind_dirs  = [os.path.join(bind_path, rdir) for \
                        rdir in 'relaxIVdriver.cpp unwcompmodule.cpp'.split()]

relax_path = os.path.join(root_path, 'thirdParty', 'Min-Cost-Flow-Class')
relax_dirs = [os.path.join(relax_path, rdir) for \
                        rdir in 'MCFClass OPTUtils RelaxIV'.split()]

aria_scripts = [os.path.join(root_path, 'bin', script) for \
                script in 'ariaPlot.py ariaDownload.py ariaExtract.py '\
                   'ariaTSsetup.py ariaAOIassist.py ariaMisclosure.py'.split()]

## setup version; adapted from MintPy
cmd = "git describe --tags"
version  = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
version  = version.decode('utf-8').strip()


# if there are new commits after the latest release
if '-' in version:
    version, num_commit = version.split('-')[:2]
    version0 = version
    version += f'-{num_commit}'
else:
    version0 = version

cmd  = f"git log -1 --date=short --format=%cd {version0}"
date = subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
date = date.decode('utf-8').strip()
print (f'ARIA-tools {version} {date}')

# If third party package is found compile as well
if os.path.isdir(relax_path):
    print('Installing ARIA-tools with support for RelaxIV')
    module1 = Extension('ARIAtools.demo', sources=bind_dirs + 
                     [os.path.join(relax_path, 'RelaxIV', 'RelaxIV.C')], 
             include_dirs=[os.path.join(root_path, 'include')] + relax_dirs)

    setup (name = 'ARIAtools',
           version = version0,
           description = 'This is the ARIA tools package with RelaxIV support',
           ext_modules = [module1],
           packages=['ARIAtools'],
           package_dir={'ARIAtools': os.path.join(root_path, 'ARIAtools')},
           scripts=aria_scripts)
else:
    # Third party package RelaxIV not found
    print('Installing ARIA-tools without support for RelaxIV')

    setup (name = 'ARIAtools',
           version = version0,
           description = 'This is the ARIA tools package without RelaxIV support',
           packages=['ARIAtools'],
           package_dir={'ARIAtools': os.path.join(root_path, 'ARIAtools')},
           scripts=aria_scripts)
