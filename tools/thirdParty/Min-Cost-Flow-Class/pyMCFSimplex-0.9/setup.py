#!/usr/bin/env python

"""
setup.py file for SWIG Wrapper for MCFSimplex
"""

from distutils.core import setup, Extension
import numpy as np # for numpy array conversion

pyMCFSimplex_module = Extension('_pyMCFSimplex',
                           sources=['pyMCFSimplex_wrap.cxx', 'MCFSimplex.cpp'],
                           )

setup (name = 'pyMCFSimplex',
       version = '0.9',
       author      = "G#.Blog - Johannes Sommer",
       author_email = "info@sommer-forst.de",
       url = r"http:\\www.sommer-forst.de\blog",
       description = "pyMCFSimplex is a Python Wrapper for MCFSimplex",
       long_description = 
"""
pyMCFSimplex is a Python Wrapper for the Minimum Cost Flow Problem Solver 
'MCFSimplex' coded and maintained at the 
Operations Research Group of the University of Pisa.
""",
       include_dirs = [np.get_include()], # Header for numpy
       ext_modules = [pyMCFSimplex_module],
       license = "LGPL 2.1",
       platforms = ["win32","linux-x86_64"],
       py_modules = ["pyMCFSimplex"],
       )
