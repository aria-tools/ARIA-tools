# ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/aria-tools/ARIA-tools/blob/master/LICENSE)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b69726c41ff9498882088afe57fe8391)](https://app.codacy.com/app/ehavazli/ARIA-tools?utm_source=github.com&utm_medium=referral&utm_content=ehavazli/ARIA-tools&utm_campaign=Badge_Grade_Settings)
[![CircleCI](https://circleci.com/gh/ehavazli/ARIA-tools.svg?style=svg)](https://circleci.com/gh/ehavazli/ARIA-tools)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. This software is open source under the terms of the GNU General Public License. Its development was funded under the NASA Sea-level Change Team (NSLCT) program and the Earth Surface and Interior (ESI) program.

For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [MintPy](https://github.com/insarlab/MintPy).
<p align="center">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/Hawaii.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/CA.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/EastCoast.png">
</p>
THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents

1.  [Software Dependencies](#software-dependencies)
2.  [Installation](#installation)
3.  [Running ARIA-tools](#running-aria-tools)
-   [Commandline download of GUNW Products](#commandline-download-of-gunw-products)
-   [Manipulating GUNW Products](#manipulating-gunw-products)
-   [Baseline and quality control plots for GUNW Products](#baseline-and-quality-control-plots-for-gunw-products)
-   [Time-series set-up of GUNW Products](#time-series-set-up-of-gunw-products)
4.  [Documentation](#documentation)
5.  [Citation](#citation)
6.  [Contributors and community contributions](#contributors)

------

## Software Dependencies
Below we list the dependencies for ARIA-tools

### Packages
```
* Python >= 3.5  (3.6 preferred)
* [PROJ 4](https://github.com/OSGeo/proj) github) >= 6.0
* [GDAL](https://www.gdal.org/) and its Python bindings >= 3.0
```

### Python dependencies
```
* [SciPy](https://www.scipy.org/)
* [netcdf4](http://unidata.github.io/netcdf4-python/netCDF4/index.html)
* [requests](https://2.python-requests.org/en/master/)
```

### Python Jupyter dependencies
```
* py3X-jupyter
* py3X-jupyter_client
* py3X-jupyter_contrib_nbextensions
* py3X-jupyter_nbextensions_configurator
* py3X-hide_code
* py3X-RISE
```

### Optional Third-party packages
```
* RelaxIV available from [Min-Cost-Flow-Class](https://github.com/frangio68/Min-Cost-Flow-Class)
```

------
## Installation
ARIA-tools package can be easily installed and used after the dependencies are installed and activated. The third-party RelaxIV package is optional (not required), and  only used when opting to minimizing phase-discontinuities. Prior to use of RelaxIV, users should conform to the RelaxIV license agreement. Easiest way of installing RelaxIV is by downloading the min-cost-flow repository in the third-party folder of the ARIAtools and use the setup.py script as outlined below. For the required dependencies, we strongly recommend using [Anaconda](https://www.anaconda.com/distribution/) package manager for easy installation of dependencies in the python environment.

Below we outline the different steps for setting up the ARIA-tools while leveraging Anaconda for installation of the requirements. Running the commands below will clone the ARIA-tools package to your local directory, create a conda environment with the name 'ARIA-tools', install dependencies to this environment and activate it.

```.tcsh
git clone https://github.com/aria-tools/ARIA-tools.git
conda env create -f ./ARIA-tools/environment.yml
conda activate ARIA-tools
```

We have included a setup.py script which allows for easy compilation and installation of third-party dependencies (c-code), as well as setting up the ARIA-tools package itself (python and command line tools).
```.tcsh
python setup.py build
python setup.py install
```

If not using the setup.py, users should compile third-party packages manually and ensure ARIA-tools and dependencies are included on their PATH and PYTHONPATH. For c-shell this can be done as follows (replace "ARIAtoolsREPO" to the location where you have cloned the ARIAtools repository):
```.tcsh
setenv PYTHONPATH $PYTHONPATH:/ARIAtoolsREPO/tools/ARIAtools
set PATH $PATH:'/ARIAtoolsREPO/tools/bin'
```

### Other installation options
The following pages might be of use to those trying to build third party packages from source.
-   [Installing dependencies from source on linux](https://github.com/aria-tools/ARIA-tools/blob/master/Linux_source_build.md)
-   [Installing dependencies from source on mac](https://github.com/aria-tools/ARIA-tools/blob/master/MacOS_source_build.md)

------
## Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/aria-tools/ARIA-tools-docs). We welcome the community to contribute other examples on how to leverage the ARIA-tools (see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md) for instructions).

### Commandline download of GUNW Products
GUNW products can be downloaded through the commandline using the *ariaDownload.py* program, which wraps around the ASF DAAC api.

### Manipulating GUNW Products
GUNW product can be manipulated (cropped, stitched, extracted) using the *ariaExtract.py* program.

### Baseline and quality control plots for GUNW Products
Quality and baseline plots for spatial-temporal contiguous interferograms can be made using the *ariaPlot.py* program.

### Time-series set-up of GUNW Products
Time-series set-up with spatial-temporal contiguous unwrapped interferograms and coherence can be done using the *ariaTSsetup.py* program.

------
## Documentation

See the [ARIA-tools-docs repository](https://github.com/aria-tools/ARIA-tools-docs) for all documentation and Jupyter Notebook Tutorials.

------
## Citation
D. Bekaert, M. Karim, L. Justin, H. Hua, P. Agram, S. Owen, G. Manipon, N. Malarout, M. Lucas, G. Sacco, L. Pan, S. Sangha, and ARIA team (2019), *Development and Dissemination of Standardized Geodetic Products by the Advanced Rapid Imaging and Analysis (ARIA) Center for Natural Hazards*, The International Union of Geodesy and Geophysics (IUGG), Montreal

------
## Contributors
-   David Bekaert
-   Simran Sangha
-   Emre Havazli
-   Brett Buzzanga
-   [_other community members_](https://github.com/aria-tools/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md).
