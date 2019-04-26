# ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/dbekaert/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. This software is open source under the terms of the GNU General Public License. Its development was funded under the NASA Sea-level Change Team (NSLCT) program and the Earth Surface and Interior (ESI) program. 


For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [PySAR](https://github.com/yunjunz/PySAR). 


THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents

1. [Download](#download)
2. [Software Dependencies](#software-dependencies)
   - [Installing dependencies on linux with Anaconda](#installing-dependencies-on-linux-with-anaconda)
   - [Installing dependencies on mac with Macports](#installing-dependencies-on-mac-with-macports)
   - [Installing dependencies on mac with Anaconda](#installing-dependencies-on-mac-with-anaconda)   
3. [Installation](#installation)
4. [Running ARIA-tools](#running-aria-tools)
   - [Manipulating GUNW Products](#manipulating-gunw-products)
   - [Baseline and quality control plots for GUNW Products](#baseline-and-quality-control-plots-for-gunw-products)
   - [Time-series set-up of GUNW Products](#time-series-set-up-of-gunw-products)
5. [Documentation](#documentation)
6. [Citation](#citation)
7. [Contributors and community contributions](#contributors)


------

## Download

To download the ARIA-tools repository:
```
git clone https://github.com/dbekaert/ARIA-tools.git
```

------

## Software Dependencies
Below we provide guidelines for Linux and Max users on how to build the required dependencies from source as not all packages are yet availble through third-party package managers. In future we will expand this section to rely on anaconda and macports for installation.

### Packages:

```
* Python >= 3.5  (3.6 preferred)
* [PROJ 4](https://github.com/OSGeo/proj) github) >= 6.0
* [GDAL](https://www.gdal.org/) and its Python bindings >= 2.5
```

### Python dependencies
```
* [SciPy](https://www.scipy.org/) 
* [netcdf4](http://unidata.github.io/netcdf4-python/netCDF4/index.html)
```

### [Installing dependencies on linux with Anaconda](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md)
### [Installing dependencies on mac with macports](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md)
### [Installing dependencies on mac with Anaconda](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_Anaconda_source_build.md) 	

------
## Installation

The ARIA-tools package can be installed by adding the tools folder on your PATH and PYTHONPATH.
This can be done by editing your private module or your favorite start-up shell.


For example, for csh do:
```
vi ~/.cshrc
```

Add the following and update *my* to the location where you cloned ARIA-tools:
```
setenv PYTHONPATH $PYTHONPATH:/my/tools
set path = ('/my/tools/python' $path)
```

------
## Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs/blob/master/README.md). We welcome the community to contribute other examples on how to leverage the ARIA-tools (see [here](https://github.com/dbekaert/ARIA-tools/blob/master/CONTRIBUTING.md) for instructions).

### Manipulating GUNW Products
GUNW product can be manipulated (cropped, stitched, extracted) using the *extractProduct.py* program. 

### Baseline and quality control plots for GUNW Products
Quality and baseline plots for spatial-temporal contiguous interferograms can be made using the *productPlot.py* program. 

### Time-series set-up of GUNW Products
Time-series set-up with spatial-temporal contiguous unwrapped interferograms and coherence can be done using the *TS_setup.py* program.


------
## Documentation

See the ARIA-tools-docs repo for all documentation:
+ [Jupyter Notebook Tutorials](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Notebooks.md)
+ [Overview of ARIA-tool modules](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Modules.md)


------
## Citation
D. Bekaert, M. Karim, L. Justin, H. Hua, P. Agram, S. Owen, G. Manipon, N. Malarout, M. Lucas, G. Sacco, L. Pan, S. Sangha, and ARIA team (2019), *Development and Dissemination of Standardized Geodetic Products by the Advanced Rapid Imaging and Analysis (ARIA) Center for Natural Hazards*, The International Union of Geodesy and Geophysics (IUGG), Montreal

------
## Contributors    

* David Bekaert
* Simran Sangha
* Emre Havazli
* Brett Buzzanga
* [_other community members_](https://github.com/dbekaert/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/dbekaert/ARIA-tools/blob/master/CONTRIBUTING.md).
