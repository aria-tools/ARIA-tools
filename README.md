# ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/dbekaert/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. This software is open source under the terms of the GNU General Public License. Its development was funded under the NASA Sea-level Change Team (NSLCT) program and the Earth Surface and Interior (ESI) program. 

THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS.
USE AT YOUR OWN RISK.

For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [PySAR](https://github.com/yunjunz/PySAR). 


## Contents

1. [Download](#download)
2. [Software Dependencies](#software-dependencies)
   - [Installing dependencies with Anaconda](#dependencies-with-anaconda)
   - [Installing dependencies with Macports](#dependencies-with-macports)
3. [Installation](#installation)
4. [Running ARIA-tools](#running-aria-tools)


------

## 1. Download

To download the ARIA-tools repository:
```
git clone https://github.com/dbekaert/ARIA-tools.git
```

## 2. Software Dependencies
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

### [Dependencies with Anaconda](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md)
### [Dependencies with macports](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md)


## 3. Installation

The ARIA-tools package can be installed by adding the tools folder on your PATH and PYTHONPATH.
This can be done by editing your private module or your favorite start-up shell.


For example, for csh do:
```
vi ~/.cshrc
```

Add the following and update ***my*** to the location where you cloned ARIA-tools:
```
setenv PYTHONPATH $PYTHONPATH:/my/tools
set path = ('/my/tools' $path)
```


## 4. Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs/blob/master/README.md).

### 4.1. Manipulating GUNW Products
### 4.2. Baseline and quality control plots for GUNW Products
### 4.3 Time-series set-up of GUNW Products


## 5. Documentation

See the ARIA-tools-docs repo for all documentation:
+ [Jupyter Notebook Tutorials](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Notebooks.md)
+ [Overview of ARIA-tool modules](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Modules.md)

## 6. Citation
**TODO**


## 7. Contributors    

* David Bekaert
* Simran Sangha
* [_other community members_](https://github.com/dbekaert/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/dbekaert/ARIA-tools/blob/master/CONTRIBUTING.md).
