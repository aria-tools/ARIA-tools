## ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/dbekaert/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [PySAR](https://github.com/yunjunz/PySAR). 


### 1. Download

To download the ARIA-tools repository:
```
git clone https://github.com/dbekaert/ARIA-tools.git
```

### 2. Installation
The ARIA-tools package can be installed by adding the tools folder on your PATH and PYTHONPATH.
This can be done by editing your private module or your favorite start-up shell.


For example, for csh do:
```
vi ~/.cshrc
```

Add the following and update ***my*** to the location where you cloned ARIA-tools:
```
setenv PYTHONPATH /my/tools
set path = ('/my/tools' $path)
```

**Dependencies:**
The following dependencies needs to be installed. Note that GDAL currently needs to be build manually from source as no pre-build is available via anaconda or macports for version 2.5. In addition PROJ 4 version 6.x is required for GDAL 2.5. Below we provide guidelines for Linux users on how to build these dependencies from source. In future we will expand this section to rely on anaconda and macports for installation.
- GDAL v2.5 (currently only available on [GDAL](https://www.gdal.org/) github)
- PROJ 6 (see [PROJ4](https://proj4.org/index.html) webpage and [PROJ4](https://github.com/OSGeo/proj) github)
- [SciPy](https://www.scipy.org/) 
- [netcdf4](http://unidata.github.io/netcdf4-python/netCDF4/index.html)
- [HDF5](https://www.h5py.org/) 


[**Linux**](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md)


### 3. Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs/blob/master/README.md).

#### 3.1. Manipulating GUNW Products
#### 3.2. Baseline and quality control plots for GUNW Products
#### 3.3 Time-series set-up of GUNW Products


### 4. Documentation

See the ARIA-tools-docs repo for all documentation:
+ [Jupyter Notebook Tutorials](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Notebooks.md)
+ [Overview of ARIA-tool modules](https://github.com/dbekaert/ARIA-tools-docs/tree/master/Modules.md)

### 5. Citation
**TODO**


### 6. Contributors    

* David Bekaert
* Simran Sangha
* [_other community members_](https://github.com/dbekaert/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/dbekaert/ARIA-tools/blob/master/CONTRIBUTING.md).
