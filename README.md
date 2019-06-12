# ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/dbekaert/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. This software is open source under the terms of the GNU General Public License. Its development was funded under the NASA Sea-level Change Team (NSLCT) program and the Earth Surface and Interior (ESI) program.


For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [PySAR](https://github.com/yunjunz/PySAR).
<p align="center">
  <img height="250" src="https://github.com/dbekaert/ARIA-tools-docs/blob/master/images/Hawaii.png">
  <img height="250" src="https://github.com/dbekaert/ARIA-tools-docs/blob/master/images/CA.png">
  <img height="250" src="https://github.com/dbekaert/ARIA-tools-docs/blob/master/images/EastCoast.png">
</p>
THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents

1. [Software Dependencies](#software-dependencies)
2. [Installation](#installation)
3. [Running ARIA-tools](#running-aria-tools)
### Commandline downloading of GUNW Products
   - [Commandline download of GUNW Products](#commandline-download-of-gunw-products)
   - [Manipulating GUNW Products](#manipulating-gunw-products)
   - [Baseline and quality control plots for GUNW Products](#baseline-and-quality-control-plots-for-gunw-products)
   - [Time-series set-up of GUNW Products](#time-series-set-up-of-gunw-products)
4. [Documentation](#documentation)
5. [Citation](#citation)
6. [Contributors and community contributions](#contributors)


------

## Software Dependencies
Below we list the dependencies for ARIA-tools

### Packages:
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

------
## Installation
ARIA-tools package can be easily installed and used after the dependencies are installed and activated.
We strongly recommend using [Anaconda](https://www.anaconda.com/distribution/) package manager for the installation of dependencies in python environment.

Below we outline the differnt steps for setting up the ARIA-tools while leveraging Anaconda for installation of the requirements. Running the commands below will clone the ARIA-tools package to your local directory, create a conda environment with the name 'ARIA-tools', install dependencies to this environment and activate it.

```
git clone https://github.com/dbekaert/ARIA-tools.git
conda config --add channels conda-forge
conda create -n ARIA-tools --file ./ARIA-tools/requirements.txt --yes
conda activate ARIA-tools
```

After the installation of ARIA-tools package and dependencies we need to update our PATH and PYTHONPATH variables and add a new PROJ_LIB variable to our shell environment.
This can be done by editing your private module or your favorite start-up shell profile.


For example, for csh do:
```
vi ~/.cshrc
```

Add the following and update *my* to the location where you cloned ARIA-tools:
```
setenv PYTHONPATH $PYTHONPATH:/my/tools/python
setenv PROJ_LIB /my/python/directory/share/proj
set path = ('/my/tools/python' $path)
```


### Other installation options
The following pages might be of use to those trying to build thrid party packages from source.
- [Installing dependencies on linux with Anaconda](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md)
- [Installing dependencies on mac with macports](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md)
- [Installing dependencies on mac with Anaconda](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_Anaconda_source_build.md) 	

------
## Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs/blob/master/README.md). We welcome the community to contribute other examples on how to leverage the ARIA-tools (see [here](https://github.com/dbekaert/ARIA-tools/blob/master/CONTRIBUTING.md) for instructions).

### Commandline download of GUNW Products
GUNW products can be downloaded through the commandline using the *productAPI.py* program, which wraps around the ASF DAAC api.

### Manipulating GUNW Products
GUNW product can be manipulated (cropped, stitched, extracted) using the *extractProduct.py* program.

### Baseline and quality control plots for GUNW Products
Quality and baseline plots for spatial-temporal contiguous interferograms can be made using the *productPlot.py* program.

### Time-series set-up of GUNW Products
Time-series set-up with spatial-temporal contiguous unwrapped interferograms and coherence can be done using the *TS_setup.py* program.


------
## Documentation

See the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs) for all documentation and Jupyter Notebook Tutorials.

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
