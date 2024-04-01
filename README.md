# ARIA-tools

[![Language](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-green.svg)](https://github.com/aria-tools/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. This software is open source under the terms of the [Apache 2.0 License](LICENSE). Its development was funded under the NASA Sea-level Change Team (NSLCT) program and the Earth Surface and Interior (ESI) program.

For a full overview of available ARIA standard products and their specification, see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be downloaded for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series. 

Actual time-series processing is not supported in ARIA-tools. However, outputs are compatible with third-party time-series InSAR packages such as the "Generic InSAR Analysis Toolbox" ([GIAnT](http://earthdef.caltech.edu/projects/giant/wiki)) and the "Miami INsar Time-series software in PYthon" ([MintPy](https://github.com/insarlab/MintPy)).
<p align="center">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/Hawaii.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/CA.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/EastCoast.png">
</p>
THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents
1.  [Software Dependencies](#software-dependencies)
2.  [Installation](#installation)
    -   [Conda](#conda)
    -   [Other installation options](#other-installation-options)
    -   [ARIA-tools with support for S3 virtual data access](#aria-tools-with-support-for-s3-virtual-data-access)
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
ARIA-tools package can be easily installed and used after the dependencies are installed and activated. The third-party RelaxIV package is optional (not required), and  only used when opting to minimizing phase-discontinuities. Prior to use of RelaxIV, users should conform to the RelaxIV license agreement. The easiest way of installing RelaxIV is by downloading the min-cost-flow repository in the third-party folder of the ARIAtools and using the setup.py script as outlined below.

__[Conda](https://docs.conda.io/en/latest/index.html)__ is a cross-platform way to use Python that allows you to setup and use "virtual environments," which allows for the easy installation and management of all of the required dependencies. We recommend using the [Miniforge](https://github.com/conda-forge/miniforge) conda environment manager, which uses conda-forge as its default code repo. Alternatively, see __[here](https://docs.anaconda.com/anaconda/install/)__ for help installing Anaconda and __[here](https://docs.conda.io/en/latest/miniconda.html)__ for installing Miniconda.

### Conda
Below we outline the different steps for setting up the ARIA-tools while leveraging Anaconda for installation of the requirements. Run the commands below to download/clone the ARIA-tools package to your local directory.:

```.tcsh
conda config --add channels conda-forge
conda install mamba
git clone https://github.com/aria-tools/ARIA-tools.git
cd ARIA-tools
```

Run the commands below to install dependencies to a new conda environment `ARIA-tools` and activate it:

```.tcsh
mamba env create -f environment.yml
conda activate ARIA-tools
```

Or run the commands below to install dependencies to an existing conda environment (`base` by default):

```.tcsh
mamba install -c conda-forge --yes --file requirements.txt
```

We have included a `setup.py` script which allows for easy compilation and installation of third-party dependencies (c-code), as well as for setting up the ARIA-tools package itself (python and command line tools).
```.tcsh
python -m pip install -e .
```

If not using the setup.py, users should compile third-party packages manually and ensure ARIA-tools and dependencies are included on their PATH and PYTHONPATH. For c-shell this can be done as follows (replace "ARIAtoolsREPO" to the location where you have cloned the ARIAtools repository):
```.tcsh
setenv PYTHONPATH ${PYTHONPATH}:{$PWD}/tools/ARIAtools
setenv PATH ${PATH}:${PWD}/tools/ARIAtools
```

To avoid potential issues associated with dependencies when cloning new ARIA-tools commits, it is advised to regularly maintain your conda environment as so (making sure to adjust the conda environment argument name `--name ARIA-tools` as appropriate):
```.tcsh
conda env update --name ARIA-tools --file environment.yml --prune
```

### Other installation options
The following pages might be of use to those trying to build third party packages from source.
-   [Installing dependencies from source on linux](https://github.com/aria-tools/ARIA-tools/blob/master/LinuxSourceBuild.md)
-   [Installing dependencies from source on mac](https://github.com/aria-tools/ARIA-tools/blob/master/MacOSSourceBuild.md)


### ARIA-tools with support for S3 virtual data access
GDAL Virtual File Systems capabilities (vsicurl) can be leveraged in ARIA-tools to avoid download of product during processing. 

Minimum requirements:
```
* [GDAL](https://www.gdal.org/) and its Python bindings >= 3.0
* Linux kernel >=4.3 
* libnetcdf >=4.5 
```

A '~/.netrc' file with earthdata credential included
```
echo "machine urs.earthdata.nasa.gov login myUsername password myPassword" > ~/.netrc
chmod 600 ~/.netrc
```

In addition, users should set the following environment variables:
```.bash
export GDAL_HTTP_COOKIEFILE=/tmp/cookies.txt
export GDAL_HTTP_COOKIEJAR=/tmp/cookies.txt
export VSI_CACHE=YES
```

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
Buzzanga, B., Bekaert, D. P. S., Hamlington, B. D., & Sangha, S. S. (2020). Towards Sustained Monitoring of Subsidence at the Coast using InSAR and GPS: An Application in Hampton Roads, Virginia. Geophysical Research Letters, 47, e2020GL090013. [https://doi.org/10.1029/2020GL090013](https://doi.org/10.1029/2020GL090013)

------
## Contributors
-   David Bekaert
-   Simran Sangha
-   Emre Havazli
-   Brett Buzzanga
-   [other community members](https://github.com/aria-tools/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md).
