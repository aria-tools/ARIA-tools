# ARIA-tools
[![Language](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-green.svg)](https://github.com/aria-tools/ARIA-tools/blob/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate standard InSAR products from Sentinel-1 (ARIA GUNW-S1) and NISAR (NISAR_L2_GUNW). This software is open source under the terms of the [Apache 2.0 License](LICENSE). Its development was funded under the ROSES awards from the NASA Sea-level Change Team (NSLCT) program, the Earth Surface and Interior (ESI) program, and under the NISAR Science Team (NISAR-ST) program.


> [!IMPORTANT]
> ARIA-tools now supports **NISAR GUNW** products created routinely by the NASA-ISRO SAR mission. We welcome the community for reporting bugs using the issue ticket functionality

> [!NOTE]
> **Standardized Date and Phase Convention**
>
> To prevent confusion during time-series analysis, **all** ARIA-tools outputs enforce a consistent date and phase convention, regardless of the native format of the input product:
> * **Date Order:** ARIA-tools strictly uses a `ReferenceDate_SecondaryDate` naming convention, where the reference scene is the more recent pass and the secondary scene is the earlier pass. 
> * **Phase Sign:** Because the more recent acquisition is used as the reference, the unwrapped phase convention is standardized across all outputs. Negative phase differences indicate movement away from the sensor, and positive phase differences indicate movement towards the sensor.
> * **NISAR vs. S1:** While ARIA-S1-GUNW products natively follow this convention out of the box, NISAR GUNW products natively use the opposite date order. ARIA-tools automatically flips the dates and mathematical signs for NISAR inputs during extraction to ensure all your downstream time-series outputs are 100% consistent!

------

ARIA tools includes support for:
- ARIA Geocoded Unwrapped Interferogram from Sentinel-1 (GUNW-S1) can be downloaded for free from the [ASF DAAC vertex page](https://search.asf.alaska.edu/#/?dataset=SENTINEL-1%20INTERFEROGRAM%20(BETA)) selecting "ARIA S1 GUNW" under "datasets".  Users can also request free on-demand products through the [ASF on-demand system](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/). These products are added to the standard product archive. A log-on using the NASA Earthdata credentials is required to order new or download data from the archive. [Product Specification Document](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/#product-packaging).
-  NISAR Gecoded Unwrapped Interferogram products (NISAR_L2_GUNW) can be downloaded for free from the [ASF DAAC vertex page](https://search.asf.alaska.edu/#/?dataset=NISAR&prodConfig=PR&sciProducts=GUNW) selecting "NISAR" under "datasets" and selecting "GUNW" under the filter criteria. [Product Specification document](https://nisar-docs.asf.alaska.edu/gunw/). 


The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series. 

Actual time-series processing is not supported in ARIA-tools. However, outputs are compatible with third-party time-series InSAR packages such the "Miami INsar Time-series software in PYthon" ([MintPy](https://github.com/insarlab/MintPy)).
<p align="center">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/Hawaii.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/CA.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/EastCoast.png">
</p>

> [!CAUTION]
> THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents
1.  [Software Dependencies](#software-dependencies)
2.  [Installation](#installation)
    -   [Conda](#conda)
    -   [Installing a stable release from conda](#installing-a-stable-release-from-conda)
    -   [Installing the latest development branch](#installing-the-latest-development-branch)
    -   [Other installation options](#other-installation-options)
    -   [ARIA-tools with support for virtual data access](#aria-tools-with-support-for-virtual-data-access)
3.  [Running ARIA-tools](#running-aria-tools)
    -   [Commandline download of GUNW Products](#commandline-download-of-gunw-products)
    -   [Manipulating GUNW Products](#manipulating-gunw-products)
    -   [Baseline and quality control plots for GUNW Products](#baseline-and-quality-control-plots-for-gunw-products)
    -   [Time-series set-up of GUNW Products](#time-series-set-up-of-gunw-products)
    -   [Ordering ARIA S1 GUNW Products on demand](#ordering-aria-s1-gunw-products-on-demand)
4.  [Documentation](#documentation)
5.  [Citation](#citation)
6.  [Contributors and community contributions](#contributors)

------

## Software Dependencies
Below we list the key dependencies for ARIA-tools. See environment.yml for complete list.

### Packages
```
* Python >= 3.8  (3.9 preferred)
* [PROJ 4](https://github.com/OSGeo/proj) github) >= 6.0
* [GDAL](https://www.gdal.org/) and its Python bindings >= 3.7.0
```

### Python dependencies
```
* [SciPy](https://www.scipy.org/)
* [netcdf4](http://unidata.github.io/netcdf4-python/netCDF4/index.html)
* [requests](https://2.python-requests.org/en/master/)
* [asf_search](https://github.com/asfadmin/Discovery-asf_search) >=12.0.1
* [asf_search[asf-enumeration]](https://github.com/ASFHyP3/asf-enumeration) 
* [hyp3_sdk](https://github.com/ASFHyP3/hyp3-sdk)
```

### Optional Python Jupyter dependencies
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

ARIA-tools has been tested on the following system:
- Linux v.7 and up

Below we demonstrate how to build and setup an environment from scratch through Linux through a `TCSH` shell.

ARIA-tools package can be easily installed and used after the dependencies are installed and activated.

__[Conda](https://docs.conda.io/en/latest/index.html)__ is a cross-platform way to use Python that allows you to setup and use "virtual environments," which allows for the easy installation and management of all of the required dependencies. We recommend using the [Miniforge](https://github.com/conda-forge/miniforge) conda environment manager, which uses conda-forge as its default code repo. Alternatively, see __[here](https://docs.anaconda.com/anaconda/install/)__ for help installing Anaconda and __[here](https://docs.conda.io/en/latest/miniconda.html)__ for installing Miniconda.

### Conda
Below we outline the different steps for setting up the ARIA-tools while leveraging Miniforge for installation of the requirements.

Run the commands below to download and setup your Miniforge environment manager:

```.tcsh
cd ~/tools
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

# specify a new directory when prompted
bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p miniforge
miniforge/bin/mamba init tcsh

# reset shell
csh
```

### Installing a stable release from conda
To install a stable release of aria-tools from conda, refer to these instructions and disregard the rest of the installation section.
```.tcsh
# OPTIONAL: create/activate new env
# highly recommended to avoid potential conflicts with other packages
conda create --name ARIA-tools-conda
conda activate ARIA-tools-conda

# install aria tools
mamba install aria-tools
```


### Installing the latest development branch
To access and install the latest development branch from github, refer instead to these installation instructions.

Run the commands below to download/clone the ARIA-tools package to your local directory:

```.tcsh
cd ~/tools
git clone https://github.com/aria-tools/ARIA-tools.git
cd ARIA-tools
```

Run the commands below to install dependencies to a new conda environment `ARIA-tools` and activate it.
Make sure to activate your environment each time you open a new session:

```.tcsh
mamba env create -f environment.yml
conda activate ARIA-tools
```

We have included a `setup.py` script which allows for easy installation and setting up the ARIA-tools package itself (python and command line tools).
```.tcsh
python -m pip install -e .
```

If not using the setup.py, users should compile third-party packages manually and ensure ARIA-tools and dependencies are included on their PATH and PYTHONPATH. For `TCSH` shell this can be done as follows (replace `{$PWD}/tools/ARIAtools` to the location where you have cloned the ARIAtools repository):
```.tcsh
setenv PYTHONPATH ${PYTHONPATH}:{$PWD}/tools/ARIAtools
setenv PATH ${PATH}:${PWD}/tools/bin
```

To avoid potential issues associated with dependencies when cloning new ARIA-tools commits, it is advised to regularly maintain your conda environment as so (making sure to adjust the conda environment argument name `--name ARIA-tools` as appropriate):
```.tcsh
mamba env update --name ARIA-tools --file environment.yml --prune
```

GNU Parallel (https://www.gnu.org/software/parallel/) will write output to stdout that requests that the user cite their paper. We can use this command to suppress this output:
```
echo 'will cite' | parallel --citation
```

### Other installation options
The following pages might be of use to those trying to build third party packages from source.
-   [Installing dependencies from source on linux](https://github.com/aria-tools/ARIA-tools/blob/master/LinuxSourceBuild.md)
-   [Installing dependencies from source on mac](https://github.com/aria-tools/ARIA-tools/blob/master/MacOSSourceBuild.md)


### ARIA-tools with support for virtual data access
GDAL Virtual File Systems capabilities (vsicurl) can be leveraged in ARIA-tools to avoid downloading GUNW products during processing. This is supported for both **Sentinel-1 GUNW** and **NISAR GUNW** products.

To use virtual access, generate a URL list using `ariaDownload.py -o url`, which produces a `.txt` file of product URLs. Pass this `.txt` file directly to `ariaExtract.py`, `ariaTSsetup.py`, or `ariaPlot.py` via the `-f` option — no download step is needed. ARIA-tools will stream data on-the-fly from the ASF archive using GDAL's `/vsicurl/` driver.

For virtual processing, a local metadata cache (`aria_meta_cache.json`) is automatically created on first run to avoid repeated remote reads of product headers. Subsequent runs read from this cache for faster initialization. When expanding a time series with additional products, the existing cache is reused — only new URLs trigger remote metadata reads.

Minimum requirements:
```
* [GDAL](https://www.gdal.org/) and its Python bindings >= 3.8
* Linux kernel >=4.3 
* libnetcdf >=4.5 
```

A '~/.netrc' file with earthdata credential included
```
echo "machine urs.earthdata.nasa.gov login myUsername password myPassword" > ~/.netrc
chmod 600 ~/.netrc
```

------
## Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/aria-tools/ARIA-tools-docs). We welcome the community to contribute other examples on how to leverage the ARIA-tools (see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md) for instructions).


### Commandline download of GUNW Products
ARIA GUNW-S1/NISAR_L2_GUNW products can be downloaded through the commandline using the *ariaDownload.py* program, which wraps around the ASF DAAC api. There is the option for virtual data access by using `ariaDownload.py -o url`, which creates a `.txt` file with HTTPS paths to archived products instead of downloading them. This URL file can be passed directly to subsequent programs (*ariaExtract.py*, *ariaTSsetup.py*, *ariaPlot.py*) for streaming access without local downloads (see [ARIA-tools with support for virtual data access](#aria-tools-with-support-for-virtual-data-access)).

### Manipulating GUNW Products
ARIA GUNW-S1/NISAR_L2_GUNW products can be manipulated (cropped, stitched, extracted) using the *ariaExtract.py* program. 

### Baseline and quality control plots for GUNW Products
ARIA GUNW-S1/ NISAR_L2_GUNW quality and baseline plots for spatial-temporal contiguous interferograms can be made using the *ariaPlot.py* program.

### Time-series set-up of GUNW Products
ARIA GUNW-S1/NISAR_L2_GUNW time-series set-up with spatial-temporal contiguous unwrapped interferograms and coherence can be done using the *ariaTSsetup.py* program.

### Ordering ARIA S1 GUNW Products on demand
The *ariaOrderASF.py* program supports on-demand ordering of **ARIA Sentinel-1 GUNW** products through the [ASF HyP3 on-demand processing system](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/). This tool is for GUNW-S1 products only (not NISAR) and allows users to build and order additional interferometric pairs that are not yet in the ASF archive, using each user's monthly HyP3 credit quota.  A `~/.netrc` file with NASA Earthdata credentials is required for authentication.

> [!NOTE]  
> We support extraction of correction layers (e.g. Troposphere, Ionosphere, Solid Earth Tides) as well as geometry information (e.g. incidence angle, look angle, baselines, etc) embeded within the GUNW products 

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
-   Alexander Fore
-   Marin Govorcin
-   Charles Marshak
-   Joseph Kennedy
-   [other community members](https://github.com/aria-tools/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md).

