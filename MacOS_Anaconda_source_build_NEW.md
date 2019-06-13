# MacOSX
Here we provide guidelines on how to build GDAL 3.0+ and PROJ 4 (6.X.X) for **MacOS** machines using **Anaconda**.  For a variant to this using Macports click [here](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md).

For **Linux** installation instructions click [here](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md).

------
## Contents

0. [Anaconda3 Setup](#anaconda3-setup)
1. [Install with Anaconda](#install-with-anaconda)
2. [Jupyter Notebooks SETUP](#jupyter-notebooks-setup)
3. [Setting of environment variables](#setting-of-environment-variables)
4. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
------
## 0. Anaconda3 Setup
First install **python3** using either [Anaconda3](https://www.anaconda.com/distribution/) or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).

Below we use a clean installation of Miniconda3. First we will download Miniconda3:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
Next execute the installer script and follow the instructions as provided by the installer.

The miniconda installation does not contain all the python module we need.
Use the **conda** excecutable to install the compiler tools that are needed for PROJ 4 installation and ARIA-tools.
```
conda install autoconf automake libtool numpy netcdf4 matplotlib pandas sqlite pkg-config shapely postgresql libcxx lapack --yes
```

## 1. Install with Anaconda
ARIA-tools require GDAL3 environtment and GDAL3 can be installed using Anaconda package manager.
The lines below will create and environment called ARIA-tools and install required packages

```
conda config --add channels conda-forge
conda create -n ARIA-tools python=3.7 libtool numpy netcdf4 matplotlib pandas sqlite pkg-config shapely postgresql libcxx lapack
source activate ARIA-tools
conda install gdal
```

------
## 2.  Jupyter Notebooks Setup
Instructions to install jupyter notebooks in a conda environment

```
conda install -c conda-forge jupyterlab --yes
```

Conda will install all required jupyter packages.
------
## 3. Setting of environment variables
Edit your private module or start-up shell and add the PROJ and GDAL environment variables.

For example for csh do:
```
vi ~/.cshrc
```

Add the following and update ***my*** to the location where you installed the packages.
```
setenv LD_LIBRARY_PATH /my/python/directory/lib:/my/gdal/install/directory/lib
setenv GDAL_DATA /my/gdal/install/directory/share/gdal
setenv PYTHONPATH /my/gdal/install/directory/lib/python3.7/site-packages
set path = ('/my/gdal/install/directory/bin' $path)
```

------
## 4. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
