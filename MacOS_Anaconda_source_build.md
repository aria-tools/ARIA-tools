# MacOSX
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines using **Anaconda**.  For a variant to this using Macports click [here](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md).


For **Linux** installation instructions click [here](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md).



------
## Contents

0. [XCode and Command Line Tools](#xcode-and-command-line-tools)
1. [Install with Anaconda](#install-with-anaconda)
2. [Building from source](#building-from-source)
    1. [Anaconda3 Setup](#anaconda3-setup)
    2. [PROJ4 Setup](#proj4-setup)
    3. [GDAL Setup](#gdal-setup)
3. [Jupyter Notebooks SETUP](#Jupyter-Notebooks-Setup)
4. [Setting of environment variables](#setting-of-environment-variables)
5. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
------
## 0. XCode and Command Line Tools
Make sure that you have XCode and Command Line Tools installed. XCode 10.X versions do not install the header files under /usr/include but we need that directory for the installation.
For that purpose, after installing XCode and Command Line Tools run:
```
open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
```
this command will install header files under /usr/include

------

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
## 2. Building from source
You can build gdal and proj4 from source. Other required packages can be installed by using conda package manager.

### i. Anaconda3 Setup
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

### ii. PROJ4 Setup
Clone the **PROJ 4** repository from github and install at least the version 6 release (i.e. main branch).
```
git clone https://github.com/OSGeo/proj.4 proj
```

Build the PROJ package.
```
cd proj
./autogen.sh
setenv CXXFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
setenv CFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
./configure --prefix=/my/proj/install/directory --disable-static
make -j4
make install

cd /my/proj/install/directory/lib
rm libproj.dylib
mv libproj.15.dylib libinternalproj.15.dylib
ln -s libinternalproj.15.dylib libinternalproj.dylib

```

### iii. GDAL Setup
Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with
```
cd gdal/gdal/
./autogen.sh
./configure --with-proj=/my/proj/install/directory --prefix=/my/gdal/install/directory --with-sqlite3  --with-libjson-c --with-geos --with-netcdf
make -j4

cd apps/
make test_ogrsf

cd ../swig/python/
setenv CPATH /Library/Developer/CommandLineTools/usr/include/c++/v1
python setup.py build

cd ../../
make install
```

------
## 3.  Jupyter Notebooks Setup
Instructions to install jupyter notebooks in a conda environment

```
conda install -c conda-forge jupyterlab --yes
```

Conda will install all required jupyter packages.
------
## 4. Setting of environment variables
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
## 5. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
