# MacOSX
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines using **Anaconda**.  For a variant to this using Macports click [here](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_source_build.md).


For **Linux** installation instructions click [here](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md).



------
## Contents

1. [XCode and Command Line Tools](#xcode-and-command-line-tools)
2. [Anaconda3 and PROJ 4 SETUP](#anaconda3-and-proj-4-setup) 
3. [GDAL SETUP](#gdal-setup)
4. [Setting of environment variables](#setting-of-environment-variables)
5. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)


------
## XCode and Command Line Tools
Make sure that you have XCode and Command Line Tools installed. XCode 10.X versions do not install the header files under /usr/include but we need that directory for the installation.
For that purpose, after installing XCode and Command Line Tools run:
```
open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
```
this command will install some header files under /usr/include


------
## Anaconda3 and PROJ 4 setup
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

### 2. PROJ 4 SETUP
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

------
## GDAL SETUP
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
## Setting of environment variables
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


## [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)

