# MacOSx
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines using **Macports**. For a variant to this using Anaconda click [here](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_Anaconda_source_build.md). 


For **Linux** installation instructions see [here](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md). 

------
## Contents

1. [MacPorts](#macports)
2. [PROJ 4 SETUP](#proj-4-setup) 
3. [GDAL SETUP](#gdal-setup)
4. [Setting of environment variables](#setting-of-environment-variables)
5. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)


------
## MacPorts
First install **Macports** following the instructions at https://www.macports.org/.

Use the **ports** package manager to install python 3.X.X and associated packages, and compiler tools that are needed for PROJ 4 installation and ARIA-tools. 
The following instructions are tested using python 3.6.8. Note that below we default the python executable to python3.

```
sudo port -N install autoconf automake libtool json-c
sudo port -N install python3X py3X-scipy py3X-matplotlib py3X-pandas py3X-shapely py3X-netcdf4 py3X-pkgconfig
sudo port select python python3X
```
Place macports install directories (/opt/local/...) ahead of system directories (/usr/...) so that compilers find the right libraries.
```
export PATH="/opt/local/bin:/opt/local/lib:/opt/local/sbin:$(getconf PATH)"
```

------
## PROJ 4 SETUP
```
cd /my/proj
```

Clone the **PROJ 4** repository from github and install at least the version 6 release (i.e. main branch).

```
git clone https://github.com/OSGeo/proj.4 proj
```

Build the PROJ package.
```
mkdir install
cd proj
./autogen.sh
./configure --prefix=/my/proj/install 
make -j4
make install
```

------
## GDAL SETUP
```
cd /my/gdal
mkdir -p install
```

Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings and add them to the location of the python library. Note python bindings are built for the **python** executable. If you are having both python2 and 3 installed make sure you are installing it for python3. We have experienced conflicts it a pre-existing GDAL version was already installed with Macports. 
```
cd /my/gdal/gdal/
./configure --with-proj=/my/proj/install --prefix=/my/gdal/install 
make -j4 
make install
```
if configure fails, the following options may help:
```
./configure --without-libtool --with-proj=/my/proj/install --prefix=/my/gdal/install \
            --with-libjson-c=internal --with-local=/opt/local --with-python=yes \
            --with-kea=no
```

------
## Setting of environment variables
Edit your private module or start-up shell and add the PROJ and GDAL environment variables.

For example for bash do:
```
vi ~/.bash_profile
```

Add the following and update ***my*** to the location where you installed the packages.
```
export PROJ_LIB=/my/proj/install/share/proj
export PATH="/my/gdal/install/bin:/opt/local/bin:/opt/local/lib:/opt/local/sbin:$(getconf PATH)"
```

------
## [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
