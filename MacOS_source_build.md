## MacOSx
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines. For **Linux** see the respective installation instructions. 

### 1. MacPorts
First install **Macports** following the instructions at https://www.macports.org/.

Use the **ports** package manager to install python 3.X.X and associated packages, and compiler tools that are needed for PROJ 4 installation and ARIA-tools. 
The following instructions are tested using python 3.6.8. Note that below we default the python executable to python3.

```
sudo port -N install autoconf automake libtool json-c
sudo port -N install python3X py3X-scipy py3X-matplotlib py3X-pandas py3X-shapely
sudo port select python python3X
sudo port -N install netcdf hdf5 pkgconfig
```
Place macports install directories (/opt/local/...) ahead of system directories (/usr/...) so that compilers find the right libraries.
```
export PATH="/opt/local/bin:/opt/local/lib:/opt/local/sbin:$(getconf PATH)"
```

### 2. PROJ 4 SETUP
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

### 3. GDAL SETUP
```
cd /my/gdal
mkdir -p install
```

Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings. Note python bindings are build for the **python** executable. If you are having both python2 and 3 installed make sure you are installing it for python3.
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

### 4. Setting of environment variables:
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
