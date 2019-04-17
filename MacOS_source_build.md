## MacOSx
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines. For **Linux** see the respective installation instructions. 

** Note this is untested on a clean install as of 4/17/2019, but will be conducted shortly. **

### 1. MacPorts
First install **Macports** following the instructions at https://www.macports.org/.

Use the **ports** package manager to install python 3.X.X and associated packages, and compiler tools that are needed for PROJ 4 installation and ARIA-tools.
```
sudo port install py3X-scipy py3X-matplotlib py3X-pandas py3X-shapely
sudo port install autoconf automake libtool numpy hdf5 netcdf4 
```
Place macports install directories (/opt/local/...) ahead of system directories (/usr/...) so that compilers find the right libraries.
```
export PATH="/opt/local/bin:/opt/local/lib:/opt/local/sbin:$(getconf PATH)"
```

### 2. PROJ 4 SETUP
Clone the **PROJ 4** repository from github and install at least the version 6 release (i.e. main branch).

```
git clone https://github.com/OSGeo/proj.4 /my/proj
```

Build the PROJ package.
```
cd /my/proj
mkdir install
./autogen.sh
./configure --prefix=./install --disable-static
make -j4
make install
```

### 3. GDAL SETUP
```
cd /my/gdal; mkdir -p install
```

Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings:
```
cd /my/gdal/gdal/
./configure --without-libtool --with-proj=/my/proj/install --prefix=/my/gdal/install \
            --with-libjson-c=internal --with-local=/opt/local --with-python=yes \
            --with-kea=no

make -j4
make install
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
