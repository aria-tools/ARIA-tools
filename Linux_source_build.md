## Linux
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **Linux** machines. For **Mac** see the respective installation instructions. 


### 1. Anaconda3
First install **python3** using either [Anaconda3](https://www.anaconda.com/distribution/) or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html). 

Below we use a clean installation of Miniconda3. First we will download Miniconda3:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Next execute the installer script and follow the instructions as provided by the installer.

The miniconda installation does not contain all the python module we need.
Use the **conda** excecutable to install the compiler tools that are needed for PROJ 4 installation and ARIA-tools.
```
conda install autoconf automake libtool patchelf numpy netcdf4 matplotlib pandas shapely --yes
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
```

You will need to manualy fix the linking of the PROJ libraries.
```
cd /my/proj/install/directory/lib
mv libproj.so.15.0.0 libinternalproj.so.15.0.0
ln -s libinternalproj.so.15.0.0 libinternalproj.so.15
ln -s libinternalproj.so.15.0.0 libinternalproj.so
rm -f libproj.*
patchelf --set-soname libinternalproj.so libinternalproj.so
```

### 3. GDAL SETUP
Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings:
```
cd gdal/gdal/
./configure --without-libtool --with-proj=/my/proj/install/directory --prefix=/my/gdal/install/directory --with-python
make -j4
make install
```


### 4. Setting of environment variables:
Edit your private module or start-up shell and add the PROJ and GDAL environment variables.

For example for csh do:
```
vi ~/.cshrc
```

Add the following and update ***my*** to the location where you installed the packages.
```
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:/my/proj/install/directory/lib:/my/gdal/install/directory/lib
setenv PROJ_LIB /my/proj/install/directory/share/proj
setenv PYTHONPATH $PYTHONPATH:/my/gdal/install/directory/lib/python3.7/site-packages
set path = ('/my/gdal/install/directory/bin' $path)
```
