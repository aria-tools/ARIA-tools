## LINUX
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 6 from source for **LINUX** machines. For MAC see the respective installation instructions. 


### 1. Anaconda3
First install **python3** using either [Anaconda3](https://www.anaconda.com/distribution/) or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html). Below we use a clean installation of Miniconda3.

First download Miniconda3 and then execute the installer script. Follow the instructions as provided by the installer.
'''
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
'''

Use the **conda** excecutable to install the compiler tools.
'''
conda install autoconf automake libtool patchelf
'''

### 2. PROJ 6 SETUP
Clone the PROJ4 repository from github or use at least version 6 release (i.e. main branch).
'''
git clone https://github.com/OSGeo/proj.4 proj
'''

Build the PROJ package.
'''
cd proj
./autogen.sh
setenv CXXFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
setenv CFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
./configure --prefix=/my/proj/install/directory --disable-static
make -j4
make install
'''

You will need to manualy fix the linking of the PROJ libraries.
'''
cd /my/proj/install/directory/lib
mv libproj.so.15.0.0 libinternalproj.so.15.0.0
ln -s libinternalproj.so.15.0.0 libinternalproj.so.15
ln -s libinternalproj.so.15.0.0 libinternalproj.so
rm -f libproj.*
patchelf --set-soname libinternalproj.so libinternalproj.so
'''

### 3. GDAL SETUP
Clone the GDAL respository from github with a version of at least 2.5 (i.e. main branch).
'''
git clone https://github.com/OSGeo/gdal
'''

Build the GDAL package:
'''
cd gdal/gdal/
./configure --without-libtool --with-proj=/my/proj/install/directory --prefix=/my/gdal/install/directory --with-python
make -j4
make install
'''


### 4. Setting of environment variables:
Edit your private module or start-up shell and add the PROJ and GDAL environment variables.

For example for csh do:
'''
vi ~/.cshrc
'''

Add the following and update **my** to the location where you installed the packages.
'''
setenv LD_LIBRARY_PATH /my/proj/install/directory/lib:/my/gdal/install/directory/lib
setenv PROJ_LIB /my/proj/install/directory/share/proj
setenv PYTHONPATH /my/gdal/install/directory/lib/python3.7/site-packages
set path = ('/my/gdal/install/directory/bin' $path)
'''
