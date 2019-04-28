# MacOSx
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines using **Macports**. For a variant to this using Anaconda click [here](https://github.com/dbekaert/ARIA-tools/blob/master/MacOS_Anaconda_source_build.md). 


For **Linux** installation instructions see [here](https://github.com/dbekaert/ARIA-tools/blob/master/Linux_source_build.md). 

------
## Contents
Make sure to run these in the specified order!
1. [MacPorts](#macports)
2. [GDAL COMMANDLINE SETUP](#gdal-commandline-setup)
3. [GDAL PYTHON bindings](#gdal-python-bindings)
4. [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)


------
## MacPorts
First install **Macports** following the instructions at https://www.macports.org/.

Use the **ports** package manager to install python 3.X.X and associated packages, and compiler tools that are needed for PROJ 4 installation and ARIA-tools. 
The following instructions are tested using python 3.6.8. Note that below we default the python executable to python3.

```
sudo port -N install autoconf automake libtool json-c netcdf proj6
sudo port -N install python3X py3X-scipy py3X-matplotlib py3X-pandas py3X-shapely py3X-netcdf4 py3X-pkgconfig
sudo port select python python3X
```
Place macports install directories (/opt/local/...) ahead of system directories (/usr/...) so that compilers find the right libraries. e.g. for *bash*:
```
export PATH="/opt/local/bin:/opt/local/lib:/opt/local/sbin:$(getconf PATH)"
```

------
## GDAL COMMANDLINE SETUP
```
cd /my/gdal
mkdir install
```

Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```


We will first build the GDAL package and command line tools. Note it is better not to build the python tools at the same time given library and compiler conflicts with NumPy. 

```
cd /my/gdal/gdal/
./configure --prefix=/my/gdal/install --with-proj=/opt/local/lib/proj6 --with-sqlite3  --with-libjson-c=internal
```
if configure fails, you might want to try and specify additional options, such as:
```
--with-kea=no
```
Next make and install GDAL:
```
make -j4 
make install
```

Next you can inspect if the following passes without errors:
```
cd apps
make test_ogrsf
cd ..
```


Before proceeding to the python bindings installation make sure your *LD_LIBRARY_PATH*, *DYLD_LIBRARY_PATH* and *PATH* in your shell to point to the new GDAL installation lib and bin folders. E.g for *bash*:
```
export LD_LIBRARY_PATH="/my/gdal/install/lib"
export DYLD_LIBRARY_PATH="/my/gdal/install/lib"
export PATH="/my/gdal/install/bin:$(getconf PATH)"
```


------
## GDAL PYTHON bindings

Make sure your *LD_LIBRARY_PATH*, *DYLD_LIBRARY_PATH* and *PATH*  are pointing to your GDAL installation.
```
echo $LD_LIBRARY_PATH
echo $DYLD_LIBRARY_PATH
echo $PATH
```

To install the python bindings:
```
cd swig/python
```
Edit the **setup.cfg** file such that *gdal-config* is pointing to your install folder of GDAL.
```
vim setup.cfg
#update gdal_config=/my/gdal/install/bin/gdal-config
```


The GDAL python bindings are built here for the **python** executable. If you are having both python2 and 3 installed make sure you are installing it for python3. Here we assume python is an alias for python3 and specify the use of **gcc** and **g++**. Note the later should be the same as used for building **NumPy**. If the test provided below fail you could remove **gcc** and **g++** specification and try again.
```
CC=gcc CXX=g++ python setup.py build >& blog
```
Note if you re-run the build command make sure to delete the build folder.


Next, we will extract the library build paths:
```
grep lgdal blog > linkcmds
```
Edit this file and ***remove all*** '-L/opt/local/lib' occurrences and 'blog:' occurrences at beginning of line.

Once done you can execute the commands using:
```
chmod +x linkcmds
./linkcmds
```
To set the location the installation for install add the python bindings to your *PYTHONPATH* and source the environment. E.g for *bash*:
```
export PYTHONPATH="/my/gdal/install/lib/python:${PYTHONPATH}"
mkdir /my/gdal/install/lib/python
```

Now we will install them at the same location as out GDAL command line tools using:
```
python setup.py install --home=/my/gdal/install/
```

To test your GDAL python bindings are correctly working with **NumPy** try:
```
python
> from osgeo import gdal_array
> import gdalnumeric
```

------
## [Return to back to ARIA-tools page](https://github.com/dbekaert/ARIA-tools)
