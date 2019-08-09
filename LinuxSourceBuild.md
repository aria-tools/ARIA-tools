# Linux
Here we provide guidelines on how to build GDAL 3.0+ and PROJ 4 (6.X.X) from source for **Linux** machines.

For **Mac** see the respective installation instructions provided on the [ARIA-tools](https://github.com/aria-tools/ARIA-tools) page.

------
## Contents

1.  [Anaconda3](#anaconda3)
2.  [PROJ 4 SETUP](#proj-4-setup)
3.  [GDAL SETUP](#gdal-setup)
4.  [Jupyter Notebooks SETUP](#jupyter-notebooks-setup)
5.  [Setting of environment variables](#setting-of-environment-variables)
6.  [Return to back to ARIA-tools page](https://github.com/aria-tools/ARIA-tools)

------
## Anaconda3
First install **python3** using either [Anaconda3](https://www.anaconda.com/distribution/) or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).

Below we use a clean installation of Miniconda3. First we will download Miniconda3:
```.tcsh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Next execute the installer script and follow the instructions as provided by the installer.

The miniconda installation does not contain all the python module we need.
Use the **conda** excecutable to install the compiler tools that are needed for PROJ 4 installation and ARIA-tools.
We can add the conda-forge channel to our installation so the packages in conda-forge are indexed.

```.tcsh
conda config --add channels conda-forge
```
After adding conda-forge we can make the installation of required packages

```.tcsh
conda install autoconf automake libtool numpy netcdf4 matplotlib pandas sqlite pkg-config shapely postgresql libcxx lapack --yes
```
As soon as conda is done with the installation you can run update command below and finish conda installation
```.tcsh
conda update --all
```

------
## PROJ 4 SETUP
Clone the **PROJ 4** repository from github and install at least the version 6 release (i.e. main branch).
```.tcsh
git clone https://github.com/OSGeo/proj.4 proj
```

Build the PROJ package.
```.tcsh
cd proj
./autogen.sh
setenv CXXFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
setenv CFLAGS "-DPROJ_RENAME_SYMBOLS -O2"
./configure --prefix=/my/proj/install/directory --disable-static
make -j4
make install
```

------
## GDAL SETUP
Clone the GDAL repository from github with a version of at least 3.0 (i.e. main branch).
```.tcsh
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings:
```.tcsh
cd gdal/gdal/
./configure --without-libtool --with-proj=/my/proj/install/directory --prefix=/my/gdal/install/directory --with-python
make -j4
make install
```

------
## Jupyter Notebooks Setup
Instructions to install jupyter notebooks in a conda environment

```.tcsh
conda install -c conda-forge jupyterlab jupyter_contrib_nbextensions rise --yes
```

Conda will install all required jupyter packages plus contributed notebook extensions, configurator, hide code and RISE to turn notebooks into slideshow.

------
## Setting of environment variables
Edit your private module or start-up shell and add the PROJ and GDAL environment variables.

For example for csh do:
```.tcsh
vi ~/.cshrc
```

Add the following and update ***my*** to the location where you installed the packages.
```.tcsh
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH:/my/proj/install/directory/lib:/my/gdal/install/directory/lib
setenv PROJ_LIB /my/proj/install/directory/share/proj
setenv PYTHONPATH $PYTHONPATH:/my/gdal/install/directory/lib/python3.7/site-packages
set path = ('/my/gdal/install/directory/bin' $path)
```

------
## [Return to back to ARIA-tools page](https://github.com/aria-tools/ARIA-tools)
