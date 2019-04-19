## MacOSX
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines. For **Linux** see the respective installation instructions.

### 0. XCode and Command Line Tools
Make sure that you have XCode and Command Line Tools installed.

### 1. Anaconda3
First install **python3** using either [Anaconda3](https://www.anaconda.com/distribution/) or [Miniconda3](https://docs.conda.io/en/latest/miniconda.html).

Below we use a clean installation of Miniconda3. First we will download Miniconda3:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
```
Next execute the installer script and follow the instructions as provided by the installer.

The miniconda installation does not contain all the python module we need.
Use the **conda** excecutable to install the compiler tools that are needed for PROJ 4 installation and ARIA-tools.
```
conda install autoconf automake libtool numpy netcdf4 matplotlib pandas sqlite pkg-config shapely --yes
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

### 3. GDAL SETUP
Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Build the GDAL package with the python bindings:
```
cd gdal/gdal/
./configure --with-proj=/my/proj/install/directory --prefix=/my/gdal/install/directory --with-python
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
setenv PROJ_LIB /my/proj/install/directory/share/proj
setenv PYTHONPATH $PYTHONPATH:/my/gdal/install/directory/lib/python3.X/site-packages
set path = ('/my/gdal/install/directory/bin' $path)
```
