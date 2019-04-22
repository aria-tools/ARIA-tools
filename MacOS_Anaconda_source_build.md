## MacOSX
Here we provide guidelines on how to build GDAL 2.5+ and PROJ 4 (6.X.X) from source for **MacOS** machines. For **Linux** see the respective installation instructions.

### 0. XCode and Command Line Tools
Make sure that you have XCode and Command Line Tools installed. XCode 10.X versions do not install the header files under /usr/include but we need that directory for the installation.
For that purpose, after installing XCode and Command Line Tools run:
```
open /Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg
```
this command will install some header files under /usr/include

Since some required header files are missing we will need to point our environment to C++ header files under Xcode
```
setenv CPATH /Library/Developer/CommandLineTools/usr/include/c++/v1
```
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
conda install autoconf automake libtool numpy netcdf4 matplotlib pandas sqlite pkg-config shapely postgresql libcxx lapack proj4=6.0 --yes
```

Setup your shell environment to include python libraries before starting GDAL installation.
An example of C-shell would look like below
```
setenv LD_LIBRARY_PATH /my/python/directory/lib
```

### 2. GDAL SETUP
Clone the GDAL repository from github with a version of at least 2.5 (i.e. main branch).
```
git clone https://github.com/OSGeo/gdal
```

Please make sure you followed the instructions in section 0
```
setenv CPATH /Library/Developer/CommandLineTools/usr/include/c++/v1
```

Build the GDAL package with the python bindings:
```
cd gdal/gdal/
./configure --with-proj=/my/python/install/directory --prefix=/my/gdal/install/directory --with-python
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
setenv LD_LIBRARY_PATH /my/python/directory/lib
setenv GDAL_DATA /my/gdal/install/directory/share/gdal
set path = ('/my/gdal/install/directory/bin' $path)
```
