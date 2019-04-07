## ARIA-tools

[![Language](https://img.shields.io/badge/python-3.5%2B-blue.svg)](https://www.python.org/)
[![Latest version](https://img.shields.io/badge/latest%20version-v1.0.0--dev-green.svg)](https://github.com/dbekaert/ARIA-tools/blob/master/docs/download.md)
[![License](https://img.shields.io/badge/license-GPL-yellow.svg)](https://github.com/dbekaert/ARIA-tools/master/LICENSE)

ARIA-tools is an open-source package in Python which contains tools to manipulate ARIA standard InSAR products. For a full overview of available ARIA standard products and their specification see the products page on the [ARIA website](https://aria.jpl.nasa.gov). Currently, support for the ARIA Geocoded Unwrapped Interferogram (GUNW) product is included. Products can be download for free from the [ARIA-products page](https://aria-products.jpl.nasa.gov) and the [ASF DAAC vertex page](https://vertex.daac.asf.alaska.edu/#) under missions and beta-products, but require log-on using the NASA Earthdata credentials.
The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series programs such as GIAnT and [PySAR](https://github.com/yunjunz/PySAR). 


### 1. Download

To download the ARIA-tools repository:
```
git clone https://github.com/dbekaert/ARIA-tools.git
```

To download the ARIA-tools documentation repository:
```
git clone https://github.com/dbekaert/ARIA-tools-docs.git
```

### 2. Installation

TODO: PythonPATH setup, dependencies and requirements (e.g. GDAL >2.4, RELAX), split in example for macports and anaconda with the required pacakges and portfiles.

### 3. Running ARIA-tools

The ARIA-tools scripts are highly modulized in Python and therefore allows for building your own processing workflow. Below, we show how to call some of the functionality. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/dbekaert/ARIA-tools-docs/README.md).

#### Manipulating GUNW Products

TODO: list the commandline parsing of the main script which allows for cropping/stiching and layer-extraction, link to the documentation page for a note-book example

#### Baseline and quality control plots for GUNW Products

TODO: list the commandline parsing of the main script, link to the documentation page for a note-book example

#### Time-series set-up of GUNW Products

TODO: list the commandline parsing of the main script, link to the documentation page for a note-book example


### 4. Documentation

See the ARIA-tools-docs repo for all documentation:
+ [Jupyter Notebook Tutorials](https://github.com/dbekaert/ARIA-tools-docs/tree/master/docs/Notebooks.md)
+ [Overview of tool functions](https://github.com/dbekaert/ARIA-tools-docs/tree/master/docs/Programs.md)   

### 5. Citation

TODO


### Contributors    

* David Bekaert
* Simran Sangha
* [_other community members_](https://github.com/dbekaert/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/dbekaert/ARIA-tools/CONTRIBUTING.md)

