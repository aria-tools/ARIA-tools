# ARIA-tools
[![Language](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache_2.0-green.svg)](https://github.com/aria-tools/ARIA-tools/blob/master/LICENSE)
[![Version](https://img.shields.io/badge/version-2.0-orange.svg)](https://github.com/aria-tools/ARIA-tools/releases)

ARIA-tools is an open-source package in Python which contains tools to manipulate standard InSAR products from Sentinel-1 (ARIA GUNW-S1) and NISAR (NISAR_L2_GUNW). This software is open source under the terms of the [Apache 2.0 License](LICENSE). Its development was funded under the ROSES awards from the NASA Sea-level Change Team (NSLCT) program, the Earth Surface and Interior (ESI) program, and under the NISAR Science Team (NISAR-ST) program.

## Highlights

- Unified command surface through `aria-tools <command> [options]`
- Support for Sentinel-1 ARIA GUNW and NISAR GUNW products
- Virtual-access workflows through GDAL `/vsicurl/` URL lists
- Native Python process execution for extract and timeseries workflows
- Modern `pyproject.toml` packaging with one full install surface
- Install-first test and contributor workflow with `pytest`, `ruff`, and `pre-commit`


> [!IMPORTANT]
> ARIA-tools now supports **NISAR GUNW** products created routinely by the NASA-ISRO SAR mission. We welcome the community for reporting bugs using the issue ticket functionality

> [!NOTE]
> **Standardized Date and Phase Convention**
>
> To prevent confusion during time-series analysis, **all** ARIA-tools outputs enforce a consistent date and phase convention, regardless of the native format of the input product:
> * **Date Order:** ARIA-tools strictly uses a `ReferenceDate_SecondaryDate` naming convention, where the reference scene is the more recent pass and the secondary scene is the earlier pass. 
> * **Phase Sign:** Because the more recent acquisition is used as the reference, the unwrapped phase convention is standardized across all outputs. Negative phase differences indicate movement away from the sensor, and positive phase differences indicate movement towards the sensor.
> * **NISAR vs. S1:** While ARIA-S1-GUNW products natively follow this convention out of the box, NISAR GUNW products natively use the opposite date order. ARIA-tools automatically flips the dates and mathematical signs for NISAR inputs during extraction to ensure all your downstream time-series outputs are 100% consistent!

------

ARIA tools includes support for:
- ARIA Geocoded Unwrapped Interferogram from Sentinel-1 (GUNW-S1) can be downloaded for free from the [ASF DAAC vertex page](https://search.asf.alaska.edu/#/?dataset=SENTINEL-1%20INTERFEROGRAM%20(BETA)) selecting "ARIA S1 GUNW" under "datasets".  Users can also request free on-demand products through the [ASF on-demand system](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/). These products are added to the standard product archive. A log-on using the NASA Earthdata credentials is required to order new or download data from the archive. [Product Specification Document](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/#product-packaging).
-  NISAR Gecoded Unwrapped Interferogram products (NISAR_L2_GUNW) can be downloaded for free from the [ASF DAAC vertex page](https://search.asf.alaska.edu/#/?dataset=NISAR&prodConfig=PR&sciProducts=GUNW) selecting "NISAR" under "datasets" and selecting "GUNW" under the filter criteria. [Product Specification document](https://nisar-docs.asf.alaska.edu/gunw/). ARIA-tools also supports the **NISAR Private Ephemeral Archive**, which contains pre-public NISAR L2 products available to the NISAR Cal/Val team. Access requires valid and approved access through NASA Earthdata credentials configured in `~/.netrc`; when present, ephemeral archive products are automatically included in NISAR search results.


The ARIA-tools package includes functionality to crop/merge data and meta-data layers for multiple standard products, extraction of data and meta-data layers from these products, and the set-up and the preparation for time-series. 

Actual time-series processing is not supported in ARIA-tools. However, outputs are compatible with third-party time-series InSAR packages such the "Miami INsar Time-series software in PYthon" ([MintPy](https://github.com/insarlab/MintPy)).
<p align="center">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/Hawaii.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/CA.png">
<img height="250" src="https://github.com/aria-tools/ARIA-tools-docs/blob/master/images/EastCoast.png">
</p>

> [!CAUTION]
> THIS IS RESEARCH CODE PROVIDED TO YOU "AS IS" WITH NO WARRANTIES OF CORRECTNESS. USE AT YOUR OWN RISK.

## Contents
1.  [Requirements](#requirements)
2.  [Installation](#installation)
3.  [Running ARIA-tools](#running-aria-tools)
4.  [Virtual Access](#virtual-access)
5.  [Developer Workflow](#developer-workflow)
6.  [Documentation](#documentation)
7.  [Citation](#citation)
8.  [Contributors and community contributions](#contributors)

------

## Requirements
`pyproject.toml` defines the authoritative pip install surface. The repo-level
`environment.yml` provides the matching all-in conda environment used for local
development and CI.

Minimum requirements:

- Python >= 3.8
- [GDAL](https://www.gdal.org/) and its Python bindings >= 3.7.0
- [PROJ](https://github.com/OSGeo/proj) >= 6.0

The standard install includes runtime dependencies plus the contributor tools
used in this repository, including `pytest`, `ruff`, and `pre-commit`.

------
## Installation

### Installing a stable release from conda
Install the published package from conda-forge:

```bash
conda create --name ARIA-tools
conda activate ARIA-tools
mamba install aria-tools
```


### Installing the latest development branch
Clone the repo and create the contributor environment:

```bash
cd ~/tools
git clone https://github.com/aria-tools/ARIA-tools.git
cd ARIA-tools
mamba env create -f environment.yml
conda activate ARIA-tools
```

Install the package in editable mode:

```bash
python -m pip install -e .
```

To keep the conda environment aligned with the repo over time:

```bash
mamba env update --name ARIA-tools --file environment.yml --prune
```

The primary command surface is `aria-tools <command> [options]`. Legacy script
entry points remain installed during the migration period.

```bash
aria-tools download --track 004 --output count
aria-tools extract -f "products/*.nc" -d Download
aria-tools timeseries -f "products/*.nc"
```

### Other installation options
If you need to build third-party geospatial dependencies from source, see:

- [LinuxSourceBuild.md](LinuxSourceBuild.md)
- [MacOSSourceBuild.md](MacOSSourceBuild.md)

To use Earthdata-backed archive access, configure `~/.netrc`:

```bash
echo "machine urs.earthdata.nasa.gov login myUsername password myPassword" > ~/.netrc
chmod 600 ~/.netrc
```

------
## Running ARIA-tools

ARIA-tools provides a unified command-line interface through the `aria-tools` command. The package is highly modular and allows for building custom processing workflows. Below, we show how to use the main commands. For detailed documentation, examples, and Jupyter notebooks see the [ARIA-tools-docs repository](https://github.com/aria-tools/ARIA-tools-docs). We welcome the community to contribute examples (see [CONTRIBUTING.md](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md) for instructions).

### Quick Start

The primary interface is `aria-tools <command> [options]`. Available commands:

```bash
# Get top-level help
aria-tools --help
aria-tools --version

# Get subcommand help
aria-tools extract --help
aria-tools timeseries --help

# Download products (or generate URL list)
aria-tools download --track 004 --output count

# Extract layers from products
aria-tools extract -f "products/*.nc" -w workdir -l unwrappedPhase,coherence

# Prepare time-series stack
aria-tools timeseries -f "products/*.nc" -w workdir

# Generate quality plots
aria-tools plot -f "products/*.nc" -w workdir
```

### Command Overview

| Command | Purpose |
| --- | --- |
| `aria-tools download` | Search, count, download, or emit URL lists for products |
| `aria-tools extract` | Crop, stitch, and export product layers |
| `aria-tools timeseries` | Prepare stacks and supporting outputs for time-series workflows |
| `aria-tools plot` | Generate QC and baseline plots |
| `aria-tools order` | Request on-demand Sentinel-1 products through HyP3 |
| `aria-tools misclosure` | Analyze phase-triplet misclosure |
| `aria-tools aoi` | Assist AOI generation from ASF metadata |
| `aria-tools kml2box` | Convert KML/KMZ polygons to GeoJSON bounding boxes |

### Commandline download of GUNW Products
ARIA GUNW-S1/NISAR_L2_GUNW products can be downloaded through the command line using `aria-tools download`, which wraps the ASF DAAC API. 

**Virtual Access Mode**: Use `aria-tools download -o url` to create a `.txt`
file with HTTPS URLs to archived products instead of downloading them. This
URL file can be passed directly to `aria-tools extract`, `aria-tools
timeseries`, and `aria-tools plot` for streaming access without local
downloads.

Example:
```bash
# Download products
aria-tools download --track 004 --start 20200101

# Or generate URL list for virtual access
aria-tools download --track 004 --start 20200101 --output url
```

### Manipulating GUNW Products
ARIA GUNW-S1/NISAR_L2_GUNW products can be manipulated (cropped, stitched, extracted) using `aria-tools extract`.

Example:
```bash
# Extract specific layers with DEM
aria-tools extract -f "products/*.nc" -w workdir \
  -l unwrappedPhase,coherence,amplitude \
  -d Download -b "33 35 -118 -116"
```

### Baseline and quality control plots for GUNW Products
Quality and baseline plots for spatial-temporal contiguous interferograms can be generated using `aria-tools plot`.

Example:
```bash
aria-tools plot -f "products/*.nc" -w workdir
```

### Time-series set-up of GUNW Products
Time-series preparation with spatial-temporal contiguous unwrapped interferograms and coherence can be done using `aria-tools timeseries`.

Example:
```bash
aria-tools timeseries -f "products/*.nc" -w workdir
```

### Ordering ARIA S1 GUNW Products on demand
On-demand ordering of **ARIA Sentinel-1 GUNW** products is supported through `aria-tools order`, which interfaces with the [ASF HyP3 on-demand processing system](https://hyp3-docs.asf.alaska.edu/guides/gunw_product_guide/). This is for GUNW-S1 products only (not NISAR) and allows users to build and order additional interferometric pairs not yet in the ASF archive, using your monthly HyP3 credit quota. A `~/.netrc` file with NASA Earthdata credentials is required for authentication.

Example:
```bash
aria-tools order --track 004 --start 20200101 --end 20200201
```

> [!NOTE]
> **Legacy Script Compatibility**: For users migrating from ARIA-tools v1,
> legacy script entry points (`ariaDownload.py`, `ariaExtract.py`, etc.)
> remain available during the transition period. See
> [docs/MIGRATION.md](docs/MIGRATION.md) for command mapping details.

> [!NOTE]  
> We support extraction of correction layers (e.g. Troposphere, Ionosphere, Solid Earth Tides) as well as geometry information (e.g. incidence angle, look angle, baselines, etc) embeded within the GUNW products 

------
## Virtual Access
GDAL Virtual File System capabilities can be used to process archived products
without downloading local `.nc` or `.h5` files first. This workflow is
supported for both Sentinel-1 GUNW and NISAR GUNW products.

Typical flow:

```bash
aria-tools download --track 004 --output url
aria-tools extract -f products.txt -w workdir
aria-tools timeseries -f products.txt -w workdir
```

Notes:

- URL-list workflows rely on GDAL `/vsicurl/` access to the ASF archive
- ARIA-tools creates a local `aria_meta_cache.json` sidecar on first use to
  cache remote product metadata and speed up subsequent runs
- The first remote run is expected to be slower while the metadata cache is
  being populated
- Full live validation still depends on network access and appropriate
  Earthdata credentials

------
## Developer Workflow
For local development and CI-style checks, install the package in the active
`ARIA-tools` environment:

```bash
python -m pip install -e .
```

Common validation commands:

```bash
python -m ruff check src/aria_tools tests/unit tests/integration
python -m ruff format --check src/aria_tools tests/unit tests/integration
python -m pytest tests/unit -q
python -m pytest tests/integration -q
python -m pytest tests/regression/validate_test.py -q
```

The test suite uses markers to separate quick offline checks from slower or
environment-dependent workflows:

- `offline`
- `slow`
- `network_required`
- `credentialed`

Examples:

```bash
python -m pytest tests/regression -q --run-slow --run-network --run-credentialed
pre-commit run --all-files
```

------
## Documentation
See the [ARIA-tools-docs repository](https://github.com/aria-tools/ARIA-tools-docs) for tutorials and notebook material.

For repo-local modernization notes and implementation details, see:

- [docs/MIGRATION.md](docs/MIGRATION.md)
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)

------
## Citation
Buzzanga, B., Bekaert, D. P. S., Hamlington, B. D., & Sangha, S. S. (2020). Towards Sustained Monitoring of Subsidence at the Coast using InSAR and GPS: An Application in Hampton Roads, Virginia. Geophysical Research Letters, 47, e2020GL090013. [https://doi.org/10.1029/2020GL090013](https://doi.org/10.1029/2020GL090013)

------
## Contributors
-   David Bekaert
-   Simran Sangha
-   Emre Havazli
-   Brett Buzzanga
-   Alexander Fore
-   Marin Govorcin
-   Charles Marshak
-   Joseph Kennedy
-   [other community members](https://github.com/aria-tools/ARIA-tools/graphs/contributors)

We welcome community contributions. For instructions see [here](https://github.com/aria-tools/ARIA-tools/blob/master/CONTRIBUTING.md).
