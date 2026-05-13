#! /usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha, David Bekaert, Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import logging

import ARIAtools.extractProduct
import ARIAtools.util.dem
import ARIAtools.util.log
import ARIAtools.util.vrt
import tile_mate
from ARIAtools.constants import ARIA_LAYERS

from aria_tools.core.workflows import (
    run_extract_workflow,
)
from aria_tools.errors import AriaToolsError

LOGGER = logging.getLogger("ariaExtract.py")


def createParser():
    """
    Extract specified product layers. The default will export all layers.
    """
    parser = argparse.ArgumentParser(
        description="Program to extract data and meta-data layers from "
        "standard GUNW products. Program will handle cropping/"
        "stitching when needed. By default, the program will crop "
        "all IFGs to bounds determined by the common intersection "
        "and bbox (if specified)"
    )

    parser.add_argument(
        "-f",
        "--file",
        dest="imgfile",
        type=str,
        required=True,
        help="List of Sentinel-1 GUNW or NISAR GUNW products "
        "(wildcards supported) or txt file with product urls "
        "for virtual access without downloading. For virtual "
        "processing a local metadata cache is created on "
        "first run; subsequent runs read from the cache for "
        "faster initialization.",
    )
    parser.add_argument(
        "-w",
        "--workdir",
        dest="workdir",
        default="./",
        help="Specify directory to deposit all outputs. Default is local "
        "directory where script is launched.",
    )
    parser.add_argument(
        "-l",
        "--layers",
        dest="layers",
        default=None,
        help="Specify the layers to extract as a comma-separated list enclosed "
        'in single quotes. Allowed values include: "unwrappedPhase", '
        '"coherence", "amplitude", "bPerpendicular", "bParallel", '
        '"incidenceAngle", "lookAngle", "azimuthAngle", "ionosphere", '
        '"troposphereWet", "troposphereHydrostatic", "troposphereTotal", '
        '"solidEarthTide". If left blank, only the bounding box will be '
        "extracted.",
    )
    parser.add_argument(
        "-tm",
        "--tropo_models",
        dest="tropo_models",
        type=str,
        default="all",
        help="Specify the weather model(s) to extract. "
        'The default is "all", which extracts all models in the product.',
    )
    parser.add_argument(
        "-d",
        "--demfile",
        dest="demfile",
        type=str,
        default=None,
        help='DEM file. To download new DEM, specify "Download".',
    )
    parser.add_argument(
        "-p",
        "--projection",
        dest="projection",
        default="4326",
        type=str,
        help="EPSG projection code for DEM. By default 4326. "
        'Specify "native" to pass most common '
        "projection from stack.",
    )
    parser.add_argument(
        "-b",
        "--bbox",
        dest="bbox",
        type=str,
        default=None,
        help="Provide either valid shapefile or Lat/Lon Bounding SNWE. -- "
        "Example : '19 20 -99.5 -98.5'",
    )
    parser.add_argument(
        "-m",
        "--mask",
        dest="mask",
        type=str,
        default=None,
        help="Specify either path to valid water mask, or "
        "download using one of the following "
        f"data sources: {tile_mate.stitcher.DATASET_SHORTNAMES}",
    )
    parser.add_argument(
        "-at",
        "--amp_thresh",
        dest="amp_thresh",
        default=None,
        type=str,
        help='Amplitude threshold below which to mask. Specify "None" to not '
        'use amplitude mask. By default "None".',
    )
    parser.add_argument(
        "-nt",
        "--num_threads",
        dest="num_threads",
        default="2",
        type=str,
        help="Specify number of threads for multiprocessing operation in "
        'gdal. By default "2". Can also specify "All" to use all '
        "available threads.",
    )
    parser.add_argument(
        "-of",
        "--outputFormat",
        dest="outputFormat",
        type=str,
        default="VRT",
        help='GDAL compatible output format (e.g., "ENVI", "GTiff"). By '
        "default files are generated virtually except for "
        '"bPerpendicular", "bParallel", "incidenceAngle", "lookAngle", '
        '"azimuthAngle", "unwrappedPhase" as these are require either '
        "DEM intersection or corrections to be applied",
    )
    parser.add_argument(
        "-croptounion",
        "--croptounion",
        action="store_true",
        dest="croptounion",
        help="If turned on, IFGs cropped to bounds based off of union and "
        "bbox (if specified). Program defaults to crop all IFGs to "
        "bounds based off of common intersection and bbox (if "
        "specified).",
    )
    parser.add_argument(
        "-ml",
        "--multilooking",
        dest="multilooking",
        type=int,
        default=None,
        help="Multilooking factor is an integer multiple of standard "
        "resolution. E.g. 2 = 90m*2 = 180m",
    )
    parser.add_argument(
        "-rr",
        "--rankedResampling",
        action="store_true",
        dest="rankedResampling",
        help="If turned on, IFGs resampled based off of the average of pixels "
        "in a given resampling window corresponding to the connected "
        "component mode (if multilooking specified). Program defaults "
        "to lanczos resampling algorithm through gdal (if multilooking "
        "specified).",
    )
    parser.add_argument(
        "-mo",
        "--minimumOverlap",
        dest="minimumOverlap",
        type=float,
        default=0.0081,
        help="Minimum km\u00b2 area of overlap of scenes wrt specified "
        "bounding box. Default 0.0081 = 0.0081km\u00b2 = area of single "
        "pixel at standard 90m resolution",
    )
    parser.add_argument(
        "-if",
        "--iono_filter",
        action="store_true",
        dest="iono_filter",
        help="Enable spatial filtering and quadratic surface approximation "
        "of the NISAR ionosphere layer. Caution: This may smooth out "
        "valid short-wavelength signals. (Note: This filter is always "
        "enforced for S1 GUNWs to mitigate large, unreliable artifacts).",
    )
    parser.add_argument(
        "--version",
        dest="version",
        default=None,
        help="Specify version as str, e.g. 2_0_4 or all prods; default: all",
    )
    parser.add_argument(
        "--nc_version",
        dest="nc_version",
        default="1b",
        help="Specify netcdf version as str, e.g. 1c or all prods; default: 1b",
    )
    parser.add_argument(
        "-verbose",
        "--verbose",
        action="store_true",
        dest="verbose",
        help="Toggle verbose mode on.",
    )
    parser.add_argument(
        "--log-level",
        choices=["debug", "info", "warning", "error"],
        default="info",
        help="Logger log level. Default: info.",
    )
    return parser


def main():
    """Main workflow for extracting layers from ARIA products."""
    # Parse command line args
    parser = createParser()
    args = parser.parse_args()

    log_level = {
        "debug": logging.DEBUG,
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
    }[args.log_level]
    logging.basicConfig(level=log_level, format=ARIAtools.util.log.FORMAT)
    LOGGER.info(f"ARIAtools version: {ARIAtools.__version__}")
    print("*****************************************************************")
    LOGGER.info("*** Extract Product Function ***")
    print("*****************************************************************")
    # Switch tropo models to None if troposphere isn't specified
    if args.layers is None or "tropo" not in args.layers:
        args.tropo_models = None
    # Check whether all necessary inputs were specified.
    # some products require a DEM to extract -- if any of those are requested,
    # ensure that a valid DEM is specified
    if args.layers is not None:
        # format list of layers
        layers = list(args.layers.split(","))
        layers = [i.replace(" ", "") for i in layers]
        layers = ["all" if layer.lower() == "all" else layer for layer in layers]
        layers = [
            "troposphere*" if layer.startswith("troposphere") else layer
            for layer in layers
        ]

        # list of layers requiring DEM for extraction
        LAYERS_REQUIRING_DEM = {
            "all",
            "bPerpendicular",
            "bParallel",
            "incidenceAngle",
            "lookAngle",
            "azimuthAngle",
            "solidEarthTide",
            "troposphere*",
        }

        # check that DEM is specified depending on layers requested
        if len(LAYERS_REQUIRING_DEM.intersection(layers)) > 0:
            if args.demfile is None:
                error_msg = (
                    "A valid DEM must be specified when extracting any of "
                    f"{', '.join(LAYERS_REQUIRING_DEM)}"
                )
                LOGGER.error(error_msg)
                raise AriaToolsError(error_msg)

    try:
        run_extract_workflow(args, logger=LOGGER, valid_layers=ARIA_LAYERS)
        print("*****************************************************************")
        LOGGER.info("*** Product extraction completed successfully ***")
        print("*****************************************************************")
    except Exception as e:
        print("*****************************************************************")
        LOGGER.error("*** Product extraction FAILED ***")
        LOGGER.error(f"Error: {e}")
        print("*****************************************************************")
        raise


if __name__ == "__main__":
    main()
