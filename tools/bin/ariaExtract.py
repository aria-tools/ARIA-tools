#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import argparse

# Import ARIA-tools functions
from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.extractProduct import cmdLineParse,export_products,prep_dem,merged_productbbox,finalize_metadata

if __name__ == '__main__':
    """
    Main driver
    """
    inps = cmdLineParse()

    print("***Extract Product Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    if not inps.layers:
        print ('No layers specified; only creating bounding box shapes')

    elif inps.layers.lower()=='all':
        print('All layers are to be extracted, pass all keys.')
        inps.layers=list(standardproduct_info.products[1][0].keys())
        # Must remove productBoundingBoxes, because it's not a raster layer
        inps.layers.remove('productBoundingBox')
        inps.layers.remove('productBoundingBoxFrames')
        # Must remove pair_name, because it's not a raster layer
        inps.layers.remove('pair_name')

    else:
        inps.layers=list(inps.layers.split(','))


    # extract/merge productBoundingBox layers for each pair and update dict,
    # report common track bbox (default is to take common intersection, but user may specify union), and expected shape for DEM.
    standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, arrshape, proj = merged_productbbox(standardproduct_info.products[1], os.path.join(inps.workdir,'productBoundingBox'), standardproduct_info.bbox_file, inps.croptounion)
    # Load mask (if specified).
    if inps.mask is not None:
        inps.mask=gdal.Warp('', inps.mask, options=gdal.WarpOptions(format="MEM", cutlineDSName=prods_TOTbbox, outputBounds=open_shapefile(standardproduct_info.bbox_file, 0, 0).bounds, dstNodata=0))
        inps.mask.SetProjection(proj)
        # If no data value
        if inps.mask.GetRasterBand(1).GetNoDataValue():
            inps.mask=np.ma.masked_where(inps.mask.ReadAsArray() == inps.mask.GetRasterBand(1).GetNoDataValue(), inps.mask.ReadAsArray())
        else:
            inps.mask=inps.mask.ReadAsArray()


    # Download/Load DEM & Lat/Lon arrays, providing bbox, expected DEM shape, and output dir as input.
    if inps.demfile is not None:
        # Pass DEM-filename, loaded DEM array, and lat/lon arrays
        inps.demfile, demfile, Latitude, Longitude = prep_dem(inps.demfile, standardproduct_info.bbox_file, prods_TOTbbox, proj, arrshape=arrshape, workdir=inps.workdir, outputFormat=inps.outputFormat)
    else:
        demfile=None ; Latitude=None ; Longitude=None

    # Extract user expected layers
    export_products(standardproduct_info.products[1], standardproduct_info.bbox_file, prods_TOTbbox, inps.layers, dem=demfile, lat=Latitude, lon=Longitude, mask=inps.mask, outDir=inps.workdir,outputFormat=inps.outputFormat, stichMethodType = inps.stichMethodType, verbose=inps.verbose)
