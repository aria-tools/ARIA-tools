#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import sys
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

# Suppress warnings
gdal.PushErrorHandler('CPLQuietErrorHandler')

# Import functions
from ARIAtools.ARIAProduct import ARIA_standardproduct
from ARIAtools.shapefile_util import open_shapefile
from ARIAtools.productPlot import cmdLineParse,plot_class


if __name__ == '__main__':
    '''
        Main driver.
    '''
    inps = cmdLineParse()

    print("***Plotting Function:***")
    # if user bbox was specified, file(s) not meeting imposed spatial criteria are rejected.
    # Outputs = arrays ['standardproduct_info.products'] containing grouped “radarmetadata info” and “data layer keys+paths” dictionaries for each standard product
    # In addition, path to bbox file ['standardproduct_info.bbox_file'] (if bbox specified)
    standardproduct_info = ARIA_standardproduct(inps.imgfile, bbox=inps.bbox, workdir=inps.workdir, verbose=inps.verbose)

    # If user requests to generate all plots.
    if inps.plotall:
        print('"-plotall"==True. All plots will be made.')
        inps.plottracks=True
        inps.plotbperp=True
        inps.plotcoh=True
        inps.plotbperpcoh=True
        inps.makeavgoh=True


    if inps.plottracks or inps.plotcoh or inps.makeavgoh or inps.plotbperpcoh:
        # Import functions
        from ARIAtools.extractProduct import merged_productbbox

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


    # Make spatial extent plot
    if inps.plottracks:
        print("- Make plot of track latitude extents vs bounding bbox/common track extent.")
        make_plot=plot_class([[j['productBoundingBox'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'productBoundingBox'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, croptounion=inps.croptounion)
        make_plot.plot_extents()


    # Make pbaseline plot
    if inps.plotbperp:
        print("- Make baseline plot and histogram.")
        make_plot=plot_class([[j['bPerpendicular'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'bPerpendicular'))
        make_plot.plot_pbaselines()


    # Make average land coherence plot
    if inps.plotcoh:
        print("- Make average IFG coherence plot in time, and histogram of average IFG coherence.")
        make_plot=plot_class([[j['coherence'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'coherence'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask)
        make_plot.plot_coherence()


    # Generate average land coherence raster
    if inps.makeavgoh:
        print("- Generate 2D raster of average coherence.")
        make_plot=plot_class([[j['coherence'] for j in standardproduct_info.products[1]], [j["pair_name"] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'coherence'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask, outputFormat=inps.outputFormat)
        make_plot.plot_avgcoherence()


    # Make pbaseline/coherence combo plot
    if inps.plotbperpcoh:
        print("- Make baseline plot that is color-coded with respect to mean IFG coherence.")
        make_plot=plot_class([[j['bPerpendicular'] for j in standardproduct_info.products[1]],  [j["pair_name"] for j in standardproduct_info.products[1]], [j['coherence'] for j in standardproduct_info.products[1]]], workdir=os.path.join(inps.workdir,'bPerpendicular'), bbox_file=standardproduct_info.bbox_file, prods_TOTbbox=prods_TOTbbox, mask=inps.mask)
        make_plot.plotbperpcoh()
