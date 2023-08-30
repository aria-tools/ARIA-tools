#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Marin Govorcin
# Copyright 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
import xarray as xr
from numpy.typing import NDArray

from typing import List, Optional, Tuple
from osgeo import gdal
from pathlib import Path


# import util modules
from ARIAtools.stitiching_util import (get_GUNW_array, get_GUNW_attr,
                                       frame_overlap, combine_data_to_single,
                                       write_GUNW_array, snwe_to_extent,
                                       _nan_filled_array)


GUNW_LAYERS = {'unwrappedPhase': 'NETCDF:"%s":/science/grids/data/unwrappedPhase',
               'coherence': 'NETCDF:"%s":/science/grids/data/coherence',
               'connectedComponents': 'NETCDF:"%s":/science/grids/data/connectedComponents',
               'ionosphere': 'NETCDF:"%s":/science/grids/corrections/derived/ionosphere/ionosphere'}


def fit_surface(data, order=2):
    dshape = data.shape
    length, width = dshape[-2:]
    mask_in = np.ones((length, width), dtype=np.float32)
    mask = (mask_in != 0).flatten()

    # for 2d only
    data = data.reshape(-1, 1)

    mask *= ~np.isnan(data.flatten())
    mask *= (data.flatten() != 0.)

    # design matrix
    xx, yy = np.meshgrid(np.arange(0, width),
                         np.arange(0, length))
    xx = np.array(xx, dtype=np.float32).reshape(-1, 1)
    yy = np.array(yy, dtype=np.float32).reshape(-1, 1)
    ones = np.ones(xx.shape, dtype=np.float32)
    
    # Bilinear
    if order==1.5:
        G = np.hstack((yy, xx, yy*xx, ones))
    elif order==2:
        # Quadratic
        G = np.hstack((yy**2, xx**2, yy*xx, yy, xx, ones))

    # estimate ramp
    X = np.dot(np.linalg.pinv(G[mask, :], rcond=1e-15), data[mask, :])
    surface = np.dot(G, X)
    surface = np.array(surface, dtype=data.dtype)

    return surface.reshape(dshape)


def _get_overlay(xr_ds1, xr_ds2):
    S = max(xr_ds1.y.min().values, xr_ds2.y.min().values)
    N = min(xr_ds1.y.max().values, xr_ds2.y.max().values)
    W = max(xr_ds1.x.min().values, xr_ds2.x.min().values)
    E = min(xr_ds1.x.max().values, xr_ds2.x.max().values)
    return S, N, W, E


def _get_median_offsets2frames(xr_data_list, xr_mask_list, ix1, ix2):
    S, N, W, E = _get_overlay(xr_data_list[ix1], xr_data_list[ix2])
    
    # Get overlap
    cropped_ds1 = xr_data_list[ix1].ionosphere.sel(y=slice(N, S), x=slice(W, E)).copy()
    cropped_mask1 = xr_mask_list[ix1].mask.sel(y=slice(N, S), x=slice(W, E)).copy()
    ds1 = np.ma.masked_array(cropped_ds1.values, mask=~cropped_mask1.values)

    cropped_ds2 = xr_data_list[ix2].ionosphere.sel(y=slice(N, S), x=slice(W, E)).copy()
    cropped_mask2 = xr_mask_list[ix2].mask.sel(y=slice(N, S), x=slice(W, E)).copy()
    ds2 = np.ma.masked_array(cropped_ds2.values, mask=~cropped_mask2.values)

    return np.nanmedian((ds1 - ds2).filled(fill_value=np.nan))

def stitch_ionosphere_frames(input_iono_files: List[str],
                             direction_N_S: Optional[bool] = True):
                             
    # Initalize variables for raster attributes
    iono_attr_list = []  # ionosphere raster metadata
    iono_xr_list = []
    mask_xr_list =[]

    # Loop through files
    for iono_file in input_iono_files:
        filename = iono_file.split(':')[1]
        iono_attr_list.append(get_GUNW_attr(iono_file))
        iono_xr = xr.open_dataset(iono_file, engine='rasterio').squeeze()

        # Generate mask using unwrapPhase connectedComponents
        mask_xr = xr.open_dataset(GUNW_LAYERS['connectedComponents'] % filename, 
                                  engine='rasterio').squeeze()
        mask = np.bool_(mask_xr.connectedComponents.data != 0)
        mask_xr['connectedComponents'].values = mask
        mask_xr = mask_xr.rename_vars({'connectedComponents':'mask'})
        # Interpolate to iono grid
        mask_xr = mask_xr.interp_like(iono_xr)
        
        iono_xr_list.append(iono_xr)
        mask_xr_list.append(mask_xr)

    # Remove intermidate variables
    del iono_xr, mask_xr, mask
    
    # Get SNWE and LATLON_SPACING
    SNWE = np.vstack([d['SNWE'] for d in iono_attr_list])
    LATLON = np.vstack([[d['LAT_SPACING'], d['LON_SPACING']] for d in iono_attr_list])

    # get sorted indices for frame bounds, from South to North
    sorted_ix = np.argsort(SNWE[:, 0], axis=0)

    if direction_N_S:
        sorted_ix = sorted_ix[::-1]

    # Step 1: adjusted frames using the median offset in the overlap region
    #         by using only reliable areas (connctedComponents != 0)
    for ix1, ix2 in zip(sorted_ix[:-1], sorted_ix[1:]):
        diff = _get_median_offsets2frames(iono_xr_list, mask_xr_list, ix1, ix2)
        iono_xr_list[ix2]['ionosphere'] += diff

    # Step 2: Merged ionosphere and mask datasets
    data_list = [d.ionosphere.data for d in iono_xr_list]
    mask_list = [d.mask.data for d in mask_xr_list]

    combined_iono = combine_data_to_single(data_list,
                                           SNWE.tolist(),
                                           LATLON.tolist(),
                                           method = 'mean',
                                           latlon_step=LATLON[0,:].tolist())

    combined_mask = combine_data_to_single(mask_list,
                                           SNWE.tolist(),
                                           LATLON.tolist(),
                                           method = 'min',
                                           latlon_step=LATLON[0,:].tolist())


    # Step 3: Fit quadratic surface 
    # Mask combined_iono before surface fitting
    combined_iono_msk = combined_iono[0].copy()
    mask = ~np.nan_to_num(combined_mask[0], 0).astype(np.bool_)
    combined_iono_msk[mask] = np.nan

    # Get surface
    surface = fit_surface(combined_iono_msk)
    surface = np.ma.masked_array(surface, mask=np.isnan(combined_iono[0]))
    surface = surface.filled(fill_value=0.)

    return surface, combined_iono[1], combined_iono[2]

## MAIN

def export_ionosphere(input_iono_files: List[str],
                      arrres: List[float],
                      output_iono: Optional[str] = './ionosphere',
                      output_format: Optional[str] = 'ISCE',
                      bounds: Optional[tuple] = None,
                      clip_json: Optional[str] = None,
                      mask_file: Optional[str] = None,
                      verbose: Optional[bool] = False,
                      overwrite: Optional[bool] = True) -> None:


        # Outputs
        output_iono = Path(output_iono).absolute()
        if not output_iono.parent.exists():
            output_iono.parent.mkdir()

        # create temp files
        temp_iono_out = output_iono.parent / ('temp_' + output_iono.name)

        # Create VRT and exit early if only one frame passed,
        # and therefore no stitching needed
        if len(input_ionp_files) == 1:
            gdal.BuildVRT(str(temp_iono_out.with_suffix('.vrt')),
                          input_ionp_files)
        
        else:
            (combined_iono, 
            snwe, latlon_spacing) = stitch_ionosphere_frames(input_ionp_files,
                                                             direction_N_S=True)

            # Write
            # write stitched ionosphere
            write_GUNW_array(
                temp_iono_out, combined_iono, snwe,
                format=output_format, verbose=verbose,
                update_mode=overwrite, add_vrt=True, nodata=0.0)
        
        # Crop
        [print(f'Cropping to {bounds}') if verbose and bounds else None]
        if overwrite:
            [print(f'Removing {output_iono}') if verbose else None]
            output_iono.unlink(missing_ok=True)
        
        
        # Crop if selected
        ds = gdal.Warp(str(output_iono),
                       str(temp_iono_out.with_suffix('.vrt')),
                       format=output_format,
                       cutlineDSName=clip_json,
                       xRes=arrres[0], yRes=arrres[1],
                       targetAlignedPixels=True,
                       # cropToCutline = True,
                       outputBounds=bounds
                       )
        ds = None
        # Update VRT
        [print(f'Writing {output_iono}, {output.with_suffix(".vrt")}')
         if verbose else None]
        gdal.Translate(str(output_iono.with_suffix('.vrt')),
                       str(output_iono), format="VRT")
        # Remove temp files
        [ii.unlink() for ii in [temp_iono_out, temp_iono_out.with_suffix('.vrt'),
                                temp_iono_out.with_suffix('.xml'),
                                temp_iono_out.with_suffix('.hdr'),
                                temp_iono_out.with_suffix('.aux.xml')] if ii.exists()]

        # Mask
        if mask_file:
            if isinstance(mask_file, str):
                mask = gdal.Open(mask_file)
            else:
                # for gdal instance, from prep_mask
                mask = mask_file

            mask_array = mask.ReadAsArray()
            array = get_GUNW_array(str(output_iono.with_suffix('.vrt')))
            update_array = mask_array * array

            update_file = gdal.Open(str(output_iono), gdal.GA_Update)
            update_file = update_file.GetRasterBand(1).WriteArray(update_array)
            update_file = None

 