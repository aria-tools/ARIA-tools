# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & Brett Buzzanga & David Bekaert
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import glob
import logging
import os
import shutil

import affine
import copy
import osgeo.gdal
import pyproj
import rasterio
import tile_mate

import ARIAtools.util.shp
import ARIAtools.util.vrt


LOGGER = logging.getLogger(__name__)


def prep_mask(
        product_dict, maskfilename, bbox_file, prods_TOTbbox, proj,
        amp_thresh=None, arrres=None, workdir='./', outputFormat='ENVI',
        num_threads='2', multilooking=None, rankedResampling=False):
    """
    Function to load and export mask file with tile_mate
    """
    LOGGER.debug("prep_mask")

    # If specified DEM subdirectory exists, delete contents
    workdir = os.path.join(workdir, 'mask')
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)

    # Get bounds of user bbox_file
    bounds = ARIAtools.util.shp.open_shp(bbox_file).bounds

    # File must be physically extracted, cannot proceed with VRT format
    # Defaulting to ENVI format
    if outputFormat == 'VRT':
        outputFormat = 'ENVI'

    # Set output res
    if multilooking is not None:
        arrres = [arrres[0] * multilooking, arrres[1] * multilooking]

    # Download mask
    if maskfilename.lower() == 'download' or \
            maskfilename.lower() in tile_mate.stitcher.DATASET_SHORTNAMES:
        # if download specified, default to esa world cover mask
        if maskfilename.lower() == 'download':
            maskfilename = 'esa_world_cover_2021'
        lyr_name = copy.deepcopy(maskfilename)
        LOGGER.debug('Downloading water mask: %s', lyr_name)

        # set file names
        uncropped_maskfilename = os.path.join(workdir,
                                              f'{maskfilename}_uncropped.tif')
        maskfilename = os.path.join(workdir, f'{maskfilename}.msk')
        ref_file = os.path.join(workdir, 'tmp_referencefile')

        # download mask
        dat_arr, dat_prof = tile_mate.get_raster_from_tiles(
            bounds, tile_shortname=lyr_name)

        # fill permanent water body
        if lyr_name == 'esa_world_cover_2021':
            dat_arr[dat_arr == 80] = 0
            dat_arr[dat_arr != 0] = 1

        # assign datatype and set resampling mode
        dat_arr = dat_arr.astype('byte')
        f_dtype = 'uint8'
        resampling_mode = rasterio.warp.Resampling.nearest

        # get output parameters from temp file
        crs = pyproj.CRS.from_wkt(proj)
        osgeo.gdal.Warp(
            ref_file, product_dict[0], format=outputFormat,
            outputBounds=bounds, xRes=arrres[0], yRes=arrres[1],
            targetAlignedPixels=True, multithread=True,
            options=['NUM_THREADS=%s' % (num_threads)])
        with rasterio.open(ref_file) as src:
            reference_gt = src.transform
            resize_col = src.width
            resize_row = src.height

        # remove temporary file
        for j in glob.glob(ref_file + '*'):
            if os.path.isfile(j):
                os.remove(j)

        # save uncropped raster to file
        with rasterio.open(uncropped_maskfilename, 'w',
                           height=resize_row, width=resize_col, count=1,
                           dtype=f_dtype, crs=crs,
                           transform=affine.Affine(*reference_gt)) as dst:
            rasterio.warp.reproject(
                source=dat_arr, destination=rasterio.band(dst, 1),
                src_transform=dat_prof['transform'], src_crs=dat_prof['crs'],
                dst_transform=reference_gt, dst_crs=crs,
                resampling=resampling_mode)

        # save cropped mask with precise spacing
        osgeo.gdal.Warp(
            maskfilename, uncropped_maskfilename, format=outputFormat,
            outputBounds=bounds, outputType=osgeo.gdal.GDT_Byte,
            xRes=arrres[0], yRes=arrres[1], targetAlignedPixels=True,
            multithread=True, options=['NUM_THREADS=%s' % (num_threads)])

        update_file = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
        update_file.SetProjection(proj)
        update_file.GetRasterBand(1).SetNoDataValue(0.)
        osgeo.gdal.Translate(f'{maskfilename}.vrt',
                             maskfilename, format='VRT')

    # User specified mask
    else:
        LOGGER.debug("Using user specified mask %s" % maskfilename)
        # Path to local version of user specified mask
        user_mask = os.path.abspath(maskfilename)  # for clarity
        user_mask_n = os.path.basename(os.path.splitext(user_mask)[0])
        local_mask = os.path.join(workdir, f'{user_mask_n}.msk')
        local_mask_unc = os.path.join(workdir, f'{user_mask_n}_uncropped.msk')
        if user_mask == local_mask:
            LOGGER.debug(
                'The mask you specified already exists in %s, '
                'using the existing one...' % os.path.dirname(local_mask))

            # move all original files to temp path to circumvent gdal issues
            temp_workdir = os.path.join(workdir, 'tmp_dir')
            os.makedirs(temp_workdir, exist_ok=True)
            local_mask_noext = os.path.join(workdir, '%s.' %(user_mask_n))
            for j in glob.glob(local_mask_noext + '*'):
                shutil.move(j, temp_workdir)

            temp_local_mask = os.path.join(
                temp_workdir, '%s.msk' %(user_mask_n))
            ds = osgeo.gdal.Open(temp_local_mask)
     
            # remove temporary file
            shutil.rmtree(temp_workdir)

        else:
            # move the mask to the local directory and built a VRT for it
            osgeo.gdal.UseExceptions()

            # shutil.copy(user_mask, local_mask_unc)
            ds = osgeo.gdal.BuildVRT(
                f'{local_mask_unc}.vrt', user_mask, outputBounds=bounds)
            assert ds is not None, f'Could not open user mask: {user_mask}'

        # crop the user mask and write
        osgeo.gdal.Warp(
            f'{local_mask}', ds, format=outputFormat,
            cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrres[0],
            yRes=arrres[1], targetAlignedPixels=True, multithread=True,
            options=[f'NUM_THREADS={num_threads}'])

        # set projection of the local mask
        mask_file = osgeo.gdal.Open(local_mask, osgeo.gdal.GA_Update)
        mask_file.SetProjection(proj)

        # create vrt for local cropped mask
        osgeo.gdal.Translate(f'{local_mask}.vrt', local_mask, format='VRT')

        # assign new local mask for amp thresh
        maskfilename = local_mask

    # Make average amplitude mask
    if amp_thresh is not None:
        amp_file = ARIAtools.util.vrt.rasterAverage(
            os.path.join(workdir, 'avgamplitude'), product_dict, bounds,
            prods_TOTbbox, arrres, outputFormat=outputFormat,
            thresh=amp_thresh)

        # Update mask with average amplitude
        mask_file = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
        mask_arr = mask_file.ReadAsArray()
        mask_file.GetRasterBand(1).WriteArray(mask_arr * amp_file)

    # crop/expand mask to DEM size?
    osgeo.gdal.Warp(
        maskfilename, maskfilename, format=outputFormat,
        cutlineDSName=prods_TOTbbox, outputBounds=bounds, xRes=arrres[0],
        yRes=arrres[1], targetAlignedPixels=True, multithread=True,
        options=[f'NUM_THREADS={num_threads} -overwrite'])

    mask = osgeo.gdal.Open(maskfilename, osgeo.gdal.GA_Update)
    mask.SetProjection(proj)
    mask.SetDescription(maskfilename)
    mask_array = mask.ReadAsArray()
    mask_array[mask_array != 1] = 0
    mask.GetRasterBand(1).WriteArray(mask_array)

    # Update VRT
    translate_options = osgeo.gdal.TranslateOptions(format="VRT")
    osgeo.gdal.Translate(
        maskfilename + '.vrt', maskfilename, options=translate_options)

    # return filename of mask
    return maskfilename
