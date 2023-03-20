#!/usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Marin Govorcin
# Copyright 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import warnings
import numpy as np
from numpy.typing import NDArray, ArrayLike

from typing import List, Optional, Tuple, Union
from osgeo import gdal, osr, gdal_array
from pathlib import Path

from matplotlib import pyplot as plt
import matplotlib as mpl

'''
NOTE:
Sequential stitcher relies on connected components [region of pixels] in the overlap between two frames with 
assumption that each component is unwrapped correctly by SNAPHU. This might not always be the case. If there 
are unwrapping errors within the overlapping component, they will propage to the stitched image.
Potential fix could be to manually create connected component around the unwrapping error and try to
re-run stitching.

Connected Component 0 represents masked pixels by SNAPHU. Stitcher merges 0 overlapping components by default 
for the sake of vizualization. However, these unwrapped phase pixels are unreliable and will often be 
misaligned as 2-pi integer cycles shift does not apply here. User is advised to mask these pixes for further processing.

TODO: Re-enumeration of connected components in the stitched image. Due the merging of overlapping components, there
are gaps in enumeration of connected component labels. This should not affect the futher processing, but for the sake
of consistency, add function to re-enumerate components.

DISCLAIMER : This is development script. Requires some additional clean-up and restructuring
'''

################### STITCHING SUBROUTINES ##################################
def stitch_2frames(unw_data1 : NDArray, conncomp_data1 : NDArray, rdict1 : dict,  
                   unw_data2 : NDArray, conncomp_data2 : NDArray, rdict2 : dict,
                   correction_method : Optional[str] = 'cycle2pi',
                   range_correction : Optional[bool] = False, 
                   verbose: Optional[bool] = False) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Function to sequentially stitch two frames along the same track. Mean offset or 2pi integer cycles are 
    estimated for each overlapping connected components, which after correction are merged into the same 
    component. Code runs forward and backward, first correcting Northern frame with respect to Southern,
    and then the opposite if corrected components overlap with multiple components in Southern frame.

    Parameters
    ----------
    unw_data1 : array
        array containing unwrapped phase in Frame-1 (South)
    conncomp_data1 : array
        array containing connected components in Frame-1 (South)
    rdict1 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..] for Frame-1 (South)
    unw_data2 : array
        array containing unwrapped phase in Frame-2 (North)
    conncomp_data2 : array
        array containing connected components in Frame-2 (North)
    rdict2 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..] for Frame-2 (North)
    correction_method : str 
        method to correct offset between overlapping connected components. Available options:
         "meanoff" - mean offset between components
         "cycle2pi" - 2pi integer cycles between components 
    verbose : bool
        print info messages [True/False]
                                            
    Returns
    -------
    unw_data1 : array
        corrected array containing unwrapped phase in Frame-1 (South)
    conncomp_data1 : array
        corrected array containing connected components in Frame-1 (South)
    unw_data2 : array
        corrected array containing unwrapped phase in Frame-2 (North)
    conncomp_data2 : array
        corrected array containing connected components in Frame-2 (North)

    TODO: combine corrected array in this function, return only the stitch unwrapped Phase and connectedComponent
    """

    # Adjust connected Component in Frame 2 to start with last component number in Frame-1
    conncomp_data2 = conncomp_data2 + np.nanmax(conncomp_data1)
    conncomp_data2[conncomp_data2==np.nanmax(conncomp_data1)] = 0.0

    ###### GET FRAME OVERLAP ##########
    box_1, box_2  = frame_overlap(rdict1['SNWE'], 
                                  rdict2['SNWE'],
                                  [rdict1['LAT_SPACING'], rdict1['LON_SPACING']],
                                  [rdict2['LAT_SPACING'], rdict2['LON_SPACING']])

    ###### LOOP OVER COMPONENTS WITHIN THE OVERLAP ########## 
    # Get connected component pairs
    if verbose:
        print('\nGetting overlapping components')
        
    conn_comp_pairs, conn_comp_reverse = get_overlaping_conncomp(conncomp_data1[box_1],
                                                                 conncomp_data2[box_2])

    # Forward correction
    for pair in conn_comp_pairs:
        diff, cycles2pi, range_corr = _integer_2pi_cycles(unw_data1[box_1], conncomp_data1[box_1], np.float32(pair[1]),
                                                         unw_data2[box_2], conncomp_data2[box_2], np.float32(pair[0]), 
                                                         range_correction=range_correction, print_msg=verbose)

        # Correction methods: mean difference, 2pi integer cycles
        if correction_method == 'cycle2pi':
            correction = cycles2pi
        elif correction_method == 'meandiff':
            correction = diff 
        else:
            raise ValueError(f'Wrong correction method {correction_method},',
                                f'Select one of available: "cycle2pi", "meandiff"')
        
        if correction:
            if range_correction: correction += range_corr # add range correction 

            unw_data2[conncomp_data2 == np.float32(pair[0])] += correction
            conncomp_data2[conncomp_data2 == np.float32(pair[0])] = np.float32(pair[1])

            #  Correct reverse pairs as we updated connected components
            if conn_comp_reverse.any():
                conn_comp_reverse[conn_comp_reverse[:,0] == pair[0], 0] = np.float32(pair[1]) 
        else:
            print('SKIPPED!:', correction, pair, '\n')
    
    # Backward correction
    for pair in conn_comp_reverse:
        diff, cycles2pi, range_corr = _integer_2pi_cycles(unw_data1[box_1], conncomp_data1[box_1], np.float32(pair[1]),
                                                          unw_data2[box_2], conncomp_data2[box_2], np.float32(pair[0]), 
                                                          range_correction=range_correction, print_msg=verbose)
        
        # Correction methods: mean difference, 2pi integer cycles
        if correction_method == 'cycle2pi':
            correction = cycles2pi
        elif correction_method == 'meanoff':
            correction = diff 
        else:
            raise ValueError(f'Wrong correction method {correction_method},',
                                f'Select one of available: "cycle2pi", "meanoff"')
        if correction:
            if range_correction: correction += range_corr # add range correction 

            unw_data1[conncomp_data1 == np.float32(pair[1])] -= correction
            conncomp_data1[conncomp_data1 == np.float32(pair[1])] = np.float32(pair[0])
                
    return unw_data1, unw_data2, conncomp_data1, conncomp_data2


def stitch_2frames_metadata(unw_data1 : NDArray, rdict1 : dict,  
                   unw_data2 : NDArray, rdict2 : dict,
                   verbose: Optional[bool] = False) -> Tuple[NDArray, NDArray]:
    """Sequential stitching function implementation from `stitch_2frames`
    for metadata layers

    Parameters
    ----------
    unw_data1 : array
        array containing unwrapped phase in Frame-1 (South)
    rdict1 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..] for Frame-1 (South)
    unw_data2 : array
        array containing unwrapped phase in Frame-2 (North)
    rdict2 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..] for Frame-2 (North)
    verbose : bool
        print info messages [True/False]
                                            
    Returns
    -------
    unw_data1 : array
        corrected array containing unwrapped phase in Frame-1 (South)
    unw_data2 : array
        corrected array containing unwrapped phase in Frame-2 (North)

    TODO: combine corrected array in this function, return only the stitch unwrapped Phase
    """

    ###### GET FRAME OVERLAP ##########
    box_1, box_2  = frame_overlap(rdict1['SNWE'], 
                                  rdict2['SNWE'],
                                  [rdict1['LAT_SPACING'], rdict1['LON_SPACING']],
                                  [rdict2['LAT_SPACING'], rdict2['LON_SPACING']],
                                  [-0.1,0.1])

    ###### EXAMINE THE OVERLAP ##########
    for i in range(unw_data1.shape[0]):
        diff = _metadata_offset(unw_data1[i][box_1],
                            unw_data2[i][box_2], 
                            print_msg=verbose)
        
        unw_data2[i] += diff
                
    return unw_data1, unw_data2


def get_overlaping_conncomp(conn_comp1 : NDArray, 
                            conn_comp2 : NDArray) -> Tuple[NDArray, NDArray]:
    """
    Get forward and backward pairs of overlapping connected components with the 
    number of overlaping pixels.
    Return [connectComponent-Fram2, connectComponent-Frame1, number of pixels]

    Parameters
    ----------
    conn_comp1 : array
        array containing connected components in Frame-1 (South)
    conn_comp2 : array
        array containing connected components in Frame-2 (North)

    Returns
    -------
    conn_pairs : array
        forward pairs of overlapping components
    conn_pairs_reverse : array
        backward pairs of overlapping components
    """

    conn_comp_union = np.empty((0,3), dtype=np.int32)

    # Get unique components
    concomp1 = np.unique(conn_comp1)
    concomp2 = np.unique(conn_comp2)

    # Loop through them and connect size and number of overlapping data
    for ix2 in concomp2:
        for ix1 in concomp1:
            if ~np.isnan(ix1) and ~np.isnan(ix2):
                data1 = conn_comp1.copy()
                data2 = conn_comp2.copy() 

                # Use selected component data, mask other components
                data1[data1!=ix1] =np.nan
                data2[data2!=ix2] =np.nan
                # Number of overlapping points
                npoints_overlap = np.nansum(~np.isnan(data1 - data2))

                if npoints_overlap > 0: 
                    # Array : conn-label Frame-2, coon-label Frame-1, # points overlap
                    connecting_array = np.array([ix2, ix1, npoints_overlap], dtype=np.int32)
                    # Skip 0 component combination with other components
                    if not ix1==0.0 and not ix2==0.0:
                        conn_comp_union = np.concatenate((conn_comp_union, 
                                                          np.atleast_2d(connecting_array)), axis=0)
                    # Get 0 components in both frames 
                    elif ix1==0.0 and ix2==0.0:
                        conn_comp_union = np.concatenate((conn_comp_union, 
                                                          np.atleast_2d(connecting_array)), axis=0)


    # Find components to correct in Frame 2 
    comp2correct = np.unique(conn_comp_union[:,0])

    conn_pairs = np.empty((0,3), dtype=np.int32)
    conn_pairs_reverse = np.empty((0,3), dtype=np.int32)
    
    for component in comp2correct:
        # find number of times components is referenced
        count = np.sum(~np.isnan(np.where(conn_comp_union[:,0] == component)[0]))
        
        if count > 1:
        #find the one with most points in overlap
            max_points = np.max(conn_comp_union[np.where(conn_comp_union[:,0] == component)[0], 2])
            # store it for forward stitching
            ix = np.where(conn_comp_union[:,2] == max_points)[0]
            conn_pairs = np.concatenate((conn_pairs, 
                                         np.atleast_2d(conn_comp_union[ix,:])), axis=0)
            
            # store it for backward stitching
            iy = np.where((conn_comp_union[:,0] == component) & (conn_comp_union[:,2] != max_points))[0]
            conn_pairs_reverse = np.concatenate((conn_pairs_reverse, 
                                                 np.atleast_2d(conn_comp_union[iy, :])), axis=0)

        # If count is 1 store the connecting component pair
        else:
            ik = np.where(conn_comp_union[:,0] == component)[0]
            conn_pairs = np.concatenate((conn_pairs, np.atleast_2d(conn_comp_union[ik,:])), axis=0) 

    return conn_pairs, conn_pairs_reverse 

def _integer_2pi_cycles(unw1 : NDArray, concom1 : NDArray, ix1 : np.float32,
                        unw2 : NDArray, concom2 : NDArray, ix2 : np.float32,
                        range_correction : Optional[bool] = False,
                        print_msg : Optional[bool] = True,) -> Tuple[np.float32, np.float32, np.float32]:
    """
    Get mean difference of unwrapped Phase values for overlaping connected components as 2pi int cycles
    
    Parameters
    ----------
    unw1 : array
        array containing unwrapped phase in Frame-1 (South)
    concom1 : array
        array containing connected components in Frame-1 (South)
    ix1 : float
        connected component lablein Frame-1 (South), (1...n, 0 is unrealiable)
        
    unw2 : array
        array containing unwrapped phase in Frame-2 (North)
    concom2 : array
        array containing connected components in Frame-2 (North)
    ix2 : float
        connected component label in Frame-2 (North), (1...n, 0 is unrealiable)
    range_correction: bool
        calculate range correction due to small non 2-pi shift caused by ESD different
        between frames

    Returns
    -------
    diff_value : float
        mean offset between overlapping comoponent in two frames
    correction2pi : float
        2pi interger cycles between overlapping comoponent in two frames
    range_corr : float
        correction for non 2-pi shift between overlapping comoponent in two frames 
    """  
    # Not sure if copy() is needed, left it here to avoid overwriting array due to numpy referencing
    # TODO: use np.where to find the component in the data and keep overlap dimensions for comparison
    data1 = unw1.copy()
    data2 = unw2.copy()
    
    data1[concom1 != ix1] = np.nan
    data2[concom2 != ix2] = np.nan

    diff_value = np.nanmean(data1 - data2)
    std_value = np.nanstd(data1 - data2)
    n_points = np.sum(~np.isnan(data1 - data2))

    # Number of 2pi integer jumps
    num_jump = (np.abs(diff_value) + np.pi) // (2.*np.pi)

    if diff_value < 0:
        num_jump *= -1

    correction2pi = 2.* np.pi *num_jump

    # Get range_correction if selected
    if range_correction:
        range_corr = _range_correction(data1, data2)
    else:
       range_corr = 0 

    if n_points != 0:
        if print_msg:
            print(f'  Frame-1 component: {ix1} - Frame-2 component: {ix2}\n'
                f'   Number of points: {n_points}\n'
                f'   Mean diff: {diff_value:.2f}, std: {std_value:.2f} rad\n'
                f'   Number of 2pi cycles: {num_jump}\n'
                f'   Correction2pi: {correction2pi:.2f}')
            if range_correction: print(f'   Range Correction: {range_corr:.2f}\n',
                                       f'  2piCorr + RangeCorr: {correction2pi + range_corr:.2f}\n')

        return diff_value, correction2pi, range_corr
    else:
        print(f' {np.sum(~np.isnan(data1))}: {np.sum(~np.isnan(data2))}\n'
              f' Number of points: {n_points}')
        return None, None, None

def _range_correction(unw1 : NDArray,
                      unw2 : NDArray,) -> np.float32:
    """
    Calculate range correction due to small non 2-pi shift caused by ESD different
    between frames. If ESD is not used, this correction can be skipped

    Parameters
    ----------
    unw1 : array
        array containing unwrapped phase in Frame-1 (South)
    unw2 : array
        array containing unwrapped phase in Frame-2 (North)
    
    Returns
    -------
    range_corr : float
        correction for non 2-pi shift 
    """

    # Wrap unwrapped Phase in Frame-1 and Frame-2 
    unw1_wrapped = np.mod(unw1,(2*np.pi))-np.pi
    unw2_wrapped = np.mod(unw2,(2*np.pi))-np.pi

    # Get the difference between wrapped images
    arr =unw1_wrapped - unw2_wrapped
    arr -= np.round(arr/(2*np.pi))*2*np.pi
    range_corr = np.angle(np.nanmean(np.exp(1j*arr)))

    return range_corr


def _metadata_offset(unw1 : NDArray, unw2 : NDArray,
                     print_msg : Optional[bool] = True) -> Tuple[np.float32]:
    """
    Get mean difference of metadata layers
    
    Parameters
    ----------
    unw1 : array
        array containing unwrapped phase in Frame-1 (South)
        
    unw2 : array
        array containing unwrapped phase in Frame-2 (North)

    Returns
    -------
    diff_value : float
        mean offset between overlapping area in two frames
    """  
    # Not sure if copy() is needed, left it here to avoid overwriting array due to numpy referencing
    # TODO: use np.where to find the component in the data and keep overlap dimensions for comparison
    data1 = unw1.copy()
    data2 = unw2.copy()
    
    diff_value = np.nanmean(data1 - data2)
    n_points = np.sum(~np.isnan(data1 - data2))

    if n_points != 0:
        if print_msg:
            print(f'   Number of points: {n_points}\n'
                f'   Mean diff: {diff_value:.2f}\n')

        return diff_value
    else:
        print(f' {np.sum(~np.isnan(data1))}: {np.sum(~np.isnan(data2))}\n'
              f' Number of points: {n_points}')
        return None


########################### MAIN #############################################

def product_stitch_sequential(input_unw_files : List[str],
                              input_conncomp_files : List[str],
                              output_unw : Optional[str]  = './unwMerged',
                              output_conn : Optional[str] = './connCompMerged',
                              output_format : Optional[str] = 'ENVI',
                              bounds : Optional[tuple] = None, 
                              clip_json : Optional[str] = None,
                              mask_file : Optional[str] = None,
                              correction_method : Optional[str] = 'cycle2pi', # [meandiff, cycle2pi]
                              range_correction : Optional[bool] = True,
                              verbose : Optional[bool] = False,
                              save_fig : Optional[bool] = False,
                              overwrite : Optional[bool] = True) -> None:
    """
    Sequential stitching of frames along the track. Starts from the Southern frame and goes towards 
    he North. Stitching is perform with forward and backward corrections of overlapping components
    between two neighboring frames. The stitched track is stored locally, and then cropped to meet 
    the defined ARIAtools bounding box and masked with a water-mask if selected.

    Parameters
    ----------
    input_unw_files : list
        list of S1 files with Unwrapped Phase (multiple frames along same track)
        ARIA GUNW: 'NETCDF:"%path/S1-GUNW-*.nc":/science/grids/data/unwrappedPhase'
    input_conncomp_files : list
        list of S1 files with Connected Components of unwrapped phase (multiple frames along same track)
        ARIA GUNW: 'NETCDF:"%path/S1-GUNW-*.nc":/science/grids/data/connectedComponents'
    output_unw : str
        str pointing to path and filename of output stitched Unwrapped Phase
    output_conn : str
        str pointing to path and filename of output stitched Connected Components
    output_format : str
        output format used for gdal writer [Gtiff, ENVI, VRT], default is ENVI
        NOTE: did not test how it works with formats other than default, should be ok
    bounds : tuple
        (West, South, East, North) bounds obtained in ariaExtract.py
    clip_json : str
        path to /productBoundingBox.json producted by ariaExtract.py
    mask_file : str
        path to water mask file, example: %aria_extract_path/mask/watermask.msk.vrt
    correction_method : str
        correction method for overlapping components,available options:
         meanoff - mean offset between components
         cycle2pi - 2pi integer cycles between components 
    range_correction : bool
        use correction for non 2-pi shift in overlapping components [True/False]
    verbose : bool
        print info messages [True/False]
    save_fig : bool
        save figure with stitched outputs [True/False]
    overwrite : bool
        overwrite stitched products [True/False]
                                       
    NOTE: Move cropping (gdal.Warp to bounds and clip_json) and masking to ariaExtract.py to make this function
          modular for other use
    """   

    # Outputs
    output_unw = Path(output_unw).absolute()
    output_conn = Path(output_conn).absolute()
    
    # Get raster attributes [SNWE, latlon_spacing, length, width. nodata]
    # from each input file

    # Initalize variables
    unw_attr_dicts = []  # unwrappedPhase
    conncomp_attr_dicts = [] # connectedComponents
    
    # We assume that the filename order is the same for unwrappedPhase and
    # connectedComponent lists 
    temp_snwe_list = []
    temp_latlon_spacing_list = []

    for unw_file, conn_file in zip(input_unw_files, input_conncomp_files):
        unw_attr_dicts.append(get_GUNW_attr(unw_file))
        conncomp_attr_dicts.append(get_GUNW_attr(conn_file))
        # get frame bounds, assume are the same for unw and conncomp
        temp_snwe_list.append(unw_attr_dicts[-1]['SNWE'])
        temp_latlon_spacing_list.append([unw_attr_dicts[-1]['LAT_SPACING'], 
                                        unw_attr_dicts[-1]['LON_SPACING']])

    # get sorted indices for frame bounds, from South to North
    # Sequential stitching starts from the most south frame and moves 
    # forward to next one in the North direction
    # TODO: add option to reverse direction of stitching
    sorted_ix = np.argsort(np.array(temp_snwe_list)[:,0], axis=0)

    # Loop through attributes
    snwe_list = [temp_snwe_list[ii] for ii in sorted_ix]
    latlon_spacing_list = [temp_latlon_spacing_list[ii] for ii in sorted_ix]
    
    # Loop through sorted frames, and stitch neighboring frames
    for i, (ix1, ix2) in enumerate(zip(sorted_ix[:-1], sorted_ix[1:])):
        if verbose:
            print('Frame-1:', unw_attr_dicts[ix1]['PATH'].split('"')[1].split('/')[-1])
            print('Frame-2:', unw_attr_dicts[ix2]['PATH'].split('"')[1].split('/')[-1])
        # Frame1
        frame1_unw_array = get_GUNW_array(unw_attr_dicts[ix1]['PATH']) 
        frame1_conn_array = get_GUNW_array(conncomp_attr_dicts[ix1]['PATH']) 
        # Frame2
        frame2_unw_array = get_GUNW_array(unw_attr_dicts[ix2]['PATH']) 
        frame2_conn_array = get_GUNW_array(conncomp_attr_dicts[ix2]['PATH']) 

        # Mask nodata values
        frame1_unw_array[frame1_unw_array == unw_attr_dicts[ix1]['NODATA']] =np.nan
        frame1_conn_array[frame1_conn_array == conncomp_attr_dicts[ix1]['NODATA']] =np.nan
        frame2_unw_array[frame2_unw_array == unw_attr_dicts[ix2]['NODATA']] =np.nan
        frame2_conn_array[frame2_conn_array == conncomp_attr_dicts[ix2]['NODATA']] =np.nan

        if i == 0:
            (corr_uw1, corr_unw2,
             corr_conn1, corr_conn2) = stitch_2frames(frame1_unw_array, frame1_conn_array, unw_attr_dicts[ix1],  
                                                      frame2_unw_array, frame2_conn_array, unw_attr_dicts[ix2],
                                                      correction_method=correction_method, range_correction=range_correction,
                                                      verbose=verbose)
            # Store corrected values
            corrected_unw_arrays = [corr_uw1, corr_unw2]
            corrected_conn_arrays = [corr_conn1, corr_conn2]
        else:
            (corr_uw1, corr_unw2,
             corr_conn1, corr_conn2) = stitch_2frames(corrected_unw_arrays[-1], corrected_conn_arrays[-1], unw_attr_dicts[ix1],  
                                                      frame2_unw_array, frame2_conn_array, unw_attr_dicts[ix2], 
                                                      correction_method=correction_method, range_correction=range_correction,
                                                      verbose=verbose)

            # Overwrite the last element in corrected arrays
            # TODO: check how to do this without using del
            del corrected_unw_arrays[-1], corrected_conn_arrays[-1] 
            corrected_unw_arrays.extend([corr_uw1, corr_unw2])
            corrected_conn_arrays.extend([corr_conn1, corr_conn2])      

    # Combine corrected unwrappedPhase and connectedComponents arrays
    combined_unwrap, combined_snwe, _ = combine_data_to_single(corrected_unw_arrays, snwe_list, 
                                                               latlon_spacing_list, method='mean')
    combined_conn, _, _ = combine_data_to_single(corrected_conn_arrays, snwe_list,
                                                 latlon_spacing_list, method='min')

    # replace nan with 0.0
    combined_unwrap = np.nan_to_num(combined_unwrap, nan=0.0)
    combined_conn = np.nan_to_num(combined_conn, nan=-1.0)

    # Write 
    # create temp files
    temp_unw_out = output_unw.parent / ('temp_' + output_unw.name)
    temp_conn_out = output_conn.parent / ('temp_' + output_conn.name)
    # write stitched unwrappedPhase
    write_GUNW_array(temp_unw_out, combined_unwrap, combined_snwe, 
                     format=output_format, verbose=verbose, update_mode=overwrite, 
                     add_vrt=True, nodata=0.0)
    
    # write stitched connectedComponents 
    write_GUNW_array(temp_conn_out, combined_conn, combined_snwe, 
                     format=output_format, verbose=verbose, update_mode=overwrite, 
                     add_vrt=True, nodata=-1.0)

    # Crop 
    [print(f'Cropping to {bounds}') if verbose and bounds else None]
    if overwrite:
        [print(f'Removing {output_unw}, {output_conn}') if verbose else None]
        output_unw.unlink(missing_ok=True)
        output_conn.unlink(missing_ok=True)

    # NOTE: Run gdal.Warp on temp file, if input and output are the same
    #       warp creates empty raster, investigate why
    #       Also, it looks like it is important to close gdal.Warp
    #       gdal.Warp/Translate add 6 seconds to runtime

    for output, input in zip([output_unw, output_conn], [temp_unw_out, temp_conn_out]):
        # Crop if selected 
        ds = gdal.Warp(str(output), 
                       str(input.with_suffix('.vrt')),
                       format=output_format,
                       cutlineDSName=clip_json,
                       outputBounds=bounds,
                       #cropToCutline = True,
                    )
        ds = None
        # Update VRT
        [print(f'Writing {output}, {output.with_suffix(".vrt")}') if verbose else None]
        gdal.Translate(str(output.with_suffix('.vrt')), str(output), format="VRT")
        # Remove temp files
        [ii.unlink() for ii in [input, input.with_suffix('.vrt'), 
                                input.with_suffix('.hdr'), input.with_suffix('.aux.xml')]]

        # Mask
        if mask_file:
            [print(f'Mask {output} with {mask_file}') if verbose else None]

            mask_array = get_GUNW_array(mask_file).astype(np.float32)
            array = get_GUNW_array(str(output.with_suffix('.vrt')))

            if output == output_conn:
                # Mask connected components
                mask_array[mask_array==0.0] = np.nan
                array[array==-1.0] = np.nan
                update_array = mask_array * array
                update_array = np.nan_to_num(update_array, nan=-1.0)
                #update_array[np.isnan(update_array)] = -1.0
            else:
                update_array = mask_array * array

            update_file=gdal.Open(str(output), gdal.GA_Update)
            update_file=update_file.GetRasterBand(1).WriteArray(update_array)
            update_file=None

    # Plot stitched
    # NOTE: saving output figure adds 4 seconds 
    if save_fig:
        plot_GUNW_stitched(str(output_unw.with_suffix('.vrt')), str(output_conn.with_suffix('.vrt')))                

    return combined_unwrap, combined_conn, combined_snwe


def product_stitch_sequential_metadata(input_unw_files : List[str],
                              output_unw : Optional[str]  = './unwMerged',
                              output_format : Optional[str] = 'ENVI',
                              verbose : Optional[bool] = False) -> None:
    """
    Sequential stitching of frames implementation from `product_stitch_sequential`
    for metadata layers

    Parameters
    ----------
    input_unw_files : list
        list of S1 files with Unwrapped Phase (multiple frames along same track)
        ARIA GUNW: 'NETCDF:"%path/S1-GUNW-*.nc":/science/grids/data/unwrappedPhase'
    output_unw : str
        str pointing to path and filename of output stitched Unwrapped Phase
    output_format : str
        output format used for gdal writer [Gtiff, ENVI, VRT], default is ENVI
        NOTE: did not test how it works with formats other than default, should be ok
    verbose : bool
        print info messages [True/False]
                                       
    NOTE: Move cropping (gdal.Warp to bounds and clip_json) and masking to ariaExtract.py to make this function
          modular for other use
    """

    # Create VRT and exit early if only one frame passed,
    # and therefore no stitching needed
    if len(input_unw_files) == 1:
        gdal.BuildVRT(output_unw+'.vrt', input_unw_files)
        return

    # Outputs
    output_unw = Path(output_unw).absolute()
    
    # Get raster attributes [SNWE, latlon_spacing, length, width. nodata]
    # from each input file

    # Initalize variables
    unw_attr_dicts = []  # metadata layers
    temp_snwe_list = []
    temp_latlon_spacing_list = []

    for unw_file in input_unw_files:
        unw_attr_dicts.append(get_GUNW_attr(unw_file))
        temp_snwe_list.append(unw_attr_dicts[-1]['SNWE'])
        temp_latlon_spacing_list.append([unw_attr_dicts[-1]['LAT_SPACING'], 
                                        unw_attr_dicts[-1]['LON_SPACING']])

    # get sorted indices for frame bounds, from South to North
    # Sequential stitching starts from the most south frame and moves 
    # forward to next one in the North direction
    # TODO: add option to reverse direction of stitching
    sorted_ix = np.argsort(np.array(temp_snwe_list)[:,0], axis=0)

    # Loop through attributes
    snwe_list = [temp_snwe_list[ii] for ii in sorted_ix]
    latlon_spacing_list = [temp_latlon_spacing_list[ii] for ii in sorted_ix]
    
    # Loop through sorted frames, and stitch neighboring frames
    for i, (ix1, ix2) in enumerate(zip(sorted_ix[:-1], sorted_ix[1:])):
        if verbose:
            print('Frame-1:', unw_attr_dicts[ix1]['PATH'].split('"')[1].split('/')[-1])
            print('Frame-2:', unw_attr_dicts[ix2]['PATH'].split('"')[1].split('/')[-1])
        # Frame1
        frame1_unw_array = get_GUNW_array(unw_attr_dicts[ix1]['PATH']) 
        # Frame2
        frame2_unw_array = get_GUNW_array(unw_attr_dicts[ix2]['PATH']) 

        # Mask nodata values
        frame1_unw_array[frame1_unw_array == unw_attr_dicts[ix1]['NODATA']] =np.nan
        frame2_unw_array[frame2_unw_array == unw_attr_dicts[ix2]['NODATA']] =np.nan

        if i == 0:
            (corr_uw1, corr_unw2) = stitch_2frames_metadata(frame1_unw_array, unw_attr_dicts[ix1],
                                                      frame2_unw_array, unw_attr_dicts[ix2],
                                                      verbose=verbose)
            # Store corrected values
            corrected_unw_arrays = [corr_uw1, corr_unw2]
        else:
            (corr_uw1, corr_unw2) = stitch_2frames_metadata(corrected_unw_arrays[-1], unw_attr_dicts[ix1],
                                                      frame2_unw_array, unw_attr_dicts[ix2],
                                                      verbose=verbose)

            # Overwrite the last element in corrected arrays
            # TODO: check how to do this without using del
            del corrected_unw_arrays[-1]
            corrected_unw_arrays.extend([corr_uw1, corr_unw2])

    # Combine corrected unwrappedPhase arrays
    combined_unwrap, combined_snwe, _ = combine_data_to_single(corrected_unw_arrays, snwe_list, 
                                                               latlon_spacing_list, method='mean',
                                                               latlon_step=[-0.1,0.1])

    # replace nan with 0.0
    combined_unwrap = np.nan_to_num(combined_unwrap, nan=0.0)

    # Write 
    # create temp files
    unw_out = output_unw.parent / (output_unw.name)
    # write stitched unwrappedPhase
    write_GUNW_array(unw_out, combined_unwrap, combined_snwe, 
                     format=output_format, verbose=verbose, 
                     add_vrt=True, nodata=0.0)


    return


############################################################################
#################### READ/WRITE GDAL UTILITIES #############################
def get_GUNW_attr(filename : Union[str, Path]) -> dict:
    """
    Use GDAl to get raster metadata
    
    Parameters
    ----------
    filename : str
        path to raster 

    Returns
    -------
    raster_attr : dict
        raster attribute dict 
        [path, nodata, length, width, snwe, lon_spacing, lat_spacing, projection]
    """

    # Use GDAL to read GUNW netcdf
    ds = gdal.Open(filename, gdal.GA_ReadOnly)

    # Get GUNW Raster attributes
    nodata = ds.GetRasterBand(1).GetNoDataValue()

    # Get raster geographical information
    transform = ds.GetGeoTransform()
    xsize = ds.RasterXSize
    ysize = ds.RasterYSize
    snwe = [transform[3] + ysize * transform[5], transform[3],
            transform[0], transform[0] + xsize * transform[1]]
    lon_spacing = transform[1]
    lat_spacing = transform[5]
    
    projection = ds.GetProjection()

    # wrap raster info in dict
    raster_attr = {'PATH'   : filename,
                   'NODATA' : nodata,
                   'LENGTH' : ysize,
                   'WIDTH'  : xsize,
                   'SNWE'   : snwe,
                   'LON_SPACING' : lon_spacing,
                   'LAT_SPACING' : lat_spacing,
                   'PROJECTION' : projection}

    # close
    ds = None

    return raster_attr

def get_GUNW_array(filename : Union[str, Path],
                   subset : Optional[tuple] = None) -> NDArray:
    """
    Use GDAl to get raster data [as array]

    Parameters
    ----------
    filename : str
        path to raster
    subset : slice
        subset created with np.s_ TODO: insert tuple of (x1,x2, y1,y2)
        and then convert to slice with np.s_

    Returns
    -------
    data : array
        raster data 2D array [length x width]
    """
    
    # Use GDAL to read GUNW netcdf
    ds = gdal.Open(filename, gdal.GA_ReadOnly)

    data =  ds.ReadAsArray()
    # close
    ds = None

    # Subset array
    if subset:
        data = data[subset]

    return data

def write_GUNW_array(output_filename: Union[str, Path], 
                     array:np.ndarray, 
                     snwe:list,
                     nodata: Optional[str] = 'NAN',
                     format : Optional[str] ='ENVI',
                     epsg : Optional[int]= 4326,
                     add_vrt : Optional[bool] = True, 
                     verbose : Optional[bool] = False,
                     update_mode : Optional[bool] = True) -> None:
    """
    Use GDAl to write raster 

    Parameters
    ----------
    output_filename : str
        path to raster
    array : array
        numpy array of raster to be written
    snwe  : list
        (South, North, West, East) bounds of raster to be written
    nodata : str
        value or nan for NODATA (used for VRT creation)
    format : str
        output raster format, default is ENVI
    epsg : int
        projection epsg, default is 4326 for WGS84
    add_vrt : bool
        flag to create VRT for output raster [True/False] 
    verbose : bool
        print info messages [True/False]
    update_mode : bool
        flag to overwrite the existing files [True/False]
    """

    array_type = gdal_array.NumericTypeCodeToGDALTypeCode(array.dtype)

    # Output path
    output = Path(output_filename).absolute()
    output_vrt = output.with_suffix('.vrt') 

    if update_mode:
        [print(f'Remove {output}') if verbose else None]
        output.unlink(missing_ok=True)
        if add_vrt:
            [print(f'Remove {output_vrt}') if verbose else None]
            output_vrt.unlink(missing_ok=True)
    
    # Get lat, lon pixel spacing
    num_bands = 1
    if len(array.shape) > 2:
        num_bands = array.shape[0]
        x_step = (snwe[3] - snwe[2]) / array.shape[2]
        y_step = (snwe[0] - snwe[1]) / array.shape[1]
    else:
        x_step = (snwe[3] - snwe[2]) / array.shape[1]
        y_step = (snwe[0] - snwe[1]) / array.shape[0]

    # Geotransform  
    geo = (snwe[2], x_step, 0, snwe[1], 0, y_step)   
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg) # get projection
    
    # Write
    driver = gdal.GetDriverByName(format)
    if len(array.shape) > 2:
        out_ds = driver.Create(str(output),
                           array.shape[2],
                           array.shape[1], 
                           array.shape[0], 
                           array_type)
    else:
        out_ds = driver.Create(str(output), 
                           array.shape[1], 
                           array.shape[0], 
                           1, 
                           array_type)
    
    out_ds.SetProjection(srs.ExportToWkt())
    out_ds.SetGeoTransform(geo)

    if verbose:
        print(f'Writing {output}')

    for i in range(num_bands):
        band = out_ds.GetRasterBand(i+1)
        if num_bands > 1:
            band.WriteArray(array[i])
        else:
            band.WriteArray(array)
        band.FlushCache()
        band.ComputeStatistics(False)

    # Close
    out_ds = None

    if add_vrt:
        # Build virtual VRT  
        vrt = gdal.BuildVRT(str(output_vrt), str(output), srcNodata=nodata)
        vrt.FlushCache()
        vrt = None      

# Extract overlap bounds
def frame_overlap(snwe1 : list,
                  snwe2 : list,
                  latlon_step1 : list,
                  latlon_step2 : list,
                  latlon_step: Optional[list] = [-0.000833334, 0.000833334] \
                  ) -> Tuple[tuple, tuple]:
    """
    Parameters
    ----------
    Use raster metadata to find overlap between two Images
    snwe1 : list
        [South, North, West, East] bounds of Image-1
    snwe2 : list
        [South, North, West, East] bounds of Image-2
    latlon_step1 : list
        latitude and longitude pixel spacing of Image-1
    latlon_step2 : list
        latitude and longitude pixel spacing of Image-2

    Returns
    -------
    subset1 : np.slice
        Intersection subset [y1:y2, x1:x2] for the Image-1
    subset2 : np.slice
        Intersection subset [y1:y2, x1:x2] for the Image-2
    """ 

    snwe = np.vstack([snwe1, snwe2])
    #Find overlap bounds
    overlap_snwe = [np.max(snwe[:,0]), np.min(snwe[:,1]), 
                    np.max(snwe[:,2]), np.min(snwe[:,3])]
    
    # Georeferenced space to image coordinate space
    x1, y1 = lalo2xy(overlap_snwe[1], overlap_snwe[2], snwe1, latlon_step1) # Frame-1 
    x2, y2 = lalo2xy(overlap_snwe[1], overlap_snwe[2], snwe2, latlon_step2) # Frame-2

    # Overlap bounds - force overlaps to have same dimensions
    # latlon_spacing sometimes diff at 13th decimal
    length = int(round((overlap_snwe[0] - overlap_snwe[1]) / latlon_step[0])) 
    width  = int(round((overlap_snwe[3] - overlap_snwe[2]) / latlon_step[1]))

    subset1 = np.s_[y1:y1 + length, x1:x1 + width]
    subset2 = np.s_[y2:y2 + length, x2:x2 + width]

    return subset1, subset2

def lalo2xy(lat : np.float32, 
            lon: np.float32, 
            data_snwe : list, 
            latlon_step :list, 
            rounding_method : Optional[str] = 'floor') -> Tuple[np.int16, np.int16]:
    """
    Georeferenced coordinates to image space coordinates. GDAL raster starting point
    is the upper left corner.

    Parameters
    ----------
    lat : float
        search latitude 
    lon : float
        search longitude
    data_snwe : list
        [South, North, West, East] bounds of raster
        North (y0) and West(x0) as reference point
    latlon_step : list
        pixel spacing [latitude_spacing, longitude_spacing]                              
    rounding_method : str
        rounding method, default is 'floor' other option is 'around'
        Read notes below. TODO. test different rounding routines

    Returns
    -------
    x : int
        coordinate in image space (x-axis/columns, direction of width) from x0 (top up column)
    y : int
        coordinate in image space (y-axis/rows, direction of length) from y0 (top left row)
    """

    # np.floor works better with points and raster - Need to check why
    # but with two rasters sometimes one pixel is missing or is redundant
    if rounding_method == 'floor':
        x = int(np.floor((lon - data_snwe[2]) / latlon_step[1] + 0.01))
        y = int(np.floor((lat - data_snwe[1]) / latlon_step[0] + 0.01))
    
    # np.around works better with two rasters
    # test it out, I think it has something to how numpy floor is rounding negative values
    # example np.around(-125.2) = -125 np.floor(-125.2) = -126
    # np.around(125.6) = 126, np.floor(125.6) = 125
    elif rounding_method == 'around': 
        x = int(np.around((lon - data_snwe[2]) / latlon_step[1] + 0.01))
        y = int(np.around((lat - data_snwe[1]) / latlon_step[0] + 0.01))

    return x, y

def combine_data_to_single(data_list : list, 
                           snwe_list : list, 
                           latlon_step_list : list,
                           method : Optional[str] = 'mean',
                           latlon_step: Optional[list] = [-0.000833334, 0.000833334] \
                           ) -> Tuple[NDArray, NDArray, list]:
    """
    Merge multiple arrays to one array. Combine them in ndarray, then apply function along the n_layers axis

    Parameters
    ----------
    data_list : list
        list of arrays containing raster values
    snwe_list : list
        list of arrays containing snwe (extent) values
    latlon_step_list : list
        list of arrays containing pixel spacing in lat and lon direction for each dataset
    method : str
        method to merge overlapping pixes, use mean, min, max erc..
        TODO: need to refine this part of code
                                        
    Returns
    -------
    comb_data : ndarray 
        combined data [n_frames, length, width] 
    SNWE : array
        extent of the combined data
    latlon_step : array
        pixel spacing in lat, lon of combined data
    
    """
    # Get the maximum extent of all data
    n = len(data_list)
    snwe_all = np.squeeze([snwe for snwe in snwe_list])
    
    SNWE = np.array([np.min(snwe_all[:,0]), np.max(snwe_all[:,1]),
                     np.min(snwe_all[:,2]), np.max(snwe_all[:,3])]).T
           
    length = abs(int(np.around((SNWE[1] - SNWE[0]) / latlon_step[0] + 0.01)))
    width = abs(int(np.around((SNWE[2] - SNWE[3]) / latlon_step[1] + 0.01)))
    
    # create combined data array
    # handle if 3D metadata layer
    if len(data_list[0].shape) > 2:
        comb_data = np.empty((n, data_list[0].shape[0], length, width), dtype=np.float64) * np.nan
    else:
        comb_data = np.empty((n, length, width), dtype=np.float64) * np.nan
    for i, data in enumerate(data_list):
        x, y = np.abs(lalo2xy(SNWE[1], SNWE[2], snwe_list[i], latlon_step_list[i], 'around'))
        # handle if 3D metadata layer
        if len(data.shape) > 2:
            comb_data[i, 0:data.shape[0], y:y+data.shape[1], x: x+data.shape[2]] = data 
        else:
            comb_data[i, y:y+data.shape[0], x: x+data.shape[1]] = data 

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        # combine using numpy
        if method =='mean':
            comb_data = np.nanmean(comb_data, axis=0)
        elif method =='median': 
            comb_data = np.nanmedian(comb_data, axis=0)
        elif method =='min': 
            comb_data = np.nanmin(comb_data, axis=0)
        elif method =='max': 
            comb_data = np.nanmax(comb_data, axis=0)    

    return comb_data, SNWE,  latlon_step

#################### PLOTTING UTILITIES ########################################

def snwe_to_extent(snwe: list) -> list:
    '''
    Convert SNWE to extent for matplotlib plotting
    '''
    extent = [snwe[2], snwe[3], snwe[0], snwe[1]]
    
    return extent

def plot_GUNW_stitched(stiched_unw_filename: str,
                       stiched_conn_filename: str) -> None:
    '''
    Plotting function for stitched outputs
    '''
    # no display
    mpl.use('Agg')
    # Save plot
    cmap = plt.cm.cividis_r  # define the colormap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be red
    cmaplist[0] = (.9, .1, .1, 1.0)

    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list(
    'Custom cmap', cmaplist, cmap.N)
    #bounds = np.linspace(0, 256, 257)

    output_dir = Path(stiched_unw_filename).absolute()
    output_fig = output_dir.parent / 'stitched.png'
    output_fig.unlink(missing_ok=True) 

    # Load Data
    stitched_unw = get_GUNW_array(stiched_unw_filename)
    stitched_conn = get_GUNW_array(stiched_conn_filename)
    stitched_attr = get_GUNW_attr(stiched_unw_filename)

    # ConnComp discrete colormap
    #bounds = np.linspace(0, int(np.nanmax(stitched_conn)), int(np.nanmax(stitched_conn)+1))
    bounds = np.linspace(0, 30, 31)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    # Mask
    stitched_unw[stitched_unw==0.0] = np.nan
    stitched_conn[stitched_conn==-1.0] = np.nan

    # Common plot options
    plot_kwargs = {'extent' : snwe_to_extent(stitched_attr['SNWE']),
                   'interpolation' : 'nearest'}

    # Figure
    fig, axs = plt.subplots(1,3, dpi=300, sharey=True)
    im1=axs[0].imshow(np.mod(stitched_unw, 4*np.pi), cmap='jet', **plot_kwargs) # Re-wrapped
    im2=axs[1].imshow(stitched_unw * (0.0556 / (6*np.pi)), cmap='jet', 
                      clim=[-0.2, 0.2], **plot_kwargs) # Unwrapped
    im3=axs[2].imshow(stitched_conn, cmap=cmap, norm=norm, **plot_kwargs) # Connected Components

    for im, ax, label, txt, in zip([im1, im2, im3],
                                    axs, 
                                    ['rad','m','#'],
                                    ['Re-wrapped phase w 20 rad', 'Unwrapped phase [m]', 'Connected Components']):
        fig.colorbar(im, ax=ax, location='bottom', shrink=0.8, label=label)
        ax.set_title(txt, fontsize=10)

    fig.tight_layout()
    fig.savefig(str(output_fig))
