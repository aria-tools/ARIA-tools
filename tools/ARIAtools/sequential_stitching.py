#!/usr/bin/env python3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Marin Govorcin
# Copyright 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import numpy as np
from numpy.typing import NDArray

from typing import List, Optional, Tuple
from osgeo import gdal, osr
from pathlib import Path


# import util modules
from ARIAtools.stitiching_util import (get_GUNW_array, get_GUNW_attr,
                                       frame_overlap, combine_data_to_single,
                                       write_GUNW_array, snwe_to_extent,
                                       _nan_filled_array)


from matplotlib import pyplot as plt
import matplotlib as mpl

'''
NOTE:
Sequential stitcher relies on connected components [region of pixels] in the
overlap between two frames with assumption that each component is unwrapped
correctly by SNAPHU. This might not always be the case. If there are
unwrapping errors within the overlapping component, they will propage to the
stitched image. Potential fix could be to manually create connected component
around the unwrapping error and try to re-run stitching.

Connected Component 0 represents masked pixels by SNAPHU. Stitcher merges 0
overlapping components by default for the sake of vizualization. However,
these unwrapped phase pixels are unreliable and will often be misaligned as
2-pi integer cycles shift does not apply here. User is advised to mask these
pixels for further processing.

TODO: Re-enumeration of connected components in the stitched image. Due the
merging of overlapping components, there are gaps in enumeration of connected
component labels. This should not affect the futher processing, but for the
sake of consistency, add function to re-enumerate components.

DISCLAIMER : This is development script. Requires some additional clean-up
and restructuring
'''

#  STITCHING SUBROUTINES


def stitch_unwrapped_frames(input_unw_files: List[str],
                            input_conncomp_files: List[str],
                            proj: Optional[str] = 'EPSG:4326',
                            correction_method: Optional[str] = 'cycle2pi',
                            range_correction: Optional[bool] = True,
                            direction_N_S: Optional[bool] = False,
                            verbose: Optional[bool] = False):

    # Get raster attributes [SNWE, latlon_spacing, length, width, nodata]
    # from each input file

    # Initalize variables
    unw_attr_dicts = []  # unwrappedPhase
    conncomp_attr_dicts = []  # connectedComponents

    # We assume that the filename order is the same for unwrappedPhase and
    # connectedComponent lists
    temp_snwe_list = []

    for unw_file, conn_file in zip(input_unw_files, input_conncomp_files):
        unw_attr_dicts.append(get_GUNW_attr(unw_file, proj=proj))
        conncomp_attr_dicts.append(get_GUNW_attr(conn_file, proj=proj))
        # get frame bounds, assume are the same for unw and conncomp
        temp_snwe_list.append(unw_attr_dicts[-1]['SNWE'])

    # get sorted indices for frame bounds, from South to North
    sorted_ix = np.argsort(np.array(temp_snwe_list)[:, 0], axis=0)
    # reverse direction of stitching, move from North to South
    # NOTE: should not make any difference
    if direction_N_S:
        sorted_ix = sorted_ix[::-1]

    # Wrap stitching parameters into dict
    stitching_dict = dict(correction_method=correction_method,
                          range_correction=range_correction,
                          verbose=verbose)

    # Loop through sorted frames, and stitch neighboring frames
    for i, (ix1, ix2) in enumerate(zip(sorted_ix[:-1], sorted_ix[1:])):
        if verbose:
            print(50*'*')
            print(
                'Frame-1:',
                unw_attr_dicts[ix1]['PATH'].split('"')[1].split('/')[-1])
            print(
                'Frame-2:',
                unw_attr_dicts[ix2]['PATH'].split('"')[1].split('/')[-1])

        # Get numpy masked arrays
        # Frame1
        frame1_unw_array = get_GUNW_array(
            filename=unw_attr_dicts[ix1]['PATH'],
            proj=proj,
            nodata=unw_attr_dicts[ix1]['NODATA']
        )
        frame1_conn_array = get_GUNW_array(
            filename=conncomp_attr_dicts[ix1]['PATH'],
            proj=proj,
            nodata=conncomp_attr_dicts[ix1]['NODATA']
        )
        # Frame2
        frame2_unw_array = get_GUNW_array(
            filename=unw_attr_dicts[ix2]['PATH'],
            proj=proj,
            nodata=unw_attr_dicts[ix2]['NODATA']
        )
        frame2_conn_array = get_GUNW_array(
            filename=conncomp_attr_dicts[ix2]['PATH'],
            proj=proj,
            nodata=conncomp_attr_dicts[ix2]['NODATA']
        )
        # capture all lyr nodata values
        conncomp_nodata_values = [0, -1, conncomp_attr_dicts[ix1]['NODATA'],
                                                      conncomp_attr_dicts[ix2]['NODATA']]
        conncomp_nodata_values = list(set(conncomp_nodata_values))
        unw_nodata_values = [0, unw_attr_dicts[ix1]['NODATA'],
                                             unw_attr_dicts[ix2]['NODATA']]
        unw_nodata_values = list(set(unw_nodata_values))
        # create mask field for all arrays
        conncomp_mask_ix1 = np.isin(frame1_conn_array, conncomp_nodata_values)
        conncomp_mask_ix2 = np.isin(frame2_conn_array, conncomp_nodata_values)
        unw_mask_ix1 = np.isin(frame1_unw_array, conncomp_nodata_values)
        unw_mask_ix2 = np.isin(frame2_unw_array, conncomp_nodata_values)
        # apply mask
        frame1_unw_array = np.ma.masked_where(unw_mask_ix1,
            frame1_unw_array)
        frame1_conn_array = np.ma.masked_where(conncomp_mask_ix1,
            frame1_conn_array)
        frame2_unw_array = np.ma.masked_where(unw_mask_ix2,
            frame2_unw_array)
        frame2_conn_array = np.ma.masked_where(conncomp_mask_ix2,
            frame2_conn_array)

        if i == 0:
            (corr_unw, corr_conn, corr_dict) = stitch_unw2frames(
                frame1_unw_array,
                frame1_conn_array,
                unw_attr_dicts[ix1],
                frame2_unw_array,
                frame2_conn_array,
                unw_attr_dicts[ix2],
                **stitching_dict)
        else:
            (corr_unw, corr_conn, corr_dict) = stitch_unw2frames(
                corr_unw,
                corr_conn,
                corr_dict,
                frame2_unw_array,
                frame2_conn_array,
                unw_attr_dicts[ix2],
                **stitching_dict)

    # replace nan with 0.0
    corr_unw = np.nan_to_num(corr_unw.data, nan=0.0)
    corr_conn = np.nan_to_num(corr_conn.data, nan=-1.0)

    return corr_unw, corr_conn, corr_dict['SNWE']


def stitch_unw2frames(unw_data1: NDArray, conn_data1: NDArray, rdict1: dict,
                      unw_data2: NDArray, conn_data2: NDArray, rdict2: dict,
                      correction_method: Optional[str] = 'cycle2pi',
                      range_correction: Optional[bool] = False,
                      verbose: Optional[bool] = False) -> \
        Tuple[NDArray, NDArray, dict]:
    """
    Function to sequentially stitch two frames along the same track. Mean
    offset or 2pi integer cycles are estimated for each overlapping connected
    components, which after correction are merged into the same component.
    Code runs forward and backward, first correcting Northern frame with
    respect to Southern, and then the opposite if corrected components
    overlap with multiple components in Southern frame.

    Parameters
    ----------
    unw_data1 : array
        array containing unwrapped phase in Frame-1 (South)
    conn_data1 : array
        array containing connected components in Frame-1 (South)
    rdict1 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..]
        for Frame-1 (South)
    unw_data2 : array
        array containing unwrapped phase in Frame-2 (North)
    conn_data2 : array
        array containing connected components in Frame-2 (North)
    rdict2 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..]
        for Frame-2 (North)
    correction_method : str
        method to correct offset between overlapping connected components.
        Available options:
         "meanoff" - mean offset between components
         "cycle2pi" - 2pi integer cycles between components
    verbose : bool
        print info messages [True/False]

    Returns
    -------
    unw_data1 : array
        corrected array containing unwrapped phase in Frame-1 (South)
    conn_data1 : array
        corrected array containing connected components in Frame-1 (South)
    unw_data2 : array
        corrected array containing unwrapped phase in Frame-2 (North)
    conn_data2 : array
        corrected array containing connected components in Frame-2 (North)

    """

    # Adjust connected Component in Frame 2 to start
    # with last component number in Frame-1

    conn_data2 = conn_data2 + np.nanmax(conn_data1)
    conn_data2[conn_data2 == np.nanmax(conn_data1)] = 0.0

    # GET FRAME OVERLAP
    box_1, box_2 = frame_overlap(
        rdict1['SNWE'], rdict2['SNWE'],
        [rdict1['LAT_SPACING'], rdict1['LON_SPACING']],
        [rdict2['LAT_SPACING'], rdict2['LON_SPACING']],
        [rdict1['LAT_SPACING'], rdict1['LON_SPACING']])

    # LOOP OVER COMPONENTS WITHIN THE OVERLAP
    # Get connected component pairs
    if verbose:
        print('\nGetting overlapping components')

    # Forward correction
    conn_pairs = get_overlapping_conn(conn_data1[box_1],
                                      conn_data2[box_2])

    for pair in conn_pairs:
        diff, cycles2pi, range_corr = _integer_2pi_cycles(
            unw_data1[box_1],
            conn_data1[box_1],
            np.float32(pair[1]),
            unw_data2[box_2],
            conn_data2[box_2],
            np.float32(pair[0]),
            range_correction=range_correction,
            print_msg=verbose)

        # Correction methods: mean difference, 2pi integer cycles
        if correction_method == 'cycle2pi':
            correction = cycles2pi
        elif correction_method == 'meanoff':
            correction = diff
            range_correction = False
        else:
            raise ValueError(f'Wrong correction method {correction_method}, ',
                             'Select one of available: "cycle2pi", "meanoff"')

        # add range correction
        if range_correction:
            correction += range_corr
        ik = conn_data2 == np.float32(pair[0])
        unw_data2[ik] += correction
        conn_data2[ik] = np.float32(pair[1])

    # Backward correction
    conn_reverse = get_overlapping_conn(conn_data2[box_2],
                                        conn_data1[box_1])

    # Keep only different componentes in pairing
    ik = np.where(conn_reverse[:, 0] != conn_reverse[:, 1])
    conn_reverse = conn_reverse[ik]

    for pair in conn_reverse:
        print('Going backward!') if verbose else None
        diff, cycles2pi, range_corr = _integer_2pi_cycles(
            unw1=unw_data1[box_1],
            concom1=conn_data1[box_1],
            ix1=np.float32(pair[0]),
            unw2=unw_data2[box_2],
            concom2=conn_data2[box_2],
            ix2=np.float32(pair[1]),
            range_correction=range_correction,
            print_msg=verbose)

        # Correction methods: mean difference, 2pi integer cycles
        if correction_method == 'cycle2pi':
            correction = cycles2pi
        elif correction_method == 'meanoff':
            correction = diff
            range_correction = False
        else:
            raise ValueError(f'Wrong correction method {correction_method}, ',
                             'Select one of available: "cycle2pi", "meanoff"')

        # add range correction
        if range_correction:
            correction += range_corr

        ik = conn_data1 == np.float32(pair[0])
        unw_data1[ik] -= correction
        conn_data1[ik] = np.float32(pair[1])

    # Update connected component frame 2 naming
    idx1 = np.max(conn_data1)
    idx = np.unique(conn_data2[conn_data2 > idx1]).compressed()
    conn_data2 = update_connect_components(conn_data2, idx, idx1+1)

    # Combine corrected unwrappedPhase and connectedComponents arrays
    comb_snwe = [rdict1['SNWE'], rdict2['SNWE']]
    comb_latlon = [[rdict1['LAT_SPACING'], rdict1['LON_SPACING']],
                   [rdict2['LAT_SPACING'], rdict2['LON_SPACING']]]

    (combined_unwrap,
     combined_snwe,
     combined_latlon_spacing) = combine_data_to_single(
        [_nan_filled_array(unw_data1),
         _nan_filled_array(unw_data2)],
        comb_snwe, comb_latlon,
        method='mean',
        latlon_step=comb_latlon[0])
    combined_conn, _, _ = combine_data_to_single(
        [_nan_filled_array(conn_data1),
         _nan_filled_array(conn_data2)],
        comb_snwe, comb_latlon,
        method='min',
        latlon_step=comb_latlon[0])
    # combined dict
    combined_dict = dict(SNWE=combined_snwe,
                         LAT_SPACING=combined_latlon_spacing[0],
                         LON_SPACING=combined_latlon_spacing[1])

    combined_unwrap = np.ma.masked_invalid(combined_unwrap)
    np.ma.set_fill_value(combined_unwrap, 0.)
    combined_conn = np.ma.masked_invalid(combined_conn)
    np.ma.set_fill_value(combined_conn, -1.)

    return combined_unwrap, combined_conn, combined_dict


def get_overlapping_conn(conn1: NDArray,
                         conn2: NDArray) -> Tuple[NDArray, NDArray]:
    """
    Get forward and backward pairs of overlapping connected components with
    the number of overlaping pixels.
    Return [connectComponent-Frame2, connectComponent-Frame1, number of pixels]

    Parameters
    ----------
    conn1 : array
        array containing connected components in Frame-1 (South)
    conn2 : array
        array containing connected components in Frame-2 (North)

    Returns
    -------
    conn_pairs : array
        forward pairs of overlapping components
    conn_pairs_reverse : array
        backward pairs of overlapping components
    """

    conn_union = np.empty((0, 3), dtype=np.int32)

    # Get unique components
    if np.ma.is_masked(conn1):
        concomp1 = np.unique(conn1).compressed()
    else:
        concomp1 = np.unique(conn1)

    if np.ma.is_masked(conn2):
        concomp2 = np.unique(conn2).compressed()
    else:
        concomp2 = np.unique(conn2)

    # Loop through them and connect size and number of overlapping data
    for ix2 in concomp2:
        for ix1 in concomp1:
            # Skip 0 component combination with other components
            if not ix1 == 0 and not ix2 == 0:
                idx = np.where((conn1 == ix1) & (conn2 == ix2))[0]
                if np.count_nonzero(idx) > 0:
                    carray = np.array([ix2, ix1, np.count_nonzero(idx)],
                                      dtype=np.int32, ndmin=2)

                    conn_union = np.concatenate((conn_union, carray), axis=0)

            # Get 0 components in both frames
            elif ix1 == 0 and ix2 == 0:
                idx = np.where((conn1 == ix2) & (conn2 == ix1))[0]
                if np.count_nonzero(idx) > 0:
                    carray = np.array([ix2, ix1, np.count_nonzero(idx)],
                                      dtype=np.int32, ndmin=2)

                    conn_union = np.concatenate((conn_union, carray), axis=0)

    # Find components to correct in Frame 2
    conn_pairs = np.empty((0, 3), dtype=np.int32)

    for k in np.unique(conn_union[:, 0]):
        ik = conn_union[:, 0] == k
        # find number of times components is referenced
        count = np.sum(conn_union[:, 0] == k)

        if count > 1:
            max_points = np.max(conn_union[ik][:, 2])
            # Select the one with the most points
            ik = np.where((conn_union[:, 0] == k) &
                          (conn_union[:, 2] == max_points))[0]
            # Select first if there are more pairs with same num of points
            ik = np.array(ik[0], ndmin=1) if ik.shape[0] > 1 else ik

        conn_pairs = np.concatenate((conn_pairs, conn_union[ik]), axis=0)

    return conn_pairs


def update_connect_components(conncomp_array: NDArray,
                              unique_components: NDArray,
                              renumber_from: int = 21) -> NDArray:
    c = conncomp_array.copy()
    for ix, component in enumerate(unique_components):
        c[conncomp_array == component] = renumber_from + ix
    return conncomp_array


def _integer_2pi_cycles(unw1: NDArray, concom1: NDArray, ix1: np.float32,
                        unw2: NDArray, concom2: NDArray, ix2: np.float32,
                        range_correction: Optional[bool] = False,
                        print_msg: Optional[bool] = False) -> \
        Tuple[np.float32, np.float32, np.float32]:
    """
    Get mean difference of unwrapped Phase values for overlapping
    connected components as 2pi int cycles

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
        connected component label in Frame-2 (North),
        (1...n, 0 is unrealiable)
    range_correction: bool
        calculate range correction due to small non 2-pi shift caused by
        ESD different between frames

    Returns
    -------
    diff_value : float
        mean offset between overlapping component in two frames
    correction2pi : float
        2pi interger cycles between overlapping component in two frames
    range_corr : float
        correction for non 2-pi shift between overlapping component
        in two frames
    """
    # find the component in the data and keep overlap
    # dimensions for comparison
    idx = np.where((concom1 == ix1) & (concom2 == ix2))

    diff = unw1[idx] - unw2[idx]

    # Masked array to array
    if np.ma.isMaskedArray(diff):
        diff = diff.data

    median_diff = np.nanmedian(diff)
    std_value = np.nanstd(diff)
    n_points = np.count_nonzero(diff)

    # Number of 2pi integer jumps
    num_jump = (np.abs(median_diff) + np.pi) // (2.*np.pi)

    if median_diff < 0:
        num_jump *= -1

    correction2pi = 2. * np.pi * num_jump

    # Get range_correction if selected
    if range_correction:
        range_corr = _range_correction(unw1[idx], unw2[idx])
    else:
        range_corr = 0

    # Note: range correctio sometimes gives oposite sign of
    #        correction, and add half or one cycle more.
    #        not sure, why that happens?? below is a hardcoded solution
    if np.abs(median_diff - (correction2pi + range_corr)) > 3.14:
        range_corr *= -1

    if print_msg:
        print(
            f' Frame-1 component: {ix1} - Frame-2 component: {ix2}\n'
            f'   Number of points: {n_points}\n'
            f'   Median diff: {median_diff:.2f}, std: {std_value:.2f} rad\n'
            f'   Number of 2pi cycles: {num_jump}\n'
            f'   Correction2pi: {correction2pi:.2f}')
        if range_correction:
            print(
                f'   Range Corr: {range_corr:.2f} \n',
                f' 2piCorr + RangeCorr: {correction2pi + range_corr:.2f}\n')

    return median_diff, correction2pi, range_corr


def _range_correction(unw1: NDArray,
                      unw2: NDArray) -> np.float32:
    """
    Calculate range correction due to small non 2-pi shift caused by ESD
    different between frames. If ESD is not used, this correction
    can be skipped

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
    unw1_wrapped = np.mod(unw1, (2*np.pi)) - np.pi
    unw2_wrapped = np.mod(unw2, (2*np.pi)) - np.pi

    # Get the difference between wrapped images
    arr = unw1_wrapped - unw2_wrapped
    arr -= np.round(arr/(2*np.pi))*2*np.pi
    range_corr = np.angle(np.nanmean(np.exp(1j*arr)))

    return range_corr


# Stitching routines for another layers

def stitch_2frames_metadata(
        unw_data1: NDArray, rdict1: dict,
        unw_data2: NDArray, rdict2: dict,
        verbose: Optional[bool] = False) -> Tuple[NDArray, NDArray]:
    """Sequential stitching function implementation from `stitch_2frames`
    for metadata layers

    Parameters
    ----------
    unw_data1 : array
        array containing unwrapped phase in Frame-1 (South)
    rdict1 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..]
        for Frame-1 (South)
    unw_data2 : array
        array containing unwrapped phase in Frame-2 (North)
    rdict2 : dict
        dict containing raster metadata [SNWE, latlon_spacing, nodata etc..]
        for Frame-2 (North)
    verbose : bool
        print info messages [True/False]

    Returns
    -------
    unw_data1 : array
        corrected array containing unwrapped phase in Frame-1 (South)
    unw_data2 : array
        corrected array containing unwrapped phase in Frame-2 (North)

    TODO: combine corrected array in this function, return only the
          stitched unwrapped Phase
    """

    # GET FRAME OVERLAP
    box_1, box_2 = frame_overlap(
        rdict1['SNWE'], rdict2['SNWE'],
        [rdict1['LAT_SPACING'], rdict1['LON_SPACING']],
        [rdict2['LAT_SPACING'], rdict2['LON_SPACING']],
        [-0.1, 0.1])

    # EXAMINE THE OVERLAP
    for i in range(unw_data1.shape[0]):
        diff = _metadata_offset(unw_data1[i][box_1],
                                unw_data2[i][box_2],
                                print_msg=verbose)

        unw_data2[i] += diff

    return unw_data1, unw_data2


def _metadata_offset(unw1: NDArray, unw2: NDArray,
                     print_msg: Optional[bool] = False) -> Tuple[np.float32]:
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
    # Not sure if copy() is needed, left it here to avoid overwriting array
    # due to numpy referencing
    # TODO: use np.where to find the component in the data and keep overlap
    # dimensions for comparison
    data1 = unw1.copy()
    data2 = unw2.copy()

    diff_value = np.nanmean(data1 - data2)
    n_points = np.sum(~np.isnan(data1 - data2))

    if n_points != 0:
        if print_msg:
            print(
                f'   Number of points: {n_points}\n'
                f'   Mean diff: {diff_value:.2f}\n')

        return diff_value
    else:
        print(f' {np.sum(~np.isnan(data1))}: {np.sum(~np.isnan(data2))}\n'
              f' Number of points: {n_points}')
        return None

# MAIN


def product_stitch_sequential(input_unw_files: List[str],
                              input_conncomp_files: List[str],
                              arrres: List[float],
                              epsg: Optional[str] = '4326',
                              output_unw: Optional[str] = './unwMerged',
                              output_conn: Optional[str] = './connCompMerged',
                              output_format: Optional[str] = 'ENVI',
                              bounds: Optional[tuple] = None,
                              clip_json: Optional[str] = None,
                              mask_file: Optional[str] = None,
                              # [meandiff, cycle2pi]
                              correction_method: Optional[str] = 'cycle2pi',
                              range_correction: Optional[bool] = True,
                              verbose: Optional[bool] = False,
                              save_fig: Optional[bool] = False,
                              overwrite: Optional[bool] = True) -> None:
    """
    Sequential stitching of frames along the track. Starts from the Southern
    frame and goes towards the North. Stitching is perform with forward and
    backward corrections of overlapping components between two neighboring
    frames. The stitched track is stored locally, and then cropped to meet
    the defined ARIAtools bounding box and masked with a water-mask
    if selected.

    Parameters
    ----------
    input_unw_files : list
        list of S1 files with Unw Phase (multiple frames along same track)
        ARIA GUNW:
        'NETCDF:"%path/S1-GUNW-*.nc":/science/grids/data/unwrappedPhase'
    input_conncomp_files : list
        list of S1 files with Connected Components of unwrapped phase
        (multiple frames along same track)
        ARIA GUNW:
        'NETCDF:"%path/S1-GUNW-*.nc":/science/grids/data/connectedComponents'
    output_unw : str
        str pointing to path and filename of output stitched Unwrapped Phase
    output_conn : str
        str pointing to path and filename of output stitched
        Connected Components
    output_format : str
        output format used for gdal writer [e.g., Gtiff ENVI], default is ENVI
    bounds : tuple
        (West, South, East, North) bounds obtained in ariaExtract.py
    clip_json : str
        path to /productBoundingBox.json producted by ariaExtract.py
    mask_file : str
        path to water mask file, example:
        %aria_extract_path/mask/watermask.msk.vrt
    correction_method : str
        correction method for overlapping components, available options:
         meanoff - mean offset between components
         cycle2pi - 2pi integer cycles between components
    range_correction : bool
        use correction for non 2-pi shift in overlapping components
        [True/False]
    verbose : bool
        print info messages [True/False]
    save_fig : bool
        save figure with stitched outputs [True/False]
    overwrite : bool
        overwrite stitched products [True/False]

    NOTE: Move cropping (gdal.Warp to bounds and clip_json) and masking
          to ariaExtract.py to make this function modular for other use
    """

    # Outputs
    output_unw = Path(output_unw).absolute()
    if not output_unw.parent.exists():
        output_unw.parent.mkdir()
    output_conn = Path(output_conn).absolute()
    if not output_conn.parent.exists():
        output_conn.parent.mkdir()
    # create temp files
    temp_unw_out = output_unw.parent / ('temp_' + output_unw.name)
    temp_conn_out = output_conn.parent / ('temp_' + output_conn.name)
    # obtain reference epsg code to assign to intermediate outputs
    ref_proj_str = get_GUNW_attr(input_unw_files[0])['PROJECTION']
    srs = osr.SpatialReference(wkt=ref_proj_str)
    srs.AutoIdentifyEPSG()
    ref_proj = srs.GetAuthorityCode(None)

    # Create VRT and exit early if only one frame passed,
    # and therefore no stitching needed
    if len(input_unw_files) == 1:
        gdal.BuildVRT(str(temp_unw_out.with_suffix('.vrt')),
                      input_unw_files)
        gdal.BuildVRT(str(temp_conn_out.with_suffix('.vrt')),
                      input_conncomp_files)

    else:
        (combined_unwrap,
         combined_conn,
         combined_snwe) = stitch_unwrapped_frames(
            input_unw_files,
            input_conncomp_files,
            proj=f'EPSG:{ref_proj}',
            correction_method=correction_method,
            range_correction=range_correction,
            direction_N_S=True,
            verbose=verbose)

        # Write
        # write stitched unwrappedPhase
        write_GUNW_array(
            temp_unw_out, combined_unwrap, combined_snwe,
            format=output_format, epsg=int(ref_proj), verbose=verbose,
            update_mode=overwrite, add_vrt=True, nodata=0.0)

        # write stitched connectedComponents
        write_GUNW_array(
            temp_conn_out, combined_conn, combined_snwe,
            format=output_format, epsg=int(ref_proj), verbose=verbose,
            update_mode=overwrite, add_vrt=True, nodata=-1.0)

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

    for output, input in zip([output_conn, output_unw],
                             [temp_conn_out, temp_unw_out]):
        # Crop if selected
        ds = gdal.Warp(str(output),
                       str(input.with_suffix('.vrt')),
                       format=output_format,
                       cutlineDSName=clip_json,
                       xRes=arrres[0], yRes=arrres[1],
                       targetAlignedPixels=True,
                       # cropToCutline = True,
                       dstSRS=f'EPSG:{epsg}',
                       outputBounds=bounds,
                       outputType=gdal.GDT_Float32
                       )
        ds = None
        # Update VRT
        [print(f'Writing {output}, {output.with_suffix(".vrt")}')
         if verbose else None]
        gdal.Translate(str(output.with_suffix('.vrt')),
                       str(output), format="VRT")
        # Remove temp files
        [ii.unlink() for ii in [input, input.with_suffix('.vrt'),
                                input.with_suffix('.xml'),
                                input.with_suffix('.hdr'),
                                input.with_suffix('.aux.xml')] if ii.exists()]

        # Mask
        if mask_file:
            if isinstance(mask_file, str):
                mask = gdal.Open(mask_file)
            else:
                # for gdal instance, from prep_mask
                mask = mask_file

            mask_array = mask.ReadAsArray()
            array = get_GUNW_array(str(output.with_suffix('.vrt')),
                proj=f'EPSG:{epsg}')

            if output == output_conn:
                # Mask connected components
                array[array == -1.0] = np.nan
                update_array = mask_array * array
                update_array = np.nan_to_num(update_array, nan=-1.0)
                # update_array[np.isnan(update_array)] = -1.0
            else:
                update_array = mask_array * array

            update_file = gdal.Open(str(output), gdal.GA_Update)
            update_file = update_file.GetRasterBand(1).WriteArray(update_array)
            update_file = None
            
        # mask out zeros
        array = get_GUNW_array(str(output.with_suffix('.vrt')),
            proj=f'EPSG:{epsg}')
        if output == output_conn:
            # Mask connected components
            array[(array == -1) | (array == 0)] = np.nan
        else:
            concomp_array = \
                get_GUNW_array(str(output_conn.with_suffix('.vrt')),
                    proj=f'EPSG:{epsg}')
            concomp_array[(concomp_array == -1) \
                | (concomp_array == 0)] = np.nan
            concomp_array = np.isnan(concomp_array)
            array[array == 0] = np.nan
            array[concomp_array] = np.nan
        update_file = gdal.Open(str(output), gdal.GA_Update)
        update_file = update_file.GetRasterBand(1).WriteArray(array)
        update_file = None

    # Plot stitched
    # NOTE: saving output figure adds 4 seconds
    if save_fig:
        plot_GUNW_stitched(str(output_unw.with_suffix('.vrt')),
            str(output_conn.with_suffix('.vrt')))

    # Remove temp files
    [ii.unlink() for ii in [temp_unw_out,
                            temp_unw_out] if ii.exists()]


def plot_GUNW_stitched(stiched_unw_filename: str,
                       stiched_conn_filename: str, epsg: str) -> None:
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
    # bounds = np.linspace(0, 256, 257)

    output_dir = Path(stiched_unw_filename).absolute()
    output_fig = output_dir.parent / 'stitched.png'
    output_fig.unlink(missing_ok=True)

    # Load Data
    stitched_unw = get_GUNW_array(stiched_unw_filename, proj=f'EPSG:{epsg}')
    stitched_conn = get_GUNW_array(stiched_conn_filename, proj=f'EPSG:{epsg}')
    stitched_attr = get_GUNW_attr(stiched_unw_filename, proj=f'EPSG:{epsg}')

    # ConnComp discrete colormap
    # bounds = np.linspace(0, int(np.nanmax(stitched_conn)),
    #                      int(np.nanmax(stitched_conn)+1))
    bounds = np.linspace(0, 30, 31)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    # Mask
    stitched_unw[stitched_unw == 0.0] = np.nan
    stitched_conn[stitched_conn == -1.0] = np.nan

    # Common plot options
    plot_kwargs = {'extent': snwe_to_extent(stitched_attr['SNWE']),
                   'interpolation': 'nearest'}

    # Figure
    fig, axs = plt.subplots(1, 3, dpi=300, sharey=True)
    # Re-wrapped
    im1 = axs[0].imshow(np.mod(stitched_unw, 4*np.pi),
                        cmap='jet', **plot_kwargs)
    im2 = axs[1].imshow(stitched_unw * (0.0556 / (6*np.pi)), cmap='jet',
                        clim=[-0.2, 0.2], **plot_kwargs)  # Unwrapped
    # Connected Components
    im3 = axs[2].imshow(stitched_conn, cmap=cmap, norm=norm, **plot_kwargs)

    for im, ax, label, txt, in zip([im1, im2, im3],
                                   axs,
                                   ['rad', 'm', '#'],
                                   ['Re-wrapped phase w 20 rad',
                                    'Unwrapped phase [m]',
                                    'Connected Components']):
        fig.colorbar(im, ax=ax, location='bottom', shrink=0.8, label=label)
        ax.set_title(txt, fontsize=10)

    fig.tight_layout()
    fig.savefig(str(output_fig))
    fig.clear()
    plt.close(fig)
