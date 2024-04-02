#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Alex Fore
# Copyright (c) 2024, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import argparse
import numpy as np
import logging
import tarfile

import osgeo.gdal

ENVI_FILES = {
    'azimuthAngle': 'azimuthAngle/20230829_20230724',
    'troposphereTotal_HRRR': 'troposphereTotal/HRRR/20230829_20230724',
    'troposphereTotal_HRRR_dates_20230724': \
        'troposphereTotal/HRRR/dates/20230724',
    'troposphereTotal_HRRR_dates_20230829': \
        'troposphereTotal/HRRR/dates/20230829',
    'coherence': 'coherence/20230829_20230724',
    'unwrappedPhase': 'unwrappedPhase/20230829_20230724',
    'connectedComponents': 'connectedComponents/20230829_20230724',
    'ionosphere': 'ionosphere/20230829_20230724'
}

ENVI_FILES_TSSETUP = {
    'azimuthAngle': 'azimuthAngle/20230829_20230724',
    'bPerpendicular_20230829_20230724': 'bPerpendicular/20230829_20230724',
    'bPerpendicular_20230829_20230817': 'bPerpendicular/20230829_20230817',
    'connectedComponents_20230829_20230724': \
        'connectedComponents/20230829_20230724',
    'connectedComponents_20230829_20230817': \
        'connectedComponents/20230829_20230817',
    'incidenceAngle': 'incidenceAngle/20230829_20230724',
    'ionosphere_20230829_20230724': 'ionosphere/20230829_20230724',
    'ionosphere_20230829_20230817': 'ionosphere/20230829_20230817',
    'lookAngle': 'lookAngle/20230829_20230724',
    'solidEarthTide_20230829_20230724': 'solidEarthTide/20230829_20230724',
    'solidEarthTide_20230829_20230817': 'solidEarthTide/20230829_20230817',
    'troposphereTotal_HRRR_20230829_20230817': \
        'troposphereTotal/HRRR/20230829_20230817',
    'troposphereTotal_HRRR_20230829_20230724': \
        'troposphereTotal/HRRR/20230829_20230724',
    'troposphereTotal_HRRR_dates_20230724': \
        'troposphereTotal/HRRR/dates/20230724',
    'troposphereTotal_HRRR_dates_20230817': \
        'troposphereTotal/HRRR/dates/20230817',
    'troposphereTotal_HRRR_dates_20230829': \
        'troposphereTotal/HRRR/dates/20230829',
    'unwrappedPhase_20230829_20230724': 'unwrappedPhase/20230829_20230724',
    'unwrappedPhase_20230829_20230817': 'unwrappedPhase/20230829_20230817',
}

LOGGER = logging.getLogger('validate_test.py')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'flavor', choices=['extract', 'tssetup'],
        help="Flavor: either 'extract or tssetup'")
    parser.add_argument(
        '-l', '--log-level', default='info', help='Logger log level')
    args = parser.parse_args()

    format = '%(levelname)s %(name)s %(message)s'

    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=format)

    test_dir = os.path.join('test_outputs', args.flavor)
    ref_dir = os.path.join('golden_test_outputs', args.flavor)

    if args.flavor == 'extract':
        iter_dict = ENVI_FILES

    else:
        iter_dict = ENVI_FILES_TSSETUP

    with tarfile.open(os.path.join(
            'golden_test_outputs', args.flavor+'.tar.gz')) as tar:
        tar.extractall(os.path.join('golden_test_outputs'))

    any_fail = False
    for filetype, filepath in iter_dict.items():
        test_file = os.path.join(test_dir, filepath)
        ref_file = os.path.join(ref_dir, filepath)

        if not os.path.isfile(test_file):
            LOGGER.warning(
                'For key: {}, output file does not exist!'.format(filetype))
            any_file = True
            continue

        test = osgeo.gdal.Open(test_file)
        ref = osgeo.gdal.Open(ref_file)

        test_data = test.ReadAsArray()
        ref_data = ref.ReadAsArray()

        # Compare mean/std of data arrays
        fill_value = ref.GetRasterBand(1).GetNoDataValue()

        mask_ref = ref_data != fill_value
        mask_test = test_data != fill_value

        # Verify all values are finite
        num_nonfinite_test = (~np.isfinite(test_data[mask_test])).sum()
        num_nonfinite_ref = (~np.isfinite(ref_data[mask_ref])).sum()

        if num_nonfinite_test > 0:
            LOGGER.warning((
                'For test data key: {}, {} non-finite value(s)!').format(
                    filetype, num_nonfinite_test))
            any_fail = True

        if num_nonfinite_ref > 0:
            LOGGER.warning((
                'For ref data key: {}, {} non-finite value(s)!').format(
                    filetype, num_nonfinite_ref))
            any_fail = True

        # Check that they have the same shape
        if not mask_ref.shape == mask_test.shape:
            LOGGER.warning((
                'For key: {}, ref and test data have different shape, '
                'ref: {}, test: {}!').format(
                    filetype, mask_ref.shape, mask_test.shape))
            any_fail = True

        else:
            # Check that the same pixels have fill value in both arrays
            if not (mask_ref == mask_test).all():
                LOGGER.warning((
                    'For key: {}, ref and test data have fill values in '
                    'different places!').format(filetype))
                any_fail = True

            mask_joint = np.logical_and.reduce((
                mask_ref, np.isfinite(test_data), np.isfinite(ref_data)))

            # Compute mean/std of the data values
            delta = (test_data[mask_joint] - ref_data[mask_joint]) / np.median(
                test_data[mask_joint])

            if not (np.abs(delta) < 10**-8).all():
                LOGGER.warning((
                    'For key: {}, ref and test data differ by more than '
                    'threshold!').format(filetype))
                any_fail = True

            mean_diff = (test_data[mask_joint]-ref_data[mask_joint]).mean()
            std_diff = (test_data[mask_joint]-ref_data[mask_joint]).std()
            LOGGER.info('For key %s, mean/std diff: %e %e' % (
                filetype, mean_diff, std_diff))

        if not any_fail:
            LOGGER.info("For key: {}, all tests passed.".format(filetype))

if __name__ == "__main__":
    main()
