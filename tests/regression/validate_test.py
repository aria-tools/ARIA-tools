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
import pytest
import subprocess

import osgeo.gdal

ENVI_FILES = {}
ENVI_FILES['extract'] = {
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
ENVI_FILES['tssetup'] = {
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

LOGGER = logging.getLogger(__name__)

class ENVIDataTester(object):
    def __init__(self, test_dir, ref_dir, flavor):
        self.test_data = {}
        self.ref_data = {}

        for filetype, filepath in ENVI_FILES[flavor].items():

            test_file = os.path.join(test_dir, filepath)
            test = osgeo.gdal.Open(test_file)
            self.test_data[filetype] = {
                'data': test.ReadAsArray(),
                'fill_value': test.GetRasterBand(1).GetNoDataValue(),
                'filename': test_file}

            ref_file = os.path.join(ref_dir, filepath)
            ref = osgeo.gdal.Open(ref_file)
            self.ref_data[filetype] = {
                'data': ref.ReadAsArray(),
                'fill_value': ref.GetRasterBand(1).GetNoDataValue(),
                'filename': ref_file}

    def test_nonfinite(self):
        any_fail = False
        for filetype in self.ref_data:
            assert filetype in self.test_data

            ref_data = self.ref_data[filetype]['data']
            test_data = self.test_data[filetype]['data']
            fill_value = self.ref_data[filetype]['fill_value']

            test_file = self.test_data[filetype]['filename']
            ref_file = self.ref_data[filetype]['filename']

            mask_ref = ref_data != fill_value
            mask_test = test_data != fill_value

            num_nonfinite_test = (~np.isfinite(test_data[mask_test])).sum()
            num_nonfinite_ref = (~np.isfinite(ref_data[mask_ref])).sum()

            if num_nonfinite_test > 0:
                LOGGER.warning((
                    'For test data file: {}, {} non-finite value(s)!').format(
                        test_file, num_nonfinite_test))
                any_fail = True

            if num_nonfinite_ref > 0:
                LOGGER.warning((
                    'For ref data file: {}, {} non-finite value(s)!').format(
                        ref_file, num_nonfinite_ref))

        if any_fail:
            raise AssertionError('Non-finite values found in test data!')

    def test_shape(self):
        any_fail = False
        for filetype in self.ref_data:
            assert filetype in self.test_data

            ref_data = self.ref_data[filetype]['data']
            test_data = self.test_data[filetype]['data']
            test_file = self.test_data[filetype]['filename']

            if not ref_data.shape == test_data.shape:
                LOGGER.warning((
                    'For file: {}, ref and test data have different shape, '
                    'ref: {}, test: {}!').format(
                        test_file, mask_ref.shape, mask_test.shape))
                any_fail = True

        if any_fail:
            raise AssertionError(
                'Different array shapes found in test and ref data!')

    def test_fill(self):
        any_fail = False
        for filetype in self.ref_data:
            assert filetype in self.test_data

            ref_data = self.ref_data[filetype]['data']
            test_data = self.test_data[filetype]['data']
            fill_value = self.ref_data[filetype]['fill_value']

            test_file = self.test_data[filetype]['filename']
            ref_file = self.ref_data[filetype]['filename']

            mask_ref = ref_data != fill_value
            mask_test = test_data != fill_value

            if not (mask_ref == mask_test).all():
                LOGGER.warning((
                    'For file: {}, ref and test data have fill values in '
                    'different places!').format(test_file))
                any_fail = True

        if any_fail:
            raise AssertionError((
                'Fill values found in different array locations in test '
                'and ref data!'))

    def test_values(self):
        any_fail = False
        for filetype in self.ref_data:
            assert filetype in self.test_data

            ref_data = self.ref_data[filetype]['data']
            test_data = self.test_data[filetype]['data']
            fill_value = self.ref_data[filetype]['fill_value']

            test_file = self.test_data[filetype]['filename']
            ref_file = self.ref_data[filetype]['filename']

            mask_ref = ref_data != fill_value
            mask_test = test_data != fill_value

            if not np.allclose(test_data, ref_data, equal_nan=True):
                LOGGER.warning((
                    'For file: {}, ref and test data differ by more than '
                    'threshold!').format(test_file))
                any_fail = True

        if any_fail:
            raise AssertionError(
                'Test and ref data differ by more than threshold')

class AriaToolsScriptTester():
    @pytest.fixture(scope='class')
    def tester(self):
        flavor = 'tssetup'
        with tarfile.open(os.path.join(
                'golden_test_outputs', self.FLAVOR+'.tar.gz')) as tar:
            tar.extractall(os.path.join('golden_test_outputs'))

        test_dir = os.path.join('test_outputs', self.FLAVOR)
        ref_dir = os.path.join('golden_test_outputs', self.FLAVOR)
        return ENVIDataTester(test_dir, ref_dir, self.FLAVOR)

    def test_nonfinite(self, tester):
        return tester.test_nonfinite()

    def test_shape(self, tester):
        return tester.test_shape()

    def test_fill(self, tester):
        return tester.test_fill()

    def test_values(self, tester):
        return tester.test_values()

@pytest.mark.usefixtures('run_tssetup_test')
class TestAriaTSsetup(AriaToolsScriptTester):
    FLAVOR='tssetup'

@pytest.mark.usefixtures('run_extract_test')
class TestAriaExtract(AriaToolsScriptTester):
    FLAVOR='extract'

@pytest.fixture(scope='session')
def sync_golden_data():
    return_code = subprocess.call(
        'aws s3 sync s3://aria-tools/tests/regression/ . --no-sign-request',
        shell=True)
    if return_code != 0:
        LOGGER.error("Error syncing golden test data!")
        raise AssertionError("Error syncing golden test data!")

@pytest.fixture(scope='session')
def run_extract_test(sync_golden_data):
    return_code = subprocess.call('./run_extract_test.py', shell=True)
    if return_code != 0:
        LOGGER.error("Error running ariaExtract test case!")
        raise AssertionError("Error running ariaExtract test case!")

@pytest.fixture(scope='session')
def run_tssetup_test(sync_golden_data):
    return_code = subprocess.call('./run_tssetup_test.py', shell=True)
    if return_code != 0:
        LOGGER.error("Error running ariaTSsetup test case!")
        raise AssertionError("Error running ariaTSsetup test case!")
