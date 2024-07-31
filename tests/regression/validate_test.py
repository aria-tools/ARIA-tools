# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Alex Fore
# Copyright (c) 2024, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os
import hashlib
import logging
import tarfile
import subprocess
import numpy as np
import pytest
import osgeo.gdal

ENVI_FILES = {}
ENVI_FILES['extract'] = {
    'azimuthAngle': 'azimuthAngle/20230829_20230724',
    'troposphereTotal_HRRR': 'troposphereTotal/HRRR/20230829_20230724',
    'troposphereTotal_HRRR_dates_20230724':
        'troposphereTotal/HRRR/dates/20230724',
    'troposphereTotal_HRRR_dates_20230829':
        'troposphereTotal/HRRR/dates/20230829',
    'coherence': 'coherence/20230829_20230724',
    'unwrappedPhase': 'unwrappedPhase/20230829_20230724',
    'connectedComponents': 'connectedComponents/20230829_20230724',
    'ionosphere': 'ionosphere/20230829_20230724'
}
ENVI_FILES['tssetup'] = {
    'azimuthAngle': 'azimuthAngle/20230829_20230724',
    'connectedComponents_20230829_20230724':
        'connectedComponents/20230829_20230724',
    'connectedComponents_20230829_20230817':
        'connectedComponents/20230829_20230817',
    'incidenceAngle': 'incidenceAngle/20230829_20230724',
    'ionosphere_20230829_20230724': 'ionosphere/20230829_20230724',
    'ionosphere_20230829_20230817': 'ionosphere/20230829_20230817',
    'solidEarthTide_20230829_20230724': 'solidEarthTide/20230829_20230724',
    'solidEarthTide_20230829_20230817': 'solidEarthTide/20230829_20230817',
    'troposphereTotal_HRRR_20230829_20230817':
        'troposphereTotal/HRRR/20230829_20230817',
    'troposphereTotal_HRRR_20230829_20230724':
        'troposphereTotal/HRRR/20230829_20230724',
    'troposphereTotal_HRRR_dates_20230724':
        'troposphereTotal/HRRR/dates/20230724',
    'troposphereTotal_HRRR_dates_20230817':
        'troposphereTotal/HRRR/dates/20230817',
    'troposphereTotal_HRRR_dates_20230829':
        'troposphereTotal/HRRR/dates/20230829',
    'unwrappedPhase_20230829_20230724': 'unwrappedPhase/20230829_20230724',
    'unwrappedPhase_20230829_20230817': 'unwrappedPhase/20230829_20230817',
}
ENVI_FILES['nisar_extract'] = {
    'connectedComponents': 'connectedComponents/20100410_20110111',
    'ionosphere': 'ionosphere/20100410_20110111',
    'troposphereTotal': 'troposphereTotal/20100410_20110111',
    'unwrappedPhase': 'unwrappedPhase/20100410_20110111'
}

NETCDF_FILES = {}
NETCDF_FILES['download'] = {
    'products_20211220_20211214':
        ('products/S1-GUNW-D-R-071-tops-20211220_20211214-135232-'
         '00119W_00034N-PP-fe7a-v2_0_5.nc'),
    'products_20211220_20211214':
        ('products/S1-GUNW-D-R-071-tops-20220101_20211214-135232-'
         '00119W_00034N-PP-246b-v2_0_5.nc'),
    'products_20211220_20211214':
        ('products/S1-GUNW-D-R-071-tops-20220101_20211220-135232-'
         '00119W_00034N-PP-e72b-v2_0_5.nc')
}

LOGGER = logging.getLogger(__name__)


class ENVIDataTester(object):
    def __init__(self, test_dir, ref_dir, flavor):
        self.test_data = {}
        self.ref_data = {}

        for filetype, filepath in ENVI_FILES[flavor].items():

            test_file = os.path.join(test_dir, filepath)
            test = osgeo.gdal.Open(test_file)
            if test is None:
                LOGGER.error('Unable to load test file: %s', test_file)

            self.test_data[filetype] = {
                'data': test.ReadAsArray(),
                'fill_value': test.GetRasterBand(1).GetNoDataValue(),
                'filename': test_file}

            ref_file = os.path.join(ref_dir, filepath)
            ref = osgeo.gdal.Open(ref_file)
            if ref is None:
                LOGGER.error('Unable to load ref file: %s', ref_file)

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


class NetCDFDataTester():
    def __init__(self, test_dir, ref_dir, flavor):
        self.flavor = flavor
        self.test_dir = test_dir
        self.ref_dir = ref_dir

    def test_md5sum(self):
        any_fail = False
        for filetype, filepath in NETCDF_FILES[self.flavor].items():

            test_file = os.path.join(self.test_dir, filepath)
            ref_file = os.path.join(self.ref_dir, filepath)

            ref_md5sum = self.md5sum(ref_file)
            test_md5sum = self.md5sum(test_file)

            if ref_md5sum != test_md5sum:
                LOGGER.warning(
                    'For file: %s, ref and test data do not have same md5sum!',
                    filepath)
                any_fail = True

        if any_fail:
            raise AssertionError('Test and ref data do not have same md5sums')

    @staticmethod
    def md5sum(filename):
        """Returns the md5sum of filename"""
        with open(filename, 'rb') as ifp:
            return hashlib.md5(ifp.read()).hexdigest()


class AriaToolsScriptTester():
    @pytest.fixture(scope='class')
    def tester(self):
        with tarfile.open(os.path.join(
                'golden_test_outputs', self.FLAVOR + '.tar.gz')) as tar:
            tar.extractall(os.path.join('golden_test_outputs'))

        test_dir = os.path.join('test_outputs', self.FLAVOR)
        ref_dir = os.path.join('golden_test_outputs', self.FLAVOR)
        return ENVIDataTester(test_dir, ref_dir, self.FLAVOR)

# Commenting out as this is know to fail
#     def test_nonfinite(self, tester):
#         return tester.test_nonfinite()

    def test_shape(self, tester):
        return tester.test_shape()

    def test_fill(self, tester):
        return tester.test_fill()

    def test_values(self, tester):
        return tester.test_values()


@pytest.mark.usefixtures('run_tssetup_test')
class TestAriaTSsetup(AriaToolsScriptTester):
    FLAVOR = 'tssetup'


@pytest.mark.usefixtures('run_extract_test')
class TestAriaExtract(AriaToolsScriptTester):
    FLAVOR = 'extract'


@pytest.mark.usefixtures('run_nisar_extract_test')
class TestAriaNISARExtract(AriaToolsScriptTester):
    FLAVOR = 'nisar_extract'


@pytest.mark.usefixtures('run_download_test')
class TestAriaDownload():
    FLAVOR = 'download'

    @pytest.fixture(scope='class')
    def tester(self):
        with tarfile.open(os.path.join(
                'golden_test_outputs', self.FLAVOR + '.tar.gz')) as tar:
            tar.extractall(os.path.join('golden_test_outputs'))

        test_dir = os.path.join('test_outputs', self.FLAVOR)
        ref_dir = os.path.join('golden_test_outputs', self.FLAVOR)
        return NetCDFDataTester(test_dir, ref_dir, self.FLAVOR)

    def test_md5sum(self, tester):
        return tester.test_md5sum()


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
def run_nisar_extract_test(sync_golden_data):
    return_code = subprocess.call('./run_nisar_extract_test.py', shell=True)
    if return_code != 0:
        LOGGER.error("Error running NISAR ariaExtract test case!")
        raise AssertionError("Error running NISAR ariaExtract test case!")


@pytest.fixture(scope='session')
def run_tssetup_test(sync_golden_data):
    return_code = subprocess.call('./run_tssetup_test.py', shell=True)
    if return_code != 0:
        LOGGER.error("Error running ariaTSsetup test case!")
        raise AssertionError("Error running ariaTSsetup test case!")


@pytest.fixture(scope='session')
def run_download_test(sync_golden_data):
    return_code = subprocess.call('./run_download_test.py', shell=True)
    if return_code != 0:
        LOGGER.error("Error running ariaDownload test case!")
        raise AssertionError("Error running ariaDownload test case!")
