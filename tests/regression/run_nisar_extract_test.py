#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Alex Fore
# Copyright (c) 2024, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import os
import shutil
import argparse
import subprocess
import contextlib
import tarfile
import logging

LOGGER = logging.getLogger('run_nisar_extract_test.py')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-l', '--log-level', default='info', help='Logger log level')
    parser.add_argument(
        '--old', default=False, action='store_true',
        help='Use old command line interface')
    args = parser.parse_args()

    FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    log_level = {
        'debug': logging.DEBUG, 'info': logging.INFO,
        'warning': logging.WARNING, 'error': logging.ERROR}[args.log_level]
    logging.basicConfig(level=log_level, format=FORMAT)

    def run_subproc(exec_string, message, raise_exception=False):
        """Helper function to run subprocess and handle return values"""
        return_code = subprocess.call(exec_string, shell=True)
        if return_code != 0:
            LOGGER.error('%s failed!' % message)
            if raise_exception:
                raise Exception('%s failed!' % message)

    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree('golden_test_inputs/nisar_extract')

    # Uncompress .tar.gz file with input data in it
    with tarfile.open(
            os.path.join('golden_test_inputs/nisar_extract.tar.gz')) as tar:
        tar.extractall('golden_test_inputs')

    with contextlib.suppress(FileNotFoundError):
        shutil.rmtree('test_outputs/nisar_extract')
    os.makedirs('test_outputs/nisar_extract')

    # extract azimuth angle
    exec_string = (
        'ariaExtract.py -f "golden_test_inputs/nisar_extract/products/*h5" '
        '-l "unwrappedPhase,troposphereTotal,ionosphere" '
        '-w test_outputs/nisar_extract '
        '-d golden_test_inputs/nisar_extract/DEM/glo_90.dem')
    if not args.old:
        exec_string += ' --log-level %s' % args.log_level
    run_subproc(exec_string, 'ariaExtract', raise_exception=True)


if __name__ == "__main__":
    main()
