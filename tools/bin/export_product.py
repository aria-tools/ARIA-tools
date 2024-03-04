#!/usr/bin/env python
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Alex Fore
# Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import argparse
import json
import os

import ARIAtools.extractProduct

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'export_products_json', help='Export products json input file')
    args = parser.parse_args()

    with open(args.export_products_json) as ifp:
        export_product_args = json.load(ifp)

    ii, ilayer, outname, prod_arr = \
        ARIAtools.extractProduct.export_product_worker(*export_product_args)

    outputs = {
        'ii': ii, 'ilayer': ilayer, 'outname': outname, 'prod_arr': prod_arr}

    base_dir = os.path.abspath(os.path.dirname(args.export_products_json))
    outfile = os.path.join(base_dir, 'outputs_%2.2d_%2.2d.json' % (
        ii, ilayer))

    with open(outfile, 'w') as ofp:
        json.dump(outputs, ofp)

if __name__ == "__main__":
    main()
