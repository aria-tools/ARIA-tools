#! /usr/bin/env python3
"""ARIA-tool to run time series preparation.
Copyright 2019, by the California Institute of Technology. ALL RIGHTS
RESERVED. United States Government Sponsorship acknowledged.

Author(s): Simran Sangha, David Bekaert, & Emre Havazli
"""
# Import functions
from ARIAtools.tsSetup import cmdLineParse,main
from pkg_resources import get_distribution

if __name__ == '__main__':
    print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    inps = cmdLineParse()
    main(inps)
