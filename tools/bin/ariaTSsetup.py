#! /usr/bin/env python3
"""ARIA-tool to run time series preparation.
Copyright (c) 2023, by the California Institute of Technology. ALL RIGHTS
RESERVED. United States Government Sponsorship acknowledged.

Author(s): Simran Sangha, David Bekaert, & Emre Havazli
"""
# Import functions
from ARIAtools.tsSetup import cmd_line_parse, main
from pkg_resources import get_distribution

if __name__ == '__main__':
    try:
        print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except:
        pass
    inps = cmd_line_parse()
    main(inps)
