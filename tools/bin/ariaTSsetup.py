#! /usr/bin/env python3
"""
ARIA-tool to run time series preparation.

Copyright 2019, by the California Institute of Technology. ALL RIGHTS
RESERVED. United States Government Sponsorship acknowledged.
Author(s): Simran Sangha, David Bekaert, & Emre Havazli
"""
# Import functions
from ARIAtools.ts_setup import cmd_line_parse, main

if __name__ == '__main__':
    inps = cmd_line_parse()
    main(inps)
