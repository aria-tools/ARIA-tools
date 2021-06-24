#! /usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Simran Sangha & David Bekaert
# Copyright 2019, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Import ARIA-tools functions
from ARIAtools.extractProduct import cmdLineParse,main
from pkg_resources import get_distribution

if __name__ == '__main__':
    try:
        print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except:
        pass
    inps = cmdLineParse()
    main(inps)
