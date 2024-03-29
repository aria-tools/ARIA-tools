#! /usr/bin/env python3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Compute and analyze the phase triplet misclosure for a series of
# interferograms.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import functions
from ARIAtools.computeMisclosure import cmdLineParse,main
from pkg_resources import get_distribution

if __name__ == '__main__':
    try:
        print ('ARIA-tools Version:', get_distribution('ARIAtools').version)
    except:
        pass
    inps = cmdLineParse()
    main(inps)
