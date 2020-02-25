Introduction
------------

This file is a short user manual for our distribution of CS2, a C++
implementation of the latest variants of Relaxation algorithms for linear
Min Cost Flow problems.

Standard Disclaimer
-------------------

This code is provided "as is", without any explicit or implicit warranty that
it will properly behave or it will suit you needs. Although codes reaching the
distribution phase have usually been extensively tested, we cannot guarantee
that they are absolutely bug-free (who can?). Any use of the codes is at you
own risk: in no case we could be considered liable for any damage or loss you
may suffer, either directly or indirectly, for having used this code. More
details about the non-warranty attached to this code are available in the
license description file.

The code also comes with a "good will only" support: feel free to contact us
for any comments/critics/bug report/request help you may have, we will be
happy to try to answer and help you. But we cannot spend much time solving
your problems, just the time to read a couple of e-mails and send you fast
suggestions if the problem is easily solvable. Apart from that, we can’t offer
you any support.

License
-------

This code is provided free of charge for academic purposes under the "academic
license": see the file academicl.txt. Commercial use of this code is only
allowed at a fee: please contact the original author dimitrib@mit.edu for
further information.

How to use it
-------------

This solver is a part of the MCFClass project. To use it:

- Download the current version of the MCFClass project from

    https://github.com/frangio68/Min-Cost-Flow-Class

- Copy the RelaxIV directory in the root directory of the MCFClass project
  (at the same level as lib/, MCFClass/, ...)

- Comment out the two lines

    #MCFR4DIR = $(libMCFClDIR)RelaxIV/
    #include $(MCFR4DIR)makefile

  in lib/makefile

- Build the library by typing

    make -f makefile-lib

  into the lib/ folder

- Alternatively, build the library and test executables by typing "make"
  in the test/ folder.
