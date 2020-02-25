# Min-Cost-Flow-Class

ReadMe for the MCFClass project, a set of C++ solvers for (Linear or Convex
Quadratic Separable) Min Cost Flow Problem solvers under the same interface
deriving from the base class MCFClass.

The aim of MCFClass is to provide an abstraction layer between practitioners
who need to solve MCF problems within complex applications and developers of
MCF software. The idea is to provide an interface which caters for all the
needs that a practitioner can have, thereby allowing him/her to use whichever
algorithm - among those that have an implementation conforming to this
interface - without bothering with the details of the implementation, and to
easily switch between different algorithms.

MCFClass defines and "exports" the types of "flows" (`MCFClass::FNumber`),
"costs" (`MCFClass::CNumber`) and so on, together with a set of comparison
operators (ETZ, GTZ, ...) which automatically detect whether or not the
underlying types are integers or floats, inserting appropriate "epsilons" in
the latter case and avoiding them (for speed) in the former; these things are
sorted out at compile time without user intervention.

This release comprises:

-  [`docs/`](docs): HTML doxygen documentation, also available at

    https://frangio68.github.io/Min-Cost-Flow-Class/

-  [`doxygen/`](doxygen): files to produce the documentation

-  [`License.md`](License.md): the text of the "GNU Lesser General Public License",
   Version 3.0, under which most of this code is distributed
   (but not all of it, see RelaxIV below)

-  [`MCFClass/`](MCFClass): definition of the base class

-  [`MCFClone/`](MCFClone): implements a "fake" MCF solver that takes two "real" ones
   and does everything on both; useful for testing the solvers (either for
   correctness or for efficiency) when used within "complex" approaches

-  [`MCFCplex/`](MCFCplex): implements a MCF solver conforming to the MCFClass interface
   based on calls to the commercial (but now free for academic purposes) Cplex solver from IBM

-  [`MCFSimplex/`](MCFSimplex): implements a MCF solver conforming to the MCFClass interface
   based on the primal and dual revised network simplex algorithm

-  [`OPTUtils/`](OPTUtils): contains the `OPTUtils.h` file with a few minor utility functions

-  [`ReadMe.md`](ReadMe.md): this file

-  [`RelaxIV/`](RelaxIV): implements a MCF solver conforming to the `MCFClass` interface
   based on the RELAXIV code by D. Bertsekas and P. Tseng, as described in
       Bertsekas, Dimitri P., and Paul Tseng.
       "RELAX-IV: A faster version of the RELAX code for solving minimum
       cost flow problems." (1994), Report LIDS-P-2276, MIT.

   be aware that RelaxIV is distributed under a less permissive academic
   license than the rest of the code (see [`RelaxIV/academicl.txt`](RelaxIV/academicl.txt) for details),
   which only applies to researchers of noncommercial and academic
   institutions, e.g., universities; if the license does not apply to you,
   you must not download the source or delete it immediately

-  [`SPTree/`](SPTree): implements a MCF solver partly conforming to the MCFClass
   interface, in the sense that is is only able to solve MCF instances that
   are in fact Shortest Path Tree ones (that is, only one source node and
   no arc capacities), but then does so using SPT algorithms (both
   label-setting and label-correcting variants can be used) that are much
   faster than complete MCF ones

-  [`pyMCFSimplex-0.9/`](pyMCFSimplex-0.9): a Python-Wrapper for the MCFSimplex solver by Johannes
   from the G#.Blog, check `README.txt` for details

-  [`test/`](test): contains two example Main files to use the library. One solves
   a given MCF instance with any one MCF solver, which can be chosen by
   just changing two lines of code. The other compares the results of two
   solvers in order to verify that they agree. See the comments in both
   files for more details

There are two more complete solvers available under the `MCFClass` interface,
namely CS2 and MCFZIB. These are, however, distributed under a more
restrictive academic license, which has to be explicitly accepted before
getting hold of the code. Request forms are available at

  http://www.di.unipi.it/optimize/Software/MCF.html

## Build and install

You can either use [CMake](https://cmake.org) or plain makefiles to build the
library, your choice.

### Using CMake

- Clone the project from the repository and navigate inside its main directory.

- If you installed the requirements you should be fine. Configure the project with:
```sh
mkdir build
cd build
cmake ..
```

- You can now build the library:
```sh
make
```

- Optionally, you can install the library with:
```sh
sudo make install
```

- After the library is configured and built, you can use it in your CMake project with:
```cmake
find_package(MCFClass)
target_link_libraries(<my_target> MCFClass::MCFClass)
```

### Using makefiles

- To create the library, go into `lib/` and type
```sh
make -f makefile-lib
```

- To test the library, go into `test` and type `make`.

- If you want to use the MCFCplex class, which comes commented out by default,
  uncomment the two lines in `lib/makefile`:
```makefile
MCFCxDIR = $(libMCFClDIR)MCFCplex/
include $(MCFCxDIR)makefile
```

  Then edit `extlib/makefile-libCPLEX` to insert the right Cplex
  path libraries. You can similarly enable (or disable) any solver, both the
  LGPL ones and those under the academic license, if you have obtained them,
  by commenting out (or commenting) the corresponding two lines in `lib/makefile`.


## Other stuff

More information about (some of) the implemented algorithms can be found at

  http://pages.di.unipi.it/frangio/abstracts.html#JOC06

A further solver, MCFIntPnt (based on Interior-Point algorithms) has been
developed, but it has not yet reached a sufficient maturuty to be distributed;
its principles are discussed at

  http://pages.di.unipi.it/frangio/abstracts.html#SIOPT04
  http://pages.di.unipi.it/frangio/abstracts.html#COAP06
  http://pages.di.unipi.it/frangio/abstracts.html#OMS06
