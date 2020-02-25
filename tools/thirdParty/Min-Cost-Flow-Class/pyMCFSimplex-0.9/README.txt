*******************************************************************************
*			          pyMCFSimplex                                *
*******************************************************************************
* Version 0.9          *
************************
* Johannes Sommer, 2013*
************************

1. What?
--------
pyMCFSimplex is a Python-Wrapper for the C/C++ MCFSimplex Solver Class from 
the University of Pisa.
MCFSimplex is a Class that solves reasonable sized Minimum Cost Flow Problems
very fast.
See also [1] for a comparison.

2. How?
-------
pyMCFSimplex was being made through SWIG. Don't ask for the time I spent on
figuring out how SWIG works. With more knowledge in C/C++ I would have been
faster - but I'm not a C/C++ guy! I want it in Python!

3. Who?
-------
The authors of MCFSimplex are Alessandro Bertolini and Antonio Frangioni from
the Operations Research Group at the Dipartimento di Informatica of the 
University of Pisa.
http://www.di.unipi.it/optimize/Software/MCF.html#MCFSimplex

pyMCFSimplex is brought to you by Johannes from the G#.Blog 
http://www.sommer-forst.de/blog. 
Feel free to contact me: info(at)sommer-forst.de

4. Installation
---------------
Installation prerequisites are:
- Python 2.7 or Python 2.6 (only Windows)
- numpy (tested with 1.6.1)
- a build environment for source distribution install.

4.1 Windows
-----------
Select the appropriate MSI package for your installed Python version (2.6 or 2.7) an simply execute the installer.

4.2 Linux
---------
Untar the binary dist package pyMCFSimplex-0.9.linux-x86_64.tar.gz with

tar xfvz pyMCFSimplex-0.9.linux-x86_64.tar.gz

It will install into /usr/local/lib/python2.7/dist-packages/.

4.3 Source Distribution
-----------------------
Grab the pyMCFSimplex-0.9_src_dist.zip file, extract it and run

a) linux:
sudo python setup.py install

b) windows:
start a command line as Administrator and run

python setup.py install


5. Usage
--------
Hey! Wait a minute! Please note that pyMCFSimplex's version is 0.9. So this is not productive! I tested all the relevant methods and all of my tests went fine, but if you run into an error that occurs in the underlying C++ library Python might crash, because at the moment I don't know yet how to deal with C++ exceptions and SWIG. Feel free to contact me, if you know more on this!

Here is a first start. "sample.dmx" must be in the same location of your python script.
With these lines of code you can parse a minimum cost flow problem in DIMACS file format and solve it.

from pyMCFSimplex import *
print "pyMCFSimplex Version '%s' successfully imported." % version()
mcf = MCFSimplex()
print "MCFSimplex Class successfully instantiated."
FILENAME = 'sample.dmx'
print "Loading network from DIMACS file %s.." % FILENAME
f = open(FILENAME,'r')
inputStr = f.read()
f.close()
mcf.LoadDMX(inputStr)

print "Setting time.."
mcf.SetMCFTime()
print "Solving problem.."
mcf.SolveMCF()
if mcf.MCFGetStatus() == 0:
    print "Optimal solution: %s" %mcf.MCFGetFO()
    print "Time elapsed: %s sec " %(mcf.TimeMCF())
else:
    print "Problem unfeasible!"
    print "Time elapsed: %s sec " %(mcf.TimeMCF())

If you want to load a network not from a DIMACS file, you'll have to call LoadNet() while passing C-arrays to the method.
C arrays in Python? Yes - don't worry. There are helper methods in pyMCFSimplex, that'll do this for you.
Look at the following piece of code.

mcf = MCFSimplex()
print "MCFSimplex Class successfully instantiated."
print "Reading sample data.."

'''
Problem data of a MCFP in DIMACS notation

c Problem line (nodes, links)
p min 4 5
c
c Node descriptor lines (supply+ or demand-)
n 1 4
n 4 -4
c
c Arc descriptor lines (from, to, minflow, maxflow, cost)
a 1 2 0 4 2
a 1 3 0 2 2
a 2 3 0 2 1
a 2 4 0 3 3
a 3 4 0 5 1
'''

# MCFP problem transformed to integers and lists
nmx     = 4 # max number of nodes
mmx     = 5 # max number of arcs
pn      = 4 # current number of nodes
pm      = 5 # current number of arcs
pU      = [4,2,2,3,5] # column maxflow
pC      = [2,2,1,3,1] # column cost
pDfct   = [-4,0,0,4]  # node deficit (supply/demand)
pSn     = [1,1,2,2,3] # column from
pEn     = [2,3,3,4,4] # column to

# call LoadNet() with the return values of the helper methods
# e.g. CreateDoubleArrayFromList(pU) takes a python list and returns a pointer to a 
# corresponding C array, that is passed as an argument to the method LoadNet()
mcf.LoadNet(nmx, mmx, pn, pm, CreateDoubleArrayFromList(pU), CreateDoubleArrayFromList(pC),
            CreateDoubleArrayFromList(pDfct), CreateUIntArrayFromList(pSn),
            CreateUIntArrayFromList(pEn))

print "Setting time.."
mcf.SetMCFTime()
mcf.SolveMCF()
if mcf.MCFGetStatus() == 0:
    print "Optimal solution: %s" %mcf.MCFGetFO()
    print "Time elapsed: %s sec " %(mcf.TimeMCF())
else:
    print "Problem unfeasible!"
    print "Time elapsed: %s sec " %(mcf.TimeMCF())


Please check out the sample script gsharpblog_mcfsolve_test.py for more information.


6. Good to know
---------------
I changed the original source code of MCFClass.h a little bit for SWIG compatibility.
All changes are marked by the following comment line "//Johannes Sommer".
This included:
- LoadDMX() accepts in pyMCFSimplex a string value (original: c iostream). The original LoadDMX method is omitted.
- as SWIG cannot deal with nested classes, I pulled the classes Inf, MCFState and MCFException out of the main class MCFClass.

Perhaps the above mentioned changes to the original source is not necessary, if you know SWIG very well.
But I could not figure out how to get these things to work in the SWIG interface file.
Useful hints are very welcome.

[1] http://bolyai.cs.elte.hu/egres/tr/egres-13-04.ps

