##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of libMCFCl						     #
#   all modules are given in input                                           #
#									     #
#   Input:  $(CC)        = compiler command				     #
#           $(SW)        = compiler options                                  #
#           $(OPTH)      = the OPT****.h files to include		     #
#           $(OPTINC)    = the -I$(OPT****.h directory)			     #
#           $(libMCFClDIR) = the basic directory for the library	     #
#									     #
#   Output: $(libMCFClOBJ) = the library (that must be built)		     #
#           $(libMCFClLIB) = the external libreries + -L<libdirs>            #
#           $(libMCFClH)   = all the .h files of the MCF library	     #
#           $(libMCFClINC) = all the -I$(include dirs) of the MCF library    #
#									     #
#                                VERSION 2.10				     #
#                	        24 - 04 - 2017				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# MCFClass
MCFClDIR = $(libMCFClDIR)MCFClass/
include $(MCFClDIR)makefile

# CS2
#MCFC2DIR = $(libMCFClDIR)CS2/
#include $(MCFC2DIR)makefile

# MCFClone
MCFCeDIR = $(libMCFClDIR)MCFClone/
include $(MCFCeDIR)makefile

# MCFCplex
#MCFCxDIR = $(libMCFClDIR)MCFCplex/
#include $(MCFCxDIR)makefile

# MCFSimplex
MCFSxDIR = $(libMCFClDIR)MCFSimplex/
include $(MCFSxDIR)makefile

# MCFZIB
#MCFZBDIR = $(libMCFClDIR)MCFZIB/
#include $(MCFZBDIR)makefile

# RelaxIV
MCFR4DIR = $(libMCFClDIR)RelaxIV/
include $(MCFR4DIR)makefile

# SPTree
MCFSTDIR = $(libMCFClDIR)SPTree/
include $(MCFSTDIR)makefile

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

# object files (this is local)
libMCFClOBJ =	$(MCFC2OBJ)\
		$(MCFCxOBJ)\
		$(MCFSxOBJ)\
		$(MCFZBOBJ)\
		$(MCFSTOBJ)\
		$(MCFR4OBJ)
# libreries
libMCFClLIB =	$(MCFC2LIB)\
		$(MCFCxLIB)\
		$(MCFSxLIB)\
		$(MCFZBLIB)\
		$(MCFSTLIB)\
		$(MCFR4LIB)
# includes
libMCFClH   =	$(MCFClH)\
		$(MCFC2H)\
		$(MCFCeH)\
		$(MCFCxH)\
		$(MCFSxH)\
		$(MCFZBH)\
		$(MCFR4H)\
		$(MCFSTH)
# include dirs
libMCFClINC =	$(MCFClINC)\
		$(MCFC2INC)\
		$(MCFCeINC)\
		$(MCFCxINC)\
		$(MCFSxINC)\
		$(MCFZBINC)\
		$(MCFR4INC)\
		$(MCFSTINC)

############################ End of makefile #################################
