##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of MCFSolve						     #
#									     #
#   'make clean' cleans up						     #
#   'make' or 'make MCFSolve' builds the module				     #
#									     #
#                                VERSION 3.10				     #
#                	        16 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# basic directory
DIR = ./

# basic libreries
LIB = -lm

# debug switches
SW = -Wall -g -std=c++14

# production switches
#SW = -O3 -DNDEBUG -std=c++14

# compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: MCFSolve TestMCF

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libreries + -L<libdirs>
# *H   is the list of all include files
# *INC is the -I<include directories>

# the whole libMCF, with all required external libreries
libMCFClDIR = ../
include $(libMCFClDIR)lib/makefile-inc

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

MCFSolve: Main.o $(libMCFClOBJ)
	$(CC) -o MCFSolve Main.o $(libMCFClOBJ) $(libMCFClLIB) $(LIB) $(SW)

TestMCF: TestMain.o $(libMCFOBJ)
	$(CC) -o TestMCF TestMain.o $(libMCFClOBJ) $(libMCFClLIB) $(LIB) $(SW)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(DIR)Main.o: $(DIR)Main.C $(libMCFClH)
	$(CC) -c $*.C -o $@ $(libMCFClINC) $(SW)

$(DIR)TestMain.o: $(DIR)TestMain.C $(libMCFClH)
	$(CC) -c $*.C -o $@ $(libMCFClINC) $(SW)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -Rf MCFSolve TestMCF *.o *~ *.dSYM

############################ End of makefile #################################
