##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of MCFSimplex 						     #
#									     #
#   Input:  $(CC)       = compiler command				     #
#           $(SW)       = compiler options                                   #
#           $(MCFClH)   = the include files for MCFClass		     #
#           $(MCFClINC) = the -I$(MCFClass.h directory)			     #
#           $(MCFSxDIR) = the directory where the source is		     #
#									     #
#   Output: $(MCFSxOBJ) = the final object(s) / library			     #
#           $(MCFSxLIB) = external libreries + -L<libdirs>		     #
#           $(MCFSxH)   = the .h files to include			     #
#           $(MCFSxINC) = the -I$(include directories)                       #
#									     #
#                                VERSION 2.00				     #
#                	        13 - 11 - 2018				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFSxDIR)*.o $(MCFSxDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFSxOBJ = $(MCFSxDIR)MCFSimplex.o
MCFSxLIB =
MCFSxH   = $(MCFSxDIR)MCFSimplex.h
MCFSxINC = -I$(MCFSxDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCFSxOBJ): $(MCFSxDIR)MCFSimplex.C $(MCFSxH) $(MCFClH)
	$(CC) -c $*.C -o $@ $(MCFClINC) $(SW)

############################# End of makefile ################################
