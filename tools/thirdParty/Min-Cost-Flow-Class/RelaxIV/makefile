##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of RelaxIV 						     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(MCFClH)   = the include files for MCFClass		     #
#           $(MCFClINC) = the -I$(MCFClass.h directory)			     #
#           $(MCFR4DIR) = the directory where the source is		     #
#									     #
#   Output: $(MCFR4OBJ) = the final object(s) / library			     #
#           $(MCFR4LIB) = external libreries + -L<libdirs>		     #
#           $(MCFR4H)   = the .h files to include			     #
#           $(MCFR4INC) = the -I$(include directories)                       #
#									     #
#                                VERSION 4.00				     #
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
	rm -f $(MCFR4DIR)*.o $(MCFR4DIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFR4OBJ = $(MCFR4DIR)RelaxIV.o
MCFR4LIB =
MCFR4H   = $(MCFR4DIR)RelaxIV.h
MCFR4INC = -I$(MCFR4DIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCFR4OBJ): $(MCFR4DIR)RelaxIV.C $(MCFR4H) $(MCFClH)
	$(CC) -c $*.C -o $@ $(MCFClINC) $(SW)

############################# End of makefile ################################
