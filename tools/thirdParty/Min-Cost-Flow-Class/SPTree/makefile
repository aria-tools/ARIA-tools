##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of SPTree 							     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(MCFClH)   = the include files for MCFClass		     #
#           $(MCFClINC) = the -I$(MCFClass.h directory)			     #
#           $(MCFSTDIR) = the directory where the source is		     #
#									     #
#   Output: $(MCFSTOBJ) = the final object(s) / library			     #
#           $(MCFSTLIB) = external libreries + -L<libdirs>		     #
#           $(MCFSTH)   = the .h files to include			     #
#           $(MCFSTINC) = the -I$(include directories)                       #
#									     #
#                                VERSION 3.00				     #
#                	        07 - 05 - 2012				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFSTDIR)*.o $(MCFSTDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFSTOBJ = $(MCFSTDIR)SPTree.o
MCFSTLIB =
MCFSTH   = $(MCFSTDIR)SPTree.h
MCFSTINC = -I$(MCFSTDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(MCFSTOBJ): $(MCFSTDIR)SPTree.C $(MCFClH)
	$(CC) -c $*.C -o $@ $(MCFClINC) $(SW)

############################# End of makefile ################################
