##############################################################################
################################## makefile ##################################
##############################################################################
#									     #
#   makefile of MCFClone 						     #
#									     #
#   Input:  $(MCFCeDIR) = the directory where the source is                  #
#									     #
#   Output: $(MCFCeH)   = the .h files to include			     #
#           $(MCFCeINC) = the -I$(MCFClone directory)			     #
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
	rm -f $(MCFCeDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFCeH   = $(MCFCeDIR)MCFClone.h
MCFCeINC = -I$(MCFCeDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

############################# End of makefile ################################
