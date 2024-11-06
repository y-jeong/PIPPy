FC = gfortran
FCFLAGS = -g -free-form -O3 

CPP = cpp -P -traditional
CPPFLAGS =

LD = $(FC)
LIB =

BINROOT=./
EX=pip.x
SHELL=/bin/sh

ARCH ?= polaris
include arch/$(ARCH)

$(info ARCH is [${ARCH}])
$(info PIPPY_ROOT_DIR is [${PIPPY_ROOT_DIR}])

# -- all src files
SRC = $(wildcard *.f)
MOD = $(wildcard *.mod)
OBJ = $(SRC:.f=.o)

#
# -- target : 	Dependencies
# --		Rule to create target

$(EX): 	$(OBJ)
	$(LD) $(LDFLAGS) -o $@ $(FCFLAGS) $(OBJ) $(LIB)

####################################################################

$(OBJ): %.o: %.f
	$(FC) $(FCFLAGS) -c $<

#
# -- Remove *.o and *~ from the directory
clean:
	rm -f *.o *.mod *~
#
# -- Remove *.o, *~, and executable from the directory
realclean:
	rm -f *.o *.mod *~ ./$(EX)

#
# -- Simple dependencies

$(OBJ) : $(MOD)

