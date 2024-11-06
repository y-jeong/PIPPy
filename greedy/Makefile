MAGMADIR=/soft/libraries/math_libs/magma-2.8.0/PrgEnv-gnu
CUDADIR=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/cuda

FC         = nvfortran
FCFLAGS    = -g -O3 -Wall

CC         = nvc
CFLAGS     = -Wall -I$(MAGMADIR)/include -I$(CUDADIR)/include
#CFLAGS     = -Wall

#CPP = cpp -P -traditional
#CPPFLAGS =

#LD         = $(FC)
LIBS        = -L$(MAGMADIR)/lib -lmagma -Wl,-rpath,$(MAGMADIR)/lib
#LIBS        =  -L$(MAGMADIR)/lib -lmagma

# BLAS and LAPACK libraries
LIBS       += -L$(NVIDIA_PATH)/compilers/lib -lblas -llapack
#LDFLAGS    = -lblas -llapack
LDFLAGS    =

BINROOT=./
EX=forward.x
SHELL=/bin/sh

ARCH ?= polaris
include arch/$(ARCH)

$(info ARCH is [${ARCH}])

# -- all src files
MOD = $(wildcard *.mod)
OBJ = \
	forward.o		\
	offset.o		\
	magma2.o		\
	magma2_sfortran.o	\
	magma2_dfortran.o	\
	magma2_cfortran.o	\
	magma2_zfortran.o	\
	magma2_common.o		\

all: $(EX)

# generate precisions
magma2_sfortran.F90: magma2_zfortran.stamp
magma2_dfortran.F90: magma2_zfortran.stamp
magma2_cfortran.F90: magma2_zfortran.stamp
magma2_zfortran.stamp: magma2_zfortran.F90
#	$(MAGMADIR)/tools/codegen.py magma2_zfortran.F90
	./codegen.py magma2_zfortran.F90
	touch $@

# force modules to compile in correct order
forward.o: magma2.o
magma2.o: magma2_sfortran.o magma2_dfortran.o magma2_cfortran.o magma2_zfortran.o
magma2_sfortran.o: magma2_common.o
magma2_dfortran.o: magma2_common.o
magma2_cfortran.o: magma2_common.o
magma2_zfortran.o: magma2_common.o
#
# -- target : 	Dependencies
# --		Rule to create target

#$(EX): 	$(OBJ)
#	$(LD) $(LDFLAGS) -o $@ $(FCFLAGS) $(OBJ) $(LIB)

$(EX): $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $^ \
		$(LIBS)

####################################################################

#$(OBJ): %.o: %.F90
#	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

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

