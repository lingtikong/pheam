.SUFFIXES : .o .cpp
# compiler and flags
CC     = g++ -Wno-unused-result
#CC     = icc
LINK   = $(CC) -static
CFLAGS = -O3 $(DEBUG)

#
OFLAGS = -O3 $(DEBUG)
INC    = $(FFTINC) $(LPKINC) $(USRINC) $(SPGINC) $(GSLINC) $(OPMINC)
LIB    = $(FFTLIB) $(LPKLIB) $(USRLIB) $(SPGLIB) $(GSLLIB) $(OPMLIB) $(SYSLIB)

# icc specific
#SYSLIB = -lstdc++ -lpthread -lguide

# fftw 3 library
#FFTINC    = -I/opt/fftw/fftw3/include
#FFTLIB    = -L/opt/fftw/fftw3/lib -lfftw3

# parallelize part of the code via OpenMP
# g++
OPMINC = -fopenmp -DOMP
OPMLIB = -lgomp -lpthread
# icc
#OPMINC = -openmp -parallel -fast -DOMP
#OPMLIB = -openmp -parallel -fast

# Lapack library
# cLapack
LPKINC = -I/opt/clapack/3.2.1/include
LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm
# MKL Lapack
LPKINC = -I/opt/intel/cmkl/10.2.5.035/include -DMKL
LPKLIB = -L/opt/intel/cmkl/10.2.5.035/lib/em64t -lmkl_lapack95_lp64 -lmkl_intel_lp64 -Wl,--start-group -lmkl_intel_thread  -lmkl_core -Wl,--end-group -liomp5 -lpthread

# spglib, used to get the irreducible q-points
SPGINC = -I/opt/spglib/0.7.1/include
SPGLIB = -L/opt/spglib/0.7.1/lib -lsymspg

# gsllib, used to get the spline interpolation
GSLINC = -I/opt/gsl/include
GSLLIB = -L/opt/gsl/lib -lgslcblas -lgsl

# Debug flags
#DEBUG = -g -DDEBUG
#====================================================================
# executable name
ROOT   = pheam
EXE    = $(ROOT)
#====================================================================
# source and rules
SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

#====================================================================
all:  ${EXE}

${EXE}: $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod

tar:
	rm -f ${ROOT}.tar; tar -czvf ${ROOT}.tar.gz *.cpp  *.h Makefile README

.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
