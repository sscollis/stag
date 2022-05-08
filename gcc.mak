#=============================================================================
#
#  Makefile for stag (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 12-08-99
#
#=============================================================================
NAME    = stag
DEBUG   = -O2 -fopenmp
#FLAGS  = -cpp -freal-4-real-8 -fdefault-real-8 -DR8 $(DEBUG)
FLAGS   = -cpp -fdefault-real-8 -DR8 $(DEBUG)
OFLAGS  = $(DEBUG) 
ifeq ($(OPENBLAS_DIR),)
  OPENBLAS_DIR = /usr/local/opt/openblas
endif
LIB     = -L$(OPENBLAS_DIR)/lib -lopenblas slatec/slatec.a
FC      = gfortran 
F77     = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All objects are listed here
#
OBJS = stag.o getver.o

all: $(NAME) slatec.a mkgrid d2s
#
# SSC:  Versions 1-3 and 5 not yet updated
#
#all: $(NAME) slatec.a mkgrid d2s stag_v1 stag_v2 stag_v3 stag_v4 \
#stag_v5 stag_v6

$(NAME): $(OBJS) slatec.a
	$(FC) $(OFLAGS) $(OBJS) $(LIB) -o $(NAME)

slatec.a:
	cd slatec && make F77=$(F77) FC=$(FC)

stag_v1: stag_v1.o
	$(FC) $(OFLAGS) stag_v1.o getver.o $(LIB) -o stag_v1

stag_v2: stag_v2.o
	$(FC) $(OFLAGS) stag_v2.o getver.o $(LIB) -o stag_v2

stag_v3: stag_v3.o
	$(FC) $(OFLAGS) stag_v3.o getver.o $(LIB) -o stag_v3

stag_v4: stag_v4.o
	$(FC) $(OFLAGS) stag_v4.o getver.o $(LIB) -o stag_v4

stag_v5: stag_v5.o
	$(FC) $(OFLAGS) stag_v5.o getver.o $(LIB) -o stag_v5

stag_v6: stag_v6.o
	$(FC) $(OFLAGS) stag_v6.o getver.o $(LIB) -o stag_v6

mkgrid: mkgrid.o
	$(FC) $(FLAGS) $(DEBUG) -o mkgrid mkgrid.o

d2s: d2s.o
	$(FC) $(FLAGS) $(DEBUG) -o d2s d2s.o

clean:
	$(RM) -f *.o *.mod $(NAME) stag_v4 stag_v6 mkgrid d2s
	cd slatec && make clean

.f90.o:
	$(FC) $(FLAGS) -c $*.f90

.f.o:
	$(F77) $(FLAGS) -c $*.f

.c.o:
	$(CC) -O2 -c $*.c
