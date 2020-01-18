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
#FLAGS   = -cpp -freal-4-real-8 -fdefault-real-8 -DR8 $(DEBUG)
FLAGS   = -cpp -fdefault-real-8 -DR8 $(DEBUG)
OFLAGS  = $(DEBUG) 
LIB     = -L$(HOME)/local/OpenBLAS/lib -lopenblas slatec/slatec.a
COMP    = gfortran 
F77     = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All objects are listed here
#
OBJS = stag.o getver.o

all: $(NAME) mkgrid d2s

$(NAME): $(OBJS)
	$(COMP) $(OFLAGS) $(OBJS) $(LIB) -o $(NAME)

mkgrid: mkgrid.o
	$(COMP) $(FLAGS) $(DEBUG) -o mkgrid mkgrid.o

d2s: d2s.o
	$(COMP) $(FLAGS) $(DEBUG) -o d2s d2s.o

clean:
	/bin/rm -f *.o *.mod $(NAME) mkgrid d2s

.f90.o:
	$(COMP) $(FLAGS) -c $*.f90

.f.o:
	$(F77) $(FLAGS) -c $*.f

.c.o:
	$(CC) -O2 -c $*.c
