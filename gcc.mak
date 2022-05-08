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

all: $(NAME) slatec.a mkgrid d2s
#
# SSC:  Versions 1-3 and 5 not yet updated
#
#all: $(NAME) slatec.a mkgrid d2s stag_v1 stag_v2 stag_v3 stag_v4 \
#stag_v5 stag_v6

$(NAME): $(OBJS) slatec.a
	$(COMP) $(OFLAGS) $(OBJS) $(LIB) -o $(NAME)

slatec.a:
	cd slatec && make

stag_v1: stag_v1.o
	$(COMP) $(OFLAGS) stag_v1.o getver.o $(LIB) -o stag_v1

stag_v2: stag_v2.o
	$(COMP) $(OFLAGS) stag_v2.o getver.o $(LIB) -o stag_v2

stag_v3: stag_v3.o
	$(COMP) $(OFLAGS) stag_v3.o getver.o $(LIB) -o stag_v3

stag_v4: stag_v4.o
	$(COMP) $(OFLAGS) stag_v4.o getver.o $(LIB) -o stag_v4

stag_v5: stag_v5.o
	$(COMP) $(OFLAGS) stag_v5.o getver.o $(LIB) -o stag_v5

stag_v6: stag_v6.o
	$(COMP) $(OFLAGS) stag_v6.o getver.o $(LIB) -o stag_v6

mkgrid: mkgrid.o
	$(COMP) $(FLAGS) $(DEBUG) -o mkgrid mkgrid.o

d2s: d2s.o
	$(COMP) $(FLAGS) $(DEBUG) -o d2s d2s.o

clean:
	$(RM) -f *.o *.mod $(NAME) stag_v4 stag_v6 mkgrid d2s
	cd slatec && make clean

.f90.o:
	$(COMP) $(FLAGS) -c $*.f90

.f.o:
	$(F77) $(FLAGS) -c $*.f

.c.o:
	$(CC) -O2 -c $*.c
