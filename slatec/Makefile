#=============================================================================
#
#  Makefile for stag (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 12-08-99
#
#=============================================================================
NAME    = slatec.a
DEBUG   = -O2
FLAGS   = -cpp -std=legacy $(DEBUG)
OFLAGS  = $(DEBUG) -o $(NAME)
LIB     =
FC      = gfortran 
F77     = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 .f
#
#  All objects are listed here
#
SRCS = $(wildcard *.f)
OBJS = $(SRCS:.f=.o)

$(NAME): $(OBJS)
	$(AR) -cr $(NAME) $(OBJS)

clean:
	$(RM) -f *.o *.mod *.a

.f90.o:
	$(FC) $(FLAGS) -c $*.f90

.f.o:
	$(F77) $(FLAGS) -c $*.f

.c.o:
	$(CC) -O2 -c $*.c
