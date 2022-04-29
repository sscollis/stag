# Stag:  Staggered0grid 2d Incompressible Navier--Stokes

Solves the two-dimensional incompressible Navier-Stokes 
equations on a staggered mesh system.

I have updated this version to build on Darwin with gfortran

## Build

    make
    cd test
    ../mkgrid
      128 128
      0 0
    ../stag < test.inp 

Visualize the results using Paraview.

## Notes

This code was developed primarily for educational purposes and comes
in several versions that demonstrate different capabilities and levels
of complexity.

Version | Notes
--------|-------------------------------------
  3     | Stable solver	
  4     | Multiple pressure solvers
  5     | REAL defines added to version 3
  6     | Immersed boundary conditions added to 3

S. Scott Collis\
Sat Jan 18 07:55:57 MST 2020
