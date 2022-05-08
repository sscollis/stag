# Stag:  Staggered-grid 2d Incompressible Navier--Stokes

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

1. Uses OpenBLAS which you can install using homebrew on MacOS
2. Automatically builds a subset of SLATEC that is needed

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

---

<p align=center>
<img src=https://github.com/sscollis/stag/blob/master/omega.gif>
<br>Vortex rebound, contours of vorticity.</p>

---

S. Scott Collis\
Sat Jan 18 07:55:57 MST 2020
