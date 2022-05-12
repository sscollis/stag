# STAG: Staggered-grid 2d Incompressible Navier-Stokes Solver

Solves the two-dimensional incompressible Navier-Stokes 
equations on a staggered mesh system.

I have updated this version to build on MacOS and Linux with GCC gfortran.

## Build

To build using GCC on Mac or Linux:

    ln -s gcc.mak Makefile
    make
    cd test
    ../mkgrid
      128 128
      0 0
    ../stag < test.inp 

### Notes
1. Uses OpenBLAS which you can install using homebrew on MacOS
2. Automatically builds a subset of SLATEC that is needed
3. Visualize the results using Paraview -- STAG generates both a `grid.xyz`
   and `output.q.*` files in Plot3d format.  

## Vortex test case

The `test` directory contains a simple vortex rebound problem used to
test STAG.  The `README.md` file therein describes the test case and how
to run and visualize it. 

## Versions 

This code was developed primarily for educational purposes and comes
in several versions that demonstrate different capabilities and levels
of complexity.

Version | Notes
--------|------------------------------------------------------------------
  3     | Stable solver	
  4     | Multiple pressure solvers
  5     | REAL defines added to version 3
  6     | Immersed boundary conditions added to version 3

---

<p align=center>
<img src=https://github.com/sscollis/stag/blob/master/omega.gif>
<br>Vortex rebound, contours of vorticity.</p>

---

S. Scott Collis
