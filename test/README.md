# Basic STAG test case:  Vortex Rebound 

## Make the grid

    ./mk-grid.sh

## Run the code

    ../stag < test.inp

## Visualize

Use Paraview to visualize the results with the gridfile `grid.xyz`
and Q-files `output.q.*`

Note that the Q-file variables are:

    p     = pressure
    u     = velocity in x-direction
    v     = velocity in y-direction
    omega = vorticity in the z-direction
    ke    = kinetic energy

so exercise caution in using computing qualtities that assume
default Q-file inputs.

S. Scott Collis\
