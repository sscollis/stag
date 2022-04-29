# Basic Stag run

## Make the grid

    ../mkgrid
    128 128
    0 0

## Run the code

    ../stag < ic.inp
    ../stag < test.inp

# Visualize

Use Paraview to visualize the results with the gridfile `grid.xyz`
and Q-files `output.q.*`

Note that the Q-file variables are:

    p = pressure
    u = velocity in x-direction
    v = velocity in y-direction
    omega = vorticity in the z-direction
    ke = kinetic energy

so exercise caution in using computing qualtities that assume
default q-file inputs.

S. Scott Collis
