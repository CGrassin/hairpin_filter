# Hairpin filter

This repository contains the source-code for my hairpin filter generator and openEMS model.

Hosted here, with technical details and actual measurements: https://charleslabs.fr/en/project-Hairpin+filter+design

License: MIT

## openEMS simulation

To run the openEMS simulation, you need to install openEMS through octave or MATLAB, and Paraview (optional).
You can change the simulation parameters in the setup block. Run the file to launch the openEMS analysis (may take a while!).

The ouput is a S-parameter plot and the "vtp" and "vtr" files to do more post-processing or display in Paraview, for instance.