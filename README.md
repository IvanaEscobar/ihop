%%

Package ''ihop'' is a ray-trace code, for calculating acoustics.
The package should interact with the MITgcm kernel on various levels 
(initialisation, time-stepping, post-processing).

Development is up to date with MITgcm checkpoint68s

It starts with taking hydrography at a single time-step along a track and 
calculates arrival times. Lastly, it outputs the travel time and ray angle as
diagnostics.


ADD useIHOP to PARAMS.h in the code modifications to use this package

Package is dependent on the following packages:
- *cal* for storing times of sound transmissions
- *cost* for aggregation of acoustic cost function contributions


## Tips
For input, you will be asked to generate range points along a 2D plane between 
a source and receiver. The number of range points can vary from 2 to N. The 
position of a receiver _must_ be contained within the ranges specified, e.g. 
`ihop_rr < ihop_ranges(N) - <step size>`. In general, the step size is 10% the 
maximum ocean depth defined in your `.bty` file. It's recommended to place an 
`ihop_ranges` point at the recevier lat, lon position. 

# TO-DO
- PYTHON: add simple input file generation
- FORTRAN77: add simple verification problem
- PYTHON: add synthetic observation data file generation 
