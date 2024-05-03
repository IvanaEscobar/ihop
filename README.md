%%
%%

Package ''ihop'' is a ray-trace code, for calculating acoustics.
The package should interact with the MITgcm kernel on various levels 
(initialisation, time-stepping, post-processing).

Development is up to date with MITgcm checkpoint68s

It starts with taking hydrography at a single time-step along a track and 
calculates arrival times. Lastly, it outputs the travel time and ray angle as
diagnostics.


ADD useIHOP to PARAMS.h in the code modifications to use this package
