#ROMS-BEC analytical version for testing

Uses analytically set (no forcing/grid files) grid, forcing, and 
boundary conditions. 

FLAT\_BOTTOM is the current active cpp key.
Use command 'grep FLAT\_BOTTOM \*' to see FLAT\_BOTTOM configurations.

MPI is deactivated in set\_global\_definitions.h because this can run on 1 core.

cppdefs.h controls the active cpp keys in all scripts.

param.h sets the grid size (currently 10x10x10).

ana\_grid.F sets the domain size and resolution (currently 10 m), land masses
(MASKING) and bathymetry.

ana\_init.F sets the initial conditions of the water column.

analytical.F adds any point sources if ANA\_PSOURCE is defined.

