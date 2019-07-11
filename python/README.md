# MinVel Python routines
This program is used to calculate anharmonic p- and s-wave velocity and density for zero-porosity mineral aggregates. It is based on the work of Hacker and Abers (2004, updated in 2016) with additional minerals and optimized for related work. Running the code requires a mineral physics database, MineralPhysicsDatabase.nc. Combined with a petrologic database (Sowers and Boyd, 2019) and thermal and geologic models, results from this code are used in as part of an effort to produce a three dimensional national crustal model for use with seismic hazard studies. (https://www.sciencebase.gov/catalog/item/5b046c6ce4b0d8682b96353a).

Installation:
The directory in which MinVel.py and MinVelWrapper.py exist must be in the PYTHONPATH.
The program requires the numpy and netCDF4 modules. To install them, run 'make init'.
To test the program, run 'make test'

Please see or query the code for help (e.g. prompt>> python3 MinVelWrapper -h).

Any questions regarding the code can be directed to olboyd@usgs.gov.

References:

Abers, G.A., and Hacker, B.R., 2016, A MATLAB toolbox and Excel workbook for calculating the densities, seismic wave speeds, and major element composition of minerals and rocks at pressure and temperature: Geochem. Geophys. Geosys., v. 16, p. 616–624.

Hacker, B.R., and Abers, G.A., 2004, Subduction Factory 3: An Excel worksheet and macro for calculating the densities, seismic wave speeds, and H2O contents of minerals and rocks at pressure and temperature: Geochem. Geophys. Geosys., v. 5, no. Q01005.

Sowers, T., and Boyd, O.S., 2019, Petrologic and mineral physics database for use with the U.S. Geological Survey National Crustal Model: U.S. Geological Survey Open-File Report 2019–1035, 17 p., https://doi.org/10.3133/ofr20191035.
