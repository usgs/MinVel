# MinVel
This program is used to calculate anharmonic p- and s-wave velocity and density for zero-porosity mineral aggregates. It is based on the work of Hacker and Abers (2004, updated in 2016) with additional minerals and optimized for related work. Running the code requires a mineral physics database, MineralPhysicsDatabase_181107.mat. Combined with a petrologic database (Sowers and Boyd, in press) and thermal and geologic models, results from this code are used in the program RockProp.m, which accounts for rock porosity and seismic attenuation. These codes are part of an effort to produce a three dimensional national crustal model for use with seismic hazard studies. (https://www.sciencebase.gov/catalog/item/5b046c6ce4b0d8682b96353a).

Please see or query the code for help (e.g. Matlab>> help MinVel).

Any questions regarding the code can be directed to olboyd@usgs.gov.

References:

Abers, G.A., and Hacker, B.R., 2016, A MATLAB toolbox and Excel workbook for calculating the densities, seismic wave speeds, and major element composition of minerals and rocks at pressure and temperature: Geochem. Geophys. Geosys., v. 16, p. 616–624.

Hacker, B.R., and Abers, G.A., 2004, Subduction Factory 3: An Excel worksheet and macro for calculating the densities, seismic wave speeds, and H2O contents of minerals and rocks at pressure and temperature: Geochem. Geophys. Geosys., v. 5, no. Q01005.

Sowers, T. and O. S. Boyd, in press, Petrologic and Mineral Physics Database for use with the USGS National Crustal Model (ver. 1.0): U.S. Geological Survey Open File Report.

Disclaimer, once approved for release:
This software has been approved for release by the U.S. Geological Survey (USGS). Although the software has been subjected to rigorous review, the USGS reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.
