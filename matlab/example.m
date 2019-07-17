% Example
% Geophysical parameters for 20% Quartz, 20% low Albite, 30% Forsterite, and 30% Fayalite
% and 0.1, 0.3, and 0.5 MPa, and 300, 400, and 500 K
Comp.Min = [1 5 12 13];
Comp.Fr = [0.2 0.2 0.3 0.3];
Pressure = [0.1e6 0.3e6 0.5e6];
Temperature = [300 400 500];
MV = MinVel(Comp, Pressure, Temperature)
