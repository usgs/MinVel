# Example 1, values on the command line
# Geophysical parameters for 20% Quartz, 20% low Albite, 30% Forsterite, and 30% Fayalite at 
# 300, 400, and 500K and 0.1, 0.3, and 0.5 MPa
python MinVelWrapper.py -t 300,400,500 -p 0.1e6,0.3e6,0.5e6 -cm 1,5,12,13 -cv 0.2,0.2,0.3,0.3

# Example 2, values in a file
# Geophysical parameters for 20% Quartz, 20% low Albite, 30% Forsterite, and 30% Fayalite at 
# 300, 400, and 500K and 0.1, 0.3, and 0.5 MPa
python MinVelWrapper.py -f example.Input
