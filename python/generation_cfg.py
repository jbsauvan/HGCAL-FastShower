
import math as m
from geometry_cfg import geometry_cell_side as a

generation_energy = 100.
generation_fluctuation = True

# nbr hits per GeV, used only if fluctuation is false
generation_number_of_hits_per_gev = 1000
# mip for 200 microns Si
generation_mip_energy = 0.000055
# sampling fraction in mips
# taken from 600mips in layer max for 75 GeV e-
generation_sampling = 0.0055
# electronic noise
generation_noise = True
# noise in mips
generation_noise_sigma = 0.3

# Incident position and angle
# below we shoot at center of the (3,4) cell in the parameterised triangular geometry
incident_i = 3
incident_j = 4
generation_incident_x = incident_i*a/2.+incident_j*a/2.
generation_incident_y = incident_j*a*m.sqrt(3.)+(incident_i+1)%2*a*m.sqrt(3.)/6.
generation_incident_eta = 2.


