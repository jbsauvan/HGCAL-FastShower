
import math as m
from geometry_cfg import geometry_cell_side as a

generation_energy = 100.
generation_fluctuation = False

generation_number_of_hits_per_gev = 1000
generation_mip_energy = 0.000055
generation_sampling = 0.0055
generation_noise = True
generation_noise_sigma = 1.

# Incident position and angle
incident_i = 3
incident_j = 4
generation_incident_x = incident_i*a/2.+incident_j*a/2.
generation_incident_y = incident_j*a*m.sqrt(3.)+(incident_i+1)%2*a*m.sqrt(3.)/6.
generation_incident_eta = 2.


