
import math as m
from geometry_cfg import geometry_side_length as a

generation_energy = 100.
generation_fluctuation = False
generation_number_of_hits_per_gev = 1000
generation_alpha = 0.25
generation_mip_energy = 0.000055
generation_layers_energy = [40.0,69.8,119.6,178.9,248.8,315.1,382.0,431.6,477.7,
                              498.7,533.6,514.8,490.0,435.1,386.7,325.4,277.9,224.4,186.5,
                              145.3,108.7,73.7,52.1,33.0,22.5,13.1,8.6,4.8]
incident_i = 3
incident_j = 4
generation_incident_x = incident_i*a/2.+incident_j*a/2.
generation_incident_y = incident_j*a*m.sqrt(3.)+(incident_i+1)%2*a*m.sqrt(3.)/6.
generation_incident_eta = 2.

generation_r0layer15 = 2.3/m.log(10.)
generation_shower_transverse_parameters = [9.-(18./63.), 135./630., 45./630.]
