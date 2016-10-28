
import math as m

geometry_type = 'Triangles'
# layer:
# if between 0 and 27, the generated energy is weighted according to the layer profile (elayers[i])
# the incident position is also shifted according to layer z position (zlayers[i]) and incident direction (etainc)
# if set to -1, the full generated energy is deposited without longitudinal weighting and the
# incident position is not shifted, that is it corresponds to layer 0 with the full energy deposited in it
geometry_layer = -1
# layers' z positions
# from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
# the values are the silicon (centre) z positions of the 28 layers wrt HGCAL z entry position in cm
geometry_layers_z = [.765,1.515,2.745,3.495,4.725,5.475,6.705,7.455,8.685,9.435,
                       10.745,11.615,12.925,13.795,15.105,15.975,17.285,18.155,19.465,20.335,
                       21.785,22.855,24.305,25.375,26.825,27.895,29.345,30.415]
# entry point in HGCal
z0 = 320.
# absolute HGCal layer positions
geometry_layers_z[:] = [z+z0 for z in geometry_layers_z]

# multiply cell side by sqrt(6) for triangles to get equal area
geometry_cell_side = 1.190

# Define (eta,phi) window to build the geometry
geometry_eta_min = 1.9
geometry_eta_max = 2.1
geometry_phi_min = -m.pi/24. # -7.5 deg 
geometry_phi_max = m.pi/24. # +7.5 deg 
#
geometry_file = 'data/AC_11.json'
