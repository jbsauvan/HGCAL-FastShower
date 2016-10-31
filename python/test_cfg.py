
# FIXME: find a better way to import
try:
    # Standalone imports
    from geometry_cfg import *
    from shower_cfg import *
    from generation_cfg import *
    from display_cfg import *
    from generation_utils import shoot_center_cell, shoot_vertex_cell
except ImportError:
    # CMSSW imports
    from HGCalSimulation.FastShower.geometry_cfg import *
    from HGCalSimulation.FastShower.shower_cfg import *
    from HGCalSimulation.FastShower.generation_cfg import *
    from HGCalSimulation.FastShower.display_cfg import *

events = 1
debug = False
output_file = 'test.root'

geometry_type = 'Triangles'
generation_incident_eta, generation_incident_phi = \
shoot_center_cell(2.0, 0.,
                  geometry_eta_min, geometry_eta_max, 
                  geometry_phi_min, geometry_phi_max, 
                  geometry_layers_z[geometry_layer if geometry_layer!=-1 else 0],
                  geometry_cell_side, 
                  geometry_type
                 )
