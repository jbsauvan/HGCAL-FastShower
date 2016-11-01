
# FIXME: find a better way to import
try:
    # Standalone imports
    from geometry_cfg import *
    from shower_cfg import *
    from generation_cfg import *
    from display_cfg import *
    from generation_utils import shoot_cell_center, shoot_cell_vertex, shoot_cell_edge
except ImportError:
    # CMSSW imports
    from HGCalSimulation.FastShower.geometry_cfg import *
    from HGCalSimulation.FastShower.shower_cfg import *
    from HGCalSimulation.FastShower.generation_cfg import *
    from HGCalSimulation.FastShower.display_cfg import *
    from HGCalSimulation.FastShower.generation_utils import shoot_cell_center, shoot_cell_vertex, shoot_cell_edge

events = 1
debug = False
output_file = 'test.root'

geometry_type = 'Triangles'
generation_incident_eta, generation_incident_phi = \
shoot_cell_edge(eta=2.0, phi=0., edge_number=0,
                  # geometry window
                  eta_min=geometry_eta_min, eta_max=geometry_eta_max, 
                  phi_min=geometry_phi_min, phi_max=geometry_phi_max, 
                  z=geometry_layers_z[geometry_layer if geometry_layer!=-1 else 0],
                  cell_side=geometry_cell_side, 
                  type=geometry_type
                 )
