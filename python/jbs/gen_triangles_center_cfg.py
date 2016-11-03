try:
    # Standalone imports
    from geometry_triangle_small_cfg import *
    from shower_cfg import *
    from generation_cfg import *
    from display_cfg import *
    from generation_utils import shoot_cell_center
except ImportError:
    # CMSSW imports
    from HGCalSimulation.FastShower.geometry_triangle_small_cfg import *
    from HGCalSimulation.FastShower.shower_cfg import *
    from HGCalSimulation.FastShower.generation_cfg import *
    from HGCalSimulation.FastShower.display_cfg import *
    from HGCalSimulation.FastShower.generation_utils import shoot_cell_center

events = 1000
debug = False
output_file = 'output/triangles_small_center.root'

generation_noise = False
generation_incident_eta, generation_incident_phi = \
shoot_cell_center(eta=2.0, phi=0.,
                  # geometry window
                  eta_min=geometry_eta_min, eta_max=geometry_eta_max, 
                  phi_min=geometry_phi_min, phi_max=geometry_phi_max, 
                  z=geometry_layers_z[geometry_layer if geometry_layer!=-1 else 0],
                  cell_side=geometry_cell_side, 
                  type=geometry_type
                 )
