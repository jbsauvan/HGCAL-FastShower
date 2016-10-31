import math as m
from ROOT import TVector2

# FIXME: the cell side can be used as a global multiplicative factor

def vertices(type, i, a):
  if type is 'Triangles':
    if i%2: 
      return [(a/2., a/(2.*m.sqrt(3.))),
              (-a/2., a/(2.*m.sqrt(3.))),
              (0., -a/m.sqrt(3.))
             ]
    else:
      return [(a/2.,-a/(2.*m.sqrt(3.))),
              (0., a/m.sqrt(3.)),
              (-a/2.,-a/(2.*m.sqrt(3.)))
             ]

  elif type is 'Hexagons':
      return [(a*m.sqrt(3)/2., -a/2.),
              (a*m.sqrt(3)/2., a/2.),
              (0., a),
              (-a*m.sqrt(3)/2., a/2.),
              (-a*m.sqrt(3)/2., -a/2.),
              (0., -a)
             ]
  return []

# find the center of the cell closest to (eta,phi)
# it assumes that the cell 0 is located at (eta_max,phi_min)
def shoot_center_cell(eta, phi, eta_min, eta_max, phi_min, phi_max, z, cell_side, type):
  if eta<eta_min or eta>eta_max or\
     TVector2.Phi_mpi_pi(phi-phi_min)<0 or TVector2.Phi_mpi_pi(phi-phi_max)>0:
    raise Exception("Error: trying to shoot particle outside geometry window")
  theta0 = 2.*m.atan(m.exp(-eta_max))
  theta = 2.*m.atan(m.exp(-eta))
  r0 = z*m.tan(theta0)
  r = z*m.tan(theta)
  x0 = r0*m.cos(phi_min)
  y0 = r0*m.sin(phi_min)
  x = r*m.cos(phi)
  y = r*m.sin(phi)
  if type is 'Triangles':
    dxdi = cell_side/2.
    dxdj = cell_side/2.
    dydj = cell_side*m.sqrt(3.)/2.
  else:
    dxdi = cell_side*m.sqrt(3.)
    dxdj = cell_side*m.sqrt(3.)/2.
    dydj = cell_side*3./2.
  j = round((y-y0)/dydj)
  i = round((x-x0 - j*dxdj)/dxdi)
  x_center = x0 + i*dxdi + j*dxdj
  y_center = y0 + j*dydj
  if type is 'Triangles' and int(i)%2: y_center += cell_side*m.sqrt(3)/6.
  r_center = m.sqrt(x_center**2 + y_center**2)
  theta_center = m.atan(r_center/z)
  eta_center = -m.log(m.tan(theta_center/2.))
  phi_center = m.copysign(m.acos(x_center/r_center), y_center)
  return eta_center, phi_center


# FIXME: avoid duplicating code
# find vertices of the cell closest to (eta,phi)
# it assumes that the cell 0 is located at (eta_max,phi_min)
def shoot_vertex_cell(eta, phi, vertex_number, 
                      eta_min, eta_max, phi_min, phi_max, z, cell_side, type):
  if eta<eta_min or eta>eta_max or\
     TVector2.Phi_mpi_pi(phi-phi_min)<0 or TVector2.Phi_mpi_pi(phi-phi_max)>0:
    raise Exception("Error: trying to shoot particle outside geometry window")
  # Find center
  theta0 = 2.*m.atan(m.exp(-eta_max))
  theta = 2.*m.atan(m.exp(-eta))
  r0 = z*m.tan(theta0)
  r = z*m.tan(theta)
  x0 = r0*m.cos(phi_min)
  y0 = r0*m.sin(phi_min)
  x = r*m.cos(phi)
  y = r*m.sin(phi)
  if type is 'Triangles':
    dxdi = cell_side/2.
    dxdj = cell_side/2.
    dydj = cell_side*m.sqrt(3.)/2.
  else:
    dxdi = cell_side*m.sqrt(3.)
    dxdj = cell_side*m.sqrt(3.)/2.
    dydj = cell_side*3./2.
  j = round((y-y0)/dydj)
  i = round((x-x0 - j*dxdj)/dxdi)
  x_center = x0 + i*dxdi + j*dxdj
  y_center = y0 + j*dydj
  if type is 'Triangles' and int(i)%2: y_center += cell_side*m.sqrt(3)/6.
  # Get position of vertices
  vertex = vertices(type, int(i), cell_side)
  if vertex_number>=len(vertex):
    raise Exception('Vertex {0} doesn\'t exist for type {1}'.format(vertex_number, type))
  x_vertex, y_vertex = vertex[vertex_number]
  x_vertex += x_center
  y_vertex += y_center
  r_vertex = m.sqrt(x_vertex**2 + y_vertex**2)
  theta_vertex = m.atan(r_vertex/z)
  eta_vertex = -m.log(m.tan(theta_vertex/2.))
  phi_vertex = m.copysign(m.acos(x_vertex/r_vertex), y_vertex)
  return eta_vertex, phi_vertex


