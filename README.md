# Fast shower integration tool for HGCAL geometry studies


## Installation
### Standalone installation
The program can be used in a standalone way, but it needs a few dependencies installed on the system:
```
g++ version >=4.8
python version 2.7
boost_python
ROOT
```
First installation
```bash
git clone https://github.com/LLRCMS/HGCAL-FastShower.git FastShower
cd FastShower
# setup environment variables
source setupenv
# compile and create executable
make
```
Once it is installed and compiled, it is only needed to setup the environment each time:
```bash
source setupenv
```

### CMSSW installation
The program can also be installed as a CMSSW package.
```bash
# First install a recent CMSSW release
cmsenv
git clone https://github.com/LLRCMS/HGCAL-FastShower.git HGCalSimulation/FastShower
scram b
```

## Running and configuration
After the installation steps an executable is produced, `shower_simulation.exe`. It takes one argument, which is a python file containing parameter definitions.  
 A test file is available in the `python` directory.
```bash
shower_simulation.exe python/test_cfg.py
```

The configuration is written in python, generally in several files containing various groups of parameters. The parameters can be separated in the following groups:  

* General, for general parameters
* Geometry, for the geometry definition
* Generation, to configure the particle generation
* Shower, to parametrize the shower shape
* Display, to configure the display of events

The list of parameters in each group is as follows.

### General parameters

| Name | Values | Definition |  
| -------- | -------- | ------------- |  
| `events` | `int` |  Number of events to generate |  
| `debug` | `bool` |  Flag to activate the debug mode |  
| `output_file` | `string` |  Name of the output file |  

### Geometry parameters 

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `geometry_type` | `string` (`Hexagons`, `Triangles` or `External`) | Type of the geometry |  
| `geometry_file` | `string` | Input JSON file for external geometries (`geometry_type = External`) |  
| `geometry_layer` | `int in [-1,28]` | Index of the simulated layer. `-1` will sum all the layers |  
|  `geometry_layers_z` | `list(float)` |  z position of the layers (cm) |  
| `geometry_cell_side` | `float` | Size parameter of the cells (cm) (`geometry_type != External`) |
| `geometry_eta/phi_min/max` | `float` | Simulated (eta,phi) window | 

### Generation parameters

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `generation_energy` |  `float` | Energy of the incident particle (GeV) |  
| `generation_incident_eta/phi` | `float`  | (eta,phi) of the incident particle |  
| `generation_fluctuation` | `bool`  | Flag to activate fluctuations. If `False`, the mean energy profile is generated |  
| `generation_number_of_hits_per_gev` | `int` | Number of generated hits / GeV (`generation_fluctuation = False`)  |  
| `generation_mip_energy` | `float`  | MIP energy value (GeV) |  
| `generation_sampling` | `float`  | Sampling fraction |  
| `generation_noise` | `bool` | Flag to activate electronic noise |  
| `generation_noise_sigma` | `float`  | Noise in MIPs (`generation_noise = True`) |  

### Shower parameters

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `shower_moliere_radius` | `float` |  Moliere radius (cm)|  
| `shower_radiation_length` | `float`  |  Radiation length (cm) |  
| `shower_critical_energy` | `float`  | Critical energy (MeV) |  
| `shower_layers_energy ` | `dict(float)`  | Mean layer energy profile |  
| `shower_alpha` |  `float` | Sampling parameter. Determines the number of hits / GeV (`generation_fluctuation = True`) |  
| `shower_transverse_parameters` | `dict(float)`  | Parameters for the transverse profile (parabolic function)  |   
| `shower_longitudinal_parameters` | `dict(float)`  | Parameters for the longitudinal profile |  

### Display parameters

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `display_events` | `int` |  Number of events to display |  

## Predefined configurations and utilities

### Predefined geometries
Configuration files with predefined geometries are available:

* `geometry_hexagon_large_cfg`, for hexagonal cells with the same size as large cells in CMSSW
* `geometry_hexagon_small_cfg`, for hexagonal cells with the same size as small cells in CMSSW
* `geometry_triangle_large_cfg`, for triangle cells with the same area as large cells in CMSSW
* `geometry_triangle_small_cfg`, for triangle cells with the same area as small cells in CMSSW

### Particle position utility
The position of the generated particle is given in (eta,phi), which requires some calculation if ones wants to generate the particle at a particular position with respect to the cell grid.   
Three functions are available to retrieve the (eta,phi) at the center, vertex or edge of a cell. They take an initial (eta,phi) as input and returns the (eta,phi) at the center of the closest cell, or at a given vertex or edge of the closest cell. They are defined in the file `python/generation_utils.py`. Currently they can only be used for the `Triangles` and `Hexagons` geometries, but not for external geometries.

`shoot_cell_center` will return the (eta,phi) at the center of the closest cell.  
```python
def shoot_cell_center(
	eta, phi, # initial (eta,phi)
	eta_min, eta_max, # Geometry window
	phi_min, phi_max,  # ...
	z, # z of the layer
	cell_side, # cell side
	type # type of geometry
)
```
`shoot_cell_vertex` will return the (eta,phi) at a specified vertex of the closest cell.   
```python
def shoot_cell_vertex(
 	eta, phi, # initial (eta,phi)
	vertex_number, # vertex index
	eta_min, eta_max, # Geometry window
	phi_min, phi_max,  # ...
	z, # z of the layer
	cell_side, # cell side
	type # type of geometry
)
```
`shoot_cell_edge` will return the (eta,phi) at the middle of a specified edge of the closest cell.   
```python
def shoot_cell_edge(
 	eta, phi, # initial (eta,phi)
	edge_number, # edge index
	eta_min, eta_max, # Geometry window
	phi_min, phi_max,  # ...
	z, # z of the layer
	cell_side, # cell side
	type # type of geometry
)
``` 
 
 

