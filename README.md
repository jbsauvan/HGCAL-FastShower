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
### Configuration
The configuration is written in python, generally in several files containing various groups of parameters. The parameters can be separated in the following groups:  

* General, for general parameters
* Geometry, for the geometry definition
* Generation, to configure the particle generation
* Shower, to parametrize the shower shape
* Display, to configure the display of events

The list of parameters in each group is as follows.

#### General parameters

| Name | Values | Definition |  
| -------- | -------- | ------------- |  
| `events` | `int` |  Number of events to generate |  
| `debug` | `bool` |  Flag to activate the debug mode |  
| `output_file` | `string` |  Name of the output file |  

#### Geometry parameters 

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `geometry_type` | `string` (`Hexagons`, `Triangles` or `External`) | Type of the geometry |  
| `geometry_file` | `string` | Input JSON file for external geometries (`geometry_type = External`) |  
| `geometry_layer` | `int in [-1,28]` | Index of the simulated layer. `-1` will sum all the layers |  
|  `geometry_layers_z` | `list(float)` |  z position of the layers (cm) |  
| `geometry_cell_side` | `float` | Size parameter of the cells (cm) (`geometry_type != External`) |
| `geometry_eta/phi_min/max` | `float` | Simulated (eta,phi) window | 

#### Generation parameters

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

#### Shower parameters

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `shower_moliere_radius` | `float` |  Moliere radius (cm)|  
| `shower_radiation_length` | `float`  |  Radiation length (cm) |  
| `shower_critical_energy` | `float`  | Critical energy (MeV) |  
| `shower_layers_energy ` | `dict(float)`  | Mean layer energy profile |  
| `shower_alpha` |  `float` | Sampling parameter. Determines the number of hits / GeV (`generation_fluctuation = True`) |  
| `shower_transverse_parameters` | `dict(float)`  | Parameters for the transverse profile (parabolic function)  |   
| `shower_longitudinal_parameters` | `dict(float)`  | Parameters for the longitudinal profile |  

#### Display parameters

| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `display_events` | `int` |  Number of events to display |  


