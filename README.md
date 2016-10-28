# Fast shower integration tool for HGCAL geometry studies.


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

* General, for general parameters that doesn't fit in any of the other categories
* Geometry, for the geometry definition
* Generation, to configure the particle generation
* Shower, to parametrize the shower shape
* Display, to configure the display of events

The list of parameters in each group is the following:

* General  
| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `events` | `int` |  number of events to generate |  
| `debug` | `bool` |  flag to activate the debug mode |  
| `output_file` | `string` |  name of the output file |  

* Geometry  
| Name | Values | Definition |  
| ------ | ------ | ------------ |  
| `geometry_type` | `string` (`Hexagons`, `Triangles` or `External`) | Type of the geometry |  
| `geometry_file` | `string` | Input JSON file for external geometries (`geometry_type = External`) |  
| `geometry_layer` | `int` in `[-1,28]` | Index of the simulated layer. `-1` will sum all the layers |  
|  `geometry_layers_z` | `list(float)` |  z position of the layers |  
| `geometry_cell_side` | `float` | Size parameter of the cells (`geometry_type != External`) |
| `geometry_eta/phi_min/max` | `float` | Simulated (eta,phi) window | 

