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

##Running
After the installation steps an executable is produced, `shower_simulation.exe`. It takes one argument, which is a python file containing parameter definitions.  
 A test file is available in the `python` directory.
```bash
shower_simulation.exe python/test_cfg.py
```
