#ifndef Constants_h__
#define Constants_h__

#include <string>

namespace Constants {
  extern const bool debug;
  extern const bool fluctuation;

  extern const double energy;
  extern const double nhitspergev; // nbr hits per GeV, used only if fluctuation is false 
  extern const double alpha; // the sampling parameter, determines nhitspergev if fluctuation is true 
  extern const double mipenergy; // mip energy for 200 microns Si

  // geometry
  extern const bool readgeom;

  // all parameters bellow correspond to simplified parameterised geometry, ignored if readgeom=true
  extern const int nx; // number of cells along x
  extern const int ny; // number of cells along y
  extern const double a; // the hexagon side in cm for large cells from SiHexGeometry.pdf 
  extern const double asqrt3;
  extern const double asqrt3over2;
  extern const double aover2;
  extern const int itype;

  // importing a full geometry from json file, used if readgeom=true
  extern const std::string geomfile; // aligned centered, d=11

  // coordinates system for parameterised (infinite) geometry
  extern const int ioffsetparam; // <0 starting cell along x, this is to fully fill the display window

  // layer depth
  // if between 0 and 27, the generated energy is weighted according to the layer profile (elayers[i])
  // the incident position is also shifted according to layer z position (zlayers[i]) and incident direction (etainc)
  // if set to -1, the full generated energy is deposited without longitudinal weighting and the
  // incident position is not shifted, that is it corresponds to layer 0 with the full energy deposited in it
  extern const int klayer; 

  // display
  extern const double nxdisplay; // the number of hexagons to display along the x axis (parameterised case)
  extern const double xdisplayoffsetfull; // center of cell (0,0) x-shift (full geometry from json)
  extern const double ydisplayoffsetfull; // center of cell (0,0) y-shift (full geometry from json)
  extern const int nevtdisplay; // number of individual events to display

  // shooting position and direction
  // shooting position is the transverse position at HGCAL entry, ie at layer 0
  // below we shoot at center of the (4,8) cell in the parameterised geometry
  extern const int iinc; 
  extern const int jinc;              
  extern const double xinc;
  extern const double yinc;

  extern const double etainc;

  // layers' z positions
  // from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
  // the values are the silicon (centre) z positions of the 28 layers wrt HGCAL z entry position in cm 
  extern const double zlayers[28];
  // layers energy profile
  // from TP geometry, average for 35 geV Pt electrons
  // to simulate 28-layers V7 geometry, group the first two layers and drop the last one according to:
  // https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
  extern const double elayers[28];

  // shower parameters:
  // transverse profile described by an exponential at a given depth
  // exponential parameter set from TP studies, 90% containment in 2.3cm at layer 15
  extern const double r0layer15; 
  // evolution vs depth described by a parabolic function
  // set from TP studies (AMM), accuracy better than 5%
  extern const double a0; 
  extern const double a1;
  extern const double a2;


}

#endif
