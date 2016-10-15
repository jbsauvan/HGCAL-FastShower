#include <string>
#include <cmath>


namespace Constants {
  extern const bool debug = false;
  extern const bool fluctuation = true;

  extern const double energy = 100.;
  extern const double nhitspergev = 1000.; // nbr hits per GeV, used only if fluctuation is false 
  extern const double alpha = 0.25; // the sampling parameter, determines nhitspergev if fluctuation is true 
  extern const double mipenergy = 0.000055; // for 200 microns Si
  extern const double sampling=0.0055; // taken from 600mips in layer max for 75 GeV e-

  // electronic noise
  extern const bool noise=true;
  extern const double sigmanoise=1.; // noise in mips

  // geometry
  extern const bool readgeom=false;

  // all parameters bellow correspond to simplified parameterised geometry, ignored if readgeom=true
  extern const int nx=31; // number of cells along x
  extern const int ny=31; // number of cells along y
  //extern const double a = 0.649; // the cell side in cm for large cells from SiHexGeometry.pdf 
  extern const double a = 1.590; // multiply cell side by sqrt(6) for triangles to get equal area
  //extern const double a = 0.22727; // the cell side in cm for large cells from SiHexGeometry.pdf 
  //extern const double a = 0.476; // the cell side in cm for small cells from SiHexGeometry.pdf 

  extern const double asqrt3 = a*std::sqrt(3.);
  extern const double asqrt3over2 = asqrt3/2.;
  extern const double aover2 = a/2.;
  //extern const int itype=0; // polygon type (0=hexagons, 1=triangles)
  extern const int itype=1; // polygon type (0=hexagons, 1=triangles)

  // importing a full geometry from json file, used if readgeom=true
  extern const std::string geomfile="./AC_11.json"; // aligned centered, d=11

  // coordinates system for parameterised (infinite) geometry
  extern const int ioffsetparam=-9; // <0 starting cell along x, this is to fully fill the display window

  // layer depth
  // if between 0 and 27, the generated energy is weighted according to the layer profile (elayers[i])
  // the incident position is also shifted according to layer z position (zlayers[i]) and incident direction (etainc)
  // if set to -1, the full generated energy is deposited without longitudinal weighting and the
  // incident position is not shifted, that is it corresponds to layer 0 with the full energy deposited in it
  extern const int klayer=-1; 

  // display
  extern const double nxdisplay=15; // the number of hexagons to display along the x axis (parameterised case)
  extern const double xdisplayoffsetfull=0.5; // center of cell (0,0) x-shift (full geometry from json)
  extern const double ydisplayoffsetfull=0.5; // center of cell (0,0) y-shift (full geometry from json)
  extern const int nevtdisplay=5; // number of individual events to display

  // shooting position and direction
  // shooting position is the transverse position at HGCAL entry, ie at layer 0
  // below we shoot at center of the (4,8) cell in the parameterised hexagonal geometry
  //extern const int iinc = 4; 
  //extern const int jinc = 8; 
  //extern const double xinc = iinc*asqrt3+jinc*asqrt3over2;
  //extern const double yinc = jinc*asqrt3*asqrt3over2/a;
  // below we shoot at center of the (3,4) cell in the parameterised triangular geometry
  extern const int iinc = 3; 
  extern const int jinc = 4;              
  extern const double xinc = iinc*aover2+jinc*aover2;
  extern const double yinc = jinc*asqrt3over2+(iinc+1)%2*asqrt3/6.;
  // below we shoot at the corner od cell(3,4) in the parameterised triangular geometry
  //extern const double xinc = iinc*aover2+jinc*aover2 + aover2;
  //extern const double yinc = jinc*asqrt3over2 -asqrt3over2/3.;
  // below we shoot at the center of the module for full geometry 
  //extern const double xinc = 0.;
  //extern const double yinc = 0.;

  extern const double etainc=2.;

  // layers' z positions, this should move in the geometry
  // from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
  // the values are the silicon (centre) z positions of the 28 layers wrt HGCAL z entry position in cm 
  extern const double zlayers[28] = {.765,1.515,2.745,3.495,4.725,5.475,6.705,7.455,8.685,9.435,
    10.745,11.615,12.925,13.795,15.105,15.975,17.285,18.155,19.465,20.335,
    21.785,22.855,24.305,25.375,26.825,27.895,29.345,30.415};
 
  // shower parameters:
  extern const double radlen=0.968; // radiation length in cm
  extern const double moliere=2.3; // moliere radius in cm
  extern const double criten=5.36; // critical energy in MeV


}
