

#include "TMath.h"

#ifdef STANDALONE
#include "ShowerShapeHexagon.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShapeHexagon.h"
#endif



double ShowerShapeHexagon::firstNeighboors(int i, int j) {

  // todo: add protections for the grid boarders
  int index;
  double sum=0.;
  Cell *cell;
  index = 1000*(i+1)+j;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];
  index = 1000*(i)+j+1;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];  
  index = 1000*(i-1)+j+1;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];
  index = 1000*(i-1)+j;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];
  index = 1000*(i)+j-1;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];
  index = 1000*(i+1)+j-1;
  cell = cellMap_[index];
  sum += (*enrjMap_)[cell];
  return sum;
  
}
