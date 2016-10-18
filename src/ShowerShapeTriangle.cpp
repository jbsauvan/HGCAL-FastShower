
#include "TMath.h"

#ifdef STANDALONE
#include "ShowerShapeTriangle.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShapeTriangle.h"
#endif



double ShowerShapeTriangle::firstNeighboors(int i, int j) {

  // todo: add protections for the grid boarders
  int index;
  double sum=0.;
  Cell *cell;
  if (imax_%2==0) { // upward triangles
    index = 1000*(i+1)+j;
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
    index = 1000*(i-1)+j;
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
    index = 1000*(i+1)+j-1;
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
  } else { // downward triangles
    index = 1000*(i-1)+j;
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
    index = 1000*(i-1)+(j+1);
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
    index = 1000*(i+1)+j;
    cell = cellMap_[index];
    sum += (*enrjMap_)[cell];
  }
  return sum;

}
