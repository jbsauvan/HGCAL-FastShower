
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
  if (imax_%2==0) { // upward triangles
    sum += enrjMap_.at(Cell::id(i+1,j));
    sum += enrjMap_.at(Cell::id(i-1,j));
    sum += enrjMap_.at(Cell::id(i+1,j-1));
  } else { // downward triangles
    sum += enrjMap_.at(Cell::id(i-1,j));
    sum += enrjMap_.at(Cell::id(i-1,j+1));
    sum += enrjMap_.at(Cell::id(i+1,j));
  }
  return sum;

}
