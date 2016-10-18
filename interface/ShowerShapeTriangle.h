
#ifndef __HGCalSimulation_FastShower_ShowerShapeTriangle_h__
#define __HGCalSimulation_FastShower_ShowerShapeTriangle_h__

#ifdef STANDALONE
#include "ShowerShape.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShape.h"
#endif

class ShowerShapeTriangle: public ShowerShape {

  // a class for groups of cells
  
  public:
  
  ShowerShapeTriangle() {} // default constructor
  ShowerShapeTriangle(std::map<Cell*,double,CellComp>* enrjMap): ShowerShape(enrjMap) {}
  ~ShowerShapeTriangle() {}
  
  double firstNeighboors(int i, int j);
  
  private:
  
};

#endif
