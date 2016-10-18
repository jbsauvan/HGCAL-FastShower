
#ifndef __HGCalSimulation_FastShower_ShowerShapeHexagon_h__
#define __HGCalSimulation_FastShower_ShowerShapeHexagon_h__

#ifdef STANDALONE
#include "ShowerShape.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShape.h"
#endif

class ShowerShapeHexagon: public ShowerShape {

  // a class for groups of cells
  
  public:
  
    ShowerShapeHexagon() {} // default constructor
    ShowerShapeHexagon(std::map<Cell*,double,CellComp>* enrjMap): ShowerShape(enrjMap) {}
    ~ShowerShapeHexagon() {}

    double firstNeighboors(int i, int j);

  private:
  
};

#endif
