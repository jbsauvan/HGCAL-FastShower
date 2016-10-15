
#ifndef ShowerShapeHexagon_h__
#define ShowerShapeHexagon_h__

#include "ShowerShape.h"

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
