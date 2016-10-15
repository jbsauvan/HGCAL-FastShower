
#ifndef ShowerShapeTriangle_h__
#define ShowerShapeTriangle_h__

#include "ShowerShape.h"

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
