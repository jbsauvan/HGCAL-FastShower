
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
  
  ShowerShapeTriangle(const std::unordered_map<uint32_t,double>& enrjMap,
        const std::unordered_map<uint32_t,Cell>& cellMap): ShowerShape(enrjMap, cellMap) {}
  ~ShowerShapeTriangle() {}
  
  double firstNeighboors(int i, int j);
  
  private:
  
};

#endif
