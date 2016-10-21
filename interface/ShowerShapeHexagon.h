
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
  
    ShowerShapeHexagon(const std::unordered_map<uint32_t,double>& enrjMap,
        const std::unordered_map<uint32_t,Cell>& cellMap): ShowerShape(enrjMap, cellMap) {}
    ~ShowerShapeHexagon() {}

    double firstNeighboors(int i, int j);

  private:
  
};

#endif
