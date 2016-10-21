
#ifndef __HGCalSimulation_FastShower_ShowerShape_h__
#define __HGCalSimulation_FastShower_ShowerShape_h__

#include <unordered_map>
#ifdef STANDALONE
#include "Cell.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#endif

class ShowerShape {

  // a class for groups of cells
  
  public:
  
    ShowerShape(const std::unordered_map<uint32_t,double>&,
        const std::unordered_map<uint32_t,Cell>&);
    ~ShowerShape() {}

    const Cell* maxCell() {return maxCell_;}
    double maxE1() {return maxE1_;} // hotest cell energy
    int imax() {return imax_;}
    int jmax() {return jmax_;}
    double firstNeighboors();
    virtual double firstNeighboors(int i, int j)=0;
    double energySum() {return energySum_;} // total energy

  protected:
  
    const std::unordered_map<uint32_t,double>& enrjMap_;  
    const std::unordered_map<uint32_t,Cell>& cellMap_;  // cell 1D indexing
    const Cell* maxCell_; // hotest cell
    double maxE1_; // hotest cell energy
    int imax_; // hotest cell i index
    int jmax_; // hotest cell j index
    double energySum_;
  
};

#endif
