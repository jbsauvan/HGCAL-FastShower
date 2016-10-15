
#ifndef ShowerShape_h__
#define ShowerShape_h__

#include <map>
#include "Cell.h"

class ShowerShape {

  // a class for groups of cells
  
  public:
  
    ShowerShape() {} // default constructor
    ShowerShape(std::map<Cell*,double,CellComp>* enrjMap);
    ~ShowerShape() {if (enrjMap_) delete enrjMap_;}

    Cell *maxCell() {return maxCell_;}
    double maxE1() {return maxE1_;} // hotest cell energy
    int imax() {return imax_;}
    int jmax() {return jmax_;}
    double firstNeighboors();
    virtual double firstNeighboors(int i, int j)=0;
    double energySum() {return energySum_;} // total energy

  protected:
  
    std::map<Cell*,double,CellComp>* enrjMap_;  
    std::map<int,Cell*> cellMap_;  // cell 1D indexing
    Cell *maxCell_; // hotest cell
    double maxE1_; // hotest cell energy
    int imax_; // hotest cell i index
    int jmax_; // hotest cell j index
    double energySum_;
  
};

#endif
