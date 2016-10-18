#include "ShowerShape.h"
#include "TMath.h"
#include <iostream>


ShowerShape::ShowerShape(std::map<Cell*,double,CellComp>* enrjMap):
 enrjMap_(enrjMap) {
  
  energySum_ = 0.;
  maxE1_=-999.;
  imax_=-999; 
  jmax_=-999;
  std::map<Cell*,double,CellComp>::iterator ic;
  for (ic=enrjMap_->begin();ic!=enrjMap_->end();ic++) { 
    energySum_ = energySum_ + ic->second;
    if (ic->second>maxE1_) {
      maxCell_ = ic->first;
      maxE1_ = ic->second;
      imax_ = maxCell_->getIIndex();
      jmax_ = maxCell_->getJIndex();
    }
  } 
  //std::cout << "ShowerShape::ShowerShape imax, jmax " << imax_ << " " << jmax_ << std::endl;   
  
  // now prepare structure for easy navigation
  int index=-999; //assume grid cannot be larger than 1000 in j
  for (ic=enrjMap_->begin();ic!=enrjMap_->end();ic++) { 
    if (ic->first->getJIndex()>999) {
      std::cout << "[ShowerShape::ShowerShape] error, j index greater than 1000!! " << std::endl;
      return;
    } else {
      index = 1000*ic->first->getIIndex() + ic->first->getJIndex();
      cellMap_[index] = ic->first;   
    }
  }
  
}

double ShowerShape::firstNeighboors() {

  return firstNeighboors(imax_,jmax_);

}
