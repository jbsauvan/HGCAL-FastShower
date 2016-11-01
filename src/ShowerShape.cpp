
#include <iostream>
#include "TMath.h"

#ifdef STANDALONE
#include "ShowerShape.h"
#else
#include "HGCalSimulation/FastShower/interface/ShowerShape.h"
#endif


ShowerShape::ShowerShape(const std::unordered_map<uint32_t,double>& enrjMap,
    const std::unordered_map<uint32_t,Cell>& cellMap):
  enrjMap_(enrjMap),
  cellMap_(cellMap) {
  
  energySum_ = 0.;
  maxE1_ = std::numeric_limits<double>::lowest();
  imax_ = std::numeric_limits<int16_t>::lowest(); 
  jmax_ = std::numeric_limits<int16_t>::lowest();
  for (const auto& id_cell : cellMap) { 
    // FIXME: should avoid searching into enrjMap even if unordered_map search is constant time
    double energy = 0.;
    try {
      energy = enrjMap.at(id_cell.first);
    } catch(const std::out_of_range&) {
      // do nothing (energy = 0)
    }

    energySum_ = energySum_ + energy;
    if (energy>maxE1_) {
      maxCell_ = &id_cell.second;
      maxE1_ = energy;
      imax_ = maxCell_->getIIndex();
      jmax_ = maxCell_->getJIndex();
    }
  } 
  //std::cout << "ShowerShape::ShowerShape imax, jmax " << imax_ << " " << jmax_ << std::endl;   
  
  // now prepare structure for easy navigation
  //int index=-999; //assume grid cannot be larger than 1000 in j
  //for (ic=enrjMap_->begin();ic!=enrjMap_->end();ic++) { 
    //if (ic->first->getJIndex()>999) {
      //std::cout << "[ShowerShape::ShowerShape] error, j index greater than 1000!! " << std::endl;
      //return;
    //} else {
      //index = 1000*ic->first->getIIndex() + ic->first->getJIndex();
      //cellMap_[index] = ic->first;   
    //}
  //}
  
}

double ShowerShape::firstNeighboors() {

  return firstNeighboors(imax_,jmax_);

}
