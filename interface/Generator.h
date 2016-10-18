#ifndef __HGCalSimulation_FastShower_Generator_h__
#define __HGCalSimulation_FastShower_Generator_h__

#include <map>
#include "TH1F.h"

#ifdef STANDALONE
#include "Cell.h"
#include "Geometry.h"
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif


class Generator {

  public:
    Generator(const Parameters&);
    ~Generator() {};

    void simulate();
    //void display(Geometry& geometry, std::map<Cell,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt=0);
    void display(std::map<Cell*,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt=0);


  private:
    Geometry geometry_;
    const Parameters& parameters_;
    

};


#endif
