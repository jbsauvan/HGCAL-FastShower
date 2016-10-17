#ifndef Generator_h__
#define Generator_h__

#include <map>
#include "Cell.h"
#include "Geometry.h"
#include "Parameters.h"
#include "TH1F.h"


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
