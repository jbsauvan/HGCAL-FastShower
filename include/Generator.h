#ifndef Generator_h__
#define Generator_h__

#include <map>
#include "Cell.h"
#include "Geometry.h"
#include "TH1F.h"


class Generator {

  public:
    Generator() {};
    ~Generator() {};

    void simulate(int);
    //void display(Geometry& geometry, std::map<Cell,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt=0);
    void display(Geometry& geometry, std::map<Cell*,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt=0);

};


#endif
