

#ifndef Cell_h__
#define Cell_h__

#include "TVectorD.h"

class Cell {
  // an hexagonal cell

  public:
    Cell() {} // default constructor
    // constructor for full geometries
    Cell(TVectorD *, std::vector<TVectorD *> *, double, int, int); 
    Cell(const Cell&);
    Cell& operator=(const Cell&);
    ~Cell();

    // getters
    TVectorD getPosition() const {return *position_;}
    std::vector<TVectorD *> getVertices() const {return *vertices_;}
    double getOrientation() const {return orientation_;}
    int getIIndex() const {return i_index_;}
    int getJIndex() const {return j_index_;}
    bool isFullCell() const {return orientation_==90.;}
    bool isHalfCell() const {return int(orientation_)%60==0;}

  private:
    TVectorD *position_; // centre position in absolute coordinates
    std::vector<TVectorD *> *vertices_; // vertices positions in absolute coordinates
    double orientation_; // orientation for halh cells
    int i_index_; // I index
    int j_index_;  // J index

};

struct CellComp {
  bool operator() (const Cell& c1, const Cell& c2) const {
    if (c1.getIIndex()!=c2.getIIndex()) return (c1.getIIndex()<c2.getIIndex());
    else return (c1.getJIndex()<c2.getJIndex());
  }
};


#endif
