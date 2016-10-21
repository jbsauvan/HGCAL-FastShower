
#ifndef __HGCalSimulation_FastShower_Geometry_h__
#define __HGCalSimulation_FastShower_Geometry_h__

#include <string>
#include <vector>
#include <unordered_map>
#include "TVectorD.h"
#include "TMatrixD.h"
#ifdef STANDALONE
#include "Cell.h"
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif


class Geometry {

  // a tesselation of the plane with polygonal cells
  // by default the plane is at z=0 but can be shifted to a z position that
  // corresponds to a layer position from TP geometry 
  // itype=0: hexagonal cells
  // itype=1: triangular cells

  public:

    Geometry(const Parameters::Geometry& params):
      parameters_(params){}
    ~Geometry() {}

    // FIXME: can simplify the interface by having one single public construct method
    // construct parametrised geometry, default grid 11x11, ie +-5 around the central cell
    void constructFromParameters(bool); 
    // construct geometry from json file
    void constructFromJson(bool);

    const Cell& closestCell(double x, double y) const; // the cell that contains the point
    bool isInCell(const TVectorD& position, const Cell& cell) const; // test if a point is within a cell
    //bool isInRealCell(TVectorD position, Cell cell); // to apply further mouse bite or virtual dead region 
    TVectorD positionInCell(const TVectorD& position) const; // relative position within the cell

    const TVectorD& getPosition(int i, int j) const; // position of cell i,j

    //vect<Cell> getNeighbours (int i, radius r); // not yet implemented
    //vect<Cell> getFirstNeighbours (int i); // not yet implemented
    //vect<Cell> getSecondNeighbours (int i); // not yet implemented

    // getters
    int getNumberOfRows() const {return nrows_;} // parameterized case
    int getNumberOfCols() const {return ncols_;} // parameterized case
    const std::unordered_map<uint32_t, Cell>& getCells() const {return cells_;}
    //int getIIndex(const Cell& cell) const;
    //int getJIndex(const Cell& cell) const;
    int getLayer() const {return klayer_;}
    double getZlayer() const {return zlayer_;}
    Parameters::Geometry::Type getType() const {return itype_;} 

    double a() const {return a_;}
    double asqrt3() const {return asqrt3_;}
    double aover2() const {return aover2_;}
    double asqrt3over2() const {return asqrt3over2_;}

    void draw(const Parameters::Display& params, double scale=0.1);
    void print();

  private:

    //void setCells(std::vector<Cell *> cells) {cells_ = cells;}
    void setNrows(int nrows) {nrows_ = nrows;}
    void setNcols(int ncols) {ncols_ = ncols;}
    void setLayer(int klayer);
    void setZlayer(double zlayer) {zlayer_ = zlayer;}
    void setType (Parameters::Geometry::Type itype) {itype_=itype;}

    std::unordered_map<uint32_t, Cell> cells_;
    int nrows_;
    int ncols_;
    int klayer_;
    double zlayer_;
    Parameters::Geometry::Type itype_; // cell type
    double a_;
    double asqrt3_;
    double aover2_;
    double asqrt3over2_;
    const Parameters::Geometry& parameters_;

};



#endif
