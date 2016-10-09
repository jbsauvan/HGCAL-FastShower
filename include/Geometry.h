
#ifndef Geometry_h__
#define Geometry_h__

#include <string>
#include <vector>
#include "TVectorD.h"
#include "TMatrixD.h"
#include "Cell.h"


class Geometry {

  // a tesselation of the plane with hexagonal cells
  // by default the plane is at z=0 but can be shifted to a z position that
  // corresponds to a layer position from TP geometry 

  public:

    Geometry() {}
    ~Geometry() {}

    // construct parametrised geometry
    void constructFromParameters(int nrows=11,int ncols=11,int klayer=-1); //default grid 11x11, ie +-5 around the central cell
    // construct geometry from json
    void constructFromJson(const std::string&);
    //std::vector<Cell *> getCells() {return cells_;}

    Cell closestCell(double x, double y); // the cell that contains the point
    bool isInCell(TVectorD position, Cell cell); // test if a point is within a cell
    //bool isInRealCell(TVectorD position, Cell cell); // to apply further mouse bite or virtual dead region 
    TVectorD positionInCell(TVectorD position); // relative position within the cell

    TVectorD getPosition(int i, int j); // position of cell i,j

    //vect<Cell> getNeighbours (int i, radius r); // not yet implemented
    //vect<Cell> getFirstNeighbours (int i); // not yet implemented
    //vect<Cell> getSecondNeighbours (int i); // not yet implemented

    // getters
    int getNumberOfRows() {return nrows_;} // parameterized case
    int getNumberOfCols() {return ncols_;} // parameterized case
    std::vector<Cell *> getCells() {return cells_;}
    int getIIndex(Cell cell);
    int getJIndex(Cell cell);
    int getLayer() {return klayer_;}
    double getZlayer() {return zlayer_;}

    void draw(double scale=0.1);
    void print();

  private:

    void setCells(std::vector<Cell *> cells) {cells_ = cells;}
    void setNrows(int nrows) {nrows_ = nrows;}
    void setNcols(int ncols) {ncols_ = ncols;}
    void setLayer(int klayer);
    void setZlayer(double zlayer) {zlayer_ = zlayer;}

    std::vector<Cell *> cells_;
    int nrows_;
    int ncols_;
    int klayer_;
    double zlayer_;

    TMatrixD rotation60;

};



#endif
