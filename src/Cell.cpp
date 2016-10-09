#include "Cell.h"
#include "Constants.h"
#include "TMath.h"


using namespace Constants;

Cell::Cell(const Cell& cell) {
  orientation_ = cell.getOrientation();
  i_index_ = cell.getIIndex();
  j_index_ = cell.getJIndex();
  position_ = new TVectorD(2);
  (*position_)(0)=(cell.getPosition())(0); 
  (*position_)(1)=(cell.getPosition())(1);
  vertices_ = new std::vector<TVectorD *>;
  for (unsigned int i=0;i<cell.getVertices().size();i++) {
    TVectorD *vertex = new TVectorD(2);
    (*vertex)(0)=(*cell.getVertices()[i])(0);
    (*vertex)(1)=(*cell.getVertices()[i])(1);
    vertices_->push_back(vertex);
  }

}

Cell& Cell::operator=(const Cell& cell) { 
  if (this != &cell) {
    delete position_;
    position_ = new TVectorD(2);
    (*position_)(0)=(cell.getPosition())(0); 
    (*position_)(1)=(cell.getPosition())(1); 
    vertices_ = new std::vector<TVectorD *>;
    for (unsigned int i=0;i<cell.getVertices().size();i++) {
      TVectorD *vertex = new TVectorD(2);
      (*vertex)(0)=(*cell.getVertices()[i])(0);
      (*vertex)(1)=(*cell.getVertices()[i])(1);
      vertices_->push_back(vertex);
    }  
    for (unsigned int i=0;i<cell.getVertices().size();i++) {
      delete cell.getVertices()[i];
    }
    orientation_ = cell.getOrientation();
    i_index_ = cell.getIIndex();
    j_index_ = cell.getJIndex();  
  }
  return *this;

}

Cell::Cell(TVectorD *position, std::vector<TVectorD *> *vertices, double orientation, int i_index, int j_index) {
  position_ = position;
  vertices_ = vertices;
  orientation_ = orientation;
  i_index_ = i_index;
  j_index_ = j_index;

}

Cell::Cell(double x, double y){
  position_ = new TVectorD(2);
  (*position_)(0) = x; 
  (*position_)(1) = y;
  orientation_ = 90.;
  double t = y*a/asqrt3;
  t = x - t;
  i_index_ = TMath::Nint(t/asqrt3);
  double yprime = y*a/asqrt3over2;
  j_index_ = TMath::Nint(yprime/asqrt3);
  std::vector<TVectorD *> *vertices = new std::vector<TVectorD *>;
  for (int ii=0;ii<6;ii++) {
    TVectorD *vertex = new TVectorD(2);
    (*vertex)(0) = x+offsetx[ii];
    (*vertex)(1) = y+offsety[ii];
    vertices->push_back(vertex);
  }  
  vertices_ = vertices;

}

Cell::Cell(int i, int j) {
  position_ = new TVectorD(2);
  (*position_)(0)=ioffsetparam*asqrt3+i*asqrt3+j*asqrt3over2;
  double yprime=j*asqrt3;
  (*position_)(1)=yprime*asqrt3over2/a;
  orientation_ = 90.;
  i_index_ = i;
  j_index_ = j;
  std::vector<TVectorD *> *vertices = new std::vector<TVectorD *>;
  for (int ii=0;ii<6;ii++) {
    TVectorD *vertex = new TVectorD(2);
    (*vertex)(0) = (*position_)(0)+offsetx[ii];
    (*vertex)(1) = (*position_)(1)+offsety[ii];
    vertices->push_back(vertex);
  }  
  vertices_ = vertices;
}

Cell::~Cell() {
  if (position_) delete position_;
  if (vertices_) {
    for (unsigned int i=0;i<vertices_->size();i++) {
      delete getVertices()[i];
    }
    delete vertices_;
  }

}
