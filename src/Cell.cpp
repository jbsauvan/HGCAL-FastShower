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

Cell::~Cell() {
  if (position_) delete position_;
  if (vertices_) {
    for (unsigned int i=0;i<vertices_->size();i++) {
      delete getVertices()[i];
    }
    delete vertices_;
  }

}
