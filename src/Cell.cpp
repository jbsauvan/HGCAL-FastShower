

#include "TMath.h"
#ifdef STANDALONE
#include "Cell.h"
#else
#include "HGCalSimulation/FastShower/interface/Cell.h"
#endif


uint32_t Cell::id(int16_t i, int16_t j)
{
  uint32_t id = 0;
  id |= i;
  id |= (j<<16);
  return id;
}


//Cell::Cell(const Cell& cell):
  //position_(cell.getPosition()),
  //vertices_(cell.getVertices()),
  //orientation_(cell.getOrientation()),
  //i_index_(cell.getIIndex()),
  //j_index_(cell.getJIndex()),
//{
//}

//Cell& Cell::operator=(const Cell& cell) { 
  //if (this != &cell) {
    //delete position_;
    //position_ = new TVectorD(2);
    //(*position_)(0)=(cell.getPosition())(0); 
    //(*position_)(1)=(cell.getPosition())(1); 
    //vertices_ = new std::vector<TVectorD *>;
    //for (unsigned int i=0;i<cell.getVertices().size();i++) {
      //TVectorD *vertex = new TVectorD(2);
      //(*vertex)(0)=(*cell.getVertices()[i])(0);
      //(*vertex)(1)=(*cell.getVertices()[i])(1);
      //vertices_->push_back(vertex);
    //}  
    //for (unsigned int i=0;i<cell.getVertices().size();i++) {
      //delete cell.getVertices()[i];
    //}
    //orientation_ = cell.getOrientation();
    //i_index_ = cell.getIIndex();
    //j_index_ = cell.getJIndex();  
  //}
  //return *this;

//}

Cell::Cell(TVectorD&& position, std::vector<TVectorD>&& vertices, double orientation, int i_index, int j_index):
  position_(std::move(position)),
  vertices_(std::move(vertices)),
  orientation_(orientation)
{
  if(i_index<std::numeric_limits<int16_t>::min() ||
      i_index>std::numeric_limits<int16_t>::max() ||
      j_index<std::numeric_limits<int16_t>::min() ||
      j_index>std::numeric_limits<int16_t>::max()
    )
  {
    throw std::string("Cell index outside of 16bits integer limits");
  }
  i_index_ = (int16_t)i_index;
  j_index_ = (int16_t)j_index;
  id_ = Cell::id(i_index_, j_index_);
}

