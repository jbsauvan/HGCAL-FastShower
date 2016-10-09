

#include <fstream>
#include <iostream>
#include "Geometry.h"
#include "Constants.h"
#include "json/json.h"
#include "TPolyLine.h"

using namespace Constants;

std::string fixedLength(int value, int digits = 5, std::string type="FH") {

  //  In the ex. below wants
  //  for int value 0 -> it shoud display 000
  //  for value -9 -> it should return -009
  //  for value say -50 -> it should return -050
  //  for value say -110 -> it should return -110

  if (type != "FH" && type != "VHH") {
    std::cout << "Unknown cell type name " << type << std::endl;
    return "";
  }

  unsigned int uvalue = value;
  std::string result;
  while (digits-- > 0) {
    result += ('0' + uvalue % 10);
    uvalue /= 10;
  }
  if (value < 0) {
    result += '-';
  }
  std::reverse(result.begin(), result.end());
  result = type + result;
  //std::cout << "string result " << result << std::endl;

  return result;

}


void Geometry::constructFromJson(const std::string& filename) {

  // first set klayer and zlayer according to user parameters
  setLayer(klayer);
  double zlayer;
  if (klayer == -1) zlayer = 0.;  // entry face required
  else zlayer = zlayers[klayer]; // else offset from the layer z position
  // to force being at the center for any requested layer
  //zlayer = 0.; 
  setZlayer(zlayer);

  // json format from Marina 06/2016
  // units are mm, converted in cm

  Json::Reader reader;
  Json::Value obj;
  std::ifstream ifs(filename);
  //bool success = reader.parse(ifs, obj);
  reader.parse(ifs, obj);
  const Json::Value& version = obj["Version"]; 
  std::cout << " " << std::endl;
  std::cout << "Reading geometry from JSON file: " << filename << std::endl;
  std::cout << "Geometry version " << version.asString() << std::endl;

  // module meta data
  const Json::Value& module_area = obj["Module"]["Module_area"]; 
  std::cout << "Module area : " << module_area.asDouble() << std::endl;
  const Json::Value& module_cells_and_half_cells = obj["Module"]["Module_cells_(1/1, 1/2)"]; 
  std::cout << "Module cells (1/1, 1/2) : " << module_cells_and_half_cells.asUInt() << std::endl;
  const Json::Value& module_full_cells = obj["Module"]["Module_cells_1/1"]; 
  std::cout << "Module cells_1/1 : " << module_full_cells.asUInt() << std::endl;
  const Json::Value& module_half_cells = obj["Module"]["Module_cells_1/2"]; 
  std::cout << "Module cells 1/2 : " << module_half_cells.asUInt() << std::endl;
  const Json::Value& module_third_cells = obj["Module"]["Module_cells_1/3"]; 
  std::cout << "Module cells 1/3 : " << module_third_cells.asUInt() << std::endl;
  const Json::Value& module_center_coord = obj["Module"]["Module_center_coord"]; 
  std::cout << "Module center coord : " << module_center_coord[0].asDouble() << 
    " " << module_center_coord[1].asDouble() << std::endl;
  const Json::Value& module_orientation = obj["Module"]["Module_orientation"]; 
  std::cout << "Module orientation : " << module_orientation.asDouble() << std::endl;
  //const Json::Value& module_vertex_coord = obj["Module"]["Module_vertex_coord"]; 
  std::cout << "Module vertex coord : " << std::endl;
  for (int i=0; i<6; i++) {
    const Json::Value& module_vertex_coord = obj["Module"]["Module_vertex_coord"][i]; 
    std::cout << " vertex " << i << " : " <<  module_vertex_coord[0].asDouble() << 
      " " << module_vertex_coord[1].asDouble() << std::endl;
  }   
  std::cout << " " << std::endl;

  // full hexagon cells meta data
  const Json::Value& full_hexagons_area = obj["Module"]["FH"]["FH_area"]; 
  std::cout << "Full hexagon area " << full_hexagons_area.asDouble() << std::endl;
  const Json::Value& full_hexagons_area_module = obj["Module"]["FH"]["FH_area_module"]; 
  std::cout << "Full hexagon area module " << full_hexagons_area_module.asDouble() << std::endl;
  const Json::Value& full_hexagons_count = obj["Module"]["FH"]["FH_count"]; 
  std:: cout << "Full hexagons nbr of cells " << full_hexagons_count.asUInt() << std::endl;
  std::cout << " " << std::endl;

  // construct full hexagon cells
  //for (unsigned int icell=0; icell<1; icell++) {
  for (unsigned int icell=0; icell<full_hexagons_count.asUInt(); icell++) {
    std::string cell_name = fixedLength(icell,5,"FH");
    const Json::Value& hexagon_attributes = obj["Module"]["FH"][cell_name]; 
    int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
    int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();
    std::cout << "New cell " << cell_name << " : " << std::endl;
    std::cout << " mapping coordinates : " << 
      hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
      hexagon_attributes[0]["mapping_coord"][1].asInt() << std::endl;
    std::cout << " center coordinates : " << 
      hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
      hexagon_attributes[1]["center_coord"][1].asDouble() << std::endl;
    TVectorD *position = new TVectorD(2);
    (*position)(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
    (*position)(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
    std::cout << " orientation : " << 
      hexagon_attributes[2]["orientation"].asDouble() << std::endl;
    double orientation = hexagon_attributes[2]["orientation"].asDouble(); 
    std::cout << " vertex_coordinates : " << std::endl;
    std::vector<TVectorD *> *vertices = new std::vector<TVectorD *>;
    for (int i=0; i<6; i++) {
      TVectorD *vertex = new TVectorD(2);
      std::cout << "  vertex " << i << " " << 
        hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble() << " " <<
        hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble() << std::endl;
      (*vertex)(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.; // cm
      (*vertex)(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.; // cm
      vertices->push_back(vertex);
    } 
    // create and fill cells
    Cell *aCell = new Cell(position, vertices, orientation, i_index, j_index); 
    cells_.push_back(aCell); 
  }
  std::cout << " " << std::endl;

  // half hexagon cells meta data
  const Json::Value& half_hexagons_area = obj["Module"]["Edge_VHH"]["VHH_area"]; 
  std::cout << "Half hexagon area (edges) : " << half_hexagons_area.asDouble() << std::endl;
  const Json::Value& half_hexagons_area_module = obj["Module"]["Edge_VHH"]["VHH_area_module"]; 
  std::cout << "Half hexagon area module (edges) : " << half_hexagons_area_module.asDouble() << std::endl;
  const Json::Value& half_hexagons_count = obj["Module"]["Edge_VHH"]["VHH_count"]; 
  std:: cout << "Half hexagons nbr of cells (edges) : " << half_hexagons_count.asUInt() << std::endl;

  // construct half hexagon cells
  std::cout << " " << std::endl;
  for (unsigned int icell=0; icell<half_hexagons_count.asUInt(); icell++) {
    std::string cell_name = fixedLength(icell,5,"VHH");
    const Json::Value& hexagon_attributes = obj["Module"]["Edge_VHH"][cell_name]; 
    std::cout << "New cell " << cell_name << " : " << std::endl;
    std::cout << " mapping coordinates : " << 
      hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
      hexagon_attributes[0]["mapping_coord"][1].asInt() << std::endl;
    int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
    int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();
    std::cout << " center coordinates : " << 
      hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
      hexagon_attributes[1]["center_coord"][1].asDouble() << std::endl;
    TVectorD *position = new TVectorD(2);
    (*position)(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
    (*position)(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
    std::cout << " orientation : " << 
      hexagon_attributes[2]["orientation"].asDouble() << std::endl;
    double orientation = hexagon_attributes[2]["orientation"].asDouble(); 
    std::cout << " vertex_coordinates : " << std::endl;
    std::vector<TVectorD *> *vertices = new std::vector<TVectorD *>;
    for (int i=0; i<4; i++) {
      TVectorD *vertex = new TVectorD(2);
      std::cout << "  vertex " << i << " " << 
        hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble() << " " <<
        hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble() << std::endl;
      (*vertex)(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.;// cm
      (*vertex)(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.;// cm
      vertices->push_back(vertex);
    }
    // create and fill cells
    Cell *aCell = new Cell(position, vertices, orientation, i_index, j_index); 
    cells_.push_back(aCell);     
  }  
  std::cout << " " << std::endl;

}

void Geometry::constructFromParameters(int nrows,int ncols,int klayer) {

  // a tesselation of the plane with hexagons
  std::cout << " " << std::endl;
  std::cout << "Building parametrized geometry : " << nrows << " " << ncols  << std::endl;

  setNrows(nrows);
  setNcols(ncols);
  setLayer(klayer);

  rotation60 = TMatrixD(2,2);

  std::vector<Cell *> cells;
  // reserve memory for vector of cells
  cells.reserve(nrows*ncols);

  double zlayer;
  if (klayer == -1) zlayer = 0.;  // entry face required
  else zlayer = zlayers[klayer]; // else offset from the layer z position
  // to force being at the center for any requested layer
  //zlayer = 0.; 
  setZlayer(zlayer);

  // I index run along x axis is defined such that hexagons are adjacent by side along this axis
  // J index runs along the y' axis is rotated by 60deg wrt x axis  

  // define an offset along the x-axis to fully fill the display pad with hexagons
  double xoffset = ioffsetparam*asqrt3;

  for (int i=0; i<nrows_;i++) {

    for (int j=0; j<ncols_;j++) {

      double x = xoffset + i*asqrt3 + j*asqrt3over2;
      double yprime = j*asqrt3;
      // get back to the orthogonal y axis
      double y = yprime*asqrt3over2/a;

      //std::cout << "cell i,j,x,y " << i << " " << j << " " << x << " " << y << std::endl;
      cells.push_back(new Cell(x,y));

    }

  }
  setCells(cells);

  double data[4] = {0.5, asqrt3over2/a, -asqrt3over2/a, 0.5};
  rotation60.SetMatrixArray(data);

}

void Geometry::setLayer(int klayer) {

  if (klayer<-1 || klayer>=28) {
    std::cout << "[Geometry::Geometry] error, invalid klayer " << klayer << std::endl;
    std::cout << "[Geometry::Geometry] setting klayer to entry face" << std::endl;
    klayer_ = 0;
  } else if (klayer==-1) klayer_ = 0;
  else klayer_ = klayer;

}

void Geometry::print() {

  std::cout << "Printing the geometry : " << std::endl;
  if (klayer != -1) std::cout << "the layer plane is " << klayer_ << " at z position " << zlayers[klayer_] << std::endl;
  else std::cout << "the layer plane is " << klayer_ << " at z position 0." << std::endl;
  std::vector<Cell *>::iterator ic;
  for (ic=cells_.begin();ic!=cells_.end();ic++) { 
    std::cout << "new cell with indices " <<
      "("<< (*ic)->getIIndex() << "," << (*ic)->getJIndex() << ")" <<  
      " and position " << 
      "("<<((*ic)->getPosition())(0) << "," << ((*ic)->getPosition())(1) << ")" << 
      std::endl;   
  }  

}

void Geometry::draw(double scale) {

  double summitx[7]; double summity[7];
  int nvertices=6;

  // for the case of the full geometry from json, offset the display by (0.5,0.5) such that
  // the center of the module is at the center of the pad
  double xdisplayoffset=0.;
  double ydisplayoffset=0.;
  if (readgeom) {
    xdisplayoffset=xdisplayoffsetfull;
    ydisplayoffset=ydisplayoffsetfull;
  }

  for (std::vector<Cell *>::iterator ic=cells_.begin();ic!=cells_.end();ic++) { 
    if ((*ic)->isHalfCell()) nvertices=4; // half cells   
    for (int i=0;i<nvertices;i++) summitx[i]=((*((*ic)->getVertices()[i]))(0)*scale+xdisplayoffset);
    summitx[nvertices]=((*((*ic)->getVertices()[0]))(0)*scale+xdisplayoffset);
    for (int i=0;i<nvertices;i++) summity[i]=((*((*ic)->getVertices()[i]))(1)*scale+ydisplayoffset);
    summity[nvertices]=((*((*ic)->getVertices()[0]))(1)*scale+ydisplayoffset);
    TPolyLine *hexagon = new TPolyLine(nvertices+1,summitx,summity);
    hexagon->SetFillColor(38);
    hexagon->SetLineColor(4);
    hexagon->SetLineWidth(1);
    //hexagon->Draw("f");
    hexagon->Draw();   
  } 

}

TVectorD Geometry::getPosition(int i, int j) {

  TVectorD position(2);
  if (!readgeom) { // parameterised geometry
    position(0) = ioffsetparam*asqrt3+i*asqrt3+j*asqrt3over2;
    double yprime = j*asqrt3;
    position(1) = yprime*asqrt3over2/a;
  } else { // full geometry
    std::vector<Cell *>::iterator ic;
    for (ic=cells_.begin();ic!=cells_.end();ic++) { 
      if ((*ic)->getIIndex()!=i) continue; 
      if ((*ic)->getJIndex()!=j) continue; 
    }  
    if (ic!=cells_.end()) position = (*ic)->getPosition();
  }
  return position;

}

Cell Geometry::closestCell(double x, double y) {

  // beware that this function does not explicit check that the point is within
  // the hexagon, therefore if the hexagon grid is too small it will return
  // the closest cell (as the name indicates)

  // here we implement the a temporary bruteforce function for the more general geometry read from the json file
  // this needs to be optimize for timing efficiency

  // intialize to Cel(-999,-999) to identify unfound cells
  //int ifound=-999;
  //int jfound=-999;

  double r2min=99999.;
  std::vector<Cell *>::iterator ic, icfound;
  for (ic=cells_.begin();ic!=cells_.end();ic++) { 

    double xcell = (*ic)->getPosition()(0);     
    double ycell = (*ic)->getPosition()(1);
    // now compute the distances and take the smallest one
    double r2 = (xcell-x)*(xcell-x) + (ycell-y)*(ycell-y);
    if (r2<r2min) {
      r2min=r2;
      icfound=ic;
    }

  } 

  if (icfound==cells_.end())
    std::cout << "[Geometry::closestCell] Cell not found!! x, y " << 
      x << " " << y << std::endl;

  return Cell(**icfound);

}

TVectorD Geometry::positionInCell(TVectorD position) {

  TVectorD relativeposition(2);
  relativeposition=position-closestCell(position(0),position(1)).getPosition();
  return relativeposition;

}

bool Geometry::isInCell(TVectorD position, Cell cell) {

  // implementation below works for any convex cell described by its vertices
  // assumes vertices are consecutive along the cell perimeter and ordered along direct rotation

  // loop on pair of consective vertices 
  for (unsigned int i=0;i<cell.getVertices().size()-1;i++) {

    //TVectorD *vertex = new TVectorD(2);
    double xa =(*cell.getVertices()[i])(0);
    double ya =(*cell.getVertices()[i])(1);
    double xb =(*cell.getVertices()[i+1])(0);
    double yb =(*cell.getVertices()[i+1])(1);
    double sign = (xb-xa)*(position(1)-ya) - (yb-ya)*(position(0)-xa);
    if (sign<0.) return false;
  }  

  return true;

}

int Geometry::getIIndex(Cell cell) {

  return cell.getIIndex();

}

int Geometry::getJIndex(Cell cell) {

  return cell.getJIndex();

}

