

#include <fstream>
#include <iostream>
#include <algorithm>

#include "TPolyLine.h"

#ifdef STANDALONE
#include "Geometry.h"
#include "json/json.h"
#else
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#include "HGCalSimulation/FastShower/interface/json/json.h"
#endif

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


void Geometry::constructFromJson(bool debug) {

  itype_ = Parameters::Geometry::Type::External;
  const std::string& filename = parameters_.file;

  // first set klayer and zlayer according to user parameters
  setLayer(parameters_.layer);
  double zlayer;
  if (klayer_ == -1) zlayer = 0.;  // entry face required
  else zlayer = parameters_.layers_z[klayer_]; // else offset from the layer z position
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
  if(debug) {
    std::cout << " " << std::endl;
    std::cout << "Reading geometry from JSON file: " << filename << std::endl;
    std::cout << "Geometry version " << version.asString() << std::endl;
  }

  // module meta data
  const Json::Value& module = obj["Module"];
  if(module.isNull()) throw std::string("No module information found in json file");
  const Json::Value& module_area = module["Module_area"]; 
  const Json::Value& module_cells_and_half_cells = module["Module_cells_(1/1, 1/2)"]; 
  const Json::Value& module_full_cells = module["Module_cells_1/1"]; 
  const Json::Value& module_half_cells = module["Module_cells_1/2"]; 
  const Json::Value& module_third_cells = module["Module_cells_1/3"]; 
  const Json::Value& module_center_coord = module["Module_center_coord"]; 
  const Json::Value& module_orientation = module["Module_orientation"]; 
  const Json::Value& module_vertices = module["Module_vertex_coord"];
  //bool module_vertices_ok = false;
  std::vector<std::pair<double,double>> module_vertex_coordinates;
  if(!module_vertices.isNull() && module_vertices.isArray()) {
    //module_vertices_ok = true;
    for (unsigned i=0; i<module_vertices.size(); i++) {
      const Json::Value& coord = module_vertices[i]; 
      if(coord.isNull() || !coord.isArray() || coord.size()!=2) {
        //module_vertices_ok = false;
        break;
      }
      module_vertex_coordinates.emplace_back(coord[0].asDouble(), coord[1].asDouble());
    }
  }
  if(!module_cells_and_half_cells.isNull() &&
      !module_full_cells.isNull() && !module_half_cells.isNull()) {
    if(module_cells_and_half_cells.asUInt()!=module_half_cells.asUInt()+module_full_cells.asUInt()) {
      std::cout<<"Inconsistency in the number of full and half cells\n";
    }
  }
  // FIXME: should add more checks
  if(debug) {
    std::cout << "Module area : " << module_area.asDouble() << std::endl;
    std::cout << "Module cells (1/1, 1/2) : " << module_cells_and_half_cells.asUInt() << std::endl;
    std::cout << "Module cells_1/1 : " << module_full_cells.asUInt() << std::endl;
    std::cout << "Module cells 1/2 : " << module_half_cells.asUInt() << std::endl;
    std::cout << "Module cells 1/3 : " << module_third_cells.asUInt() << std::endl;
    std::cout << "Module center coord : " << module_center_coord[0].asDouble() << 
      " " << module_center_coord[1].asDouble() << std::endl;
    std::cout << "Module orientation : " << module_orientation.asDouble() << std::endl;
    std::cout << "Module vertex coord : " << std::endl;
    std::cout << " " << std::endl;
  }

  // full hexagon cells meta data
  const Json::Value& hexagons = module["FH"];
  if(hexagons.isNull() || !hexagons.isObject()) throw std::string("Cannot find full hexagons information");
  const Json::Value& full_hexagons_area = hexagons["FH_area"]; 
  const Json::Value& full_hexagons_area_module = hexagons["FH_area_module"]; 
  const Json::Value& full_hexagons_count = hexagons["FH_count"]; 
  // FIXME: should add more checks
  if(debug) {
    std::cout << "Full hexagon area " << full_hexagons_area.asDouble() << std::endl;
    std::cout << "Full hexagon area module " << full_hexagons_area_module.asDouble() << std::endl;
    std:: cout << "Full hexagons nbr of cells " << full_hexagons_count.asUInt() << std::endl;
    std::cout << " " << std::endl;
  }

  // construct full hexagon cells
  for (unsigned int icell=0; icell<full_hexagons_count.asUInt(); icell++) {
    std::string cell_name = fixedLength(icell,5,"FH");
    // FIXME: add more checks
    const Json::Value& hexagon_attributes = hexagons[cell_name]; 
    int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
    int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();
    TVectorD position(2);
    position(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
    position(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
    double orientation = hexagon_attributes[2]["orientation"].asDouble();
    std::vector<TVectorD> vertices;
    for (unsigned i=0; i<hexagon_attributes[3]["vertex_coord_abs"].size(); i++) {
      vertices.emplace_back(2);
      vertices.back()(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.; // cm
      vertices.back()(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.; // cm

    }
    cells_.emplace(
        Cell::id(i_index, j_index),
        Cell(std::move(position), std::move(vertices), orientation, i_index, j_index)
        );

    if(debug) {
      std::cout << "New cell " << cell_name << " : " << std::endl;
      std::cout << " mapping coordinates : " << 
        hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
        hexagon_attributes[0]["mapping_coord"][1].asInt() << std::endl;
      std::cout << " center coordinates : " << 
        hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
        hexagon_attributes[1]["center_coord"][1].asDouble() << std::endl;
      std::cout << " orientation : " << 
        hexagon_attributes[2]["orientation"].asDouble() << std::endl;
      std::cout << " vertex_coordinates : " << std::endl;
      for(const auto& vertex : vertices) {
        std::cout << "  vertex " << 
          vertex(0) << " " <<
          vertex(1) << std::endl;
      }
    }
  }
  if(debug) std::cout << " " << std::endl;

  // half hexagon cells meta data
  const Json::Value& half_hexagons = module["Edge_VHH"];
  const Json::Value& half_hexagons_area = half_hexagons["VHH_area"]; 
  const Json::Value& half_hexagons_area_module = half_hexagons["VHH_area_module"]; 
  const Json::Value& half_hexagons_count = half_hexagons["VHH_count"]; 
  // FIXME: should add more checks
  if(debug) {
    std::cout << "Half hexagon area (edges) : " << half_hexagons_area.asDouble() << std::endl;
    std::cout << "Half hexagon area module (edges) : " << half_hexagons_area_module.asDouble() << std::endl;
    std:: cout << "Half hexagons nbr of cells (edges) : " << half_hexagons_count.asUInt() << std::endl;
    std::cout << " " << std::endl;
  }

  // construct half hexagon cells
  for (unsigned int icell=0; icell<half_hexagons_count.asUInt(); icell++) {
    std::string cell_name = fixedLength(icell,5,"VHH");
    // FIXME: add more checks
    const Json::Value& hexagon_attributes = half_hexagons[cell_name]; 
    int i_index = hexagon_attributes[0]["mapping_coord"][0].asInt();
    int j_index = hexagon_attributes[0]["mapping_coord"][1].asInt();
    TVectorD position(2);
    position(0) = hexagon_attributes[1]["center_coord"][0].asDouble()/10.; // cm
    position(1) = hexagon_attributes[1]["center_coord"][1].asDouble()/10.; // cm
    double orientation = hexagon_attributes[2]["orientation"].asDouble(); 
    std::vector<TVectorD> vertices;
    for (unsigned i=0; i<hexagon_attributes[3]["vertex_coord_abs"].size(); i++) {
      vertices.emplace_back(2);
      vertices.back()(0) = hexagon_attributes[3]["vertex_coord_abs"][i][0].asDouble()/10.;// cm
      vertices.back()(1) = hexagon_attributes[3]["vertex_coord_abs"][i][1].asDouble()/10.;// cm
    }
    cells_.emplace(
        Cell::id(i_index, j_index),
        Cell(std::move(position), std::move(vertices), orientation, i_index, j_index)
        ); 

    if(debug) {
      std::cout << "New cell " << cell_name << " : " << std::endl;
      std::cout << " mapping coordinates : " << 
        hexagon_attributes[0]["mapping_coord"][0].asInt() << " " <<
        hexagon_attributes[0]["mapping_coord"][1].asInt() << std::endl;
      std::cout << " center coordinates : " << 
        hexagon_attributes[1]["center_coord"][0].asDouble() << " " <<
        hexagon_attributes[1]["center_coord"][1].asDouble() << std::endl;

      std::cout << " orientation : " << 
        hexagon_attributes[2]["orientation"].asDouble() << std::endl;
      std::cout << " vertex_coordinates : " << std::endl;
      for(const auto& vertex : vertices) {
        std::cout << "  vertex " << 
          vertex(0) << " " <<
          vertex(1) << std::endl;
      }
    }
  }  
  if(debug) std::cout << " " << std::endl;

}

void Geometry::constructFromParameters(bool debug) {
  // a tesselation of the plane with polygons
  double a(parameters_.cell_side);
  int nrows(parameters_.cells_nx);
  int ncols(parameters_.cells_ny);
  int klayer(parameters_.layer);
  Parameters::Geometry::Type itype(parameters_.type);

  if(debug) {
    std::cout << " " << std::endl;
    std::cout << "Building parametrized geometry : " << nrows << " " << ncols  << std::endl;
    // itype=0: heaxgonal cells
    // itype=1: triangular cells
    if (itype==Parameters::Geometry::Type::Hexagons) std:: cout << "with hexagonal cells " << std::endl;
    else if (itype==Parameters::Geometry::Type::Triangles) std:: cout << "with triangular cells " << std::endl;

    // itype=0: heaxgonal cells
    // itype=1: triangular cells
    if (itype==Parameters::Geometry::Type::Hexagons) std:: cout << "with hexagonal cells " << std::endl;
    else if (itype==Parameters::Geometry::Type::Triangles) std:: cout << "with triangular cells " << std::endl;
  }
  
  setNrows(nrows);
  setNcols(ncols);
  setLayer(klayer);
  setType(itype);
  a_ = a;
  asqrt3_ = a*std::sqrt(3.);
  asqrt3over2_ = asqrt3_/2.;
  aover2_ = a/2.;

  // vertices coordinates wrt cell center
  const int nverticeshexagon=6; // hexagons
  const int nverticestriangle=3; // hexagons  
  const std::array<double, nverticeshexagon> hexagonoffsetx = {asqrt3over2_,asqrt3over2_,0.,-asqrt3over2_,-asqrt3over2_,0};
  const std::array<double,nverticeshexagon> hexagonoffsety = {-aover2_,aover2_,a_,aover2_,-aover2_,-a_};
  const std::array<double,nverticestriangle> uptriangleoffsetx = {aover2_,0.,-aover2_};
  const std::array<double, nverticestriangle> uptriangleoffsety = {-asqrt3over2_/3.,asqrt3_/3.,-asqrt3over2_/3.};
  const std::array<double, nverticestriangle> downtriangleoffsetx = {aover2_,-aover2_,0.};
  const std::array<double,nverticestriangle> downtriangleoffsety = {asqrt3over2_/3.,asqrt3over2_/3.,-asqrt3_/3.};
  int nvertices=nverticeshexagon;
  if (itype==Parameters::Geometry::Type::Triangles) nvertices=nverticestriangle;
  cells_.clear();
  //cells_.reserve(nrows*ncols);

  double zlayer;
  if (klayer == -1) zlayer = 0.;  // entry face required
  else zlayer = parameters_.layers_z[klayer]; // else offset from the layer z position
  // to force being at the center for any requested layer
  //zlayer = 0.; 
  setZlayer(zlayer);

  // I index run along x axis is defined such that hexagons are adjacent by side along this axis
  // J index runs along the y' axis is rotated by 60deg wrt x axis  

  // define an offset along the x-axis to fully fill the display pad with cells
  double xoffset = parameters_.offset*asqrt3_;
  if (itype==Parameters::Geometry::Type::Triangles) xoffset = parameters_.offset*aover2_;

 for (int i=0; i<nrows_;i++) {
    for (int j=0; j<ncols_;j++) {
      if (debug) {
        std::cout << "Creating new cell of type " << static_cast<std::underlying_type<Parameters::Geometry::Type>::type>(itype) << " : " << std::endl;
        std::cout << " mapping coordinates : " << 
          i << " " <<
          j << std::endl;
      }
      double x=0., y=0.;
      switch(itype) {
        case Parameters::Geometry::Type::Hexagons:
        {
          x = xoffset + i*asqrt3_ + j*asqrt3over2_;
          double yprime = j*asqrt3_;
          // get back to the orthogonal y axis
          y = yprime*asqrt3over2_/a_;
          break;
        }
        case Parameters::Geometry::Type::Triangles:
        {
          x = xoffset + i*aover2_ + j*aover2_;
          y = j*asqrt3over2_;
          if (i%2 == 1) y = y + asqrt3_/6.; // cell center is shifted in y for downward triangles
          break;
        }
      };

      TVectorD position(2);
      position(0) = x;
      position(1) = y;
      if (debug) {
        std::cout << " center coordinates : " << 
          x << " " <<
          y << std::endl;
      }
      double orientation = 90.;
      if (itype == Parameters::Geometry::Type::Triangles && i%2 != 0) orientation =  -90.; // for downward triangles
      if (debug) {
        std::cout << " orientation : " << 
        orientation << std::endl;
      }
      std::vector<TVectorD> vertices;
      for (int iv=0; iv<nvertices; iv++) {
        vertices.emplace_back(2);
        switch(itype) {
          case Parameters::Geometry::Type::Hexagons:
          {
            vertices.back()(0) = x+hexagonoffsetx[iv];
            vertices.back()(1) = y+hexagonoffsety[iv];
            break;
          }
          case Parameters::Geometry::Type::Triangles:
          {
            switch(i%2) {
              case 0:
              {
                vertices.back()(0) = x+uptriangleoffsetx[iv];
                vertices.back()(1) = y+uptriangleoffsety[iv];
                break;
              }
              default:
              {
                vertices.back()(0) = x+downtriangleoffsetx[iv];
                vertices.back()(1) = y+downtriangleoffsety[iv];
                break;
              }
            }
            break;
          }
        };
        if (debug) {
          std::cout << "  vertex " << i << " " << 
            vertices.back()(0) << " " <<
            vertices.back()(1) << std::endl;
        }
      }

      cells_.emplace(
          Cell::id(i, j),
          Cell(std::move(position),std::move(vertices),orientation,i,j)
          );

    }
  }
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
  if (klayer_ != -1) std::cout << "the layer plane is " << klayer_ << " at z position " << parameters_.layers_z[klayer_] << std::endl;
  else std::cout << "the layer plane is " << klayer_ << " at z position 0." << std::endl;
  for (const auto& id_cell : cells_) { 
    const auto& cell = id_cell.second;
    std::cout << "new cell with indices " <<
      "("<< cell.getIIndex() << "," << cell.getJIndex() << ")" <<  
      " and position " << 
      "("<<cell.getPosition()(0) << "," << cell.getPosition()(1) << ")" << 
      std::endl;   
  }  

}

void Geometry::draw(const Parameters::Display& params, double scale) {

  std::array<double, 7> summitx, summity;

  // for the case of the full geometry from json, offset the display by (0.5,0.5) such that
  // the center of the module is at the center of the pad
  double xdisplayoffset=0.;
  double ydisplayoffset=0.;
  if (itype_==Parameters::Geometry::Type::External) {
    xdisplayoffset=params.offset_x;
    ydisplayoffset=params.offset_y;
  }

  for (const auto& id_cell : cells_) { 
    const auto& cell = id_cell.second;
    unsigned i=0;
    for (const auto& vertex : cell.getVertices()) {
      summitx[i]=vertex(0)*scale+xdisplayoffset;
      summity[i]=vertex(1)*scale+ydisplayoffset;
      i++;
    }
    unsigned nvertices = cell.getVertices().size();
    summitx[nvertices]=cell.getVertices()[0](0)*scale+xdisplayoffset;
    summity[nvertices]=cell.getVertices()[0](1)*scale+ydisplayoffset;
    TPolyLine polygon(nvertices+1,summitx.data(),summity.data());
    polygon.SetFillColor(38);
    polygon.SetLineColor(4);
    polygon.SetLineWidth(1);
    polygon.Draw();
  } 

}

const TVectorD& Geometry::getPosition(int i, int j) const {

  //TVectorD position(2);
  //// for parameterised geometries uses indices
  //// this is error prone (code duplication), do we really gain time?
  //if (getType()!=Parameters::Geometry::Type::External) { 
    //if (getType()==Parameters::Geometry::Type::Hexagons) { // hexagons
      //position(0) = parameters_.offset*asqrt3_+i*asqrt3_+j*asqrt3over2_;
      //double yprime = j*asqrt3_;
      //position(1) = yprime*asqrt3over2_/a_;
    //} else { // triangles
      //position(0) = parameters_.offset*asqrt3_+i*aover2_+j*aover2_;
      //position(1) = j*asqrt3over2_;    
      //if (i%2 == 1) position(1) = position(1) + asqrt3_/6.; // cell center is shifted in y for downward triangles
    //}
  //} else { // full geometry
  //for (const auto& id_cell : cells_) { 
    //if (cell.getIIndex()!=i) continue; 
    //if (cell.getJIndex()!=j) continue; 
    //return cell.getPosition();
  //}  
  //throw std::string("Didn't find cell");
  // This will throw an exception if the cell is not found
  return cells_.at(Cell::id(i,j)).getPosition();
}

const Cell& Geometry::closestCell(double x, double y) const {

  // beware that this function does not explicit check that the point is within
  // the cell, therefore if the cell grid is too small it will return
  // the closest cell (as the name indicates)

  // here we implement the a temporary bruteforce function for the more general geometry read from the json file
  // this needs to be optimize for timing efficiency

  // to test time spent here
  //return *cells_.begin();
 
  double r2min = std::numeric_limits<double>::max();
  auto icfound = cells_.cend();
  for (auto ic=cells_.cbegin();ic!=cells_.cend();ic++) { 
    double xcell = ic->second.getPosition()(0);     
    double ycell = ic->second.getPosition()(1);
    // now compute the distances and take the smallest one
    double r2 = (xcell-x)*(xcell-x) + (ycell-y)*(ycell-y);
    if (r2<r2min) {
      r2min=r2;
      icfound=ic;
    }
  } 

  if (icfound==cells_.cend()) {
    // FIXME: do we want to throw an exception?
    throw std::string("Didn't find closest cell");
    //std::cout << "[Geometry::closestCell] Cell not found!! x, y " << 
      //x << " " << y << std::endl;
  }

  return icfound->second;
}

TVectorD Geometry::positionInCell(const TVectorD& position) const {

  TVectorD relativeposition(2);
  relativeposition=position-closestCell(position(0),position(1)).getPosition();
  return relativeposition;

}

bool Geometry::isInCell(const TVectorD& position, const Cell& cell) const {

  // implementation below works for any convex cell described by its vertices
  // assumes vertices are consecutive along the cell perimeter and ordered along direct rotation

  // loop on pair of consective vertices 
  const auto& vertices = cell.getVertices();
  for (unsigned int i=0;i<vertices.size()-1;i++) {
    double xa = vertices[i](0);
    double ya = vertices[i](1);
    double xb = vertices[i+1](0);
    double yb = vertices[i+1](1);
    double sign = (xb-xa)*(position(1)-ya) - (yb-ya)*(position(0)-xa);
    if (sign<0.) return false;
  }  

  return true;
}

//int Geometry::getIIndex(const Cell& cell) const {

  //return cell.getIIndex();

//}

//int Geometry::getJIndex(const Cell& cell) const {

  //return cell.getJIndex();

//}

