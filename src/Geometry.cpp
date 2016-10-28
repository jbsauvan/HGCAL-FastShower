

#include <fstream>
#include <iostream>
#include <algorithm>

#include "TPolyLine.h"
#include "TVector2.h"


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
  int klayer(parameters_.layer);
  Parameters::Geometry::Type itype(parameters_.type);

  if(debug) {
    std::cout << " " << std::endl;
    std::cout << "Building parametrized geometry :\n";
    std::cout << "  " << parameters_.eta_min << " < eta < " << parameters_.eta_max << "\n";
    std::cout << "  " << parameters_.phi_min << " < phi < " << parameters_.phi_max << "\n";
    if (itype==Parameters::Geometry::Type::Hexagons) std:: cout << "with hexagonal cells " << std::endl;
    else if (itype==Parameters::Geometry::Type::Triangles) std:: cout << "with triangular cells " << std::endl;
  }
  
  setLayer(klayer);
  setType(itype);
  a_ = a;
  asqrt3_ = a*std::sqrt(3.);
  asqrt3over2_ = asqrt3_/2.;
  aover2_ = a/2.;
  a3over2_ = aover2_*3.;

  // vertices coordinates wrt cell center
  const int nverticeshexagon=6; // hexagons
  const int nverticestriangle=3; // hexagons  
  const std::array<double, nverticeshexagon> hexagonoffsetx = {{asqrt3over2_,asqrt3over2_,0.,-asqrt3over2_,-asqrt3over2_,0}};
  const std::array<double,nverticeshexagon> hexagonoffsety = {{-aover2_,aover2_,a_,aover2_,-aover2_,-a_}};
  const std::array<double,nverticestriangle> uptriangleoffsetx = {{aover2_,0.,-aover2_}};
  const std::array<double, nverticestriangle> uptriangleoffsety = {{-asqrt3over2_/3.,asqrt3_/3.,-asqrt3over2_/3.}};
  const std::array<double, nverticestriangle> downtriangleoffsetx = {{aover2_,-aover2_,0.}};
  const std::array<double,nverticestriangle> downtriangleoffsety = {{asqrt3over2_/3.,asqrt3over2_/3.,-asqrt3_/3.}};
  int nvertices=nverticeshexagon;
  if (itype==Parameters::Geometry::Type::Triangles) nvertices=nverticestriangle;
  cells_.clear();

  double zlayer;
  if (klayer == -1) zlayer = 0.;  // entry face required
  else zlayer = parameters_.layers_z[klayer]; // else offset from the layer z position
  // to force being at the center for any requested layer
  //zlayer = 0.; 
  setZlayer(zlayer);

  // I index run along x axis is defined such that hexagons are adjacent by side along this axis
  // J index runs along the y' axis is rotated by 60deg wrt x axis  

  // compute x,y positions of the geometry window
  double z = zlayer + 320.; // FIXME: remove hardcoded z0
  double theta_min = 2.*std::atan(std::exp(-parameters_.eta_max));
  double theta_max = 2.*std::atan(std::exp(-parameters_.eta_min));
  double r_min = z*tan(theta_min);
  double r_max = z*tan(theta_max);
  double phi_min = parameters_.phi_min;
  double phi_max = parameters_.phi_max;
  std::array<double, 6> xs = {
    r_min*cos(parameters_.phi_min),
    r_min*cos(parameters_.phi_max),
    r_max*cos(parameters_.phi_max),
    r_max*cos(parameters_.phi_min),
    r_min*cos((parameters_.phi_min+parameters_.phi_max)/2.),
    r_max*cos((parameters_.phi_min+parameters_.phi_max)/2.)
  };
  std::array<double, 6> ys = {
    r_min*sin(parameters_.phi_min),
    r_min*sin(parameters_.phi_max),
    r_max*sin(parameters_.phi_max),
    r_max*sin(parameters_.phi_min),
    r_min*sin((parameters_.phi_min+parameters_.phi_max)/2.),
    r_max*sin((parameters_.phi_min+parameters_.phi_max)/2.)
  };



  // x0,y0 is the origine of the tesselation
  // Which means that the center of cell (i,j)=(0,0) 
  // will be at (x,y)=(x0,y0)
  std::array<double, 5> dx = {
    xs[1]-xs[0],
    xs[2]-xs[0],
    xs[3]-xs[0],
    xs[4]-xs[0],
    xs[5]-xs[0]
  };
  std::array<double, 5> dy = {
    ys[1]-ys[0],
    ys[2]-ys[0],
    ys[3]-ys[0],
    ys[4]-ys[0],
    ys[5]-ys[0]
  };

  // partial derivatives of x,y vs i,j
  double dxdi = 0.;
  double dxdj = 0.;
  double dydj = 0.;
  switch(itype) {
    case Parameters::Geometry::Type::Hexagons:
      {
        dxdi = asqrt3_;
        dxdj = asqrt3over2_;
        dydj = a3over2_;
        break;
      }
    case Parameters::Geometry::Type::Triangles:
      {
        dxdi = aover2_;
        dxdj = aover2_;
        dydj = asqrt3over2_;
        break;
      }
    default:
      break;
  };
  // compute i,j window needed to cover the x,y window
  std::array<double,6> js = {
    0.,
    dy[0]/dydj,
    dy[1]/dydj,
    dy[2]/dydj,
    dy[3]/dydj,
    dy[4]/dydj
  };
  std::array<double,6> is = {
    0.,
    (dx[0]-js[1]*dxdj)/dxdi,
    (dx[1]-js[2]*dxdj)/dxdi,
    (dx[2]-js[3]*dxdj)/dxdi,
    (dx[3]-js[4]*dxdj)/dxdi,
    (dx[4]-js[5]*dxdj)/dxdi
  };
  int imin = std::round(*(std::min_element(is.begin(), is.end())));
  int jmin = std::round(*(std::min_element(js.begin(), js.end())));
  int imax = std::round(*(std::max_element(is.begin(), is.end())));
  int jmax = std::round(*(std::max_element(js.begin(), js.end())));

  double x_min = std::numeric_limits<double>::max();
  double x_max = std::numeric_limits<double>::lowest();
  double y_min = std::numeric_limits<double>::max();
  double y_max = std::numeric_limits<double>::lowest();
  // build cells inside the requested window
  for (int i=imin; i<=imax;i++) {
    for (int j=jmin; j<=jmax;j++) {
      double x = xs[0] + i*dxdi + j*dxdj;
      double y = ys[0] + j*dydj;
      // up and down triangle barycenters are not aligned
      if(itype==Parameters::Geometry::Type::Triangles && i%2) y += asqrt3_/6.;
      double r = std::sqrt(x*x + y*y);
      double phi = std::acos(x/r);
      // check if cell is inside boundaries. If not, skip it
      if(!(r>=r_min && r<=r_max &&
           TVector2::Phi_mpi_pi(phi-phi_min)>=0 && TVector2::Phi_mpi_pi(phi-phi_max)<=0
          ))
      {
        continue;
      }
      if(x>x_max) x_max = x;
      if(x<x_min) x_min = x;
      if(y>y_max) y_max = y;
      if(y<y_min) y_min = y;

      if (debug) {
        std::cout << "Creating new cell of type " << static_cast<std::underlying_type<Parameters::Geometry::Type>::type>(itype) << " : " << std::endl;
        std::cout << " mapping coordinates : " << 
          i << " " <<
          j << std::endl;
      }
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
          default:
            break;
        };
        if (debug) {
          std::cout << "  vertex " << iv << " " << 
            vertices.back()(0) << " " <<
            vertices.back()(1) << std::endl;
        }
      }

      auto found = cells_.emplace(
          Cell::id(i, j),
          Cell(std::move(position),std::move(vertices),orientation,i,j)
          );
      if(!found.second)
      {
        std::cout<<"Warning: Cell with indices"<<i<<" "<<j<<" already exists (id="<<Cell::id(i,j)<<")\n";
        std::cout<<"Warning: This may indicate a bug in the id definition\n";
      }

    }
  }
  // build histogram of cells
  cell_histogram_ = new TH2Poly("cells", "cells",
      x_min>0?x_min*0.9:x_min*1.1, x_max>0?x_max*1.1:x_max*0.9,
      y_min>0?y_min*0.9:y_min*1.1, y_max>0?y_max*1.1:y_max*0.9
      );
  for (const auto& id_cell : cells_) { 
    const auto& cell = id_cell.second;
    std::vector<double> binsx, binsy;
    for (const auto& vertex : cell.getVertices()) {
      binsx.emplace_back(vertex(0));
      binsy.emplace_back(vertex(1));
    }
    binsx.emplace_back(cell.getVertices()[0](0));
    binsy.emplace_back(cell.getVertices()[0](1));
    cell_histogram_->AddBin(binsx.size(), binsx.data(), binsy.data());
  }
}

void Geometry::setLayer(int klayer) {

  if (klayer<-1 || klayer>=28) {
    std::stringstream error;
    error << "[Geometry] error, invalid klayer " << klayer;
    throw error.str();
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

void Geometry::draw(const Parameters::Display& params) {

  cell_histogram_->Draw();

  //std::array<double, 7> summitx, summity;
  //for (const auto& id_cell : cells_) { 
    //const auto& cell = id_cell.second;
    //unsigned i=0;
    //for (const auto& vertex : cell.getVertices()) {
      //summitx[i]=vertex(0);
      //summity[i]=vertex(1);
      //i++;
    //}
    //unsigned nvertices = cell.getVertices().size();
    //summitx[nvertices] = cell.getVertices()[0](0);
    //summity[nvertices] = cell.getVertices()[0](1);
    //// Calling Draw makes the current pad take the ownership of the object
    //// So raw pointers are used, and the objects are deleted when the pad is deleted
    //TPolyLine* polygon = new TPolyLine(nvertices+1,summitx.data(),summity.data());
    //polygon->SetFillColor(38);
    //polygon->SetLineColor(4);
    //polygon->SetLineWidth(1);
    ////polygon->Draw();
  //} 

}

const TVectorD& Geometry::getPosition(int i, int j) const {
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


