////////////////////////////////////////////////////////////////////////////////
//
// Releases notes:
//
// 30/03: initial version with exponential profile and no fluctuations
// 31/03: option for small cells 
// 01/04: transverse parameter set to 10mm, corresponding to 90% containment in 23mm (Rm) 
// 06/04: layer structure in depth, no rotation between the different planes 
// 06/04: transverse profile varying with depth  
// 06/04: added incidence angle 
// 08/04: drawing of the mean energy sharing
// 08/04: drawing of a module, aligned case
// 11/04: small mods to make it compilable with aClic
// 13/04: added sampling fluctuations
// 13/04: added display of a individual events
// 13/04: added calculation of reconstructed asymetry, barycenter, shower width in x and y
// 20/04: random generation of incident position
// 03/05: energy fraction in each layer according to TP profile 
// 28/06: added possibility of reading the geometry from JSON description
// 30/09: full integration of the detailed/json geometry up to display and hit generation
// 04/10: full integration of the previous parameterised version 
// 
////////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TF1.h"
#include "TSystem.h"
#include "Math/Vector2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "TText.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TStopwatch.h"

#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <iostream>
#include <fstream>

#include "json.h"

using namespace std;
using namespace TMath;

// steering 
int nevents = 10;
//int nevents = 1;
const bool debug = false;
//const bool fluctuation = true;
const bool fluctuation = false;
  
const double energy = 100.;
const double nhitspergev = 1000.; // nbr hits per GeV, used only if fluctuation is false 
const double alpha = 0.25; // the sampling parameter, determines nhitspergev if fluctuation is true 
const double mipenergy = 0.000055; // mip energy for 200 microns Si

// geometry
const bool readgeom=true;
//const bool readgeom=false;

// all parameters bellow correspond to simplified parameterised geometry, ignored if readgeom=true
const int nx=31; // number of cells along x
const int ny=31; // number of cells along y
//const double a = 0.649; // the hexagon side in cm for large cells from SiHexGeometry.pdf 
const double a = 0.22727; // the hexagon side in cm for large cells from SiHexGeometry.pdf 
//const double a = 0.476; // the hexagon side in cm for small cells from SiHexGeometry.pdf 
const double asqrt3 = a*sqrt(3.);
const double asqrt3over2 = asqrt3/2.;
const double aover2 = a/2.;
const double offsetx[6] = {asqrt3over2,asqrt3over2,0.,-asqrt3over2,-asqrt3over2,0};
const double offsety[6] = {-aover2,aover2,a,aover2,-aover2,-a};

// importing a full geometry from json file, used if readgeom=true
const string geomfile="/data_CMS/cms/charlot/CMSSW/CMSSW_7_6_3/src/AC_11.json"; // aligned centered, d=11

// coordinates system for parameterised (infinite) geometry
const int ioffsetparam=-9; // <0 starting cell along x, this is to fully fill the display window

// layer depth
// if between 0 and 27, the generated energy is weighted according to the layer profile (elayers[i])
// the incident position is also shifted according to layer z position (zlayers[i]) and incident direction (etainc)
// if set to -1, the full generated energy is deposited without longitudinal weighting and the
// incident position is not shifted, that is it corresponds to layer 0 with the full energy deposited in it
const int klayer=-1; 

// display
const double nxdisplay=15; // the number of hexagons to display along the x axis (parameterised case)
const double xdisplayoffsetfull=0.5; // center of cell (0,0) x-shift (full geometry from json)
const double ydisplayoffsetfull=0.5; // center of cell (0,0) y-shift (full geometry from json)
const int nevtdisplay=5; // number of individual events to display

// shooting position and direction
// shooting position is the transverse position at HGCAL entry, ie at layer 0
// below we shoot at center of the (4,8) cell in the parameterised geometry
//const int iinc = 4-ioffsetparam; 
const int iinc = 4; 
const int jinc = 8;              
//const double xinc = iinc*asqrt3+jinc*asqrt3over2;
//const double yinc = jinc*asqrt3*asqrt3over2/a;
//  below we shoot at the border between the (4,8) and (5,8) cells in x and y
//const double epsilon=0.000;
//const double xinc = iinc*asqrt3+jinc*asqrt3over2+asqrt3over2-epsilon;
//const double yinc = jinc*asqrt3*asqrt3over2/a;
// for full geometry the origin is at (0,0)
const double xinc = 0.;
const double yinc = 0.;

const double etainc=2.;

// layers' z positions
// from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
// the values are the silicon (centre) z positions of the 28 layers wrt HGCAL z entry position in cm 
const double zlayers[28] = {.765,1.515,2.745,3.495,4.725,5.475,6.705,7.455,8.685,9.435,
                       10.745,11.615,12.925,13.795,15.105,15.975,17.285,18.155,19.465,20.335,
                       21.785,22.855,24.305,25.375,26.825,27.895,29.345,30.415};
// layers energy profile
// from TP geometry, average for 35 geV Pt electrons
// to simulate 28-layers V7 geometry, group the first two layers and drop the last one according to:
// https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
const double elayers[28] = {40.0,69.8,119.6,178.9,248.8,315.1,382.0,431.6,477.7,
                              498.7,533.6,514.8,490.0,435.1,386.7,325.4,277.9,224.4,186.5,
                              145.3,108.7,73.7,52.1,33.0,22.5,13.1,8.6,4.8};

// shower parameters:
// transverse profile described by an exponential at a given depth
// exponential parameter set from TP studies, 90% containment in 2.3cm at layer 15
const double r0layer15=2.3/std::log(10.); 
// evolution vs depth described by a parabolic function
// set from TP studies (AMM), accuracy better than 5%
const double a0=9.-(18./63.); const double a1=135./630.; const double a2=45./630.;

  
std::string fixedLength(int value, int digits = 5, string type="FH") {

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
 
class Cell {

  // an hexagonal cell
  
  public:
  
  Cell() {} // default constructor
  // constructors for parametrized geometries
  Cell(double, double);
  Cell(int, int);
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
  bool operator() (const Cell& c1, const Cell& c2) const 
  {
    if (c1.getIIndex()!=c2.getIIndex()) return (c1.getIIndex()<c2.getIIndex());
    else return (c1.getJIndex()<c2.getJIndex());
  }
};

Cell::Cell(const Cell& cell){

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
  void constructFromJson(const string&);
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

void Geometry::constructFromJson(const string& filename) {
  
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
  ifstream ifs(filename);
  bool success = reader.parse(ifs, obj);
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
  int ifound=-999;
  int jfound=-999;
  
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

    TVectorD *vertex = new TVectorD(2);
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

void display(Geometry& geometry, std::map<Cell,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt=0) {

  double xdisplayoffset=0.;
  double ydisplayoffset=0.;
  if (readgeom) {
    xdisplayoffset=xdisplayoffsetfull;
    ydisplayoffset=ydisplayoffsetfull;
  }

  std::string title1, title2, title3, title4;
  char str[20];
  if (ievt == 0) title1 = "Mean energy profile in layer ";
  else title1 = "Event " + std::to_string(ievt) + " energy profile in layer ";
  title1 = title1 + std::to_string(klayer);
  title2 = "E = ";
  int ires = sprintf(str,"%4.1f",energy);
  std::string string=str;
  title2 = title2 + string;
  title2 = title2 + " GeV", 
  title3 = "position = (";
  ires = sprintf(str,"%3.1f",xinc);
  string=str;
  title3 = title3 + string;
  title3 = title3 + ",";
  ires = sprintf(str,"%3.1f",yinc);
  string=str;  
  title3 = title3 + string;
  title3 = title3 + ") cm";
  title4 = "eta = ";
  ires = sprintf(str,"%3.1f",etainc);
  string=str;   
  title4 = title4 + string;
  std::string title = title1 + ", ";
  title = title + title2;
  title = title + ", ";
  title = title + title3;
  title = title + ", ";
  title = title + title4;
  
  TCanvas *c1 = new TCanvas(title.c_str(),title.c_str(),40,40,700,700);
  double scale=1./(nxdisplay*asqrt3);
  double textsize=0.02;
  geometry.draw(scale);

  std::map<Cell,TH1F*,CellComp>::iterator ic;
  for (ic=hCellEnergyEvtMap.begin(); ic!=hCellEnergyEvtMap.end(); ic++) {
    // print mean energies
    double enrj = (ic->second)->GetMean();
    if (enrj<0.1) continue;
    int ires = sprintf(str,"%4.1f",(ic->second)->GetMean());
    //if (enrj<0.01) continue;
    //int ires = sprintf(str,"%5.2f",(ic->second)->GetMean());
    TText *t = new
     TText((ic->first).getPosition()(0)*scale+xdisplayoffset,(ic->first).getPosition()(1)*scale+ydisplayoffset,str);
    t->SetTextAlign(22);
    t->SetTextColor(kBlack);
    if (enrj>=1.) t->SetTextColor(kRed);
    t->SetTextFont(43);
    t->SetTextSize(20*11/nxdisplay);
    t->Draw();
    TPaveText *leg1 = new TPaveText(.05,.91,.35,.97);
    leg1->AddText(title1.c_str());
    leg1->SetFillColor(kWhite);
    leg1->SetTextSize(0.02);
    leg1->Draw();
    TPaveText *leg2 = new TPaveText(.045,.85,.18,.88);
    leg2->AddText(title2.c_str());
    leg2->SetFillColor(kWhite);
    leg2->SetTextSize(0.02);
    leg2->SetTextColor(kBlue);
    leg2->SetBorderSize(0.0);
    leg2->Draw();
    TPaveText *leg3 = new TPaveText(.06,.79,.25,.84);
    leg3->AddText(title3.c_str());
    leg3->SetFillColor(kWhite);
    leg3->SetTextSize(0.02);
    leg3->SetTextColor(kBlue);
    leg3->SetBorderSize(0.0);
    leg3->Draw();
    TPaveText *leg4 = new TPaveText(.05,.76,.13,.79);
    leg4->AddText(title4.c_str());
    leg4->SetFillColor(kWhite);
    leg4->SetTextSize(0.02);
    leg4->SetTextColor(kBlue);
    leg4->SetBorderSize(0.0);
    leg4->Draw();
  } 
 
}

void shower_simulation()
{

  // output file
  string fileName = "hgcal_shower_simulation.root";
  string histoFileName = "hist_"+fileName;
     
  // some initializations
  double energygen=0.;
  double energyrec=0.;
  
  // randome engine
  //TRandom *gun = new TRandom(); 
  TRandom *gun = new TRandom3(); 
  
  // incident direction
  // coordinate of origin of simulated geometry/module in CMS frame is given by etainc,phiinc and z=320.;
  // set phiinc to 0., can be updated when we will have seval modules
  double phiinc = 0.;
  double thetainc = 2.*std::atan(std::exp(-etainc));
  double z0 = 320.; // z ccordinate of first plane
  double rt = z0*tan(thetainc); 
  TVectorD dir(3);  
  dir(0) = rt*cos(phiinc);
  dir(1) = rt*sin(phiinc);
  dir(2) = z0;
  
  TStopwatch t;
  t.Start();   
    
  Geometry geometry;
  if (!readgeom) {
    //Geometry geometry(nx,ny); // constructor for z=0 => HGCAL front face
    geometry.constructFromParameters(nx,ny,klayer); // constructor for layer klayer
    std::cout << " " << std::endl;
   } else {
    // constructing geometry from JSON
    geometry.constructFromJson(geomfile);
  }
  geometry.print();
  
  // draw the geometry
  std::string title;
  char str[20];
  title = "Layer ";
  title = title + std::to_string(klayer);
  TCanvas *c1 = new TCanvas(title.c_str(),title.c_str(),40,40,700,700);
  double scale=1.;
  double textsize=0.02;
  geometry.draw(scale);
  
  // layer weight
  double layer_weight=1.;
  if (klayer!=-1) {
    double total_weight = 0.;
    for (int i=0;i<28;i++) total_weight = total_weight + elayers[i];
    layer_weight = elayers[klayer]/total_weight;
  }
   
  // energy map
  std::map<Cell,double,CellComp> enrjMap;
   
  // book the histograms  
  TH1F hTransverseProfile("hTransverseProfile","Generated transverse profile (cm)",100,0.,20.);
  TH1F hPhiProfile("hPhiProfile","Generated azimuthal profile (cm)",100,0.,6.3);
  TH1F hEnergyGen("hEnergyGen","Generated total energy",100,0.,300.);
  TH1F hSpotEnergy("hSpotEnergy","Generated spot energy",100,0.,0.1);
  TH1F hLandau("hLandau","Landau fluctations",100,0.,5.);
  TH1F hCellEnergyDist("hCellEnergyDist","Cell energy",1000,0.,100.);
  TH1F hEnergySum("hEnergySum","Energy sum",120,0.,120.);
  
  std::map<Cell, TH1F*, CellComp> hCellEnergyMap;
  std::map<Cell, TH1F*, CellComp> hCellEnergyEvtMap;
  
  string hName; 
  std::vector<Cell *>::iterator ic;
  std::vector<Cell *> cells=geometry.getCells();
  for (ic=cells.begin();ic!=cells.end();ic++) { 
    int i = (*ic)->getIIndex();
    int j = (*ic)->getJIndex();
    hName="hCellEnergy[";
    hName += std::to_string(i);
    hName += ",";
    hName += std::to_string(j);
    hName += "]";
    //std::cout << "hName " << hName << std::endl;
    hCellEnergyMap[**ic] = new TH1F(hName.c_str(),"Energy in cell [i,j])",100,0.,100.);
  }
  
  if (nevtdisplay>0) {
    std::vector<Cell *>::iterator ic;
    std::vector<Cell *> cells=geometry.getCells();
    for (ic=cells.begin();ic!=cells.end();ic++) { 
      int i = (*ic)->getIIndex();
      int j = (*ic)->getJIndex();
      hName="hCellEnergyEvt[";
      hName += std::to_string(i);
      hName += ",";
      hName += std::to_string(j);
      hName += "]";
      //std::cout << "hName " << hName << std::endl;
      hCellEnergyEvtMap[**ic] = new TH1F(hName.c_str(),"Event Energy in cell [i,j])",100,0.,100.);
    } 
  }
  
  std::cout << " " << std::endl;
  if (debug) std::cout << "incident position: " <<"("<<xinc<<","<<yinc<<")"<<std::endl; 
  if (debug) std::cout << "incident energy: " <<energy<<" GeV"<<std::endl; 
  if (debug) std::cout << "energy in layer: " <<energy*layer_weight<<" GeV"<<std::endl; 
  
  if (debug) std::cout<< "cell grid: " <<"("<<nx<<","<<ny<<")"<< std::endl;
  if (debug) std::cout<< "hexagon side: " <<a<< std::endl;
  
  if (debug) std::cout<< "moliere radius (at layer 15): " << 2.3*r0layer15 << " cm" << std::endl;
  if (debug) std::cout<< "nbr hits per GeV: " << nhitspergev << std::endl;
  
  if (debug) std::cout<< "requested events: " << nevents << std::endl;
  
  // start main loop on all events
  for (int iev=1; iev<=nevents; iev++) {

    // generate new event
    cout << "================ Simulating event: " << iev << " ================" << endl;    
    
    // initialize energies
    std::vector<Cell *> cells=geometry.getCells();
    for (ic=cells.begin();ic!=cells.end();ic++) enrjMap[**ic]=0.;

    energygen = 0.;
    energyrec=0.;
     	
    // generate energy spots according to transverse profile
    // transverse profile at a given layer scaled according to TP results
    double r0=r0layer15; // shower average value if layer -1 requested
    if (klayer!=-1) r0=(a0+a1*klayer+a2*klayer*klayer)*r0layer15/28.;    
    
    // energy spot 
    // no fluctuations: fixed energy = 1. / nhitspergev
    // fluctuations: alpha/sqrt(E) -> Poissonian nbr hits of energy 1/alpha^2
    // where alpha is the stochastic term of the resolution
    int nhits;
    double denrj;
    if (!fluctuation) {
      nhits = int(energy*layer_weight*nhitspergev);
      denrj = 1./nhitspergev;
    } else {
      denrj = alpha*alpha;
      nhits = gun->Poisson(energy*layer_weight/denrj); 
    }  
    
    if (debug) 
     std::cout << " number of generated hits " << nhits << " with energy " << denrj <<
     std::endl;
     
    // incident position
    double xinccor = xinc;
    //double xinccor = xinc - asqrt3over2 + asqrt3*gun->Rndm();
    if (debug) std::cout << "shooting position = ("<< xinccor <<","<<yinc<<")"<<std::endl;

    for (int i=0; i<nhits; i++) {

      double r = gun->Exp(r0); // exponential exp(-r/r0)
      double phi = gun->Rndm()*TMath::TwoPi();
      double x = r*cos(phi) + xinccor;
      double y = r*sin(phi) + yinc;

      // add here translation for the requested layer
      double z = geometry.getZlayer();
      if (z!=0.) {
        x = x + z*(dir)(0)/(dir)(2);
        y = y + z*(dir)(1)/(dir)(2);
      }

      TVectorD pos(2);
      pos(0)=x;
      pos(1)=y;
      
      if (debug) {
        std::cout << " new simulated hit with energy " << denrj << " and position(x,y) " 
        <<"("<<x<<","<<y<<")"<< " cell(i,j) " 
        <<"("<<geometry.getIIndex(geometry.closestCell(x,y))
        <<","<<geometry.getJIndex(geometry.closestCell(x,y))
        <<")"<<" cell(x,y) "
        <<"("<<geometry.closestCell(x,y).getPosition()(0)
        <<","<<geometry.closestCell(x,y).getPosition()(1)
        <<")"
        <<" isincell(cell) "<<geometry.isInCell(pos,geometry.closestCell(x,y))
        <<" position in cell " 
        <<"("<<geometry.positionInCell(pos)(0)
        <<","<<geometry.positionInCell(pos)(1)
        <<")"
        <<std::endl;
      }

      energygen += denrj;

      // map generated point into hexagonal gemetry
      Cell cell = geometry.closestCell(x,y);
      
      // for half-cell or boarder cells, check it is within the cell
      bool isincell = geometry.isInCell(pos,cell);
      
      // add energy to corresponding cell
      if (isincell) 
       enrjMap[cell] += denrj;      
      if (isincell) 
       energyrec += denrj; 
	
      // fill shower histograms
      hTransverseProfile.Fill(r,denrj);
      hPhiProfile.Fill(phi,denrj);
      hSpotEnergy.Fill(denrj);
      
    }

     std::cout << "simulated energy " << energygen << std::endl;
     std::cout << "simulated energy inside cells " << energyrec << std::endl;
    
    // fill histograms
    hEnergyGen.Fill(energygen,1.);
    hEnergySum.Fill(energyrec,1.);
       
    for (ic=cells.begin();ic!=cells.end();ic++) { 
      hCellEnergyMap[**ic]->Fill(enrjMap[**ic]);
    }

    if (debug) std::cout << " incident energy: " << energy << " simulated energy: " << energygen << std::endl;    

    // if requested display a few events
    if (iev<=nevtdisplay) {
      for (ic=cells.begin();ic!=cells.end();ic++) { 
        hCellEnergyEvtMap[**ic]->Reset();
        hCellEnergyEvtMap[**ic]->Fill(enrjMap[**ic]);
      }
      display(geometry,hCellEnergyEvtMap,iev);    

    }

  }
    
  cout << endl;
  cout << nevents << " events generated " << endl;  
  cout << endl;
  
  // Exporting histograms to file
  TFile hFile(histoFileName.c_str(),"RECREATE");
  
  hEnergyGen.Write();
  hTransverseProfile.Write();
  hPhiProfile.Write();
  hSpotEnergy.Write();
  hCellEnergyDist.Write();
  hEnergySum.Write();

   for (ic=cells.begin();ic!=cells.end();ic++) { 
     hCellEnergyMap[**ic]->Write();
   }  
  
  hFile.Write();
  hFile.Close();
  
  t.Stop();
  t.Print();
  
  display(geometry,hCellEnergyMap);    
  
}

int main(){
  shower_simulation();
}

