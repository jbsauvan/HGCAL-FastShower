////////////////////////////////////////////////////////////////////////////////
//
// Releases notes:
//
// 30/03: initial version with exponential profile and no fluctuations
// 31/03: added simulation for the small cells
// 01/04: transverse parameter set to 10mm, corresponding to 90% containment in 23mm (Rm) 
// 
////////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TSystem.h"
#include "Math/Vector2D.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom3.h"
#include "TMath.h"

#include <vector>
#include <string>


using namespace std;
using namespace TMath;

// steering
int nevents = 1000.;
//int nevents = 0;
const bool debug = false;
//const bool debug = true;
  
const double energy = 100.;
const double nhitspergev = 10.; // nbr hits per GeV

const double a = 0.649; // the hexagon side in cm, large cells from SiHexGeometry.pdf
//const double a = 0.476; // the hexagon side in cm, small cells from SiHexGeometry.pdf
const double asqrt3 = a*sqrt(3.);
const double asqrt3over2 = asqrt3/2.;

// shooting position
// center of the (5,5) cell
const double xinc = 5*asqrt3+5*asqrt3over2;
const double yinc = 5*asqrt3*asqrt3over2/a;
// shifted by 0.1cm along x and y
//const double xinc = (5*asqrt3+5*asqrt3over2)+0.1;
//const double yinc = (5*asqrt3*asqrt3over2/a)+0.1;

class Cell {

  // an hexagonal cell
  
  public:
  
  Cell(double, double);
  Cell(int, int);
  Cell(const Cell&);
  Cell& operator=(const Cell&);
  ~Cell() {delete position;}
  
  TVectorD getPosition() const {return *position;}
  
  private:
  
  TVectorD *position;
  
};

Cell::Cell(const Cell& cell){

  position = new TVectorD(2);
  (*position)(0)=(cell.getPosition())(0); 
  (*position)(1)=(cell.getPosition())(1);

}

Cell& Cell::operator=(const Cell& cell) { 

  if (this != &cell) {
    delete position;
    position = new TVectorD(2);
    (*position)(0)=(cell.getPosition())(0); 
    (*position)(1)=(cell.getPosition())(1);  
  }
  return *this;

}

Cell::Cell(double x, double y){

  position = new TVectorD(2);
  (*position)(0)=x; 
  (*position)(1)=y;

}

Cell::Cell(int i, int j) {

  position = new TVectorD(2);
  (*position)(0)=i*asqrt3+j*asqrt3over2;
  double yprime=j*asqrt3;
  (*position)(1)=yprime*asqrt3over2/a;

}

class Geometry {

  // a tesselation of the plane with hexagonal cells

  public:
  
  Geometry(int nrows=11,int ncols=11); // constructor, default grid 11x11, ie +-5 around the central cell
  ~Geometry() {}
  
  std::vector<Cell> getCells() {return cells_;}
  int getNumberOfRows() {return nrows_;}
  int getNumberOfCols() {return ncols_;}
  
  TVectorD getPosition(int i, int j); // position of cell i
  TVectorD getPosition(Cell cell) {return cell.getPosition();} // position of cell i
  //vect<Cell> getNeighbours (int i, radius r); // all neighbours within radius r
  //vect<Cell> getFirstNeighbours (int i); // all neighbousrs within first anulus
  //vect<Cell> getSecondNeighbours (int i); // all neighbousrs within second anulus
  
  Cell closestCell(double x, double y); // the cell that contains the point
  TVectorD positionInCell(TVectorD position); // relative position within the cell
  bool isInCell(TVectorD position, Cell cell); // test if a point is within a cell
  //bool isInRealCell(TVectorD position, Cell cell); // to apply further mouse bite or dead region within the the hexagon virtual cell
  int getIIndex(Cell cell);
  int getJIndex(Cell cell);
  
  void print();
  
  private:
  
  std::vector<Cell> cells_;
  int nrows_;
  int ncols_;
  
  TMatrixD rotation60;
   
};

Geometry::Geometry(int nrows,int ncols):rotation60(2,2),nrows_(nrows),ncols_(ncols) {

  // a tesselation of the plane with hexagons
  // I index run along x axis is defined such that hexagons are adjacent by side along this axis
  // J index runs along the y' axis is rotated by 60deg wrt x axis
  
  for (int i=0; i<nrows_;i++) {
  
    for (int j=0; j<ncols_;j++) {
    
      double x = i*asqrt3 + j*asqrt3over2;
      double yprime = j*asqrt3;
      // get back to the orthogonal y axis
      double y = yprime*asqrt3over2/a;
      
      //std::cout << "cell i,j,x,y " << i << " " << j << " " << x << " " << y << std::endl;
      cells_.push_back(Cell(x,y));
    
    }
  
  }

  //TArrayD data(4);
  double data[4] = {0.5, asqrt3over2/a, -asqrt3over2/a, 0.5};
  rotation60.SetMatrixArray(data);

}

void Geometry::print() {

  std::cout << "printing the geometry: " << std::endl;
  std::vector<Cell>::iterator ic;
  for (ic=cells_.begin();ic!=cells_.end();ic++) { 
     std::cout << "new cell with indices " <<
     "("<< getIIndex(*ic) << "," << getJIndex(*ic) << ")" <<  
     " and position " << 
     "("<<(ic->getPosition())(0) << "," << (ic->getPosition())(1) << ")" << 
     std::endl;   
  }  

}

TVectorD Geometry::getPosition(int i, int j) {

  TVectorD position(2);
  position(0) = i*asqrt3+j*asqrt3over2;
  double yprime = j*asqrt3;
  position(1) = yprime*asqrt3over2/a;
  return position;

}

Cell Geometry::closestCell(double x, double y) {

  // intialize to Cel(-999,-999) to identify unfound cells
  int ifound=-999;
  int jfound=-999;
  
  double r2min=99999.;
  for (int i=0; i<nrows_;i++) {
  
    for (int j=0; j<ncols_; j++) {
    
      // filter out far cells
      double xcell = getPosition(i,j)(0);     
      double ycell = getPosition(i,j)(1);
      if (fabs(xcell-x)>4*asqrt3) continue;
      if (fabs(ycell-y)>4*asqrt3) continue;
      
      // now compute the distances and take the smallest one
      double r2 = (xcell-x)*(xcell-x) + (ycell-y)*(ycell-y);
      if (r2<r2min) {
        r2min=r2;
	ifound=i;
	jfound=j;
      }
      
    }
    
  } 
  
  return Cell(ifound,jfound);   
 
}

TVectorD Geometry::positionInCell(TVectorD position) {

  TVectorD relativeposition(2);
  relativeposition=position-closestCell(position(0),position(1)).getPosition();
  return relativeposition;

}

bool Geometry::isInCell(TVectorD position, Cell cell) {

  TVectorD cellPosition = cell.getPosition();
  TVectorD relativeposition=position-cellPosition;
  double r2 = relativeposition(0)*relativeposition(0)+relativeposition(1)*relativeposition(1);
  
  // first eliminate distances bigger than the circle
  if (sqrt(r2)>a) return false;
  
  // the point is in the circle of radius a, now check the hexagon sides
  if (fabs(relativeposition(0))>asqrt3over2) return false;
  
  TVectorD relativeposition60 = rotation60*relativeposition;
  if (fabs(relativeposition60(0))>asqrt3over2) return false;
  
  TVectorD relativeposition120 = rotation60*relativeposition60;
  if (fabs(relativeposition120(0))>asqrt3over2) return false;
  
  return true;

}

int Geometry::getIIndex(Cell cell) {

  double x = cell.getPosition()(1)*a/asqrt3;
  x = cell.getPosition()(0) - x;
  return TMath::Nint(x/asqrt3);
  
}

int Geometry::getJIndex(Cell cell) {

  double yprime = cell.getPosition()(1)*a/asqrt3over2;
  return TMath::Nint(yprime/asqrt3);
  
}



void shower_simulation()
{

  // output file
  string fileName = "hgcal_shower_simulation.root";
  string histoFileName = "hist_"+fileName;
     
  // some initializations
  double energysum;
  
  // shower parameters
  // transverse profile described by an exponential
  // exponential parameter set from TP studies, 90% containment in 2.3cm
  const double r0=2.3/std::log(10.); 
  
  // geometry
  const int nx=11; // number of cells along x
  const int ny=11; // number of cells along y
  Geometry geometry(nx,ny);
  std::cout << " " << std::endl;
  geometry.print();
  
  // energy array
  double enrj[nx][ny];    
  
  // book the histograms
  
  TH1F hTransverseProfile("hTransverseProfile","Generated transverse profile (cm)",100,0.,20.);
  TH1F hPhiProfile("hPhiProfile","Generated azimuthal profile (cm)",100,0.,6.3);
  TH1F hEnergySum("hEnergySum","Generated total energy",100,0.,200.);
  TH1F hSpotEnergy("hSpotEnergy","Generated spot energy",100,0.,10.);
  
  TH1F *hCellEnergy[nx][ny];
  string hName; 
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      hName="hCellEnergy[";
      hName += std::to_string(i);
      hName += ",";
      hName += std::to_string(j);
      hName += "]";
      hCellEnergy[i][j] = new TH1F(hName.c_str(),"Energy in cell [i,j])",100,0.,100.);
    }
  }
  
  std::cout << " " << std::endl;
  if (debug) std::cout << "incident position: " <<"("<<xinc<<","<<yinc<<")"<<std::endl; 
  if (debug) std::cout << "incident energy: " <<energy<<" GeV"<<std::endl; 
  
  if (debug) std::cout<< "cell grid: " <<"("<<nx<<","<<ny<<")"<< std::endl;
  if (debug) std::cout<< "hexagon side: " <<a<< std::endl;
  
  if (debug) std::cout<< "moliere radius: " << 2.3*r0 << " cm" << std::endl;
  if (debug) std::cout<< "nbr hits per GeV: " << nhitspergev << std::endl;
  
  if (debug) std::cout<< "requested events: " << nevents << std::endl;
  
  // randome engine
  //TRandom *gun = new TRandom(); 
  //TRandom *gun = new TRandom1(); 
  //TRandom *gun = new TRandom2();  
  TRandom *gun = new TRandom3(); 
  
  // start main loop on all events
  for (int iev=1; iev<=nevents; iev++) {

    // generate new event
    if (debug) cout << "================ Generating event: " << iev << " ================" << endl;    
    
    // initialize energies
    for (int i=0;i<nx;i++) 
      for (int j=0;j<ny;j++)
        enrj[i][j] = 0.;
    energysum = 0.;
     	
    // generate energy spots according to transverse profile
    // each energy spot has fixed energy = 1. / nhitspergev
    // later implement the sampling fluctuations
    // 25%/sqrt(E) -> epsilon = 0.0625 ls /binGeV and 16 points per GeV
    int nhits = int(energy*nhitspergev);
    double denrj = energy / nhits;
    if (debug) 
     std::cout << " number of generated hits " << nhits << " with energy " << denrj <<
     std::endl;
     
    for (int i=0; i<nhits; i++) {

      double r = gun->Exp(r0); // exponential exp(-r/r0)
      double phi = gun->Rndm()*TMath::TwoPi();
      double x = r*cos(phi) + xinc;
      double y = r*sin(phi) + yinc;
      //double x = xinc;
      //double y = yinc;
      if (debug) {
        TVectorD pos(2);
	pos(0)=x;
	pos(1)=y;
        std::cout << " new generated hit with energy " << denrj << " and position(x,y) " 
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
      
      hTransverseProfile.Fill(r,denrj);

      energysum += denrj;

      // map generated point in hexagonal gemetry
      Cell cell = geometry.closestCell(x,y);

      // add energy to corresponding cell
      int iindex = geometry.getIIndex(cell);
      int jindex = geometry.getJIndex(cell);
      if (iindex<0 && iindex>=geometry.getNumberOfRows()) {
        std::cout << "cell iindex out of range " << iindex << std::endl;
	continue;
      }	
      if (jindex<0 && iindex>=geometry.getNumberOfCols()) {
        std::cout << "cell jindex out of range " << jindex << std::endl;
	continue;
      }	

      enrj[iindex][jindex] += denrj;      
    
      // fill histograms
      hTransverseProfile.Fill(r,denrj);
      hPhiProfile.Fill(phi,denrj);
      hSpotEnergy.Fill(denrj);

    }
    
     // fill histograms
     hEnergySum.Fill(energysum,1.);
     for (int i=0; i<nx; i++) 
       for (int j=0; j<ny; j++) 
          hCellEnergy[i][j]->Fill(enrj[i][j]);
     if (debug) std::cout << " incident energy: " << energy << " simulated energy: " << energysum << std::endl;    

  }
    
  cout << endl;
  cout << nevents << " events generated " << endl;  
  cout << endl;
  
  // Exporting histograms to file
  TFile hFile(histoFileName.c_str(),"RECREATE");
  
  hEnergySum.Write();
  hTransverseProfile.Write();
  hPhiProfile.Write();
  hSpotEnergy.Write();
  for (int i=0; i<nx; i++) 
    for (int j=0; j<ny; j++) 
       hCellEnergy[i][j]->Write();   
  
  hFile.Write();
  hFile.Close();
  
 
}

int main(){
  shower_simulation();
}

