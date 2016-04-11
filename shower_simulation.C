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
// 08/04: drawing of the energy sharing
// 08/04: drawing of a module, aligned case
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
#include <cstdio>

using namespace std;
using namespace TMath;

// steering 
int nevents = 1000.;
//int nevents = 0;
const bool debug = false;
//const bool debug = true;
  
const double energy = 100.;
const double nhitspergev = 10.; // nbr hits per GeV 

// geometry
const int nx=31; // number of cells along x
const int ny=31; // number of cells along y
const int klayer=-1; // layer depth
const double a = 0.649; // the hexagon side in cm for large cells from SiHexGeometry.pdf 
//const double a = 0.476; // the hexagon side in cm for small cells from SiHexGeometry.pdf 
const double asqrt3 = a*sqrt(3.);
const double asqrt3over2 = asqrt3/2.;
const double aover2 = a/2.;

// display
const double nxdisplay=15; // the number of hexagons to display along the x axis 
const int ioffset=-9; // <0 starting cell along x, this is to fully fill the display window

// shooting position and direction
// for the non projective case this is the position at HGCAL entry (z=320cm) 
// center of the (4,8) cell 
const int iinc = 4-ioffset;
const int jinc = 8;
const double xinc = iinc*asqrt3+jinc*asqrt3over2;
const double yinc = jinc*asqrt3*asqrt3over2/a;
const double etainc=2.;

// layers' z positions for the non projective case 
// from CMSSW V7 geometry: https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
// the values are the silicon (centre) z positions of the 28 layers wrt HGCAL z entry position in cm 
const double zlayers[28] = {.765,1.515,2.745,3.495,4.725,5.475,6.705,7.455,8.685,9.435,
                       10.745,11.615,12.925,13.795,15.105,15.975,17.285,18.155,19.465,20.335,
                       21.785,22.855,24.305,25.375,26.825,27.895,29.345,30.415};

// shower parameters:
// transverse profile described by an exponential at a given depth
// exponential parameter set from TP studies, 90% containment in 2.3cm at layer 15
const double r0layer15=2.3/std::log(10.); 
// evolution vs depth described by a parabolic function
// set from TP studies (AMM), accuracy better than 5%
const double a0=9.-(18./63.); const double a1=135./630.; const double a2=45./630.;
  

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
  
  TVectorD *position; // centre position
   
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
  (*position)(0)=ioffset*asqrt3+i*asqrt3+j*asqrt3over2;
  double yprime=j*asqrt3;
  (*position)(1)=yprime*asqrt3over2/a;

}

class Geometry {

  // a tesselation of the plane with hexagonal cells
  // by default the plane is at z=0 but can be shifted to a z position that
  // corresponds to a layer position from TP geometry 

  public:
  
  Geometry(int nrows=11,int ncols=11,int klayer=-1); //default grid 11x11, ie +-5 around the central cell
  ~Geometry() {}
  
  void setLayer(int klayer);
  std::vector<Cell> getCells() {return cells_;}
  int getNumberOfRows() {return nrows_;}
  int getNumberOfCols() {return ncols_;}
  
  TVectorD getPosition(int i, int j); // position of cell i,j
  TVectorD getPosition(Cell cell) {return cell.getPosition();} // position of cell
  //vect<Cell> getNeighbours (int i, radius r); // all neighbours within radius r
  //vect<Cell> getFirstNeighbours (int i); // all neighbousrs within first anulus
  //vect<Cell> getSecondNeighbours (int i); // all neighbousrs within second anulus
  
  Cell closestCell(double x, double y); // the cell that contains the point
  TVectorD positionInCell(TVectorD position); // relative position within the cell
  bool isInCell(TVectorD position, Cell cell); // test if a point is within a cell
  //bool isInRealCell(TVectorD position, Cell cell); // to apply further mouse bite or dead region within the the hexagon virtual cell
  int getIIndex(Cell cell);
  int getJIndex(Cell cell);
  int getLayer() {return klayer_;}
  double getZlayer() {return zlayer_;}
 
  void draw(double scale=0.1);
  void print();
  
  private:
  
  std::vector<Cell> cells_;
  int nrows_;
  int ncols_;
  int klayer_;
  double zlayer_;
  
  TMatrixD rotation60;
   
};

Geometry::Geometry(int nrows,int ncols,int klayer):rotation60(2,2),nrows_(nrows),ncols_(ncols) {

  // a tesselation of the plane with hexagons

  // define the layer plane
  klayer_ = klayer;
  if (klayer<-1 || klayer>=28) {
    std::cout << "[Geometry::Geometry] error, invalid klayer " << klayer << std::endl;
    std::cout << "[Geometry::Geometry] setting klayer to -1 (entry face)" << std::endl;
    klayer_ = -1;
  } 
  if (klayer_ == -1) zlayer_ = 0.;  // entry face required
  else zlayer_ = zlayers[klayer_]; // else offset from the layer z position
  
  // I index run along x axis is defined such that hexagons are adjacent by side along this axis
  // J index runs along the y' axis is rotated by 60deg wrt x axis  
  
  // define an offset along the x-axis to fully fill the display pad with hexagons
  double xoffset = ioffset*asqrt3;
  
  for (int i=0; i<nrows_;i++) {
  
    for (int j=0; j<ncols_;j++) {
    
      double x = xoffset + i*asqrt3 + j*asqrt3over2;
      double yprime = j*asqrt3;
      // get back to the orthogonal y axis
      double y = yprime*asqrt3over2/a;
      
      //std::cout << "cell i,j,x,y " << i << " " << j << " " << x << " " << y << std::endl;
      cells_.push_back(Cell(x,y));
    
    }
  
  }

  double data[4] = {0.5, asqrt3over2/a, -asqrt3over2/a, 0.5};
  rotation60.SetMatrixArray(data);

}

void Geometry::setLayer(int klayer) {

  if (klayer<-1 || klayer>=28) {
    std::cout << "[Geometry::Geometry] error, invalid klayer " << klayer << std::endl;
    std::cout << "[Geometry::Geometry] setting klayer to -1 (entry face)" << std::endl;
    klayer_ = 0;
  }

  // entry face required
  if (klayer == -1) zlayer_ = 0.;
  // else offset from the layer z position
  else zlayer_ = zlayers[klayer];
  
}
  
void Geometry::print() {

  std::cout << "printing the geometry: " << std::endl;
  if (klayer_ != -1) std::cout << "the layer plane is " << klayer_ << " at z position " << zlayers[klayer_] << std::endl;
  else std::cout << "the layer plane is " << klayer_ << " at z position 0." << std::endl;
  std::vector<Cell>::iterator ic;
  for (ic=cells_.begin();ic!=cells_.end();ic++) { 
     std::cout << "new cell with indices " <<
     "("<< getIIndex(*ic) << "," << getJIndex(*ic) << ")" <<  
     " and position " << 
     "("<<(ic->getPosition())(0) << "," << (ic->getPosition())(1) << ")" << 
     std::endl;   
  }  

}

void Geometry::draw(double scale) {

  double summitx[7]; double summity[7];
  double offsetx[7]={-asqrt3over2,0.,asqrt3over2,asqrt3over2,0.,-asqrt3over2,-asqrt3over2};
  double offsety[7]={aover2,a,aover2,-aover2,-a,-aover2,aover2};
  for (ic=cells_.begin();ic!=cells_.end();ic++) { 
    for (int i=0;i<7;i++) summitx[i]=(ic->getPosition()(0)+offsetx[i])*scale;
    for (int i=0;i<7;i++) summity[i]=(ic->getPosition()(1)+offsety[i])*scale;
    TPolyLine *hexagon = new TPolyLine(7,summitx,summity);
    hexagon->SetFillColor(38);
    hexagon->SetLineColor(4);
    hexagon->SetLineWidth(1);
    //hexagon->Draw("f");
    hexagon->Draw();   
  } 
  
  // now draw a module of 11 cells size, aligned with central cell
  // later the module properties are to be added to the geometry object
  int nxmod=11; // number of cells along module x-axis, TP geometry
  // module centred on cell closest to (0.5,0.5) in display frame
  int jcenter=TMath::Nint(nxdisplay/2);
  int icenter=jcenter-TMath::Nint(jcenter*asqrt3over2)-ioffset;
  Cell center(icenter,jcenter);
  double offsetxmod[7]={-asqrt3over2*nxmod,0.,asqrt3over2*nxmod,asqrt3over2*nxmod,0.,
   -asqrt3over2*nxmod,-asqrt3over2*nxmod};
  double offsetymod[7]={aover2*nxmod,a*nxmod,aover2*nxmod,-aover2*nxmod,
   -a*nxmod,-aover2*nxmod,aover2*nxmod};
  for (int i=0;i<7;i++) summitx[i]=(center.getPosition()(0)+offsetxmod[i])*scale;
  for (int i=0;i<7;i++) summity[i]=(center.getPosition()(1)+offsetymod[i])*scale;
  TPolyLine *hexagon = new TPolyLine(7,summitx,summity);
  hexagon->SetFillColor(38);
  hexagon->SetLineColor(4);
  hexagon->SetLineWidth(3);
  //hexagon->Draw("f");
  hexagon->Draw();   
   

}

TVectorD Geometry::getPosition(int i, int j) {

  TVectorD position(2);
  position(0) = ioffset*asqrt3+i*asqrt3+j*asqrt3over2;
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
      if (fabs(xcell-x)>4*asqrt3) continue;
      double ycell = getPosition(i,j)(1);
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

  if ((ifound==-999 || jfound==-999)) 
   std::cout << "[Geometry::closestCell] Cell not found!! x, y " << 
    x << " " << y << std::endl;
   
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

void display(Geometry& geometry, TH1F *hCellEnergy[nx][ny]) {

  std::string title1, title2, title3, title4;
  char str[20];
  title1 = "Energy profile in layer ";
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
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) { 
       // print mean energies
       double enrj = hCellEnergy[i][j]->GetMean();
       if (enrj<0.1) continue;
       int ires = sprintf(str,"%4.1f",hCellEnergy[i][j]->GetMean());
       TText *t = new
        TText(geometry.getPosition(i,j)(0)*scale,geometry.getPosition(i,j)(1)*scale,str);
       t->SetTextAlign(22);
       t->SetTextColor(kBlack);
       if (enrj>=10.) t->SetTextColor(kRed);
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
 
}

void shower_simulation()
{

  // output file
  string fileName = "hgcal_shower_simulation.root";
  string histoFileName = "hist_"+fileName;
     
  // some initializations
  double energysum=0.;
  // incident direction
  TVectorD dir(3);
  double thetainc = 2.*std::atan(std::exp(-etainc));
  dir(0) = xinc;
  dir(1) = yinc;
  double rinc = sqrt(xinc*xinc+yinc*yinc);
  dir(2) = rinc / std::tan(thetainc);
  
  //Geometry geometry(nx,ny); // constructor for z=0 => HGCAL front face
  Geometry geometry(nx,ny,klayer); // constructor for layer klayer
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
      //std::cout << "hName " << hName << std::endl;
      hCellEnergy[i][j] = new TH1F(hName.c_str(),"Energy in cell [i,j])",100,0.,100.);
    }
  }

  std::cout << " " << std::endl;
  if (debug) std::cout << "incident position: " <<"("<<xinc<<","<<yinc<<")"<<std::endl; 
  if (debug) std::cout << "incident energy: " <<energy<<" GeV"<<std::endl; 
  
  if (debug) std::cout<< "cell grid: " <<"("<<nx<<","<<ny<<")"<< std::endl;
  if (debug) std::cout<< "hexagon side: " <<a<< std::endl;
  
  if (debug) std::cout<< "moliere radius (at layer 15): " << 2.3*r0layer15 << " cm" << std::endl;
  if (debug) std::cout<< "nbr hits per GeV: " << nhitspergev << std::endl;
  
  if (debug) std::cout<< "requested events: " << nevents << std::endl;
  
  // randome engine
  //TRandom *gun = new TRandom(); 
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
    double r0=r0layer15; // shower average value if layer -1 requested
    if (klayer!=-1) r0=(a0+a1*klayer+a2*klayer*klayer)*r0layer15/28.;
    
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

      // add here z translation
      double z = geometry.getZlayer();
      if (z!=0.) {
        x = x + z*(dir)(0)/(dir)(2);
        y = y + z*(dir)(1)/(dir)(2);
      }

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

      energysum += denrj;

      // map generated point in hexagonal gemetry
      Cell cell = geometry.closestCell(x,y);

      // add energy to corresponding cell
      int iindex = geometry.getIIndex(cell);
      int jindex = geometry.getJIndex(cell);
      if (iindex<0 || iindex>=geometry.getNumberOfRows()) {
        if (debug) std::cout << "cell iindex out of range " << iindex << std::endl;
	continue;
      }	
      if (jindex<0 || iindex>=geometry.getNumberOfCols()) {
        if (debug) std::cout << "cell jindex out of range " << jindex << std::endl;
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
  
  display(geometry,hCellEnergy);
  
}

int main(){
  shower_simulation();
}

