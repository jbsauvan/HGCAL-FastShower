#include <iostream>
#include "Generator.h"
#include "ShowerParametrization.h"
#include "ShowerShapeHexagon.h"
#include "ShowerShapeTriangle.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TText.h"
#include "TPaveText.h"


Generator::
Generator(const Parameters& params):
  geometry_(params.geometry()),
  parameters_(params)
{
}

void Generator::simulate() {
  unsigned nevents = parameters_.general().events;
  bool debug = parameters_.general().debug;

  // output file
  std::string fileName = "hgcal_shower_simulation.root";
  std::string histoFileName = "hist_"+fileName;

  // some initializations
  double energygen=0.;
  double energygenincells=0.;
  double energyrec=0.;

  // randome engine
  //TRandom *gun = new TRandom(); 
  TRandom *gun = new TRandom3(); 

  // incident direction
  // coordinate of origin of simulated geometry/module in CMS frame is given by etainc,phiinc and z=320.;
  // set phiinc to 0., can be updated when we will have seval modules
  double phiinc = 0.;
  double thetainc = 2.*std::atan(std::exp(-parameters_.generation().incident_eta));
  double z0 = 320.; // z ccordinate of first plane
  double rt = z0*tan(thetainc); 
  TVectorD dir(3);  
  dir(0) = rt*cos(phiinc);
  dir(1) = rt*sin(phiinc);
  dir(2) = z0;

  TStopwatch t;
  t.Start();   

  if (parameters_.geometry().type!=Parameters::Geometry::Type::External) {
    geometry_.constructFromParameters(parameters_.general().debug); // constructor for layer klayer
    std::cout << " " << std::endl;
  } else {
    // construct geometry from JSON
    geometry_.constructFromJson(parameters_.general().debug);
  }
  geometry_.print();

  // draw the geometry
  std::string title;
  //char str[20];
  title = "Layer ";
  title = title + std::to_string(parameters_.geometry().layer);
  TCanvas *c1 = new TCanvas(title.c_str(),title.c_str(),40,40,700,700);
  double scale=1.;
  //double textsize=0.02;
  geometry_.draw(parameters_.display(), scale);

  // shower parametrization
  ShowerParametrization aShowerParametrization(parameters_);
  
//   // layer weight
//   double layer_weight=1.;
//   if (klayer!=-1) {
//     double total_weight = 0.;
//     for (int i=0;i<28;i++) total_weight = total_weight + elayers[i];
//     layer_weight = elayers[klayer]/total_weight;
//   }
  // moved inside the event loop
  //if (klayer!=-1) layer_weight = aShowerParametrization.getLayerProfile()[klayer];

  // energy map
  std::map<Cell*,double,CellComp> enrjMap;

  // book the histograms  
  TH1F hTransverseProfile("hTransverseProfile","Generated transverse profile (cm)",100,0.,20.);
  TH1F hPhiProfile("hPhiProfile","Generated azimuthal profile (cm)",100,0.,6.3);
  TH1F hEnergyGen("hEnergyGen","Generated total energy",100,0.,300.);
  TH1F hSpotEnergy("hSpotEnergy","Generated spot energy",100,0.,0.1);
  TH1F hLandau("hLandau","Landau fluctations",100,0.,5.);
  TH1F hCellEnergyDist("hCellEnergyDist","Cell energy",1000,0.,100.);
  TH1F hEnergySum("hEnergySum","Energy sum",120,0.,120.);

  std::map<Cell*, TH1F*, CellComp> hCellEnergyMap;
  std::map<Cell*, TH1F*, CellComp> hCellEnergyEvtMap;

  std::string hName; 
  std::vector<Cell *>::iterator ic;
  std::vector<Cell *> cells=geometry_.getCells();
  for (ic=cells.begin();ic!=cells.end();ic++) { 
    int i = (*ic)->getIIndex();
    int j = (*ic)->getJIndex();
    hName="hCellEnergy[";
    hName += std::to_string(i);
    hName += ",";
    hName += std::to_string(j);
    hName += "]";
    hCellEnergyMap[*ic] = new TH1F(hName.c_str(),"Energy in cell [i,j])",100,0.,100.);
  }

  if (parameters_.display().events>0) {
    std::vector<Cell *>::iterator ic;
    std::vector<Cell *> cells=geometry_.getCells();
    for (ic=cells.begin();ic!=cells.end();ic++) { 
      int i = (*ic)->getIIndex();
      int j = (*ic)->getJIndex();
      hName="hCellEnergyEvt[";
      hName += std::to_string(i);
      hName += ",";
      hName += std::to_string(j);
      hName += "]";
      hCellEnergyEvtMap[*ic] = new TH1F(hName.c_str(),"Event Energy in cell [i,j])",100,0.,100.);
    } 
  }

  std::cout << " " << std::endl;
  if (debug)
  {
    std::cout << "incident position: " <<"("<<
      parameters_.generation().incident_x<<","<<
      parameters_.generation().incident_y<<")"<<std::endl; 
  }
  if (debug) std::cout << "incident energy: " <<parameters_.generation().energy<<" GeV"<<std::endl; 
  //if (debug) std::cout << "energy in layer: " <<energy*layer_weight<<" GeV"<<std::endl; 
  if (debug) std::cout << "requested layer: " <<parameters_.geometry().layer<<std::endl; 

  if (debug)
  {
    std::cout<< "cell grid: " <<"("<<
      parameters_.geometry().cells_nx<<","<<
      parameters_.geometry().cells_ny<<")"<< std::endl;
  }
  if (debug) std::cout<< "hexagon side: " <<parameters_.geometry().cell_side<< std::endl;

  //if (debug) std::cout<< "moliere radius (at layer 15): " << 2.3*r0layer15 << " cm" << std::endl;
  if (debug) std::cout<< "moliere radius: " << aShowerParametrization.getMoliereRadius()  << " cm" << std::endl;
  if (debug) std::cout<< "nbr hits per GeV: " << parameters_.generation().number_of_hits_per_gev << std::endl;

  if (debug) std::cout<< "requested events: " << nevents << std::endl;

  // start main loop on all events
  for (unsigned iev=1; iev<=nevents; iev++) {

    // generate new event
    std::cout << "================ Simulating event: " << iev << " ================" << std::endl;    

    // initialize energies
    std::vector<Cell *> cells=geometry_.getCells();
    //for (ic=cells.begin();ic!=cells.end();ic++) enrjMap[**ic]=0.;
    for (ic=cells.begin();ic!=cells.end();ic++) enrjMap[*ic]=0.;

    energygen = 0.;
    energygenincells = 0.;
    energyrec = 0.;
    
    if (parameters_.generation().noise) {
      double calibratednoise = parameters_.generation().noise_sigma*parameters_.generation().mip_energy/parameters_.generation().sampling;
      for (ic=cells.begin();ic!=cells.end();ic++) {
        double enoise = gun->Gaus(0.,calibratednoise);
	enrjMap[*ic] += enoise;
	energyrec += enoise;
	if (debug) std::cout << "adding " << enoise << " GeV noise in cell " << std::endl;      
      }
    }
      
    // generate energy spots according to transverse profile
    double r0 = aShowerParametrization.r0(parameters_.geometry().layer);
    // take longitudinal profile as mean energy per layer for fixed energy
    double layer_weight = 1.;
    if (parameters_.geometry().layer!=-1) layer_weight = aShowerParametrization.getLayerProfile()[parameters_.geometry().layer];

    // energy spot 
    // no fluctuations: fixed energy = 1. / nhitspergev
    // fluctuations: alpha/sqrt(E) -> Poissonian nbr hits of energy 1/alpha^2
    // where alpha is the stochastic term of the resolution
    int nhits;
    double denrj;
    if (!parameters_.generation().fluctuation) {
      nhits = int(parameters_.generation().energy*layer_weight*parameters_.generation().number_of_hits_per_gev);
      denrj = 1./parameters_.generation().number_of_hits_per_gev;
    } else {
      denrj = aShowerParametrization.spotEnergy();
      nhits = gun->Poisson(parameters_.generation().energy*layer_weight/denrj); 
    }  

    if (debug) 
      std::cout << " number of generated hits " << nhits << " with energy " << denrj <<
        std::endl;

    // incident position
    double xinccor = parameters_.generation().incident_x;
    //double xinccor = xinc - asqrt3over2 + asqrt3*gun->Rndm();
    if (debug) std::cout << "shooting position = ("<< xinccor <<","<<parameters_.generation().incident_y<<")"<<std::endl;

    for (int i=0; i<nhits; i++) {

      double r = gun->Exp(r0); // exponential exp(-r/r0)
      double phi = gun->Rndm()*TMath::TwoPi();
      double x = r*cos(phi) + xinccor;
      double y = r*sin(phi) + parameters_.generation().incident_y;

      // add here translation for the requested layer
      double z = geometry_.getZlayer();
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
          <<"("<<geometry_.getIIndex(*geometry_.closestCell(x,y))
          <<","<<geometry_.getJIndex(*geometry_.closestCell(x,y))
          <<")"<<" cell(x,y) "
          <<"("<<geometry_.closestCell(x,y)->getPosition()(0)
          <<","<<geometry_.closestCell(x,y)->getPosition()(1)
          <<")"
          <<" isincell(cell) "<<geometry_.isInCell(pos,*geometry_.closestCell(x,y))
          <<" position in cell " 
          <<"("<<geometry_.positionInCell(pos)(0)
          <<","<<geometry_.positionInCell(pos)(1)
          <<")"
          <<std::endl;
      }

      energygen += denrj;

      // map generated point into geometry
      Cell *cell = geometry_.closestCell(x,y);

      // for half-cell or boarder cells, check it is within the cell
      bool isincell = geometry_.isInCell(pos,*cell);
      if (!isincell) { 
        std::cout << "[main] point is not inside the closest cell!  x,y" << x << " " << y << 
          " cell position " << cell->getPosition()(0) << " " << cell->getPosition()(1) << 
          " closest cell indices " << geometry_.getIIndex(*cell) << " " <<  geometry_.getJIndex(*cell)<< std::endl;
      }

      // add energy to corresponding cell
      if (isincell) {
        enrjMap[cell] += denrj;    
	energygenincells += denrj;  
        energyrec += denrj; 
      }	

     // fill shower histograms
      hTransverseProfile.Fill(r,denrj);
      hPhiProfile.Fill(phi,denrj);
      hSpotEnergy.Fill(denrj);

    }

    std::cout << "simulated energy " << energygen << std::endl;
    std::cout << "simulated energy inside cells " << energygenincells << std::endl;
    std::cout << "reconstructed energy inside cells (includes noise) " << energyrec << std::endl;
    
    ShowerShape *aShowerShape;
    if (parameters_.geometry().type!=Parameters::Geometry::Type::Triangles) { // hexagons
      aShowerShape = new ShowerShapeHexagon(&enrjMap);
    } else { // triangles
      aShowerShape = new ShowerShapeTriangle(&enrjMap);   
    }  
    std::cout << "cell max i,j " << aShowerShape->maxCell()->getIIndex() << " " << aShowerShape->maxCell()->getJIndex()
    << " with energy " << aShowerShape->maxE1() << std::endl;
    std::cout << "energy in first neighboors " << aShowerShape->firstNeighboors() << std::endl;

    // fill histograms
    hEnergyGen.Fill(energygen,1.);
    hEnergySum.Fill(energyrec,1.);

    for (ic=cells.begin();ic!=cells.end();ic++) { 
      hCellEnergyMap[*ic]->Fill(enrjMap[*ic]);
    }

    if (debug) std::cout << " incident energy: " << parameters_.generation().energy << " simulated energy: " << energygen << std::endl;    

    // if requested display a few events
    if (iev<=nevents) {
      for (ic=cells.begin();ic!=cells.end();ic++) { 
        hCellEnergyEvtMap[*ic]->Reset();
        hCellEnergyEvtMap[*ic]->Fill(enrjMap[*ic]);
      }
      display(hCellEnergyEvtMap,iev);    

    }

  }

  std::cout << std::endl;
  std::cout << nevents << " events generated " << std::endl;  
  std::cout << std::endl;

  // Exporting histograms to file
  TFile hFile(histoFileName.c_str(),"RECREATE");

  hEnergyGen.Write();
  hTransverseProfile.Write();
  hPhiProfile.Write();
  hSpotEnergy.Write();
  hCellEnergyDist.Write();
  hEnergySum.Write();

  for (ic=cells.begin();ic!=cells.end();ic++) { 
    hCellEnergyMap[*ic]->Write();
  }  

  hFile.Write();
  hFile.Close();

  t.Stop();
  t.Print();

  display(hCellEnergyMap);    

}



void Generator::display(std::map<Cell*,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt) {
  double xdisplayoffset=0.;
  double ydisplayoffset=0.;
  if (parameters_.geometry().type==Parameters::Geometry::Type::External) {
    xdisplayoffset=parameters_.display().offset_x;
    ydisplayoffset=parameters_.display().offset_y;
  }

  std::string title1, title2, title3, title4;
  char str[20];
  if (ievt == 0) title1 = "Mean energy profile in layer ";
  else title1 = "Event " + std::to_string(ievt) + " energy profile in layer ";
  title1 = title1 + std::to_string(parameters_.geometry().layer);
  title2 = "E = ";
  sprintf(str,"%4.1f",parameters_.generation().energy);
  std::string string=str;
  title2 = title2 + string;
  title2 = title2 + " GeV", 
         title3 = "position = (";
  sprintf(str,"%3.1f",parameters_.generation().incident_x);
  string=str;
  title3 = title3 + string;
  title3 = title3 + ",";
  sprintf(str,"%3.1f",parameters_.generation().incident_y);
  string=str;  
  title3 = title3 + string;
  title3 = title3 + ") cm";
  title4 = "eta = ";
  sprintf(str,"%3.1f",parameters_.generation().incident_eta);
  string=str;   
  title4 = title4 + string;
  std::string title = title1 + ", ";
  title = title + title2;
  title = title + ", ";
  title = title + title3;
  title = title + ", ";
  title = title + title4;

  TCanvas *c1 = new TCanvas(title.c_str(),title.c_str(),40,40,700,700);
  double scale=1./(parameters_.display().size*geometry_.asqrt3());
  if (parameters_.geometry().type==Parameters::Geometry::Type::Triangles) scale=1./(parameters_.display().size*geometry_.aover2());
  if (parameters_.geometry().type==Parameters::Geometry::Type::External) scale=1./(parameters_.display().size*0.22727*sqrt(3.)); //hard coded a size for json geometry, to be improved
  //double textsize=0.02;
  geometry_.draw(parameters_.display(), scale);

  std::map<Cell*,TH1F*,CellComp>::iterator ic;
  for (ic=hCellEnergyEvtMap.begin(); ic!=hCellEnergyEvtMap.end(); ic++) {
    // print mean energies
    double enrj = (ic->second)->GetMean();
    if (enrj<0.1) continue;
    sprintf(str,"%4.1f",(ic->second)->GetMean());
    //if (enrj<0.01) continue;
    //int ires = sprintf(str,"%5.2f",(ic->second)->GetMean());
    TText *t = new
      TText((ic->first)->getPosition()(0)*scale+parameters_.display().offset_x,(ic->first)->getPosition()(1)*scale+parameters_.display().offset_y,str);
    t->SetTextAlign(22);
    t->SetTextColor(kBlack);
    if (enrj>=1.) t->SetTextColor(kRed);
    t->SetTextFont(43);
    t->SetTextSize(20*11/parameters_.display().size);
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
