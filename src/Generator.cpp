#include <iostream>
#include "Generator.h"
#include "TRandom3.h"
#include "Constants.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFile.h"
#include "TText.h"
#include "TPaveText.h"

using namespace Constants;

void Generator::simulate(int nevents) {
  // output file
  std::string fileName = "hgcal_shower_simulation.root";
  std::string histoFileName = "hist_"+fileName;

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
    geometry.constructFromParameters(nx,ny,klayer,itype); // constructor for layer klayer
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
  //double textsize=0.02;
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

  std::string hName; 
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
    std::cout << "================ Simulating event: " << iev << " ================" << std::endl;    

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
      if (!isincell) { 
        std::cout << "[main] point is not inside the closest cell!  x,y" << x << " " << y << 
          " cell position " << cell.getPosition()(0) << " " << cell.getPosition()(1) << 
          " closest cell indices " << geometry.getIIndex(cell) << " " <<  geometry.getJIndex(cell)<< std::endl;
      }

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
    hCellEnergyMap[**ic]->Write();
  }  

  hFile.Write();
  hFile.Close();

  t.Stop();
  t.Print();

  display(geometry,hCellEnergyMap);    

}



void Generator::display(Geometry& geometry, std::map<Cell,TH1F*,CellComp>& hCellEnergyEvtMap, int ievt) {
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
  if (geometry.getType()==1) scale=1./(nxdisplay*aover2);
  //double textsize=0.02;
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
