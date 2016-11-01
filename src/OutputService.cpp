
#ifdef STANDALONE
#include "OutputService.h"
#else
#include "HGCalSimulation/FastShower/interface/OutputService.h"
#endif


OutputService::
OutputService(const std::string& file_name):
  file_(TFile::Open(file_name.c_str(), "recreate")),
  tree_(new TTree("tree", "tree")),
  run_(0),
  event_(0),
  cell_n_(0)
{
  tree_->Branch("run", &run_, "run/i");
  tree_->Branch("event", &event_, "event/i");
  tree_->Branch("gen_energy", &gen_energy_, "gen_energy/F");
  tree_->Branch("gen_eta", &gen_eta_, "gen_eta/F");
  tree_->Branch("gen_phi", &gen_phi_, "gen_phi/F");
  tree_->Branch("cell_n", &cell_n_, "cell_n/i");
  tree_->Branch("cell_energy", &cell_energy_);
  tree_->Branch("cell_x", &cell_x_);
  tree_->Branch("cell_y", &cell_y_);
  tree_->Branch("cell_z", &cell_z_);
  tree_->Branch("cell_eta", &cell_eta_);
  tree_->Branch("cell_phi", &cell_phi_);
}

OutputService::
~OutputService()
{
  tree_->Write();
  file_->Close();
}


void 
OutputService::
fillTree(const Event& event, const Geometry& geometry)
{
  clear();
  run_ = event.run();
  event_ = event.event();
  gen_energy_ = event.generatedEnergy();
  gen_eta_ = event.generatedEta();
  gen_phi_ = event.generatedPhi();
  for(const auto& id_hit : event.hits())
  {
    // skip zero and negative energies
    if(id_hit.second<=0.) continue;
    const auto& cell = geometry.getCells().at(id_hit.first);
    double x = cell.getPosition()(0);
    double y = cell.getPosition()(1);
    double z = cell.getPosition()(2);
    double r = std::sqrt(x*x + y*y);
    double theta = std::atan(r/z);
    double eta = -std::log(std::tan(theta/2.));
    double phi = std::copysign(std::acos(x/r),y);
    cell_energy_.emplace_back(id_hit.second);
    cell_x_.emplace_back(x);
    cell_y_.emplace_back(y);
    cell_z_.emplace_back(z);
    cell_eta_.emplace_back(eta);
    cell_phi_.emplace_back(phi);
  }
  cell_n_ = cell_energy_.size();
  tree_->Fill();
}


void
OutputService::
clear()
{
  run_ = 0;
  event_ = 0;
  gen_energy_ = 0.;
  gen_eta_ = 0.;
  gen_phi_ = 0.;
  cell_n_ = 0;
  cell_energy_.clear();
  cell_x_.clear();
  cell_y_.clear();
  cell_z_.clear();
  cell_eta_.clear();
  cell_phi_.clear();
}



