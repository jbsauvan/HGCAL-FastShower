
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
  cell_n_ = event.hits().size();
  for(const auto& id_hit : event.hits())
  {
    cell_energy_.emplace_back(id_hit.second);
    const auto& cell = geometry.getCells().at(id_hit.first);
    cell_x_.emplace_back(cell.getPosition()(0));
    cell_y_.emplace_back(cell.getPosition()(1));
  }
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
}



