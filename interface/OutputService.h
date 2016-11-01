#ifndef __HGCalSimulation_FastShower_OutputService_h__
#define __HGCalSimulation_FastShower_OutputService_h__


#include <vector>

#include "TFile.h"
#include "TTree.h"

#ifdef STANDALONE
#include "Event.h"
#include "Geometry.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#include "HGCalSimulation/FastShower/interface/Geometry.h"
#endif


class OutputService
{
  public:
    OutputService(const std::string&);
    ~OutputService();

    void fillTree(const Event&, const Geometry&);

  private:
    std::unique_ptr<TFile> file_;
    TTree* tree_; // the tree is owned by file_ and deleted when the file is closed

    // tree branches
    unsigned run_;
    unsigned event_;
    float gen_energy_;
    float gen_eta_;
    float gen_phi_;
    unsigned cell_n_;
    std::vector<float> cell_energy_;
    std::vector<float> cell_x_;
    std::vector<float> cell_y_;
    std::vector<float> cell_z_;
    std::vector<float> cell_eta_;
    std::vector<float> cell_phi_;

    void clear();


};

#endif
