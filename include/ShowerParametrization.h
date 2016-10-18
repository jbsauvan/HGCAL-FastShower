
#ifndef ShowerParametrization_h__
#define ShowerParametrization_h__

#include "Parameters.h"


class ShowerParametrization {

   // Electromagnetic shower parametrization 
   // values from G. Grindhammer et al., hep-ex/0001020
   // adapted for HGCAL TP geometry

  public:

    ShowerParametrization() {}
    ShowerParametrization(const Parameters& params):
     radiationLength_(params.shower().radiation_length),
     moliereRadius_(params.shower().moliere_radius),
     criticalEnergy_(params.shower().critical_energy),
     alpha_(params.shower().alpha) 
    {

      const auto& longitudinal_params = params.shower().longitudinal_parameters;
      // FIXME: throw excetion if parameter doesn't exist. Should make sure that it is properly caught
      meant0_ = longitudinal_params.at("meant0");
      meant1_ = longitudinal_params.at("meant1");
      meanalpha0_ = longitudinal_params.at("meanalpha0");
      meanalpha1_ = longitudinal_params.at("meanalpha1");
      sigmalnt0_ = longitudinal_params.at("sigmalnt0");
      sigmalnt1_ = longitudinal_params.at("sigmalnt1");
      sigmalnalpha0_ = longitudinal_params.at("sigmalnalpha0");
      sigmalnalpha1_ = longitudinal_params.at("sigmalnalpha1");
      corrlnalphalnt0_ = longitudinal_params.at("corrlnalphalnt0");
      corrlnalphalnt1_ = longitudinal_params.at("corrlnalphalnt1");  

      // FIXME: can compute the normalized profile in python
      if(params.shower().layers_energy.size()!=28) throw std::string("The size of shower_layers_energy should be 28");
      std::copy_n(params.shower().layers_energy.begin(), layerProfile_.size(), layerProfile_.begin());
      double total_weight=0.;
      for (unsigned i=0;i<layerProfile_.size();i++) total_weight = total_weight + layerProfile_[i];
      for (unsigned i=0;i<layerProfile_.size();i++) layerProfile_[i] = layerProfile_[i]/total_weight;
      // transverse parameters
      // exponential parameter set from TP studies, 90% containment in 2.3cm at layer 15
      //r0layer15_ = 2.3/std::log(10.);  
      r0layer15_ = moliereRadius_/std::log(10.);  
      const auto& transverse_params = params.shower().transverse_parameters;
      // FIXME: throw excetion if parameter doesn't exist. Should make sure that it is properly caught
      a0_ = transverse_params.at("a0"); 
      a1_ = transverse_params.at("a1");
      a2_ = transverse_params.at("a2");
    }
    ~ShowerParametrization() {}

    // average medium
    double getRadiationLength() const {return radiationLength_;}
    double getMoliereRadius() const {return moliereRadius_;}
    double getCriticalEnergy() const {return criticalEnergy_;}

    // longitudinal
    double meanT(double lny) const {return meant0_+meant1_*lny;}
    double meanAlpha(double lny) const {return meanalpha0_+meanalpha1_*lny;}
    double sigmaLnT(double lny) const {return 1./(sigmalnt0_+sigmalnt1_*lny);}
    double sigmaLnAlpha(double lny) const {return 1./(sigmalnalpha0_+sigmalnalpha1_*lny);}
    double correlationAlphaT(double lny) const {return corrlnalphalnt0_+corrlnalphalnt1_*lny;}        
    const std::array<double,28>& getLayerProfile() const {return layerProfile_;}
    // transversal
    double r0(int klayer) {return klayer==-1 ? r0layer15_ : (a0_ + a1_*klayer + a2_*klayer*klayer)*r0layer15_/28.;}

    // fluctuations
    double spotEnergy() {return alpha_*alpha_;}
    
  private:
  
    // average medium parameters
    double radiationLength_; // in cm
    double moliereRadius_;   // in cm 
    double criticalEnergy_;  // in MeV 

    // longitudinal parametrisation
    double meant0_;
    double meant1_;
    double meanalpha0_;
    double meanalpha1_;
    double sigmalnt0_;
    double sigmalnt1_;
    double sigmalnalpha0_;
    double sigmalnalpha1_;
    double corrlnalphalnt0_;
    double corrlnalphalnt1_;
    std::array<double,28> layerProfile_;
    
    // transverse parametrisation
    double r0layer15_;
    double a0_;
    double a1_;
    double a2_;

    // for fluctuation
    double alpha_; // the stochastic coefiscient in GeV^1/2
    
};

#endif
