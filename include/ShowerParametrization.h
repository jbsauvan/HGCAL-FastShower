
#ifndef ShowerParametrization_h__
#define ShowerParametrization_h__


class ShowerParametrization {

   // Electromagnetic shower parametrization 
   // values from G. Grindhammer et al., hep-ex/0001020
   // adapted for HGCAL TP geometry

  public:

    ShowerParametrization() {}
    ShowerParametrization(double radiationLength, double moliereRadius, double criticalEnergy, double alpha):
     radiationLength_(radiationLength), moliereRadius_(moliereRadius), criticalEnergy_(criticalEnergy),
     alpha_(alpha) 
    {
      // longitudinal parameters, from TP studies https://indico.cern.ch/event/353342/contributions/830862/attachments/699977/961072/20141118CharlotHGCAL.ppt.pdf
      meant0_=-1.396;
      meant1_=1.007;
      meanalpha0_=-0.0433;
      meanalpha1_=0.540;
      sigmalnt0_=-2.506;
      sigmalnt1_=1.245;
      sigmalnalpha0_=-0.08442;
      sigmalnalpha1_=0.7904;
      corrlnalphalnt0_=0.7858;
      corrlnalphalnt1_=-0.0232;  
      // mean layer energy profile, taken from TP studies for35 GeV Pt electrons
      // to simulate 28-layers V7 geometry, group the first two layers and drop the last one according to:
      // https://indico.cern.ch/event/458374/contribution/9/attachments/1179028/1828217/Andreev_29Oct2015.pdf 
      double elayers[28] = {40.0,69.8,119.6,178.9,248.8,315.1,382.0,431.6,477.7,
                              498.7,533.6,514.8,490.0,435.1,386.7,325.4,277.9,224.4,186.5,
                              145.3,108.7,73.7,52.1,33.0,22.5,13.1,8.6,4.8};
      double total_weight=0.;
      for (int i=0;i<28;i++) total_weight = total_weight + elayers[i];
      for (int i=0;i<28;i++) layerProfile_[i] = elayers[i]/total_weight;

      // transverse parameters
      // exponential parameter set from TP studies, 90% containment in 2.3cm at layer 15
      //r0layer15_ = 2.3/std::log(10.);  
      r0layer15_ = moliereRadius_/std::log(10.);  
      // evolution vs depth described by a parabolic function
      // set from TP studies (AMM), accuracy better than 5%
      a0_=9.-(18./63.); 
      a1_=135./630.; 
      a2_=45./630.;
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
    const double *getLayerProfile() const {return layerProfile_;}
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
    double layerProfile_[28];
    
    // transverse parametrisation
    double r0layer15_;
    double a0_;
    double a1_;
    double a2_;

    // for fluctuation
    double alpha_; // the stochastic coefiscient in GeV^1/2
    
};

#endif
