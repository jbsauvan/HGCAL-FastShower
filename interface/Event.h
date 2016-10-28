#ifndef __HGCalSimulation_FastShower_Event_h__
#define __HGCalSimulation_FastShower_Event_h__


#include <unordered_map>

class Event 
{
  private:
    Event() {}; // No default constructor, use Event(run, event)

  public:
    Event(uint32_t, uint32_t);
    ~Event() {};

    void fillHit(uint32_t, double);
    void setGenerated(double, double, double);
    void clear();

    const uint32_t run() const {return run_;}
    const uint32_t event() const {return event_;}
    const double generatedEnergy() const {return generated_energy_;}
    const double generatedEta() const {return generated_eta_;}
    const double generatedPhi() const {return generated_phi_;}
    const std::unordered_map<uint32_t, double>& hits() const {return hits_;}

  private:
    uint32_t run_;
    uint32_t event_;
    // FIXME: improve generated information
    double generated_energy_;
    double generated_eta_;
    double generated_phi_;
    std::unordered_map<uint32_t, double> hits_;



};

#endif
