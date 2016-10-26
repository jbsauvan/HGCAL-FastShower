

#ifdef STANDALONE
#include "Event.h"
#else
#include "HGCalSimulation/FastShower/interface/Event.h"
#endif


Event::
Event(uint32_t run, uint32_t event):
  run_(run),
  event_(event)
{
}

void
Event::
fillHit(uint32_t id, double energy)
{
  auto itr_insert = hits_.emplace(id, energy);
  // if id already inserted, add the energy
  if(!itr_insert.second) itr_insert.first->second += energy;
}

void
Event::
setGenerated(double energy, double x, double y, double eta)
{
  generated_energy_ = energy;
  generated_x_ = x;
  generated_y_  = y;
  generated_eta_ = eta;
}

void
Event::
clear()
{
  hits_.clear();
}
