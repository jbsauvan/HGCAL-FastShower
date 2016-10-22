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
    void clear();

    const uint32_t run() const {return run_;}
    const uint32_t event() const {return event_;}
    const std::unordered_map<uint32_t, double>& hits() const {return hits_;}

  private:
    uint32_t run_;
    uint32_t event_;
    std::unordered_map<uint32_t, double> hits_;



};

#endif
