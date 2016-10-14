#ifndef Parameters_h__
#define Parameters_h__

#include <boost/python.hpp>

class Parameters
{
  public:
    struct General
    {
      unsigned events;
      bool debug;
    };
    struct Geometry
    {
      std::string type;
      double side_length;
      int cells_x;
      int cells_y;
      int layer;
    };
    struct Generation
    {
      bool fluctuation;
      double energy;
      int number_of_hits_per_gev;
      double alpha;
      double mip_energy;
      double position;

    };
    struct Display
    {
      unsigned events;
      int size;
    };

  public:
    Parameters();
    ~Parameters() {};

    void read(const std::string&);
    void print() const;

    const General& general() const {return general_;}
    const Geometry& geometry() const {return geometry_;}
    const Generation& generation() const {return generation_;}
    const Display& display() const {return display_;}

  private:
    void fillGeneral(boost::python::dict& dict);
    void fillGeometry(boost::python::dict& dict);
    void fillGeneration(boost::python::dict& dict);
    void fillDisplay(boost::python::dict& dict);


    General general_;
    Geometry geometry_;
    Generation generation_;
    Display display_;

};

#endif
