#ifndef Parameters_h__
#define Parameters_h__

#include <map>
#include <string>
#include <boost/python.hpp>

class Parameters
{
  public:
    struct General
    {
      General();
      unsigned events;
      bool debug;
    };
    struct Geometry
    {
      enum class Type {Hexagons, Triangles, External, Undefined};
      const static std::map<std::string, Type> type_map_;
      Geometry();
      Type type;
      int layer;
      std::vector<double> layers_z;
      // internal infinit geometries
      double side_length;
      int cells_x;
      int cells_y;
      // external geometries
      std::string file;
    };
    struct Generation
    {
      Generation();
      bool fluctuation;
      double energy;
      int number_of_hits_per_gev;
      double alpha;
      double mip_energy;
      double position;
      std::vector<double> layers_energy;

    };
    struct Display
    {
      Display();
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
