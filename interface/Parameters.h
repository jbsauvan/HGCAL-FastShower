#ifndef __HGCalSimulation_FastShower_Parameters_h__
#define __HGCalSimulation_FastShower_Parameters_h__

#include <map>
#include <string>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>

class Parameters
{
  public:
    struct General
    {
      General();
      unsigned events;
      bool debug;
      std::string output_file;
    };
    struct Geometry
    {
      enum class Type {Hexagons, Triangles, External, Undefined};
      const static std::map<std::string, Type> type_map_;
      Geometry();
      Type type;
      int layer;
      std::vector<double> layers_z;
      // internal infinite geometries
      double cell_side;
      double eta_min;
      double eta_max;
      double phi_min;
      double phi_max;
      // external geometries
      std::string file;
    };
    struct Shower 
    {
      Shower();
      double moliere_radius;
      double radiation_length;
      double critical_energy;
      std::map<std::string, double> transverse_parameters;
      std::map<std::string, double> longitudinal_parameters;
      std::vector<double> layers_energy;
      double alpha;
    };
    struct Generation
    {
      Generation();
      bool fluctuation;
      double energy;
      int number_of_hits_per_gev;
      double mip_energy;
      double sampling;
      bool noise;
      double noise_sigma;
      double incident_eta;
      double incident_phi;
    };
    struct Display
    {
      Display();
      unsigned events;
    };

  public:
    Parameters();
    ~Parameters() {};

    void read(const std::string&);
    void print() const;

    const General& general() const {return general_;}
    const Geometry& geometry() const {return geometry_;}
    const Shower& shower() const {return shower_;}
    const Generation& generation() const {return generation_;}
    const Display& display() const {return display_;}

  private:
    void fillGeneral(const boost::python::dict& dict);
    void fillGeometry(const boost::python::dict& dict);
    void fillShower(const boost::python::dict& dict);
    void fillGeneration(const boost::python::dict& dict);
    void fillDisplay(const boost::python::dict& dict);

    // Converters to std objects
    template<typename T>
    std::vector<T>
    toStdVector(const boost::python::list& plist)
    {
      boost::python::stl_input_iterator<T> begin(plist), end;
      return std::vector<T>(begin,end);
    }
    template<typename Key, typename Value>
    std::map<Key,Value>
    toStdMap(const boost::python::dict& pdict)
    {
      std::map<Key, Value> output;
      boost::python::list keys = pdict.keys();
      for(unsigned i=0; i<len(keys);++i)
      {
        Key key = boost::python::extract<Key>(keys[i]);
        Value value = boost::python::extract<Value>(pdict[key]);
        output.emplace(key, value);
      }
      return output;
    }


    General general_;
    Geometry geometry_;
    Shower shower_;
    Generation generation_;
    Display display_;

};

#endif
