#include <iostream>
#include <string>

#ifdef STANDALONE
#include "Generator.h"
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Generator.h"
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif




int main(int argc, char** argv) {

  if(argc!=2)
  {
    std::cout<<"Usage:\n";
    std::cout<<" shower_simulation.exe config_file_name\n";
    return 1;
  }

  std::string config_file(argv[1]);
  Parameters params;
  try
  {
    params.read(config_file);
  }
  catch(const std::string& excep)
  {
    std::cout<<"An error occured while reading the python config file:\n";
    std::cout<<excep<<"\n";
    return 2;
  }
  catch(const boost::python::error_already_set& excep)
  {
    std::cout<<"An error occured while reading the python config file:\n";
    PyErr_Print();
    PyErr_Clear();
    return 2;
  }
  catch(...)
  {
    std::cout<<"An undefined error occured while reading the python config file\n";
    return 2;
  }

  if(params.general().debug) params.print();

  try
  {
    Generator generator(params);
    generator.simulate();
  }
  // FIXME: improve exception management
  catch(const std::string& excep)
  {
    std::cout<<"An error occured while initializing the generator:\n";
    std::cout<<excep<<"\n";
    return 3;
  }
  catch(...)
  {
    std::cout<<"An error occured while initializing the generator:\n";
    return 3;
  }
  return 0;
}
