
#include <iostream>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>
#include "Parameters.h"


using namespace boost;

Parameters::
Parameters()
{
}

void
Parameters::
read(const std::string& file)
{
  Py_Initialize();
  python::object main_module = python::import("__main__");
  python::dict main_namespace = python::extract<python::dict>(main_module.attr("__dict__"));
  python::object py_module = python::exec_file(file.c_str(), main_namespace);
  fillGeneral(main_namespace);
  fillGeometry(main_namespace);
  fillGeneration(main_namespace);
  fillDisplay(main_namespace);
}

void 
Parameters::
fillGeneral(python::dict& dict)
{
  general_.events = python::extract<int>(dict["events"]);
  general_.debug = python::extract<bool>(dict["debug"]);
}

void 
Parameters::
fillGeometry(python::dict& dict)
{
  geometry_.type = python::extract<std::string>(dict["geometry_type"]);
  geometry_.side_length = python::extract<double>(dict["geometry_side_length"]);
  geometry_.cells_x = python::extract<int>(dict["geometry_cells_x"]);
  geometry_.cells_y = python::extract<int>(dict["geometry_cells_y"]);
  geometry_.layer = python::extract<int>(dict["geometry_layer"]);
  // Read vector of layers z
  python::list layers_z = python::extract<python::list>(dict["geometry_layers_z"]);
  python::stl_input_iterator<double> begin(layers_z), end;
  geometry_.layers_z = std::vector<double>(begin,end);
}

void 
Parameters::
fillGeneration(python::dict& dict)
{
  generation_.energy = python::extract<double>(dict["generation_energy"]);
  generation_.fluctuation = python::extract<bool>(dict["generation_fluctuation"]);
  generation_.number_of_hits_per_gev = python::extract<int>(dict["generation_number_of_hits_per_gev"]);
  generation_.alpha = python::extract<double>(dict["generation_alpha"]);
  generation_.mip_energy = python::extract<double>(dict["generation_mip_energy"]);
  // Read vector of layers energies
  python::list layers_energy = python::extract<python::list>(dict["generation_layers_energy"]);
  python::stl_input_iterator<double> begin(layers_energy), end;
  generation_.layers_energy = std::vector<double>(begin,end);
}

void 
Parameters::
fillDisplay(python::dict& dict)
{
  display_.events = python::extract<int>(dict["display_events"]);
  display_.size = python::extract<int>(dict["display_size"]);
}


void
Parameters::
print() const
{
  std::cout<<"|Configuration parameters\n";
  std::cout<<"|- General\n";
  std::cout<<"|-- Events = "<<general_.events<<"\n";
  std::cout<<"|-- Debug = "<<general_.debug<<"\n";
  std::cout<<"|- Geometry\n";
  std::cout<<"|-- Type = "<<geometry_.type<<"\n";
  std::cout<<"|-- SideLength = "<<geometry_.side_length<<"\n";
  std::cout<<"|-- Cells_X = "<<geometry_.cells_x<<"\n";
  std::cout<<"|-- Cells_Y = "<<geometry_.cells_y<<"\n";
  std::cout<<"|-- Layer = "<<geometry_.layer<<"\n";
  std::cout<<"|-- Layers z = [";
  for(const auto& z : geometry_.layers_z)
  {
    std::cout<<z<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|- Generation\n";
  std::cout<<"|-- Energy = "<<generation_.energy<<"\n";
  std::cout<<"|-- Fluctuation = "<<generation_.fluctuation<<"\n";
  std::cout<<"|-- Hits/GeV = "<<generation_.number_of_hits_per_gev<<"\n";
  std::cout<<"|-- Alpha = "<<generation_.alpha<<"\n";
  std::cout<<"|-- Mip energy = "<<generation_.mip_energy<<"\n";
  std::cout<<"|-- Position = "<<generation_.position<<"\n";
  std::cout<<"|-- Layers energy = [";
  for(const auto& e : generation_.layers_energy)
  {
    std::cout<<e<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|- Display\n";
  std::cout<<"|-- Events = "<<display_.events<<"\n";
  std::cout<<"|-- Size = "<<display_.size<<"\n";
}
