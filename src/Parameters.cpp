
#include <iostream>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>
#include "Parameters.h"


using namespace boost;

Parameters::General::
General():
  events(0),
  debug(false)
{
}

const std::map<std::string, Parameters::Geometry::Type>
Parameters::Geometry::type_map_ = {
  {"Hexagons", Parameters::Geometry::Type::Hexagons},
  {"Triangles", Parameters::Geometry::Type::Triangles},
  {"External", Parameters::Geometry::Type::External}
};

Parameters::Geometry::
Geometry():
  type(Type::Undefined),
  layer(-1),
  layers_z(0),
  side_length(0),
  cells_x(0),
  cells_y(0),
  file("")
{
}

Parameters::Generation::
Generation():
  fluctuation(false),
  energy(0.),
  number_of_hits_per_gev(0),
  alpha(0.),
  mip_energy(0.),
  position(0.),
  layers_energy(0),
  incident_x(0.),
  incident_y(0.),
  incident_eta(0.),
  r0layer15(0.),
  shower_transverse_parameters(0)
{
}

Parameters::Display::
Display():
  events(0),
  size(0)
{
}

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
  // Read parameters common to all geometry types
  std::string type = python::extract<std::string>(dict["geometry_type"]);
  // FIXME: not needed, map::at will throw an exception if not defined
  if(Geometry::type_map_.find(type)==Geometry::type_map_.end()) throw std::string("Unknown type of geometry");
  geometry_.type = Geometry::type_map_.at(type);
  geometry_.layer = python::extract<int>(dict["geometry_layer"]);
  // 
  python::list layers_z = python::extract<python::list>(dict["geometry_layers_z"]);
  python::stl_input_iterator<double> begin(layers_z), end;
  geometry_.layers_z = std::vector<double>(begin,end);
  // Read parameters for internal geometries (infinite hexagons or triangles)
  if(geometry_.type!=Geometry::Type::External)
  {
    geometry_.side_length = python::extract<double>(dict["geometry_side_length"]);
    geometry_.cells_x = python::extract<int>(dict["geometry_cells_x"]);
    geometry_.cells_y = python::extract<int>(dict["geometry_cells_y"]);
  }
  // Read parameters for external json geometries
  else
  {
    geometry_.file = python::extract<std::string>(dict["geometry_file"]);
  }

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
  python::stl_input_iterator<double> begin_energy(layers_energy), end_energy;
  generation_.layers_energy = std::vector<double>(begin_energy,end_energy);
  generation_.incident_x = python::extract<double>(dict["generation_incident_x"]);
  generation_.incident_y = python::extract<double>(dict["generation_incident_y"]);
  generation_.incident_eta = python::extract<double>(dict["generation_incident_eta"]);
  generation_.r0layer15 = python::extract<double>(dict["generation_r0layer15"]);
  //
  python::list shower_transverse_parameters = python::extract<python::list>(dict["generation_shower_transverse_parameters"]);
  python::stl_input_iterator<double> begin_profile(shower_transverse_parameters), end_profile;
  generation_.shower_transverse_parameters = std::vector<double>(begin_profile,end_profile);
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
  std::cout<<"|-- Type = "<<static_cast<std::underlying_type<Geometry::Type>::type>(geometry_.type)<<"\n";
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
  std::cout<<"|-- File = "<<geometry_.file<<"\n";
  std::cout<<"|- Generation\n";
  std::cout<<"|-- Energy = "<<generation_.energy<<"\n";
  std::cout<<"|-- Fluctuation = "<<generation_.fluctuation<<"\n";
  std::cout<<"|-- Hits/GeV = "<<generation_.number_of_hits_per_gev<<"\n";
  std::cout<<"|-- Alpha = "<<generation_.alpha<<"\n";
  std::cout<<"|-- Mip energy = "<<generation_.mip_energy<<"\n";
  std::cout<<"|-- Position = ("<<generation_.incident_x<<" "<<generation_.incident_y<<")\n";
  std::cout<<"|-- Angle (eta) = "<<generation_.incident_eta<<"\n";
  std::cout<<"|-- Layers energy = [";
  for(const auto& e : generation_.layers_energy)
  {
    std::cout<<e<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|-- R_0(15) = "<<generation_.r0layer15<<"\n";
  std::cout<<"|-- Transverse parameters = [";
  for(const auto& p : generation_.shower_transverse_parameters)
  {
    std::cout<<p<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|- Display\n";
  std::cout<<"|-- Events = "<<display_.events<<"\n";
  std::cout<<"|-- Size = "<<display_.size<<"\n";
}
