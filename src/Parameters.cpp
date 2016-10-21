
#include <iostream>
#ifdef STANDALONE
#include "Parameters.h"
#else
#include "HGCalSimulation/FastShower/interface/Parameters.h"
#endif


using namespace boost;

Parameters::General::
General():
  events(0),
  debug(false),
  output_file("")
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
  cell_side(0),
  offset(0),
  cells_nx(0),
  cells_ny(0),
  file("")
{
}

Parameters::Shower::
Shower():
  moliere_radius(0.),
  radiation_length(0.),
  critical_energy(0.),
  transverse_parameters(),
  longitudinal_parameters(),
  layers_energy(0),
  alpha(0.)
{
}

Parameters::Generation::
Generation():
  fluctuation(false),
  energy(0.),
  number_of_hits_per_gev(0),
  mip_energy(0.),
  sampling(0.),
  noise(false),
  noise_sigma(0.),
  incident_x(0.),
  incident_y(0.),
  incident_eta(0.)
{
}

Parameters::Display::
Display():
  events(0),
  size(0),
  offset_x(0.),
  offset_y(0.)
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
  fillShower(main_namespace);
  fillGeneration(main_namespace);
  fillDisplay(main_namespace);
}

void 
Parameters::
fillGeneral(const python::dict& dict)
{
  general_.events = python::extract<int>(dict["events"]);
  general_.debug = python::extract<bool>(dict["debug"]);
  general_.output_file = python::extract<std::string>(dict["output_file"]);
}

void 
Parameters::
fillGeometry(const python::dict& dict)
{
  // Read parameters common to all geometry types
  std::string type = python::extract<std::string>(dict["geometry_type"]);
  // FIXME: not needed, map::at will throw an exception if not defined
  if(Geometry::type_map_.find(type)==Geometry::type_map_.end()) throw std::string("Unknown type of geometry");
  geometry_.type = Geometry::type_map_.at(type);
  geometry_.layer = python::extract<int>(dict["geometry_layer"]);
  // 
  python::list layers_z = python::extract<python::list>(dict["geometry_layers_z"]);
  geometry_.layers_z = toStdVector<double>(layers_z);
  // Read parameters for internal geometries (infinite hexagons or triangles)
  if(geometry_.type!=Geometry::Type::External)
  {
    geometry_.cell_side = python::extract<double>(dict["geometry_cell_side"]);
    geometry_.offset = python::extract<int>(dict["geometry_offset"]);
    geometry_.cells_nx = python::extract<int>(dict["geometry_cells_nx"]);
    geometry_.cells_ny = python::extract<int>(dict["geometry_cells_ny"]);
  }
  // Read parameters for external json geometries
  else
  {
    geometry_.file = python::extract<std::string>(dict["geometry_file"]);
  }

}

void 
Parameters::
fillShower(const python::dict& dict)
{
  shower_.moliere_radius = python::extract<double>(dict["shower_moliere_radius"]);
  shower_.radiation_length = python::extract<double>(dict["shower_radiation_length"]);
  shower_.critical_energy = python::extract<double>(dict["shower_critical_energy"]);
  python::dict transverse_parameters = python::extract<python::dict>(dict["shower_transverse_parameters"]);
  shower_.transverse_parameters = toStdMap<std::string,double>(transverse_parameters);
  python::dict longitudinal_parameters = python::extract<python::dict>(dict["shower_longitudinal_parameters"]);
  shower_.longitudinal_parameters = toStdMap<std::string,double>(longitudinal_parameters);
  python::list layers_energy = python::extract<python::list>(dict["shower_layers_energy"]);
  shower_.layers_energy = toStdVector<double>(layers_energy);
  shower_.alpha = python::extract<double>(dict["shower_alpha"]);
}

void 
Parameters::
fillGeneration(const python::dict& dict)
{
  generation_.fluctuation = python::extract<bool>(dict["generation_fluctuation"]);
  generation_.energy = python::extract<double>(dict["generation_energy"]);
  generation_.number_of_hits_per_gev = python::extract<int>(dict["generation_number_of_hits_per_gev"]);
  generation_.mip_energy = python::extract<double>(dict["generation_mip_energy"]);
  generation_.sampling = python::extract<double>(dict["generation_sampling"]);
  generation_.noise = python::extract<bool>(dict["generation_noise"]);
  generation_.noise_sigma = python::extract<double>(dict["generation_noise_sigma"]);
  generation_.incident_x = python::extract<double>(dict["generation_incident_x"]);
  generation_.incident_y = python::extract<double>(dict["generation_incident_y"]);
  generation_.incident_eta = python::extract<double>(dict["generation_incident_eta"]);
}

void 
Parameters::
fillDisplay(const python::dict& dict)
{
  display_.events = python::extract<int>(dict["display_events"]);
  display_.size = python::extract<int>(dict["display_size"]);
  display_.offset_x = python::extract<double>(dict["display_offset_x"]);
  display_.offset_y = python::extract<double>(dict["display_offset_y"]);
}



void
Parameters::
print() const
{
  std::cout<<"|Configuration parameters\n";
  std::cout<<"|- General\n";
  std::cout<<"|-- Events = "<<general_.events<<"\n";
  std::cout<<"|-- Debug = "<<general_.debug<<"\n";
  std::cout<<"|-- Output file = "<<general_.output_file<<"\n";
  std::cout<<"|- Geometry\n";
  std::cout<<"|-- Type = "<<static_cast<std::underlying_type<Geometry::Type>::type>(geometry_.type)<<"\n";
  std::cout<<"|-- SideLength = "<<geometry_.cell_side<<"\n";
  std::cout<<"|-- Cells_NX = "<<geometry_.cells_nx<<"\n";
  std::cout<<"|-- Cells_NY = "<<geometry_.cells_ny<<"\n";
  std::cout<<"|-- Offset = "<<geometry_.offset<<"\n";
  std::cout<<"|-- Layer = "<<geometry_.layer<<"\n";
  std::cout<<"|-- Layers z = [";
  for(const auto& z : geometry_.layers_z)
  {
    std::cout<<z<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|-- File = "<<geometry_.file<<"\n";
  std::cout<<"|- Shower\n";
  std::cout<<"|-- Moliere radius = "<<shower_.moliere_radius<<"\n";
  std::cout<<"|-- Radiation length = "<<shower_.radiation_length<<"\n";
  std::cout<<"|-- Critical energy = "<<shower_.critical_energy<<"\n";
  std::cout<<"|-- Transverse parameters = [";
  for(const auto& name_value : shower_.transverse_parameters)
  {
    std::cout<<name_value.first<<"("<<name_value.second<<") ";
  }
  std::cout<<"]\n";
  std::cout<<"|-- Longitudinal parameters = [";
  for(const auto& name_value : shower_.longitudinal_parameters)
  {
    std::cout<<name_value.first<<"("<<name_value.second<<") ";
  }
  std::cout<<"]\n";
  std::cout<<"|-- Layers energy = [";
  for(const auto& e : shower_.layers_energy)
  {
    std::cout<<e<<" ";
  }
  std::cout<<"]\n";
  std::cout<<"|-- Alpha = "<<shower_.alpha<<"\n";
  std::cout<<"|- Generation\n";
  std::cout<<"|-- Energy = "<<generation_.energy<<"\n";
  std::cout<<"|-- Fluctuation = "<<generation_.fluctuation<<"\n";
  std::cout<<"|-- Hits/GeV = "<<generation_.number_of_hits_per_gev<<"\n";
  std::cout<<"|-- Mip energy = "<<generation_.mip_energy<<"\n";
  std::cout<<"|-- Sampling = "<<generation_.sampling<<"\n";
  std::cout<<"|-- Noise = "<<generation_.noise<<"\n";
  std::cout<<"|-- Noise sigma = "<<generation_.noise_sigma<<"\n";
  std::cout<<"|-- Position = ("<<generation_.incident_x<<" "<<generation_.incident_y<<")\n";
  std::cout<<"|-- Angle (eta) = "<<generation_.incident_eta<<"\n";
  std::cout<<"|- Display\n";
  std::cout<<"|-- Events = "<<display_.events<<"\n";
  std::cout<<"|-- Size = "<<display_.size<<"\n";
  std::cout<<"|-- Offset = "<<display_.offset_x<<" "<<display_.offset_y<<"\n";
}
