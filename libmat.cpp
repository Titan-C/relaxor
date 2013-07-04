#include "material.h"
#include "material.cpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#define _2pi 8*atan(1)

double simpson_int(const std::vector<double>& f_array, const std::vector<double>& weight){
  unsigned int length=weight.size();

  double Integral = f_array[0]*weight[0];

  for(unsigned int i=1; i<length; i+=2)
    Integral+=4*f_array[i]*weight[i];

  for(unsigned int i=2; i<length-1; i+=2)
    Integral+=2*f_array[i]*weight[i];

  length--;
  Integral+=f_array[length]*weight[length];

  return Integral/3;
}

std::vector<double> wavearray(double amplitude, unsigned int tau, unsigned int length, double phase){
  std::vector<double> wave;
  wave.resize(length);
  unsigned int wavetop = (tau>=length)? length : tau;
    for(unsigned int i=0; i<wavetop; i++)
      wave[i]=amplitude*cos(_2pi*(i-phase)/tau);

  unsigned int periods = length/tau;
  for(unsigned int i=1; i<periods;i++){
    for(unsigned int j = 0; j<tau ;j++)
      wave[i*tau+j]=wave[j];
  }
  return wave;
}

BOOST_PYTHON_MODULE(libmat)
{
    using namespace boost::python;
    
    def("simpson_int", simpson_int);
    def("wavearray", wavearray);
    
    class_<std::vector<double> >("double_vector")
        .def(vector_indexing_suite<std::vector<double> >())
    ;
    class_<Material>("Material", init<unsigned int,double, std::string, optional<bool, bool, bool> >())
      .def("state", &Material::state)
      .def("oven", &Material::oven)
      .def("inicio", &Material::set_interaction_dipole_config)
      .def("set_rho", &Material::set_rho)
      .def("set_ExpID", &Material::set_ExpID)
      .def("__repr__", &Material::repr)
    ;
}
