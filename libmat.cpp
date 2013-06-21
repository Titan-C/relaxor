#include "material.h"
#include "material.cpp"
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

BOOST_PYTHON_MODULE(libmat)
{
    using namespace boost::python;
    class_<std::vector<double> >("double_vector")
        .def(vector_indexing_suite<std::vector<double> >())
    ;
    class_<Material>("Material", init<unsigned int,double, std::string, optional<bool, bool, bool> >())
      .def("state", &Material::state)
      .def("inicio", &Material::set_interaction_dipole_config)
      .def("set_rho", &Material::set_rho)
      .def("set_ExpID", &Material::set_ExpID)
      .def("__repr__", &Material::repr)
    ;
}
