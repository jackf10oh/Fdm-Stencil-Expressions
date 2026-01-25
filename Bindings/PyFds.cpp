// PyFds.cpp
//
// pybind11 bindsing for LinOps mesh + discretization
//
// JAF 1/25/2026 

#include<pybind11/pybind11.h> 
#include<pybind11/stl.h>
#include<pybind11/functional.h> 

#include "../LinOps/All.hpp"
#include "../LinOpsXD/All.hpp"

namespace py = pybind11; 

PYBIND11_MODULE(PyFds, m)
{
  py::class_<LinOps::Mesh1D, std::shared_ptr<LinOps::Mesh1D>>(m, "Mesh1D")
    .def(py::init<double,double,std::size_t>(), 
           py::arg("start")=0.0, 
           py::arg("stop")=1.0, 
           py::arg("nsteps")=11
    )
    .def(py::init<std::vector<double>>(), 
           py::arg("vals")
    );    

  py::class_<LinOps::MeshXD, std::shared_ptr<LinOps::MeshXD>>(m,"MeshXD")
    .def(py::init<double,double,std::size_t,std::size_t>(), 
           py::arg("start")=0.0, 
           py::arg("stop")=1.0, 
           py::arg("nsteps")=11, 
           py::arg("ndims")=2
    )
    .def(py::init<std::vector<std::pair<double,double>>, std::vector<std::size_t>>(), 
           py::arg("endpoints_list"), 
           py::arg("nsteps_list")
    )
    .def(py::init<std::shared_ptr<const LinOps::Mesh1D>>(),
      py::arg("mesh1d")
    )
    .def(py::init<std::vector<std::shared_ptr<const LinOps::Mesh1D>>>(),
           py::arg("mesh1d_list")
    );
}

