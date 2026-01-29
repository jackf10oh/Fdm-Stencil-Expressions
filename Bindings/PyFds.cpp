// PyFds.cpp
//
// pybind11 bindsing for LinOps mesh + discretization
//
// JAF 1/25/2026 

#include<pybind11/pybind11.h> 
#include<pybind11/stl.h>
#include<pybind11/functional.h> 
#include<pybind11/eigen.h>

#include<LinOps/All.hpp>
#include<LinOpsXD/All.hpp>

namespace py = pybind11; 

PYBIND11_MODULE(PyFds, m)
{
  // Mesh1D ================================================================================ 
  py::class_<LinOps::Mesh1D, std::shared_ptr<LinOps::Mesh1D>>(m, "Mesh1D")
    .def(py::init<double,double,std::size_t>(), 
           py::arg("start")=0.0, 
           py::arg("stop")=1.0, 
           py::arg("nsteps")=11
    )
    .def(py::init<std::vector<double>>(), 
           py::arg("vals")
    )
    .def("size",&LinOps::Mesh1D::size, "size of current mesh");    

  // MeshXD ================================================================================ 
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
    )
    .def("GetMesh", 
           &LinOps::MeshXD::GetMesh,
           py::arg("dim")=1)
    .def("dims", 
           &LinOps::MeshXD::dims)
    .def("dim_size", 
           &LinOps::MeshXD::dim_size,
           py::arg("dim")=1)
    .def("sizes_product", 
           &LinOps::MeshXD::sizes_product)
    .def("sizes_middle_product", 
           &LinOps::MeshXD::sizes_middle_product,
           py::arg("start"),
           py::arg("stop"));
  
  // Discretization1D ================================================================================ 
  py::class_<LinOps::Discretization1D, std::unique_ptr<LinOps::Discretization1D>>(m,"Discretization1D")
    .def(py::init<std::size_t>(), 
           py::arg("size")=0
    )
    .def(py::init<std::shared_ptr<const LinOps::Mesh1D>>(), 
           py::arg("mesh")
    )
    .def(py::init<Eigen::VectorXd&&>(),
           py::arg("arr")
    )
    .def(py::init<const Eigen::VectorXd&>(),
           py::arg("arr")
    )
    .def("values",
            [](const LinOps::Discretization1D& self){return self.values();}, 
            py::return_value_policy::reference_internal
    )
    .def("size", 
            &LinOps::Discretization1D::size, 
            "size of current discretization"
    )
    .def("resize",
            &LinOps::Discretization1D::resize, 
            py::arg("mesh"),"resize Discretization1D to fit on a Mesh1D"
    )
    .def("set_init",
            [](LinOps::Discretization1D& self, double c){self.set_init(c);}, 
            py::arg("constant")
    )
    .def("set_init",
            [](LinOps::Discretization1D& self, const LinOps::Mesh1D_SPtr_t& m, double c){self.set_init(m,c);}, 
            py::arg("mesh"),
            py::arg("constant")
    )
    .def("set_init",
      [](LinOps::Discretization1D& self, std::function<double(double)> f){self.set_init(f); },
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::Discretization1D& self, const LinOps::Mesh1D_SPtr_t& m, std::function<double(double)> f){self.set_init(m,f); },
      py::arg("mesh"),
      py::arg("func")
    ); 

  // DiscretizationXD ================================================================================ 
  py::class_<LinOps::DiscretizationXD, std::unique_ptr<LinOps::DiscretizationXD>>(m,"DiscretizationXD")
    .def(py::init<std::size_t>(), 
           py::arg("size")=0
    )
    .def(py::init<std::shared_ptr<const LinOps::MeshXD>>(), 
           py::arg("mesh")
    )
    .def(py::init<Eigen::VectorXd&&>(),
           py::arg("arr")
    )
    .def(py::init<const Eigen::VectorXd&>(),
           py::arg("arr")
    )
    .def("values",
            [](const LinOps::DiscretizationXD& self){return self.values();}, 
            py::return_value_policy::reference_internal
    )
    .def("sizes_product", 
            &LinOps::DiscretizationXD::sizes_product, 
            "product of each dimensions size"
    )
    .def("dims", 
            &LinOps::DiscretizationXD::dims, 
            "number of dimensions in DiscretizationXD"
    )
    .def("dim_size", 
            &LinOps::DiscretizationXD::dim_size,
            py::arg("dim")=0, 
            "size of specific dimension in DiscretizationXD"
    )
    .def("sizes_middle_product", 
            &LinOps::DiscretizationXD::sizes_middle_product,
            py::arg("start"), 
            py::arg("stop"), 
            "product of sizes of dimensions in [start,stop)" 
    )
    .def("resize",
            &LinOps::DiscretizationXD::resize, 
            py::arg("mesh"),"resize DiscretizationXD to fit on a MeshXD"
    )
    .def("set_init",
            [](LinOps::DiscretizationXD& self, double c){self.set_init(c);}, 
            py::arg("constant")
    )
    .def("set_init",
            [](LinOps::DiscretizationXD& self, const LinOps::MeshXD_SPtr_t& m, double c){self.set_init(m,c);}, 
            py::arg("mesh"),
            py::arg("constant")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, std::function<double(double)> f){self.set_init(f); },
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, const LinOps::MeshXD_SPtr_t& m, std::function<double(double)> f){self.set_init(m,f); },
      py::arg("mesh"),
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, std::function<double(double,double)> f){self.set_init(f); },
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, const LinOps::MeshXD_SPtr_t& m, std::function<double(double,double)> f){self.set_init(m,f); },
      py::arg("mesh"),
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, std::function<double(double,double,double)> f){self.set_init(f); },
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, const LinOps::MeshXD_SPtr_t& m, std::function<double(double,double,double)> f){self.set_init(m,f); },
      py::arg("mesh"),
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, std::function<double(double,double,double,double)> f){self.set_init(f); },
      py::arg("func")
    )
    .def("set_init",
      [](LinOps::DiscretizationXD& self, const LinOps::MeshXD_SPtr_t& m, std::function<double(double,double,double,double)> f){self.set_init(m,f); },
      py::arg("mesh"),
      py::arg("func")
    ); 
}

