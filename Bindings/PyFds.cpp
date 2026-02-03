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
#include "Wave2D.hpp" 

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
    .def("at",[](const LinOps::Mesh1D& self, std::size_t i){return self.at(i);}, py::arg("idx"), "values of mesh at an index")    
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
           [](LinOps::MeshXD& self, std::size_t i){ return self.GetMeshAt(i); },
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
  // Mesh1D ================================================================================ 
  py::class_<Wave2D>(m, "Wave2D")
    .def(py::init<>())
    .def("SetDomain", [](Wave2D& self, LinOps::MeshXD_SPtr_t m) -> Wave2D& 
      {auto a = self.Args(); a.domain_mesh_ptr=m; self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new domain mesh inside of solver") 
    .def("SetTime", [](Wave2D& self, LinOps::Mesh1D_SPtr_t m) -> Wave2D& 
      {auto a = self.Args(); a.time_mesh_ptr=m; self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new time mesh inside of solver") 
    .def("SetIC", [](Wave2D& self, std::vector<Eigen::VectorXd> v) -> Wave2D& 
      {auto a = self.Args(); a.ICs=std::move(v); self.SetArgs(std::move(a)); return self; }, py::arg("mesh"), "Sets a new initial condition inside of solver") 
    .def("SetHeight",[](Wave2D& self, double h) -> Wave2D& 
      { self.set_bump_height(h); self.Reset(); return self;}, py::arg("h")=1.0, "Sets new height for oscilation force term at origin")
    .def("SetDamping", [](Wave2D& self, double d) -> Wave2D& 
      {self.set_damping(d); self.Reset(); return self;}, py::arg("d")=2.0, "Sets a new damping rate at boundaries of the domain")
    .def("StoredData", &Wave2D::StoredData, "returns a vector of solutions. solutions are flattened in dimensional order")
    .def("Compute", &Wave2D::FillVals, "Computes solution at each entry in time. available in StoredData")
    .def("SolAt", [](Wave2D& self, double t, double x, double y){ return self.SolAt(t,x,y); }, py::arg("t"), py::arg("x"), py::arg("y"), "Returns value of solution at time t at coords (x,y)");    
}

