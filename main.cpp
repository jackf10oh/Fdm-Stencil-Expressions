// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

// #include<FDStencils/All.hpp> // must be first for plugin macro.
#include<LinOps/All.hpp> 
// #include<FDStencilsXD/All.hpp> // likewise ...
// #include<LinOpsXD/All.hpp>
// #include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

void foobar(const std::shared_ptr<const LinOps::Mesh1D>& m){
  cout << "use count: " << m.use_count() << endl; 
}

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto m = LinOps::make_mesh(0.0,4.0,5); 

  LinOps::TimeDepCoeff c = [](double t, double x){return t + x; }; 
  // LinOps::TimeDepCoeff c = [](double t, double x, double y){return t + x; }; // 2 dims fails on set_mesh(1D)
  LinOps::TimeDepCoeff c02 = [](double t){return t*t; }; 

  cout << c.Time() << endl; 
  c.set_mesh(m); 
  c.SetTime(1.0); 
  cout << c.Time() << endl; 
  cout << c.GetMat() << endl; 

  cout << c02.Time() << endl; 
  c02.SetTime(2.0); 
  cout << c02.Time() << endl; 
  cout << c02.GetMat() << endl; // SCALAR! 

  auto m2d = LinOps::make_meshXD(0.0,3.0,4,2); 
  c.set_mesh(m2d); 
  c.SetTime(4.0); 
  cout << c.GetMat() << endl; 
};

// #include "Bindings/PyPdeBase1D.hpp"

// struct HeatPDE_impl : public PDE_Base_Impl<HeatPDE_impl>
// { 
//   double diffusion = 1.0; 
//   double convection = 0.0; 
//   double reaction = 0.0; 

//   LinOps::IOp U = LinOps::IOp(); 
//   Fds::NthDerivOp Ux = Fds::NthDerivOp(1); 
//   Fds::NthDerivOp Uxx = Fds::NthDerivOp(2); 

//   // rhs expr in space 
//   decltype(diffusion*Uxx + convection*Ux + reaction*U) rhs_expr = diffusion*Uxx + convection*Ux + reaction*U;  
  
//   // lhs expr in time 
//   TExprs::NthTimeDeriv Ut = TExprs::NthTimeDeriv(1); 

//   auto& GetLhs(){ return Ut; } 
//   auto& GetRhs(){ return rhs_expr; } 
// }; 

// using HeatPDE = Concrete_PDE_1D<HeatPDE_impl>; 

//   HeatPDE pde; 

//   std::cout << "Enter Diffusion:"; 
//   std::cin >> pde.diffusion; 
//   std::cout << "Enter Convection:"; 
//   std::cin >> pde.convection; 
//   std::cout << "Enter Reaction:"; 
//   std::cin >> pde.reaction; 

//   auto my_mesh = LinOps::make_mesh(0.0, 10.0, 17);
  
//   // domain 
//   pde.Args().domain_mesh_ptr = my_mesh; 
//   // time steps 
//   pde.Args().time_mesh_ptr = LinOps::make_mesh(0.0,5.0,21); 
//   // boundary conditions 
//   pde.Args().bcs = std::make_shared<Fds::BCPair>(Fds::make_dirichlet(0.0), Fds::make_dirichlet(0.0)); 
//   // initial conditions 
//   BumpFunc f{.L=0.0, .R=10.0, .c = 5.0, .h = 3.0}; 
//   pde.Args().ICs = std::vector<Eigen::VectorXd>{ std::move(LinOps::Discretization1D().set_init(my_mesh,f).values()) };
//   // time dep flag 
//   pde.Args().time_dep_flag = false;

//   // calculate 
//   pde.Reset(); 
//   pde.FillVals(); 

//   // print 
//   print_mat(pde.StoredData(), "Solution"); 