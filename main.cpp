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

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto m1 = LinOps::make_mesh(0.0,5.0, 6); 
  auto m2 = LinOps::make_meshXD(0.0,3.0,3,2); 

  cout << LinOps::IOp(m1).GetMat() << endl; 

  cout << LinOps::IOp(m2).GetMat() << endl; 

  cout << LinOps::DirectionalRandOp(m2,0).GetMat() << endl; 
  cout << LinOps::DirectionalRandOp(m2,1).GetMat() << endl; 

  auto expr_xd = 20.0 * LinOps::DirectionalRandOp(); // should not be 1D
  auto expr_1d = 20.0 * LinOps::RandLinOp(); // should not be XD 
  cout << "expr_xd is 1D? " << LinOps::traits::is_1dim_linop_crtp<decltype(expr_xd)>::value << endl; 
  cout << "expr_xd is XD? " << LinOps::traits::is_xdim_linop_crtp<decltype(expr_xd)>::value << endl; 
  cout << "expr_1d is XD? " << LinOps::traits::is_xdim_linop_crtp<decltype(expr_1d)>::value << endl; 

  auto expr = LinOps::DirectionalRandOp(0) + 20.0 * LinOps::IOp(); // works .... 
  // auto expr = LinOps::RandLinOp() + 20.0 * LinOps::DirectionalRandOp(); // will not compile: static_assert() fails .... 
  expr.set_mesh(m2); 
  cout << expr.GetMat() << endl; 

  auto m3 = LinOps::make_mesh(0.0,10.0,11); 
  std::cout << LinOps::NthDerivOp(m3,2).GetMat() << std::endl; 

  auto m4 = LinOps::make_meshXD(0.0,3.0,4, 3); 
  auto D = LinOps::DirectionalNthDerivOp(m4,1,0); 
  cout << D.GetMat() << endl; 


  cout << "------------" << endl; 
  LinOps::IOp I; 
  I.set_mesh(m1); 
  cout << "I.get_mesh1d == nullptr" << (I.get_mesh1d() == m1) << endl; 

  auto sum_i = I + LinOps::IOp(); 
  std::cout << "rhs mesh == nullptr " << (sum_i.Rhs().get_mesh1d()==nullptr) << std::endl; 

  sum_i.set_mesh(m1); 
  std::cout << "post set_mesh. rhs mesh == nullptr " << (sum_i.Rhs().get_mesh1d()==nullptr) << std::endl; 
  cout << sum_i.Lhs().GetMat() << endl; 
  cout << sum_i.Rhs().GetMat() << endl; 
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