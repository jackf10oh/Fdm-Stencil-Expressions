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

#include<LinOps/All.hpp> 
#include<OutsideSteps/All.hpp> 
#include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 
#include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(7); 

  // defining Domain Mesh --------------------------------------
  auto r = 3.14159; // pi 
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = LinOps::make_mesh(0.0, 2.0, 51); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::Discretization1D my_vals;
  auto sin_lam = [](double x){return 2.0 * std::sin(x); }; 
  my_vals.set_init(my_mesh, sin_lam); 
  print_vec(my_vals.values(), "ICs"); 

  // LHS time derivs ----------------------------------------------------------------
  auto time_expr = TExprs::NthTimeDeriv(1); 

  // building RHS expression -----------------------------------------------------
  using D = LinOps::NthDerivOp;
  auto space_expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  auto left = OSteps::DirichletBC(0.0); 
  auto right = OSteps::NeumannBC(0.0); 

  auto bcs = OSteps::BCPair(left,right); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = std::tie(bcs), 
    .ICs = std::vector<Eigen::VectorXd>(1, my_vals.values()), 
  }; 

  // TExprs::GenSolver s(time_expr, space_expr, TExprs::PrintWrite{}); 
  // auto v = s.Calculate(args); 

  TExprs::GenInterp interp(time_expr, space_expr, args); 
  interp.FillVals(); 
  print_mat(interp.StoredData(), "Solutions through time");


  std::cout << time_mesh->size() << std::endl; 
  std::cout << interp.StoredData().size() << std::endl; 
};