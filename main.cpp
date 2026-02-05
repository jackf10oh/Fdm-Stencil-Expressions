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

#include<LinOps/Operators/DiffOps/experimental_NthDerivOp.hpp> 
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
  std::cout << std::setprecision(4); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0;
  int n_gridpoints = 41;

  // mesh in space [0, 10] 
  auto my_mesh = LinOps::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time [0, 2]
  auto time_mesh = LinOps::make_mesh(0.0, 2.0, 11); 

  // Initial values 
  BumpFunc f{.L = 4, .R = 6, .c = 5, .h = 3}; 
  Eigen::VectorXd v = LinOps::Discretization1D().set_init(my_mesh, f).values(); 
  std::vector<Eigen::VectorXd> ics = {v}; 

  // defining LHS + RHS equations -------------------
  using D = LinOps::NthDerivOp;
  auto Rhs = 0.2 * D(2) + 0.5 * D(1); 
  auto Lhs = TExprs::NthTimeDeriv(1); 

  // defining Boundary Conditions 
  auto bcs = OSteps::BCPair( OSteps::DirichletBC(), OSteps::DirichletBC()); 

  // Solving -------------------------
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh, 
    .time_mesh_ptr = time_mesh, 
    .ICs = ics 
  }; 

  TExprs::GenSolver solver(Lhs, Rhs, std::tie(bcs), TExprs::PrintWrite{}); 
  
  solver.Calculate(args); 

};