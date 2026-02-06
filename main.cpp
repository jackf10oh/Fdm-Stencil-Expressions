// main.cpp
//
//
//
// JAF 12/8/2025

// #include "DiffOps/DiffOpTraits.hpp"
// #include "DiffOps/NthDerivOp.hpp"

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

#include<DiffOps/DiffOpTraits.hpp>
#include<DiffOps/DiffOpBase.hpp>
#include<DiffOps/NthDerivOp.hpp>

using std::cout, std::endl, std::cin;

template<std::size_t N>
using D = DiffOps::NthDerivOp<N>;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  double diffusion, convection; 

  cout << "Enter Diffusion: ";
  cin >> diffusion; 
  cout << "Enter Convection: "; 
  cin >> convection; 

  // Defining LHS Expression in time + RHS expression in space 
  auto Lhs = TExprs::NthTimeDeriv(1); 
  auto Rhs = diffusion * D<2>{} + convection * D<1>{}; 

  // defining boundary conditions 
  auto bcs = OSteps::BCPair(OSteps::DirichletBC(0.0), OSteps::DirichletBC(0.0)); 

  // defining Args Mesh --------------------------------------

  // mesh in space [0, 10]
  auto r = 10.0;
  int n_gridpoints = 21; 
  auto domain = LinOps::make_mesh(0.0,r,n_gridpoints); 

  BumpFunc f{.L = 4, .R = 6, .c = 5, .h = 3}; 
  Eigen::VectorXd v = LinOps::Discretization1D().set_init(domain, f).values(); 
  
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = domain, 
    .time_mesh_ptr = LinOps::make_mesh(0.0,4.0, 30),
    .ICs = {std::move(v)}
  }; 

  // Solving 
  TExprs::GenSolver solver(Lhs, Rhs, std::tie(bcs), TExprs::PrintWrite{}); 

  solver.Calculate(args); 



}