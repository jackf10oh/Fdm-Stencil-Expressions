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

#include "Utilities/PrintVec.hpp"

#include "FDStencils/All.hpp"

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 101;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 1.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::Discretization1D my_vals;
  // Bump centered at r/2. Zero at 0.0 and r. 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);};  
  my_vals.set_init(my_mesh, func); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  auto left_bc = Fds::make_dirichlet(0.0); 
  auto right_bc = Fds::make_dirichlet(0.0); 

  // Solving --------------------------------------------------------------------- 
  Fds::SolverArgs1D args 
  {
    .domain_mesh_ptr   = my_mesh,  
    .time_mesh_ptr     = time_mesh, 
    .bcs_pair          = { left_bc, right_bc},  
    .ICs               = my_vals, 
    .time_dep_flag     = false 
  };

  Fds::Solver1D s(expr); 
  // auto result = s.Calculate( args );
  auto result = s.CalculateImp( args );

  // Printing ---------------------------------------------------------------- 
  print_vec(my_vals, "ICs");  
  print_vec(result, "solution"); 

  // expr.set_mesh(my_mesh); 
  // std::cout << expr.GetMat() << std::endl; 
};

// explicit step 
// [U(t+1) - U(t)] / dt = D * U(t) 
// U(t+1) = U(t) + dt * D * U(t) 


// implicit step 
// [U(t+1) - U(t)] / dt = D * U(t+1) 
// U(t+1) - U(t) = dt * D * U(t+1) 
// I * U(t+1) - dt * D * U(n+1) = U(t) 
// ( I - dt * D) * U(t+1) = U(t) 
// U(t+1) = (I - dt * D).inv(U(t)) 

