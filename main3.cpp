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
#include "Utilities/BumpFunc.hpp"

#include "FDStencils/All.hpp" // must be first for plugin macro.
#include "FDStencils/Interp1D.hpp" // Interp1D class 
#include "LinOps/All.hpp" 
#include "FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/All.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(4); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 101;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 4.0, 401); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::Discretization1D my_vals;
  // Bump centered at r/2. Zero at 0.0 and r. 
  BumpFunc f{.L=0.0, .R=r, .c=r/2, .h=2}; 
  my_vals.set_init(my_mesh, f); 
  print_vec(my_vals, "ICs");

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_dirichlet(0.0); 
  auto bcs = std::make_shared<Fds::BCPair>(left, right); 

  // Solving --------------------------------------------------------------------- 
  Fds::SolverArgs1D args 
  {
    .domain_mesh_ptr   = my_mesh,  
    .time_mesh_ptr     = time_mesh, 
    .bcs               = bcs,  
    .ICs               = my_vals, 
    .time_dep_flag     = false 
  };

  Interp1D calc(expr, args); 
  // Printing ---------------------------------------------------------------- 
  while(true){
    double t; 
    cout << "enter t: "; 
    std::cin >> t; 

    auto x_it = my_mesh->cbegin(); 
    auto end = std::prev(my_mesh->cend()); 
    
    std::cout << "["; 
    for(; x_it!=end; x_it++)
    {
      std::cout << calc.Value(t,*x_it) << ", "; 
    }; 
    std::cout << calc.Value(t,*x_it) << "]" << endl; 
  }
};
 