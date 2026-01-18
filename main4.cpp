// main4.cpp 
//
// Testing out GenSolver on PDEs 
//
// JAF 1/17/2026 

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include "Utilities/PrintVec.hpp"
#include "Utilities/BumpFunc.hpp"

#include "FDStencils/All.hpp" // must be first for plugin macro. 
#include "LinOps/All.hpp" 
// #include "FDStencilsXD/All.hpp" // likewise ...
// #include "LinOpsXD/All.hpp"

#include "LhsExpressions/All.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 101;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 5.0, 8); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::Discretization1D my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  BumpFunc f{.L=0.0, .R=r, .c=r/2, .h=1.0, .focus=20.0};
  my_vals.set_init(my_mesh, f); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto expr = 0.2 * D(2) - 0.5 * D(1); 
  // auto expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<BCPair>(left,right); 

  // LHS time derivs ----------------------------------------------------------------
  // auto lhs_expr = NthTimeDeriv(1); 
  auto lhs_expr = NthTimeDeriv(1); 

  // Solving --------------------------------------------------------------------- 
  GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(1, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  GenSolver s(lhs_expr, expr); 
  // auto v = s.Calculate(args); 
  auto v = s.CalculateImp(args); 

  // Printing ---------------------------------------------------------------- 
  print_vec(my_vals, "ICs"); 
  print_vec(v,"Sol"); 
};