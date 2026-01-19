// main4.cpp 
//
// Testing out GenSolver on 2D PDEs 
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
#include "FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/All.hpp"
#include "TExprs/All.hpp"

using std::cout, std::endl;

// using namespace Fds; 
// using namespace LinOps; 
// using namespace TExprs

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = LinOps::make_meshXD(0.0,r,n_gridpoints,2); 
  // mesh in time 
  auto time_mesh = LinOps::make_mesh(0.0, 3.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  LinOps::DiscretizationXD my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  BumpFunc bump_1d{.L=0.0, .R=r, .c=r/2, .h=5.0, .focus=25.0};
  auto f_lambda = [&](double x, double y){ return bump_1d(x) * bump_1d(y); }; 
  my_vals.set_init(my_mesh, f_lambda); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::DirectionalNthDerivOp;
  auto expr = 0.2 * D(2,0) + 0.2 * D(2,1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_neumann(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCListXD>();
  bcs->list.emplace_back(left,right); 
  bcs->list.emplace_back(left,right); 

  // LHS time derivs ----------------------------------------------------------------
  auto lhs_expr = TExprs::NthTimeDeriv(2); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(2, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenSolver s(lhs_expr, expr); 
  auto v = s.Calculate(args); 
  // auto v = s.CalculateImp(args); 

  // Printing ---------------------------------------------------------------- 
  // print_vec(my_vals, "ICs"); 
  // print_vec(v,"Sol"); 

  print_mat(my_mesh->OneDim_views(my_vals.values(),0), "ICs"); 
  print_mat(my_mesh->OneDim_views(v.values(),0), "Sol"); 
};