// main6.cpp 
//
// Testing out GenInterp (wraps GenSolver) on 1D PDEs 
//
// JAF 1/120/2026 

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

#include "TExprs/GenInterp.hpp" // GenInterp 

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
  int n_gridpoints = 11;
  // mesh in space 
  auto my_mesh = LinOps::make_meshXD(0.0,r,n_gridpoints,2); 
  // auto my_mesh = LinOps::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = LinOps::make_mesh(0.0, 3.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  // LinOps::Discretization1D my_vals;
  LinOps::DiscretizationXD my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  BumpFunc bump_1d{.L=0.0, .R=r, .c=r/2, .h=5.0, .focus=15.0};
  auto f_lambda = [&](double x, double y){ return bump_1d(x) * bump_1d(y); }; 
  my_vals.set_init(my_mesh, f_lambda); 

  // building RHS expression -----------------------------------------------------
  // using D = Fds::NthDerivOp; 
  using D = Fds::DirectionalNthDerivOp;
  auto expr = 0.2 * D(2,0) + 0.2 * D(2,1);  

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_neumann(-10.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(10.0); 

  // auto bcs = std::make_shared<Fds::BCPair>(left, right);

  auto bcs = std::make_shared<Fds::BCListXD>();
  bcs->list.emplace_back(left,right); 
  bcs->list.emplace_back(left,right); 

  // LHS time derivs ----------------------------------------------------------------
  auto lhs_expr = TExprs::NthTimeDeriv(1); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(1, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenInterp my_interp(lhs_expr, expr, args); 

  // manual printing ---------------------------------------------------------------- 
  my_interp.FillVals(); 
  auto views = my_mesh->OneDim_views(my_interp.m_data.at(100));
  print_mat(views, "Sol at t=3.0"); 

  // Interpolating! ---------------------------------------------------------------
  double t; 
  std::cout << "Enter t: "; 
  std::cin >> t; 

  std::cout <<"sol at " << t << "-----------------" << std::endl; 

  auto row_print = [&](auto x_it)
  {
    std::cout <<"["; 
    for(auto it = my_mesh->GetMeshAt(1)->cbegin(); it != my_mesh->GetMeshAt(1)->cend(); it++)
    {
      std::cout << my_interp.SolAt(t,*x_it,*it); 
      std::cout << ((it != std::prev(my_mesh->GetMeshAt(1)->cend())) ? "," : "]"); 
    } 
      std::cout << ((x_it != std::prev(my_mesh->GetMeshAt(0)->cend())) ? ",\n" : ""); 
  };

  std::cout << "["; 
  for(auto row_it = my_mesh->GetMeshAt(0)->cbegin(); row_it != my_mesh->GetMeshAt(0)->cend(); row_it++)
  {
    row_print(row_it);
  }
  std::cout << "]" << std::endl; 

};