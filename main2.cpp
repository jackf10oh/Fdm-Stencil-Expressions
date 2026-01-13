// main2.cpp
//
//
//
// JAF 1/12/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include "Utilities/PrintVec.hpp"

#include "FDStencils/All.hpp" // must be first for plugin macro. 
#include "LinOps/All.hpp" 
#include "FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/All.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 21;
  // mesh in space 
  auto my_meshXD = Fds::make_meshXD(0.0,r,n_gridpoints,2); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 5.0, 101); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::DiscretizationXD my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  auto func = [r, smush=10](double x, double y){
    auto bump1D = [&](double x){ return std::pow(x*(r-x)*(4.0/(r*r)),smush); }; 
    return bump1D(x)*bump1D(y);
  };  

  my_vals.set_init(my_meshXD, func); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::DirectionalNthDerivOp;
  auto expr = 0.2 * D(2,0) + 0.2 * D(2,1) - 0.5 * D(1,0); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_dirichlet(0.0); 

  auto bcs = std::make_shared<BCListXD>(); 
  bcs->list.emplace_back(left,right); 
  bcs->list.emplace_back(left,right); 

  // Solving --------------------------------------------------------------------- 
  Fds::SolverArgsXD args 
  {
    .domain_mesh_ptr   = my_meshXD,  
    .time_mesh_ptr     = time_mesh, 
    .bcs               = bcs,  
    .ICs               = my_vals, 
    .time_dep_flag     = true 
  };

  Fds::SolverXD s(expr); 
  // auto result = s.Calculate( args );
  auto result = s.CalculateImp( args );

  // Printing ---------------------------------------------------------------- 
  print_mat(my_vals.OneDim_views(), "ICs");  
  print_mat(result.OneDim_views(), "solution"); 

  // expr.set_mesh(my_mesh); 
  // std::cout << expr.GetMat() << std::endl; 
};