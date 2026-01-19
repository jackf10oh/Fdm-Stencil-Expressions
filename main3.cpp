// main3.cpp 
//
// Testing out LhsExecutor on PDEs 
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
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 5.0, 501); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::Discretization1D my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  BumpFunc f{.L=0.0, .R=r, .c=r/2, .h=1.0, .focus=20.0};
  my_vals.set_init(my_mesh, f); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto expr = 0.2 * D(2); 
  // auto expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_dirichlet(0.0); 

  auto bcs = std::make_shared<BCPair>(left,right); 
  // bcs->pair.first=left; 
  // bcs->pair.second=right; 

  // LHS time derivs ----------------------------------------------------------------
  // auto lhs_expr = NthTimeDeriv(1); 
  auto lhs_expr = NthTimeDeriv(2); 

  // Solving --------------------------------------------------------------------- 
  LhsExecutor exec(lhs_expr); 
  
  cout << exec.m_num_nodes << endl; 
  cout << exec.m_order << endl; 
  cout << exec.m_stored_sols.size() << endl; 

  auto it = time_mesh->cbegin(); 
  auto end = time_mesh->cend();

  exec.ConsumeSolution(my_vals.values()); 
  exec.ConsumeSolution(my_vals.values()); 
  exec.ConsumeTime(*it++);
  exec.ConsumeTime(*it++);

  expr.set_mesh(my_mesh); 

  for(; it!= end; it++)
  {
    Eigen::VectorXd rhs = exec.BuildRhs(*it); 
    // print_vec(exec.m_weights_calc.m_arr, "Weights"); 
    
    // Explicit Step 
    Eigen::VectorXd next_sol = exec.inv_coeff_util()*expr.GetMat()*exec.MostRecentSol() + rhs; 
    
    // Apply BCs 
    bcs->SetSol(next_sol, my_mesh); 
    exec.ConsumeSolution(next_sol);
    exec.ConsumeTime(*it); 
  }

  // Printing ---------------------------------------------------------------- 
  print_vec(my_vals, "ICs"); 
  print_vec(exec.m_stored_sols[0],"Sol"); 
  // auto foo = exec.BuildRhs(*it); 
};
