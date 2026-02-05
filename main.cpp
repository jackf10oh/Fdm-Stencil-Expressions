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
  int n_gridpoints = 11;

  // mesh in space [0, 10] 
  auto my_mesh = LinOps::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time [0, 2]
  auto time_mesh = LinOps::make_mesh(0.0, 2.0, 17); 

  using LinOps::NthDerivOp;
  NthDerivOp fdm_stencil; 
  fdm_stencil.set_mesh(my_mesh); 
  cout << fdm_stencil.GetMat() << endl; 
};