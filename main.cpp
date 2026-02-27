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
  int n_gridpoints = 21;
  auto domain = LinOps::make_mesh(0,r, n_gridpoints); 

  auto domain2 = LinOps::make_meshXD(0.0, r, n_gridpoints, 2); 

  // Making discretization -------------------------------
  BumpFunc my_func{.L=0.0, .R=10.0, .c=4.0, .h=4.0}; 
  auto v = LinOps::make_Discretization(domain, my_func); 

  auto v2 = LinOps::make_Discretization(domain2, [&](double x, double y){return my_func(x) * my_func(y);}); 

  // print -------------------- 
  print_vec(v, "discretization");

  print_mat(domain2->OneDim_views(v2.values()), "sol2d"); 
};