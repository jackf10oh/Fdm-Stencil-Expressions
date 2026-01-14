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

  // defining Domain + Time  Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 21;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 5.0, 101); 

  // Making Discretization + setting ICs ----------------------------- 
  Discretization1D my_vals; 
  BumpFunc f{.L = 0.0, .R=r, .c=7.5, .h=4.0}; 
  my_vals.set_init(my_mesh, f); 

  // Printing ----------------------------------------------
  print_vec(my_vals, "Bump"); 

};

// auto get_interval = [](const auto& v, double c){
//   if(v.size() < 2) throw std::runtime_error("size of v < 2"); 
//   if(c < v.cbegin()[0]) throw std::runtime_error("c < v[0]"); 
//   auto after = std::lower_bound(v.cbegin(), v.cend(), c);
//   if(after == v.cbegin()) return std::pair(after, after+1); 
//   if(after == v.cend()) throw std::runtime_error("right bound == v.cend()") 
//   return std::pair(after-1, after); 
// } 