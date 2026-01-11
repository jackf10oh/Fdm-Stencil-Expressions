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

#include "FDStencils/All.hpp" // must include first for plugin to take effect over linops?
#include "LinOps/All.hpp" 
#include "FDStencilsXD/FdmPluginXD.hpp"
#include "FDStencilsXD/All.hpp"
#include "LinOpsXD/All.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

auto lam00 = [](){return 0.0;}; 
auto lam01 = [](double x){return x*x - x + 1.5;}; 
auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto my_meshes = LinOps::make_meshXD(0.0,1.0, 3, 3);
  // auto my_meshes = LinOps::make_mesh(0.0,3.0, 4);

  IOpXD I(my_meshes); 
  DirectionalRandOp R(my_meshes, 0);  

  // I.SetTime(1.0); 
  // R.SetTime(2.0); 
  
  auto expr = I + R; 

  I.SetTime(1.0); 
  R.SetTime(2.0); 

  cout << I.Time() << endl; 
  cout << R.Time() << endl; 
  cout << expr.Time() << endl; 

  cout << internal::is_exprxd_crtp<decltype(expr)>::value << endl; 

  // cout << R.GetMat().toDense() << endl; 

};

