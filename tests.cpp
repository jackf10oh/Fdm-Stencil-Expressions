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
#include "LinOpsXD/All.hpp"
#include "FDStencilsXD/DiffOps/DirectionalNthDerivOp.hpp"

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto my_meshes = LinOps::make_meshXD(0.0,2.0, 3, 3);
  // auto my_meshes = LinOps::make_mesh(0.0,3.0, 4);
  
  LinOps::IOpXD I; 
  I.set_mesh(my_meshes); 
  // cout << I.GetMat() << endl; 

  auto D = Fds::DirectionalNthDerivOp(1, 2); // first order , X direction 
  D.set_mesh(my_meshes); 
  cout << D.GetMat().toDense() << endl; 

  LinOps::DirectionalRandOp R(2); 
  R.set_mesh(my_meshes); 
  cout << "Random ----------------------" << endl << R.GetMat().toDense() << endl; 

};

// auto lam00 = [](){return 0.0;}; 
// auto lam01 = [](double x){return x*x - x + 1.5;}; 
// auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
// auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 
