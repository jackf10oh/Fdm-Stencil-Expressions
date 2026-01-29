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

#include "FDStencils/include/FDStencils/All.hpp" // must be first for plugin macro.
#include "LinOps/include/LinOps/All.hpp" 
#include "FDStencilsXD/include/FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/include/LinOpsXD/All.hpp"
#include "TExprs/include/TExprs/All.hpp" 

#include "Utilities/include/Utilities/PrintVec.hpp"
#include "Utilities/include/Utilities/BumpFunc.hpp"

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto my_mesh = LinOps::make_mesh(0.0,4.0,5); 
  double c = 2.0; 
  LinOps::IOp I(my_mesh);
  auto expr = c * I; 
  cout << expr.GetMat() << endl; 
  c = 4.0; 
  cout << expr.GetMat() << endl; 



  auto my_meshXD = LinOps::make_meshXD(0.0,4.0,5); 
  double c2 = 2.0; 
  LinOps::IOpXD I2(my_meshXD);
  auto expr2 = c2 * I2; 
  cout << expr2.GetMat() << endl; 
  c2 = 4.0; 
  cout << expr2.GetMat() << endl; 
  
};