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
  
};