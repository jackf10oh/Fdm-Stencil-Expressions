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

#include "DiffOps/All.hpp" // must include first for plugin to take effect over linops?
#include "LinOps/All.hpp" 
// #include "LinOpsXD/All.hpp"

#include "DiffOps/CoeffOps/CoeffOpBase.hpp"
#include "DiffOps/CoeffOps/TimeDepCoeff.hpp"
#include "DiffOps/CoeffOps/AutonomousCoeff.hpp"

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto my_mesh = make_mesh(0.0,5.0,6); 

};

// auto lam00 = [](){return 0.0;}; 
// auto lam01 = [](double x){return x*x - x + 1.5;}; 
// auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
// auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 
