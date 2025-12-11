// test.cpp
//
//
//
// JAF 12/5/2025

#include<iostream>
#include<iomanip>
#include<vector>
#include<eigen3/Eigen/Dense>
#include<eigen3/unsupported/Eigen/KroneckerProduct>

#include "incl/All.hpp"
#include "LinOps/All.hpp"

using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 
  auto my_mesh = make_mesh(0.0,4.0,5); 

  Discretization1D my_vals;
  // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  auto func = [](double x){return x*x;}; // x^2 
  my_vals.set_init(my_mesh, func); 

  TOp t(my_mesh); // TOp inherits from coeffopbase 
  t.SetTime(104.0);
  NthDerivOp D(my_mesh,1); // 1st derivative operator

  // cout << D.GetMat() << endl; 

  auto expr = D.compose(t.compose(D)); 
  cout << expr.GetMat() << endl; 

  // auto D2 = D.compose(D); 
  // cout << D2.GetMat() << endl; 

  auto scaled = 40.0 * D; 
  cout << D.compose(scaled).GetMat() << endl; 

}