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


  TOp t(my_mesh), t2(my_mesh); 

  cout << "--------------" << endl;
  t.SetTime(15.0);  
  cout << t.GetMat() << endl; 

  cout << "--------------" << endl; 
  auto expr = t*t; 
  expr.SetTime(20.0); 
  cout << expr.GetMat() << endl; 

  cout << "--------------" << endl; 
  cout << t.Time() << endl; 
  t.SetTime(10.0); 
  cout << t.GetMat() << endl; 
  cout << t.Time() << endl; 
  expr.SetTime(1.2);
  cout << t.GetMat() << endl; 
  cout << t.Time() << endl; 


  cout << "--------------" << endl; 
  NthDerivOp D1(my_mesh), D2(my_mesh); 
  auto expr3 = D1.compose(D2); 
  expr3.SetTime(10.0);
  cout << "expr is expr? " << is_expr_crtp<decltype(expr3)>::value << endl;
  cout << "L is op? " << is_linop_crtp<decltype(D1)>::value << endl;
  cout << "R is op? " << is_linop_crtp<decltype(D2)>::value << endl;

  cout << "--------------" << endl; 
  RandLinOp L1, L2; 
  auto expr2 = L1.compose(L1); 
  expr2.SetTime(10.0);
  cout << "expr is expr? " << is_expr_crtp<decltype(expr2)>::value << endl;
  cout << "L is op? " << is_linop_crtp<decltype(L1)>::value << endl;
  cout << "R is op? " << is_linop_crtp<decltype(L2)>::value << endl;

}