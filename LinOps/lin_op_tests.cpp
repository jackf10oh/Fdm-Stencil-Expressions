// test.cpp
//
// 
//
// JAF 12/5/2025

#include <cstdint>
#include<iostream>
#include<vector>
#include<chrono>
#include<thread>

#include<eigen3/Eigen/Core> 
// #include<eigen3/unsupported/Eigen/KroneckerProduct>

#include "All.hpp"
// #include "LinOps/Mesh.hpp"
// #include "incl/Discretization.hpp"
// // #include "incl/Utilities/StopWatch.hpp"
// // #include "incl/Operators/RandLinOp.hpp"

using std::cout, std::endl;

int main()
{
  std::shared_ptr<Mesh1D> my_mesh = std::make_shared<Mesh1D>(0.0,10,11); 

  Discretization1D my_vals;

  // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  auto func = [](double x){return x*x;}; // x^2 
  my_vals.set_init(my_mesh, func); 

  // RandLinOp op_01(2); 
  // RandLinOp op_02(2); 

  // cout << "------------" << endl;
  // cout << op_01.GetMat() << endl;

  // cout << "------------" << endl;
  // cout << op_02.GetMat() << endl;

  // cout << "------------" << endl;
  // auto s = op_01.compose(op_02); 

  // cout << "composition is composition? " << is_compose_expr<decltype(s)>::value << endl; 
  // cout << "LinOp is composition? " << is_compose_expr<decltype(op_01)>::value << endl; 
  // cout << "int is composition? " << is_compose_expr<int>::value << endl; 

  // cout << s.GetMat() << endl; 

  // cout << "------------" << endl;
  // // auto sum = op_01 + s; 
  // auto sum = op_01.compose(s + op_02); 
  // op_01+op_02;
  // cout << (op_01 + op_02).GetMat() << endl;
  // cout << (RandLinOp(3) + RandLinOp(3)).GetMat() << endl;
  // cout << "sum is sum? " << is_add_expr<decltype(sum)>::value << endl; 
  // cout << "composition is sum? " << is_add_expr<decltype(s)>::value << endl; 
  // cout << "int is sum? " << is_add_expr<int>::value << endl; 
  // cout << "linop is sum? " << is_add_expr<decltype(op_01)>::value << endl; 
  // cout << sum.GetMat() << endl;

  // cout << "------------" << endl;
  // cout << "prod is product? " << is_scalar_multiply_expr<decltype(15.0*sum)>::value << endl;
  // auto prod = 15.0 * sum; 
  // auto fancy = RandLinOp(2).compose(prod); 
  // cout << fancy.GetMat() << endl;

  cout << "------------" << endl;
  auto I = IOp(); 
  I.set_mesh(my_mesh); 
  cout << I.GetMat() << endl;
  cout << "------------" << endl;
  cout << my_vals.values() << endl; 
  cout << "------------" << endl;
  cout << I.apply(my_vals).values().transpose() << endl;
}