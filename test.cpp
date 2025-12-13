// test.cpp
//
//
//
// JAF 12/5/2025

#include<iostream>
#include<iomanip>
#include<vector>
#include<chrono>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>
#include<eigen3/Eigen/SparseLU>
#include<eigen3/unsupported/Eigen/KroneckerProduct>

#include "incl/All.hpp"
#include "LinOps/All.hpp"

#include "incl/Utilities/FornbergCalc.hpp"

using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 
  // auto my_mesh = make_mesh(0.0,5.0,6); 

  // Discretization1D my_vals;
  // // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  // auto func = [](double x){return x*x;}; // x^2 
  // my_vals.set_init(my_mesh, func); 

  using std::chrono::system_clock;
  using std::chrono::duration_cast;
  using Time_t = std::chrono::microseconds; 
  auto time_lambda = [](const MeshPtr_t& m){
    // TEST 1 
    std::cout << "------------" << std::endl; 
    NthDerivOp D(1); 
    auto start = system_clock::now(); 
    D.set_mesh(m);
    auto end = system_clock::now(); 
    auto d = duration_cast<Time_t>(end-start); 
    std::cout << "Mesh size:" << m->size() << " time:"  << d.count() << std::endl; 
  };

  auto mesh01 = make_mesh(0.0,1e3,1e3 + 1); 
  auto mesh02 = make_mesh(0.0,1e4,1e4 + 1); 
  auto mesh03 = make_mesh(0.0,1e5,1e5 + 1); 
  auto mesh04 = make_mesh(0.0,1e6,1e6 + 1); 

  time_lambda(mesh01);
  time_lambda(mesh02);
  time_lambda(mesh03);
  time_lambda(mesh04);

}