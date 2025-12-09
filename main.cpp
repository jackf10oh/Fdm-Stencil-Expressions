// main.cpp
//
//
//
// JAF 12/8/2025

#include<iostream>
#include<iomanip>
#include<vector>
#include<eigen3/Eigen/Dense>
// #include<eigen3/unsupported/Eigen/KroneckerProduct>
#include "incl/All.hpp"
#include "LinOps/All.hpp"

typedef double Real;
using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 
  auto r = 10.0; 
  auto my_mesh = std::make_shared<Mesh1D>(0.0,r,101); 

  Discretization1D my_vals;
  // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  // auto func = [](double x){return x*x;}; // x^2 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; // - x(x-10) / 10 
  my_vals.set_init(my_mesh, func); 
  print_vec(my_vals.begin(),my_vals.end(), "t=t0");



  NthDerivOp D(1); // 1st order derivative 
  auto left = std::make_shared<DirichletBC>(0.0);
  auto right = left; 
  
  double dt = 0.01;
  auto DeltaX = -0.2*D.compose(D) + 0.5*D; 
  auto Explicit_Step = IOp(my_mesh) + (dt) * DeltaX;
  Explicit_Step.lbc_ptr = left; 
  Explicit_Step.rbc_ptr = right; 
  Explicit_Step.set_mesh(my_mesh);

  double T = 10;
  int NSteps = T/dt; 
  // cout << Explicit_Step.Rhs().GetMat() << endl;
  for(int n=0; n<NSteps; n++)
  {
    my_vals = Explicit_Step.solve_implicit(my_vals); 
  }

  print_vec(my_vals.begin(),my_vals.end(), "t=T");

};
