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
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

typedef double Real;
using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 
  auto r = 10.0; 
  int n_gridpoints = 101;
  MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints); 
  
  Discretization1D my_vals;
  // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  // auto func = [](double x){return x*x;}; // x^2 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; // Bump centered at r/2. Zero at 0.0 and r 
  my_vals.set_init(my_mesh, func); 
  print_vec(my_vals.begin(),my_vals.end(), "t=t0");  

  // assembling scheme
  using D = NthDerivOp;
  double dt = 0.01;

  auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2) + 0.5*D(1));
  auto Diffusion = IOp() + (dt)*(-0.2*D(2)); 
  auto Convection = IOp() + (dt)*(-0.5*D(1)); 

  // setting boundary condition 
  auto left = std::make_shared<DirichletBC>(0.0); 
  auto right = left; 
  fdm_scheme.lbc_ptr = left; 
  fdm_scheme.rbc_ptr = right;
  Diffusion.lbc_ptr = left; 
  Diffusion.rbc_ptr = right;
  Convection.lbc_ptr = left; 
  Convection.rbc_ptr = right;
 
  // set mesh 
  fdm_scheme.set_mesh(my_mesh);
  Diffusion.set_mesh(my_mesh); 
  Convection.set_mesh(my_mesh); 

  double T = 10.0;
  int NSteps = T/dt; 
  for(int n=0; n<NSteps; n++)
  {
    // explcit euler 
    // my_vals = Explicit_Step.explicit_step(my_vals); 

    //implicit euler
    // my_vals = fdm_scheme.solve_implicit(my_vals);

    // Operator splitting methods
    auto temp = Convection.apply(my_vals); 
    my_vals = Diffusion.solve_implicit(temp); 
  }

  print_vec(my_vals.begin(),my_vals.end(), "t=T");
  
};
