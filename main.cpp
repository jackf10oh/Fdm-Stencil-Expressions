// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<Eigen/Dense>
// #include<unsupported/Eigen/KroneckerProduct>
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

typedef double Real;
using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 101;
  MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints); 
  
  // Initializing IC discretizations -------------------------------------------------------
  Discretization1D my_vals;
  // Bump centered at r/2. Zero at 0.0 and r. 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);};  
  my_vals.set_init(my_mesh, func); 
  // print to cout 
  print_vec(my_vals.begin(),my_vals.end(), "t=t0");  

  // defining boundary of T + meshing in time  ----------------------------------------------
  double T = 5.0;
  int NSteps = 300; 
  double dt = T/NSteps;

  // assembling scheme -----------------------------------------------------------
  using D = NthDerivOp;
  auto fdm_scheme = IOp(my_mesh) + (dt) * (-0.2*D(2)); // + 0.5*D(1)

  // setting boundary condition --------------------------------------------------------
  auto left = std::make_shared<DirichletBC>(0.0);; 
  auto right = std::make_shared<RobinBC>(1.0,1.0,0.0); 
  fdm_scheme.lbc_ptr = left; 
  fdm_scheme.rbc_ptr = right;

  // set mesh -----------------------------------------------------------------------
  fdm_scheme.set_mesh(my_mesh);


  // execute -----------------------------------------
  for(int n=0; n<NSteps; n++)
  {
    // explcit euler 
    // my_vals = fdm_scheme.explicit_step(my_vals); 

    //implicit euler
    my_vals = fdm_scheme.solve_implicit(my_vals);
  }

  // output -------------------------------------------------
  print_vec(my_vals.begin(),my_vals.end(), "t=T");
  
};

  // auto Diffusion = IOp() + (dt)*(-0.2*D(2)); 
  // auto Convection = IOp() + (dt)*(-0.5*D(1)); 

  // Diffusion.lbc_ptr = left; 
  // Diffusion.rbc_ptr = right;
  // Convection.lbc_ptr = left; 
  // Convection.rbc_ptr = right;
 

  // Diffusion.set_mesh(my_mesh); 
  // Convection.set_mesh(my_mesh); 

  //implicit euler
  // my_vals = fdm_scheme.solve_implicit(my_vals);

  // Operator splitting methods
  // auto temp = Convection.apply(my_vals); 
  // my_vals = Diffusion.solve_implicit(temp); 