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

#include<LinOps/All.hpp> 
#include<OutsideSteps/All.hpp> 
// #include<FDStencilsXD/All.hpp> // likewise ...
// #include<LinOpsXD/All.hpp>
// #include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  double t = 0.0; 
  auto mesh = LinOps::make_mesh(0.0,6.0,7); 

  LinOps::MatrixStorage_t R = LinOps::RandLinOp(mesh).GetMat().sparseView(); 

  auto v = LinOps::Discretization1D().set_init(mesh, [](double x){return 20.0 + x*x;}).values(); 

  auto bcs = OSteps::BCPair( OSteps::DirichletBC(0.0), OSteps::NeumannBC(2.0));

  cout << R << endl << "-------" << endl << v.transpose() << endl; 
  
  bcs.BeforeImpStep(t,mesh,R,v); 
  cout << R << endl << "-------" << endl << v.transpose() << endl; 

};