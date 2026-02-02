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
#include<OutsideSteps/BoundaryCondsXD/BCList.hpp> 

// #include<FDStencilsXD/All.hpp> // likewise ...
// #include<LinOpsXD/All.hpp>
// #include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

using std::cout, std::endl;

template<typename L, typename R>
using p = OSteps::BCPair<L,R>; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  double t = 0.0; 
  std::vector<std::pair<double,double>> ends = {{0.0,1.0},{0.0,1.0},{0.0,1.0}}; 
  std::vector<std::size_t> steps = {6,6,4};  
  auto mesh = LinOps::make_meshXD(ends, steps); 

  // auto v = LinOps::DiscretizationXD().set_init(mesh, [](double x, double y){ return x + y*y; }).values(); 

  // print_mat(mesh->OneDim_views(v), "OneDim Views(0)"); 

  auto bc1 = OSteps::DirichletBC(1.0); 
  auto bc2 = OSteps::NeumannBC(2.0); 
  auto bc3 = OSteps::DirichletBC(3.0); 
  // auto bcs_list = OSteps::BCList( p(bc1,bc1), p(bc2,bc2) ); 
  auto bcs_list = OSteps::BCList( p(bc1,bc1), p(bc2,bc2), p(bc3,bc3) ); 

  // bcs_list.SolAfterStep<OSteps::FDStep_Type::EXPLICIT>(t,mesh,v); 
  // bcs_list.SolBeforeStep<OSteps::FDStep_Type::IMPLICIT>(t,mesh,v); 
  // print_mat(mesh->OneDim_views(v), "OneDim Views(0)"); 


  LinOps::MatrixStorage_t R = (22 * LinOps::IOp(mesh)).GetMat(); 
  // cout << "R" << endl << R << endl; 

  bcs_list.MatBeforeStep(t,mesh,R); 
  cout << "R" << endl << R << endl; 

};