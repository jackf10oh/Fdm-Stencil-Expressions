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

typedef double Real;
using std::cout, std::endl;

struct foo
{
  Discretization1D apply(const Discretization1D& d) const { return d; };
};

struct bar{};

int main()
{
  std::cout << std::setprecision(2); 
  auto my_mesh = std::make_shared<Mesh1D>(0.0,4.0,5); 

  Discretization1D my_vals;
  // auto func = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // 2x^3 - 5x^2 + 3x -1 
  auto func = [](double x){return x*x;}; // x^2 
  my_vals.set_init(my_mesh, func); 

  // IOp I(my_mesh); 
  // RandLinOp L1(my_mesh); 
  // RandLinOp L2(my_mesh);  

  // cout << "Identity" << endl << I.GetMat() << endl;
  
  // cout << "-------------" << endl;
  // cout << "L1" << endl << L1.GetMat() << endl; 
  
  // cout << "-------------" << endl;
  // cout << "L2" << endl << L2.GetMat() << endl; 

  // cout << "-------------" << endl << "plugins!" << endl;
  // I.print();
  // L1.print();
  // L2.print();

  // cout << "---------------" << endl << "DiffOps" << endl; 
  // auto D = NthDerivOp(my_mesh); 
  // // cout << D.GetMat() << endl; 
  // auto Explicit_Op = IOp() + 0.1 * D.compose(D); 
  // Explicit_Op.set_mesh(my_mesh); 

  // cout << Explicit_Op.GetMat() << endl; 

  // auto c = Explicit_Op.apply(my_vals);

  // std::cout << "result" << c.values() << endl; 

  CoeffOp coeff_op(my_mesh); 
  coeff_op.SetTime(2.0);
  // cout << coeff_op.GetMat() << endl; 

  // auto comp_op = 2.0*(2.0*(2.0*(2.0*coeff_op))); 
  // // cout << comp_op.GetMat() << endl; 
  // comp_op.SetTime(5.0); 
  // cout << comp_op.Time() << endl; 
  // std::cout << "comp_op is expr? " << is_expr_crtp<decltype(comp_op)>::value << endl; 
  // std::cout << "comp_op is op? " << is_linop_crtp<decltype(comp_op)>::value << endl; 

  cout << "foo has apply? " << has_apply<foo>::value << endl;
  cout << "bar has apply? " << has_apply<bar>::value << endl;
}