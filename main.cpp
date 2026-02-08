// main.cpp
//
//
//
// JAF 12/8/2025

// #include "DiffOps/DiffOpTraits.hpp"
// #include "DiffOps/NthDerivOp.hpp"

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

// #include<LinOps/Operators/DiffOps/experimental_NthDerivOp.hpp>
#include<LinOps/All.hpp> 
#include<OutsideSteps/All.hpp> 
#include<TExprs/All.hpp> 

#include<Utilities/PrintVec.hpp>
#include<Utilities/BumpFunc.hpp>

#include<DiffOps/DiffOpTraits.hpp>
#include<DiffOps/DiffOpBase.hpp>
#include<DiffOps/NthDerivOp.hpp>

#include<Utilities/BlockDiagExpr.hpp>

using std::cout, std::endl, std::cin;

template<std::size_t N>
using D = DiffOps::NthDerivOp<N>;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // mesh in space [0, 10]
  auto r = 5.0;
  int n_gridpoints = 6; 
  auto domain = LinOps::make_meshXD(0.0,r,n_gridpoints, 2); 

  // defining matrices 
  auto A = LinOps::DirectionalNthDerivOp(domain, 1,1).GetMat(); 
  cout << A << endl; 
}