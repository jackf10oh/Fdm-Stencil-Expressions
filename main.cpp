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

// #include<DiffOps/DiffOpTraits.hpp>
// #include<DiffOps/DiffOpBase.hpp>
// #include<DiffOps/NthDerivOp.hpp>

#include<DiffOps/NthPartialDeriv.hpp>

#include<Utilities/BlockDiagExpr.hpp>
#include<Utilities/HighDimExpr.hpp> 

using std::cout, std::endl, std::cin;

template<std::size_t N, std::size_t Dir>
using D = DiffOps::NthPartialDeriv<N,Dir>;

int main()
{
  // iomanip 
  std::cout << std::setprecision(6); 

  auto domain = LinOps::make_meshXD(0.0,5.0,6,2); 

  auto my_op = D<1,0>{};   

  my_op.set_mesh(domain); 

  cout << my_op.GetMat() << endl;

}