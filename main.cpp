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
#include<Utilities/HighDimExpr.hpp> 

using std::cout, std::endl, std::cin;

template<std::size_t N>
using D = DiffOps::NthDerivOp<N>;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  auto domain = LinOps::make_meshXD(0.0,3.0,4,3); 

  auto R = LinOps::DirectionalRandOp(domain, 1).GetMat(); 

  cout << R << endl; 

  cout << endl << "==========================" << endl; 

  LinOps::MatrixStorage_t Rsmall = LinOps::RandLinOp(domain->GetMesh(0)).GetMat().sparseView();
  cout << make_BlockDiag(make_HighDim(Rsmall, 4), 4);  
}