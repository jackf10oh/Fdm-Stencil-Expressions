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

#include "FDStencils/All.hpp" // must be first for plugin macro.
// #include "FDStencils/Interp1D.hpp" // Interp1D class 
#include "LinOps/All.hpp" 
// #include "FDStencilsXD/All.hpp" // likewise ...
// #include "LinOpsXD/All.hpp"

#include "Utilities/PrintVec.hpp"
#include "Utilities/BumpFunc.hpp"

#include "LhsExpressions/All.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

template<typename TUP_T>
void tup_print_impl(const TUP_T& tup){std::cout << std::endl;}; 

template<typename TUP_T, std::size_t I, std::size_t... Trailing>
void tup_print_impl(const TUP_T& tup){ std::cout << std::get<I>(tup).toString() << ", "; tup_print_impl<TUP_T,Trailing...>(tup);};

template<typename TUP_T, std::size_t... Is>
void tup_unpack(const TUP_T& tup, const std::index_sequence<Is...>& idx_seq)
{
  tup_print_impl<TUP_T,Is...>(tup);
}

template<typename TUP_T>
void tup_print(const TUP_T& tup)
{
  constexpr std::size_t N = std::tuple_size_v<TUP_T>; 
  tup_unpack(tup, std::make_index_sequence<N>{}); 
}

int main()
{
  // iomanip 
  std::cout << std::setprecision(4); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 11;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 4.0, 401); 

  // Coeff Operator ------------------------------------------
  AutonomousCoeff c = [](double x ){return x*x; }; 
  c.set_mesh(my_mesh); 

  // FornCalc -----------------------------------------------------------
  FornCalc calc(3,2); 

  std::vector<double> v = {0.0,1.0,2.0}; 

  // auto weights = calc.GetWeights(0.0, v.begin(), v.end(), 2); 
  calc.Calculate(0.0, v.begin(), v.end(), 2); 
  print_vec(calc.m_arr, "weights"); 
  
  // LhsExpr -------------------------------------------------------- 
  using D = NthTimeDeriv; 
  auto Ut = D(1); 
  auto Utt = -1.0 * D(2); 
  // auto Utt = -1.0 * D(2); 
  auto Uttt = D(3);  

  // auto sum_expr01 = Utt + Ut; 
  // auto sum_expr02 = Uttt + sum_expr01;
  auto sum_expr03 = (Uttt + Utt) + Ut; 

  // cout << std::tuple_size_v<std::remove_reference_t<decltype(sum_expr03.toTuple())>> << endl;

  tup_print(sum_expr03.toTuple());

  cout << Ut.CoeffAt(calc.m_arr, 3, 0) << endl; 
  cout << std::get<1>(sum_expr03.toTuple()).CoeffAt(calc.m_arr, 3,1) << endl; 
  // cout << sum_expr03.CoeffAt(calc.m_arr, 3, 0); // deleted member function  
};
 