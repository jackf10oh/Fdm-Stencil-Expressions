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


// ========================================================
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

// ========================================================================
// take a tuple and map each element to an singelton or empty tuple 
template<typename TUP_T, typename TUP_U>
auto transform_tup_impl(TUP_T& tup, TUP_U result)
{
  return result; 
}; 

template<typename TUP_T, typename TUP_U, std::size_t I, std::size_t... Trailing>
auto transform_tup_impl(TUP_T& tup, TUP_U result)
{ 
  auto result_appended = std::tuple_cat(
    result,
    std::make_tuple(std::make_tuple(std::get<I>(tup)))
  ); 
  return transform_tup_impl<TUP_T,decltype(result_appended), Trailing...>(tup, result_appended);
}; 

template<typename TUP_T, std::size_t... Is>
auto transform_tup_unpack(TUP_T& tup, const std::index_sequence<Is...>& idx_seq)
{
  return transform_tup_impl<TUP_T, std::tuple<>, Is...>(tup, std::tuple<>{}); 
}

template<typename TUP_T>
auto transform_tup(TUP_T& tup)
{
  return transform_tup_unpack(tup, std::make_index_sequence<std::tuple_size_v<TUP_T>>{}); 
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
  auto Utt = D(2); 
  // auto cUt = 2.0 * Ut;   
  auto cUt = c * Ut;  
  auto ccUt = 3.0 * cUt; 

  std::cout << "left mult -------------------------------------" << std::endl; 
  cout << Ut.CoeffAt(calc.m_arr, 3, 0) << endl; 
  // cout << cUt.CoeffAt(calc.m_arr, 3, 0) << endl; 
  // cout << ccUt.CoeffAt(calc.m_arr, 3, 0) << endl; 

  std::cout << "unary negate -------------------------------------" << std::endl; 
  auto negation = - D(1); 
  cout << negation.CoeffAt(calc.m_arr, 3, 0) << endl; 


  std::cout << "binary subtract -------------------------------------" << std::endl; 
  auto sum_expr01 = Utt + Ut; 
  auto sum_expr02 = Utt - Ut; 
  auto sum_expr03 = 2.0*Ut - Ut; 

  // cout << std::tuple_size_v<std::remove_reference_t<decltype(sum_expr02.toTuple())>> << endl;
  // cout << std::tuple_size_v<std::remove_reference_t<decltype(sum_expr01.toTuple())>> << endl;
  cout << std::get<1>(sum_expr01.toTuple()).CoeffAt(calc.m_arr, 3,0) << endl; 
  cout << std::get<1>(sum_expr02.toTuple()).CoeffAt(calc.m_arr, 3,0) << endl; 
  // cout << std::get<1>(sum_expr02.toTuple()).toString() << endl; 
  // cout << sum_expr03.CoeffAt(calc.m_arr, 3, 0); // deleted member function 

  std::cout << "Lhs Executor -------------------------------------" << std::endl; 
  LhsExecutor executor(sum_expr02); 
  cout << "executor order: " << executor.m_order << endl; 
  cout << "executor num_nodes: " << executor.m_num_nodes << endl; 

  /* std::cout << "BOQ Tuple -------------------------------------" << std::endl; 
  auto tup = std::make_tuple( "foo", "bar", 1.0, 1, 2.0 ); 

  auto singletons = transform_tup(tup); 

  cout << std::get<0>(tup) << endl; 

  cout << std::get<0>(std::get<1>(singletons)) << std::endl; 

  auto recombined = std::apply([](auto&&... args){ return std::tuple_cat(args...);}, singletons);
  
  cout << std::get<0>(recombined) << endl; 

  */ 

  /* std::cout << "Return type traits --------------------------------------------" << endl; 
  cout << typeid(TimeDerivTraits<decltype(Ut)>::CoeffAtReturnType).name() << endl;  
  cout << typeid(TimeDerivTraits<decltype(c*Ut)>::CoeffAtReturnType).name() << endl;
  
  */ 

  std::cout << "Lhs Executor II ----------------------------------------" << std::endl; 
  // auto sum_expr04 = Utt + Ut + cUt + ccUt; 
  auto sum_expr04 = Ut + cUt + ccUt; 
  LhsExecutor exec(sum_expr04); 

  std::cout << std::tuple_size_v<decltype(exec.m_scalar_coeff_sum_partition)> << std::endl; 
  std::cout << std::tuple_size_v<decltype(exec.m_mat_coeff_sum_partition)> << std::endl; 

  exec.ConsumeTime(0.0); 
  exec.ConsumeTime(1.0); 
  // exec.BuildRhs(2.0); 
  exec.SetWeightsFromTime(2.0); 

  print_vec(exec.m_weights_calc.m_arr, "exec weights");

  // cout << exec.inv_coeff_util() << endl; 

  std::cout << "CoeffMultExpr II ----------------------------------------" << std::endl; 
  AutonomousCoeff c01 = [](double x){return x*x + 10.0; }; 
  TimeDepCoeff c02 = [](double t, double x){return t + x; };
  auto expr01 = c01 * Ut; 
  auto expr02 = c02 * Utt;  

  expr01.set_mesh(my_mesh); 
  std::cout << expr01.Lhs().GetMat() << std::endl; 

  expr02.set_mesh(my_mesh); 
  expr02.SetTime(20.0); 
  std::cout << expr02.Lhs().GetMat() << std::endl; 

  
};
