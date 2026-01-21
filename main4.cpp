// main4.cpp 
//
// Testing out GenSolver on 1D PDEs 
//
// JAF 1/17/2026 

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include "Utilities/PrintVec.hpp"
#include "Utilities/BumpFunc.hpp"

#include "FDStencils/All.hpp" // must be first for plugin macro. 
#include "LinOps/All.hpp" 
// #include "FDStencilsXD/All.hpp" // likewise ...
// #include "LinOpsXD/All.hpp"
#include "TExprs/All.hpp"

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(3); 

  // defining Domain Mesh --------------------------------------
  auto r = 10.0; 
  int n_gridpoints = 16;
  // mesh in space 
  auto my_mesh = Fds::make_mesh(0.0,r,n_gridpoints); 
  // mesh in time 
  auto time_mesh = Fds::make_mesh(0.0, 5.0, 8); 

  // Initializing IC discretizations -------------------------------------------------------
  Fds::Discretization1D my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  BumpFunc f{.L=0.0, .R=r, .c=r/2, .h=1.0, .focus=20.0};
  my_vals.set_init(my_mesh, f); 

  // building RHS expression -----------------------------------------------------
  using D = Fds::NthDerivOp;
  auto expr = 0.2 * D(2) - 0.5 * D(1); 
  // auto expr = 0.2 * D(2) - 0.5 * D(1); 

  // Boundary Conditions + --------------------------------------------------------------------- 
  std::shared_ptr<Fds::IBCLeft> left = Fds::make_dirichlet(0.0); 
  std::shared_ptr<Fds::IBCRight> right = Fds::make_neumann(0.0); 

  auto bcs = std::make_shared<Fds::BCPair>(left,right); 

  // LHS time derivs ----------------------------------------------------------------
  // auto lhs_expr = NthTimeDeriv(1); 
  auto lhs_expr = TExprs::NthTimeDeriv(1); 

  // Solving --------------------------------------------------------------------- 
  TExprs::GenSolverArgs args{
    .domain_mesh_ptr = my_mesh,
    .time_mesh_ptr = time_mesh,
    .bcs = bcs, 
    .ICs = std::vector<Eigen::VectorXd>(1, my_vals.values()), 
    .time_dep_flag = false 
  }; 

  TExprs::GenSolver s(lhs_expr, expr, TExprs::PrintWrite{}); 
  // auto v = s.Calculate(args); 
  // auto v = s.CalculateImp(args); 
  s.CalculateImp(args); 

  // Printing ---------------------------------------------------------------- 
  print_vec(my_vals, "ICs"); 
  // print_vec(v,"Sol"); 
};

/* 

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

// std::cout << "Return type traits --------------------------------------------" << endl; 
// cout << typeid(TimeDerivTraits<decltype(Ut)>::CoeffAtReturnType).name() << endl;  
// cout << typeid(TimeDerivTraits<decltype(c*Ut)>::CoeffAtReturnType).name() << endl;

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

*/ 