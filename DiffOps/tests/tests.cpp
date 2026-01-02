// tests.cpp
//
// FdmStencilExpressions Test suite with GTest
//
// JAF 12/5/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../All.hpp"

// testing dirichlet boundary conditions ===========================================================================
// can be constructed from args, default args
TEST(DirichletBcSuite, DirichletConstructible)
{
  // dirichlet bc == 0.0
  DirichletBC bc; 
  // bc == 2.0
  DirichletBC bc2(2.0);
  ASSERT_EQ(bc.boundary_val,0.0); 
  ASSERT_EQ(bc2.boundary_val,2.0); 
}

// test overrides of pure virtual funcs
TEST(DirichletBcSuite, DirichletOverrides)
{
  // mesh, discretization, and operator
  auto my_mesh = make_mesh(); 
  auto s = my_mesh->size(); 
  MatrixStorage_t A = RandLinOp(my_mesh).GetMat().sparseView(); 
  Discretization1D my_disc; 
  my_disc.set_init(my_mesh, [](double x){return x*x + 2*x + 2;});

  // boundary condition 
  auto bval = 3.0; 
  auto bc_ptr = std::make_shared<DirichletBC>(bval);

  // explicit methods 
  // sets left/right end of discretization 
  bc_ptr->SetSolL(my_disc, my_mesh);
  bc_ptr->SetSolR(my_disc, my_mesh);
  ASSERT_EQ(my_disc.at(0),bval); 
  ASSERT_EQ(my_disc.at(s-1),bval); 
  
  // implicit methods 
  // change disc values to some garbage value before implicit sets them back
  my_disc.at(0)=1.2345; 
  my_disc.at(s-1)=6.789; 
  // set left/right end of discretization in implicit scheme
  bc_ptr->SetImpSolL(my_disc, my_mesh);
  bc_ptr->SetImpSolR(my_disc, my_mesh);
  // sets first/last row in a stencil
  bc_ptr->SetStencilL(A, my_mesh); 
  bc_ptr->SetStencilR(A, my_mesh);

  // first/last entry of discretization 
  ASSERT_EQ(my_disc.at(0),bval); 
  ASSERT_EQ(my_disc.at(s-1),bval);

  // first row of A 
  ASSERT_EQ(A.coeff(0,0),1.0); 
  ASSERT_EQ(A.coeff(0,1),0.0); 
  ASSERT_EQ(A.coeff(0,s-1),0.0);
  // last row of A  
  ASSERT_EQ(A.coeff(s-1,0),0.0); 
  ASSERT_EQ(A.coeff(s-1,s-2),0.0); 
  ASSERT_EQ(A.coeff(s-1,s-1),1.0); 
}

// Testing TCoeff ===========================================================================  
// Testing TCoeff can be constructed
TEST(CoeffOpTestSuite, TCoeffConstructible)
{
  TCoeff t; 
};

// Testing TCoeff methods
TEST(CoeffOpTestSuite, TCoeffSettable)
{
  // make a mesh 
  MeshPtr_t my_mesh = make_mesh(0.0,4.0,5); 
  
  // make a TCoeff that acts on functions of mesh
  TCoeff t(my_mesh); 

  // lambda to check 4 corners + middle of a matrix expr
  auto check_lambda = [s = my_mesh->size()-1](const auto& expr, double val) -> void
  {
    MatrixStorage_t Mat = expr.GetMat(); 
    // Check that entries on diag are 2
    ASSERT_EQ(Mat.coeff(0,0),val); 
    ASSERT_EQ(Mat.coeff(s,s),val); 
    ASSERT_EQ(Mat.coeff(s/2,s/2),val);

    // of diag are zero
    ASSERT_EQ(Mat.coeff(0,s),0); 
    ASSERT_EQ(Mat.coeff(s,0),0); 
  };

  // check diag is 0.0, off diag is 0.0
  check_lambda(t,0.0); 

  // set diag to 2.0 
  t.SetTime(2.0);
  // check diag is 2.0, off diag is 0.0
  check_lambda(t,2.0);

  // set diag to 4.0
  t.SetTime(4.0);
  // check diag is 2.0, off diag is 0.0
  check_lambda(t,4.0);

  ASSERT_EQ(4.0, t.Time()); 
}

// testing TimeDepCoeff class
TEST(CoeffOpTestSuite, TimeDepCoeffTest)
{
  // make a mesh 
  MeshPtr_t my_mesh = make_mesh(0.0,4.0,5); 

  // some lambdas to test out 
  auto lam01 = [](double x){return x*x;}; // x^2
  auto lam02 = [](double x){return std::sin(x);}; // sin(x)
  auto lam03 = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // some polynomial in x 
  
  // take two vectors and check each |ui-vi| < eps 
  auto check_lamda = [](Eigen::VectorXd u, double val){
    double tol = 1e-4;  
    for(int i=0; i<u.size(); i++){
      ASSERT_NEAR(u[i],val,tol);
    }
  };

  // make coeff
  TimeDepCoeff coeff(lam01,my_mesh); 

  // check they have the same values 
  check_lamda(coeff.GetDiag(), lam01(coeff.Time())); 

  // set to new functions / discretization 
  coeff = lam02; 
  coeff.SetTime(2.0); 
  // check again
  check_lamda(coeff.GetDiag(), lam02(coeff.Time())); 

  // set to new functions / discretization 
  coeff = lam03; 
  // check again
  check_lamda(coeff.GetDiag(), lam03(coeff.Time())); 
}

// testing AutonomousCoeff class
TEST(CoeffOpTestSuite, AutonomousCoeffTest)
{
  // make a mesh 
  MeshPtr_t my_mesh = make_mesh(0.0,4.0,5); 

  // some lambdas to test out 
  auto lam01 = [](double x){return x*x;}; // x^2
  auto lam02 = [](double x){return std::sin(x);}; // sin(x)
  auto lam03 = [](double x){return 2*x*x*x-5*x*x+3*x-1;}; // some polynomial in x 
  
  // take two vectors and check each |ui-vi| < eps 
  auto check_lamda = [](Eigen::VectorXd u, Eigen::VectorXd v){
    int s=u.size(), s2=v.size();
    ASSERT_EQ(s,s2);
    double tol = 1e-4;  
    for(int i=0; i<s; i++){
      ASSERT_NEAR(u[i],v[i],tol);
    }
  };

  // make discretization + coeff
  Discretization1D disc;
  disc.set_init(my_mesh,lam01);
  AutonomousCoeff coeff(lam01,my_mesh); 

  // check they have the same values 
  check_lamda(disc.values(), coeff.GetDiag()); 

  // set to new functions / discretization 
  disc.set_init(lam02); 
  coeff = lam02; 
  // check again
  check_lamda(disc.values(), coeff.GetDiag()); 

  // set to new functions / discretization 
  disc.set_init(lam03); 
  coeff = lam03; 
  // check again
  check_lamda(disc.values(), coeff.GetDiag()); 
}

// // NthDerivOp Tests =========================================================================== 

// // // testing Fornberg algorithm 
// // TEST(NthDerivOpSuite, FornbergAlgo)
// // {}

// // // testing stateful Fornberg calculator 
// // TEST(NthDerivOpSuite, FornbergCalc)
// // {}

// testing NthDerivOp is constructible 
TEST(NthDerivOpSuite, NthDerivOpConstructible)
{
  NthDerivOp my_deriv; 
}

TEST(NthDerivOpSuite, Method_set_mesh_completing)
{
  auto my_mesh_01 = make_mesh(0.0,10.0,11);
  auto my_mesh_02 = make_mesh(0.0,10.0,101);
  auto my_mesh_03 = make_mesh(0.0,10.0,1001);
  auto my_mesh_04 = make_mesh(0.0,10.0,10001);

  auto D1 = NthDerivOp(1); 
  auto D2 = NthDerivOp(1); 
  auto D3 = NthDerivOp(1); 
  auto D4 = NthDerivOp(1); 

  auto test_mesh_lam = [](auto mesh_ptr){
    auto D1 = NthDerivOp(1); 
    auto D2 = NthDerivOp(1); 
    auto D3 = NthDerivOp(1); 
    auto D4 = NthDerivOp(1); 
    D1.set_mesh(mesh_ptr); 
    D1.GetMat(); 
    D2.set_mesh(mesh_ptr); 
    D2.GetMat(); 
    D3.set_mesh(mesh_ptr); 
    D3.GetMat(); 
    D4.set_mesh(mesh_ptr);
    D1.GetMat(); 
  };

  test_mesh_lam(my_mesh_01);
  test_mesh_lam(my_mesh_02);
  test_mesh_lam(my_mesh_03);
  test_mesh_lam(my_mesh_04);
}

// testing NthDerivOp custom .compose() method  
TEST(NthDerivOpSuite, NthDerivOpCompose)
{
  using D = NthDerivOp; 
  // constructible with a certain order 
  ASSERT_EQ(D(3).Order(), 3); 

  // orders for derivative 
  std::size_t n=4, m=7; 

  // testing composition splits sums ((f+g)' = f' + g')
  auto add_expr = D(n) + D(m); 
  ASSERT_EQ(D(n).compose(add_expr).Lhs().Order(), n+n); 
  ASSERT_EQ(D(n).compose(add_expr).Rhs().Order(), n+m);
  
  // testing composition parses scalar multiply ((c*f)' = c*f')
  auto mult_expr = 4.0 * D(n); 
  // ASSERT_EQ(D(n).compose(mult_expr).Rhs().Order(), n+n);
}

// Plugin Tests ===========================================================================
// make sure functionality of standard linops still works 
TEST(FdmPluginSuite, StandardLinOps)
{
  // /*TEST(LinearOperatorSuite, IdentityConstructible)*/
  {
    // default construct uses nullptr
    IOp Identity01;
    ASSERT_EQ(Identity01.mesh(),nullptr);

    // construct with ptr arg
    auto my_mesh = make_mesh(); 
    IOp Identity02(my_mesh);
    ASSERT_EQ(Identity02.mesh(), my_mesh);   

    // Check that entries on diag are 1
    int s = my_mesh->size()-1; 
    ASSERT_EQ(Identity02.GetMat().coeff(0,0),1); 
    ASSERT_EQ(Identity02.GetMat().coeff(s,s),1); 
    ASSERT_EQ(Identity02.GetMat().coeff(s/2,s/2),1);

    // of diag are zero
    ASSERT_EQ(Identity02.GetMat().coeff(0,s),0); 
    ASSERT_EQ(Identity02.GetMat().coeff(s,0),0); 
  };

  // /*TEST(LinearOperatorSuite, RandLinOpConstructible)*/
  {
    // default construct uses nullptr
    RandLinOp Rand01;
    ASSERT_EQ(Rand01.mesh(),nullptr);

    // construct with ptr arg
    auto my_mesh = make_mesh(); 
    RandLinOp Rand02(my_mesh);
    ASSERT_EQ(Rand02.mesh(), my_mesh);   
  };

  // /*TEST(LinearOperatorSuite, RandLinOpGetMat)*/ 
  {
    // construct with ptr arg
    auto my_mesh = make_mesh(); 
    RandLinOp Rand01(my_mesh);

    ASSERT_EQ(Rand01.mesh(), my_mesh);   

    Eigen::MatrixXd result = Rand01.GetMat(); 
  };

  // /* TEST(LinearOperatorSuite, Composition)*/ 
  {
    // get a mesh
    auto my_mesh = make_mesh(); 

    // vector of [1,1,...,1]
    Discretization1D my_vals(my_mesh); 
    my_vals.set_init(1.0);

    // construct without ptr arg
    RandLinOp L1(my_mesh), L2(my_mesh);

    // composition L1( L2(.) ) 
    auto Expr = L1.compose(L2); 

    // get underlying Eigen::VectorXd results of .apply() 
    auto expression_result = Expr.apply(my_vals).values(); 
    auto manual_result = L1.apply(L2.apply(my_vals)).values(); 

    ASSERT_EQ(expression_result.size(), manual_result.size()); 
    for(int i=0; i< expression_result.size(); i++){
      ASSERT_NEAR(expression_result[i], manual_result[i], 1e-4);
    }
  };

  // /*TEST(LinearOperatorSuite, ExpressionChaining)*/
  {
    auto my_mesh = make_mesh(); 
    // just a messy expression 
    auto my_expr = (2.0*(2.0*(2.0*(2.0*IOp())))).compose(50*IOp(my_mesh) + RandLinOp() + IOp() - RandLinOp(my_mesh).compose(IOp(my_mesh)));

    // build it step by step too!
    auto tmp1 = 2.0*IOp();
    auto tmp2 = 2.0*tmp1; 
    auto tmp3 = 2.0*tmp2; 
    auto tmp4 = 2.0*tmp3;
    auto rhs1 = 50*IOp(my_mesh); 
    auto rhs2 = rhs1 + RandLinOp(); 
    auto rhs3 = rhs2 + IOp(); 
    auto rhs4 = rhs3 - RandLinOp(my_mesh).compose(IOp(my_mesh));

    auto my_expr2 = tmp4.compose(rhs4); // all temporaries still alive

    // not all lhs/rhs had a mesh in expression construction 
    my_expr2.set_mesh(my_mesh); 
    my_expr.set_mesh(my_mesh); 

    // we should be able to make into eigen Matrix no matter what
    Eigen::MatrixXd resulting_mat = my_expr.GetMat().eval(); 
  }

  // /* TEST(LinearOperatorSuite, Method_set_mesh_ExprHooking)*/ 
  {
    // Calling set_mesh on an expression E = L1 + L2 
    // should pass the mesh to L1.set_mesh() and L2.set_mesh() 
    
    // construct without mesh ptrs 
    IOp I_lval(nullptr);
    auto Expr = I_lval + IOp(nullptr);

    // make mesh and give it to expression
    auto my_mesh = make_mesh();
    Expr.set_mesh(my_mesh); 

    // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
    ASSERT_EQ(I_lval.mesh(), my_mesh);
    ASSERT_EQ(Expr.Rhs().mesh(), my_mesh);

    // test again for scalar multiply  // construct without mesh ptrs 
    IOp I2_lval(nullptr);
    double c=2.0; 
    auto Expr2 = 2.0* I2_lval;
    auto Expr3 = c*IOp(nullptr);
    Expr2.set_mesh(my_mesh); 
    Expr3.set_mesh(my_mesh); 
    // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
    ASSERT_EQ(I2_lval.mesh(), my_mesh);
    ASSERT_EQ(Expr2.Rhs().mesh(), my_mesh);
    ASSERT_EQ(Expr3.Rhs().mesh(), my_mesh);

    // test again for composition 
    RandLinOp I3_lval(nullptr); 
    IOp I4_lval(nullptr);
    auto Expr4 = I3_lval.compose(I3_lval); 
    auto Expr5 = I3_lval.compose(IOp(nullptr));
    Expr4.set_mesh(my_mesh); 
    Expr5.set_mesh(my_mesh); 
    // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
    ASSERT_EQ(I3_lval.mesh(), my_mesh);
    ASSERT_EQ(Expr4.Rhs().mesh(), my_mesh);
    ASSERT_EQ(Expr5.Rhs().mesh(), my_mesh);
  }
}

// testing out SetTime() hooking
TEST(FdmPluginSuite, Method_SetTime_Hooking)
{
  // make a mesh 
  MeshPtr_t my_mesh = make_mesh(0.0,4.0,5); 
  
  // make a TCoeff that acts on functions of mesh
  TCoeff t(my_mesh); 

  // lambda to check 4 corners + middle of a matrix expr
  auto check_lambda = [s = my_mesh->size()-1](const auto& expr, double val) -> void
  {
    MatrixStorage_t Mat = expr.GetMat(); 
    // Check that entries on diag are 2
    ASSERT_EQ(Mat.coeff(0,0),val); 
    ASSERT_EQ(Mat.coeff(s,s),val); 
    ASSERT_EQ(Mat.coeff(s/2,s/2),val);

    // of diag are zero
    ASSERT_EQ(Mat.coeff(0,s),0); 
    ASSERT_EQ(Mat.coeff(s,0),0); 
  };

  auto expr01 = t+t;
  auto expr02 = t*t; 
  auto expr03 = 3.0*t; 
  auto expr04 = t.compose(IOp()); 

  // set time, make sure it is equal 
  t.SetTime(1.0);
  ASSERT_EQ(1.0,t.Time()); 

  // set time of expression, make sure it propagates to t 
  expr01.SetTime(2.0);
  ASSERT_EQ(2.0,t.Time()); 

  // set time of expression, make sure it propagates to t 
  expr02.SetTime(3.0);
  ASSERT_EQ(3.0,t.Time()); 

  // set time of expression, make sure it propagates to t 
  expr03.SetTime(4.0);
  ASSERT_EQ(4.0,t.Time()); 

  // set time of expression, make sure it propagates to t 
  expr04.SetTime(5.0);
  ASSERT_EQ(5.0,t.Time()); 
}; 

// Taking 1 explicit step with a convection diffusion scheme 
TEST(FdmPluginSuite, Method_explicit_step_returning)
{
  auto r = 10.0; 
  int n_gridpoints = 10001;

  // uniform mesh from 0.0 to r with n_gridpoints
  MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints);

  Discretization1D my_vals;
  
  // Bump centered at r/2. Zero at 0.0 and r 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
  
  // fill my_vals with func(x)
  my_vals.set_init(my_mesh, func); 
  
  // create an fdm_scheme for convection diffusion equation 
  using D = NthDerivOp; 
  double dt = 0.01;
  // auto fdm_scheme = IOp() + (dt) * (-0.2*D(2) + 0.5*D(1));
  auto fdm_scheme = D(2);

  // set the boundary conditions to Dirichlet 0
  auto left = std::make_shared<DirichletBC>(0.0);
  auto right = std::make_shared<DirichletBC>(0.0); 
  fdm_scheme.lbc_ptr = left; 
  fdm_scheme.rbc_ptr = right; 

  // // fill out the stencil matrix
  fdm_scheme.set_mesh(my_mesh);

  // apply the stencil. this uses lbc_ptr->SetSolL & ...->...R 
  auto result = fdm_scheme.explicit_step(my_vals);

  ASSERT_EQ(result.at(0),left->boundary_val); 
  ASSERT_EQ(result.at(result.size()-1),right->boundary_val); 
}

// Taking 1 implicit step with a convection diffusion scheme 
TEST(FdmPluginSuite, Method_solve_implicit_returning)
{
  auto r = 10.0; 
  int n_gridpoints = 101;

  // uniform mesh from 0.0 to r with n_gridpoints
  MeshPtr_t my_mesh = make_mesh(0.0,r,n_gridpoints);

  Discretization1D my_vals;
  
  // Bump centered at r/2. Zero at 0.0 and r 
  auto func = [r, smush=10](double x){return std::pow(x*(r-x)*(4.0/(r*r)),smush);}; 
  
  // fill my_vals with func(x)
  my_vals.set_init(my_mesh, func); 
  
  // create an fdm_scheme for convection diffusion equation 
  using D = NthDerivOp; 
  double dt = 0.01;
  auto fdm_scheme = IOp() + (dt) * (-0.2*D(2) + 0.5*D(1));

  // set the boundary conditions to Dirichlet 0
  auto left = std::make_shared<DirichletBC>(0.0);
  auto right = left; 
  fdm_scheme.lbc_ptr = left; 
  fdm_scheme.rbc_ptr = right; 

  // // fill out the stencil matrix
  fdm_scheme.set_mesh(my_mesh);

  // // apply the stencil. this uses lbc_ptr->SetSolL & ...->...R 
  auto result = fdm_scheme.solve_implicit(my_vals);

  ASSERT_EQ(result.at(0),left->boundary_val); 
  ASSERT_EQ(result.at(result.size()-1),right->boundary_val); 
}




