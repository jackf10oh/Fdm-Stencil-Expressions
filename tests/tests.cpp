// tests.cpp
//
// FdmStencilExpressions Test suite with GTest
//
// JAF 12/5/2025

#include<iostream>
#include<vector>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../incl/All.hpp"

// Testing TCoeff can be constructed
TEST(CoeffOpTestSuite, TCoeffConstructible)
{
  TCoeff t; 
}

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

// testing out SetTime() hooking
TEST(CoeffOpTestSuite, Method_SetTime_Hooking)
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

// Fornberg tests 







