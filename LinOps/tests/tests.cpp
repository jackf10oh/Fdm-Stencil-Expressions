// tests.cpp
//
// LinOps Test suite with GTest
//
// JAF 12/5/2025

#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../All.hpp"

// Mesh Suite ---------------------------------------- 
TEST(MeshSuite1D, Mesh1DConstructible){
  // simply make a mesh and do nothing with it
  Mesh1D my_mesh; 

  // make a mesh with custom endpoints, # of steps
  Mesh1D my_custom_mesh(-10.0,20.0,31); 
  ASSERT_EQ(my_custom_mesh.size(), 31); 

  // mesh throws errors if # of grid points < 2
  EXPECT_ANY_THROW(Mesh1D(0.0, 1.0, 1));

  // mesh throws errors if x1 > x2 
  EXPECT_ANY_THROW(Mesh1D(1.0, -5.0, 6));
};

TEST(MeshSuite1D, Mesh1DIndexing){
  // make a mesh 
  Mesh1D my_mesh(0.0,10.0,11);  

  ASSERT_EQ(my_mesh[0], 0.0); 
  ASSERT_EQ(my_mesh[10], 10.0); 
  ASSERT_EQ(my_mesh.at(0), 0.0); 
  ASSERT_EQ(my_mesh.at(10), 10.0); 
  
  // operator[] doesn't check size 
  my_mesh[20]; // way out of bounds. but no error
  EXPECT_ANY_THROW(my_mesh.at(30)); 
};

TEST(MeshSuite1D, Mesh1DIterators){
  // simply make a mesh and do nothing with it
  int n_steps = 11;
  Mesh1D my_mesh(0.0,10.0,n_steps);

  // give all iterators as std::vec 
  my_mesh.begin(); 
  my_mesh.end(); 
  my_mesh.cbegin(); 
  my_mesh.cend(); 
  my_mesh.rbegin(); 
  my_mesh.rend(); 
  my_mesh.crbegin(); 
  my_mesh.crend(); 
  
  int count=0; 
  std::vector<double>::iterator it = my_mesh.begin(); 
  while(it!= my_mesh.end()){ count++; it++;}; 
  ASSERT_EQ(count,n_steps); 

  int reverse_count=0; 
  std::vector<double>::reverse_iterator rit = my_mesh.rbegin(); 
  while(rit!= my_mesh.rend()){ reverse_count++; rit++;}; 
  ASSERT_EQ(reverse_count,n_steps); 
};



// Discretization Suite ---------------------------------------- 
TEST(DiscretizationSuite1d, Disc1DConstructible)
{
  Discretization1D my_vals; 
}

TEST(DiscretizationSuite1d, Disc1DMovable)
{
  int n_steps=11; 
  auto my_mesh = make_mesh(0.0,10.0,n_steps); 
  Discretization1D moved_from(my_mesh);
  Discretization1D moved_to(std::move(moved_from));  
  ASSERT_TRUE(moved_from.mesh().expired()); 
  ASSERT_EQ(moved_from.values().data(),nullptr); // moved from now has invalid eigen::vectorxd 
  ASSERT_EQ(moved_to.size(), n_steps); 
  ASSERT_EQ(moved_to.mesh().lock(), my_mesh); 
}

TEST(DiscretizationSuite1d, Disc1DSetMesh)
{
  auto my_mesh = make_mesh(); 
  Discretization1D my_vals; 
  ASSERT_TRUE(my_vals.mesh().expired()); 
  Discretization1D discretization_w_stored_mesh(my_mesh); 
  ASSERT_FALSE(discretization_w_stored_mesh.mesh().expired());
}

TEST(DiscretizationSuite1d, Disc1DSetCosntant)
{
  auto my_mesh = make_mesh(); 
  Discretization1D my_vals(my_mesh); 

  double val_set = 0.0; 
  my_vals.set_init(val_set);
  
  ASSERT_EQ(my_vals.at(0), val_set); 
  ASSERT_EQ(my_vals.at(my_vals.size()-1), val_set); 
}

TEST(DiscretizationSuite1d, Disc1DSetByCallable)
{
  auto my_lambda = [](const double& x){return x*x;}; 
  struct Callable_t 
  {
    double operator()(const double& x){return x*x*x;}; 
  }; 
  Callable_t my_callable; 

  int n_steps = 101; 
  double left=0, right=100; 
  auto my_mesh = make_mesh(left,right,n_steps); 
  Discretization1D my_vals; 

  // with lambda 
  my_vals.set_init(my_mesh, my_lambda); 

  ASSERT_EQ(my_vals.at(0), my_lambda(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_lambda(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps-1), my_lambda(my_mesh->at(n_steps-1))); 

  my_vals.set_init(my_callable);

  ASSERT_EQ(my_vals.at(0), my_callable(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_callable(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps-1), my_callable(my_mesh->at(n_steps-1))); 

}

TEST(DiscretizationSuite1d, Disc1DIterators)
{
    // simply make a mesh and do nothing with it
  int n_steps = 11;
  auto my_mesh = make_mesh(0.0,10.0,n_steps);

  Discretization1D my_vals(my_mesh); 

  // give all iterators as std::vec 
  my_vals.begin(); 
  my_vals.end(); 
  my_vals.cbegin(); 
  my_vals.cend(); 
  // Eigen::VectorXd has no reverse iterators. 
  // my_vals.rbegin(); 
  // my_vals.rend(); 
  // my_vals.crbegin(); 
  // my_vals.crend(); 
  
  int count=0; 
  Eigen::VectorXd::iterator it = my_vals.begin(); 
  while(it!= my_vals.end()){ count++; it++;}; 
  ASSERT_EQ(count,n_steps); 
};



// Linear Operators Suite ---------------------------------------------
// Just the LinOpBase CRTP class. 
TEST(LinearOperatorSuite, IdentityConstructible)
{
  // default construct uses nullptr
  IOp Identity01;
  ASSERT_EQ(Identity01.mesh().lock(),nullptr);

  // construct with ptr arg
  auto my_mesh = make_mesh(); 
  IOp Identity02(my_mesh);
  ASSERT_EQ(Identity02.mesh().lock(), my_mesh);   

  // Check that entries on diag are 1
  int s = my_mesh->size()-1; 
  ASSERT_EQ(Identity02.GetMat().coeff(0,0),1); 
  ASSERT_EQ(Identity02.GetMat().coeff(s,s),1); 
  ASSERT_EQ(Identity02.GetMat().coeff(s/2,s/2),1);

  // of diag are zero
  ASSERT_EQ(Identity02.GetMat().coeff(0,s),0); 
  ASSERT_EQ(Identity02.GetMat().coeff(s,0),0); 
};

TEST(LinearOperatorSuite, RandLinOpConstructible)
{
  // default construct uses nullptr
  RandLinOp Rand01;
  ASSERT_EQ(Rand01.mesh().lock(),nullptr);

  // construct with ptr arg
  auto my_mesh = make_mesh(); 
  RandLinOp Rand02(my_mesh);
  ASSERT_EQ(Rand02.mesh().lock(), my_mesh);   
};

TEST(LinearOperatorSuite, RandLinOpGetMat)
{
  // construct with ptr arg
  auto my_mesh = make_mesh(); 
  RandLinOp Rand01(my_mesh);

  ASSERT_EQ(Rand01.mesh().lock(), my_mesh);   

  Eigen::MatrixXd result = Rand01.GetMat(); 
};

TEST(LinearOperatorSuite, RandLinOpApply)
{
  // setup mesh + discretization 
  auto my_mesh = make_mesh(); 

  Discretization1D my_vals;
  auto func = [](double x){return x*x;}; // x^2 
  my_vals.set_init(my_mesh, func); 

  // get a random linear operator 
  RandLinOp Rand01(my_mesh);
  // get its underlying matrix representation 
  Eigen::MatrixXd matrix_rep = Rand01.GetMat(); 

  // make sure .apply() gives the same as A*v 
  Discretization1D apply_method_result = Rand01.apply(my_vals); 
  Eigen::VectorXd manual_linalg_result = matrix_rep * my_vals.values(); 
  for(std::size_t i=0; i<apply_method_result.size(); i++){
    ASSERT_EQ(apply_method_result[i], manual_linalg_result[i]); 
  }
};

/* // Using LinOpExpr. 
TEST(LinearOperatorSuite, BasicAddition)
{
  auto my_mesh = make_mesh(); 
  // construct without ptr arg
  auto make_I_rval = [my_mesh](){return IOp(my_mesh);}; 
  IOp I_lval(my_mesh);

  auto sum01 = I_lval + I_lval;  
  auto sum02 = I_lval + make_I_rval();
  auto sum03 = make_I_rval() + I_lval;
  auto sum04 = make_I_rval() + make_I_rval(); 

  // lambda to check 4 corners + middle of a matrix expr
  auto check_lambda = [s = my_mesh->size()-1](const auto& expr) -> void
  {
    CustomStorage_t Mat = expr.GetMat(); 
    // Check that entries on diag are 2
    ASSERT_EQ(Mat.coeff(0,0),2); 
    ASSERT_EQ(Mat.coeff(s,s),2); 
    ASSERT_EQ(Mat.coeff(s/2,s/2),2);

    // of diag are zero
    ASSERT_EQ(Mat.coeff(0,s),0); 
    ASSERT_EQ(Mat.coeff(s,0),0); 
  };

  check_lambda(sum01); 
  check_lambda(sum02); 
  check_lambda(sum03); 
  check_lambda(sum04); 

  // eigen sparse matrices have no operator==() .... 
  auto to_dense = [](auto expr) -> Eigen::MatrixXd {
    Eigen::MatrixXd result;
    result = expr.GetMat();
    return result; 
  };

  // compar expr 1=2, 1=3, 1=4
  ASSERT_EQ(to_dense(sum01),to_dense(sum02));
  ASSERT_EQ(to_dense(sum01),to_dense(sum03));
  ASSERT_EQ(to_dense(sum01),to_dense(sum04));
};

TEST(LinearOperatorSuite, ScalarMultiplication)
{
  auto my_mesh = make_mesh(); 
  RandLinOp L1(my_mesh);

  // test multiplication with rvals
  2.0* RandLinOp(my_mesh);

  // multiply by scalar 
  double c = 3.0; 
  auto Expr = c * L1; 

  // evalues the expression object to a dense matrix
  auto to_dense = [](auto& expr) -> Eigen::MatrixXd {
    Eigen::MatrixXd result;
    result = expr.GetMat().eval();
    return result; 
  };

  // compar expr 1=2, 1=3, 1=4
  ASSERT_EQ(to_dense(Expr), (c*to_dense(L1)).eval());
};

TEST(LinearOperatorSuite, Composition)
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

TEST(LinearOperatorSuite, ExpressionChaining)
{
  auto my_mesh = make_mesh(); 
  // just a messy expression 
  auto my_expr = (2.0*(2.0*(2.0*(2.0*IOp())))).compose(50*IOp(my_mesh) + RandLinOp() + IOp() - RandLinOp(my_mesh).compose(IOp(my_mesh)));
  // not all lhs/rhs had a mesh in expression construction 
  my_expr.set_mesh(my_mesh); 
  // we should be able to make into eigen Matrix no matter what
  Eigen::MatrixXd resulting_mat = my_expr.GetMat(); 

  // building the same expression step by step 
  auto tmp1 = 2.0*IOp();
  auto tmp2 = 2.0*tmp1; 
  auto tmp3 = 2.0*tmp2; 
  auto tmp4 = 2.0*tmp3;
  auto rhs1 = 50*IOp(my_mesh); 
  auto rhs2 = rhs1 + RandLinOp(); 
  auto rhs3 = rhs2 + IOp(); 
  auto rhs4 = rhs3 - RandLinOp(my_mesh).compose(IOp(my_mesh));

  auto my_expr2 = tmp4.compose(rhs4); // all temporaries still alive

  // test the expression still has a .apply() method 
  Discretization1D disc; 
  disc.set_init(my_mesh, 1.0); 

  my_expr2.set_mesh(my_mesh); 
  Discretization1D result = my_expr2.apply(disc); 
  
}

TEST(LinearOperatorSuite, Method_set_mesh_ExprHooking)
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
  auto Expr5 = Expr4.compose(IOp(nullptr));
  // Expr4.set_mesh(my_mesh); 
  Expr5.set_mesh(my_mesh); 
  // both Lhs and Rhs should now have m_mesh_ptr == my_mesh
  ASSERT_EQ(I3_lval.mesh(), my_mesh);
  // ASSERT_EQ(Expr4.Rhs().mesh(), my_mesh);
  ASSERT_EQ(Expr5.Rhs().mesh(), my_mesh);
}

TEST(LinearOperatorSuite, LinOpTraits)
{
  struct foo {
    Discretization1D apply(const Discretization1D& d) const { return d; };
  };        
  struct bar{};

  // given a potential crtp mixin. see if it has .apply(discretization1d) -> const discretization1d method
  ASSERT_TRUE(has_apply<foo>::value); 
  ASSERT_FALSE(has_apply<bar>::value); 

  // see if a given type is a derived from the LinOpBase<> CRTP class
  ASSERT_TRUE(is_linop_crtp<RandLinOp>::value);
  ASSERT_FALSE(is_linop_crtp<int>::value); 
}

*/ 
