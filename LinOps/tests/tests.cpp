// tests.cpp
//
// Test suite with GTest
//
// JAF 12/5/2025

#include<vector>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>
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
  auto my_mesh = std::make_shared<Mesh1D>(0.0,10.0, n_steps); 
  Discretization1D moved_from(my_mesh);
  Discretization1D moved_to(std::move(moved_from));  
  ASSERT_EQ(moved_from.mesh(), nullptr); 
  ASSERT_EQ(moved_from.values().data(),nullptr); // moved from now has invalid eigen::vectorxd 
  ASSERT_EQ(moved_to.size(), n_steps); 
  ASSERT_EQ(moved_to.mesh(), my_mesh); 
}

TEST(DiscretizationSuite1d, Disc1DSetMesh)
{
  auto my_mesh = std::make_shared<Mesh1D>(); 
  Discretization1D my_vals; 
  ASSERT_EQ(my_vals.mesh(), nullptr); 
  Discretization1D discretization_w_stored_mesh(my_mesh); 
  ASSERT_FALSE(discretization_w_stored_mesh.mesh()==nullptr);
}

TEST(DiscretizationSuite1d, Disc1DSetCosntant)
{
  auto my_mesh = std::make_shared<Mesh1D>(); 
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
  auto my_mesh = std::make_shared<Mesh1D>(left,right,n_steps); 
  Discretization1D my_vals; 

  // with lambda 
  my_vals.set_init(my_mesh, my_lambda); 

  ASSERT_EQ(my_vals.at(0), my_lambda(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_lambda(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps), my_lambda(my_mesh->at(n_steps))); 

  my_vals.set_init(my_callable);

  ASSERT_EQ(my_vals.at(0), my_callable(my_mesh->at(0))); 
  ASSERT_EQ(my_vals.at(n_steps/2), my_callable(my_mesh->at(n_steps/2))); 
  ASSERT_EQ(my_vals.at(n_steps), my_callable(my_mesh->at(n_steps))); 

}

TEST(DiscretizationSuite1d, Disc1DIterators)
{
    // simply make a mesh and do nothing with it
  int n_steps = 11;
  auto my_mesh = std::make_shared<Mesh1D>(0.0,10.0,n_steps);

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

// Traits Suite ------------------------------------------------------- 
struct foo {
  Discretization1D apply(const Discretization1D& d) const { return d; };
};        

struct bar{};
// Linear Operators Suite ---------------------------------------------

// Expressions Suite ---------------------------------------------

