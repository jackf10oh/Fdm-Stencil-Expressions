// tests.cpp (LinOpsXD)
//
// Gtest suites for LinOpsXD
//
// JAF 1/5/2026 

#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include<Eigen/Sparse>
#include<gtest/gtest.h>
#include<gmock/gmock.h>

#include "../All.hpp"

// Mesh Suite ---------------------------------------- 
TEST(MeshSuiteXD, MeshXDConstructible){
  // simply make a mesh and do nothing with it
  MeshXD my_mesh; 

  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 
  
  // make a mesh with custom endpoints, # of steps, # of dimension 
  MeshXD from_ends_steps_dims(-10.0,20.0,31, 3); 

  // from std::vec of [left,right] and n_steps 
  MeshXD from_ends_list_steps_list01({{0.0,4.0}}, {8}); 
  MeshXD from_ends_list_steps_list02({{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}, {8,16,32}); 


  // copy 
  MeshXD from_copy(from_dims_only); 

  // mesh throws errors if # of grid points < 2
  EXPECT_ANY_THROW(MeshXD(0.0, 1.0, 1, 4));

  // mesh throws errors if x1 > x2 
  EXPECT_ANY_THROW(MeshXD(1.0, -5.0, 6, 1));
}; 

TEST(MeshSuiteXD, MeshXDMesh1DGetters){
  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 
  const MeshXD from_dims_only_const(X); 

  MeshPtr_t resulting_mesh1d = from_dims_only.GetMesh(0);
  // !!! not safe to index 
  from_dims_only.GetMesh(X+5); 
  EXPECT_ANY_THROW(resulting_mesh1d = from_dims_only.GetMeshAt(X+5));  

 
  resulting_mesh1d = from_dims_only.GetMeshAt(0); 
  // safe to index 
  EXPECT_ANY_THROW(from_dims_only.GetMeshAt(X+5));  
}; 

TEST(MeshSuiteXD, MeshXDSizesGetters){
  // simply make a mesh and do nothing with it
  MeshXD my_mesh; 

  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  MeshXD from_dims_only(X); 

  // make a mesh with custom endpoints, # of steps, # of dimension 
  std::size_t n_steps=31; 
  MeshXD from_ends_steps_dims(-10.0,20.0,n_steps, 3); 

  // from std::vec of [left,right] and n_steps 
  std::size_t X1{8}, X2{16}, X3{32}; 
  MeshXD from_ends_list_steps_list01({{0.0,4.0}}, {X1}); 
  MeshXD from_ends_list_steps_list02({{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}, {X1,X2,X3});

  // returns # of dimensions 
  ASSERT_EQ(from_dims_only.dims(), X); 
  
  // returns # of steps in a specific dimension 
  ASSERT_EQ(from_ends_steps_dims.dim_size(0), n_steps);
  ASSERT_EQ(from_ends_steps_dims.dim_size(1), n_steps);
  EXPECT_ANY_THROW(from_ends_steps_dims.dim_size(10)); // 10 > # of dims 

  // returns product of (# of steps per dim)
  ASSERT_EQ(from_ends_list_steps_list01.sizes_product(),X1); 
  ASSERT_EQ(from_ends_list_steps_list02.sizes_product(),X1*X2*X3);
  
  // returns product of (# of steps per dim) for dim [start,end)
  ASSERT_EQ(from_ends_list_steps_list02.sizes_middle_product(0,3), X1*X2*X3); 
  ASSERT_EQ(from_ends_list_steps_list02.sizes_middle_product(0,2), X1*X2); 
  ASSERT_EQ(from_ends_list_steps_list02.sizes_middle_product(1,3), X2*X3); 
  EXPECT_ANY_THROW(from_ends_list_steps_list02.sizes_middle_product(0,10)); // 10 > # of dims 
  EXPECT_ANY_THROW(from_ends_list_steps_list02.sizes_middle_product(4,0)); // start > end  
}; 

// Discretization Suite ---------------------------------------- 

// Linear Operators XD Suite ---------------------------------------------

