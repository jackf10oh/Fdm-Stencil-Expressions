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

// Mesh Suite ==================================================== 
// Test all constructors 
TEST(MeshXDSuite, MeshXDConstructible){
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

  // from std::vector<> of std::shared_ptr<Mesh1D>  
  auto mesh_1d_01 = make_mesh(); 
  auto mesh_1d_02 = make_mesh(-1.0,1.0,21); 
  auto mesh_1d_03 = make_mesh(0,10.0, 101); 
  MeshXD from_vec(std::vector<MeshPtr_t>({mesh_1d_01,mesh_1d_02,mesh_1d_03})); 

  // copy 
  MeshXD from_copy(from_dims_only); 

  // mesh throws errors if # of grid points < 2
  EXPECT_ANY_THROW(MeshXD(0.0, 1.0, 1, 4));

  // mesh throws errors if x1 > x2 
  EXPECT_ANY_THROW(MeshXD(1.0, -5.0, 6, 1));
}; 

// Testing all method that return std::shared_ptr<Mesh1D> 
TEST(MeshXDSuite, MeshXDMesh1DGetters){
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

// Testing all method that return sizes 
TEST(MeshXDSuite, MeshXDSizesGetters){
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

// Discretization Suite ============================================= 
// Test all constructors 
TEST(DiscretizationXDSuite, DiscretizationXDConstructible){
  // default constructor leaves m_vals, m_dims, m_mesh_ptr all empty 
  DiscretizationXD empty_disc; 

  // 1 dim from single size
  std::size_t N=67; 
  DiscretizationXD from_size_1d(N); 

  // from X dimensional MeshXDPtr_t
  MeshXDPtr_t meshes_1d = std::make_shared<MeshXD>(std::size_t{1}); 
  MeshXDPtr_t meshes_3d = std::make_shared<MeshXD>(std::size_t{3}); 
  DiscretizationXD from_meshes_1d(meshes_1d); 
  DiscretizationXD from_meshes_3d(meshes_3d); 

  // copy 
  DiscretizationXD from_other(from_meshes_1d); 

}

// testing move constructors ???????????????

// Testing all method that return sizes. should match original MeshXDPtr_t 
TEST(DiscretizationXDSuite, DiscretizationXDSizesGetters){
  // simply make a mesh and do nothing with it
  auto my_mesh = make_meshes(); 
  // make a unit mesh, 5 steps per dim, X dimensions 
  std::size_t X = 4; 
  auto from_dims_only = make_meshes(X); 
  // make a mesh with custom endpoints, # of steps, # of dimension 
  std::size_t n_steps=31; 
  auto from_ends_steps_dims = make_meshes(-10.0,20.0,n_steps, 3); 

  // from std::vec of [left,right] and n_steps 
  std::size_t X1{8}, X2{16}, X3{32}; 
  using bounds_vec_t = std::vector<std::pair<double,double>>; 
  bounds_vec_t bound_vals01{{0.0,4.0}}, bound_vals02{{0.0,4.0},{-2.0,2.0},{-4.0,0.0}}; 
  using nsteps_vec_t = std::vector<std::size_t>; 
  nsteps_vec_t nsteps01{X1}, nsteps02{X1, X2, X3}; 
  auto from_ends_list_steps_list01 = make_meshes(bound_vals01, nsteps01); 
  auto from_ends_list_steps_list02 = make_meshes(bound_vals02, nsteps02); 

  // size getters should give same result across DiscretizationXD / MeshXDPtr_t 
  // returns # of dimensions 
  ASSERT_EQ(from_dims_only->dims(), \
            DiscretizationXD(from_dims_only).dims()); 
  
  // returns # of steps in a specific dimension 
  ASSERT_EQ(from_ends_steps_dims->dim_size(0), \
            DiscretizationXD(from_ends_steps_dims).dim_size(0)); 
            
  ASSERT_EQ(from_ends_steps_dims->dim_size(1), \
            DiscretizationXD(from_ends_steps_dims).dim_size(1));  

  EXPECT_ANY_THROW(DiscretizationXD(from_ends_steps_dims).dim_size(10)); // 10 > # of dims 

  // // returns product of (# of steps per dim)
  ASSERT_EQ(from_ends_list_steps_list01->sizes_product(), \
            DiscretizationXD(from_ends_list_steps_list01).sizes_product());  

  ASSERT_EQ(from_ends_list_steps_list02->sizes_product(), \
            DiscretizationXD(from_ends_list_steps_list02).sizes_product());  
  
  // returns product of (# of steps per dim) for dim [start,end)
  ASSERT_EQ(from_ends_list_steps_list02->sizes_middle_product(0,3), \
            DiscretizationXD(from_ends_list_steps_list02).sizes_middle_product(0,3));  

  ASSERT_EQ(from_ends_list_steps_list02->sizes_middle_product(0,2), \
            DiscretizationXD(from_ends_list_steps_list02).sizes_middle_product(0,2));

  ASSERT_EQ(from_ends_list_steps_list02->sizes_middle_product(1,3), \
            DiscretizationXD(from_ends_list_steps_list02).sizes_middle_product(1,3));
 
  EXPECT_ANY_THROW(DiscretizationXD(from_ends_list_steps_list02).sizes_middle_product(0,10)); // 10 > # of dims 
  EXPECT_ANY_THROW(DiscretizationXD(from_ends_list_steps_list02).sizes_middle_product(4,0)); // start > end 
}

// Testing set_init() for constant ????????????????????
TEST(DiscretizationXDSuite, DiscretizationXDSetByConstant){

  MeshXDPtr_t my_meshes = make_meshes(); 
  DiscretizationXD my_disc(my_meshes); 
  
  double val_init = 6.7; 
  my_disc.set_init(val_init); 

  for(auto& val : my_disc.values()){
    ASSERT_EQ(val, val_init); 
  }
}; 

// Testing set_init() for callable  
TEST(DiscretizationXDSuite, DiscretizationXDSetByCallable){
  auto lam00 = [](){return 67.0;}; 
  // Euclidean norms in 1D, 2D, 3D 
  auto lam01 = [](double x){return std::sqrt(x*x);}; 
  auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
  auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

  auto my_mesh_1d = make_meshes(0.0,10.0,21, 1); 
  auto my_mesh_2d = make_meshes(0.0,10.0,21, 2); 
  auto my_mesh_3d = make_meshes(0.0,10.0,21, 3); 

  DiscretizationXD my_disc; 

  // set init 1D case 
  my_disc.set_init(my_mesh_1d, lam01); 
  for(std::size_t i=0; i<my_disc.dim_size(0); i++){
    double x = my_mesh_1d->GetMeshAt(0)->at(i); 
    double lam_val = lam01(x); 
    double disc_val = my_disc.values()[i];  

    ASSERT_NEAR(disc_val, lam_val, 1e-4); 
  }

  // set init 2D case 
  my_disc.set_init(my_mesh_2d, lam02); 
  for(std::size_t i=0; i<my_disc.dim_size(0); i++){

    double x = my_mesh_2d->GetMeshAt(0)->at(i); 

    for(std::size_t j=0; j<my_disc.dim_size(1); j++){

      double y = my_mesh_2d->GetMeshAt(1)->at(j); 

      double lam_val = lam02(x,y); 

      std::size_t flat_idx = i + j*my_disc.sizes_middle_product(0,1); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 

    }
  }

  // set init 3D case 
  my_disc.set_init(my_mesh_3d, lam03); 
  for(std::size_t i=0; i<my_disc.dim_size(0); i++){

    double x = my_mesh_3d->GetMeshAt(0)->at(i); 

    for(std::size_t j=0; j<my_disc.dim_size(1); j++){

      double y = my_mesh_3d->GetMeshAt(1)->at(j); 

      for(std::size_t k=0; k<my_disc.dim_size(1); k++){

      double z = my_mesh_3d->GetMeshAt(2)->at(k); 

      double lam_val = lam03(x,y,z); 

      std::size_t flat_idx = i + j*my_disc.sizes_middle_product(0,1) + k*my_disc.sizes_middle_product(0,2); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 
      }
    }
  }

  // last case: we can use lower dimensional function on higher dimension mesh
  // i.e. lam02 takes 2 dims but my_mesh3d has 3 dims
  my_disc.set_init(my_mesh_3d, lam02); 
  for(std::size_t i=0; i<my_disc.dim_size(0); i++){

    double x = my_mesh_3d->GetMeshAt(0)->at(i); 

    for(std::size_t j=0; j<my_disc.dim_size(1); j++){

      double y = my_mesh_3d->GetMeshAt(1)->at(j); 

      for(std::size_t k=0; k<my_disc.dim_size(1); k++){

      double z = my_mesh_3d->GetMeshAt(2)->at(k); 

      double lam_val = lam02(x,y); 

      std::size_t flat_idx = i + j*my_disc.sizes_middle_product(0,1) + k*my_disc.sizes_middle_product(0,2); 
      double disc_val = my_disc.values()[flat_idx];

      ASSERT_NEAR(disc_val, lam_val, 1e-4); 
      }
    }
  }
}; 

// testing move assignment from Eigen::VectorXd ?????????????????

// Linear Operators XD Suite ====================================================

