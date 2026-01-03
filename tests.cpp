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
#include<unsupported/Eigen/KroneckerProduct>
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

#include "LinOpsXD/MeshXD.hpp"
#include "LinOpsXD/BoundaryCondXD.hpp" 
#include "DiffOps/Utilities/FillStencil.hpp"
// #include "LinOpsXD/LinOpXDTraits.hpp"

using std::cout, std::endl;

auto random_banded = [](std::size_t N, int bands=1) -> MatrixStorage_t {
  if(bands<=0) throw std::invalid_argument("bands must be >= 1"); 
  auto rand = Eigen::MatrixXd::Random(N,N);
  Eigen::MatrixXd temp(rand.rows(), rand.cols()); 
  temp.setZero();
  for(int i=-bands; i<= bands; i++){
    temp.diagonal(i) = rand.diagonal(i); 
  }
  MatrixStorage_t A = temp.sparseView(); 
  return A; 
};

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  // mesh assembly. 2 dims 
  MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 5, 2); 

  // BoundaryCondXD. holds vector of pair of BcPtrs. 
  BoundaryCondXD bc_list; 
  bc_list.m_bc_list.push_back({make_neumann(0.0), make_neumann(1.0)});
  bc_list.m_bc_list.push_back({make_robin(1.0,1.0,5.0), make_robin(1.0,1.0,5.0)});
  // bc_list.m_bc_list.push_back({make_dirichlet(0.0), make_dirichlet(1.0)});

  std::size_t prod = my_meshes->sizes_product(); 
  MatrixStorage_t bc_mask(prod,prod); 

  // create an empty stencil A and set its boundarys for 1 dimensional case 
  std::size_t s = my_meshes->sizes_product();
  MatrixStorage_t A(s,s);
  bc_list.SetStencilImp(A, my_meshes); 
  cout << A << endl; 

  // // create an empty stencil A and set its boundarys for 1 dimensional case 
  // std::size_t s1 = my_meshes->GetMesh(0)->size();
  // MatrixStorage_t A(s1,s1);
  // bc_list.SetStencilImp(A, my_meshes); 
  // cout << A << endl; 

  // make a mask of this dimensions bcs 
  // std::size_t s2 = my_meshes->GetMesh(1)->size();
  // auto mat_row = flat_stencil(bc_list[1], my_meshes->GetMesh(1)); 
  // auto diag_expr = make_SparseDiag(mat_row);  
  // cout << Eigen::KroneckerProductSparse(diag_expr,I) << endl; 
};

// auto random_banded = [](std::size_t N, int bands=1) -> MatrixStorage_t {
//   if(bands<=0) throw std::invalid_argument("bands must be >= 1"); 
//   auto rand = Eigen::MatrixXd::Random(N,N);
//   Eigen::MatrixXd temp(rand.rows(), rand.cols()); 
//   temp.setZero();
//   for(int i=-bands; i<= bands; i++){
//     temp.diagonal(i) = rand.diagonal(i); 
//   }
//   MatrixStorage_t A = temp.sparseView(); 
//   return A; 
// };
