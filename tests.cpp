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

  // vector of pair of BcPtrs. what BoundaryCondXD will look like. 
  std::vector<std::pair<BcPtr_t,BcPtr_t>> bc_list; 
  bc_list.push_back({make_neumann(0.0), make_neumann(1.0)});
  bc_list.push_back({make_dirichlet(0.0), make_dirichlet(1.0)});

  std::size_t s = my_meshes->GetMesh(0)->size();
  MatrixStorage_t A(s,s);
  bc_list[0].first->SetStencilL(A,my_meshes->GetMesh(0));
  bc_list[0].first->SetStencilR(A,my_meshes->GetMesh(0));

  cout << A << endl; 

  IOp I; 
  I.set_mesh(my_meshes->GetMesh(1)); 
  cout << Eigen::KroneckerProductSparse(I.GetMat(), A) << endl; 
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
