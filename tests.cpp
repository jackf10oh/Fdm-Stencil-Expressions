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


#include "LinOpsXD/All.hpp"
// #include "LinOpsXD/MeshXD.hpp"
// #include "LinOpsXD/BoundaryCondXD.hpp" 
// #include "DiffOps/Utilities/FillStencil.hpp"
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
  
  // block random operator 
  DirectionalRandOp L(my_meshes,0); 

  // assembling BoundaryCondXD 
  BoundaryCondXD bc_list; 
  bc_list.m_bc_list.resize(0); 
  bc_list.m_bc_list.push_back({make_neumann(1.0),make_neumann(1.0)});
  bc_list.m_bc_list.push_back({make_neumann(1.0),make_neumann(1.0)});

  bc_list.SetStencilImp(L.GetMat(), L.mesh());


  Eigen::MatrixXd foo = L.GetMat(); 
  cout << foo << endl; 
};