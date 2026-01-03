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

#include "DiffOps/All.hpp" // must include first for plugin to take effect over linops 

#include "LinOps/All.hpp"

#include "LinOpsXD/All.hpp"

#include "Utilities/StrideIterator.hpp"

using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  std::vector<double> v(10); 
  for(auto i=0; i<v.size(); ++i) v[i] = i; 
  print_vec(v, "v"); 

  auto stride_view = make_strided_MemView(v, 2, 3); 
  print_vec(stride_view, "strided v");

};

  // // mesh assembly. 2 dims 
  // MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 5, 2);
  
  // // block random operator 
  // DirectionalRandOp L(my_meshes,0); 

  // // assembling BoundaryCondXD 
  // BoundaryCondXD bc_list; 
  // bc_list.m_bc_list.resize(0); 
  // bc_list.m_bc_list.push_back({make_neumann(1.0),make_neumann(1.0)});
  // bc_list.m_bc_list.push_back({make_neumann(1.0),make_neumann(1.0)});

  // bc_list.SetStencilImp(L.GetMat(), L.mesh());


  // Eigen::MatrixXd foo = L.GetMat(); 
  // cout << foo << endl; 