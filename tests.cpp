// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<Eigen/Dense>
#include<unsupported/Eigen/KroneckerProduct>
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

// #include "LinOpsXD/MeshXD.hpp"
// #include "LinOpsXD/LinOpXDTraits.hpp"
#include "DiffOps/Utilities/CirculantExpr.hpp"
#include "DiffOps/Utilities/SparseDiagExpr.hpp"

using std::cout, std::endl;

auto rand_banded = [](std::size_t N, int bands=1){
  auto rand = Eigen::MatrixXd::Random(N,N);
  Eigen::MatrixXd temp(rand.rows(), rand.cols()); 
  temp.setZero();
  for(int i=-bands; i<=bands; i++){
    temp.diagonal(i) = rand.diagonal(i); 
  }
  MatrixStorage_t A = temp.sparseView(); 
  return A; 
}; 

int main()
{
  std::cout << std::setprecision(2); 

  // create mesh + boundary condition 
  auto my_mesh = make_mesh(0,4,6); 
  BcPtr_t bc = std::make_shared<NeumannBC>(); 

  // first deriv stencil 
  auto D = MatrixStorage_t(my_mesh->size(), my_mesh->size());
  bc->SetStencilL(D,my_mesh); 
  bc->SetStencilR(D,my_mesh); 

  // create flat 1 row matrix 
  MatrixStorage_t A(1,my_mesh->size());

  // use a BC to set left/right end of A 
  // set left side of row 
  bc->SetStencilL(A, my_mesh); 
  // store into a temp vector 
  std::vector<double> temp(A.valuePtr(), A.valuePtr()+A.nonZeros());
  // set right side of row. this zeros left side 
  bc->SetStencilR(A, my_mesh);
  // repopulate left side 
  for(auto i=0; i<temp.size(); i++) A.coeffRef(0,i) = temp[i]; 

  cout << "------------" << endl; 
  // cout << A << endl; 

  cout << "------------" << endl; 
  auto foo = make_SparseDiag(A); 
  cout << foo.eval() << endl; 

  cout << "------------" << endl; 
  MatrixStorage_t I(my_mesh->size(), my_mesh->size()); 
  I.setIdentity();  
  // put boundary conditions on outter dimension 
  cout <<  Eigen::KroneckerProductSparse(foo,I) << endl;

  cout << "------------" << endl; 
  // put the boundary conditions on inner dimension 
  cout << Eigen::KroneckerProductSparse(I,D) << endl; 

};

  // auto my_mesh = make_mesh(0.0,3.0,4); 
  // NthDerivOp D(my_mesh); 
  // cout << D.GetMat() << endl; 

  // cout << "---------------------" << endl; 
  // MatrixStorage_t I(5,5); 
  // I.setIdentity(); 
  // auto prod = Eigen::kroneckerProduct(I, D.GetMat()); 
  // cout << prod << endl; 

  //   cout << "---------------" << endl; 
  // auto start = A.outerIndexPtr(); 
  // for(auto i=0; i<A.outerSize(); i++){
  //   cout << start[i] << endl; 
  // }

  // cout << "---------------" << endl; 
  // auto inner_start = A.innerIndexPtr(); 
  // for(auto i=0; i<A.nonZeros(); i++){
  //   cout << inner_start[i] << endl; 
  // }

  // cout << "---------------" << endl; 
  // auto val_start = A.valuePtr(); 
  // for(auto i=0; i<A.nonZeros(); i++){
  //   cout << val_start[i] << endl; 
  // }

  // std::vector<std::pair<double,double>> axes = {{0.0,1.0},{0.0,2.0},{0.0,3.0}}; 
  // std::vector<std::size_t> nsteps = {3,5,7}; 
  // MeshXDPtr_t my_mesh = std::make_shared<MeshXD>(axes,nsteps); 

  // cout << "my_mesh product_size: " << my_mesh->sizes_product() << endl; 

  // // how to iterate through dynamic multi dim meshes?
  // std::size_t end = my_mesh->sizes_product();  
  // std::vector<double> coords(my_mesh->dims()); 
  // for(std::size_t flat_i=0; flat_i<end; flat_i++){
  //   for(std::size_t dim=0; dim<my_mesh->dims(); dim++){
  //     std::size_t dim_i = flat_i;
  //     dim_i /= my_mesh->sizes_partial_product(dim); 
  //     dim_i %= my_mesh->dim_size(dim); 
  //     cout << dim_i; 
  //     cout << ", ";
  //   }
  //   cout << endl; 
  // }