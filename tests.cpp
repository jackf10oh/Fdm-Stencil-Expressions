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
// #include<eigen3/unsupported/Eigen/KroneckerProduct>
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

// #include "LinOpsXD/MeshXD.hpp"
// #include "LinOpsXD/LinOpXDTraits.hpp"
#include "DiffOps/Utilities/make_diagonal.hpp"

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

  // create mesh 
  auto my_mesh = make_mesh(); 

  // create flat 1 row matrix 
  // MatrixStorage_t A(1,my_mesh->size());

  // cout << "------------" << endl; 
  // cout << A << endl; 


  // // use a BC to set left/right end of A 
  // BcPtr_t bc = std::make_shared<NeumannBC>(); 

  // // set left side of row 
  // bc->SetStencilL(A, my_mesh); 

  // // store into a temp vector 
  // std::vector<double> temp(A.valuePtr(), A.valuePtr()+A.nonZeros());

  // // set right side of row. this zeros left side 
  // bc->SetStencilR(A, my_mesh);
  // // repopulate left side 
  // for(auto i=0; i<temp.size(); i++) A.coeffRef(0,i) = temp[i]; 

  // cout << "------------" << endl; 
  // cout << A << endl; 

  Eigen::VectorXd vec(4);
  vec << 1, 2, 4, 8;
  Eigen::MatrixXd mat;
  mat = makeCirculant(vec);
  std::cout << mat << std::endl;
};

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


// printed output 
/*
0, 0, 0, 
1, 0, 0, 
2, 0, 0, 
0, 1, 0, 
1, 1, 0, 
2, 1, 0, 
0, 2, 0, 
1, 2, 0, 
2, 2, 0, 
0, 3, 0, 
1, 3, 0, 
2, 3, 0, 
0, 4, 0, 
1, 4, 0, 
2, 4, 0, 
0, 0, 1, 
1, 0, 1, 
2, 0, 1, 
0, 1, 1, 
1, 1, 1, 
2, 1, 1, 
0, 2, 1, 
1, 2, 1, 
2, 2, 1, 
0, 3, 1, 
1, 3, 1, 
2, 3, 1, 
0, 4, 1, 
1, 4, 1, 
2, 4, 1, 
0, 0, 2, 
1, 0, 2, 
2, 0, 2, 
0, 1, 2, 
1, 1, 2, 
2, 1, 2, 
0, 2, 2, 
1, 2, 2, 
2, 2, 2, 
0, 3, 2, 
1, 3, 2, 
2, 3, 2, 
0, 4, 2, 
1, 4, 2, 
2, 4, 2, 
0, 0, 3, 
1, 0, 3, 
2, 0, 3, 
0, 1, 3, 
1, 1, 3, 
2, 1, 3, 
0, 2, 3, 
1, 2, 3, 
2, 2, 3, 
0, 3, 3, 
1, 3, 3, 
2, 3, 3, 
0, 4, 3, 
1, 4, 3, 
2, 4, 3, 
0, 0, 4, 
1, 0, 4, 
2, 0, 4, 
0, 1, 4, 
1, 1, 4, 
2, 1, 4, 
0, 2, 4, 
1, 2, 4, 
2, 2, 4, 
0, 3, 4, 
1, 3, 4, 
2, 3, 4, 
0, 4, 4, 
1, 4, 4, 
2, 4, 4, 
0, 0, 5, 
1, 0, 5, 
2, 0, 5, 
0, 1, 5, 
1, 1, 5, 
2, 1, 5, 
0, 2, 5, 
1, 2, 5, 
2, 2, 5, 
0, 3, 5, 
1, 3, 5, 
2, 3, 5, 
0, 4, 5, 
1, 4, 5, 
2, 4, 5, 
0, 0, 6, 
1, 0, 6, 
2, 0, 6, 
0, 1, 6, 
1, 1, 6, 
2, 1, 6, 
0, 2, 6, 
1, 2, 6, 
2, 2, 6, 
0, 3, 6, 
1, 3, 6, 
2, 3, 6, 
0, 4, 6, 
1, 4, 6, 
2, 4, 6, 
*/
