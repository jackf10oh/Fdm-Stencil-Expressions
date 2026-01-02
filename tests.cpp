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

using std::cout, std::endl;

int main()
{
  std::cout << std::setprecision(2); 
  auto rand = Eigen::MatrixXd::Random(11,11);
  Eigen::MatrixXd temp(rand.rows(), rand.cols()); 
  temp.setZero();
  int n = 2; 
  for(int i=-n; i<= n; i++){
    temp.diagonal(i) = rand.diagonal(i); 
  }
  MatrixStorage_t A = temp.sparseView(); 

  // MatrixStorage_t A(rand.rows(), rand.cols()); 
  // rand.diagonal(0).asDiagonal().evalTo(A);   
  cout << A << endl; 
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
