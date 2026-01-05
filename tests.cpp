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

#include "DiffOps/All.hpp" // must include first for plugin to take effect over linops?
#include "LinOps/All.hpp" 
#include "LinOpsXD/All.hpp"
#include "Utilities/PrintVec.hpp"

using std::cout, std::endl;

auto lam00 = [](){return 1.0;}; 
auto lam01 = [](double x){return std::sqrt(x*x);}; 
auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

int main()
{
  // // iomanip 
  std::cout << std::setprecision(2); 

  // mesh assembly. 2 dims 
  MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 7, 2);

  // discretization.
  DiscretizationXD my_vals; 

  my_vals.set_init(my_meshes, lam02); 

  // boundary conditions 
  // auto bc = make_dirichlet(4.0); 
  
  // auto views = my_vals.OneDim_views(0); 
  // for(auto& v : views){
  //   bc->SetSolL(v, my_meshes->GetMesh(0)); 
  //   bc->SetSolR(v, my_meshes->GetMesh(0)); 
  // }

  // print as flat vector 
  // cout << my_vals.values() << endl; 

  // print as a 2d matrix 
  cout << Eigen::Map<Eigen::MatrixXd>(my_vals.values().data(), my_meshes->dim_size(0), my_meshes->dim_size(1)) << endl; 

  // from a given "slice" that looks like 2d. print the matrix 
  // std::size_t ith_slice = 20; 
  // std::size_t offset = ith_slice * my_meshes->sizes_middle_product(0,2); 
  // cout << Eigen::Map<Eigen::MatrixXd>(my_vals.values().data() + offset, my_meshes->dim_size(0), my_meshes->dim_size(1)) << endl; 

};

  // // how to iterate through dynamic multi dim meshes?
  // std::size_t end = my_meshes->sizes_product();  
  // std::vector<double> coords(my_meshes->dims()); 
  // for(std::size_t flat_i=0; flat_i<end; flat_i++){
  //   for(std::size_t dim=0; dim<my_meshes->dims(); dim++){
  //     std::size_t dim_i = flat_i;
  //     dim_i /= my_meshes->sizes_middle_product(0, dim); 
  //     dim_i %= my_meshes->dim_size(dim); 
  //     cout << dim_i; 
  //     cout << ", ";
  //   }
  //   cout << endl; 
  // } 

// 1st dim offset : 4*n 
// 2nd dim offset : (n / my_meshes->sizes_middle_product(0,ith_dim))

/*
0 -> 0,4,8,12 
1 -> 1,5,9,13
2 -> 2,6,10,14 
3 -> 3,7,11,15  

4 -> 0+16, 4+16, 8+16, 12+16
5 -> 1+16,5+16,9+16,13+16
6 -> 2+16,6+16,10+16,14+16 
7 -> 3+16,7+16,11+16,15+16 

*/
// printed output 
/*
0, 0, 0, 
1, 0, 0, 
2, 0, 0, 
3, 0, 0, 
0, 1, 0, 
1, 1, 0, 
2, 1, 0, 
3, 1, 0, 
0, 2, 0, 
1, 2, 0, 
2, 2, 0, 
3, 2, 0, 
0, 3, 0, 
1, 3, 0, 
2, 3, 0, 
3, 3, 0, 
0, 0, 1, 
1, 0, 1, 
2, 0, 1, 
3, 0, 1, 
0, 1, 1, 
1, 1, 1, 
2, 1, 1, 
3, 1, 1, 
0, 2, 1, 
1, 2, 1, 
2, 2, 1, 
3, 2, 1, 
0, 3, 1, 
1, 3, 1, 
2, 3, 1, 
3, 3, 1, 
0, 0, 2, 
1, 0, 2, 
2, 0, 2, 
3, 0, 2, 
0, 1, 2, 
1, 1, 2, 
2, 1, 2, 
3, 1, 2, 
0, 2, 2, 
1, 2, 2, 
2, 2, 2, 
3, 2, 2, 
0, 3, 2, 
1, 3, 2, 
2, 3, 2, 
3, 3, 2, 
0, 0, 3, 
1, 0, 3, 
2, 0, 3, 
3, 0, 3, 
0, 1, 3, 
1, 1, 3, 
2, 1, 3, 
3, 1, 3, 
0, 2, 3, 
1, 2, 3, 
2, 2, 3, 
3, 2, 3, 
0, 3, 3, 
1, 3, 3, 
2, 3, 3, 
3, 3, 3, 
*/ 