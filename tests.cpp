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

// #include "LinOpsXD/All.hpp"
#include "LinOpsXD/MeshXD.hpp" 
#include "LinOpsXD/DiscretizationXD.hpp" 
#include "LinOpsXD/OperatorsXD/IOpXD.hpp" 
#include "LinOpsXD/OperatorsXD/DirectionalRandOp.hpp" 

#include "Utilities/StrideIterator.hpp"

using std::cout, std::endl;

int main()
{
  // // iomanip 
  std::cout << std::setprecision(2); 

  // mesh assembly. 2 dims 
  MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 4, 3);
  
  // given an ith dimension from meshes 
  std::size_t ith_dim = 1; 

  // there are (sizes_product) / (ith dim size) many copies that look like 1D meshes of it 
  std::size_t num_copies = my_meshes->sizes_product() / my_meshes->dim_size(ith_dim); 

  // iterate through the copies 
  for(std::size_t n=0; n<num_copies; n++)
  {
    std::size_t offset = n * my_meshes->sizes_middle_product(0, ith_dim+1); 
    // std::size_t step = 
    cout << "offset: " << offset << endl; 
    
  }

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
0 -> 0 
1 -> 1 
2 -> 2 
3 -> 3 
4 -> 16 + 1 

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