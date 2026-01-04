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
  MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 3, 3);

  DiscretizationXD my_vals(my_meshes); 
  my_vals.set_init(1.5); 
  
  std::vector<MemView<StrideIterator<Eigen::VectorXd::iterator>>> result; 
  result.reserve()

  // given an ith dimension from meshes 
  std::size_t ith_dim = 0; 

  // there are (sizes_product) / (ith dim size) many copies that look like 1D meshes of it 
  std::size_t num_copies = my_meshes->sizes_product() / my_meshes->dim_size(ith_dim); 

  // iterate through the copies 
  for(std::size_t n=0; n<num_copies; n++)
  {
    std::size_t mod = my_meshes->sizes_middle_product(0, ith_dim); 
    std::size_t scale = mod * my_meshes->dim_size(ith_dim); 
    std::size_t offset = (mod ? n % mod : n) + (scale * (n/mod));  
    // cout << "offset: " << offset << " stride: " << mod << endl;

    StrideIterator begin (my_vals.values().begin()+offset, mod);  
    auto end = begin + my_meshes->dim_size(ith_dim); 

    MemView<StrideIterator<Eigen::VectorXd::iterator>> view(begin, end);
       
    for(auto it=view.begin(); it!= view.end(); it++){
      *it = n; 
    }
  }

  cout << my_vals.values() << endl; 
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