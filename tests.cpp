// main.cpp
//
//
//
// JAF 12/8/2025

#include<iostream>
#include<iomanip>
#include<vector>
#include<eigen3/Eigen/Dense>
#include<eigen3/unsupported/Eigen/KroneckerProduct>
#include "DiffOps/All.hpp"
#include "LinOps/All.hpp"

#include "LinOpsXD/MeshXD.hpp"
#include "LinOpsXD/LinOpXDTraits.hpp"

typedef double Real;
using std::cout, std::endl;


template<typename F, typename Tup_t, typename = void>
struct has_apply_result_double_impl : public std::false_type{}; 

template<typename F,typename Tup_t>
struct has_apply_result_double_impl<
  F,  
  Tup_t,
  std::void_t<decltype(std::apply(std::declval<F>(),std::declval<Tup_t>()))>> 
  : public std::is_same<double, decltype(std::apply(std::declval<F>(),std::declval<Tup_t>()))>
{}; 

template<typename F, typename Tup_t> 
using has_apply_result_double = has_apply_result_double_impl<F,Tup_t>; 

int main()
{
  // given a function F(double ...) unwrap a vector and apply it 
  auto lam_03 = [](double x, double y, double z){return x*x + y*y + z*z;}; 

  std::tuple<double,double,double> my_point{ 1.0,1.0,1.0 }; 
  std::tuple<double,double,double> my_point2{ 1.0,1.0, 0.0 }; 
  std::tuple<double,double,double, double> my_point3{ 1.0,1.0, 0.0, 1.0}; 

  cout << std::apply(lam_03, my_point) << endl; 
  // cout << std::apply(lam_03, my_point3) << endl; 

  cout << has_apply_result_double<decltype(lam_03),decltype(my_point)>::value << endl; 
  // cout << has_apply_result_double<decltype(lam_03),decltype(my_point3)>::value << endl; 

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
