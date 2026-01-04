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

using std::cout, std::endl;

// test if callable type F can be invoked on N doubles 
template<typename F, std::size_t N=20, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
struct callable_traits
{
  constexpr static bool is_callable = std::is_invocable<F,Args...>::value;

  constexpr static std::size_t num_args(){ 
    if constexpr (is_callable){
      return sizeof...(Args); 
    }
     else{
      return callable_traits<F,N-1,double, Args...>::num_args(); 
    }
  } 
}; 

template<typename F, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
struct callable_traits<F,0,Args...>
{
  constexpr static bool is_callable = std::is_invocable<F,Args...>::value; 

  using result_type = int; 

  constexpr static std::size_t num_args(){ 
    if constexpr (is_callable){
      return sizeof...(Args); 
    }
     else{
      static_assert(false, "maximum length of args reached"); 
      // return callable_traits<F,N-1,double, Args...>::num_args(); 
    }
  } 
}; 

template<typename F, std::size_t N, typename... Args> 
struct result_traits
{
  using result_type = typename result_traits<F,N-1, double, Args...>::result_type; 
}; 

template<typename F, typename... Args> 
struct result_traits<F,0,Args...>
{
  using result_type = typename std::invoke_result<F,Args...>::type; 
}; 

int main()
{
  // // iomanip 
  std::cout << std::setprecision(2); 

  // mesh assembly. 2 dims 
  MeshXDPtr_t my_meshes = std::make_shared<MeshXD>(0.0,1.0, 3, 3);

  // discretization.
  DiscretizationXD my_vals(my_meshes); 
  my_vals.set_init(0.0); 

  auto lam00 = [](){return 1.0;}; 
  auto lam01 = [](double x){return x*x + 1.0;}; 
  auto lam02 = [](double x, double y){return x*x + y*y;}; 
  auto lam03 = [](double x, double y, double z){return "foobar";}; 

  cout << callable_traits<decltype(lam00)>::num_args() << endl; 
  cout << callable_traits<decltype(lam01)>::num_args() << endl; 
  cout << callable_traits<decltype(lam02)>::num_args() << endl; 
  cout << callable_traits<decltype(lam03)>::num_args() << endl; 

  cout << typeid(result_traits<decltype(lam00),0>::result_type).name() << endl; 
  cout << typeid(result_traits<decltype(lam01),1>::result_type).name() << endl; 
  cout << typeid(result_traits<decltype(lam03),3>::result_type).name() << endl; 

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