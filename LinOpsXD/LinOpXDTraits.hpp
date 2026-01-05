// LinOpXDTraits.hpp
//
//
//
// JAF 12/29/25 

#ifndef LINOPXDTRAITS_H
#define LINOPXDTRAITS_H

#include<tuple>
#include<functional>
#include<type_traits> 

template<typename F, std::size_t N_doubles, typename... Args> 
struct result_traits
{
  using result_type = typename result_traits<F,N_doubles-1, double, Args...>::result_type; 
}; 

template<typename F, typename... Args> 
struct result_traits<F,0,Args...>
{
  using result_type = typename std::invoke_result<F,Args...>::type; 
}; 

// test if callable type F can be invoked on N doubles up to max_n_args
static constexpr std::size_t max_N_doubles = 20; 
template<typename F, std::size_t N_doubles=max_N_doubles, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
struct arg_traits
{
  constexpr static bool is_callable = std::is_invocable<F,Args...>::value;

  constexpr static std::size_t num_args(){ 
    if constexpr (is_callable){
      return sizeof...(Args); 
    }
     else{
      return arg_traits<F,N_doubles-1,double, Args...>::num_args(); 
    }
  }
  
  using result_type = typename result_traits<F,num_args()>::result_type; 
}; 

// terminating case 
template<typename F, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
struct arg_traits<F,0,Args...>
{
  constexpr static bool is_callable = std::is_invocable<F,Args...>::value; 

  constexpr static std::size_t num_args(){ 
    if constexpr (is_callable){
      return sizeof...(Args); 
    }
     else{
      static_assert(false, "maximum length of args reached"); 
      // return callable_traits<F,N-1,double, Args...>::num_args(); 
    }
  } 

  using result_type = typename result_traits<F,num_args()>::result_type; 
}; 

template<typename F>
struct callable_traits
{
  constexpr static std::size_t num_args = arg_traits<F>::num_args(); 
  using result_type = typename result_traits<F, num_args>::result_type; 
}; 

#endif // LinOpXDTraits.hpp