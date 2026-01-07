// Traits.hpp
//
//
//
// JAF 12/11/2025

#ifndef FDMSTENCIL_TRAITS_H
#define   FDMSTENCIL_TRAITS_H

#include<type_traits>
#include "FdmPlugin.hpp"
#include "CoeffOps/CoeffOpBase.hpp"

// Given a type, see if it is derived from CoeffOpBase's crtp scheme ------------------------------------
template<typename T, typename = void> 
struct is_coeffop_crtp : std::false_type{}; 

template<typename T>
struct is_coeffop_crtp<T, std::void_t<typename std::remove_cv_t<std::remove_reference_t<T>>::is_coeff_flag>>: std::true_type{}; 

// given a type, see if it has a .sparseView() method ---------------------------------------------------
template<typename T, typename = void>
struct has_sparseview_method : std::false_type{}; 

template<typename T>
struct has_sparseview_method<T, std::void_t<decltype(std::declval<T>().sparseView())>> : std::true_type{}; 

// minimum # of doubles + return type of callable type F --------------------------------
template<typename F>
class callable_traits
{
  private: 
  // ----------------------------------------------------------
  // get return type of callable type G invoked on N doubles 
  template<typename G, std::size_t N_doubles, typename... Args> 
  struct result_traits
  {
    using result_type = typename result_traits<G,N_doubles-1, double, Args...>::result_type; 
  }; 

  template<typename G, typename... Args> 
  struct result_traits<G,0,Args...>
  {
    using result_type = typename std::invoke_result<G,Args...>::type; 
  }; 

  // ---------------------------------------------------------------
  // test if callable type G can be invoked on N doubles up to max_n_args
  static constexpr std::size_t max_N_doubles = 20; 
  template<typename G, std::size_t N_doubles=max_N_doubles, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
  struct arg_traits
  {
    constexpr static bool is_callable = std::is_invocable<G,Args...>::value;

    constexpr static std::size_t num_args(){ 
      if constexpr (is_callable){
        return sizeof...(Args); 
      }
      else{
        return arg_traits<G,N_doubles-1,double, Args...>::num_args(); 
      }
    }
    
    using result_type = typename result_traits<G,num_args()>::result_type; 
  }; 

  // terminating case 
  template<typename G, typename... Args> //  typename = std::enable_if<N!= std::size_t{-1}> 
  struct arg_traits<G,0,Args...>
  {
    constexpr static bool is_callable = std::is_invocable<G,Args...>::value; 

    constexpr static std::size_t num_args(){ 
      if constexpr (is_callable){
        return sizeof...(Args); 
      }
      else{
        static_assert(false, "maximum length of args reached"); 
      }
    } 

    using result_type = typename result_traits<F,num_args()>::result_type; 
  }; 

  public:
  constexpr static std::size_t num_args = arg_traits<F>::num_args(); 
  using result_type = typename result_traits<F, num_args>::result_type; 
}; // end callable_traits<F> 

#endif