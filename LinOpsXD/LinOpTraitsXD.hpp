// LinOpXDTraits.hpp
//
//
//
// JAF 12/29/25 

#ifndef LINOPXDTRAITS_H
#define LINOPXDTRAITS_H

#include<type_traits> 
#include "../LinOps/LinOpTraits.hpp"

// Forward declarations ------------------------------------------------- 
class DiscretizationXD; 

template<typename T>
class LinOpBaseXD; 
 
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExprXD; 

// Traits ==============================================================================
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

// given a type, detect if it is derived from linopbase<> --------------------------------
template<typename T, typename = void>
struct is_linopxd_crtp_impl : std::false_type {};

template<typename T>
struct is_linopxd_crtp_impl<T, std::void_t<typename T::is_linopxd_tag>> : std::true_type {};

template<typename T>
using is_linopxd_crtp = is_linopxd_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

// // given a type T see if it is an expression ---------------------------------------------
template<typename T>
struct is_exprxd_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_exprxd_crtp_impl<LinOpExprXD<L,R,OP>>: public std::true_type
{};

template<typename T>
using is_exprxd_crtp = is_exprxd_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>;

// // given a linop expression. detect if it is a composition L1( L2( . )) -----------------
// in LinOpTraits.hpp .... 
// template<typename T>
// struct is_compose_expr: public std::false_type{};

template<typename L, typename R, typename OP>
struct is_compose_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<OperatorComposition_t, std::remove_reference_t<OP>>
{};

// // given a linop expression. detect if it is a binary addition L1+L2 ----------------------
// in LinOpTraits.hpp .... 
// template<typename T>
// struct is_add_expr: public std::false_type{};

template<typename L, typename R, typename OP>
struct is_add_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<OperatorAddition_t, std::remove_reference_t<OP>>
{};

// // given a linop expression. detect if it is a scalar multiply c*L -------------------------
// in LinOpTraits.hpp .... 
// template<typename T>
// struct is_scalar_multiply_expr: public std::false_type{};

template<typename L, typename R, typename OP>
struct is_scalar_multiply_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<ScalarMultiply_t, std::remove_reference_t<OP>>
{};

// // given a base class and flags, attach flags to base class ----------------------------------
// in LinOpTraits.hpp ...... 
/*
template<typename Base, typename... Flags>
struct make_flagged{ ... };

template<typename Base, typename... Flags>
using make_flagged_t = typename make_flagged<Base, Flags...>::type; 
*/

// // given T we may need to store it as a reference or a value in a binary expression ----------
// in LinOpTraits.hpp ..... 
// template<typename T, typename = void>
// struct Storage_t {  using type = T; }; 

template<typename LINOPXD_T>
struct Storage_t<LINOPXD_T, std::enable_if_t< is_linopxd_crtp<LINOPXD_T>::value > >
{
  using type = std::conditional_t<
    std::is_lvalue_reference<LINOPXD_T>::value && !(is_exprxd_crtp<LINOPXD_T>::value), // if T is an lvalue LinOpBase + not expression 
    LINOPXD_T, // store by lvalue LINOP_T& 
    std::remove_reference_t<LINOPXD_T> // else store rvalue by value
  >; 
}; 

// // given a type T that will be a mixin, see if it has .apply() const method ---------------------
// in LinOpTraits.hpp ..... 
// template<typename T, typename = void>
// struct has_apply : public std::false_type {};

template<typename T>
struct has_apply<T, std::void_t<decltype(std::declval<const T>().apply(std::declval<const DiscretizationXD&>()))>> : public std::true_type{};

#endif // LinOpXDTraits.hpp