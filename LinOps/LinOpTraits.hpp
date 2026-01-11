// LinOpTraits.hpp
//
// List of TMP functions for LinOps 
//
// JAF 12/7/2025

#ifndef LINOPTRAITS_H
#define LINOPTRAITS_H

#include<type_traits>

namespace LinOps{

// Forward Declarations -------------------------------------------------
class Discretization1D; 

template<typename Derived>
class LinOpBase;

template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr; 

namespace internal{

// Flags for binary + unary ops -----------------------------------------------------------------------
// flag for composition. i.e. L1( L2( . ) )
struct OperatorComposition_t{}; 

// flag type for composition. i.e. L1 + L2 
struct OperatorAddition_t{}; 

// flag type for composition. i.e. L1 - L2 
struct OperatorSubtraction_t{}; 

// flag for scalar multiplication. i.e c*L
struct ScalarMultiply_t{}; 

// flag type for negation. i.e. -L1
struct OperatorNegation_t{}; 

// Structs for binary operations f(L1,L2) to get matrix of expression -----------------------------------------------------------------------
// L1 + L2 
struct linop_bin_add_op : public internal::OperatorAddition_t
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return (A.GetMat()) + (B.GetMat()); }
}; 
// L1 - L2 
struct linop_bin_subtract_op : public internal::OperatorSubtraction_t
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return (A.GetMat()) - (B.GetMat()); }
}; 
// c * L
struct scalar_left_mult_op : public internal::ScalarMultiply_t
{
  template<typename L2>
  auto operator()(const double& c, const L2& B) const { return  c*(B.GetMat()); }
}; 
// -L 
struct unary_negate_op : public internal::OperatorNegation_t
{
  template<typename L1>
  auto operator()(const L1& B) const { return  -(B.GetMat()); }
}; 
// composition: L1( L2( . ) )
struct linopXlinop_mult_op : public internal::OperatorComposition_t
{
  template<typename L1, typename L2>
  auto operator()(const L1& A, const L2& B) const { return  (A.GetMat())*(B.GetMat()); }
}; 


// Traits =====================================================================
// given a type, detect if it is derived from linopbase<> ---------------------
template<typename T, typename = void>
struct is_linop_crtp_impl : std::false_type {};

template<typename T>
struct is_linop_crtp_impl<T, std::void_t<typename T::is_linop_tag>> : std::true_type {};

template<typename T>
using is_linop_crtp = is_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

// given a type T see if it is an expression ----------------------------------------
template<typename T>
struct is_expr_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_expr_crtp_impl<LinOpExpr<L,R,OP>>: public std::true_type
{};

template<typename T>
using is_expr_crtp = is_expr_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>;

// given a linop expression. detect if it is a composition L1( L2( . )) ---------------
template<typename T>
struct is_compose_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_compose_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorComposition_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is a binary addition L1+L2 ---------------------
template<typename T>
struct is_add_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_add_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorAddition_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is a scalar multiply c*L ---------------------
template<typename T>
struct is_scalar_multiply_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_scalar_multiply_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<ScalarMultiply_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is unary negation  -L ---------------------
template<typename T>
struct is_negation_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_negation_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorNegation_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is binary subtraction L1-L2 ---------------------
template<typename T>
struct is_subtraction_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_subtraction_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorSubtraction_t, std::remove_reference_t<OP>>
{};

// given a type T we may need to store it  ----------------------------------------
// as a reference or a value in a binary expression
template<typename T, typename = void>
struct Storage_t 
{
  using type = T; 
}; 

template<typename LINOP_T>
struct Storage_t<LINOP_T, std::enable_if_t< is_linop_crtp<LINOP_T>::value > >
{
  using type = std::conditional_t<
    std::is_lvalue_reference<LINOP_T>::value && !(is_expr_crtp<LINOP_T>::value), // if T is an lvalue LinOpBase + not expression 
    LINOP_T, // store by lvalue LINOP_T& 
    std::remove_reference_t<LINOP_T> // else store rvalue by value
  >; 
}; 

// given type T, see if instances have .left_scalar_mult_impl(double) -------------------------
template<typename T, typename = void>
struct supports_left_scalar_mult : public std::false_type{}; 

template<typename T>
struct supports_left_scalar_mult<T, 
  std::void_t<
    decltype(std::declval<T>().left_scalar_mult_impl( std::declval<double>() ))
  >
> : public std::true_type{}; 

// given a type T that will be a mixin, see if it has .apply() const method ----------------------------------
template<typename T, typename = void>
struct has_apply : public std::false_type {};

template<typename T>
struct has_apply<T, std::void_t<decltype(std::declval<const T>().apply(std::declval<const Discretization1D&>()))>> : public std::true_type{};

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

} // end namespace internal 
} // end namespace LinOps 

#endif // LinOpTraits.hpp
