// LinOpTraits.hpp
//
//
//
// JAF 12/7/2025

#ifndef LINOPTRAITS_H
#define LINOPTRAITS_H

#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp"

// Utility funcs --------------------------------------------- 
template <typename T>
constexpr bool is_lvalue(T&&) {
  return std::is_lvalue_reference<T>{};
}

// Forward Declarations -------------------------------------------------
template<typename Derived>
class LinOpBase;

template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr; 

// flag for composition. i.e. L1( L2( . ) )
struct OperatorComposition_t{}; 

// flag type for composition. i.e. L1( L2( . ) )
struct OperatorAddition_t{}; 

// flag for scalar multiplication. i.e c*L
struct ScalarMultiply_t{}; 

// Traits ---------------------------------------------------------------
// given a type, detect if it is derived from linopbase<> 
template<typename T, typename = void>
struct is_linop_crtp_impl : std::false_type {};

template<typename T>
struct is_linop_crtp_impl<T, std::void_t<typename T::is_linop_tag>> : std::true_type {};

template<typename T>
using is_linop_crtp = is_linop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 

// given a linop expression. detect if it is a composition L1( L2( . ))
template<typename T>
struct is_compose_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_compose_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorComposition_t, OP>
{};

// given a linop expression. detect if it is a binary addition L1+L2
template<typename T>
struct is_add_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_add_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<OperatorAddition_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is a scalar multiply c*L
template<typename T>
struct is_scalar_multiply_expr: public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_scalar_multiply_expr<LinOpExpr<L,R,OP>>: public std::is_base_of<ScalarMultiply_t, std::remove_reference_t<OP>>
{};

// given a base class and flags, attach flags to base class
template<typename Base, typename... Flags>
struct make_flagged
{
  // the new class
  struct result : public Base, public Flags...
  {
    // only need to take constructor from base
    result(Base b) : Base(b){};
  }; 
  using type = result;
};

template<typename Base, typename... Flags>
using make_flagged_t = typename make_flagged<Base, Flags...>::type; 

// given a type T we may need to store it as a reference or a value in a binary expression
template<typename T>
struct Storage_t
{
  using type = typename std::conditional<
  is_linop_crtp<T>::value,  // if T derives from linopbase
  std::conditional_t<
    std::is_lvalue_reference<T>::value, // if T is an lvalue 
    T, // store by reference
    std::remove_reference_t<T> // else store rvalue by value
  >, 
  T // else store Non linops as is 
  >::type; 
}; 

// given a type T see if it is an expression
template<typename T>
struct is_expr_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_expr_crtp_impl<LinOpExpr<L,R,OP>>: public std::true_type
{};

template<typename T>
struct is_expr_crtp 
  : is_expr_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>> {};

// given a type T that will be a mixin, see if it has .apply() const method 
template<typename T, typename = void>
struct has_apply : public std::false_type {};

template<typename T>
struct has_apply<T, std::void_t<decltype(std::declval<const T>().apply(std::declval<const Discretization1D&>()))>> : public std::true_type{};

#endif // LinOpTraits.hpp
