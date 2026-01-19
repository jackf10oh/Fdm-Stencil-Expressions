// LinOpXDTraits.hpp
//
//
//
// JAF 12/29/25 

#ifndef LINOPXDTRAITS_H
#define LINOPXDTRAITS_H

#include<type_traits> 
#include "../LinOps/LinOpTraits.hpp"

namespace LinOps{

// Forward Declarations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class DiscretizationXD; 

template<typename Derived>
class LinOpBaseXD;

template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExprXD; 

// Traits =====================================================================
// given a type, detect if it is derived from linopbase<>  - - - - - - - - - - 
namespace internal{
template<typename T, typename = void>
struct is_linopxd_crtp_impl : std::false_type {};

template<typename T>
struct is_linopxd_crtp_impl<T, std::void_t<typename T::is_linopxd_tag>> : std::true_type {};
} // end namespace internal 

namespace traits{
template<typename T>
using is_linopxd_crtp = internal::is_linopxd_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
}

// given a type T see if it is an expression - - - - - - - - - - - - - - - - - - - - 
namespace internal{
template<typename T>
struct is_exprxd_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_exprxd_crtp_impl<LinOpExprXD<L,R,OP>>: public std::true_type
{};
} // end namespace internal 

namespace traits{
template<typename T>
using is_exprxd_crtp = internal::is_exprxd_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>;
}

// (Different Binary Operator Detections) ===================================================
namespace traits{
// given a linop expression. detect if it is a composition L1( L2( . )) - - - - - - - - - - 
// compose<T> : public false_type{} in LinOpTraits...

template<typename L, typename R, typename OP>
struct is_compose_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<internal::OperatorComposition_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is a binary addition L1+L2 - - - - - - - - - - 
// is_add_expr<T>: public std::false_type in LinOpTraits... 

template<typename L, typename R, typename OP>
struct is_add_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<internal::OperatorAddition_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is a scalar multiply c*L - - - - - - - - - - 
// is_scalar_multiply_expr<T>: public std::false_type in LinOpTraits...

template<typename L, typename R, typename OP>
struct is_scalar_multiply_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<internal::ScalarMultiply_t, std::remove_reference_t<OP>>
{}; 

// given a linop expression. detect if it is unary negation  -L - - - - - - - - - - 
// is_negation_expr<T> : public std::false_type in LinOpTraits ...

template<typename L, typename R, typename OP>
struct is_negation_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<internal::OperatorNegation_t, std::remove_reference_t<OP>>
{};

// given a linop expression. detect if it is binary subtraction L1-L2 - - - - - - - - - - 
// is_subtraction_expr<T>: public std::false_type in LinOpTraits...

template<typename L, typename R, typename OP>
struct is_subtraction_expr<LinOpExprXD<L,R,OP>>: public std::is_base_of<internal::OperatorSubtraction_t, std::remove_reference_t<OP>>
{};

}

// given a type T we may need to store it  - - - - - - - - - - - - - - - - - - - - 
namespace traits{
// as a reference or a value in a binary expression
// template<typename T, typename = void> in LinOpTraits....
// struct Storage_t 
// {
//   using type = T; 
// }; 

template<typename LINOPXD_T>
struct Storage_t<LINOPXD_T, std::enable_if_t< is_linopxd_crtp<LINOPXD_T>::value > >
{
  using type = std::conditional_t<
    std::is_lvalue_reference<LINOPXD_T>::value && !(is_expr_crtp<LINOPXD_T>::value), // if T is an lvalue LinOpBaseXD + not expression 
    LINOPXD_T, // store by lvalue LINOP_T& 
    std::remove_reference_t<LINOPXD_T> // else store rvalue by value
  >; 
}; 

} // end namespace traits 

// given a type T that will be a mixin, see if it has .apply() const method - - - - - - - - - - - - - - - - - -
namespace traits{
template<typename T, typename = void>
struct has_applyxd : public std::false_type {};

template<typename T>
struct has_applyxd<T, std::void_t<decltype(std::declval<const T>().apply(std::declval<const DiscretizationXD&>()))>> : public std::true_type{};
} // end namespace traits

} // end namespace LinOps 

#endif // LinOpXDTraits.hpp