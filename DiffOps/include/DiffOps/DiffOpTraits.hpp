// DiffOpTraits.hpp
//
// List of TMP functions for LinOps 
//
// JAF 2/5/2026 

#ifndef DIFFOPTRAITS_H
#define DIFFOPTRAITS_H

#include<type_traits>

namespace DiffOps{

// Forward Declarations - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<typename DERIVED, std::size_t ORDER>
class DiffOpBase;

template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class DiffOpExpr; 

// (Binary Operator Function Objects) ==================================================
namespace internal{

// Structs for binary operations f(L1,L2) to combine .CoeffAt() - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// L1 + L2 
template<typename L, typename R>
struct diffop_bin_add_op
{
  L& m_lhs; 
  R& m_rhs; 
  diffop_bin_add_op(L& l, R& r) : m_lhs(l), m_rhs(r){}; 
  template<typename Cont>
  double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
  {
    return m_lhs.CoeffAt(weights, n_nodes_per_row, ith_node) + m_rhs.CoeffAt(weights, n_nodes_per_row, ith_node); 
  } 
}; 

// L1 - L2 
template<typename L, typename R>
struct diffop_bin_subtract_op
{
  L& m_lhs; 
  R& m_rhs; 
  diffop_bin_subtract_op(L& l, R& r) : m_lhs(l), m_rhs(r){}; 
  template<typename Cont>
  double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
  {
    return m_lhs.CoeffAt(weights, n_nodes_per_row, ith_node) - m_rhs.CoeffAt(weights, n_nodes_per_row, ith_node); 
  } 
}; 

// c * L
template<typename L, typename R>
struct diffop_left_mult_op
{
  L& m_lhs; 
  R& m_rhs; 
  diffop_left_mult_op(L& l, R& r) : m_lhs(l), m_rhs(r){}; 
  template<typename Cont>
  double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
  {
    return m_lhs * m_rhs.CoeffAt(weights, n_nodes_per_row, ith_node); 
  } 
}; 

// -L 
template<typename L>
struct diffop_negate_op
{
  L& m_lhs; 
  diffop_negate_op(L& l) : m_lhs(l){}; 
  template<typename Cont>
  double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
  {
    return -(m_lhs.CoeffAt(weights, n_nodes_per_row, ith_node)); 
  } 
}; 

// // composition: L1( L2( . ) ) // .compose() uses LinOpMixin expression templates... 
// struct linopXlinop_mult_op
// {
//   template<typename L1, typename L2>
//   auto operator()(const L1& A, const L2& B) const { return  (A.GetMat())*(B.GetMat()); }
// }; 

} // end namespace internal 

// Traits =====================================================================
// given a type, detect if it is derived from DiffOpBase<> - - - - - - - - - - - - 
namespace internal{
template<typename T, typename = void>
struct is_diffop_crtp_impl : std::false_type {};

template<typename T>
struct is_diffop_crtp_impl<T, std::void_t<typename T::is_diffop_tag>> : std::true_type {};

} // end namespace internal 

namespace traits{
template<typename T>
using is_diffop_crtp = internal::is_diffop_crtp_impl<std::remove_reference_t<std::remove_cv_t<T>>>; 
} // end namespace traits 

// given a type T see if it is an expression - - - - - - - - - - - - - - - - - - - - - - - - 
namespace internal{
template<typename T>
struct is_diffop_expr_crtp_impl : public std::false_type
{};

template<typename L, typename R, typename OP>
struct is_diffop_expr_crtp_impl<DiffOpExpr<L,R,OP>>: public std::true_type
{};
} // end namespace internal 

namespace traits{
template<typename T>
using is_diffop_expr_crtp = internal::is_diffop_expr_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>;
} // end namespace traits

// given a type T we may need to store it  - - - - - - - - - - - - - - - - - - - - - - - - - 
namespace traits{
// as a reference or a value in a binary expression
template<typename T, typename = void>
struct Storage_t 
{
  using type = T; 
}; 

template<typename DIFFOP_T>
struct Storage_t<DIFFOP_T, std::enable_if_t< is_diffop_crtp<DIFFOP_T>::value > >
{
  using type = std::conditional_t<
    std::is_lvalue_reference<DIFFOP_T>::value && !(is_diffop_expr_crtp<DIFFOP_T>::value), // if T is an lvalue DiffOPBase + not expression 
    DIFFOP_T, // store by lvalue DIFFOP_T& 
    std::remove_reference_t<DIFFOP_T> // else store rvalue by value
  >; 
}; 

} // end namespace traits 

} // end namespace DiffOps 

#endif // DiffOpTraits.hpp 