// TimeDerivBase.hpp
//
//
//
// JAF 1/15/2026 

#ifndef TIMEDERIVBASE_H
#define TIMEDERIVBASE_H 

#include "../FDStencils/CoeffOpBase.hpp"

// traits =====================================================
template<typename T, typename = void> 
struct is_timederiv_crtp_impl : public std::false_type{}; 

template<typename T> 
struct is_timederiv_crtp_impl<T, std::void_t<typename T::is_timederiv_tag>>: public std::true_type{}; 

template<typename T>
using is_timederiv_crtp = is_timederiv_crtp_impl<std::remove_cv_t<std::remove_reference_t<T>>>; 

// Base Definition =====================================================
template<typename Derived>
class TimeDerivBase
{
  public:
    // Type Defs --------------------------
    struct is_timederiv_tag{}; 
  public:
    // Member Data ----------------------------------
    std::size_t m_order; // nth derivative in time 
    std::size_t m_n_nodes; // # of nodes required for nth deriv 
  public:
    // Constructors + Destructor ============================================
    TimeDerivBase()=default; 
    TimeDerivBase(std::size_t order_init, std::size_t n_nodes_init)
      : m_order(order_init), m_n_nodes(n_nodes_init)
    {}
    TimeDerivBase(const TimeDerivBase& other)=default; 
    ~TimeDerivBase()=default;
    // Member Funcs ========================================================== 
    // Current Order of expression  
    std::size_t Order() const {return m_order; }; 

    template<typename Cont>
    decltype(auto) CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return static_cast<const Derived*>(this)->CoeffAt(v, n_nodes_per_row, ith_node); 
    } 

    // going to need set time and set mesh so that coeff ops can interac with Lhs... 

    // packs this pointer into iterable tuple object. note: SumExpr will override it
    auto toTuple() &
    {
      // std::cout << "Lvalue Base .toTuple()" << std::endl; 
      return std::forward_as_tuple(static_cast<Derived&>(*this)); // lvalue -> reference
    }
    auto toTuple() &&
    {
      // std::cout << "Rvalue Base .toTuple()" << std::endl; 
      return std::make_tuple(*this); // rvalue -> rvalue ref
    }

    std::string toString() const {return "hi from base"; }; 

    // Operators ================================
    template<typename RHS, typename = std::enable_if_t<is_timederiv_crtp<RHS>::value>>
    auto operator+(RHS&& rhs) & 
    {
      return make_lhssumexpr_helper(static_cast<Derived&>(*this), std::forward<RHS>(rhs)); 
    }

    template<typename RHS, typename = std::enable_if_t<is_timederiv_crtp<RHS>::value>>
    auto operator+(RHS&& rhs) &&
    {
      return make_lhssumexpr_helper(std::move(static_cast<Derived&>(*this)), std::forward<RHS>(rhs)); 
    }

}; 

#endif