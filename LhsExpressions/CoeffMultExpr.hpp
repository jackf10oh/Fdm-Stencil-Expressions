// CoeffMultExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSCOEFFMULTEXPR_H
#define LHSCOEFFMULTEXPR_H

#include "TimeDerivBase.hpp"
#include "NthTimeDeriv.hpp" 
#include "../LinOps/LinOpTraits.hpp"

template<typename COEFF_T, typename RHS_T>
class CoeffMultExpr : public TimeDerivBase<CoeffMultExpr<COEFF_T,RHS_T>> 
{
  public:
    // Type Defs --------------------- 
    using Lhs_t = typename LinOps::internal::Storage_t<COEFF_T>::type; 
    using Rhs_t = typename std::conditional_t<
      std::conjunction_v<
        std::is_same<std::remove_cv_t<std::remove_reference_t<RHS_T>>, NthTimeDeriv>, 
        std::is_lvalue_reference<RHS_T>
      >, 
      RHS_T, 
      std::remove_reference_t<RHS_T>
    >; 

    // Member Data ---------------
    Lhs_t m_coeff; 
    Rhs_t m_rhs; 
    /* store the right hand side somehow... 
    if its an lvalue store a reference. 
    if its a rvalue just store a copy.*/

  public:
    // Constructors + Destructor ====================================
    // CoeffMultExpr()=delete; // necessary?  
    CoeffMultExpr(Lhs_t c_init, Rhs_t rhs_init)
      : m_coeff(c_init), m_rhs(rhs_init), TimeDerivBase<CoeffMultExpr>(rhs_init.Order())
    {}
    CoeffMultExpr(const CoeffMultExpr& other)=default; 
    // destructor 
    ~CoeffMultExpr()=default; 

    // Member Funcs =================================================== 
    // L/R Getters
    Lhs_t& Lhs(){return m_coeff; } 
    Rhs_t& Rhs(){return m_rhs; } 
    const Lhs_t& Lhs() const {return m_coeff; } 
    const Rhs_t& Rhs() const {return m_rhs; } 

    // Must be implemented for TimeDerivBase.  Get coeff * rhs's CoeffAt() value  
    template<typename Cont>
    auto CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      if constexpr(LinOps::internal::is_coeffop_crtp<Lhs_t>::value){
        return m_coeff.GetMat() * m_rhs.CoeffAt(v,n_nodes_per_row,ith_node);  
      }
      else if constexpr(std::is_same<double, std::remove_cv_t<std::remove_reference_t<Lhs_t>>>::value){
        return m_coeff * m_rhs.CoeffAt(v,n_nodes_per_row,ith_node); 
      }
      else{
        throw std::runtime_error("LhsCoeffMultExpr error: COEFF_T is not a double or CoeffOp."); 
      }
    } 

    // using LhsBase<LhsCoeffMultExpr<COEFF_T>>::toTuple; 
    std::string toString() const {return "hi from coeffMult!"; }; 
}; 

// should be deleted for SumExprs...... 
template<
  typename Lhs, 
  typename Rhs,  
  typename = std::enable_if_t<
    std::conjunction_v<
    LinOps::internal::is_coeffop_crtp<Lhs>,
    is_timederiv_crtp<Rhs>
    >
  >
>
auto operator*(Lhs&& c, Rhs&& rhs)
{
  return CoeffMultExpr<Lhs,Rhs>(std::forward<Lhs>(c), std::forward<Rhs>(rhs)); 
}

// Operator for double c, TimeDeriv Ut making expression c*Ut 
template<
  typename Rhs, 
  typename = std::enable_if_t<
    is_timederiv_crtp<Rhs>::value
  >
>
auto operator*(double c, Rhs&& rhs)
{
  return CoeffMultExpr<double,Rhs>(c, std::forward<Rhs>(rhs)); 
}

#endif // LhsCoeffMultExpr.hpp 