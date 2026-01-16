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

template<typename COEFF_T>
class CoeffMultExpr : public TimeDerivBase<CoeffMultExpr<COEFF_T>> 
{
  public:
    // Member Data ---------------
    using Lhs_t = typename LinOps::internal::Storage_t<COEFF_T>::type; 
    Lhs_t m_coeff; 

  public:
    // Constructors + Destructor ====================================
    // CoeffMultExpr()=delete; // necessary?  
    CoeffMultExpr(Lhs_t c_init, std::size_t order=1)
      : m_coeff(c_init), TimeDerivBase<CoeffMultExpr>(order, order+1)
    {
      // std::cout << "is linop? " << LinOps::internal::is_linop_crtp<COEFF_T>::value << std::endl; 
      // std::cout << "is linop? " << LinOps::internal::is_linop_crtp<Lhs_t>::value << std::endl; 
      // std::cout << "is coeff? " << LinOps::internal::is_coeffop_crtp<COEFF_T>::value << std::endl; 
      // std::cout << "is coeff? " << LinOps::internal::is_coeffop_crtp<Lhs_t>::value << std::endl; 
      // std::cout << "is lval? " << std::is_lvalue_reference<COEFF_T>::value << std::endl; 
      // std::cout << "is lval? " << std::is_lvalue_reference<Lhs_t>::value << std::endl; 
    }
    // LhsCoeffMultExpr(const LhsCoeffMultExpr& other)=delete; 
    // destructor 
    ~CoeffMultExpr()=default; 

    // Member Funcs =================================================== 
    template<typename Cont>
    auto CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      std::size_t offset = this->m_order * n_nodes_per_row; 
      if constexpr(LinOps::internal::is_coeffop_crtp<Lhs_t>::value){
        return v[ith_node + offset] * m_coeff.GetMat();  
      }
      else if constexpr(std::is_same<double, std::remove_cv_t<std::remove_reference_t<Lhs_t>>>::value){
        return v[ith_node + offset] * m_coeff; 
      }
      else{
        throw std::runtime_error("LhsCoeffMultExpr error: COEFF_T is not a double or CoeffOp."); 
      }
    } 

    // using LhsBase<LhsCoeffMultExpr<COEFF_T>>::toTuple; 
    std::string toString() const {return "hi from coeffMult"; }; 
}; 

template<
  typename Lhs,  
  typename = std::enable_if_t<
    LinOps::internal::is_coeffop_crtp<Lhs>::value
  >
>
auto operator*(Lhs&& c, NthTimeDeriv rhs)
{
  return CoeffMultExpr<Lhs>(std::forward<Lhs>(c), rhs.Order()); 
}

auto operator*(double c, NthTimeDeriv rhs)
{
  return CoeffMultExpr<double>(c, rhs.Order()); 
}

#endif // LhsCoeffMultExpr.hpp 