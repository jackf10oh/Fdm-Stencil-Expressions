// LhsSumExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSSUMEXPR_H
#define LHSSUMEXPR_H 

#include "LhsBase.hpp"

template<typename... Args>
class LhsSumExpr : public LhsBase<LhsSumExpr<Args...>> 
{
  public:
    std::tuple<Args...> m_args; 
  public:
    // Constructors + Destructor ================================
    LhsSumExpr()=delete; 
    LhsSumExpr(std::tuple<Args...> tup_init, std::size_t order) 
      : m_args( tup_init ), LhsBase<LhsSumExpr<Args...>>(order, order+1) 
    {} 
    // destructor 
    ~LhsSumExpr()=default; 

    // Member Funcs ====================================================
    // auto toTuple() { return m_args; }
    auto& toTuple() &
    {
      // std::cout << "Lvalue SumExpr .toTuple()" << std::endl; 
      return m_args; // lvalue -> reference
    }
    auto toTuple() &&
    {
      // std::cout << "Rvalue SumExpr .toTuple()" << std::endl; 
      return std::move(m_args); // rvalue -> rvalue ref
    }
    std::string toString() const {return "hi from sum"; }; 
};

template<typename LHS, typename RHS>
auto make_lhssumexpr_helper(LHS&& lhs, RHS&& rhs)
{
  std::size_t m = std::max(lhs.Order(), rhs.Order()); 
  auto left_tup = std::forward<LHS>(lhs).toTuple(); 
  auto right_tup = std::forward<RHS>(rhs).toTuple(); 
  auto cat = std::tuple_cat(left_tup, right_tup); 
  return LhsSumExpr(cat, m); 
}

#endif // LhsSumExpr.hpp 