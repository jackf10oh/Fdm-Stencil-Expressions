// DiffOpExpr.hpp
//
// Binary and Unary expressions of DiffOps. 
//
// JAF 2/5/2026 

#ifndef DIFFOPEXPR_H
#define DIFFOPEXPR_H 

#include "DiffOpBase.hpp"

namespace DiffOps{

// Expression of L,R, BinOp =====================================================
template<typename Lhs_t, typename Rhs_t, template<typename...> class BinaryOp_t>
class DiffOpExpr : public DiffOpBase<DiffOpExpr<Lhs_t, Rhs_t, BinaryOp_t<Lhs_t,Rhs_t>>>
{
  public:
    // Type Defs --------------------------------------
    using LStorage_t = typename traits::Storage_t<Lhs_t>::type;
    using RStorage_t = typename traits::Storage_t<Rhs_t>::type;
    using Operator_t = BinaryOp_t;

  private:
    // Member Data ---------------------------------------------
    LStorage_t m_Lhs;
    RStorage_t m_Rhs;
    BinaryOp_t m_BinOp; 
    
  public:
    // Constructors + Destructor =============================================================
    // default 
    DiffOps()=delete; 
    // from Lhs + Rhs 
    DiffOps(LStorage_t A, RStorage_t B)
      : m_Lhs(A), m_Rhs(B), m_BinOp(m_Lhs,m_Rhs) 
    { /* some types of static asserts might be needed here ...*/};
    // copy 
    DiffOps(const DiffOps& other)=default; 
    // destructor
    ~DiffOps()=default;

    // Member Funcs ======================================================================
    // Matrix coefficient from finite difference weights  
    template<typename Cont>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return m_BinOp(weights, n_nodes_per_row, ith_node); 
    } 

    // Expr Only ==============================================================
    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; }; 
    const LStorage_t& Lhs() const { return m_Lhs; }; 
    // rhs 
    RStorage_t& Rhs(){ return m_Rhs; }; 
    const RStorage_t& Rhs() const { return m_Rhs; }; 
};

// specialization for only L, UnaryOp =====================================================
template<typename Lhs_t, template<typename...> class UnaryOp_t>
class DiffOpExpr<Lhs_t, void, UnaryOp_t> : public DiffOpBase<DiffOpExpr<Lhs_t, void, UnaryOp_t<Lhs_t>>>
{
  public:
    // Type Defs --------------------------------------
    using LStorage_t = typename traits::Storage_t<Lhs_t>::type;
    using RStorage_t = void; // no longer storing ... 
    using Operator_t = UnaryOp_t;

  private:
    // Member Data ---------------------------------------------
    LStorage_t m_Lhs;
    // RStorage_t m_Rhs; // no longer storing 
    BinaryOp_t m_UnarOp; 
    
  public:
    // Constructors + Destructor =============================================================
    // default 
    DiffOps()=delete; 
    // from Lhs + Rhs 
    DiffOps(LStorage_t A, RStorage_t B)
      : m_Lhs(A), m_Rhs(B), m_BinOp(m_Lhs,m_Rhs) 
    { /* some types of static asserts might be needed here ...*/};
    // copy 
    DiffOps(const DiffOps& other)=default; 
    // destructor
    ~DiffOps()=default;

    // Member Funcs ======================================================================
    // Matrix coefficient from finite difference weights  
    template<typename Cont>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return m_BinOp(weights, n_nodes_per_row, ith_node); 
    } 

    // Expr Only ==============================================================
    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; }; 
    const LStorage_t& Lhs() const { return m_Lhs; }; 
    // rhs 
    RStorage_t& Rhs(){ return m_Rhs; }; 
    const RStorage_t& Rhs() const { return m_Rhs; }; 
};

} // end namespace DiffOps 

#endif // DiffOpExpr.hpp 
