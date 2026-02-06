// LinOpExpr.hpp
//
// Complie time expression template that stores Lhs, Rhs, BinOp (UnarOp), and Mesh
//
// JAF 12/7/2025

#ifndef LINOPEXPR_H
#define LINOPEXPR_H

#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinearOpBase.hpp"
#include "LinOpTraits.hpp"

namespace LinOps{

// Expression of L,R, BinOp, + Mesh =====================================================
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr : public LinOpMixIn< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >, public LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >, public LinOpBaseXD< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >
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
    LinOpExpr()=delete; 
    // from Lhs, Rhs, BinOp, + Mesh 
    LinOpExpr(LStorage_t A, RStorage_t B, BinaryOp_t bin_op)
      : m_Lhs(A), m_Rhs(B), m_BinOp(bin_op)
    {
      constexpr bool both_linop = traits::is_linop_crtp<Lhs_t>::value && traits::is_linop_crtp<Rhs_t>::value;
      constexpr bool both_1d = traits::is_1dim_linop_crtp<Lhs_t>::value && traits::is_1dim_linop_crtp<Rhs_t>::value;
      constexpr bool both_xd = traits::is_xdim_linop_crtp<Lhs_t>::value && traits::is_xdim_linop_crtp<Rhs_t>::value;
      if constexpr(both_linop){
        static_assert(both_1d || both_xd, "Error constructing LinOpExpr: tried to mix strictly 1D operator with strictly XD operator"); 
      }
    };
    // copy 
    LinOpExpr(const LinOpExpr& other)=default; 
    // destructor
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================

    // returns combination bin_op(A,B) of 2 stored LinOps ----------------
    decltype(auto) GetMat()
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    decltype(auto) GetMat() const
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };

    // fixes ambiguous .apply() .set_mesh()? 
    using LinOpBase1D<LinOpExpr<Lhs_t,Rhs_t,BinaryOp_t> >::apply; 
    using LinOpBaseXD<LinOpExpr<Lhs_t,Rhs_t,BinaryOp_t> >::apply; 
    using LinOpBase1D<LinOpExpr<Lhs_t,Rhs_t,BinaryOp_t>>::set_mesh; 
    using LinOpBaseXD<LinOpExpr<Lhs_t,Rhs_t,BinaryOp_t>>::set_mesh; 

    // Expr Only ==============================================================

    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; }; 
    const LStorage_t& Lhs() const { return m_Lhs; }; 
    // rhs 
    RStorage_t& Rhs(){ return m_Rhs; }; 
    const RStorage_t& Rhs() const { return m_Rhs; }; 
   
};

// Specialization for Unary operators --------------------------------------------
template<typename Lhs_t, typename UnaryOp_t>
class LinOpExpr<Lhs_t, void, UnaryOp_t> : public LinOpMixIn< LinOpExpr<Lhs_t, void, UnaryOp_t> >, public LinOpBase1D< LinOpExpr<Lhs_t, void, UnaryOp_t> >, public LinOpBaseXD< LinOpExpr<Lhs_t, void, UnaryOp_t> >
{
  public:
    // Type Defs ------------------------------------------------------------------
    using LStorage_t = typename traits::Storage_t<Lhs_t>::type;
    using RStorage_t = void; // not storing a second argument anymore 
    using Operator_t = UnaryOp_t;

  private:
    // Member Data -------------------------------------------------------------
    LStorage_t m_Lhs;
    // RStorage_t m_Rhs; // not storing a second argument anymore 
    UnaryOp_t m_UnarOp; 
    
  public:
    // Constructors / Destructors ===============================================
    // default 
    LinOpExpr()=delete; 
    // from lhs, unar_op, + mesh
    LinOpExpr(LStorage_t A, UnaryOp_t unar_op)
      : m_Lhs(A), m_UnarOp(unar_op)
    {
      constexpr bool A_is_linop = traits::is_linop_crtp<Lhs_t>::value;
      static_assert(A_is_linop, "Error constructing Unary LinOpExpr: single linop A is not linop!"); 
    };
    // copy 
    LinOpExpr(const LinOpExpr& other)=default; 
    // destructors
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================

    // returns Op( A ) of 1 stored LinOps --------------------------------------
    decltype(auto) GetMat()
    {
      // return m_BinOp(m_Lhs, m_Rhs); 
      return m_UnarOp(m_Lhs); 
    };
    decltype(auto) GetMat() const
    {
      // return m_BinOp(m_Lhs, m_Rhs); 
      return m_UnarOp(m_Lhs); 
    };

    // fixes ambiguous .apply() .set_mesh()? 
    using LinOpBase1D<LinOpExpr<Lhs_t,void,UnaryOp_t>>::apply; 
    using LinOpBaseXD<LinOpExpr<Lhs_t,void,UnaryOp_t>>::apply; 
    using LinOpBase1D<LinOpExpr<Lhs_t,void,UnaryOp_t>>::set_mesh; 
    using LinOpBaseXD<LinOpExpr<Lhs_t,void,UnaryOp_t>>::set_mesh; 

    // Expr Only ==============================================================

    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; };     
    const LStorage_t& Lhs() const { return m_Lhs; };     
    // since nothing is stored. but still declared so that decltype(Rhs()) is still usable
    void Rhs(){};    
    
};

} // end namespace LinOps 

#endif // LinOpExpr.hpp