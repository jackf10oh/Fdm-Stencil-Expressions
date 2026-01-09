// LinOpExpr.hpp
//
// Complie time expression template that stores Lhs, Rhs, BinOp (UnarOp), and Mesh
//
// JAF 12/7/2025

#ifndef LINOPEXPR_H
#define LINOPEXPR_H

#include<cstdint>
#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinearOpBase.hpp"
#include "LinOpTraits.hpp"

// Expression of L,R, BinOp, + Mesh =====================================================
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr : public LinOpBase<LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t>>
{
  public:
    // Type Defs --------------------------------------
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    using RStorage_t = typename Storage_t<Rhs_t>::type;
  private:
    // Member Data ---------------------------------------------
    LStorage_t m_Lhs;
    RStorage_t m_Rhs;
    BinaryOp_t m_BinOp; 
    
  public:
    // Constructors / Destructors =============================================================
    // from Lhs, Rhs, BinOp, + Mesh 
    LinOpExpr(LStorage_t A, RStorage_t B, BinaryOp_t bin_op, MeshPtr_t m)
      : m_Lhs(A), m_Rhs(B), m_BinOp(bin_op), LinOpBase<LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t>>(m) 
    {};
    // destructor
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================
    // Expr Only ==============================================================
    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; }; 
    // rhs 
    RStorage_t& Rhs(){ return m_Rhs; }; 
   
    // Must be implemented ========================================================
    // returns combination bin_op(A,B) of 2 stored LinOps ----------------
    decltype(auto) GetMat()
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    decltype(auto) GetMat() const
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };

    // sets both stored diffops to work on a mesh ------------------
    void set_mesh(MeshPtr_t m)
    {
      this->m_mesh_ptr=m; 
      if constexpr(is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      if constexpr(is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m);
    }; 
};

// Specialization for Unary operators --------------------------------------------
template<typename Lhs_t, typename UnaryOp_t>
class LinOpExpr<Lhs_t, void, UnaryOp_t> : public LinOpBase<LinOpExpr<Lhs_t, void, UnaryOp_t>>
{
  public:
    // Type Defs ------------------------------------------------------------------
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    using RStorage_t = void; // not storing a second argument anymore 
  private:
    // Member Data -------------------------------------------------------------
    LStorage_t m_Lhs;
    // RStorage_t m_Rhs; // not storing a second argument anymore 
    UnaryOp_t m_UnarOp; 
    
  public:
    // Constructors / Destructors ===============================================
    // from lhs, unar_op, + mesh
    LinOpExpr(LStorage_t A, UnaryOp_t unar_op, MeshPtr_t m)
      : m_Lhs(A), m_UnarOp(unar_op), LinOpBase<LinOpExpr<Lhs_t, void, UnaryOp_t>>(m)
    {};

    // destructors
    ~LinOpExpr()=default;

    // Member Funcs ======================================================================
    // Expr Only ==============================================================
    // Lhs/Rhs getters --------------------------------------------------
    // Lhs 
    LStorage_t& Lhs(){ return m_Lhs; };     
    // since nothing is stored. but still declared so that decltype(Rhs()) is still usable
    void Rhs(){};    
    
    // Must be implemented ========================================================
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

    // sets both stored diffops to work on a mesh ------------------------------------
    void set_mesh(MeshPtr_t m)
    {
      // store mesh into expression 
      this->m_mesh_ptr=m; 
      // if Lhs is linop, call its set_mesh() as well
      if constexpr(is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
    }; 
};

#endif // LinOpExpr.hpp