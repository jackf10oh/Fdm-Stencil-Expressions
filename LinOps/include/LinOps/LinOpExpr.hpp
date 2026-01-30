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

namespace LinOps{

// Expression of L,R, BinOp, + Mesh =====================================================
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr : public LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >
{
  public:
    // // clarify ambiguous operators from LinOpMixin; 
    // using LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >::operator+; 
    // using LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >::operator-; 
    // using LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >::left_scalar_mult_impl; 
    // using LinOpBase1D< LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t> >::compose; 
    // Type Defs --------------------------------------
    using LStorage_t = typename traits::Storage_t<Lhs_t>::type;
    using RStorage_t = typename traits::Storage_t<Rhs_t>::type;

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
    {};
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

    // sets both stored diffops to work on a mesh1d ------------------
    void set_mesh(const Mesh1D_SPtr_t& m)
    {
      // this->m_mesh_ptr=m; 
      if constexpr(traits::is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      if constexpr(traits::is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m);
    }; 

    // ... on a meshxd ------------------
    void set_mesh(const MeshXD_SPtr_t& m)
    {
      // this->m_mesh_ptr=m; 
      if constexpr(traits::is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      if constexpr(traits::is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m);
    }; 

    // get Mesh1D_SPtr_t from stored lhs first or rhs 
    Mesh1D_SPtr_t get_mesh1d() const
    {
      if constexpr(LinOps::traits::is_linop_crtp<Lhs_t>::value){
        return m_Lhs.get_mesh1d(); 
      }
      else{
        return m_Rhs.get_mesh1d(); 
      }
    }

    // get MeshXD_SPtr_t ... 
    MeshXD_SPtr_t get_meshxd() const
    {
      if constexpr(LinOps::traits::is_linop_crtp<Lhs_t>::value){
        return m_Lhs.get_meshxd(); 
      }
      else{
        return m_Rhs.get_meshxd(); 
      }
    }

    // get weak_ptr Mesh1D_WPtr_t from stored lhs first or rhs 
    Mesh1D_WPtr_t get_weak_mesh1d() const
    {
      if constexpr(LinOps::traits::is_linop_crtp<Lhs_t>::value){
        return m_Lhs.get_weak_mesh1d(); 
      }
      else{
        return m_Rhs.get_weak_mesh1d(); 
      }
    }

    // get weak_ptr MeshXD_WPtr_t ... 
    Mesh1D_WPtr_t get_weak_meshxd() const
    {
      if constexpr(LinOps::traits::is_linop_crtp<Lhs_t>::value){
        return m_Lhs.get_weak_meshxd(); 
      }
      else{
        return m_Rhs.get_weak_meshxd(); 
      }
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

// Specialization for Unary operators --------------------------------------------
template<typename Lhs_t, typename UnaryOp_t>
class LinOpExpr<Lhs_t, void, UnaryOp_t> : public LinOpBase1D<LinOpExpr<Lhs_t, void, UnaryOp_t>>
{
  public:
    // Type Defs ------------------------------------------------------------------
    using LStorage_t = typename traits::Storage_t<Lhs_t>::type;
    using RStorage_t = void; // not storing a second argument anymore 

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
    {};
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

    // sets both stored diffops to work on a mesh ------------------------------------
    void set_mesh(const Mesh1D_SPtr_t& m)
    {
      // if Lhs is linop, call its set_mesh() as well
      if constexpr(traits::is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
    }; 

    // get Mesh1D_SPtr_t from stored lhs 
    Mesh1D_SPtr_t get_mesh1d() const
    {
      static_assert(LinOps::traits::is_linop_crtp<Lhs_t>::value, "LHS must be linop in Unary Expression");
      return m_Lhs.get_mesh1d(); 
    }

    // get MeshXD_SPtr_t ... 
    MeshXD_SPtr_t get_meshxd() const
    {
      static_assert(LinOps::traits::is_linop_crtp<Lhs_t>::value, "LHS must be linop in Unary Expression");
      return m_Lhs.get_meshxd(); 
    }

    // get weak_ptr Mesh1D_WPtr_t from stored lhs
    Mesh1D_WPtr_t get_weak_mesh1d() const
    {
      static_assert(LinOps::traits::is_linop_crtp<Lhs_t>::value, "LHS must be linop in Unary Expression");
      return m_Lhs.get_weak_mesh1d(); 
    }

    // get weak_ptr MeshXD_WPtr_t ... 
    Mesh1D_WPtr_t get_weak_meshxd() const
    {
      static_assert(LinOps::traits::is_linop_crtp<Lhs_t>::value, "LHS must be linop in Unary Expression");
      return m_Lhs.get_weak_meshxd(); 
    }

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