// LinOpExpr.hpp
//
//
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

// Expression of L,R, & BinOp --------------------------------------------
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr : public LinOpBase<LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t>>
{
  public:
    // types and constexpr flags
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    using RStorage_t = typename Storage_t<Rhs_t>::type;
  private:
    // member data   
    LStorage_t m_Lhs;
    RStorage_t m_Rhs;
    BinaryOp_t m_BinOp; 
    
  public:
    // Constructors 
    LinOpExpr(LStorage_t A, RStorage_t B, BinaryOp_t bin_op, MeshPtr_t m)
      : m_Lhs(A), m_Rhs(B), m_BinOp(bin_op), LinOpBase<LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t>>(m) 
    {};

    // destructors
    ~LinOpExpr()=default;

    // member funcs
    // getters
    LStorage_t& Lhs(){ return m_Lhs; }; 
    RStorage_t& Rhs(){ return m_Rhs; }; 
    // returns combination bin_op(A,B) of 2 stored LinOps----------------
    decltype(auto) GetMat()
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    decltype(auto) GetMat() const
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    // apply bin_op(A,B) * x --------------------
    // Discretization1D apply(const Discretization1D& d_arr) const 
    // {
    //   // we might constexpr branch this to L.apply(R.apply(d.values())) later ...
    //   // potentially introduces bug if 
    //   // R(.) maps boundaries -> L(.) uses boundaries 
    // };
    // sets both stored diffops to work on a mesh ------------------
    void set_mesh(MeshPtr_t m)
    {
      this->m_mesh_ptr=m; 
      if constexpr(is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      if constexpr(is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m);
    }; 
    // MeshPtr_t mesh() const 
    // {
    //   // if LHS is from linop base. given priority over RHS 
    //   if constexpr(is_linop_crtp<LStorage_t>::value)
    //   {
    //     // and it has a mesh 
    //     auto result = m_Lhs.mesh().lock(); 
    //     if(result) return result; 
    //   }
    //   // if RHS is from linop base. 
    //   if constexpr(is_linop_crtp<RStorage_t>::value)
    //   {
    //     // and it has a mesh 
    //     auto result = m_Lhs.mesh().lock(); 
    //     if(result) return result; 
    //   }
    //   // any other other cases give the stored mesh in the expression. presumably a nullptr. 
    //   return this->m_mesh_ptr; 
    // }
};

// Specialization for Unary operators --------------------------------------------
template<typename Lhs_t, typename UnaryOp_t>
class LinOpExpr<Lhs_t, void, UnaryOp_t> : public LinOpBase<LinOpExpr<Lhs_t, void, UnaryOp_t>>
{
  public:
    // types and constexpr flags
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    // using RStorage_t = typename Storage_t<Rhs_t>::type; // not storing a second argument anymore 
  private:
    // member data   
    LStorage_t m_Lhs;
    // RStorage_t m_Rhs; // not storing a second argument anymore 
    UnaryOp_t m_UnarOp; 
    
  public:
    // Constructors 
    LinOpExpr(LStorage_t A, UnaryOp_t unar_op, MeshPtr_t m)
      : m_Lhs(A), m_UnarOp(unar_op), LinOpBase<LinOpExpr<Lhs_t, void, UnaryOp_t>>(m)
    {};

    // destructors
    ~LinOpExpr()=default;

    // member funcs
    // getters
    LStorage_t& Lhs(){ return m_Lhs; }; 
    // RStorage_t& Rhs(){ return m_Rhs; }; 
    void Rhs(){}; // since nothing is stored. but still declared so that decltype(Rhs()) is still usable   
    // returns combination bin_op(A,B) of 2 stored LinOps----------------
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
    // apply bin_op(A,B) * x --------------------
    // Discretization1D apply(const Discretization1D& d_arr) const 
    // {
    //   // we might constexpr branch this to L.apply(R.apply(d.values())) later ...
    //   // potentially introduces bug if 
    //   // R(.) maps boundaries -> L(.) uses boundaries 
    // };
    // sets both stored diffops to work on a mesh ------------------
    void set_mesh(MeshPtr_t m)
    {
      // store mesh into expression 
      this->m_mesh_ptr=m; 
      // if Lhs is linop, call its set_mesh() as well
      if constexpr(is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      // if constexpr(is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m); // no longer altering a RHS anymore 
    }; 
    // MeshPtr_t mesh() const 
    // {
    //   // if LHS is from linop base. given priority over RHS 
    //   if constexpr(is_linop_crtp<LStorage_t>::value)
    //   {
    //     // and it has a mesh 
    //     auto result = m_Lhs.mesh().lock(); 
    //     if(result) return result; 
    //   }
    //   // // if RHS is from linop base. dont need to check RHS anymore  
    //   // if constexpr(is_linop_crtp<RStorage_t>::value)
    //   // {
    //   //   // and it has a mesh 
    //   //   const MeshPtr_t& result = m_Rhs.mesh(); 
    //   //   if(result) return result; 
    //   // }
    //   // any other other cases give the stored mesh in the expression. presumably a nullptr. 
    //   return this->m_mesh_ptr; 
    // }
};

#endif // LinOpExpr.hpp