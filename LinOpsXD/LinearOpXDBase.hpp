// LinearOpXDBase.hpp
//
//
//
// JAF 12/7/2025

#ifndef LINEAROPXDBASE_H
#define LINEAROPXDBASE_H

#include<iostream>
#include<cstdint>
#include<Eigen/Core>
#include<type_traits>
#include "MeshXD.hpp"
#include "DiscretizationXD.hpp" 

// CTRP base for 1D differential operator -------------------------------
template<typename Derived>
class LinOpXDBase 
{

  // type defs 
  public:
    typedef struct{} is_linopxd_tag; // to tell if a class derived from LinOpBase<> 
    using Derived_t = Derived; // so Plugin can access grand child class

  protected:
    // member data
    MeshXDPtr_t m_mesh_ptr; // all LinOps keep pointers to mesh they operate on 
    // member data for matrix kept in derived classes

  public:
    // Constructors ---------------------------------------------------------- 
    LinOpXDBase(MeshXDPtr_t m=nullptr){set_mesh(m);}; 
    // Destructors ----------------------------------------------------------- 
    ~LinOpXDBase()=default; 
    // member functions -------------------------------------------------------
    // member functions. implemented by derived class -------------------------
    auto& GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    const auto& GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat(); 
    };

    // multiply the underlying expression with Discretization's underlying vecXd
    DiscretizationXD apply(const DiscretizationXD& d) const 
    {
        return static_cast<const Derived_t*>(this)->apply(d);
    };

    // fit operator to a mesh of rectangular domain.
    void set_mesh(MeshXDPtr_t m) 
    {
      // ensure we aren't resetting the mesh again, or setting to nullptr
      // if((m==nullptr)||(m==m_mesh_ptr)) return; 
      // m_mesh_ptr = m; // take shared ownership of mesh
      static_cast<Derived*>(this)->set_mesh(m);
    };
    const MeshXDPtr_t& mesh() const 
    { 
      return m_mesh_ptr; 
    } 

    // // operators ---------------------------------------------------
    // // composition of linear of L1(L2( . ))
    // template<typename DerivedInner> 
    // auto compose(DerivedInner&& InnerOp) &
    // {
    //   static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
    //   return compose_impl<DerivedInner,Derived_t&>(std::forward<DerivedInner>(InnerOp)); 
    // }; // end .compose(other) & lvalue overload 

    // // // composition of linear of L1(L2( . ))
    // template<typename DerivedInner>
    // auto compose(DerivedInner&& InnerOp) && 
    // {
    //   static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
    //   return compose_impl<DerivedInner,Derived_t&&>(std::forward<DerivedInner>(InnerOp)); 
    // }; // end .compose(other) && rvalue overload  

  // private:
    // not accessibles --------------------------------------------------------------------------------------------
    // composition of linear of L1(L2( . ))
    // template<typename DerivedInner, typename Lhs_t = Derived_t> 
    // auto compose_impl(DerivedInner&& InnerOp)
    // {
    //   static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
    //   using Rhs_t = std::remove_reference_t<DerivedInner>;
    //   using LStorage_t = typename Storage_t<Lhs_t>::type;
    //   using RStorage_t = typename Storage_t<Rhs_t>::type;
    //   if constexpr(is_add_expr<std::remove_reference_t<DerivedInner>>::value){
    //     return compose(InnerOp.Lhs())+compose(InnerOp.Rhs());
    //   }
    //   else if constexpr(is_scalar_multiply_expr<std::remove_reference_t<DerivedInner>>::value){
    //     return InnerOp.Lhs() * compose(InnerOp.Rhs()); 
    //   }
    //   else{
    //     auto bin_op = [](const LStorage_t& A, const RStorage_t& B){return A.GetMat()*B.GetMat(); };
    //     using Op_t = make_flagged_t<decltype(bin_op), OperatorComposition_t>;
    //     if constexpr(std::is_lvalue_reference<Lhs_t>::value){
    //       return LinOpExpr<Lhs_t, Rhs_t, Op_t>(
    //       std::forward<Lhs_t>(static_cast<Lhs_t&>(*this)), 
    //       std::forward<DerivedInner>(InnerOp), 
    //       static_cast<Op_t>(bin_op)
    //       ); 
    //     }
    //     else{
    //       return LinOpExpr<Lhs_t, Rhs_t, Op_t>(
    //       std::forward<Lhs_t>(static_cast<Lhs_t>(*this)), 
    //       std::forward<DerivedInner>(InnerOp), 
    //       static_cast<Op_t>(bin_op)
    //       ); 
    //     }
    //   } // end else 
    // }; // end .compose_impl(other) 

}; // end LinOpBase

#endif // LinearOpBase.hpp