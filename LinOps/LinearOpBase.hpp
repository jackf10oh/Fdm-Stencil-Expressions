// LinearOpBase.hpp
//
//
//
// JAF 12/7/2025

#ifndef LINEAROPBASE_H
#define LINEAROPBASE_H

#include<iostream>
#include<cstdint>
#include<eigen3/Eigen/Core>
#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinOpTraits.hpp"

#ifndef LINOP_PLUGIN
template<typename T>
using LinOpMixIn = std::null_t; 
#else
template<typename T>
using LinOpMixIn = LINOP_PLUGIN<T>; 
#endif

// CTRP base for 1D differential operator -------------------------------
template<typename Derived>
class LinOpBase : public LinOpMixIn<LinOpBase<Derived>> 
{
  using MeshPtr_t = std::shared_ptr<Mesh1D>;
  public:
    using Derived_t = Derived; 

  protected:
    // member data
    MeshPtr_t m_mesh_ptr; // all LinOps keep pointers to mesh they operate on 
    // member data for matrix kept in derived classes

  public:
    // Constructors ---------------------------------------------------------- 
    LinOpBase(MeshPtr_t m=nullptr):m_mesh_ptr(m){}; 
    // member functions -------------------------------------------------------
    // member functions. implemented by derived class -------------------------
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };

    // multiply the underlying expression with Discretization's underlying vecXd
    Discretization1D apply(const Discretization1D& d) const 
    {
      if constexpr (has_apply<LinOpMixIn<LinOpBase<Derived>>>::value)
      {
        return LinOpMixIn<LinOpBase<Derived>>::apply(d);
      }
      else
      {
        return static_cast<Derived*>(this)->apply(d);
      }
    };

    // fit operator to a mesh of rectangular domain 
    void set_mesh(MeshPtr_t m) 
    {
      // ensure we aren't resetting the mesh again, or setting to nullptr
      if((m==nullptr)||(m==m_mesh_ptr)) return; 
      m_mesh_ptr = m; // take shared ownership of mesh
      static_cast<Derived*>(this)->set_mesh(m);
    };

    // operators ---------------------------------------------------
    // composition of linear of L1(L2( . ))
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!");
      
      if constexpr(is_add_expr<std::remove_reference_t<DerivedInner>>::value){
        // std::cout << "composing L(L1+L2)" << std::endl;
        return compose(InnerOp.Lhs())+compose(InnerOp.Rhs());
      }
      else if constexpr(is_scalar_multiply_expr<std::remove_reference_t<DerivedInner>>::value){
        // std::cout << "composing L(c*L) = c*L(L)" << std::endl;
        return InnerOp.Lhs() * compose(InnerOp.Rhs()); 
      }
      else{
        // std::cout << "some other composition..." << std::endl;
        using Lhs_t = Derived;
        using Rhs_t = std::remove_reference_t<DerivedInner>;
        using LStorage_t = typename Storage_t<Lhs_t>::type;
        using RStorage_t = typename Storage_t<Rhs_t>::type;

        auto bin_op = [](const LStorage_t& A, const RStorage_t& B){return A.GetMat()*B.GetMat(); };
        using Op_t = make_flagged_t<decltype(bin_op), OperatorComposition_t>;
        return LinOpExpr<Lhs_t, Rhs_t, Op_t>(
          *static_cast<Derived*>(this), 
          std::forward<DerivedInner>(InnerOp), 
          static_cast<Op_t>(bin_op)
        );
      } // end else 
    }; // end .compose(other) 

}; // end LinOpBase

#endif // LinearOpBase.hpp