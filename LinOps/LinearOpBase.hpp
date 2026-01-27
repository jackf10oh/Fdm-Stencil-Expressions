// LinearOpBase.hpp
//
// CTRP base for 1D differential operator L
// where L works on function discretizations across x0 , ..., xN 
//
// JAF 12/7/2025

#ifndef LINEAROPBASE_H
#define LINEAROPBASE_H

#include<iostream>
#include<cstdint>
#include<Eigen/Core>
#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinOpTraits.hpp"
#include "LinOpExpr.hpp"

namespace LinOps{

namespace internal{
#ifndef LINOP_PLUGIN
// use an empty struct if no plugin given.
struct empty{}; 
template<typename T> 
using LinOpMixIn = empty; 
#else
template<typename T>
using LinOpMixIn = LINOP_PLUGIN<T>; 
#endif
} // end namespace internal 

template<typename Derived>
class LinOpBase : public internal::LinOpMixIn<LinOpBase<Derived>> 
{
  // friend classes. LINOP_PLUGIN<T> can use private members of T 
  friend internal::LinOpMixIn<LinOpBase>; 

  // Type Defs -----------------------------------------------------
  public:
    typedef struct{} is_linop_tag; // to tell if a class derived from LinOpBase<> 
    using Derived_t = Derived; // so Plugin can access grand child class

  protected:
    // Member Data -------------------------------------------------
    Mesh1D_WPtr_t m_mesh_ptr; // all LinOps keep weak pointers to mesh they operate on
    // member data for matrix kept in derived classes

  public:
    // Constructors / destructors ============================================================ 
    LinOpBase()=default; 
    LinOpBase(Mesh1D_WPtr_t m) : m_mesh_ptr(m){}; 
    LinOpBase(const LinOpBase& other)=default; 
    ~LinOpBase()=default; 

    // member functions ============================================================

    // must be implemented by derived class =====================================
    // get the matrix form of linear operator ---------------------------------
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat(); 
    };

    // defaults given. can be overriden =====================================
    // multiply the underlying expression with Discretization's underlying vecXd
    Discretization1D apply(const Discretization1D& d) const 
    {
      if constexpr (traits::has_apply<internal::LinOpMixIn<LinOpBase<Derived>>>::value)
      {
        return internal::LinOpMixIn<LinOpBase<Derived>>::apply(d);
      }
      else
      {
        // no longer checking... Discretization1D and LinOp can point to 2 differernt Mesh1D! 
        // Mesh1D_WPtr_t m = d.mesh_weak_ptr(); 
        // if(m_mesh_ptr.owner_before(m) || m.owner_before(m_mesh_ptr)){
        //   throw std::runtime_error("Linear Operator L and discretization d point to different Mesh1D!");
        // }
        Eigen::VectorXd v = GetMat() * d.values();  // calculate A*b
        Discretization1D result(std::move(v), this->m_mesh_ptr); // move A*b into result's values
        return result;
      }
    };

    // fit operator to a mesh of rectangular domain ----------------------------------------
    void set_mesh(const Mesh1D_SPtr_t& m) 
    {
      // // ensure we aren't resetting the mesh again
      // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // // do nothing on nullptr. or throw an error 
      // auto locked = m.lock(); 
      // if(!locked) return; 
      // m_mesh_ptr = m; // store the mesh  
      // // perform work on locked 
      static_cast<Derived*>(this)->set_mesh(m);
    };
    
    // return lock of Mesh1D pointed to  ----------------------------------------------
    Mesh1D_SPtr_t mesh() const 
    { 
      return this->m_mesh_ptr.lock(); 
    } 

    // Composition of Linear Ops L1(L2( . )) (lval) ---------------------------------------------------
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp) &
    {
      static_assert(traits::is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) & lvalue overload 

    // // composition of linear of L1(L2( . )) (rval)
    template<typename DerivedInner>
    auto compose(DerivedInner&& InnerOp) && 
    {
      static_assert(traits::is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) && rvalue overload  

  private:
    // not accessibles ==============================================================================
    // composition of linear Ops L1(L2( . )) --------------------------------------------------------
    template<typename DerivedInner, typename Lhs_t = Derived_t> 
    auto compose_impl(DerivedInner&& InnerOp)
    {
      static_assert(traits::is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      if constexpr(traits::is_add_expr<std::remove_reference_t<DerivedInner>>::value){
        return compose(InnerOp.Lhs())+compose(InnerOp.Rhs());
      }
      else if constexpr(traits::is_subtraction_expr<std::remove_reference_t<DerivedInner>>::value){
        return compose(InnerOp.Lhs())-compose(InnerOp.Rhs());
      }
      else if constexpr(traits::is_negation_expr<std::remove_reference_t<DerivedInner>>::value){
        return -compose(InnerOp.Lhs());
      }
      else if constexpr(traits::is_scalar_multiply_expr<std::remove_reference_t<DerivedInner>>::value){
        return InnerOp.Lhs() * compose(InnerOp.Rhs()); 
      }
      else{
        using Rhs_t = std::remove_reference_t<DerivedInner>;
        return LinOpExpr<Lhs_t, Rhs_t, internal::linopXlinop_mult_op>(
        std::forward<Lhs_t>(static_cast<Lhs_t>(*this)), // lhs
        std::forward<DerivedInner>(InnerOp), // rhs 
        internal::linopXlinop_mult_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
        ); 
      } // end else 
    }; // end .compose_impl(other) 
    
    // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    auto left_scalar_mult_impl(double c) & {
      return LinOpExpr<double, Derived&, internal::scalar_left_mult_op>(
        c, // lhs scalar
        static_cast<Derived&>(*this), // rhs 
        internal::scalar_left_mult_op{}, // unary_op
        m_mesh_ptr // rhs's mesh 
      );
    }
    
    // Left multiply by a scalar: i.e. c*L (rval)
    auto left_scalar_mult_impl(double c) && {
      return LinOpExpr<double, Derived&&, internal::scalar_left_mult_op>(
        c, // lhs scalar
        static_cast<Derived&&>(*this), // rhs 
        internal::scalar_left_mult_op{}, // unary_op,
        m_mesh_ptr // rhs's mesh
      );
    }

  public:
    // Operators ================================================================================ 
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) & {
        return LinOpExpr<Derived&, LINOP_T, internal::linop_bin_add_op>(
        static_cast<Derived&>(*this),  // lhs
        std::forward<LINOP_T>(rhs), // rhs 
        internal::linop_bin_add_op{}, // bin_op
        m_mesh_ptr
      );
    }
    
    // L1 + L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) && {
        return LinOpExpr<Derived&&, LINOP_T, internal::linop_bin_add_op>(
        static_cast<Derived&&>(*this), // lhs 
        std::forward<LINOP_T>(rhs), // rhs
        internal::linop_bin_add_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
      );
    }
    
    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) & {
        return LinOpExpr<Derived&, LINOP_T, internal::linop_bin_subtract_op>(
        static_cast<Derived&>(*this), //lhs
        std::forward<LINOP_T>(rhs), //rhs
        internal::linop_bin_subtract_op{}, //bin_op
        m_mesh_ptr
      );
    }
    
    // L1 - L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) && {
        return LinOpExpr<Derived&&, LINOP_T, internal::linop_bin_subtract_op>(
        static_cast<Derived&&>(*this), // lhs
        std::forward<LINOP_T>(rhs), // rhs 
        internal::linop_bin_subtract_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
      );
    }
    
    // riend declare c * L (lval + rval) -------------------------------------------
    template<typename LINOP_T, typename>
    friend auto operator*(double scalar, LINOP_T&& rhs); 
    template<typename T>
    friend struct traits::supports_left_scalar_mult;
    // unary operator-() (lval) ---------------------------------------------- 
    auto operator-() & {
        return LinOpExpr<Derived&, void, internal::unary_negate_op>(
        static_cast<Derived&>(*this), // rhs 
        internal::unary_negate_op{}, // unary op 
        m_mesh_ptr // rhs mesh 
      );
    }
    // unary operator-() (rval) 
    auto operator-() && {
        return LinOpExpr<Derived&&, void, internal::unary_negate_op>(
        static_cast<Derived&&>(*this), // rhs
        internal::unary_negate_op{}, // unary_op
        m_mesh_ptr
      );
    }
}; // end LinOpBase

// operator*(c,L) outside of class .... 
template<typename LINOP_T, 
  typename = std::enable_if_t<
    traits::supports_left_scalar_mult<LINOP_T>::value
  >
>
auto operator*(double c, LINOP_T&& rhs){
  return std::forward<LINOP_T>(rhs).left_scalar_mult_impl(c); 
}

} // end namespace LinOps 

#endif // LinearOpBase.hpp