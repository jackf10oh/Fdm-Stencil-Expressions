// LinearOpBaseXD.hpp
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
#include "LinOpTraitsXD.hpp"

namespace LinOps{

#ifndef LINOPXD_PLUGIN
// use an empty struct if no plugin given.
struct emptyXD{}; 
template<typename T> 
using LinOpXDMixIn = emptyXD; 
#else
template<typename T>
using LinOpXDMixIn = LINOPXD_PLUGIN<T>; 
#endif

#ifndef CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE
#define CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Eigen::SparseMatrix<double, Eigen::RowMajor>
#endif 


template<typename Derived>
class LinOpBaseXD : public LinOpXDMixIn<LinOpBaseXD<Derived>> 
{
  // friend classes. LINOP_PLUGIN<T> can use private members of T 
  // friend LinOpXDMixIn<LinOpBase>; 

  // Type Defs -----------------------------------------------------
  public:
    typedef struct{} is_linopxd_tag; // to tell if a class derived from LinOpBase<> 
    using Derived_t = Derived; // so Plugin can access grand child class

  protected:
    // Member Data -------------------------------------------------
    MeshXDPtr_t m_mesh_ptr; // all LinOps keep weak pointers to mesh they operate on
    // member data for matrix kept in derived classes

  public:
    // Constructors / destructors ============================================================ 
    LinOpBaseXD()=default; 
    LinOpBaseXD(MeshXDPtr_t m) : m_mesh_ptr(m){}; 
    LinOpBaseXD(const LinOpBaseXD& other)=default; 
    ~LinOpBaseXD()=default; 
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
    DiscretizationXD apply(const DiscretizationXD& d) const 
    {
      if constexpr (traits::has_apply<LinOpXDMixIn<LinOpBaseXD<Derived>>>::value)
      {
        return LinOpXDMixIn<LinOpBaseXD<Derived>>::apply(d);
      }
      else
      {
        MeshXDPtr_t m = d.mesh(); 
        if(m_mesh_ptr.owner_before(m) || m.owner_before(m_mesh_ptr)){
          throw std::runtime_error("Linear Operator L and discretization d point to different Mesh1D!");
        }
        Eigen::VectorXd v = GetMat() * d.values();  // calculate A*b
        DiscretizationXD result = std::move(v); // move A*b into result's values 
        result.set_mesh(m); // make result point to same mesh as d 
        return result;
      }
    };

    // fit operator to a mesh of rectangular domain ----------------------------------------
    void set_mesh(MeshXDPtr_t m) 
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
    
    // return reference to stored MeshPtr_t ----------------------------------------------
    MeshXDPtr_t mesh() const 
    { 
      return this->m_mesh_ptr; 
    } 

    // Composition of Linear Ops L1(L2( . )) (lval) ---------------------------------------------------
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp) &
    {
      static_assert(traits::is_linopxd_crtp<DerivedInner>::value, "compose only works on other  XD linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) & lvalue overload 

    // // composition of linear of L1(L2( . )) (rval)
    template<typename DerivedInner>
    auto compose(DerivedInner&& InnerOp) && 
    {
      static_assert(traits::is_linopxd_crtp<DerivedInner>::value, "compose only works on other XD linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) && rvalue overload  

  private:
    // not accessibles ==============================================================================
    // composition of linear Ops L1(L2( . )) --------------------------------------------------------
    template<typename DerivedInner, typename Lhs_t = Derived_t> 
    auto compose_impl(DerivedInner&& InnerOp)
    {
      static_assert(traits::is_linopxd_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
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
        return LinOpExprXD<Lhs_t, Rhs_t, internal::linopXlinop_mult_op>(
        std::forward<Lhs_t>(static_cast<Lhs_t>(*this)), // lhs
        std::forward<DerivedInner>(InnerOp), // rhs 
        internal::linopXlinop_mult_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
        ); 
      } // end else 
    }; // end .compose_impl(other) 
    
    // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    auto left_scalar_mult_impl(double c) & {
      return LinOpExprXD<double, Derived&, internal::scalar_left_mult_op>(
        c, // lhs scalar
        static_cast<Derived&>(*this), // rhs 
        internal::scalar_left_mult_op{}, // unary_op
        m_mesh_ptr // rhs's mesh 
      );
    }
    
    // Left multiply by a scalar: i.e. c*L (rval)
    auto left_scalar_mult_impl(double c) && {
      return LinOpExprXD<double, Derived&&, internal::scalar_left_mult_op>(
        c, // lhs scalar
        static_cast<Derived&&>(*this), // rhs 
        internal::scalar_left_mult_op{}, // unary_op,
        m_mesh_ptr // rhs's mesh
      );
    }

  public:
    // Operators ================================================================================ 
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename LINOPXD_T, typename = std::enable_if_t<traits::is_linopxd_crtp<LINOPXD_T>::value>>
    auto operator+(LINOPXD_T&& rhs) & {
        return LinOpExprXD<Derived&, LINOPXD_T, internal::linop_bin_add_op>(
        static_cast<Derived&>(*this),  // lhs
        std::forward<LINOPXD_T>(rhs), // rhs 
        internal::linop_bin_add_op{}, // bin_op
        m_mesh_ptr
      );
    }
    
    // L1 + L2 (rval)  
    template<typename LINOPXD_T, typename = std::enable_if_t<traits::is_linopxd_crtp<LINOPXD_T>::value>>
    auto operator+(LINOPXD_T&& rhs) && {
        return LinOpExprXD<Derived&&, LINOPXD_T, internal::linop_bin_add_op>(
        static_cast<Derived&&>(*this), // lhs 
        std::forward<LINOPXD_T>(rhs), // rhs
        internal::linop_bin_add_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
      );
    }
    
    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename LINOPXD_T, typename = std::enable_if_t<traits::is_linopxd_crtp<LINOPXD_T>::value>>
    auto operator-(LINOPXD_T&& rhs) & {
        return LinOpExprXD<Derived&, LINOPXD_T, internal::linop_bin_subtract_op>(
        static_cast<Derived&>(*this), //lhs
        std::forward<LINOPXD_T>(rhs), //rhs
        internal::linop_bin_subtract_op{}, //bin_op
        m_mesh_ptr
      );
    }
    
    // L1 - L2 (rval)  
    template<typename LINOPXD_T, typename = std::enable_if_t<traits::is_linopxd_crtp<LINOPXD_T>::value>>
    auto operator-(LINOPXD_T&& rhs) && {
        return LinOpExprXD<Derived&&, LINOPXD_T, internal::linop_bin_subtract_op>(
        static_cast<Derived&&>(*this), // lhs
        std::forward<LINOPXD_T>(rhs), // rhs 
        internal::linop_bin_subtract_op{}, // bin_op
        m_mesh_ptr // lhs's mesh
      );
    }
    
    // riend declare c * L (lval + rval) -------------------------------------------
    template<typename LINOPXD_T, typename>
    friend auto operator*(double scalar, LINOPXD_T&& rhs); 
    template<typename T>
    friend struct traits::supports_left_scalar_mult;

    // unary operator-() (lval) ---------------------------------------------- 
    auto operator-() & {
        return LinOpExprXD<Derived&, void, internal::unary_negate_op>(
        static_cast<Derived&>(*this), // rhs 
        internal::unary_negate_op{}, // unary op 
        m_mesh_ptr // rhs mesh 
      );
    }
    // unary operator-() (rval) 
    auto operator-() && {
        return LinOpExprXD<Derived&&, void, internal::unary_negate_op>(
        static_cast<Derived&&>(*this), // rhs
        internal::unary_negate_op{}, // unary_op
        m_mesh_ptr
      );
    }
}; // end LinOpBase

#ifndef LINEAROPBASE_H
// This operators already defined in LinOps/LinearOpBase.hpp
// this is a pretty ugly macro solution for only defining it once 
// operator*(c,L) outside of class .... 
template<typename LINOP_T, 
  typename = std::enable_if_t<
    traits::supports_left_scalar_mult<LINOP_T>::value
  >
>
auto operator*(double c, LINOP_T&& rhs){
  return std::forward<LINOP_T>(rhs).left_scalar_mult_impl(c); 
}
#endif

} // end namespace LinOps 

#endif // LinearOpBaseXD.hpp