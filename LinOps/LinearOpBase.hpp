// LinearOpBase.hpp
//
//
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

#ifndef LINOP_PLUGIN
// use an empty struct if no plugin given.
struct empty{}; 
template<typename T> 
using LinOpMixIn = empty; 
#else
template<typename T>
using LinOpMixIn = LINOP_PLUGIN<T>; 
#endif

#ifndef CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE
#define CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Eigen::SparseMatrix<double, Eigen::ColMajor>
#endif 

// CTRP base for 1D differential operator -------------------------------
template<typename Derived>
class LinOpBase : public LinOpMixIn<LinOpBase<Derived>> 
{
  // friend classes. LINOP_PLUGIN<T> can use private members of T 
  friend LinOpMixIn<LinOpBase>; 

  // type defs 
  public:
    typedef struct{} is_linop_tag; // to tell if a class derived from LinOpBase<> 
    using Derived_t = Derived; // so Plugin can access grand child class

  protected:
    // member data
    MeshPtr_t m_mesh_ptr; // all LinOps keep weak pointers to mesh they operate on
    // member data for matrix kept in derived classes

  public:
    // Constructors ---------------------------------------------------------- 
    // LinOpBase(MeshPtr_t m=nullptr): LinOpMixIn<LinOpBase<Derived>>(){set_mesh(m);}; 
    // member functions -------------------------------------------------------
    // member functions. implemented by derived class -------------------------
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat(); 
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
        Discretization1D result; 
        result = GetMat()*d.values(); // assign to values in result. 
        result.match_mesh(d.mesh()); // set mesh in result to be same as d_arr 
        return result; 
      }
    };

    // fit operator to a mesh of rectangular domain.
    void set_mesh(MeshPtr_t m) 
    {
      // ensure we aren't resetting the mesh again, or setting to nullptr
      // auto locked = m.lock(); 
      // if(!locked) return; // do nothing on nullptr. or throw an error 
      // if(locked == m_mesh_ptr.lock()) return; // do nothing if m,m_mesh_ptr both point to same mesh 
      // m_mesh_ptr = m; // store the mesh  
      // perform work on locked 
      static_cast<Derived*>(this)->set_mesh(m);
    };
    // return reference to stored MeshPtr_t
    MeshPtr_t mesh() const 
    { 
      return m_mesh_ptr; 
    } 

    // Composition of Linear Ops L1(L2( . )) ---------------------------------------------------
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp) &
    {
      static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) & lvalue overload 

    // // composition of linear of L1(L2( . ))
    template<typename DerivedInner>
    auto compose(DerivedInner&& InnerOp) && 
    {
      static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,Derived_t&&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) && rvalue overload  

  private:
    // not accessibles --------------------------------------------------------------------------------------------
    // composition of linear Ops L1(L2( . ))
    template<typename DerivedInner, typename Lhs_t = Derived_t> 
    auto compose_impl(DerivedInner&& InnerOp)
    {
      static_assert(is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      if constexpr(is_add_expr<std::remove_reference_t<DerivedInner>>::value){
        return compose(InnerOp.Lhs())+compose(InnerOp.Rhs());
      }
      else if constexpr(is_subtraction_expr<std::remove_reference_t<DerivedInner>>::value){
        return compose(InnerOp.Lhs())-compose(InnerOp.Rhs());
      }
      else if constexpr(is_negation_expr<std::remove_reference_t<DerivedInner>>::value){
        return -compose(InnerOp.Lhs());
      }
      else if constexpr(is_scalar_multiply_expr<std::remove_reference_t<DerivedInner>>::value){
        return InnerOp.Lhs() * compose(InnerOp.Rhs()); 
      }
      else{
        using Rhs_t = std::remove_reference_t<DerivedInner>;
        return LinOpExpr<Lhs_t, Rhs_t, linopXlinop_mult_op>(
        std::forward<Lhs_t>(static_cast<Lhs_t>(*this)), 
        std::forward<DerivedInner>(InnerOp), 
        linopXlinop_mult_op{}
        ); 
      } // end else 
    }; // end .compose_impl(other) 
    // Left multiply by a scalar: i.e. c*L ------------------------------------------------------------------------- 
    auto left_scalar_mult_impl(double c) & {
      return LinOpExpr<double, Derived&, scalar_left_mult_op>(
        c,
        static_cast<Derived&>(*this), 
        scalar_left_mult_op{}
      );
    }
    // Left multiply by a scalar: i.e. c*L (rval)
    auto left_scalar_mult_impl(double c) && {
      return LinOpExpr<double, Derived&&, scalar_left_mult_op>(
        c,
        static_cast<Derived&&>(*this), 
        scalar_left_mult_op{}
      );
    }

  private:  
    // Structs for binary operations f(L1,L2) to get matrix of expression
    // L1 + L2 
    struct linop_bin_add_op : public OperatorAddition_t
    {
      template<typename L1, typename L2>
      auto operator()(const L1& A, const L2& B) const { return (A.GetMat()) + (B.GetMat()); }
    }; 
    // L1 - L2 
    struct linop_bin_subtract_op : public OperatorSubtraction_t
    {
      template<typename L1, typename L2>
      auto operator()(const L1& A, const L2& B) const { return (A.GetMat()) - (B.GetMat()); }
    }; 
    // c * L
    struct scalar_left_mult_op : public ScalarMultiply_t
    {
      template<typename L2>
      auto operator()(const double& c, const L2& B) const { return  c*(B.GetMat()); }
    }; 
    // -L 
    struct unary_negate_op : public OperatorNegation_t
    {
      template<typename L1>
      auto operator()(const L1& B) const { return  -(B.GetMat()); }
    }; 
    // composition: L1( L2( . ) )
    struct linopXlinop_mult_op : public OperatorComposition_t
    {
      template<typename L1, typename L2>
      auto operator()(const L1& A, const L2& B) const { return  (A.GetMat())*(B.GetMat()); }
    }; 
    
  public:
    // Operators ================================================================================ 
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) & {
        return LinOpExpr<Derived&, LINOP_T, linop_bin_add_op>(
        static_cast<Derived&>(*this), 
        std::forward<LINOP_T>(rhs),
        linop_bin_add_op{}
      );
    }
    // L1 + L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) && {
        return LinOpExpr<Derived&&, LINOP_T, linop_bin_add_op>(
        static_cast<Derived&&>(*this), 
        std::forward<LINOP_T>(rhs),
        linop_bin_add_op{}
      );
    }
    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) & {
        return LinOpExpr<Derived&, LINOP_T, linop_bin_subtract_op>(
        static_cast<Derived&>(*this), 
        std::forward<LINOP_T>(rhs),
        linop_bin_subtract_op{}
      );
    }
    // L1 - L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) && {
        return LinOpExpr<Derived&&, LINOP_T, linop_bin_subtract_op>(
        static_cast<Derived&&>(*this), 
        std::forward<LINOP_T>(rhs),
        linop_bin_subtract_op{}
      );
    }
    // riend declare c * L (lval + rval) 
    template<typename LINOP_T, typename>
    friend auto operator*(double scalar, LINOP_T&& rhs); 

    // unary operator-() (lval) ---------------------------------------------- 
    auto operator-() & {
        return LinOpExpr<Derived&, void, unary_negate_op>(
        static_cast<Derived&>(*this), 
        unary_negate_op{}
      );
    }
    // unary operator-() (rval) 
    auto operator-() && {
        return LinOpExpr<Derived&&, void, unary_negate_op>(
        static_cast<Derived&&>(*this), 
        unary_negate_op{}
      );
    }
  }; // end LinOpBase

template<typename LINOP_T, typename = std::enable_if_t<is_linop_crtp<LINOP_T>::value>>
auto operator*(double c, LINOP_T&& rhs){
  return std::forward<LINOP_T>(rhs).left_scalar_mult_impl(c); 
}

#endif // LinearOpBase.hpp