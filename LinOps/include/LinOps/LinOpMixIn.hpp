// LinOpMixIn.hpp 
//
// CRTP mixin class for both LinOpBase1D + XD
//
// JAF 1/29/2026

#ifndef LINOPMIXIN_H
#define LINOPMIXIN_H

#include "LinOpTraits.hpp"
#include "LinOpExpr.hpp"

namespace LinOps{

template<typename DERIVED>
class LinOpMixIn
{
  // optional macro to add features into both LinOpBase1D<> + XD<>
  #ifndef LINOP_PLUGIN
  #else
  #include LINOP_PLUGIN
  #endif

  public:
    // Type Defs -------------------------------------
    typedef struct{} is_linop_tag; // to tell if a class derived from LinOpMixIn<> 

    // Member Funcs ==========================================================

    // SetTime() + SetTime_impl() ... 

    // Getters to convert a Linear Operator to it Matrix form -----------------------
    decltype(auto) GetMat()
    {
      return static_cast<DERIVED*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const DERIVED*>(this)->GetMat(); 
    };

    // Composition of Linear Ops L1(L2( . )) (lval) ---------------------------------------------------
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp) &
    {
      static_assert(traits::is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,DERIVED&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) & lvalue overload 

    // // composition of linear of L1(L2( . )) (rval)
    template<typename DerivedInner>
    auto compose(DerivedInner&& InnerOp) && 
    {
      static_assert(traits::is_linop_crtp<DerivedInner>::value, "compose only works on other linear operators!"); 
      return compose_impl<DerivedInner,DERIVED&&>(std::forward<DerivedInner>(InnerOp)); 
    }; // end .compose(other) && rvalue overload  

  private:
    // not accessibles ==============================================================================
    // composition of linear Ops L1(L2( . )) --------------------------------------------------------
    template<typename DerivedInner, typename Lhs_t> 
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
        internal::linopXlinop_mult_op{} // bin_op
        ); 
      } // end else 
    }; // end .compose_impl(other) 
    
    // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    template<typename SCALAR_T>
    auto left_scalar_mult_impl(SCALAR_T&& c) & {
      return LinOpExpr<SCALAR_T, DERIVED&, internal::scalar_left_mult_op>(
        std::forward<SCALAR_T>(c), // lhs scalar
        static_cast<DERIVED&>(*this), // rhs 
        internal::scalar_left_mult_op{} // unary_op
      );
    }
    
    // Left multiply by a scalar: i.e. c*L (rval)
    template<typename SCALAR_T>
    auto left_scalar_mult_impl(SCALAR_T&& c) && {
      return LinOpExpr<SCALAR_T, DERIVED&&, internal::scalar_left_mult_op>(
        std::forward<SCALAR_T>(c), // lhs scalar
        static_cast<DERIVED&&>(*this), // rhs 
        internal::scalar_left_mult_op{} // unary_op,
      );
    }

  public:
    // Operators ================================================================================ 
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) & {
        return LinOpExpr<DERIVED&, LINOP_T, internal::linop_bin_add_op>(
        static_cast<DERIVED&>(*this),  // lhs
        std::forward<LINOP_T>(rhs), // rhs 
        internal::linop_bin_add_op{} // bin_op
      );
    }
    
    // L1 + L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator+(LINOP_T&& rhs) && {
        return LinOpExpr<DERIVED&&, LINOP_T, internal::linop_bin_add_op>(
        static_cast<DERIVED&&>(*this), // lhs 
        std::forward<LINOP_T>(rhs), // rhs
        internal::linop_bin_add_op{} // bin_op
      );
    }
    
    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) & {
        return LinOpExpr<DERIVED&, LINOP_T, internal::linop_bin_subtract_op>(
        static_cast<DERIVED&>(*this), //lhs
        std::forward<LINOP_T>(rhs), //rhs
        internal::linop_bin_subtract_op{} //bin_op
      );
    }
    
    // L1 - L2 (rval)  
    template<typename LINOP_T, typename = std::enable_if_t<traits::is_linop_crtp<LINOP_T>::value>>
    auto operator-(LINOP_T&& rhs) && {
        return LinOpExpr<DERIVED&&, LINOP_T, internal::linop_bin_subtract_op>(
        static_cast<DERIVED&&>(*this), // lhs
        std::forward<LINOP_T>(rhs), // rhs 
        internal::linop_bin_subtract_op{} // bin_op
      );
    }
    
    // friend declare c * L (lval + rval) -------------------------------------------
    template<typename SCALAR_T, typename LINOP_T, typename>
    friend auto operator*(SCALAR_T&& scalar, LINOP_T&& rhs); 
    template<typename T>
    friend struct traits::supports_left_scalar_mult;

    // unary operator-() (lval) ---------------------------------------------- 
    auto operator-() & {
        return LinOpExpr<DERIVED&, void, internal::unary_negate_op>(
        static_cast<DERIVED&>(*this), // lhs 
        internal::unary_negate_op{} // unary op 
      );
    }
    // unary operator-() (rval) 
    auto operator-() && {
        return LinOpExpr<DERIVED&&, void, internal::unary_negate_op>(
        static_cast<DERIVED&&>(*this), // lhs
        internal::unary_negate_op{} // unary_op
      );
    }

}; // end LinOpMixIn<> 

// operator*(c,L) outside of class ....
template<
  typename SCALAR_T, 
  typename LINOP_T, 
  typename = std::enable_if_t<
    std::conjunction_v<
    traits::is_linop_crtp<LINOP_T>, 
    traits::supports_left_scalar_mult<LINOP_T>, 
    std::is_arithmetic<std::remove_cv_t<std::remove_reference_t<SCALAR_T>>>
    > // end conjuntion 
  > // end enable_if
>
auto operator*(SCALAR_T&& c, LINOP_T&& rhs){
  return std::forward<LINOP_T>(rhs).left_scalar_mult_impl( std::forward<SCALAR_T>(c) ); 
}

} // end namespace LinOps 


#endif
