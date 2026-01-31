// CoeffOpMixIn.hpp
//
// going to need some way to model equations like 
// Ut = c(t,x) * Uxx + c(t,x) * Ux
//
// CoeffOpMixIn.hpp

#ifndef COEFFOP_H
#define COEFFOP_H 

#include "LinearOpBase.hpp"

namespace LinOps{

template<typename DERIVED>
class CoeffOpMixIn : public LinOpMixIn<CoeffOpMixIn<DERIVED>>
{

  public:
    // Type Defs --------------------------------
    // So LinOpMixIn case access DERIVED type of grandchild classes
    using DERIVED_T = DERIVED; 
    // flag type for any derived class to be picked up by is_coeffop_crtp<>::value trait 
    struct is_coeff_flag{}; 

    // Operators ================================================================== 
    // composition c * L  produces composition
    template<typename LINOP_T>
    auto operator*(LINOP_T&& rhs)
    {
      static_assert(!traits::is_coeffop_crtp<LINOP_T>::value,"Coefficients are meant to multiply c*L for L linear operator. not another Coefficient. a*b*L should be written as 1 functions");
      return LinOpMixIn<CoeffOpMixIn<DERIVED>>::compose(std::forward<LINOP_T>(rhs));
    };
    
    // delting a ton of operators out of LinOpMixIn =====================================
    // so that you can't make expressions of coefficients
    auto operator-()=delete;
    template<typename RHS>
    auto operator-(RHS&& rhs)=delete;
    template<typename RHS>
    auto operator+(RHS&& rhs)=delete; 
    template<typename RHS>
    auto compose(RHS&& InnerOp)=delete; 
    auto left_scalar_mult_impl(double c)=delete; 
};

} // end namespace LinOps 

#endif // CoeffOpMixIn.hpp 