// CoeffOp.hpp
//
// going to need some way to model equations like 
// Ut = c(t,x) * Uxx + c(t,x) * Ux
//
// CoeffOp.hpp

#include "../FdmPlugin.hpp"
#include "../../LinOps/LinearOpBase.hpp"

#ifndef COEFFOP_H
#define COEFFOP_H 

template<typename Derived>
class CoeffOpBase : public LinOpBase<CoeffOpBase<Derived>>
{
  public:
    // use member types so FdmPlugin can access grandchildren 
    using Derived_t = Derived;
  public:
    // Constructors
    CoeffOpBase(): LinOpBase<CoeffOpBase<Derived>>(){}; 
    void SetTime_impl(double t)
    {
      static_cast<Derived*>(this)->SetTime_impl(t); 
    }
    auto& GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    const auto& GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat(); 
    };
    // operator for scalar multiplication c*L
    template<typename DerivedR>
    auto operator*(DerivedR&& RHS) &
    {
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<DerivedR>(RHS));
    };
    // operator for scalar multiplication c*L
    template<typename DerivedR>
    auto operator*(DerivedR&& RHS) &&
    {
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<DerivedR>(RHS));
    };
};

#endif