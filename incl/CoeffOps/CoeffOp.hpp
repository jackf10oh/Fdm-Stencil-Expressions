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
  // use member types so FdmPlugin can access grandchildren 
  public:
    using Derived_t = Derived; 

  // operator for scalar multiplication c*L
  template<typename DerivedR>
  auto operator*(DerivedR&& RHS)
  {
    return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<DerivedR>(RHS));
  };
};

#endif