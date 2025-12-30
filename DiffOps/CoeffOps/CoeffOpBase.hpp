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
    CoeffOpBase(MeshPtr_t m=nullptr): LinOpBase<CoeffOpBase<Derived>>(m){}; 
    void SetTime_impl(double t)
    {
      static_cast<Derived*>(this)->SetTime_impl(t); 
    }
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat(); 
    };
    // default functionality to cast to Eigen::VectorXd
    Eigen::VectorXd GetDiag() const {return Eigen::VectorXd(GetMat().diagonal());};
    // Discretization1D apply(const Discretization1D& d) const
    // {
    //   return static_cast<const Derived*>(this)->apply(d); 
    // };
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

// extend is_linop_crtp trait to include any class derived from CoeffOpBase
template<typename T>
struct is_linop_crtp_impl<T, std::enable_if_t<std::is_base_of_v<LinOpBase<CoeffOpBase<T>>,T>,void>> : std::true_type {};

#endif