// CoeffOp.hpp
//
// going to need some way to model equations like 
// Ut = c(t,x) * Uxx + c(t,x) * Ux
//
// CoeffOp.hpp

#ifndef COEFFOP_H
#define COEFFOP_H 

#include "../FdmPlugin.hpp"
#include "../../LinOps/LinearOpBase.hpp"
#include "../../Utilities/SparseDiagExpr.hpp"

// TRAITS ===============================================
// Given a type, see if it is derived from CoeffOpBase's crtp scheme ------------------------------------
template<typename T, typename = void> 
struct is_coeffop_crtp : std::false_type{}; 

template<typename T>
struct is_coeffop_crtp<T, std::void_t<typename std::remove_cv_t<std::remove_reference_t<T>>::is_coeff_flag>>: std::true_type{}; 

// BASE INTERFACE ===============================================
template<typename Derived>
class CoeffOpBase : public LinOpBase<CoeffOpBase<Derived>>
{
  public:
    // flag type for any derived class to be picked up by is_coeffop_crtp<>::value trait 
    struct is_coeff_flag{}; 
    // use member types so FdmPlugin can access grandchildren 
    using Derived_t = Derived;
  public:
    // Constructors / Destructor =================================================== 
    CoeffOpBase()=default;
    CoeffOpBase(const CoeffOpBase& other) = default; 
    ~CoeffOpBase()=default; 

    // Member functions ========================================================
    // must be implemented by derived classes -----------------------------------
    // Get underlying matrix.............. 
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat();
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat();
    };
    // set stencil to new mesh............
    void set_mesh(MeshPtr_t m)
    {
      static_cast<Derived*>(this)->set_mesh(m); 
      // if(m==nullptr || m==this->m_mesh_ptr) return; // do nothing on nullptr or copy of m_mesh_ptr 
      // this->m_mesh_ptr = m; // just store it. expecting coeff ops to not depend on mesh 
    }
    // set stencil to new time..............
    void SetTime_impl(double t)
    {
      static_cast<Derived*>(this)->SetTime_impl(t); 
    }
    
    // Operators ================================================================== 
    // composition c * L ( Lval overload)
    template<typename LINOP_T>
    auto operator*(LINOP_T&& rhs) &
    {
      static_assert(!is_coeffop_crtp<LINOP_T>::value,"Coefficients are meant to multiply c*L for L linear operator. not another Coefficient. a*b*L should be written as 1 functions");
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<LINOP_T>(rhs));
    };
    // operator for scalar multiplication c * L ( Rval overload)
    template<typename LINOP_T>
    auto operator*(LINOP_T&& rhs) &&
    {
      static_assert(!is_coeffop_crtp<LINOP_T>::value,"Coefficients are meant to multiply c*L for L linear operator. not another Coefficient. a*b*L should be written as 1 functions");
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<LINOP_T>(rhs));
    };
    
    // delting a ton of operators out of LinOpBase =====================================
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

#endif // CoeffOpBase.hpp 