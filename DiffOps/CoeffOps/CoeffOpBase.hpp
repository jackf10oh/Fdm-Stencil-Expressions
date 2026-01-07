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
    // Constructors
    CoeffOpBase(MeshPtr_t m=nullptr){ set_mesh(m); };
    // Member functions ========================================================
    // have default implementations --------------------------------------------
    decltype(auto) GetMat()
    {
      return static_cast<Derived*>(this)->GetMat();
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const Derived*>(this)->GetMat();
    };
    void set_mesh(const MeshPtr_t& m)
    {
      static_cast<Derived*>(this)->set_mesh(m); 
      // if(m==nullptr || m==this->m_mesh_ptr) return; // do nothing on nullptr or copy of m_mesh_ptr 
      // this->m_mesh_ptr = m; // just store it. expecting coeff ops to not depend on mesh 
    }
    // must be implemented by derived classes -----------------------------------
    void SetTime_impl(double t)
    {
      static_cast<Derived*>(this)->SetTime_impl(t); 
    }
    // composition c * L ( Lval overload)
    template<typename DerivedR>
    auto operator*(DerivedR&& RHS) &
    {
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<DerivedR>(RHS));
    };
    // operator for scalar multiplication c * L ( Rval overload)
    template<typename DerivedR>
    auto operator*(DerivedR&& RHS) &&
    {
      return LinOpBase<CoeffOpBase<Derived>>::compose(std::forward<DerivedR>(RHS));
    };
};

// operator+ should be deleted for coeffopbase LHS/RHS... 

#endif