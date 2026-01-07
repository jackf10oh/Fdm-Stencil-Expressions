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

template<typename Derived>
class CoeffOpBase : public LinOpBase<CoeffOpBase<Derived>>
{
  public:
    // use member types so FdmPlugin can access grandchildren 
    using Derived_t = Derived;
  public:
    // Constructors
    CoeffOpBase(MeshPtr_t m=nullptr){ set_mesh(m); };
    // Member functions ========================================================
    // new functionality: must be implemented ----------------------------------
    double GetScalar() const
    {
      return static_cast<const Derived*>(this)->GetScalar(); 
    }
    // have default implementations --------------------------------------------
    decltype(auto) GetMat()
    {
      auto row_expr = Eigen::RowVectorXd::Ones(this->m_mesh_ptr->size());
      return SparseDiag(GetScalar() * row_expr);  
    };
    decltype(auto) GetMat() const
    {
      auto row_expr = static_cast<const Derived*>(this)->GetScalar() * Eigen::RowVectorXd::Ones(this->m_mesh_ptr->size());
      return SparseDiag(row_expr);  
    };
    void set_mesh(const MeshPtr_t& m)
    {
       // do nothing on nullptr or copy of m_mesh_ptr 
      if(m==nullptr || m==this->m_mesh_ptr) return;
      // just store it. expecting coeff ops to not depend on mesh 
      this->m_mesh_ptr = m; 
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