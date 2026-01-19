// AutonomousCoeff.hpp
//
//
//
// JAF 12/12/2025

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include<Eigen/Core> 
#include "../CoeffOpBase.hpp"
#include "../../LinOps/LinOpTraits.hpp" // callable_traits<F> 
#include "../../LinOps/Discretization.hpp" 
#include "../../Utilities/SparseDiagExpr.hpp"

namespace Fds{
using namespace LinOps; 

template<typename FUNC_STORAGE_T = std::function<double(double,double)>>
class AutonomousCoeff : public CoeffOpBase<AutonomousCoeff<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = AutonomousCoeff; 
  public:
    FUNC_STORAGE_T m_function;  
    Discretization1D m_diag_vals;
  public:
    // constructors ==========================================================
    AutonomousCoeff()=delete; // no default constructor
    // from callable + mesh 
    AutonomousCoeff(FUNC_STORAGE_T f_init, MeshPtr_t m = MeshPtr_t{})
      : m_function(f_init), m_diag_vals(0)
    {
      static_assert(LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args==1, "In 1D, AutonomousCoeff must have form a(x)");  
      set_mesh(m);
    }
    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=default; 
    
    // destructors =============================================================
    ~AutonomousCoeff()=default;

    // member funcs =============================================================
    // Matrix getters 
    auto GetMat()
    {
      return SparseDiag(m_diag_vals.values().transpose());  
    };
    auto GetMat() const
    {
      return SparseDiag(m_diag_vals.values().transpose());
    };
    // updated m_mesh_ptr, resize m_diag_vals
    void set_mesh(MeshPtr_t m){
      m_diag_vals.set_init(m, m_function); 
    }
    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){};
}; 

} // end namespace Fds 

#endif // AutonomousCoeff.hpp 