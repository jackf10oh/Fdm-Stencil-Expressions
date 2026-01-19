// TimeDepCoeff.hpp
//
//
//
// JAF 12/26/2025 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include<Eigen/Core> 
#include "../CoeffOpBase.hpp"
#include "../../LinOps/LinOpTraits.hpp" // callable_traits<F> 
#include "../../LinOps/Discretization.hpp"
#include "../../Utilities/SparseDiagExpr.hpp"

namespace Fds{
using namespace LinOps; 

template<typename FUNC_STORAGE_T = std::function<double(double,double)>>
class TimeDepCoeff : public CoeffOpBase<TimeDepCoeff<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = TimeDepCoeff; 
  public:
    FUNC_STORAGE_T m_function;  
    Discretization1D m_diag_vals; 
  public:
    // constructors ==========================================================
    TimeDepCoeff()=delete; // no default constructor
    // from callable + mesh 
    TimeDepCoeff(FUNC_STORAGE_T f_init, MeshPtr_t m = MeshPtr_t{})
      : m_function(f_init), m_diag_vals(0)
    {
      static_assert(LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args > 0, "Assinging functions with no arguments to TimeDepCoeff not allowed"); 
      static_assert(LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args <=2, "In 1D, TimeDepCoeff must have form a(t,x) or a(t)");  
      set_mesh(m);
    }
    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other)=delete; 
    
    // destructors =============================================================
    ~TimeDepCoeff()=default;

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
      // ensure we aren't resetting the mesh again
      if(!this->m_mesh_ptr.owner_before(m) && !m.owner_before(this->m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 

      if(!locked) return; 
      this->m_mesh_ptr = m; // store the mesh  
      // perform work on locked 
      m_diag_vals.resize(m);
    }
    // update state of TimeDepCoeff from a given t 
    void SetTime_impl(double t){
      // SetTime() stores the new time... 
      // Store the result of m_function(t,x)  
      if constexpr(LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args==2){
        auto binded = std::bind(m_function, t, std::placeholders::_1); 
        m_diag_vals.set_init(this->m_mesh_ptr, binded); 
      }
      // set diag to be f(t) 
      if constexpr(LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args==1){
        m_diag_vals.set_init( m_function(t) ); 
      }
    };
}; 

} // end namespace Fds 

#endif // TimeDepCoeff.hpp 
