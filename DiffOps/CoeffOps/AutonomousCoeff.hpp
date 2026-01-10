// AutonomousCoeff.hpp
//
//
//
// JAF 12/12/2025

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include<Eigen/Core> 
#include "CoeffOpBase.hpp"
#include "../../LinOps/LinOpTraits.hpp" // callable_traits<F> 
#include "../../Utilities/SparseDiagExpr.hpp"

// AutonomousCoeff

template<typename FUNC_STORAGE_T = std::function<double(double,double)>>
class AutonomousCoeff : public CoeffOpBase<AutonomousCoeff<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = AutonomousCoeff; 
  public:
    FUNC_STORAGE_T m_function;  
    Eigen::RowVectorXd m_diag_vals; 
  public:
    // constructors ==========================================================
    AutonomousCoeff()=delete; // no default constructor
    // from callable + mesh 
    AutonomousCoeff(FUNC_STORAGE_T f_init, MeshPtr_t m = MeshPtr_t{})
      : m_function(f_init), m_diag_vals(0)
    {
      static_assert(callable_traits<FUNC_STORAGE_T>::num_args==1, "In 1D, AutonomousCoeff must have form a(x)");  
      set_mesh(m);
    }
    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=delete; 
    
    // destructors =============================================================
    ~AutonomousCoeff()=default;

    // member funcs =============================================================
    // Matrix getters 
    auto GetMat()
    {
      return SparseDiag(m_diag_vals);  
    };
    auto GetMat() const
    {
      return SparseDiag(m_diag_vals);
    };
    // updated m_mesh_ptr, resize m_diag_vals
    void set_mesh(MeshPtr_t m){
      // ensure we aren't resetting the mesh again
      if(!this->m_mesh_ptr.owner_before(m) && !m.owner_before(this->m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 

      if(!locked) return; 
      this->m_mesh_ptr = m; // store the mesh 

      m_diag_vals.resize(locked->size()); 
      const auto data = locked->cbegin(); 
      for(std::size_t i=0; i<m_diag_vals.size(); i++){
        m_diag_vals[i] = m_function(data[i]); 
      }
    }
    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){};
}; 

#endif // AutonomousCoeff.hpp 