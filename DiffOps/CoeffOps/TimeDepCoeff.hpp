// TimeDepCoeff.hpp
//
//
//
// JAF 12/26/2025 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include<Eigen/Core> 
#include "CoeffOpBase.hpp"
#include "../Traits.hpp"
#include "../../Utilities/SparseDiagExpr.hpp"

template<typename FUNC_STORAGE_T = std::function<double(double,double)>>
class TimeDepCoeff : public CoeffOpBase<TimeDepCoeff<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = TimeDepCoeff; 
  public:
    FUNC_STORAGE_T m_function;  
    Eigen::RowVectorXd m_diag_vals; 
  public:
    // constructors ==========================================================
    TimeDepCoeff()=delete; // no default constructor
    // from callable + mesh 
    TimeDepCoeff(FUNC_STORAGE_T f_init, MeshPtr_t m=nullptr)
      : m_function(f_init), m_diag_vals(0)
    {
      static_assert(callable_traits<FUNC_STORAGE_T>::num_args > 0, "Assinging functions with no arguments to TimeDepCoeff not allowed"); 
      static_assert(callable_traits<FUNC_STORAGE_T>::num_args <=2, "In 1D, TimeDepCoeff must have form a(t,x) or a(t)");  
      set_mesh(m);
    }
    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other)=delete; 
    // destructors =============================================================
    ~TimeDepCoeff()=default;
    // member funcs =============================================================
    auto GetMat()
    {
      return SparseDiag(m_diag_vals);  
    };
    auto GetMat() const
    {
      return SparseDiag(m_diag_vals);
    };
    void set_mesh(const MeshPtr_t& m){
      if(m==nullptr || m==this->m_mesh_ptr) return; 
      m_diag_vals.resize(m->size());
      this->m_mesh_ptr = m;  
    }
    // update state of TimeDepCoeff from a given t 
    void SetTime_impl(double t){
      // SetTime() stores the new time... 
      // Store the result of m_function(t,x)  
      // m_val = m_function(t); 
      if constexpr(callable_traits<FUNC_STORAGE_T>::num_args==2){
        for(std::size_t i=0; i<m_diag_vals.size(); i++){
          m_diag_vals[i] = m_function(t, this->m_mesh_ptr->operator[](i)); 
        }
      }
      if constexpr(callable_traits<FUNC_STORAGE_T>::num_args==1){
        m_diag_vals.setConstant( m_function(t) ); 
      }
    };
}; 

#endif // TimeDepCoeff.hpp 
