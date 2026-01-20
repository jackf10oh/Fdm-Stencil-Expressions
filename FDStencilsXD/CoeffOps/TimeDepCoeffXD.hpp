// TimeDepCoeffXD.hpp
//
//
//
// JAF 1/11/2026 

#ifndef TIMEDEPCOEFFXD_H
#define TIMEDEPCOEFFXD_H

#include "../CoeffOpBaseXD.hpp" 
#include "../../LinOpsXD/LinOpTraitsXD.hpp" // callable_traits<F> 
#include "../../LinOpsXD/DiscretizationXD.hpp" 
#include "../../Utilities/SparseDiagExpr.hpp"

namespace Fds{

template<typename FUNC_STORAGE_T>
class TimeDepCoeffXD : public CoeffOpBaseXD<TimeDepCoeffXD<FUNC_STORAGE_T>>
{
  public:
    // Type Defs 
    using Derived_t = TimeDepCoeffXD; 
  public:
    // Member Data -----------------------------------
    FUNC_STORAGE_T m_function;  
    LinOps::DiscretizationXD m_diag_vals;
    std::size_t m_prod_after; 

  public:
    // Constructors + Destructor ==========================================================
    TimeDepCoeffXD()=delete; // no default constructor
    // from callable + mesh 
    TimeDepCoeffXD(FUNC_STORAGE_T f_init, LinOps::MeshXDPtr_t m = LinOps::MeshXDPtr_t{})
      : m_function(f_init), m_diag_vals(0), m_prod_after(std::size_t{1}) 
    {
      static_assert(
        std::is_same<double, typename LinOps::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value,
        "Error constructing coeff: F doesn't return double"
      );  
      set_mesh(m);
      SetTime_impl( this->m_current_time ); 
    }
    // copy constructor
    TimeDepCoeffXD(const TimeDepCoeffXD& other)=delete; 
    
    // destructor 
    ~TimeDepCoeffXD()=default;

    // Member Funcs =============================================================
    // Matrix getters 
    auto GetMat()
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals.values(), m_prod_after); 
    };
    auto GetMat() const
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals.values(), m_prod_after); 
    };
    // updated m_mesh_ptr
    void set_mesh(LinOps::MeshXDPtr_t m){
      this->m_mesh_ptr = m; 
    }
    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){
      auto locked = this->m_mesh_ptr.lock(); 
      if(!locked) return; 
      constexpr std::size_t N = LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args; 
      static_assert(N>=1, "Time Dep Coeff Requires # args >= 1 i.e. a(t, ...)"); 
      if constexpr(N==1){
        m_diag_vals.resize( locked ); 
        m_diag_vals.set_init( m_function(t) ); // m_function(t) evaluates to a double. 
      }
      else{

      std::vector<std::shared_ptr<const LinOps::Mesh1D>> s(N-1);
      for(auto i=0; i<N-1; i++) s[i] = locked->GetMesh(i);  
      auto sub_dims = std::make_shared<const LinOps::MeshXD>(s); 

      using Bind_t = typename LinOps::traits::callable_traits<FUNC_STORAGE_T>::BindFirst_t; 
      Bind_t binded(m_function, t); 
      m_diag_vals.resize( sub_dims ); 
      m_diag_vals.set_init( binded );
      m_prod_after = locked->sizes_product() / m_diag_vals.values().size(); 
      }
    };
}; 

} // end namespace Fds 

#endif // TimeDepCoeffXD.hpp 