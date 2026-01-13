// AutonomousCoeffXD.hpp
//
//
//
// JAF 1/11/2026 

#ifndef AUTONOMOUSCOEFFXD_H
#define AUTONOMOUSCOEFFXD_H

#include "../CoeffOpBaseXD.hpp" 
#include "../../LinOpsXD/LinOpTraitsXD.hpp" // callable_traits<F> 
#include "../../LinOpsXD/DiscretizationXD.hpp" 
#include "../../Utilities/SparseDiagExpr.hpp"

namespace Fds{
using namespace LinOps; 

template<typename FUNC_STORAGE_T>
class AutonomousCoeffXD : public CoeffOpBaseXD<AutonomousCoeffXD<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = AutonomousCoeffXD; 
    using RowVec_t = decltype(Eigen::VectorXd().transpose()); 
  public:
    FUNC_STORAGE_T m_function;  
    DiscretizationXD m_diag_vals;
    std::size_t m_prod_after; 
  public:
    // constructors ==========================================================
    AutonomousCoeffXD()=delete; // no default constructor
    // from callable + mesh 
    AutonomousCoeffXD(FUNC_STORAGE_T f_init, MeshXDPtr_t m = MeshXDPtr_t{})
      : m_function(f_init), m_diag_vals(0), m_prod_after(std::size_t{1}) 
    {
      static_assert(
        std::is_same<double, typename internal::callable_traits<FUNC_STORAGE_T>::result_type>::value,
        "Error constructing coeff: F doesn't return double"
      );  
      set_mesh(m);
    }
    // copy constructor
    AutonomousCoeffXD(const AutonomousCoeffXD& other)=delete; 
    
    // destructors =============================================================
    ~AutonomousCoeffXD()=default;

    // member funcs =============================================================
    // Matrix getters 
    auto GetMat()
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals.values(), m_prod_after); 
    };
    auto GetMat() const
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals.values(), m_prod_after); 
    };
    // updated m_mesh_ptr, resize m_diag_vals
    void set_mesh(MeshXDPtr_t m){
      auto locked = m.lock(); 
      if(!locked) return; 
      constexpr std::size_t N = internal::callable_traits<FUNC_STORAGE_T>::num_args; 

      if constexpr(N==0){
        m_diag_vals.set_init( m, m_function() ); 
      }
      else{
      std::vector<std::shared_ptr<const Mesh1D>> s(N);
      for(auto i=0; i<N; i++) s[i] = locked->GetMesh(i);  
      auto sub_dims = std::make_shared<const MeshXD>(s); 
      m_diag_vals.set_init(sub_dims, m_function); 
      m_prod_after = locked->sizes_product() / m_diag_vals.values().size(); 
      }
    }
    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){};
}; 

} // end namespace Fds 

#endif // AutonomousCoeffXD.hpp 