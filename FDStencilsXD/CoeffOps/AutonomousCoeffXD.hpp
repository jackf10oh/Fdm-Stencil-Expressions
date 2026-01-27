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

template<typename FUNC_STORAGE_T>
class AutonomousCoeffXD : public CoeffOpBaseXD<AutonomousCoeffXD<FUNC_STORAGE_T>>
{
  public:
    using Derived_t = AutonomousCoeffXD; 
    using RowVec_t = decltype(Eigen::VectorXd().transpose()); 
  public:
    FUNC_STORAGE_T m_function;  
    LinOps::DiscretizationXD m_diag_vals;
    std::size_t m_prod_after; 
  public:
    // Constructors + Destructor ==========================================================
    AutonomousCoeffXD()=delete; // no default constructor
    // from callable + mesh 
    AutonomousCoeffXD(FUNC_STORAGE_T f_init, const LinOps::MeshXD_SPtr_t& m = nullptr)
      : m_function(f_init), m_diag_vals(0), m_prod_after(std::size_t{1}) 
    {
      static_assert(
        std::is_same<double, typename LinOps::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value,
        "Error constructing coeff: F doesn't return double"
      );  
      if(m) set_mesh(m);
    }
    // copy constructor
    AutonomousCoeffXD(const AutonomousCoeffXD& other)=delete; 
    
    // destructor
    ~AutonomousCoeffXD()=default;

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
    // updated m_mesh_ptr, resize m_diag_vals
    void set_mesh(const LinOps::MeshXD_SPtr_t& m){

      constexpr std::size_t N = LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args; 

      if constexpr(N==0){
        m_diag_vals.set_init( m, m_function() ); 
      }
      else{
        std::vector<Mesh1D_SPtr_t> s(N);
        for(auto i=0; i<N; i++) s[i] = m->GetMesh(i);  
        auto sub_dims = std::make_shared<const LinOps::MeshXD>(std::move(s)); 
        m_diag_vals.set_init(sub_dims, m_function); 
        m_prod_after = m->sizes_product() / m_diag_vals.values().size(); 
      }
    }
    // update state of AutonomousCoeff from a given t 
    void SetTime_impl(double t){};
}; 

} // end namespace Fds 

#endif // AutonomousCoeffXD.hpp 