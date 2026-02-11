// AutonomousCoeff.hpp
//
//
//
// JAF 1/11/2026 

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include<Eigen/Core> // VectorXd
#include<Utilities/SparseDiagExpr.hpp>
#include "../../CoeffOpMixIn.hpp" 
#include "../../Vector.hpp"
#include "../../VectorXD.hpp"
#include "../../LinOpTraits.hpp" // callable_traits<F>  

namespace LinOps{

template<typename FUNC_STORAGE_T>
class AutonomousCoeff : public CoeffOpMixIn<AutonomousCoeff<FUNC_STORAGE_T>>, public LinOpBase1D<AutonomousCoeff<FUNC_STORAGE_T>>, public LinOpBaseXD<AutonomousCoeff<FUNC_STORAGE_T>>
{
  private: 
    // Member Data ----------------------------------- 
    FUNC_STORAGE_T m_function; 
    std::variant<Mesh1D_WPtr_t, MeshXD_WPtr_t> m_mesh_ptr; 
    Eigen::VectorXd m_diag_vals;
    std::size_t m_prod_after; 

  public:
    // Constructors + Destructor ==========================================================
    AutonomousCoeff()=delete; // no default constructor

    // from callable + Mesh1D
    AutonomousCoeff(FUNC_STORAGE_T f_init, const Mesh1D_SPtr_t& m)
      : m_function(f_init), m_diag_vals(0), m_prod_after(1) 
    {
      constexpr bool returns_double = std::is_same<double, typename LinOps::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value; 
      constexpr bool args_more_than_zero = (traits::callable_traits<FUNC_STORAGE_T>::num_args > 0);
      static_assert(args_more_than_zero && returns_double,"Error constructing coeff: F doesn't return double or have args > 0");  
      if(m) set_mesh(m);
    }

    // from callable + MeshXD
    AutonomousCoeff(FUNC_STORAGE_T f_init, const MeshXD_SPtr_t& m  = nullptr)
      : m_function(f_init), m_diag_vals(0), m_prod_after(1) 
    {
      constexpr bool returns_double = std::is_same<double, typename LinOps::traits::callable_traits<FUNC_STORAGE_T>::result_type>::value; 
      constexpr bool args_more_than_zero = (traits::callable_traits<FUNC_STORAGE_T>::num_args > 0);
      static_assert(args_more_than_zero && returns_double,"Error constructing coeff: F doesn't return double or have args > 0");  
      if(m) set_mesh(m);
    }

    // copy constructor
    AutonomousCoeff(const AutonomousCoeff& other)=default; 
    
    // destructor
    ~AutonomousCoeff()=default;

    // Member Funcs =============================================================
    // Matrix getters 
    auto GetMat()
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals, m_prod_after); 
    };
    auto GetMat() const
    {
      return SparseDiag<Eigen::VectorXd, SparseDiagPattern::CYCLE>(m_diag_vals, m_prod_after); 
    };

    // return weak_ptr of Mesh1D pointed to
    Mesh1D_WPtr_t get_weak_mesh1d() const { 
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr); 
      return Mesh1D_WPtr_t{}; 
    }

    // return weak_ptr of MeshXD pointed to
    MeshXD_WPtr_t get_weak_meshxd() const { 
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr); 
      return MeshXD_WPtr_t{}; 
    }

    // return Mesh1D pointed to 
    Mesh1D_SPtr_t get_mesh1d() const {
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr).lock(); 
      return Mesh1D_SPtr_t{};
    } 

    // return MeshXD pointed to 
    MeshXD_SPtr_t get_meshxd() const {
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr).lock(); 
      return MeshXD_SPtr_t{};
    } 

    // updated m_mesh_ptr, recaclulate m_diag_vals
    void set_mesh(const Mesh1D_SPtr_t& m){
      // ensure we aren't resetting the mesh again
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr))
      {
        auto my_m = std::get<Mesh1D_WPtr_t>(m_mesh_ptr); 
        if(!my_m.owner_before(m) && !m.owner_before(my_m)) return;
      }
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("AutonomousCoeff.set_mesh(m) error: std::shared_ptr<const MeshXD> is expried"); 
      m_mesh_ptr.emplace<Mesh1D_WPtr_t>(m); // store the mesh  
      // perform work on m 

      constexpr std::size_t N = LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args; 
      // Mesh1D has 1 dim
      if(N > 1) throw std::runtime_error("AutonomousCoeff.set_mesh(m) error: # of args in F > 1 (# of dims in Mesh1D)"); 

      // pack first N dimensions into small MeshXD
      m_diag_vals = LinOps::make_Discretization(m, m_function).values(); 
      m_prod_after = 1; 
    }

    // updated m_mesh_ptr, recaclulate m_diag_vals
    void set_mesh(const MeshXD_SPtr_t& m){
      // ensure we aren't resetting the mesh again
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr))
      {
        auto my_m = std::get<MeshXD_WPtr_t>(m_mesh_ptr); 
        if(!my_m.owner_before(m) && !m.owner_before(my_m)) return;
      }
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("AutonomousCoeff.set_mesh(m) error: std::shared_ptr<const MeshXD> is expried"); 
      m_mesh_ptr.emplace<MeshXD_WPtr_t>(m); // store the mesh  
      // perform work on m 

      constexpr std::size_t N = LinOps::traits::callable_traits<FUNC_STORAGE_T>::num_args; 
      if(N > m->dims()) throw std::runtime_error("AutonomousCoeff.set_mesh(m) error: # of args in F > # of dims in MeshXD"); 

      // pack first N dimensions into smaller MeshXD
      std::vector<LinOps::Mesh1D_SPtr_t> s(N);
      for(auto i=0; i<N; i++) s[i] = m->GetMesh(i);  
      auto sub_dims = std::make_shared<const LinOps::MeshXD>(std::move(s)); 
      
      // recalculate + store # of cycles 
      m_diag_vals = LinOps::make_Discretization(sub_dims, m_function).values(); 
      m_prod_after = m->sizes_product() / m_diag_vals.size(); 
    }
}; 

} // end namespace LinOps 

#endif // AutonomousCoeff.hpp 