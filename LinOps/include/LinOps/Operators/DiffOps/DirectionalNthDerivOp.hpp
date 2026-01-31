// DirectionalNthDerivOp.hpp
//
// X Dimensional version of NthDerivOp.hpp
//
// JAF 1/10/2026 

#ifndef DIRECTIONALNTHDERIVOP_H
#define DIRECTIONALNTHDERIVOP_H 

#include<cstdint>
#include<Eigen/Sparse> 
#include<unsupported/Eigen/KroneckerProduct> 
#include "../../LinearOpBase.hpp" 
#include "NthDerivOp.hpp" 

namespace LinOps{

class DirectionalNthDerivOp : public LinOpMixIn<DirectionalNthDerivOp>, public LinOpBaseXD<DirectionalNthDerivOp>
{
  private:
    // Member Data ----------------------------------------------
    MeshXD_WPtr_t m_mesh_ptr; 
    MatrixStorage_t m_mat; 
    std::size_t m_dir; 
    std::size_t m_order; 
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 
    NthDerivOp m_onedim_stencil;

  public:
    // Constructors + Destructor =====================================================
    // meshxd + order + direction 
    DirectionalNthDerivOp(const MeshXD_SPtr_t& m, std::size_t order=1, std::size_t dir=0)
      : m_order(order), m_dir(dir), m_prod_before(std::size_t{1}), m_prod_after(std::size_t{1}), m_onedim_stencil(order) 
    {set_mesh(m);};

    // order + direction 
    DirectionalNthDerivOp(std::size_t order=1, std::size_t dir=0)
      : m_order(order), m_dir(dir), m_prod_before(1), m_prod_after(1), m_onedim_stencil(order) 
    {};

    // destructor 
    ~DirectionalNthDerivOp()=default; 

    // Member Functions =====================================================
    // getters to Direction + Order 
    std::size_t Direction() const {return m_dir; }
    std::size_t Order() const {return m_order; }

    // Getters to Matrix 
    auto& GetMat(){ return m_mat; } 
    const auto& GetMat() const{return m_mat;} 

    // getters to mesh 
    MeshXD_WPtr_t get_weak_meshxd() const { return m_mesh_ptr; }
    MeshXD_SPtr_t get_meshxd() const { return m_mesh_ptr.lock(); }

    // set DIrectional Derivative to operate on a domain mesh 
    void set_mesh(const MeshXD_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;

      // store the mesh   
      m_mesh_ptr = m; 

      // perform work on m ... 

      // check this->direction < mesh-> # dims 
      if(m_dir >= m->dims())
      {
        // throw an error ...
        // throw std::runtime_error("DirectionalNthDerivOp.set_mesh() error: direction >= MeshXD.dims()"); 
        // or just resize to identity 
        std::size_t s = m->sizes_product(); 
        m_mat.resize(s, s); 
        m_mat.setIdentity(); 
        return; 
      }
      m_onedim_stencil.set_mesh(m->GetMeshAt(m_dir)); 
      // std::cout << m_onedim_stencil.GetMat() << std::endl;

      m_prod_before = m->sizes_middle_product(0,m_dir); 
      m_prod_after = (m->dims() > m_dir) ? m->sizes_middle_product(m_dir+1, m->dims()) : 1; 

      MatrixStorage_t I; 
      if(m_prod_before>1){
        I.resize(m_prod_before, m_prod_before); 
        I.setIdentity(); 
        m_mat = Eigen::KroneckerProductSparse(m_onedim_stencil.GetMat(), I); 
      }
      else{
        m_mat = m_onedim_stencil.GetMat(); 
      }; 
      if(m_prod_after>1){
        I.resize(m_prod_after, m_prod_after); 
        I.setIdentity(); 
        MatrixStorage_t temp = Eigen::KroneckerProductSparse(I, m_mat);  
        m_mat = std::move(temp); 
      }; 
    }
    
}; // end class DirectionalNthDerivOp

} // end namespace LinOps 

#endif // DirectionalNthDerivOp.hpp 

