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
#include "../../LinOpsXD/LinearOpBaseXD.hpp"
#include "../../FDStencils/DiffOps/NthDerivOp.hpp"

namespace Fds{
using namespace LinOps; 

class DirectionalNthDerivOp : public LinOpBaseXD<DirectionalNthDerivOp>
{
  private:
    typedef typename CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Matrix_t; 
    // Member Data ----------------------------------------------
    Matrix_t m_mat; 
    std::size_t m_dir; 
    std::size_t m_order; 
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 
    NthDerivOp m_onedim_stencil;

  public:
    // Constructors + Destructor =====================================================
    DirectionalNthDerivOp(MeshXDPtr_t m, std::size_t order=1, std::size_t dir=0)
      : m_order(order), m_dir(dir), m_prod_before(std::size_t{1}), m_prod_after(std::size_t{1}), m_onedim_stencil(order) 
    {set_mesh(m);};
    DirectionalNthDerivOp(std::size_t order=1, std::size_t dir=0, MeshXDPtr_t m=MeshXDPtr_t{})
      : m_order(order), m_dir(dir), m_prod_before(std::size_t{1}), m_prod_after(std::size_t{1}), m_onedim_stencil(order) 
    {set_mesh(m);};
    // destructor 
    ~DirectionalNthDerivOp()=default; 

    // Member Functions =====================================================
    auto& GetMat(){ return m_mat; } 
    const auto& GetMat() const{return m_mat;} 
    void set_mesh(MeshXDPtr_t m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 
      // take ownership of mesh 
      if(!locked) return;
      // store the mesh   
      m_mesh_ptr = m; 

      // perform work on locked... 

      // check this->direction < mesh-> # dims 
      if(m_dir >= locked->dims()) throw std::runtime_error("DirectionalNthDerivOp.set_mesh() error: direction >= MeshXD.dims()"); 
      m_onedim_stencil.set_mesh(locked->GetMeshAt(m_dir)); 
      // std::cout << m_onedim_stencil.GetMat() << std::endl;

      m_prod_before = locked->sizes_middle_product(0,m_dir); 
      m_prod_after = (locked->dims() > m_dir) ? locked->sizes_middle_product(m_dir+1, locked->dims()) : 1; 

      // std::cout <<"dims before" <<  m_prod_before << std::endl; 
      // std::cout <<"dims after" <<  m_prod_after << std::endl; 
      
      Matrix_t I; 
      if(m_prod_before>1){
        I.resize(m_prod_before, m_prod_before); 
        I.setIdentity(); 
        m_mat = Eigen::KroneckerProductSparse(m_onedim_stencil.GetMat(), I); 
        // std::cout << "intermediary ---------------" << std::endl << m_mat << std::endl; 
      }
      else{
        m_mat = m_onedim_stencil.GetMat(); 
        // std::cout << "moved ---------------" << std::endl << m_mat << std::endl; 
      }; 
      if(m_prod_after>1){
        I.resize(m_prod_after, m_prod_after); 
        I.setIdentity(); 
        Matrix_t temp = Eigen::KroneckerProductSparse(I, m_mat);  
        m_mat = std::move(temp); 
      }; 
    }
    std::size_t Order() const {return m_order; }
    std::size_t Direction() const {return m_dir; }
}; 

} // end namespace Fds 

#endif // DirectionalNthDerivOp.hpp 

