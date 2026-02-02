// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H

#include<iostream>
#include<Eigen/Sparse> 
#include<unsupported/Eigen/KroneckerProduct> 

#include "../LinearOpBase.hpp" 

namespace LinOps{

class DirectionalRandOp: public LinOpMixIn<DirectionalRandOp>, public LinOpBaseXD<DirectionalRandOp> 
{
  private:
    private:
      // Member Data --------------------------- 
      MeshXD_WPtr_t m_mesh_ptr; 
      MatrixStorage_t m_mat; 
      std::size_t m_direction; // which Mesh1D the operator acts on. 
  public:
    // Constructors + Destructo =====================================
    // direction only 
    DirectionalRandOp(std::size_t dir_init=0) 
      : m_direction(dir_init)
    {};
    // mesh + direction 
    DirectionalRandOp(const MeshXD_SPtr_t& m, std::size_t dir_init=0) 
      : m_direction(dir_init) 
    {set_mesh(m);}; 
    // destructor
    ~DirectionalRandOp()=default; 

    // Member Funcs ============================================== 
    
    // matrix getters
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 

    // return weak_ptr of MeshXD pointed to
    MeshXD_WPtr_t get_weak_meshxd() const { return m_mesh_ptr; }

    // return MeshXD pointed to 
    MeshXD_SPtr_t get_meshxd() const { return m_mesh_ptr.lock(); } 

    // set operator to domain mesh 
    void set_mesh(const MeshXD_SPtr_t& m){
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("IOp.set_mesh(m) error: std::shared_ptr<const Mesh1D> is expried"); 
      m_mesh_ptr = m; // store the mesh  

      // perform work on m ...
      
      // check this->direction < mesh-> # dims 
      if(m_direction >= m->dims())
      {
        // dont allow direction >= dims ...
        // throw std::runtime_error("DirectionalNthDerivOp.set_mesh() error: direction >= MeshXD.dims()");
        // just set m_mat to identity ...  
        std::size_t s = m->sizes_product(); 
        m_mat.resize(s,s); 
        m_mat.setIdentity();  
        return; 
      }

      // get a random matrix R that matches specific direction 
      std::size_t s = m->dim_size(m_direction); 
      MatrixStorage_t R = Eigen::MatrixXd::Random(s,s).sparseView();      
      // identity matrix for calculations later 
      MatrixStorage_t I;  

      std::size_t prod_before = m->sizes_middle_product(0,m_direction); 
      std::size_t prod_after = (m->dims() > m_direction) ? m->sizes_middle_product(m_direction+1, m->dims()) : 1; 

      if(prod_before>1){
        I.resize(prod_before, prod_before); 
        I.setIdentity(); 
        m_mat = Eigen::KroneckerProductSparse(R, I); 
      }
      else{
        m_mat = std::move(R); 
      }; 
      if(prod_after>1){
        I.resize(prod_after, prod_after); 
        I.setIdentity(); 
        MatrixStorage_t temp = Eigen::KroneckerProductSparse(I, m_mat);  
        m_mat = std::move(temp); 
      }; 
    } // end set_mesh(MeshXD_SPtr_t)  
}; 

} // end namespace LinOps 

#endif 