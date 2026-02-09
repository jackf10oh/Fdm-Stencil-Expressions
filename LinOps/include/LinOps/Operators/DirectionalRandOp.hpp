// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H

#include<Utilities/BlockDiagExpr.hpp> 
#include<Utilities/HighDimExpr.hpp> 

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
      std::size_t m_prod_before; // which Mesh1D the operator acts on. 
      std::size_t m_prod_after; // which Mesh1D the operator acts on. 
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
    auto GetMat(){ return make_HighDim(make_BlockDiag( m_mat,m_prod_before),m_prod_after); }; 
    auto GetMat() const { return make_HighDim(make_BlockDiag( m_mat,m_prod_before),m_prod_after); }; 

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
      std::size_t s = m->dim_size(m_direction); 
      m_mat = Eigen::MatrixXd::Random(s,s).sparseView(); 
      m_prod_before = m->sizes_middle_product(m_direction+1,m->dims()); 
      m_prod_after = m->sizes_middle_product(0,m_direction); 
    } // end set_mesh(MeshXD_SPtr_t)  
}; 

} // end namespace LinOps 

#endif 