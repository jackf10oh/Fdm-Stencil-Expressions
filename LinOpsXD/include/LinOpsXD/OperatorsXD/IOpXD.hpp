// IOpXD.hpp
//
// Higher dimension IOp
//
// JAF 1/3/2025 

#ifndef IDENTITYOPXD_H
#define IDENTITYOPXD_H

#include<iostream>
#include<Eigen/Sparse>
#include "../LinearOpBaseXD.hpp" 

namespace LinOps{

class IOpXD: public LinOpBaseXD<IOpXD> 
{
  private:
    // private type defs
    typedef typename CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Matrix_t; 
    // member data 
    Matrix_t m_mat; 
  public:
    // Constructors + Destructors ==========================
    IOpXD() : m_mat(0,0){}; 
    IOpXD(const MeshXD_SPtr_t& m){ set_mesh(m);}; 

    // destructor 
    ~IOpXD()=default; 
    // Fember Functions ==============================================
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 
    DiscretizationXD apply(const DiscretizationXD& d){ return d; }; 
    void set_mesh(const MeshXD_SPtr_t& m){
      // ensure we aren't resetting the mesh again
      if(!this->m_mesh_ptr.owner_before(m) && !m.owner_before(this->m_mesh_ptr)) return;

      // store the mesh  
      this->m_mesh_ptr = m;

      // perform work on m 
      std::size_t s = m->sizes_product(); 
      m_mat.resize(s,s); 
      m_mat.setIdentity(); 
    } // end set_mesh()  
}; 

} // end namespace LinOps 

#endif // IOpXD.hpp 