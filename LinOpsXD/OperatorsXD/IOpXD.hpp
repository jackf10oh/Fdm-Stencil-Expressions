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
    IOpXD(MeshXDPtr_t m=MeshXDPtr_t{}){set_mesh(m);}; 

    // destructor 
    ~IOpXD()=default; 
    // Fember Functions ==============================================
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 
    DiscretizationXD apply(const DiscretizationXD& d){ return d; }; 
    void set_mesh(MeshXDPtr_t m){
      // ensure we aren't resetting the mesh again
      if(!this->m_mesh_ptr.owner_before(m) && !m.owner_before(this->m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 
      if(!locked) return; 
      this->m_mesh_ptr = m; // store the mesh  
      // perform work on locked 
      std::size_t s = locked->sizes_product(); 
      m_mat.resize(s,s); 
      m_mat.setIdentity(); 
    } // end set_mesh()  
}; 

} // end namespace LinOps 

#endif // IOpXD.hpp 