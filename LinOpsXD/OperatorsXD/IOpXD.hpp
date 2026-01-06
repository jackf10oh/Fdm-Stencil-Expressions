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

class IOpXD: public LinOpBaseXD<IOpXD> 
{
  private:
    // private type defs
    typedef typename CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Matrix_t; 
    // member data 
    Matrix_t m_mat; 
  public:
    // constructors
    IOpXD(MeshXDPtr_t m=nullptr){set_mesh(m);}; 


    // destructors 
    ~IOpXD()=default; 
    // member functions 
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 
    DiscretizationXD apply(const DiscretizationXD& d){ return d; }; 
    void set_mesh(MeshXDPtr_t m){
      // do nothing on nullptr or same ptr 
      if(m==nullptr || m==m_mesh_ptr) return; 
      m_mesh_ptr = m; 

      std::size_t s = m_mesh_ptr->sizes_product(); 
      m_mat.resize(s,s); 
      m_mat.setIdentity(); 
    } // end set_mesh()  
}; 

#endif // IOpXD.hpp 