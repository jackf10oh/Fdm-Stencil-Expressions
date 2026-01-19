// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include<Eigen/Core>
#include<Eigen/SparseCore>
#include<Eigen/Sparse>

#include "../LinearOpBase.hpp"

namespace LinOps{

#ifndef CUSTOM_IDENTITY_MATRIX_STORAGE
using CustomStorage_t = Eigen::MatrixXd;
#else
using CustomStorage_t = typename CUSTOM_IDENTITY_MATRIX_STORAGE;
#endif

class IOp : public LinOpBase<IOp>
{
  private:
    CustomStorage_t m_Mat; 
  public: 
    // Constructors --------------------------
    IOp(std::size_t s_init=0)
    {
      resize(s_init); 
    } 
    IOp(MeshPtr_t m)
    {
      set_mesh(m);
    } 
    // Destructors ----------------------------
    ~IOp()=default;
    // member funcs
    CustomStorage_t& GetMat() { return m_Mat; };
    const CustomStorage_t& GetMat() const { return m_Mat; };
    Discretization1D apply(const Discretization1D& d_arr) const 
    { 
      return d_arr;
    }; 
    void set_mesh(MeshPtr_t m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 
      if(!locked) return; 
      m_mesh_ptr = m; // store the mesh  
      // set m_Mat to I(s,s) 
      m_Mat.resize(locked->size(), locked->size()); 
      m_Mat.setIdentity(); 
    };
    void resize(std::size_t s=0)
    {
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };
    void resize(std::size_t s1,std::size_t s2) 
    {
      m_Mat.resize(s1,s2); 
      m_Mat.setIdentity(); 
    };
};

} // end namespace LinOps 

#endif