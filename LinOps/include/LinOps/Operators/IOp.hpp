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

class IOp : public LinOpBase1D<IOp>
{
  private:
    Mesh1D_WPtr_t m_mesh_ptr; 
    MatrixStorage_t m_Mat; 
  public: 
    // Constructors --------------------------
    IOp(std::size_t s_init=0)
    {
      resize(s_init); 
    } 
    IOp(const Mesh1D_SPtr_t& m)
    {
      set_mesh(m);
    } 
    // Destructors ----------------------------
    ~IOp()=default;

    // Member Funcs  ======================================================

    // matrix getters 
    MatrixStorage_t& GetMat() { return m_Mat; };
    const MatrixStorage_t& GetMat() const { return m_Mat; };

    // Identity just returns inputs as outputs
    Discretization1D apply(const Discretization1D& d_arr) const { return d_arr; } 

    // return weak_ptr of Mesh1D pointed to. default: just convert the shared_ptr 
    Mesh1D_WPtr_t get_weak_mesh1d() const{ return m_mesh_ptr; }

    // return Mesh1D pointed to 
    Mesh1D_SPtr_t get_mesh1d() const { return m_mesh_ptr.lock(); } 

    // fit operator to a domain mesh 
    void set_mesh(const Mesh1D_SPtr_t& m) 
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("IOp.set_mesh(m) error: std::shared_ptr<const Mesh1D> is expried"); 
      m_mesh_ptr = m; // store the mesh  
      // perform work on m 
      // set m_Mat to I(s,s) 
      std::size_t s = m->size(); 
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };
    
    // resize Identity to s; 
    void resize(std::size_t s=0)
    {
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };

}; // End IOp 

} // end namespace LinOps 

#endif