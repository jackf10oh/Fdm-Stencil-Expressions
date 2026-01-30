// RandLinOp.hpp
//
// header file for a linop wrapper of a random matrix
//
// JAF 12/6/2025

#ifndef RANDLINOP_H
#define RANDLINOP_H

#include<Eigen/Core>

#include "../LinearOpBase.hpp"

namespace LinOps{

class RandLinOp : public LinOpBase1D<RandLinOp>
{
  private:
    Eigen::MatrixXd m_Mat; 
    Mesh1D_WPtr_t m_mesh_ptr; 
  public: 
    // Constructors + Destructor ===========================
    RandLinOp(std::size_t s=0)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s); 
    }
    RandLinOp(const Mesh1D_SPtr_t& m)
    {
      set_mesh(m);  
    }
    // destructors
    ~RandLinOp()=default;

    // Member Funcs ==========================================

    // matrix getters 
    Eigen::MatrixXd& GetMat() { return m_Mat; };
    const Eigen::MatrixXd& GetMat() const { return m_Mat; };

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
      // fill m_Mat with random entries 
      std::size_t s = m->size(); 
      m_Mat = Eigen::MatrixXd::Random(s,s); 
    };

    void resize(std::size_t s)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s);
    };

}; // end  RandLinOp

} // end namespace LinOps 

#endif // RandLinOp.hpp