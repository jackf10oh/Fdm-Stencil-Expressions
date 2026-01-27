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

class RandLinOp : public LinOpBase<RandLinOp>
{
  private:
    Eigen::MatrixXd m_Mat; 
  public: 
    // Constructors --------------------------
    RandLinOp(std::size_t s=1)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s); 
    }
    RandLinOp(const Mesh1D_SPtr_t& m)
    {
      set_mesh(m);  
    }
    // Destructors ----------------------------
    ~RandLinOp()=default;
    // member funcs
    Eigen::MatrixXd& GetMat() { return m_Mat; };
    const Eigen::MatrixXd& GetMat() const { return m_Mat; };
    void set_mesh(const Mesh1D_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      m_mesh_ptr = m; // store the mesh  
      // fill m_Mat with random entries 
      m_Mat = Eigen::MatrixXd::Random(m->size(),m->size()); 
    };
    void resize(std::size_t s)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s);
    };
    void resize(std::size_t s1, std::size_t s2)
    {
      m_Mat = Eigen::MatrixXd::Random(s1,s2);
    }
};

} // end namespace LinOps 

#endif // RandLinOp.hpp