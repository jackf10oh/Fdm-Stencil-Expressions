// RandLinOp.hpp
//
// header file for a linop wrapper of a random matrix
//
// JAF 12/6/2025

#ifndef RANDLINOP_H
#define RANDLINOP_H

#include<Eigen/Core>

#include "../LinearOpBase.hpp"

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
    RandLinOp(MeshPtr_t m)
    {
      set_mesh(m);  
    }
    // Destructors ----------------------------
    ~RandLinOp()=default;
    // member funcs
    Eigen::MatrixXd& GetMat() { return m_Mat; };
    const Eigen::MatrixXd& GetMat() const { return m_Mat; };
    void set_mesh(MeshPtr_t m)
    {
      // ensure we aren't resetting the mesh again, or setting to nullptr
      auto locked = m.lock(); 
      if(!locked) return; // do nothing on nullptr. or throw an error 
      if(locked == m_mesh_ptr.lock()) return; // do nothing if m,m_mesh_ptr both point to same mesh 
      m_mesh_ptr = m; // store the mesh  
      // file matrix with random entries 
      m_Mat = Eigen::MatrixXd::Random(locked->size(),locked->size()); 
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

#endif // RandLinOp.hpp