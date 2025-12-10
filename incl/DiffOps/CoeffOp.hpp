// CoeffOp.hpp
//
// going to need some way to model equations like 
// Ut = c(t,x) * Uxx + c(t,x) * Ux
//
// CoeffOp.hpp

#include "../FdmPlugin.hpp"
#include "../../LinOps/LinearOpBase.hpp"
#include "../BoundaryCond.hpp"

#ifndef COEFFOP_H
#define COEFFOP_H 

class CoeffOp : public LinOpBase<CoeffOp>
{
  public:
    CoeffOp(MeshPtr_t m=nullptr){set_mesh(m);};
    void SetTime(double t)
    {
      m_current_time = t; 
      auto rows = m_stencil.rows();
      auto cols = m_stencil.cols();
      m_stencil = t * Eigen::MatrixXd::Identity(rows,cols); 
    }
    Eigen::MatrixXd& GetMat(){ return m_stencil; }; 
    const Eigen::MatrixXd& GetMat() const { return m_stencil; }; 
    using LinOpBase<CoeffOp>::apply;     
    void set_mesh(MeshPtr_t m)
    {
      if(m==m_mesh_ptr) return; 
      if(m==nullptr){m_stencil.resize(0,0); return;}; 
      m_mesh_ptr = m; 
      auto rows = m->size();
      auto cols = m->size();
      m_stencil = m_current_time * Eigen::MatrixXd::Identity(rows,cols); 
    }
};

#endif