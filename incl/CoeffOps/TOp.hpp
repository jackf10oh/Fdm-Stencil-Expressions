// TOp.hpp
//
//
//
// JAF 12/10/2025

#ifndef CURRENTTIMEOP_H
#define CURRENTTIMEOP_H

#include "../FdmPlugin.hpp"
#include "CoeffOpBase.hpp"

class TOp : public CoeffOpBase<TOp>
{
  public:
    // constructors 
    TOp(MeshPtr_t m=nullptr)
      :CoeffOpBase()
    {
      set_mesh(m);
    }
    // member functions  
    void SetTime_impl(double t)
    {
      m_stencil.setIdentity(); 
      m_stencil = m_current_time * m_stencil;  
    };
    Eigen::MatrixXd& GetMat(){ return m_stencil; }; 
    const Eigen::MatrixXd& GetMat() const { return m_stencil; };  
    void set_mesh(MeshPtr_t m)
    {
      // store m into m_mesh_ptr. checks null 
      LinOpBase::set_mesh(m);

      // do nothing on nullptr 
      if(m==nullptr) return; 
      m_stencil.resize(m->size(), m->size()); 
      m_stencil.setIdentity(); 
      m_stencil = m_current_time * m_stencil; 
    }
};

#endif
