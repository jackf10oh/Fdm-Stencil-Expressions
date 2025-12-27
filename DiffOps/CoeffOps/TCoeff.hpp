// TCoeff.hpp
//
//
//
// JAF 12/10/2025

#ifndef CURRENTTIMEOP_H
#define CURRENTTIMEOP_H

#include "../FdmPlugin.hpp"
#include "CoeffOpBase.hpp"

class TCoeff : public CoeffOpBase<TCoeff>
{
  using Derived_t = TCoeff; 
  public:
    // constructors 
    TCoeff(MeshPtr_t m=nullptr)
      :CoeffOpBase(m)
    {
      set_mesh(m);
    }
    // member functions  
    void SetTime_impl(double t)
    {
      if(m_mesh_ptr) m_stencil.resize(m_mesh_ptr->size(), m_mesh_ptr->size());
      m_stencil.setIdentity(); 
      m_stencil = m_current_time * m_stencil;  
    };
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };  
    Discretization1D apply(const Discretization1D& d){
      Discretization1D result(d.mesh()); 
      result = m_stencil * d.values(); 
      return result; 
    }
    void set_mesh(MeshPtr_t m)
    {
      // store m into m_mesh_ptr. checks null 
      LinOpBase::set_mesh(m);

      // do nothing on nullptr 
      if(m==nullptr) return; 
      m_stencil.resize(m->size(), m->size()); 
      m_stencil.setIdentity(); 
      m_stencil *= m_current_time; 
    }
};

#endif // TCoeff.hpp
