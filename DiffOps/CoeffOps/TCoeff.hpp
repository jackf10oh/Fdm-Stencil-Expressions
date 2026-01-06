// TCoeff.hpp
//
//
//
// JAF 12/10/2025

#ifndef CURRENTTIMEOP_H
#define CURRENTTIMEOP_H

#include "CoeffOpBase.hpp"

struct read_func_obj; 

class TCoeff : public CoeffOpBase<TCoeff>
{
  friend read_func_obj; 

  using Derived_t = TCoeff; 
  public:
    // constructors 
    TCoeff(MeshPtr_t m=nullptr)
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
      // do nothing on nullptr 
      if(m==nullptr || m==m_mesh_ptr) return; 
      // store m into m_mesh_ptr.
      m_mesh_ptr=m; 
      // resize current stencil 
      m_stencil.resize(m->size(), m->size()); 
      m_stencil.setIdentity(); 
      m_stencil *= m_current_time; 
    }
}; 

struct read_func_obj{
  template<typename RHS>
  void operator()(const TCoeff& c, const RHS& rhs){std::cout << "read_func_obj called on t="<<c.m_current_time << "rhs: " << rhs << std::endl; }; 
};

#endif // TCoeff.hpp
