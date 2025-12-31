// TimeDepCoeff.hpp
//
//
//
// JAF 12/26/2025 

#ifndef TIMEDEPCOEFF_H
#define TIMEDEPCOEFF_H

#include "CoeffOpBase.hpp"

class TimeDepCoeff : public CoeffOpBase<TimeDepCoeff>
{
  public:
    using Derived_t = TimeDepCoeff; 
  public:
    std::function<double(double)> m_function;  
  public:
    // constructors 
    TimeDepCoeff()=delete; // no default constructor
    // from callable + mesh 
    TimeDepCoeff(const std::function<double(double)>& f_init, MeshPtr_t m=nullptr)
    {
      if(!f_init) throw std::runtime_error("must assign function to AutonomousCoeff"); 
      m_function = f_init; 
      set_mesh(m);
    }
    // copy constructor
    TimeDepCoeff(const TimeDepCoeff& other) 
    {
      if(!other.m_function) throw std::runtime_error("must assign function to AutonomousCoeff"); 
      m_function = other.m_function; 
      set_mesh(other.m_mesh_ptr); 
    }
    // From lambdas, etc ...
    template<
    typename Func_t,
    typename = std::enable_if_t<
      std::is_same_v< 
        double,
        std::invoke_result_t<Func_t,double>
        >
      >
    >
    TimeDepCoeff(Func_t f){
      m_function=f; 
    }
    // destructors
    ~TimeDepCoeff()=default;
    // member funcs 
    void SetTime_impl(double t){
      // store the new time 
      m_current_time=t; 
      // updated current m_stencil 
      m_stencil.setIdentity(); 
      m_stencil *= m_function(t); 
    };
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };  
    Discretization1D apply(const Discretization1D& d){
      Discretization1D result(d.mesh()); 
      result = m_stencil * d.values(); 
      return result; 
    };
    void set_mesh(MeshPtr_t m)
    {
      // do nothing on nullptr or same mesh
      if(m==nullptr || m==m_mesh_ptr) return; 
      
      // store m into m_mesh_ptr. checks null 
      m_mesh_ptr = m; 

      // resize stencil and set as F(t)
      m_stencil.resize(m_mesh_ptr->size(), m_mesh_ptr->size()); 
      m_stencil.setIdentity(); 
      m_stencil *= m_function(m_current_time); 
    }
    template<
    typename Func_t,
    typename = std::enable_if_t<
      std::is_same_v< 
        double,
        std::invoke_result_t<Func_t,double>
        >
      >
    >
    TimeDepCoeff& operator=(Func_t f){
      m_function=f; 
      SetTime_impl(m_current_time);
      return *this; 
    }
 
}; 

#endif // TimeDepCoeff.hpp 
