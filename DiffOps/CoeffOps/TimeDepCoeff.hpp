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
    TimeDepCoeff(const std::function<double(double)>& f_init, MeshPtr_t m=nullptr)
      :CoeffOpBase(m)
    {
      if(!f_init) throw std::runtime_error("must assign function to AutonomousCoeff"); 
      m_function = f_init; 
      set_mesh(m);
    }
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
      set_mesh(m_mesh_ptr); 
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
      // do nothing on nullptr 
      if(m==nullptr) return; 
      
      // store m into m_mesh_ptr. checks null 
      LinOpBase::set_mesh(m);

      m_stencil.resize(m->size(), m->size()); 
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
