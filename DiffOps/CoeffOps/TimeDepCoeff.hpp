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
    double m_val; 
  public:
    // constructors ==========================================================
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
    // destructors =============================================================
    ~TimeDepCoeff()=default;
    // member funcs =============================================================
    // update state of TimeDepCoeff from a given t 
    void SetTime_impl(double t){
      // SetTime() stores the new time... 
      // Store the result of m_function(t) 
      m_val = m_function(t); 
    };
    // return current stored result 
    double GetScalar() const{
      return m_val; 
    }
    // operators ==================================================================
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
      SetTime_impl(m_current_time); // possibly remove? 
      return *this; 
    }
}; 

#endif // TimeDepCoeff.hpp 
