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
      : CoeffOpBase<TCoeff>(m)
    {}; 
    // member functions
    double GetScalar() const
    {
      return m_current_time; 
    }
    void SetTime_impl(double t){}; // m_current_time = t already handled by SetTime(); 
}; 

#endif // TCoeff.hpp
