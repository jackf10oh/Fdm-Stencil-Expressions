// FdmPlugin.hpp
//
// Inside LinOpBase. add a new member data for storing time. 
// Operators will be able to change at different time steps. 
// default implementation will do nothing. 
//
// JAF 

#ifndef FDMPLUGIN_H
#define FDMPLUGIN_H 

#define LINOP_PLUGIN Fds::FdmPlugin
#define CUSTOM_IDENTITY_MATRIX_STORAGE Eigen::SparseMatrix<double,Eigen::RowMajor>

#include<cstdint>
#include<iostream>
#include<Eigen/SparseCore>
#include "../LinOps/LinOpTraits.hpp"

namespace Fds{
using namespace LinOps; 

// Plugin Class Def ================================================= 
template<typename BaseDerived>
class FdmPlugin
{
  // Member Data --------------------------------------------------------------
  protected:
    double m_current_time; // current time of linear operator. 

  public:
    // Constructors / Destructor =================================================
    FdmPlugin() : m_current_time(0.0){};
    FdmPlugin(const FdmPlugin& other)=default; 
    ~FdmPlugin()=default; 

    // Member Funcs =================================================
    // set the current time of the Linear Operator. Do nothing by default 
    void SetTime_impl(double t){}; 
    void SetTime(double t)
    {
      // if we are an expression 
      if constexpr (internal::is_expr_crtp<typename BaseDerived::Derived_t>::value)
      {
        auto& expr = static_cast<typename BaseDerived::Derived_t&>(*this);
        // if LHS of expr is LinOp / FdmPugin 
        if constexpr(internal::is_linop_crtp<typename BaseDerived::Derived_t::LStorage_t>::value) 
        {
          // LHS sets time 
          expr.Lhs().SetTime(t);
        }
        // if RHS of expr is LinOp / FdmPugin
        if constexpr(internal::is_linop_crtp<typename BaseDerived::Derived_t::RStorage_t>::value)
        {
          // RHS sets time 
          expr.Rhs().SetTime(t);
        } 
      }
      // we aren't an expression call the actual implementor 
      else  
      {
        static_cast<typename BaseDerived::Derived_t*>(this)->SetTime_impl(t);
      }
      // store new time.   
      m_current_time = t; 
    }
    // get the time this operator is set to
    double Time() const {
      // if we are an expression 
      if constexpr (internal::is_expr_crtp<typename BaseDerived::Derived_t>::value)
      {
        const auto& expr = static_cast<const typename BaseDerived::Derived_t&>(*this);
        // if LHS of expr is LinOp / FdmPugin 
        if constexpr(internal::is_linop_crtp<typename BaseDerived::Derived_t::LStorage_t>::value) 
        {
          // LHS returns time 
          return expr.Lhs().Time();
        }
        // if RHS of expr is LinOp / FdmPugin
        else if constexpr(internal::is_linop_crtp<typename BaseDerived::Derived_t::RStorage_t>::value)
        {
          // RHS returns time 
          return expr.Rhs().Time();
        } 
      }
      // if we aren't an expression -> linops return m_current_time 
      else
      {
        return m_current_time; 
      }
    };  
};

} // end namespace Fds 

#endif // FdmPlugin.hpp