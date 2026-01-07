// AutonomousCoeff.hpp
//
//
//
// JAF 12/12/2025

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include<functional>
#include<Eigen/Core>
#include "CoeffOpBase.hpp"

class AutonomousCoeff : public CoeffOpBase<AutonomousCoeff>
{
  public:
    using Derived_t = AutonomousCoeff; 
  public:
    std::function<double(double)> m_function;  
    Eigen::RowVectorXd m_vals; 
  public:
    // constructors =================================================
    // no default constructor
    AutonomousCoeff()=delete; 
    // from callable + mesh
    AutonomousCoeff(const std::function<double(double)>& f_init, MeshPtr_t m=nullptr)
    {
      if(!f_init) throw std::runtime_error("must assign function to AutonomousCoeff"); 
      m_function = f_init; 
      set_mesh(m);
    }
    // from callable
    AutonomousCoeff(const AutonomousCoeff& other) 
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
    AutonomousCoeff(Func_t f){
      m_function=f; 
      set_mesh(m_mesh_ptr); 
    }
    // destructors =============================================================
    ~AutonomousCoeff()=default; 
    // member functions ========================================================
    // Do nothing on SetTime. Autonomous depends on mesh only. 
    void SetTime_impl(double t){}
    // GetScalar. doesn't exist since Autonomous can't be cast -> scalar value 
    double GetScalar() const =delete;  
    // Matrix getters. 
    auto GetMat(){ return SparseDiag(m_vals); }
    const auto GetMat() const { return SparseDiag(m_vals); }
    // apply to discretization. Override CoeffOpBase since no .GetScalar() 
    Discretization1D apply(const Discretization1D& d){
      Discretization1D result(d.mesh()); 
      result = GetMat() * d.values(); 
      return result; 
    };
    // AutonomousCoeff depends on Spacial mesh as well. 
    void set_mesh(MeshPtr_t m)
    {
      // do nothing on nullptr or same mesh
      if(m==nullptr) return; 
      // store m into m_mesh_ptr. checks null 
      m_mesh_ptr=m; 
      std::size_t s1 = m_mesh_ptr->size(); 
      m_vals.resize(s1); 
      for(std::size_t i=0; i<s1; i++){
        m_vals[i] = m_function(m_mesh_ptr->operator[](i)); 
      }
    }
    
    // operators ========================================================
    // assignment by callable 
    template<
    typename Func_t,
    typename = std::enable_if_t<
      std::is_same_v< 
        double,
        std::invoke_result_t<Func_t,double>
        >
      >
    >
    AutonomousCoeff& operator=(Func_t f){
      m_function=f; 
      return *this; 
    }
};

#endif // AutonomousCoeff.hpp 