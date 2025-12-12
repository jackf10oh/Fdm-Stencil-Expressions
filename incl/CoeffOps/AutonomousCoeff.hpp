// AutonomousCoeff.hpp
//
//
//
// JAF 12/12/2025

#ifndef AUTONOMOUSCOEFF_H
#define AUTONOMOUSCOEFF_H

#include<functional>
#include "../FdmPlugin.hpp"
#include "CoeffOpBase.hpp"

class AutonomousCoeff : public CoeffOpBase<AutonomousCoeff>
{
  public:
    using Derived_t = AutonomousCoeff; 
  public:
    std::function<double(double)> m_function;  
  public:
    // constructors 
    AutonomousCoeff()=delete; // no default constructor
    AutonomousCoeff(const std::function<double(double)>& f_init, MeshPtr_t m=nullptr)
      :CoeffOpBase(m)
    {
      if(!f_init) throw std::runtime_error("must assign function to AutonomousCoeff"); 
      m_function = f_init; 
      set_mesh(m);
    }
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
    // destructors
    ~AutonomousCoeff()=default; 
    // member functions  
    void SetTime_impl(double t){};
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };  
    Discretization1D apply(const Discretization1D& d){
      Discretization1D result(d.mesh()); 
      result = m_stencil * d.values(); 
      return result; 
    };
    void set_mesh(MeshPtr_t m)
    {
      // store m into m_mesh_ptr. checks null 
      LinOpBase::set_mesh(m);

      // do nothing on nullptr 
      if(m==nullptr) return; 
      m_stencil.resize(m->size(), m->size()); 
      using T = Eigen::Triplet<double>; 
      std::vector<T> triplet_list(m->size());  
      for(int i=0; i<m->size(); i++){
        triplet_list[i] = T(i,i,m_function((*m)[i])); 
      }; 
      m_stencil.setFromTriplets(triplet_list.begin(), triplet_list.end()); 
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
    AutonomousCoeff& operator=(Func_t f){
      m_function=f; 
      set_mesh(m_mesh_ptr); 
      return *this; 
    }
};

#endif // AutonomousCoeff.hpp 