// Discretization.hpp
//
// Discretization of a function F(x0), ..., F(xN) on a mesh x0, ..., xN
//
// JAF 12/5/2025

#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include<memory>
#include<cstdint>
#include<vector>
#include<Eigen/Core>
#include <type_traits>

#include "Mesh.hpp"
#include "LinOpTraits.hpp"

namespace LinOps{

class Discretization1D
{
  private:
    // Member Data --------------------------------------------------------------------
    Mesh1D_WPtr_t m_mesh_ptr; // weak pointer to mesh the function "maps" into the reals
    Eigen::VectorXd m_vals; // eigen array of functions vals at mesh points
  public:
    // Constructors / Destructors ==============================================================
    // Args -------------------------------
    // from size 
    Discretization1D(std::size_t size_init=0): m_mesh_ptr(), m_vals(size_init){};
    // from mesh
    // from const shared_ptr<Mesh1D>
    Discretization1D(const Mesh1D_SPtr_t& mesh_init) : m_mesh_ptr(mesh_init){resize(mesh_init);}
    
    // Copy -----------------------------
    Discretization1D(const Discretization1D& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {};
    
    // copy from Eigen::VectorXd
    Discretization1D(const Eigen::VectorXd& other, Mesh1D_WPtr_t mesh_init = Mesh1D_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(other)
    {}; 
    
    // Move ---------------------------
    Discretization1D(Discretization1D&& other)
      : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {}; 
    
    // move from Eigen::VectorXD
    Discretization1D(Eigen::VectorXd&& other, Mesh1D_WPtr_t mesh_init = Mesh1D_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(std::move(other))
    {}; 
    // destructors ----------------------------------
    ~Discretization1D()=default; 

    // Member Functions ==============================================================
    // return the size of the underlying vector ------------------------
    std::size_t size() const { return m_vals.size(); }

    // get underlying values ------------------
    Eigen::VectorXd& values(){ return m_vals; }
    const Eigen::VectorXd& values() const { return m_vals; }
    
    // indexing ------------------------------------
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ assert(i<size()); return m_vals[i]; }
    const double& at(std::size_t i) const { assert(i<size()); return m_vals[i]; }

    // get underlying mesh ------------------------------
    // const auto mesh(){return m_mesh_ptr.lock(); }; 
    Mesh1D_SPtr_t get_mesh1d() const{return m_mesh_ptr.lock(); }

    // set the mesh the discretization is on -----------
    auto& set_mesh(Mesh1D_WPtr_t m){ 
      m_mesh_ptr=m;
      return *this; 
    }
    // set vector to same size as mesh -----------------
    Discretization1D& resize(const Mesh1D_SPtr_t& m){
      m_mesh_ptr=m;
      m_vals.conservativeResize(m->size()); 
      return *this; 
    }

    // set vector to a constant -------------------------------------------------------------
    auto& set_init(double val){ m_vals.setConstant(val); return *this; }
    
    // resize + set constant 
    auto& set_init(const Mesh1D_SPtr_t& m, double val){ resize(m); set_init(val); return *this; }

    // set_init() TEMPLATE ----------------------------------------------------------------- 
    // set vector F(x0),...,F(xN) for x0,...,xN in stored mesh 
    template<
    typename F,
    typename = std::enable_if_t<
      std::is_same<typename traits::callable_traits<F>::result_type, double>::value && 
      std::bool_constant<traits::callable_traits<F>::num_args == 1>::value
      >
    >
    auto& set_init(F func)
    {
      auto locked = m_mesh_ptr.lock(); 
      set_init(locked, std::move(func));
      return *this; 
    }
    // resize + set vector F(x0),...,F(xN) for x0,...,xN in stored mesh 
    template<
    typename F,
    typename = std::enable_if_t<
      std::negation<std::is_arithmetic<F>>::value // placed first to short circuit? why doesn't enable_if_t< 1 && 2 && 3 > work here?  
    >, 
    typename = std::enable_if_t<
      std::is_same<typename traits::callable_traits<F>::result_type, double>::value &&
      std::bool_constant<traits::callable_traits<F>::num_args == 1>::value
      >
    >
    auto& set_init(const Mesh1D_SPtr_t& m, F func)
    {
      m_mesh_ptr = m; 

      std::size_t s1 = m->size();
      auto data = m->cbegin(); 

      m_vals.resize(s1); 
      for(auto i=0; i<s1; i++) m_vals[i] = func(data[i]); 
      return *this; 
    }
    
    // forward iterators -------------------------------------------------------
    auto begin(){ return m_vals.begin(); }
    auto end(){ return m_vals.end(); }

    // forward citerators ----------------------------------
    auto cbegin() const { return m_vals.cbegin(); }
    auto cend() const { return m_vals.cend(); }

    // no reverse citerators in eigen 

    // Operators ================================================================
    // Assignment from Discretization1D --------------------------------
    Discretization1D& operator=(const Discretization1D& other) = default;
    Discretization1D& operator=(Discretization1D&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      m_vals=std::move(other.m_vals); 
      return *this;
    }; 
    // Assignment from Eigen::VectorXd --------------------------------
    Discretization1D& operator=(const Eigen::VectorXd& other){ 
      m_vals=other;
      auto locked = m_mesh_ptr.lock(); 
      if(locked) m_vals.resize(locked ->size());  
      return *this;
    }
    Discretization1D& operator=(Eigen::VectorXd&& other){ 
      m_vals = std::move(other);
      auto locked = m_mesh_ptr.lock(); 
      if(locked) m_vals.resize(locked ->size());  
      return *this;
    }
}; // end Discretization1D 

} // end namespace LinOps 

#endif // Discretization.hpp