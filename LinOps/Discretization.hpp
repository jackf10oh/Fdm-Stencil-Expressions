// Discretization.hpp
//
// 1D array to store values of a function at discrete points
// inside a rectangular domain 
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

class Discretization1D
{
  private:
    // member data 
    MeshPtr_t m_mesh_ptr; // shared pointer to mesh the function "maps" into the reals
    Eigen::VectorXd m_vals; // eigen array of functions vals at mesh points
  public:
    // constructors ==============================================================
    // from size 
    Discretization1D(std::size_t size_init=0): m_mesh_ptr(), m_vals(size_init){};
    // from mesh
    // Discretization1D(MeshPtr_t mesh_init): m_mesh_ptr(mesh_init), m_vals(mesh_init ? mesh_init->size() : 0){};
    // from const mesh ptr
    Discretization1D(const std::shared_ptr<const Mesh1D>& mesh_init) : m_mesh_ptr(mesh_init), m_vals(mesh_init ? mesh_init->size() : 0){}
    // Copy -----------------------------
    Discretization1D(const Discretization1D& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {};
    // copy from Eigen::VectorXd
    Discretization1D(const Eigen::VectorXd& other)
      : m_mesh_ptr(), m_vals(other)
    {}; 
    // Move ---------------------------
    Discretization1D(Discretization1D&& other)
      : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {
      other.m_mesh_ptr=MeshPtr_t{};
    }; 
    // move from Eigen::VectorXD
    Discretization1D(Eigen::VectorXd&& other)
      : m_mesh_ptr(), m_vals(other)
    {}; 
    // destructors ------------------------------------------
    ~Discretization1D()=default; 

    // member functions 
    // return the size of the underlying vector
    std::size_t size() const { return m_vals.size(); }

    // get underlying values
    Eigen::VectorXd& values(){ return m_vals; }
    const Eigen::VectorXd& values() const { return m_vals; }
    
    // indexing
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ assert(i<size()); return m_vals[i]; }
    const double& at(std::size_t i) const { assert(i<size()); return m_vals[i]; }

    // get underlying mesh
    // const auto mesh(){return m_mesh_ptr.lock(); }; 
    MeshPtr_t mesh() const{return m_mesh_ptr; }

    // set vector to same size as mesh
    void set_mesh(MeshPtr_t m){ 
      m_mesh_ptr=m;
    }
    void resize(MeshPtr_t m){
      m_mesh_ptr=m;
      auto locked = m.lock(); 
      if(!locked) m_vals.resize(0);
      else m_vals.resize(locked->size()); 
    }
    // set_init() NON TEMPLATE ----------------------------------------------------------
    // set vector to a constant
    void set_init(double val){ m_vals.setConstant(val);}
    // set vector to match a mesh size and set it constant 
    void set_init(MeshPtr_t m, double val){ resize(m), m_vals.setConstant(val); }

    // set_init() TEMPLATE ----------------------------------------------------------------- 
    // set vector to same size as mesh and set as value of callable
    template<
    typename F,
    typename = std::enable_if_t<
      std::is_same<typename callable_traits<F>::result_type, double>::value && 
      std::bool_constant<callable_traits<F>::num_args == 1>::value
      >
    >
    void set_init(F func)
    {
      if(m_mesh_ptr.expired()) return; 
      const std::shared_ptr<const Mesh1D> m = m_mesh_ptr.lock(); 
      std::transform(m->cbegin(), m->cend(), 
                    m_vals.begin(), 
                    [&func](const double x){return func(x);}
      ); 
    }
    // set vector to same size as mesh and set as value of callable
    template<
    typename F,
    typename = std::enable_if_t<
      std::negation<std::is_arithmetic<F>>::value // placed first to short circuit? why doesn't enable_if_t< 1 && 2 && 3 > work here?  
    >, 
    typename = std::enable_if_t<
      std::is_same<typename callable_traits<F>::result_type, double>::value &&
      std::bool_constant<callable_traits<F>::num_args == 1>::value
      >
    >
    void set_init(MeshPtr_t m, F func)
    {
      const std::shared_ptr<const Mesh1D> locked_ptr = m.lock(); 
      if(m.expired()){
        m_vals.resize(0);
        return; 
      } 
      const auto locked = m.lock(); 
      m_mesh_ptr = locked; 
      m_vals.resize(locked->size()); 
      std::transform(locked->cbegin(), locked->cend(), 
                    m_vals.begin(), 
                    [&func](const double x){return func(x);}
      ); 
    }
    
    // forward iterators 
    auto begin(){ return m_vals.begin(); }
    auto end(){ return m_vals.end(); }
    auto begin() const { return m_vals.begin(); }
    auto end() const { return m_vals.end(); }
    // forward citerators
    auto cbegin() const { return m_vals.cbegin(); }
    auto cend() const { return m_vals.cend(); }
    // no reverse citerators in eigen 

    // Operators ----------------------------------------------------
    Discretization1D& operator=(const Discretization1D& other) = default;
    Discretization1D& operator=(Discretization1D&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      other.m_mesh_ptr=MeshPtr_t{}; 
      m_vals=std::move(other.m_vals); 
      return *this;
    }; 
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

#endif // Discretization.hpp