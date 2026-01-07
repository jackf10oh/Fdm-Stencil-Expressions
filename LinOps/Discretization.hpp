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

#include "Mesh.hpp"
#include "LinOpTraits.hpp"

class Discretization1D
{
  private:
    // member data 
    MeshPtr_t m_mesh_ptr; // shared pointer to mesh the function "maps" into the reals
    Eigen::VectorXd m_vals; // eigen array of functions vals at mesh points
  public:
    // constructors -------------------------------------------
    // from size 
    Discretization1D(std::size_t size_init=0): m_mesh_ptr(), m_vals(size_init){};
    // from mesh
    // Discretization1D(MeshPtr_t mesh_init): m_mesh_ptr(mesh_init), m_vals(mesh_init ? mesh_init->size() : 0){};
    // from const mesh ptr
    Discretization1D(const MeshPtr_t& mesh_init) : m_mesh_ptr(mesh_init), m_vals(mesh_init ? mesh_init->size() : 0){}
    // Copy 
    Discretization1D(const Discretization1D& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {};
    // Move 
    Discretization1D(Discretization1D&& other)
      : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {other.m_mesh_ptr=nullptr;};

    // destructors ------------------------------------------
    virtual ~Discretization1D()=default; 

    // member functions 
    // return the size of the underlying vector
    std::size_t size() const { return m_vals.size(); }

    // get underlying values
    Eigen::VectorXd& values() { return m_vals; }
    const Eigen::VectorXd& values() const { return m_vals; }
    
    // indexing
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ assert(i<size()); return m_vals[i]; }
    const double& at(std::size_t i) const { assert(i<size()); return m_vals[i]; }

    // get underlying mesh
    auto mesh(){return m_mesh_ptr; }; 
    const auto& mesh() const{return m_mesh_ptr; }

    // set vector to same size as mesh
    void match_mesh(MeshPtr_t m) { m_mesh_ptr=m; m_vals.resize(m->size()); }
    // set vector to a constant
    void set_init(double val){ m_vals.setConstant(val);}
    // set vector to match a mesh size and set it constant 
    void match_mesh(MeshPtr_t m, double val){ match_mesh(m), m_vals.setConstant(val); }

    // set vector to same size as mesh and set as value of callable
    template<
    typename F,
    typename = std::enable_if_t<std::conjunction_v<
        std::is_same<typename callable_traits<F>::result_type, double>, 
        std::bool_constant<callable_traits<F>::num_args == 1>
        >
      >
    >
    void set_init(F func)
    {
      if(m_mesh_ptr==nullptr) return; 
      std::transform(m_mesh_ptr->cbegin(), m_mesh_ptr->cend(), 
                    m_vals.begin(), 
                    [&func](const double x){return func(x);}
      ); 
    }
    // set vector to same size as mesh and set as value of callable
    template<
    typename F,
    typename = std::enable_if_t<std::conjunction_v<
        std::is_same<typename callable_traits<F>::result_type, double>, 
        std::bool_constant<callable_traits<F>::num_args == 1>
        >
      >
    >
    void set_init(MeshPtr_t m, F func)
    {
      match_mesh(m); 
      std::transform(m->cbegin(), m->cend(), 
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
      other.m_mesh_ptr=nullptr; 
      m_vals=std::move(other.m_vals); 
      return *this;
    }; 
    Discretization1D& operator=(const Eigen::VectorXd& other){ 
      m_vals=other;
      // if(m_mesh_ptr) m_vals.resize(m_mesh_ptr->size());  
      return *this;
    }
    Discretization1D& operator=(Eigen::VectorXd&& other){ 
      m_vals=std::move(other);       
      // if(m_mesh_ptr) m_vals.resize(m_mesh_ptr->size());  
      return *this;
    }
}; // end Discretization1D 

#endif // Discretization.hpp