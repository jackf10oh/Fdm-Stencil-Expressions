// Mesh.hpp
//
// 1D Mesh classes for rectangular domains
//
// JAF 12/5/2025

#ifndef MESH_H
#define MESH_H

#include<cstdint>
#include<vector>
#include<memory> 
#include<type_traits>

// forward declaration -> aliases
class Mesh1D; 
using MeshPtr_t = std::shared_ptr<Mesh1D>;

class Mesh1D
{
  protected:
    // member data ------------------------------------------------------------------------
    std::vector<double> m_vals;
  public:
    // constructors ------------------------------------------------------------------------
    // uniformly on interval [x1,x2] 
    Mesh1D(double x1=0.0, double x2=1.0, std::size_t n_steps=11)
      : m_vals(n_steps)
    {
      if(n_steps < 2) throw std::invalid_argument("must have # of steps >= 2 in mesh");
      if(x1 >= x2) throw std::invalid_argument("must have x1 < x2 in mesh");
      double dx = (x2 - x1) / (n_steps-1); 
      for(std::size_t i=0; i<n_steps-1; i++) m_vals[i] = x1 + i*dx; 
      m_vals[n_steps-1] = x2; 
    }
    // copy constructors
    Mesh1D(const Mesh1D& other): m_vals(other.m_vals){}; 

    // destructors  ---------------------------------------------------------------------------
    virtual ~Mesh1D()=default; 

    // member functions ------------------------------------------------------------------------
    // size of mesh
    std::size_t size() const { return m_vals.size(); }
    // indexing
    double& operator[](std::size_t i){ return m_vals[i]; }
    const double& operator[](std::size_t i) const { return m_vals[i]; } 
    double& at(std::size_t i){ return m_vals.at(i); }
    const double& at(std::size_t i) const { return m_vals.at(i); }
    // forward iterators 
    auto begin(){ return m_vals.begin(); }
    auto end(){ return m_vals.end(); }
    auto begin() const { return m_vals.begin(); }
    auto end() const { return m_vals.end(); }
    // forward citerators
    auto cbegin() const { return m_vals.cbegin(); }
    auto cend() const { return m_vals.cend(); }
    // reverse iterators 
    auto rbegin(){ return m_vals.rbegin(); }
    auto rend(){ return m_vals.rend(); }
    auto rbegin() const { return m_vals.rbegin(); }
    auto rend() const { return m_vals.rend(); }
    // reverse citerators
    auto crbegin() const { return m_vals.crbegin(); }
    auto crend() const { return m_vals.crend(); }  
    // set vals / size from std vector or iterators? 
};

// quick helper to wrap std::make_shared<Mesh1D> or other Mesh1D derived types  
template<typename Mesh_t=Mesh1D, typename... Args> 
auto make_mesh(Args... args)
{
  static_assert(std::is_base_of<Mesh1D,Mesh_t>::value, "make_mesh() requires T in shared_ptr<T> to be derived from Mesh1D.");
  return std::make_shared<Mesh_t>(args...); 
}

#endif // Mesh.hpp