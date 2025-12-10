// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Sparse>

#include "../LinearOpBase.hpp"

class IOp : public LinOpBase<IOp>
{
  using MeshPtr_t = std::shared_ptr<Mesh1D>;
  private:
    Eigen::SparseMatrix<double> m_Mat; 
  public: 
    // Constructors --------------------------
    IOp(MeshPtr_t m=nullptr)
    {
      set_mesh(m);
    } 
    // Destructors ----------------------------
    ~IOp()=default;
    // member funcs
    auto GetMat() { return m_Mat; };
    const auto& GetMat() const { return m_Mat; };
    Discretization1D apply(const Discretization1D& d_arr) const 
    { 
      return d_arr;
    }; 
    void set_mesh(MeshPtr_t m)
    {
      // store m into inherited m_mesh_ptr data member
      // checks null type
      LinOpBase::set_mesh(m);

      // check null type
      if(m_mesh_ptr==nullptr) return;
      
      // pointer isn't null -> resize m_Mat
      m_Mat.resize(m_mesh_ptr->size(), m_mesh_ptr->size()); 
      m_Mat.setIdentity(); 
    };
    void resize(std::size_t s=0)
    {
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };
    void resize(std::size_t s1,std::size_t s2) 
    {
      m_Mat.resize(s1,s2); 
      m_Mat.setIdentity(); 
    };
};

#endif