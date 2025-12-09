// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include<eigen3/Eigen/Core>

#include "../LinearOpBase.hpp"

class IOp : public LinOpBase<IOp>
{
  using MeshPtr_t = std::shared_ptr<Mesh1D>;
  private:
    Eigen::MatrixXd m_Mat; 
  public: 
    // Constructors --------------------------
    IOp(std::size_t s=1)
    {
      m_Mat = Eigen::MatrixXd::Identity(s,s); 
    }
    IOp(MeshPtr_t m)
      : LinOpBase(m) 
    {
      set_mesh(m);  
    }
    // Destructors ----------------------------
    ~IOp()=default;
    // member funcs
    Eigen::MatrixXd& GetMat() { return m_Mat; };
    const Eigen::MatrixXd& GetMat() const { return m_Mat; };
    Discretization1D apply(const Discretization1D& d_arr) const 
    { 
      Discretization1D result(d_arr.mesh()); 
      result = m_Mat * d_arr.values(); 
      return result; 
    }; 
    void set_mesh(MeshPtr_t m)
    {
      m_Mat = Eigen::MatrixXd::Identity(m->size(),m->size()); 
    };
    void resize(std::size_t s)
    {
      m_Mat = Eigen::MatrixXd::Identity(s,s);
    };
    void resize(std::size_t s1, std::size_t s2)
    {
      m_Mat = Eigen::MatrixXd::Identity(s1,s2);
    }
};

#endif