// RandLinOp.hpp
//
// header file for a linop wrapper of a random matrix
//
// JAF 12/6/2025

#ifndef RANDLINOP_H
#define RANDLINOP_H

#include<eigen3/Eigen/Core>

#include "../LinearOpBase.hpp"

class RandLinOp : public LinOpBase<RandLinOp>
{
  using MeshPtr_t = std::shared_ptr<Mesh1D>;
  private:
    Eigen::MatrixXd m_Mat; 
  public: 
    // Constructors --------------------------
    RandLinOp(std::size_t s=1)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s); 
    }
    RandLinOp(MeshPtr_t m)
      : LinOpBase(m) 
    {
      set_mesh(m);  
    }
    // Destructors ----------------------------
    ~RandLinOp()=default;
    // member funcs
    Eigen::MatrixXd& GetMat() { return m_Mat; };
    const Eigen::MatrixXd& GetMat() const { return m_Mat; };
    Discretization1D apply(const Discretization1D& d_arr) const 
    { 
      Discretization1D result(d_arr.mesh()); 
      std::cout << result.size() << ": " << result.at(0) << std::endl;
      result = m_Mat * d_arr.values(); 
      return result; 
    }; 
    void set_mesh(MeshPtr_t m)
    {
      m_Mat = Eigen::MatrixXd::Random(m->size(),m->size()); 
    };
    void resize(std::size_t s)
    {
      m_Mat = Eigen::MatrixXd::Random(s,s);
    };
    void resize(std::size_t s1, std::size_t s2)
    {
      m_Mat = Eigen::MatrixXd::Random(s1,s2);
    }
};

#endif // RandLinOp.hpp