// BoundaryCond.hpp
//
//
//
// JAF 12/8/2025

#ifndef BOUNDARYCOND_H
#define BOUNDARYCOND_H

#include<eigen3/Eigen/Core>
#include "../LinOps/Discretization.hpp"
#include "../LinOps/Mesh.hpp"

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
class BoundaryCond
{
  protected:
    // types
    using MeshPtr_t = std::shared_ptr<Mesh1D>; 
  protected:  
    // member data 
    double m_current_time;
  public:
    // Constructors ---------------------------------------------
    BoundaryCond(): m_current_time(0.0){}; 
    // Destructors ----------------------------------------------
    virtual ~BoundaryCond()=default; 
    // Member Funcs ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTime(double t){m_current_time=t;};

    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(Eigen::MatrixXd& Mat, const MeshPtr_t& mesh)const{}; 
    virtual void SetStencilR(Eigen::MatrixXd& Mat, const MeshPtr_t& mesh)const{};
    
    // change the first/last entries in an impicit solution vector 
    virtual void SetImpSolL(Discretization1D& Sol, const MeshPtr_t& mesh)const{};
    virtual void SetImpSolR(Discretization1D& Sol, const MeshPtr_t& mesh)const{};

    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol, const MeshPtr_t& mesh)const{};
    virtual void SetSolR(Discretization1D& Sol, const MeshPtr_t& mesh)const{};
};

#endif // BoundaryCond.hpp