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

// forward declaration -> type alias 
class BoundaryCond; 
using BcPtr_t = std::shared_ptr<BoundaryCond>; 
// using BcPtr_t = BoundaryCond*; // unlikely to use. Only going to assign to boundary conditions once 
using MatrixStorage_t = Eigen::MatrixXd; 

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
class BoundaryCond
{
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
    virtual void SetStencilL(MatrixStorage_t& Mat, const MeshPtr_t& mesh)const{}; 
    virtual void SetStencilR(MatrixStorage_t& Mat, const MeshPtr_t& mesh)const{};
    
    // change the first/last entries in an impicit solution vector 
    virtual void SetImpSolL(Discretization1D& Sol)const{};
    virtual void SetImpSolR(Discretization1D& Sol)const{};

    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol)const{};
    virtual void SetSolR(Discretization1D& Sol)const{};
};

#endif // BoundaryCond.hpp