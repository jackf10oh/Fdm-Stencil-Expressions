// BoundaryCond.hpp
//
//
//
// JAF 12/8/2025

#ifndef BOUNDARYCOND_H
#define BOUNDARYCOND_H

#include<memory>
#include<cstdint>
#include<Eigen/Core>
#include "../LinOps/Discretization.hpp"
#include "../LinOps/Mesh.hpp"

// forward declaration -> type alias 
class IBoundaryCond; 
using MatrixStorage_t = Eigen::SparseMatrix<double, Eigen::RowMajor>; 
using BcPtr_t = std::shared_ptr<IBoundaryCond>; 

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
class IBoundaryCond
{
  protected:  
    // member data 
    double m_current_time;
  public:
    // Constructors ---------------------------------------------
    IBoundaryCond(): m_current_time(0.0){}; 
    // Destructors ----------------------------------------------
    virtual ~IBoundaryCond()=default; 
    // Member Funcs ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTime(double t){m_current_time=t;};

    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(Eigen::Ref<MatrixStorage_t> Mat, const MeshPtr_t& mesh)const=0; 
    virtual void SetStencilR(Eigen::Ref<MatrixStorage_t> Mat, const MeshPtr_t& mesh)const=0;
    
    // change the first/last entries in an impicit solution vector 
    virtual void SetImpSolL(Discretization1D& Sol, const MeshPtr_t& mesh)const=0;
    virtual void SetImpSolR(Discretization1D& Sol, const MeshPtr_t& mesh)const=0;

    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol, const MeshPtr_t& mesh)const=0;
    virtual void SetSolR(Discretization1D& Sol, const MeshPtr_t& mesh)const=0;
};

#endif // BoundaryCond.hpp