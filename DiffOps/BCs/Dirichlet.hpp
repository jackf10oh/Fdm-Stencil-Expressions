// Dirichlet.hpp 
//
//
//
// JAF 12/8/2025

#ifndef DIRICHLETBCS_H
#define DIRICHLETBCS_H 

#include "../BoundaryCond.hpp"

class DirichletBC : public BoundaryCond
{
  public:  
    // member data 
    double boundary_val;

  public:
    // Constructors ---------------------------------------------
    DirichletBC(double val_init=0.0) : boundary_val(val_init){}; 
    // Destructors ----------------------------------------------
    virtual ~DirichletBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.topRows(1) *= 0; Mat.coeffRef(0,0)=1;
    }; 
    virtual void SetStencilR(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.bottomRows(1) *= 0; Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=1;
    };

    virtual void SetImpSolL(Discretization1D& Sol) const override
    {Sol.at(0) = boundary_val;};
    virtual void SetImpSolR(Discretization1D& Sol) const override
    {Sol.at(Sol.size()-1) = boundary_val;};
    
    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol) const override 
    { Sol.at(0) = boundary_val;};
    virtual void SetSolR(Discretization1D& Sol)const override 
    {Sol.at(Sol.size()-1) = boundary_val;};
};

#endif