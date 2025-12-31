// Dirichlet.hpp 
//
//
//
// JAF 12/8/2025

#ifndef DIRICHLETBCS_H
#define DIRICHLETBCS_H 

#include "../BoundaryCond.hpp"

class NeumannBC : public IBoundaryCond
{
  public:  
    // member data 
    double boundary_flux;

  public:
    // Constructors ---------------------------------------------
    NeumannBC(double val_init=0.0) : boundary_flux(val_init){}; 
    // Destructors ----------------------------------------------
    virtual ~NeumannBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.topRows(1) *= 0;
      // first order derivative approximation 
      double h = mesh->operator[](1) - mesh->operator[](0);  
      Mat.coeffRef(0,0)= -1.0/h;
      Mat.coeffRef(0,1)=  1.0/h;
    }; 
    virtual void SetStencilR(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.bottomRows(1) *= 0; 
      // first order derivative approximation 
      double h = mesh->operator[](mesh->size()-1) - mesh->operator[](mesh->size()-2);  
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-2)= -1.0/h;
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=  1.0/h;
    };

    virtual void SetImpSolL(Discretization1D& Sol, const MeshPtr_t& mesh) const override
    {Sol.at(0) = boundary_flux;};
    virtual void SetImpSolR(Discretization1D& Sol, const MeshPtr_t& mesh) const override
    {Sol.at(Sol.size()-1) = boundary_flux;};
    
    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol, const MeshPtr_t& mesh) const override 
    { Sol.at(0) = boundary_flux;};
    virtual void SetSolR(Discretization1D& Sol, const MeshPtr_t& mesh) const override 
    {Sol.at(Sol.size()-1) = boundary_flux;};
};

#endif