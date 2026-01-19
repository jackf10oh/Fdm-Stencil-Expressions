// Dirichlet.hpp 
//
// dirichlet boundary value conditions of the form 
// U(0) = c for some constant c 
//
// JAF 12/8/2025

#ifndef DIRICHLETBCS_H
#define DIRICHLETBCS_H 

#include "BCLeftRight.hpp" 

namespace Fds{

class DirichletBC : public IBCLeft, public IBCRight 
{
  public:  
    // member data 
    double boundary_val;

  public:
    // Constructors ---------------------------------------------
    DirichletBC(double val_init=0.0) : boundary_val(val_init), IBCLeft(), IBCRight(){}; 
    // Destructors ----------------------------------------------
    virtual ~DirichletBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(MatrixStorage_t& Mat, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override
    {
      Mat.topRows(1) *= 0; Mat.coeffRef(0,0)=1;
    }
    virtual void SetStencilR(MatrixStorage_t& Mat, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override
    {
      Mat.bottomRows(1) *= 0; Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=1;
    }

    // change the first/last (left/right boundary) entry of a vector to implicit solution   
    virtual void SetImpSolL(StridedRef Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override
    {Sol[0] = boundary_val;}
    virtual void SetImpSolR(StridedRef Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override
    {Sol[Sol.size()-1] = boundary_val;}
    
    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(StridedRef Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override 
    { Sol[0] = boundary_val;}
    virtual void SetSolR(StridedRef Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh) const override 
    {Sol[Sol.size()-1] = boundary_val;}
};

template<typename... Args>
auto make_dirichlet(Args... args){ return std::make_shared<DirichletBC>(args...); }; 

} // end namespace Fds 

#endif