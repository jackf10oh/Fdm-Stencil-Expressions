// Robin.hpp 
//
// robin(third type) boundary conditions 
// alpha*U + beta*Ux = target_val
//
// JAF 12/8/2025

#ifndef ROBINBCS_H
#define ROBINBCS_H 

#include "../BoundaryCond.hpp"
#include "../Utilities/FornbergCalc.hpp"

class RobinBC : public IBoundaryCond
{
  public:  
    // member data 
    double boundary_target;
    double val_coeff;
    double deriv_coeff;

  public:
    // Constructors ---------------------------------------------
    // a*U + b*Ux = target
    RobinBC(double a=1.0, double b=0.0, double target=0.0) 
      : val_coeff(a), deriv_coeff(b), boundary_target(target) 
    {}; 
    // Destructors ----------------------------------------------
    virtual ~RobinBC()=default; 
    // Member Funcs ----------------------------------------------
    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencilL(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.topRows(1) *= 0;
      // first order derivative approximation 
      double h = mesh->operator[](1) - mesh->operator[](0);  
      Mat.coeffRef(0,0)=  val_coeff + deriv_coeff*(-1.0/h);
      Mat.coeffRef(0,1)=  deriv_coeff*(1.0/h);
    }; 
    virtual void SetStencilR(MatrixStorage_t& Mat, const MeshPtr_t& mesh) const override
    {
      Mat.bottomRows(1) *= 0; 
      // first order derivative approximation 
      double h = mesh->operator[](mesh->size()-1) - mesh->operator[](mesh->size()-2);  
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-2)= deriv_coeff*(-1.0/h);
      Mat.coeffRef(Mat.rows()-1, Mat.cols()-1)=  val_coeff + deriv_coeff*(1.0/h);
    };

    virtual void SetImpSolL(Discretization1D& Sol, const MeshPtr_t& mesh) const override
    {Sol.at(0) = boundary_target;};
    virtual void SetImpSolR(Discretization1D& Sol, const MeshPtr_t& mesh) const override
    {Sol.at(Sol.size()-1) = boundary_target;};
    
    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSolL(Discretization1D& Sol, const MeshPtr_t& mesh) const override 
    { 
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights((*mesh)[0], mesh->begin(), mesh->begin()+3, 1); 

      // solve the equation target = a*(S[0]) + b * ( W[0]*S[0] + W[1]*S[1] + W[2]*S[2] ) 
      // for the target value S[0] 
      double target = boundary_target; 
      target -= deriv_coeff * weights[1]*Sol[1]; 
      target -= deriv_coeff * weights[2]*Sol[2];
      target /= val_coeff + deriv_coeff * weights[0]; 

      // assign to Sol reference
      Sol[0] = target;  
      // void return type
    };
    virtual void SetSolR(Discretization1D& Sol, const MeshPtr_t& mesh) const override 
    {
      // if(Sol.size()<3 || mesh->size()<3) throw std::runtime_error("Discretization1D or Mesh1D size too small!(must be >= 3)"); 

      // up to 3 nodes, up to 1st order deriv
      FornCalc calc(3,1);

      // get forward finite difference weights for Sol[0], Sol[1], Sol[2] 
      auto weights = calc.GetWeights((*mesh)[mesh->size()-1], mesh->end()-3, mesh->end(), 1); 

      // solve the equation target = a*(S[N-1]) + b * ( W[0]*S[N-3] + W[1]*S[N-2] + W[2]*S[N-1] ) 
      // for the target value S[N-1] 
      double target = boundary_target; 
      target -= deriv_coeff * weights[0]*Sol[Sol.size()-3]; 
      target -= deriv_coeff * weights[1]*Sol[Sol.size()-2]; 
      target /= val_coeff + deriv_coeff*weights[2]; 

      // assign to Sol reference
      Sol[Sol.size()-1] = target;  
      // void return type
    };
};

#endif // Robin.hpp