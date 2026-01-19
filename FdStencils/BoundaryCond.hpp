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
#include "FdmPlugin.hpp" // MatrixStorage_t
#include "../LinOps/Mesh.hpp"

namespace Fds{

// forward declaration -> type alias 
class IBoundaryCond; 
using BcPtr_t = std::shared_ptr<IBoundaryCond>; 

using StridedRef_t = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
class IBoundaryCond
{
  private:
    // Member Data --------------------- 
    double m_current_time;
    
  public:
    // Constructors + Destructor ===========================================
    IBoundaryCond(): m_current_time(0.0){}; 
    // destructors 
    virtual ~IBoundaryCond()=default; 

    // Member Funcs ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTime(double t){m_current_time=t;};

    // change first/last (left/right boundary) row of the fdm stencil matrix
    virtual void SetStencil(MatrixStorage_t& Mat, const std::shared_ptr<const LinOps::Mesh1D>& mesh)const=0;
    
    // change the first/last entries in an impicit solution vector 
    virtual void SetImpSol(StridedRef_t Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh)const=0;

    // change the first/last (left/right boundary) entry of a vector  
    virtual void SetSol(StridedRef_t Sol, const std::shared_ptr<const LinOps::Mesh1D>& mesh)const=0;
};

} // end namespace Fds 

#endif // BoundaryCond.hpp