// BCPair.hpp
//
//
//
// JAF 12/8/2025

#ifndef BCPAIR_H
#define BCPAIR_H

#include<memory>
#include<cstdint>
#include<Eigen/Core>
#include<utility> // std::pair
#include "../../LinOps/Mesh.hpp"
#include "../../LinOps/Discretization.hpp"
#include "BCLeftRight.hpp"

namespace Fds{
using namespace LinOps; 

// type alias 
using MatrixStorage_t = Eigen::SparseMatrix<double, Eigen::RowMajor>; 

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
class BCPair : public IBoundaryCond
{
  protected:  
    // TypeDefs  
    using StridedRef = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

  public:
    // Member Data -----------------------------------------------------------
    std::pair<std::shared_ptr<IBCLeft>, std::shared_ptr<IBCRight>> pair;   
    // double m_current_time;
    
  public:
    // Constructors + Destructor =================================================
    BCPair() = delete;
    BCPair(std::shared_ptr<IBCLeft> l, std::shared_ptr<IBCRight> r)
      : pair(std::move(l), std::move(r))
    {}; 
    BCPair(std::pair<std::shared_ptr<IBCLeft>, std::shared_ptr<IBCRight>> p)
      : pair(std::move(p))
    {}; 
    BCPair(const BCPair& other)=default; 
    // destructor
    virtual ~BCPair()=default; 
    // Member Functions ==================================================================
    // Default Implemented ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTime(double t) final 
    {
      pair.first->SetTimeL(t); 
      pair.second->SetTimeR(t); 
    }

    // Pure Virtual ---------------------------------------------------
    // change last (right boundary) row of the fdm stencil matrix
    virtual void SetStencil(MatrixStorage_t& Mat, const std::shared_ptr<const Mesh1D>& mesh) const final
    {
      pair.first->SetStencilL(Mat, mesh); 
      pair.second->SetStencilR(Mat, mesh); 
    }
    
    // change the last entries in an impicit solution vector 
    virtual void SetImpSol(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh) const final
    {
      pair.first->SetImpSolL(Sol, mesh); 
      pair.second->SetImpSolR(Sol, mesh); 
    }

    // change the last (right boundary) entry of a vector  
    virtual void SetSol(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh) const final
    {
      pair.first->SetSolL(Sol, mesh); 
      pair.second->SetSolR(Sol, mesh); 
    }
};

} // end namespace Fds 

#endif // BCPair.hpp