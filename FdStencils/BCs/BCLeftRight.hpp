// BCLeftRight.hpp
//
//
//
// JAF 1/11/2025 

#ifndef BCLEFTRIGHT_H
#define BCLEFTRIGHT_H 

#include<memory>
#include<cstdint>
#include<Eigen/Core>
#include<utility> // std::pair 
#include "../../LinOps/Mesh.hpp"
#include "../../LinOps/Discretization.hpp"

namespace Fds{
using namespace LinOps; 

using MatrixStorage_t = Eigen::SparseMatrix<double, Eigen::RowMajor>; 
using StridedRef = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

class IBCLeft
{
  protected:
    // Member Data ---------------------------------------------------------
    double m_current_time;
  public:
    // Constructors + Destructor =============================================
    IBCLeft() : m_current_time(0.0){};
    IBCLeft(double t) : m_current_time(t){}; 
    IBCLeft(const IBCLeft& other)=default; 
    // destructor
    virtual ~IBCLeft()=default; 

    // Member Functions ==================================================================
    // Default Implemented ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTimeL(double t){m_current_time=t;};

    // Pure Virtual ---------------------------------------------------
    // change first (left boundary) row of the fdm stencil matrix
    virtual void SetStencilL(MatrixStorage_t& Mat, const std::shared_ptr<const Mesh1D>& mesh)const=0; 
    
    // change the first entries in an impicit solution vector 
    virtual void SetImpSolL(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh)const=0;

    // change the first (left boundary) entry of a vector  
    virtual void SetSolL(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh)const=0;
};

class IBCRight
{
  protected:
    // Member Data ---------------------------------------------------------
    double m_current_time;
  public:
    // Constructors + Destructor =============================================
    IBCRight() : m_current_time(0.0){};
    IBCRight(double t) : m_current_time(t){}; 
    IBCRight(const IBCRight& other)=default; 
    // destructor
    virtual ~IBCRight()=default; 

    // Member Functions ==================================================================
    // Default Implemented ----------------------------------------------
    // Set Time of current boundary condition
    virtual void SetTimeR(double t){m_current_time=t;};

    // Pure Virtual ---------------------------------------------------
    // change last (right boundary) row of the fdm stencil matrix
    virtual void SetStencilR(MatrixStorage_t& Mat, const std::shared_ptr<const Mesh1D>& mesh)const=0; 
    
    // change the last entries in an impicit solution vector 
    virtual void SetImpSolR(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh)const=0;

    // change the last (right boundary) entry of a vector  
    virtual void SetSolR(StridedRef Sol, const std::shared_ptr<const Mesh1D>& mesh)const=0;
};

} // end namespace Fds 

#endif // BCLeftRight.hpp 



