// BCPair.hpp
//
// 
//
// JAF 12/8/2025

#ifndef BCPAIR_H
#define BCPAIR_H

#include "../OStepBase.hpp"

namespace OSteps{

// Base Class for Boundary Conditions. All operators make no changes to stencil / solution 
template<typename LBC_T,typename RBC_T>
class BCPair : public OStepBase1D<BCPair<LBC_T,RBC_T>>
{
  public:
    // Member Data -----------------------------------------------------------
    typename std::remove_reference<LBC_T>::type m_left; 
    typename std::remove_reference<RBC_T>::type m_right; 
    
  public:
    // Constructors + Destructor =================================================
    BCPair() = delete;
    BCPair(LBC_T l, RBC_T r)
      : m_left(l),m_right(r)
    {}; 
    BCPair(const BCPair& other)=default; 
    // destructor
    virtual ~BCPair()=default; 

    // Member Functions ==================================================================

    void BeforeExpStep(double t,  const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {/* boundary conditions only alter Solution on Explicit steps...*/}

    void AfterExpStep(double t,  const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      m_left.SetSolL(t, mesh, Sol); 
      m_right.SetSolR(t, mesh, Sol); 
    }

    void BeforeImpStep(double t,  const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      m_left.SetStencilL(t, mesh, Mat); 
      m_right.SetStencilR(t, mesh, Mat); 
      m_left.SetImpSolL(t, mesh, Sol); 
      m_right.SetImpSolR(t, mesh, Sol); 
    }
    void AfterImpStep(double t,  const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    { /* boundary conditions only alter Ax=b before Implicit step...*/} 
};

} // end namespace OSteps 

#endif // BCPair.hpp