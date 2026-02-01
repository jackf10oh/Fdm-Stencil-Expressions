// OStepBase.hpp 
//
// CRTP bases for 1D / XD Outside Step. 
// Boundary Conditions will be viewed as another type of OStep in the new framework 
//
// JAF 1/31/2026 

#ifndef OSTEPBASE_H
#define OSTEPBASE_H

#include<LinOps/LinearOpBase.hpp> // LinOps::MatrixStorage_t 

namespace OSteps{
using StridedRef_t = Eigen::Ref<Eigen::VectorXd, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;
using LinOps::MatrixStorage_t; 
using LinOps::Mesh1D_SPtr_t; 
using LinOps::MeshXD_SPtr_t;

template<typename DERIVED>
class OStepBase1D
{
  public:
    // Type Defs -----------------
    typedef struct{} is_1dim_ostep_tag; 

    // Member Funcs ------------------------------------------------
    void BeforeExpStep(double t, const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->BeforeExpStep(t,mesh, Mat, Sol); 
    }
    void AfterExpStep(double t, const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->AfterExpStep(t,mesh, Mat, Sol); 
    }
    void BeforeImpStep(double t, const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->BeforeImpStep(t,mesh, Mat, Sol); 
    }
    void AfterImpStep(double t, const Mesh1D_SPtr_t& mesh, LinOps::MatrixStorage_t& Mat, StridedRef_t Sol) const 
    {
      static_cast<const DERIVED*>(this)->AfterImpStep(t,mesh, Mat, Sol); 
    }
}; 

} // end namesace OSteps

#endif // OStepBase.hpp