// BoundaryCondXD.hpp
//
//
//
// JAF 1/2/2025 

#ifndef BOUNDARYCONDXD_H
#define BOUNDARYCONDXD_H 

#include<iostream>
#include<vector>
#include<tuple>
#include "../LinOpsXD/MeshXD.hpp"
#include "../LinOpsXD/DiscretizationXD.hpp"

namespace Fds{
using namespace LinOps; 

class IBoundaryCondXD
{
  public:
    // Constructors + Destructors ========================
    IBoundaryCondXD()=default; 
    IBoundaryCondXD(const IBoundaryCondXD& other)=default;
    // destructor 
    virtual ~IBoundaryCondXD()=default; 

    // Member Functions ============================================
    // Must Implement ---------------------------------------------
    // set a DiscretizationXD to an explicit solution 
    virtual void SetSol(DiscretizationXD& Sol, const std::shared_ptr<const MeshXD>& mesh) const =0; 

    // set a DiscretizationXD to an implicit solution 
    virtual void SetImpSol(DiscretizationXD& Sol, const std::shared_ptr<const MeshXD>& mesh) const =0; 

    // set a Matrixs' row according to m_bc_list. making it an implicit stencil  
    virtual void SetStencilImp(MatrixStorage_t& Mat, const std::shared_ptr<const MeshXD>& mesh) const =0; 

    // Default Implemented ---------------------------------------------
    virtual void SetTime(double t){}; 
};

} // end namespace Fds 

#endif // BoundaryCondXD.hpp 