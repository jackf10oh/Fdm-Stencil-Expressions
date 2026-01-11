// All.hpp
//
// Include all headers from FDStencilsXD 
//
// JAF 1/10/2026 

#ifndef FDSTENCILSXD_ALL_H
#define FDSTENCILSXD_ALL_H

// #include "FdmPluginXD.hpp"
// TODO: Solver in XD  

#include "BoundaryCondXD.hpp"
#include "BCs/BCListXD.hpp" 

#include "DiffOps/DirectionalNthDerivOp.hpp"
// TODO: coeffops in XD 

#include "../LinOpsXD/All.hpp"
namespace Fds{
  using namespace LinOps; 
} // end namespace Fds 

#endif // end All.hpp 

