// All.hpp
//
// Include all headers from FDStencilsXD 
//
// JAF 1/10/2026 

#ifndef FDSTENCILSXD_ALL_H
#define FDSTENCILSXD_ALL_H

#include "FdmPluginXD.hpp"
#include "SolverXD.hpp"

#include "BoundaryCondXD.hpp"
#include "BCs/BCListXD.hpp" 

#include "DiffOps/DirectionalNthDerivOp.hpp"
#include "CoeffOpBaseXD.hpp"
#include "CoeffOps/AutonomousCoeffXD.hpp"
#include "CoeffOps/TimeDepCoeffXD.hpp"

#include "../LinOpsXD/All.hpp"
namespace Fds{
  using namespace LinOps; 
} // end namespace Fds 

#endif // end All.hpp 

