// All.hpp
//
//
//
// JAF 12/8/2025

#ifndef DIFFOPS_ALL_H
#define DIFFOPS_ALL_H 

#include "FdmPlugin.hpp"
#include "Solver1D.hpp"

#include "BoundaryCond.hpp" 
#include "BCs/BCPair.hpp"
#include "BCs/BCLeftRight.hpp"
#include "BCs/Dirichlet.hpp"
#include "BCs/Neumann.hpp"
#include "BCs/Robin.hpp"

#include "DiffOps/NthDerivOp.hpp"
#include "CoeffOpBase.hpp"
#include "CoeffOps/AutonomousCoeff.hpp"
#include "CoeffOps/TimeDepCoeff.hpp"

#include "../LinOps/All.hpp"
namespace Fds{
  using namespace LinOps; 
} // end namespace Fds 

#endif 