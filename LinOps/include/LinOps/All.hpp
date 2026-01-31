// All.hpp
// 
// header file to include full linop framework
// 
// JAF 12/5/2025 

#ifndef LINOP_ALL_H
#define LINOP_ALL_H

#include "Mesh.hpp"
#include "Discretization.hpp"

#include "LinOpTraits.hpp"
#include "LinearOpBase.hpp"
#include "LinOpExpr.hpp"

#include "Operators/IOp.hpp"
#include "Operators/RandLinOp.hpp"
#include "Operators/DirectionalRandOp.hpp"
#include "Operators/DiffOps/NthDerivOp.hpp"
#include "Operators/DiffOps/DirectionalNthDerivOp.hpp" 
#include "Operators/CoeffOps/AutonomousCoeff.hpp"
#include "Operators/CoeffOps/TimeDepCoeff.hpp"

#endif // All.hpp
