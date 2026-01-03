// BoundaryCondXD.hpp
//
//
//
// JAF 1/2/2025 

#ifndef BOUNDARYCONDXD_H
#define BOUNDARYCONDXD_H 

#include<vector>
#include<tuple>
#include "../DiffOps/BoundaryCond.hpp"
#include "MeshXD.hpp"

struct BoundaryCondXD
{
  // member data. list of boundary conditions. 1 per Dimension 
  std::vector<std::pair<BcPtr_t>> m_bc_arr; 

  // set a DiscretizationXD to an explicit solution 

  // set a DiscretizationXD to an implicit solution 

  // get a Matrix that to apply as a mask over XD fdm stencils 
  // MatrixStorage_t BoundaryRowsMask()
};

#endif // BOundaryCondXD.hpp