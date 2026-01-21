// main.cpp
//
//
//
// JAF 12/8/2025

#include<cstdint>
#include<iostream>
#include<iomanip>
#include<vector>
#include<tuple>
#include<Eigen/Dense>

#include "FDStencils/All.hpp" // must be first for plugin macro.
// #include "FDStencils/Interp1D.hpp" // Interp1D class 
#include "LinOps/All.hpp" 
#include "FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/All.hpp"
#include "TExprs/All.hpp" 

#include "Utilities/PrintVec.hpp"
#include "Utilities/BumpFunc.hpp"


using std::cout, std::endl;

int main()
{
  // iomanip 
  std::cout << std::setprecision(4); 
  
  // Domain MeshXD
  auto my_mesh = LinOps::make_meshXD(0.0,4.0,5,3); 

  // DiscretizationXD + ICs
  LinOps::DiscretizationXD my_vals; 
  auto f = [](double x, double y, double z){ return std::sqrt(x*x + y*y + z*z); }; 
  my_vals.set_init(my_mesh,f); 

  // printing out flat index -> sub dim index 
  // iterate through flat index 
  for(std::size_t flat_i=0; flat_i < my_mesh->sizes_product(); flat_i++)
  {
    std::cout << flat_i << ": ";  

    // iterate through dims 0,1,2 
    for(std::size_t ith_dim=0; ith_dim < my_mesh->dims(); ith_dim++)
    {
      std::size_t s1 = (flat_i % my_mesh->sizes_middle_product(0,ith_dim+1)); 
      std::size_t s2 = s1 / my_mesh->sizes_middle_product(0,ith_dim); 
      std::cout << s2 << ", "; 
    }
    std::cout << "\n"; 
  }
};
