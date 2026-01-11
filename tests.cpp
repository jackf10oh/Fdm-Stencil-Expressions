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

#include "Utilities/PrintVec.hpp"

#include "FDStencils/All.hpp" // must include first for plugin to take effect over linops?
#include "LinOps/All.hpp" 
#include "LinOpsXD/All.hpp"
#include "FDStencilsXD/DiffOps/DirectionalNthDerivOp.hpp"
#include "FDStencilsXD/BCs/BCListXD.hpp"

using std::cout, std::endl;

using namespace Fds; 
using namespace LinOps; 

auto lam00 = [](){return 0.0;}; 
auto lam01 = [](double x){return x*x - x + 1.5;}; 
auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

int main()
{
  // iomanip 
  std::cout << std::setprecision(2); 

  auto my_meshes = LinOps::make_meshXD(0.0,10.0, 11, 3);
  // auto my_meshes = LinOps::make_mesh(0.0,3.0, 4);

  LinOps::DiscretizationXD disc; 
  disc.set_init(my_meshes, lam03);  

  BCListXD bcs; 
  bcs.bc_list.emplace_back(make_dirichlet(1.0), make_dirichlet(1.0)); 
  bcs.bc_list.emplace_back(make_dirichlet(2.0), make_dirichlet(2.0)); 
  bcs.bc_list.emplace_back(make_neumann(1.0), make_neumann(1.0)); 

  bcs.SetSol(disc, my_meshes); 

  // view a 1D slice of disc 
  auto views = disc.OneDim_views(0);  
  cout << "1D -----------" << endl << views[0].transpose() << endl; 

  // view ith slice of disc that looks like 2D matrix 
  using Stride_t = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
  using MatView_t = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Stride_t>; 
  std::size_t offset = 1 * disc.sizes_middle_product(0,2); 
  MatView_t A(disc.values().data() + offset, disc.dim_size(0), disc.dim_size(1), Stride_t(11,1)); 

  cout << "2D -----------" << endl <<A << endl; 
  
  // LinOps::IOpXD I(my_meshes);

};

