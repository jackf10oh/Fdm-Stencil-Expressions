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

// #include "DiffOps/DiffOps/experimental_NthDerivOp.hpp" 
#include "DiffOps/All.hpp" // must include first for plugin to take effect over linops?
#include "LinOps/All.hpp" 
#include "LinOpsXD/All.hpp"
#include "Utilities/SparseDiagExpr.hpp"

#include "DiffOps/CoeffOps/AutonomousCoeff.hpp"

using std::cout, std::endl;

auto lam00 = [](){return 0.0;}; 
auto lam01 = [](double x){return x*x - x + 1.5;}; 
auto lam02 = [](double x, double y){return std::sqrt(x*x + y*y);}; 
auto lam03 = [](double x, double y, double z){return std::sqrt(x*x + y*y + z*z);}; 

int main()
{
  // // iomanip 
  std::cout << std::setprecision(2); 

  // mesh assembly. 2 dims 
  auto my_mesh = make_mesh(0.0, 5, 6);
  
  AutonomousCoeff c(lam01, my_mesh); 

  cout << c.GetMat() << endl; 

  // double val = c.GetScalar(); 

};

  // auto A = DirectionalRandOp(my_meshes, 0); 
  // auto B = DirectionalRandOp(my_meshes, 1); 
  // // auto expr = A+B; 

  // cout << "---------- A --------------" << endl << A.GetMat().toDense() << endl; 
  // cout << "---------- B --------------" << endl << B.GetMat().toDense() << endl; 
  // cout << "---------- A+B --------------" << endl << (A+B).GetMat() << endl; 
  // cout << "---------- A-B --------------" << endl << (A-B).GetMat() << endl; 
  // cout << "---------- -A --------------" << endl << (-A).GetMat() << endl; 
  // cout << "---------- 2.0*A --------------" << endl << (2.0*A).GetMat() << endl; 


  // discretization.
  // DiscretizationXD my_vals; 
  // my_vals.set_init(my_meshes, lam02); 

  // print as flat vector 
  // print_vec(my_vals.values(), "flattened vals"); 

  // print as a 2d matrix 
  // print_mat(my_vals.OneDim_views(0), "Mapped 2D"); 

  // from a given "slice" that looks like 2d. print the matrix 
  // std::size_t ith_slice = 1; 
  // std::size_t offset = ith_slice * my_meshes->sizes_middle_product(0,2); 
  // cout << Eigen::Map<Eigen::MatrixXd>(my_vals.values().data() + offset, my_meshes->dim_size(0), my_meshes->dim_size(1)) << endl; 

  // auto views = my_vals.OneDim_views(2);
  // for(auto& v : views) print_vec(v); 

  // BoundaryCondXD 
  // BoundaryCondXD bc_list; 
  // using pair_t = std::pair<BcPtr_t,BcPtr_t>; 
  // pair_t p1 = {make_dirichlet(1.0), make_dirichlet(2.0)}; 
  // pair_t p2 = {make_dirichlet(3.0), make_dirichlet(4.0)};
  // // pair_t p2 = {make_robin(1,1,3.0), make_robin(1,1,4.0)};
  // bc_list.m_bc_list = {p1,p2}; 

  // bc_list.SetSol(my_vals, my_meshes); 
  // bc_list.SetImpSol(my_vals, my_meshes); 

  // // how to iterate through dynamic multi dim meshes?
  // std::size_t end = my_meshes->sizes_product();  
  // std::vector<double> coords(my_meshes->dims()); 
  // for(std::size_t flat_i=0; flat_i<end; flat_i++){
  //   for(std::size_t dim=0; dim<my_meshes->dims(); dim++){
  //     std::size_t dim_i = flat_i;
  //     dim_i /= my_meshes->sizes_middle_product(0, dim); 
  //     dim_i %= my_meshes->dim_size(dim); 
  //     cout << dim_i; 
  //     cout << ", ";
  //   }
  //   cout << endl; 
  // } 

// 1st dim offset : 4*n 
// 2nd dim offset : (n / my_meshes->sizes_middle_product(0,ith_dim))

/*
0 -> 0,4,8,12 
1 -> 1,5,9,13
2 -> 2,6,10,14 
3 -> 3,7,11,15  

4 -> 0+16, 4+16, 8+16, 12+16
5 -> 1+16,5+16,9+16,13+16
6 -> 2+16,6+16,10+16,14+16 
7 -> 3+16,7+16,11+16,15+16 

*/
// printed output 
/*
0, 0, 0, 
1, 0, 0, 
2, 0, 0, 
3, 0, 0, 
0, 1, 0, 
1, 1, 0, 
2, 1, 0, 
3, 1, 0, 
0, 2, 0, 
1, 2, 0, 
2, 2, 0, 
3, 2, 0, 
0, 3, 0, 
1, 3, 0, 
2, 3, 0, 
3, 3, 0, 
0, 0, 1, 
1, 0, 1, 
2, 0, 1, 
3, 0, 1, 
0, 1, 1, 
1, 1, 1, 
2, 1, 1, 
3, 1, 1, 
0, 2, 1, 
1, 2, 1, 
2, 2, 1, 
3, 2, 1, 
0, 3, 1, 
1, 3, 1, 
2, 3, 1, 
3, 3, 1, 
0, 0, 2, 
1, 0, 2, 
2, 0, 2, 
3, 0, 2, 
0, 1, 2, 
1, 1, 2, 
2, 1, 2, 
3, 1, 2, 
0, 2, 2, 
1, 2, 2, 
2, 2, 2, 
3, 2, 2, 
0, 3, 2, 
1, 3, 2, 
2, 3, 2, 
3, 3, 2, 
0, 0, 3, 
1, 0, 3, 
2, 0, 3, 
3, 0, 3, 
0, 1, 3, 
1, 1, 3, 
2, 1, 3, 
3, 1, 3, 
0, 2, 3, 
1, 2, 3, 
2, 2, 3, 
3, 2, 3, 
0, 3, 3, 
1, 3, 3, 
2, 3, 3, 
3, 3, 3, 
*/ 

/*   
  // boundary conditions 
  auto set_dim_boundaries = [&](std::size_t dim, BcPtr_t bc){
    auto views = my_vals.OneDim_views(dim); 
    std::size_t s1 = views.size(); 
    // iterate through views that look like Mesh1D
    for(std::size_t i=0; i<s1; i++){
      bool set_by_low_dim = false; 
      // determine if it has been set by lower dimension 
      std::size_t s2 = 1; 
      for(int dim_i=0; dim_i<dim; dim_i++){
        std::size_t s3 = my_vals.dim_size(dim_i); 
        // in first group of values
        if((i/s2)%s3 == 0) set_by_low_dim = true;  
        // in last group of values
        if((i/s2)%s3 == s3-1) set_by_low_dim = true;  
        
        // next check uses a large bucket 
        s2 *= s3; 
      }
      // if it hasn't, use BC on it. 
      if(!set_by_low_dim){
        bc->SetSolL(views[i], my_meshes->GetMesh(dim)); 
        bc->SetSolR(views[i], my_meshes->GetMesh(dim)); 
      }
    }
  }; 

  set_dim_boundaries(0, make_dirichlet(1.0)); 
  set_dim_boundaries(1, make_dirichlet(2.0)); 
  // set_dim_boundaries(2, make_dirichlet(3.0)); 
*/