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
#include "LinOps/All.hpp" 
#include "FDStencilsXD/All.hpp" // likewise ...
#include "LinOpsXD/All.hpp"
#include "TExprs/All.hpp" 

#include "Utilities/PrintVec.hpp"
#include "Utilities/BumpFunc.hpp"

using std::cout, std::endl;



// using interp_t = ; 
// struct bar : private foo, public interp_t
// {
//   bar(double a1, double b1, double c1) 
//     : foo{.a=a1, .b=b1, .c=c1},
//     interp_t(foo::GetLhsExpr(), foo::GetRhsExpr()) 
//   {};
// };

int main()
{
  // iomanip 
  // std::cout << std::setprecision(4); 
  foo my_foo(1.0 ,0.0, 0.0); 
};

  // Domain MeshXD ----------------------------------------
  // auto my_mesh = LinOps::make_mesh(0.0,4.0,5); 

  // double r = 4.0;
  // std::size_t N = 5;
  // auto my_mesh = LinOps::make_meshXD(0.0,r,N, 2);

  // std::vector<std::pair<double,double>> end_points = {{0.0,1.0}, {0.0,5.0}}; 
  // std::vector<std::size_t> n_steps = {4,9}; 
  // auto my_mesh = LinOps::make_meshXD(end_points, n_steps);

  // std::cout << "my_mesh type: " << typeid(decltype(my_mesh)).name() << std::endl; 
  
  // std::shared_ptr<const LinOps::MeshXD> copied = make_meshXD(my_mesh); 

  // std::cout << "# Dims: " << copied->dims() << std::endl; 
  // std::cout << "dim(0) size: " << copied->dim_size(0) << std::endl; 
  // std::cout << "sizes prod: " << copied->sizes_product() << std::endl; 


  // ICs --------------
  // LinOps::DiscretizationXD my_vals;

  // Bump centered at (r/2,r2). Zero at x=0, x=r, y=0, y=r. 
  // BumpFunc bump_1d{.L=0.0, .R=r, .c=r/2, .h=5.0, .focus=15.0};
  // auto f_lambda = [&](double x, double y){ return bump_1d(x) * bump_1d(y); }; 
  // auto f_lambda = [](double x, double y){ return x + y; };
  // my_vals.set_init(my_mesh, f_lambda); 


  // Print Matrix
  // auto views = my_mesh->OneDim_views(my_vals.values(), 1); 
  // print_mat(views, "Matrix"); 

  // std::vector<double> coords(2); 
  // while(true)
  // {
  //   std::cout << "enter x: "; 
  //   std::cin >> coords[0]; 
  //   std::cout << "enter y: "; 
  //   std::cin >> coords[1];
  //   std::cout << "value at (x,y) : " << LinearInterp(coords, my_vals.values(), my_mesh)<<std::endl; 
  // }


  // print interpolations 
  // std::vector<double> coords(2); 
  // coords[0] = 2.0; 
  // for(std::size_t j=0; j<N; j++)
  // {
  //   coords[1] = j * (r/(N-1)); 
  //   std::cout << coords[0] << ", " << coords[1] << std::endl; 
  //   LinearInterp_impl(coords, my_vals.values(), my_mesh); 
  //   // std::cout << coords[0] << ", " << coords[1] << ":" << LinearInterp_impl(coords, my_vals.values(), my_mesh) << std::endl; 
  // }

  // DiscretizationXD + ICs
  // LinOps::DiscretizationXD my_vals; 
  // auto f = [](double x, double y, double z){ return std::sqrt(x*x + y*y + z*z); }; 
  // my_vals.set_init(my_mesh,f); 

  // printing out flat index -> sub dim index 
  // iterate through flat index 
  // for(std::size_t flat_i=0; flat_i < my_mesh->sizes_product(); flat_i++)
  // {
  //   std::cout << flat_i << ": ";  

  //   // iterate through dims 0,1,2 
  //   for(std::size_t ith_dim=0; ith_dim < my_mesh->dims(); ith_dim++)
  //   {
  //     std::size_t s1 = (flat_i % my_mesh->sizes_middle_product(0,ith_dim+1)); 
  //     std::size_t s2 = s1 / my_mesh->sizes_middle_product(0,ith_dim); 
  //     std::cout << s2 << ", "; 
  //   }
  //   std::cout << "\n"; 
  // }