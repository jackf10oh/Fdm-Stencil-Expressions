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

template<typename Cont>
auto get_interval(const Cont& v, double c)
{
// runtime checks 
if(v.size() < 2) throw std::runtime_error("size of v < 2"); 
if(c < v.cbegin()[0]) throw std::runtime_error("c < v[0]"); 

// right side a(i+1) 
auto after = std::lower_bound(v.cbegin(), v.cend(), c);
if(after == v.cend()) throw std::runtime_error("right bound == v.cend()"); 

// if a(i+1) == v[0] bump it by 1. 
auto before = (after==v.cbegin()) ? after++ : std::prev(after); 
return std::pair(before, after); 
}; 

double LinearInterp_recursive_impl(
  const std::vector<double>& coords, 
  const Eigen::VectorXd& v, 
  const std::shared_ptr<const LinOps::MeshXD>& m,
  std::size_t ith_dim,
  std::size_t cumulative_offset = 0
)
{
  auto sub_dim_m = m->GetMeshAt(ith_dim); 
  auto bounding_interval = get_interval(*sub_dim_m, coords[ith_dim]);  

  if(ith_dim == 0)
  {
    std::size_t final_offset_01 = cumulative_offset + std::distance(sub_dim_m->cbegin(), bounding_interval.first); 
    std::size_t final_offset_02 = final_offset_01 + 1; 
    // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
    double result = v[final_offset_01] + (coords[ith_dim] - *bounding_interval.first) * (v[final_offset_02]-v[final_offset_01]) / (*bounding_interval.second - *bounding_interval.first);  
    return result; 
  }
  else
  {
    std::size_t stride_size = m->sizes_middle_product(0,ith_dim); 
    std::size_t interval_start_idx = std::distance(sub_dim_m->cbegin(), bounding_interval.first); 
    std::size_t offset_01 = cumulative_offset + stride_size * (interval_start_idx); 
    std::size_t offset_02 = cumulative_offset + stride_size * (interval_start_idx+1); 

    double endpoint_01 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, offset_01); 
    double endpoint_02 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, offset_02);
    
    // result = y1 + (c-x1) * (y2-y1) / (x2-x1)
    double result = endpoint_01 + (coords[ith_dim] - *bounding_interval.first) * (endpoint_02 - endpoint_01) / (*bounding_interval.second - *bounding_interval.first); 
    return result; 
  }
}; 

double LinearInterp(const std::vector<double>& coords, const Eigen::VectorXd& v, const std::shared_ptr<const LinOps::MeshXD>& m)
{
  return LinearInterp_recursive_impl(coords, v, m, coords.size()-1);
}; 

int main()
{
  // iomanip 
  // std::cout << std::setprecision(4); 
  
  // Domain MeshXD ----------------------------------------
  auto my_mesh = LinOps::make_mesh(0.0,4.0,5); 

  // double r = 4.0;
  // std::size_t N = 5;
  // auto my_mesh = LinOps::make_meshXD(0.0,r,N, 2);

  // std::vector<std::pair<double,double>> end_points = {{0.0,1.0}, {0.0,5.0}}; 
  // std::vector<std::size_t> n_steps = {4,9}; 
  // auto my_mesh = LinOps::make_meshXD(end_points, n_steps);

  std::cout << "my_mesh type: " << typeid(decltype(my_mesh)).name() << std::endl; 
  

  std::shared_ptr<const LinOps::MeshXD> copied = make_meshXD(my_mesh); 

  std::cout << "# Dims: " << copied->dims() << std::endl; 
  std::cout << "dim(0) size: " << copied->dim_size(0) << std::endl; 
  std::cout << "sizes prod: " << copied->sizes_product() << std::endl; 

};


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