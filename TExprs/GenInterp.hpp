// GenInterp.hpp
//
//
//
// JAF 1/19/2026 

#ifndef GENINTERP_H
#define GENINTERP_H 

#include<memory>
#include<vector>
#include<utility> // std::pair
#include<Eigen/Core>

#include "GenSolver.hpp"

namespace TExprs{

template<typename LHS_EXPR, typename RHS_EXPR, typename M, typename B, typename C>
class GenInterp
{
  public:
    // Type Defs ----------------------------- 
    struct VecSaveWrite
    {
      std::shared_ptr<std::vector<Eigen::VectorXd>> m_ptr; 
      VecSaveWrite()=delete; 
      VecSaveWrite(std::shared_ptr<std::vector<Eigen::VectorXd>> p_init) : m_ptr(p_init){}; 
      VecSaveWrite(const VecSaveWrite& other) : m_ptr(other.m_ptr){}; 
      void SaveSolution(Eigen::VectorXd&& sol){ m_ptr->emplace_back(sol); }; 
      void ConsumeLastSolution(Eigen::VectorXd&& sol){ m_ptr->emplace_back(sol); }; 
    }; 

    // Member Data ------------------------------
    std::shared_ptr<std::vector<Eigen::VectorXd>> m_data; 
    VecSaveWrite m_sol_writer = VecSaveWrite(m_data); 
    GenSolver<LHS_EXPR, RHS_EXPR, VecSaveWrite> m_solver; 
    GenSolverArgs<M,B,C> m_args; 
    bool m_calculated; 

  public: 
    // Constructors + Destructor ===================================
    GenInterp()=delete; 
    GenInterp(LHS_EXPR& lhs, RHS_EXPR& rhs, const GenSolverArgs<M,B,C>& args_init = GenSolverArgs<M,B,C>{})
      : m_data(std::make_shared<std::vector<Eigen::VectorXd>>(0)), 
      m_sol_writer(m_data), 
      m_args(args_init), 
      m_solver(lhs,rhs,m_sol_writer), 
      m_calculated(false) 
    {
      if(m_args.time_mesh_ptr) m_data->reserve(m_args.time_mesh_ptr->size());     
    }
    GenInterp(const GenInterp& other)=delete; 
    // destructor 
    ~GenInterp()=default; 

    // Member Funcs ====================================================
    // get value of Solution at any point t,x1,x2,...xn in time/space 
    double at(double t, const std::vector<double>& coords){
      // if m_data is empty... 
      if(!m_calculated) FillVals(); 

      // find index in m_args.time_mesh_ptr
      auto time_interval_pair = get_interval(*m_args.time_mesh_ptr, t);

      // find left / right value in linear interpolation 
      double val_01 =  LinearInterp_impl(coords, (*m_data)[std::distance(m_args.time_mesh_ptr->cbegin(), time_interval_pair.first)]); 
      double val_02 =  LinearInterp_impl(coords, (*m_data)[std::distance(m_args.time_mesh_ptr->cbegin(), time_interval_pair.second)]);
      
      // linear interpolation (t-t1) * (y2 - y1) / (t2 - t1) 
      return val_01 + (t - *time_interval_pair.first) * (val_02 - val_01) / (*time_interval_pair.second - *time_interval_pair.first); 
    }

    // set m_args to a new input
    void SetArgs(const GenSolverArgs<M,B,C>& args_switch)
    {
      m_args = args_switch; 
      m_data->resize(0); 
      m_data->reserve(m_args.time_mesh_ptr->size()); 
      m_calculated=false; 
    }

    // Populate m_data with solutions at each step in time 
    void FillVals(std::size_t num_iters=20)
    {
      if(!m_calculated)
      {
        std::for_each(m_args.ICs.cbegin(), m_args.ICs.cend(), [&](const auto& ic){m_data->push_back(ic);});  
        // WritePolicy moves all solutions at each time step to m_data
        m_solver.CalculateImp(m_args, num_iters); 
        m_calculated = true; 
      }
    }

  public:
    // Unreachable ----------------------------------------------------
    double LinearInterp_impl(const std::vector<double>& coords, const Eigen::VectorXd& v)
    {
      std::shared_ptr<const LinOps::MeshXD> high_dim_m = nullptr; 

      if constexpr(std::is_same<M, std::shared_ptr<const LinOps::MeshXD>>::value) high_dim_m = m_args.domain_mesh_ptr; 
      else high_dim_m = std::make_shared<const LinOps::MeshXD>( std::vector<std::shared_ptr<const LinOps::Mesh1D>>(1, m_args.domain_mesh_ptr) ); 

      if(coords.size() != high_dim_m->dims()) throw std::runtime_error(" GenInterp::LinearInterp(...) error. coords.size() != GenSolverArgs.domain_mesh_ptr->dims(). or != in case of Mesh1D"); 

      // Current Mesh1D we are working on 
      const auto sub_dim_m = high_dim_m->GetMeshAt(high_dim_m->dims()-1); 

      // find pair [x1, x2] with x1 <= x <= x2 (x = coords[ith_dim])
      auto interval_pair = get_interval(*sub_dim_m, coords[coords.size()-1]); 

      std::pair<double,double> sub_dim_vals = LinearInterp_recursive_impl(coords, v, high_dim_m, coords.size()-1); 

      return sub_dim_vals.first + (coords[coords.size()-1] - *interval_pair.first) * (sub_dim_vals.second - sub_dim_vals.first) / (*interval_pair.second - *interval_pair.first);  
    }

    std::pair<double,double> LinearInterp_recursive_impl(
      const std::vector<double>& coords, 
      const Eigen::VectorXd& v, 
      const std::shared_ptr<const LinOps::MeshXD>& m,
      std::size_t ith_dim,
      std::size_t cumulative_offset = 0
    )
    {
      // Current Mesh1D we are working on 
      const auto sub_dim_m = m->GetMeshAt(ith_dim); 

      // find pair [x1, x2] with x1 <= x <= x2 (x = coords[ith_dim])
      auto internval_pair = get_interval(*sub_dim_m, coords[ith_dim]); 

      // if this is first dimension 
      if(ith_dim == 0)
      {
        // ! requires offsets of all higher dimension. 
        // ! actually return disc(idx(...) + some_offset_from_higher_dims)
        // return discretization(cumulative_offset+idx(p.first)), discretization(cumulative_offset+idx(p.second)) 
        double val1 = v[cumulative_offset + std::distance(sub_dim_m->cbegin(), internval_pair.first)]; 
        double val2 = v[cumulative_offset + std::distance(sub_dim_m->cbegin(), internval_pair.second)]; 
        return std::pair(val1,val2); 
      }
      else if(ith_dim == (coords.size()-1)) // this is the outtermost dimension. 
      {
        // cumulative offset needs to be manually set to nonzero before it can accumulate products 
        cumulative_offset = m->sizes_middle_product(0, ith_dim) * std::distance(sub_dim_m->cbegin(), internval_pair.first); 
        
        // gets 2 linear interpretations from lower dimension 
        std::pair<double,double> val_pair = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, cumulative_offset);
        return val_pair; 
      }
      else // this is an intermediate dimension
      {
        // calculate new offset by incorporating this 
        std::size_t cumulative_offset_01 = cumulative_offset + m->sizes_middle_product(0, ith_dim) * std::distance(sub_dim_m->cbegin(), internval_pair.first); 
        std::size_t cumulative_offset_02 = cumulative_offset + m->sizes_middle_product(0, ith_dim) * std::distance(sub_dim_m->cbegin(), internval_pair.second); 

        // get 2 pairs of linear interpolant along subdimension (ith_dim-1)
        std::pair<double,double> val_pair_01 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, cumulative_offset_01);
        std::pair<double,double> val_pair_02 = LinearInterp_recursive_impl(coords, v, m, ith_dim-1, cumulative_offset_02);

        // perform linear interpolation in this dimension to produce pair of values 
        double val_01 = (coords[ith_dim] - *internval_pair.first) * (val_pair_01.second-val_pair_01.first) / (*internval_pair.second-*internval_pair.first);
        double val_02 = (coords[ith_dim] - *internval_pair.first) * (val_pair_02.second-val_pair_02.first) / (*internval_pair.second-*internval_pair.first);

        // return std::pair<double,double> 
        return {val_01, val_02}; 
      }
    }; 

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
}; 

} // end namespace TExprs 

#endif // GenInterp.hpp


// Linear interpolation scratch work 


// 1D Case: 

/* 
Find pair [x1,x2] with x1 <= x <= x2] and return (x-x1) * (val[idx(x2)] - val[idx(x1)]) / (x2-x1)
*/

// 2D Case:

/* 
find pair [y1,y2] with y1 <= y <= y2.

get idx(y1) and idx(y2)

along OneDim_view( idx(y1) ) 
get [x1,x2] with x1 <= x <= x2 and return ... 

along OneDim_view( idx(y2) )
get [x1,x2] ......

return (y-y1) * (1D_Interp(idx(y2)) - 1D_Interp(idx(y1))) / (y2-y1)
*/

// 3D Case:

/*
find pair [z1,z2] with z1 <= z <= z2 

get idx(z1) and idx(z2) 

along TwoDim_view( idx(z1) ) 
get [y1,y2] -> get [x1,x2] (twice) ... return 

along TwoDim_view( idx(z2) ) 
get [y1,y2] -> get [x1,x2] (twice) ... return 

return (z-z1) * (2D_Interp(idx(z2)) - 2D_INterp(idx(z1))) / (z2 - z1)

*/

// In General ND Case:

/*
find pair [x1, x2] with x1 <= x <= x2 


( 
  necessary data for N-1 Interpolation: 
    std::size_t N-1 : next inner dimension to get new [x1,x2] pairs in
    std::size_t k : offset to get solution value at in higher dimension   
)
along N-1 Dim View ( idx(x1) ) get val 1 
along N-1 Dim View ( idx(x2) ) get val 2

return (x - x1) * (val2 - val1) / (x2 - x1)
*/


/*
find pair [x1, x2] with x1 <= x <= x2 
auto p = get_interval(*m->GetMeshAt( current_ith_dimension ) )

if this is first dimension 
{
  ! requires offsets of all higher dimension. 
  ! actually return disc(idx(...) + some_offset_from_higher_dims)
  return discretization(idx(p.first)), discretization(idx(p.second)) 
}

if this is an intermediate dimension
{
  auto pair1 = LinearInterp_recursive_impl(...) with idx(p.first), (current dimension-1)
  auto pair2 = LinearInterp_recursive_impl(...) with idx(p.second), (current dimension-1)

  return (middle of pair1), (middle of pair2) 
}

( 
  necessary data for N-1 Interpolation: 
    std::size_t N-1 : next inner dimension to get new [x1,x2] pairs in
    std::size_t k : offset to get solution value at in higher dimension   
)
*/

