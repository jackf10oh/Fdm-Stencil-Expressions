// Interp1D.hpp
//
// Class that linearly interpolates between points on a SolverArgs1D's meshes 
//
// JAF 1/13/2026 

#ifndef INTERP1D_H
#define INTERP1D_H 

#include "Solver1D.hpp" 

namespace Fds{
using namespace LinOps; 

template<typename EXPR_T, template<typename MAT_T> class EIGENSOLVER_T=Eigen::BiCGSTAB>
class Interp1D{

  private:
    // Member Data -----------------------------------------------------
    bool m_filled; 
    std::vector<Discretization1D> m_vals; 
    EXPR_T& m_expr; 
    SolverArgs1D m_args; 
    EIGENSOLVER_T<MatrixStorage_t> m_solver; 

  public:
    // Constructors + Destructor =========================================================
    Interp1D()=delete; 
    Interp1D(EXPR_T& expr, SolverArgs1D args_init)
      : m_expr(expr), m_args(args_init), m_filled(false), m_vals(0), m_solver()
    {
      m_vals.reserve(m_args.time_mesh_ptr->size()); 
    }; 
    // Interp1D(EXPR_T& expr)
    //   : m_expr(expr), m_args(), m_filled(false), m_vals(0), m_solver()
    // {}; 
    Interp1D(const Interp1D& other)=delete; 

    // destructor
    ~Interp1D()=default; 

    // Member Functions =====================================================
    void SetArgs(SolverArgs1D args)
    {
      m_args = args; 
      m_vals.resize(0); 
      m_vals.reserve(m_args.time_mesh_ptr->size()); 
    }
    // from m_vals, get linear interpolant. 
    double Value(double t, double x)
    {
      if(!m_filled) FillVals(); 
      
      std::size_t j,k; 

      auto space_p = get_interval(*m_args.domain_mesh_ptr, x);
      k = std::distance(m_args.domain_mesh_ptr->cbegin(), space_p.first);

      auto time_p = get_interval(*m_args.time_mesh_ptr, t);
      j = std::distance(m_args.time_mesh_ptr->cbegin(), time_p.first);
    
      double dx = m_vals[j][k+1] - m_vals[j][k];  
      double x1 = m_vals[j][k] + dx * ((x - *space_p.first) / (*space_p.second - *space_p.first)); 

      dx = m_vals[j+1][k+1] - m_vals[j+1][k];
      double x2 = m_vals[j][k] + dx * ((x - *space_p.first) / (*space_p.second - *space_p.first));  

      double result = x1 + (x2 - x1) * ((t - *time_p.first) / (*time_p.second - *time_p.first)); 
      return result; 
    }; 
  private:
    // Unreachable ----------------------------------------------------------- 
    // solve PDE Implicitly, copy to m_vals -----------------------------------
    void FillVals(std::size_t max_iters=20){

      // resize m_vals to store solution at each time step 
      m_vals.reserve(m_args.time_mesh_ptr->size()); 
      m_vals.resize(0); 

      // initialize Solution to start at ICs 
      Discretization1D solution = m_args.ICs; 
      m_vals.emplace_back(solution); 

      // resize stencils  
      m_expr.set_mesh(m_args.domain_mesh_ptr); 

      // iterate t = t0, t1, ... , tn = T 
      auto it = m_args.time_mesh_ptr->cbegin()+1; 
      auto end = m_args.time_mesh_ptr->cend(); 
      double t;
      double t_prev = m_args.time_mesh_ptr->cbegin()[0];  
      MatrixStorage_t stencil; 
      auto I = IOp(m_args.domain_mesh_ptr).GetMat();  
      m_solver.setMaxIterations(max_iters); 
      switch (m_args.time_dep_flag)
      {
      // PDE is time dependent 
      case true:
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // update RHS's time 
          m_expr.SetTime(t);

          // update BCs time 
          // args.bcs_pair.first->SetTime(t);
          m_args.bcs->SetTime(t);

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          // args.bcs_pair.first->SetStencilL(stencil, args.domain_mesh_ptr);
          m_args.bcs->SetStencil(stencil, m_args.domain_mesh_ptr);

          // set solutions L/R side 
          // args.bcs_pair.first->SetImpSolL(solution.values(), args.domain_mesh_ptr);
          m_args.bcs->SetImpSol(solution.values(), m_args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          m_solver.compute(stencil);
          Eigen::VectorXd v = m_solver.solveWithGuess(solution.values(), solution.values()); 
          solution = std::move(v); 
          m_vals.push_back(solution); 


          // update into old_time 
          t_prev = t; 
        }
        break;

      // PDE is not time dependent 
      default: // false  
        MatrixStorage_t bcs_mask(m_args.domain_mesh_ptr->size(), m_args.domain_mesh_ptr->size()); 
        // args.bcs_pair.first->SetStencilL(bcs_mask, args.domain_mesh_ptr);
        m_args.bcs->SetStencil(bcs_mask, m_args.domain_mesh_ptr);
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          overwrite_stencil(stencil, bcs_mask);

          // set solutions L/R side 
          // args.bcs_pair.first->SetImpSolL(solution.values(), args.domain_mesh_ptr);
          m_args.bcs->SetImpSol(solution.values(), m_args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          m_solver.compute(stencil);
          Eigen::VectorXd v = m_solver.solveWithGuess(solution.values(), solution.values()); 
          solution = std::move(v);
          m_vals.push_back(solution); 

          // update into old_time 
          t_prev = t; 
        }
        break;
      }
      // store mesh into Discretization1D solution
      // solution.set_mesh(m_args.domain_mesh_ptr); 
      m_vals.emplace_back(std::move(solution)); 

      // solution fully calculated. 
      // return solution; 
      m_filled=true; 
    } 

    // From a point x in [a,b] get smallest pair [a(i), a(i+1)] containing x
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

} // end namespace Fds 

#endif // end Interp1D.hpp 