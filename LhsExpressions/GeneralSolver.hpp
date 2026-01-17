// GeneralSolver.hpp
//
// class to solve any PDE of the form LHS_EXPR_IN_TIME_DERIVS = RHS_EXPR_IN_SPACE_DERIVS
//
// JAF 1/17/2026 

#ifndef GENERALSOLVER_H
#define GENERALSOLVER_H 

#include "../LinOps/Mesh.hpp"

template<typename ANYMESH_SHAREDPTR_T, typename ANYBC_SHAREDPTR_T, typename CONTAINER_T>
struct GenSolverArgs
{
  // shared_ptr to const Mesh1D or const MeshXD 
  // std::variant<std::shared_ptr<const Mesh1D>, std::shared_ptr<const MeshXD>> domain_mesh_ptr; 
  ANYMESH_SHAREDPTR_T domain_mesh_ptr; 

  // shared_ptr to const Mesh1D 
  std::shared_ptr<const LinOps::Mesh1D> time_mesh_ptr; 

  // shared_ptr to BcPtr_t or BCXDPtr_t
  ANYBC_SHAREDPTR_T bcs;
  
  // List of underlying Eigen::VectorXd from Discretization1D (XD)  
  CONTAINER_T ICs;

  // boolean flag. 
  bool time_dep_flag=true; 
}; 

// CTAD Guide 
template<typename M, typename B, typename C>
GenSolverArgs(M, std::shared_ptr<const LinOps::Mesh1D>, B, C, bool)
-> GenSolverArgs<M, B, C>;

template<typename LHS_EXPR, typename RHS_EXPR>
class GenSolver
{
  private:
    // Member Data -------------------------------------
    LHS_EXPR& m_lhs; // expression of time derivatives 
    RHS_EXPR& m_rhs; // expression of spatial derivatives 

  public:
    // Constructors + Destructor ===========================
    GenSolver()=delete; 
    GenSolver(LHS_EXPR& l_init, RHS_EXPR& r_init)
      : m_lhs(l_init), m_rhs(r_init)
    {}
    GenSolver(const GenSolver& other)=delete;
    // destructor 
    ~GenSolver()=default;  

    // Member Funcs =================================================
    template<typename... Ts>
    auto Calculate(const GenSolverArgs<Ts...>& args)
    {

      LhsExecutor exec(m_lhs); 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = args.time_mesh_ptr->cbegin() + args.ICs.size(); 
      auto end = args.time_mesh_ptr->cend(); 

      exec.ConsumeSolutionList(args.ICs.cbegin(), args.ICs.cend()); 
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), it); 

      switch (args.time_dep_flag)
      {
      case true:
        for(; it!= end; it++)
        {
          m_rhs.SetTime(*it);
          exec.SetTime(*it); 

          Eigen::VectorXd rhs = exec.BuildRhs(*it); 
          
          // Explicit Step 
          Eigen::VectorXd next_sol = exec.inv_coeff_util()*m_rhs.GetMat()*exec.MostRecentSol() + rhs; 
          
          // Apply BCs 
          args.bcs->SetSol(next_sol, args.domain_mesh_ptr); 

          // push Solution, time to executor 
          exec.ConsumeSolution(next_sol);
          exec.ConsumeTime(*it); 
        }
        break;
      
      default: // false 
        for(; it!= end; it++)
        {
          Eigen::VectorXd rhs = exec.BuildRhs(*it); 
          
          // Explicit Step 
          Eigen::VectorXd next_sol = exec.inv_coeff_util()*m_rhs.GetMat()*exec.MostRecentSol() + rhs; 
          
          // Apply BCs 
          args.bcs->SetSol(next_sol, args.domain_mesh_ptr); 

          // push Solution, time to executor 
          exec.ConsumeSolution(next_sol);
          exec.ConsumeTime(*it); 
        }
        break;
      } // end switch(args.time_dep_flag)

      // give last solution as result
      return std::move(exec.MostRecentSol()); 
    }
}; 

#endif

