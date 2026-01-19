// GeneralSolver.hpp
//
// class to solve any PDE of the form LHS_EXPR_IN_TIME_DERIVS = RHS_EXPR_IN_SPACE_DERIVS
//
// JAF 1/17/2026 

#ifndef GENERALSOLVER_H
#define GENERALSOLVER_H 

#include "TExprExecutor.hpp" // TExprExecutor 
#include "../LinOps/Mesh.hpp" // Mesh1D
#include "../LinOpsXD/MeshXD.hpp" // MeshXD 
#include "../LinOps/Discretization.hpp" // Discretization1D ( unnecessary? )
#include "../LinOpsXD/DiscretizationXD.hpp" // DiscretizationXD ( unnecessary? ) 
#include "../LinOps/Operators/IOp.hpp" // IOp
#include "../LinOpsXD/OperatorsXD/IOpXD.hpp" // IOpXD
#include "../Utilities/FillStencil.hpp" // overwrite_stencil 

namespace TExprs{

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

template<typename LHS_EXPR, typename RHS_EXPR, template<typename MAT_T> class EIGENSOLVER_T=Eigen::BiCGSTAB>
class GenSolver
{
  private:
    // Member Data -------------------------------------
    LHS_EXPR& m_lhs; // expression of time derivatives 
    RHS_EXPR& m_rhs; // expression of spatial derivatives 
    EIGENSOLVER_T<TExprs::internal::MatrixStorage_t> m_solver; // Eigen sparse iterative solver

  public:
    // Constructors + Destructor ===========================
    GenSolver()=delete; 
    GenSolver(LHS_EXPR& l_init, RHS_EXPR& r_init)
      : m_lhs(l_init), m_rhs(r_init), m_solver()
    {}
    GenSolver(const GenSolver& other)=delete;
    // destructor 
    ~GenSolver()=default;  

    // Member Funcs =================================================
    // Takes GenSolverArgs and find solution U(t=T) with explicit steps 
    template<typename M, typename B, typename C>
    auto Calculate(const GenSolverArgs<M,B,C>& args)
    {

      TExprs::internal::TExprExecutor exec(m_lhs); 

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
          m_rhs.SetTime(*std::prev(it));
          exec.SetTime(*it); 
          args.bcs->SetTime(*it); 

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
      // pack last solution into Discretization1D (XD) and store mesh 
      using DISCRETIZATION_T = std::conditional_t<
        std::is_same<std::shared_ptr<const LinOps::Mesh1D>, M>::value,
        LinOps::Discretization1D, 
        LinOps::DiscretizationXD
      >; 
      DISCRETIZATION_T result; 
      result = std::move(exec.MostRecentSol()); 
      result.resize(args.domain_mesh_ptr); 
      return result; 
    }

    // Takes GenSolverArgs and find solution U(t=T) with explicit steps 
    template<typename M, typename B, typename C>
    auto CalculateImp(const GenSolverArgs<M,B,C>& args, std::size_t max_iters = 10)
    {
      m_solver.setMaxIterations(max_iters); 

      TExprs::internal::TExprExecutor exec(m_lhs); 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = args.time_mesh_ptr->cbegin() + args.ICs.size(); 
      auto end = args.time_mesh_ptr->cend(); 

      exec.ConsumeSolutionList(args.ICs.cbegin(), args.ICs.cend()); 
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), it); 

      TExprs::internal::MatrixStorage_t I; 
      if constexpr(std::is_same<std::shared_ptr<const LinOps::Mesh1D>, M>::value){
        I = LinOps::IOp(args.domain_mesh_ptr).GetMat();  
      }
      if constexpr(std::is_same<std::shared_ptr<const LinOps::MeshXD>, M>::value){
        I = LinOps::IOpXD(args.domain_mesh_ptr).GetMat();  
      }

      switch (args.time_dep_flag)
      {
      case true:
        for(; it!= end; it++)
        {
          m_rhs.SetTime(*it);
          exec.SetTime(*it); 
          args.bcs->SetTime(*it); 

          Eigen::VectorXd rhs = exec.BuildRhs(*it); 

          // Apply BCs 
          args.bcs->SetImpSol(rhs, args.domain_mesh_ptr); 

          // build implicit stencil 
          TExprs::internal::MatrixStorage_t A = I - exec.inv_coeff_util()*m_rhs.GetMat(); 
          args.bcs->SetStencil(A, args.domain_mesh_ptr); 

          // Implicit Step 
          m_solver.compute(A); 
          Eigen::VectorXd next_sol = m_solver.solveWithGuess(rhs,rhs); 

          // push Solution, time to executor 
          exec.ConsumeSolution(next_sol);
          exec.ConsumeTime(*it); 
        }
        break;
      
      default: // false 
        TExprs::internal::MatrixStorage_t bcs_mask(args.domain_mesh_ptr->size(), args.domain_mesh_ptr->size()); 
        // args.bcs_pair.first->SetStencilL(bcs_mask, args.domain_mesh_ptr);
        args.bcs->SetStencil(bcs_mask, args.domain_mesh_ptr);
        for(; it!= end; it++)
        {
          Eigen::VectorXd rhs = exec.BuildRhs(*it); 

          // Apply BCs 
          args.bcs->SetImpSol(rhs, args.domain_mesh_ptr); 

          // build implicit stencil 
          TExprs::internal::MatrixStorage_t A = I - exec.inv_coeff_util()*m_rhs.GetMat(); 
          overwrite_stencil(A, bcs_mask); // applies BCs. since the rows don't change through time.

          // Implicit Step 
          m_solver.compute(A); 
          Eigen::VectorXd next_sol = m_solver.solveWithGuess(rhs,rhs); 

          // push Solution, time to executor 
          exec.ConsumeSolution(next_sol);
          exec.ConsumeTime(*it); 
        }
        break;
      } // end switch(args.time_dep_flag)
      // pack last solution into Discretization1D (XD) and store mesh 
      using DISCRETIZATION_T = std::conditional_t<
        std::is_same<std::shared_ptr<const LinOps::Mesh1D>, M>::value,
        LinOps::Discretization1D, 
        LinOps::DiscretizationXD
      >; 
      DISCRETIZATION_T result; 
      result = std::move(exec.MostRecentSol()); 
      result.resize(args.domain_mesh_ptr); 
      return result; 
    }

}; 

} // end namespace TExprs 

#endif

