// GeneralSolver.hpp
//
// class to solve any PDE of the form LHS_EXPR_IN_TIME_DERIVS = RHS_EXPR_IN_SPACE_DERIVS
//
// JAF 1/17/2026 

#ifndef GENERALSOLVER_H
#define GENERALSOLVER_H 

#include<iostream>
#include<Eigen/Core> // VectorXd
#include<LinOps/Mesh.hpp> // Mesh1D
#include<LinOps/MeshXD.hpp> // MeshXD 
#include<LinOps/Operators/IOp.hpp> // IOp
// #include<Utilities/FillStencil.hpp> // overwrite_stencil 
#include "TExprExecutor.hpp" // TExprExecutor 

namespace TExprs{

// Generalized Args for PDEs in 1D or XD
template<typename ANYMESH_SHAREDPTR_T = LinOps::Mesh1D_SPtr_t>
struct GenSolverArgs
{
  // shared_ptr to const Mesh1D or const MeshXD 
  ANYMESH_SHAREDPTR_T domain_mesh_ptr; 

  // shared_ptr to const Mesh1D 
  LinOps::Mesh1D_SPtr_t time_mesh_ptr; 
  
  // List of underlying Eigen::VectorXd from Discretization1D (XD)  
  std::vector<Eigen::VectorXd> ICs;
}; 

// CTAD Guide 
template<typename M>
GenSolverArgs(M, std::shared_ptr<const LinOps::Mesh1D>, std::vector<Eigen::VectorXd>)
-> GenSolverArgs<M>;

// Write Policies. i.e. write previous solution to cout, write to std::vector, write to CSV 
struct EmptyWrite
{
  // no member data 
  EmptyWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol){}; 
  void ConsumeLastSolution(Eigen::VectorXd&& sol){}; 
};

struct FinalWrite
{
  // no member data 
  FinalWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol){}; 
  auto ConsumeLastSolution(Eigen::VectorXd&& sol){ return sol; }; 
};

struct PrintWrite
{
  bool first_entry=true; 
  PrintWrite()=default; 
  void SaveSolution(Eigen::VectorXd&& sol)
  {
    if(first_entry){
      std::cout << "[";
      first_entry=false;
    }
    std::cout << "["; 
    auto it=sol.cbegin(); 
    auto end = std::prev(sol.cend()); 
    for(; it!=end; it++){
      std::cout << *it << ", ";
    }
    std::cout << *it << "],\n";  
  }; 
  void ConsumeLastSolution(Eigen::VectorXd&& sol)
  {
    std::cout << "["; 
    auto it=sol.cbegin(); 
    auto end = std::prev(sol.cend()); 
    for(; it!=end; it++){
      std::cout << *it << ", ";
    }
    std::cout << *it << "]]\n";  
  }; 
};

// Generalized Solver for PDEs in 1D or XD
template<typename LHS_EXPR, typename RHS_EXPR, typename OSTEP_TUP, typename WRITE_POLICY_T=FinalWrite, template<typename MAT_T> class EIGENSOLVER_T=Eigen::BiCGSTAB>
class GenSolver
{
  private:
    // Member Data -------------------------------------
    LHS_EXPR& m_lhs; // expression of time derivatives 
    RHS_EXPR& m_rhs; // expression of spatial derivatives 
    typename std::remove_reference<OSTEP_TUP>::type m_ostep_tup; 
    WRITE_POLICY_T m_write_p; 

  public:
    // Constructors + Destructor ===========================
    GenSolver()=delete; 
    GenSolver(LHS_EXPR& l_init, RHS_EXPR& r_init, OSTEP_TUP ostep_init, WRITE_POLICY_T wp_init=WRITE_POLICY_T{})
      : m_lhs(l_init), m_rhs(r_init), m_ostep_tup(ostep_init),m_write_p(wp_init)
    {}
    GenSolver(const GenSolver& other)=delete;
    // destructor 
    ~GenSolver()=default;  

    // Member Funcs =================================================
    // Takes GenSolverArgs and find solution U(t=T) with explicit steps 
    template<typename M>
    auto Calculate(GenSolverArgs<M> args)
    {

      TExprs::internal::TExprExecutor exec(m_lhs); 
      MatrixStorage_t Mat; 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = args.time_mesh_ptr->cbegin() + args.ICs.size() - 1; 
      auto end = args.time_mesh_ptr->cend() - 1;

      exec.ConsumeSolutionList(args.ICs.begin(), args.ICs.end()); // no more .cbegin() -> should use move semantics
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), std::next(it)); 

      double t = *it; 
      while(it!= end)
      {
        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ 
            ((lam_args.template BeforeLinAlgebra<decltype(it),decltype(exec),decltype(m_rhs),OSteps::FDStep_Type::IMPLICIT>(it, args.domain_mesh_ptr, exec, m_rhs)), ...); 
          }, 
          m_ostep_tup
        ); 

        Mat = m_rhs.GetMat(); 

        // working on right end of [t(n-1), t(n)]
        ++it;
        t = *it;  

        exec.BuildNextTime(t); 

        // scale Mat according to 1 / dt ... 
        std::visit([&](const auto& inv){ Mat = inv * Mat; }, exec.inv_coeff()); 
        // outside steps matrix before step(Mat) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, Mat)), ...); }, 
          m_ostep_tup
        ); 

        Eigen::VectorXd rhs = std::move(exec.RhsVector()); 
        // outside steps solution before step (rhs) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolBeforeStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, rhs)), ...); }, 
          m_ostep_tup
        ); 

        // Explicit Step 
        Eigen::VectorXd next_sol = Mat * exec.MostRecentSol() + rhs; 
        
        // ourside steps solution after step(next_sol) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolAfterStep<OSteps::FDStep_Type::EXPLICIT>(t, args.domain_mesh_ptr, next_sol)), ...); }, 
          m_ostep_tup
        ); 

        // give expiring solution to WRITE_POLICY_T
        m_write_p.SaveSolution(std::move(exec.ExpiringSol())); 

        // push Solution, time to executor 
        exec.ConsumeSolution(next_sol);
        exec.ConsumeTime(t); 
      }    

      // Write remaining solutions to m_write_p 
      for(auto i=0; i<exec.StoredSols().size()-1; i++) m_write_p.SaveSolution( std::move(exec.StoredSols()[i])); 

      // write policy also determines return type / how to handle last solution
      return m_write_p.ConsumeLastSolution(std::move( exec.MostRecentSol() )); 
    } // end .Calculate(args) 

    // Takes GenSolverArgs and find solution U(t=T) with implicit steps 
    template<typename M>
    auto CalculateImp(GenSolverArgs<M> args, std::size_t max_iters = 20)
    {
      TExprs::internal::TExprExecutor exec(m_lhs); 
      EIGENSOLVER_T<MatrixStorage_t> iterative_solver; // Eigen sparse iterative solver
      iterative_solver.setMaxIterations(max_iters); 
      MatrixStorage_t Mat; 
      // // If I wanted a truly adaptive Domain Mesh, We need to move this elsewhere. i.e. iterate along diagonal....  
      // MatrixStorage_t I = LinOps::IOp(args.domain_mesh_ptr).GetMat(); 

      exec.set_mesh(args.domain_mesh_ptr);
      m_rhs.set_mesh(args.domain_mesh_ptr); 

      auto it = args.time_mesh_ptr->cbegin() + args.ICs.size() - 1; 
      auto end = args.time_mesh_ptr->cend() - 1;

      exec.ConsumeSolutionList(args.ICs.cbegin(), args.ICs.cend()); 
      exec.ConsumeTimeList(args.time_mesh_ptr->cbegin(), std::next(it)); 

      double t = *it; 
      while(it!= end)
      {
        // outside steps before any type of linear algebra is performed... 
        std::apply(
          [&](auto&... lam_args){ ((lam_args.template BeforeLinAlgebra<decltype(it),decltype(exec),decltype(m_rhs),OSteps::FDStep_Type::IMPLICIT>(it, args.domain_mesh_ptr, exec, m_rhs)), ...); }, 
          m_ostep_tup
        ); 

        // working on right end of [t(n-1), t(n)]
        ++it;
        t = *it; 

        // moved into an ostep that sets time + mesh of lhs executor / rhs expression 
        // m_rhs.SetTime(t);
        // exec.SetTime(t); 

        exec.BuildNextTime(t); 
        
        // Mat = I - inv_coeff * FDStencil ; 
        Mat = m_rhs.GetMat(); 
        std::visit([&](const auto& inv){ Mat = -inv * Mat; }, exec.inv_coeff()); 
        for(std::size_t i=0; i<Mat.rows(); ++i) Mat.coeffRef(i,i) += 1.0; 

        // outside steps matrix before step(Mat) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template MatBeforeStep<OSteps::FDStep_Type::IMPLICIT>(t, args.domain_mesh_ptr, Mat)), ...); }, 
          m_ostep_tup
        ); 

        Eigen::VectorXd rhs = std::move(exec.RhsVector()); 
        // outside steps solution before step (rhs) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolBeforeStep<OSteps::FDStep_Type::IMPLICIT>(t, args.domain_mesh_ptr, rhs)), ...); }, 
          m_ostep_tup
        ); 

        // Explicit Step ( U(n+1) = D*U(n) + U(n) )
        // Implicit Step (I - D(t+1))*U(n+1) = rhs 
        iterative_solver.compute(Mat); 
        Eigen::VectorXd next_sol = iterative_solver.solveWithGuess(rhs,rhs); 
        
        // ourside steps solution after step(next_sol) 
        std::apply(
          [&](const auto&... lam_args){ ((lam_args.template SolAfterStep<OSteps::FDStep_Type::IMPLICIT>(t, args.domain_mesh_ptr, next_sol)), ...); }, 
          m_ostep_tup
        ); 

        // give expiring solution to WRITE_POLICY_T
        m_write_p.SaveSolution(std::move(exec.ExpiringSol())); 

        // push Solution, time to executor 
        exec.ConsumeSolution(next_sol);
        exec.ConsumeTime(t); 
      }    

      // Write remaining solutions to m_write_p 
      for(auto i=0; i<exec.StoredSols().size()-1; i++) m_write_p.SaveSolution( std::move(exec.StoredSols()[i])); 

      // write policy also determines return type / how to handle last solution
      return m_write_p.ConsumeLastSolution(std::move( exec.MostRecentSol() )); 
    } // end .CalculateImp(args) 

}; 

} // end namespace TExprs 

#endif

// explicit step 
// [U(t+1) - U(t)] / dt = D * U(t) 
// U(t+1) = U(t) + dt * D * U(t) 


// implicit step 
// [U(t+1) - U(t)] / dt = D * U(t+1) 
// U(t+1) - U(t) = dt * D * U(t+1) 
// I * U(t+1) - dt * D * U(n+1) = U(t) 
// ( I - dt * D) * U(t+1) = U(t) 
// U(t+1) = (I - dt * D).inv(U(t)) 


// With Forcing Terms f(t,u,x,...) 
// explicit step 
// [U(t+1) - U(t)] / dt = D * U(t) + f(t,u,x,...)
// U(t+1) = U(t) + dt * D * U(t) + dt*f(t,u,x,...) 

// IMEX step (implicit step in D(U), explicit in f(t,u,x...))   
// [U(t+1) - U(t)] / dt = D * U(t+1) + f(t,u,x,...)
// U(t+1) - U(t) = dt * D * U(t+1) + dt*f(t,u,x,...)
// I * U(t+1) - dt * D * U(n+1) = U(t) + dt*f(t,u,x,...)
// ( I - dt * D) * U(t+1) = U(t) + dt*f(t,u,x,...)
// U(t+1) = (I - dt * D).inv(   U(t)+dt*f(t,u,x,...)   ) 

// proper signature on forcing term f? 
// in 1D f(double t, const Eigen::VectorXd& U, const Mesh1D_SPtr_t& m)
// in XD f(double t, const Eigen::VectorXd& U, const MeshXD_SPtr_t& m)