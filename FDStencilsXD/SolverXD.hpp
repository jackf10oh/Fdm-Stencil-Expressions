// SolverXD.hpp
//
//
//
// JAF 1/11/2026 

#ifndef SOLVERXD_H
#define SOLVERXD_H

#include<cstdint>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/IterativeLinearSolvers> // BICGSTAB

#include "FdmPluginXD.hpp"
#include "BoundaryCondXD.hpp"
#include "../LinOpsXD/MeshXD.hpp"
#include "../LinOpsXD/DiscretizationXD.hpp"
#include "../LinOpsXD/OperatorsXD/IOpXD.hpp"
#include "../Utilities/FillStencil.hpp"

namespace Fds{

struct SolverArgsXD{
  // Member Data -------------------
  std::shared_ptr<const LinOps::MeshXD> domain_mesh_ptr; // XD ! 
  std::shared_ptr<const LinOps::Mesh1D> time_mesh_ptr; // 1D ! 
  BcXDPtr_t bcs; 
  LinOps::DiscretizationXD ICs; // XD ! 
  bool time_dep_flag=true; 
}; 

template<typename EXPR_T, template<typename MAT_T> class EIGENSOLVER_T=Eigen::BiCGSTAB>
class SolverXD{

  private:
    // member data -------------------------------------------
    EXPR_T& m_expr; 
    EIGENSOLVER_T<MatrixStorage_t> m_solver; 

  public:
    // Constructors / Destructors =========================================
    SolverXD()=delete;
    SolverXD(EXPR_T& expr_init)
      : m_expr(expr_init), m_solver()
    {};
    SolverXD(const SolverXD& other)=delete; 
    ~SolverXD()=default;

    // Member Functions =====================================
    // solve PDE Explicitly -----------------------------------
    LinOps::DiscretizationXD Calculate(const SolverArgsXD& args){

      // initialize Solution to start at ICs 
      LinOps::DiscretizationXD solution = args.ICs; 

      // resize stencils  
      m_expr.set_mesh(args.domain_mesh_ptr); 

      // iterate t = t0, t1, ... , tn = T 
      auto it = args.time_mesh_ptr->cbegin()+1; 
      auto end = args.time_mesh_ptr->cend(); 
      double t;
      double t_prev = args.time_mesh_ptr->cbegin()[0];  
      switch (args.time_dep_flag)
      {
      // PDE is time dependent 
      case true:
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // update RHS's time 
          m_expr.SetTime(t_prev);

          // update BCs time 
          // args.bcs_pair.first->SetTime(t);
          args.bcs->SetTime(t);

          // explicit step;
          Eigen::VectorXd v = solution.values() + (t-t_prev) * m_expr.GetMat() * solution.values(); 
          solution = std::move(v); 
          args.bcs->SetSol(solution.values(), args.domain_mesh_ptr); 


          // update into old_time 
          t_prev = t; 
        }
        break;
      // PDE is not time dependent 
      default: // false 
        MatrixStorage_t stencil = m_expr.GetMat(); 
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // explicit step;
          Eigen::VectorXd v = solution.values() + (t-t_prev) * stencil * solution.values(); 
          solution = std::move(v); 

          // update into old_time 
          t_prev = t; 
        }
        break;
      }
      // store mesh into Discretization1D solution
      solution.set_mesh(args.domain_mesh_ptr); 
      // solution fully calculated. 
      return solution; 
    }

    // solve PDE Implicitly -----------------------------------
    LinOps::DiscretizationXD CalculateImp(const SolverArgsXD& args, std::size_t max_iters=20){

      // initialize Solution to start at ICs 
      LinOps::DiscretizationXD solution = args.ICs; 

      // resize stencils  
      m_expr.set_mesh(args.domain_mesh_ptr); 

      // iterate t = t0, t1, ... , tn = T 
      auto it = args.time_mesh_ptr->cbegin()+1; 
      auto end = args.time_mesh_ptr->cend(); 
      double t;
      double t_prev = args.time_mesh_ptr->cbegin()[0];  
      MatrixStorage_t stencil; 
      MatrixStorage_t I = LinOps::IOpXD(args.domain_mesh_ptr).GetMat();  
      m_solver.setMaxIterations(max_iters); 
      switch (args.time_dep_flag)
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
          args.bcs->SetTime(t);

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          // args.bcs_pair.first->SetStencilL(stencil, args.domain_mesh_ptr);
          args.bcs->SetStencil(stencil, args.domain_mesh_ptr);

          // set solutions L/R side 
          // args.bcs_pair.first->SetImpSolL(solution.values(), args.domain_mesh_ptr);
          args.bcs->SetImpSol(solution.values(), args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          m_solver.compute(stencil);
          Eigen::VectorXd v = m_solver.solveWithGuess(solution.values(), solution.values()); 
          solution = std::move(v); 

          // update into old_time 
          t_prev = t; 
        }
        break;

      // PDE is not time dependent 
      default: // false  
        MatrixStorage_t bcs_mask(args.domain_mesh_ptr->sizes_product(), args.domain_mesh_ptr->sizes_product()); 
        args.bcs->SetStencil(bcs_mask, args.domain_mesh_ptr);
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          overwrite_stencil(stencil, bcs_mask);

          // set solutions L/R side 
          args.bcs->SetImpSol(solution.values(), args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          m_solver.compute(stencil);
          Eigen::VectorXd v = m_solver.solveWithGuess(solution.values(), solution.values()); 
          solution = std::move(v); 

          // update into old_time 
          t_prev = t; 
        }
        break;
      }
      // store mesh into Discretization1D solution
      solution.set_mesh(args.domain_mesh_ptr); 
      // solution fully calculated. 
      return solution; 
    }
};

} // end namespace Fds 

#endif // SolverXD.hpp    