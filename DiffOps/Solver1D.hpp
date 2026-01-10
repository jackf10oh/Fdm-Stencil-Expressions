// Solver.hpp
//
// 
//
// JAF 1/8/2025 

#ifndef SOLVER_H
#define SOLVER_H

#include<cstdint>
#include<iostream>
#include<Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/LU>
#include<Eigen/SparseCore>
#include<Eigen/SparseLU> 
#include<Eigen/IterativeLinearSolvers> // BICGSTAB

#include "FdmPlugin.hpp"
#include "BoundaryCond.hpp"
#include "../LinOps/Mesh.hpp"
#include "../LinOps/Discretization.hpp"
#include "../LinOps/Operators/IOp.hpp"
#include "../Utilities/FillStencil.hpp"

// Type Defs -------- 
using MatrixStorage_t = Eigen::SparseMatrix<double, Eigen::RowMajor>; 

struct SolverArgs1D{
  // Member Data -------------------
  std::shared_ptr<const Mesh1D> domain_mesh_ptr; 
  std::shared_ptr<const Mesh1D> time_mesh_ptr; 
  std::pair<BcPtr_t, BcPtr_t> bcs_pair; 
  Discretization1D ICs;
  bool time_dep_flag=true; 
}; 

template<typename EXPR_T>
class Solver1D{

  private:
    // member data -------------------------------------------
    EXPR_T& m_expr; 

  public:
    // Constructors / Destructors =========================================
    Solver1D()=delete;
    Solver1D(EXPR_T& expr_init)
      : m_expr(expr_init)
    {};
    Solver1D(const Solver1D& other)=delete; 
    ~Solver1D()=default;

    // Member Functions =====================================
    // solve PDE Explicitly -----------------------------------
    Discretization1D Calculate(const SolverArgs1D& args){

      // initialize Solution to start at ICs 
      Discretization1D solution = args.ICs; 

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
          args.bcs_pair.first->SetTime(t);
          args.bcs_pair.second->SetTime(t);

          // explicit step;
          Eigen::VectorXd v = solution.values() + (t-t_prev) * m_expr.GetMat() * solution.values(); 
          args.bcs_pair.first->SetSolL(v, args.domain_mesh_ptr); 
          args.bcs_pair.second->SetSolR(v, args.domain_mesh_ptr); 
          solution = std::move(v); 

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
    Discretization1D CalculateImp(const SolverArgs1D& args, std::size_t max_iters=20){

      // initialize Solution to start at ICs 
      Discretization1D solution = args.ICs; 

      // resize stencils  
      m_expr.set_mesh(args.domain_mesh_ptr); 

      // iterate t = t0, t1, ... , tn = T 
      auto it = args.time_mesh_ptr->cbegin()+1; 
      auto end = args.time_mesh_ptr->cend(); 
      double t;
      double t_prev = args.time_mesh_ptr->cbegin()[0];  
      MatrixStorage_t stencil; 
      auto I = IOp(args.domain_mesh_ptr).GetMat();  
      Eigen::BiCGSTAB<MatrixStorage_t> mat_solver(stencil);
      mat_solver.setMaxIterations(max_iters); 
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
          args.bcs_pair.first->SetTime(t);
          args.bcs_pair.second->SetTime(t);

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          args.bcs_pair.first->SetStencilL(stencil, args.domain_mesh_ptr);
          args.bcs_pair.first->SetStencilR(stencil, args.domain_mesh_ptr);

          // set solutions L/R side 
          args.bcs_pair.first->SetImpSolL(solution.values(), args.domain_mesh_ptr);
          args.bcs_pair.second->SetImpSolR(solution.values(), args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          mat_solver.compute(stencil);
          Eigen::VectorXd v = mat_solver.solveWithGuess(solution.values(), solution.values()); 
          solution = std::move(v); 

          // update into old_time 
          t_prev = t; 
        }
        break;

      // PDE is not time dependent 
      default: // false  
        MatrixStorage_t bcs_mask(args.domain_mesh_ptr->size(), args.domain_mesh_ptr->size()); 
        args.bcs_pair.first->SetStencilL(bcs_mask, args.domain_mesh_ptr);
        args.bcs_pair.first->SetStencilR(bcs_mask, args.domain_mesh_ptr);
        for(; it!=end; it++){
          // get new time, and dt 
          t = *it;

          // calculate matrix A 
          stencil = I - (t-t_prev) * m_expr.GetMat(); 

          // set first / last row of A 
          overwrite_stencil(stencil, bcs_mask);

          // set solutions L/R side 
          args.bcs_pair.first->SetImpSolL(solution.values(), args.domain_mesh_ptr);
          args.bcs_pair.second->SetImpSolR(solution.values(), args.domain_mesh_ptr);

          // implicit step. modifies solution in place 
          mat_solver.compute(stencil);
          Eigen::VectorXd v = mat_solver.solveWithGuess(solution.values(), solution.values()); 
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

#endif // Solver1D.hpp