// GeneralSolver.hpp
//
// class to solve any PDE of the form LHS_EXPR_IN_TIME_DERIVS = RHS_EXPR_IN_SPACE_DERIVS
//
// JAF 1/17/2026 

#ifndef GENERALSOLVER_H
#define GENERALSOLVER_H 

template<typename ANYMESH_SHAREDPTR_T, typename ANYBC_SHAREDPTR_T, typename CONTAINER_T>
struct GenSolverArgs
{
  // shared_ptr to const Mesh1D or const MeshXD 
  // std::variant<std::shared_ptr<const Mesh1D>, std::shared_ptr<const MeshXD>> domain_mesh_ptr; 
  ANYMESH_SHAREDPTR_T domain_mesh_ptr; 

  // shared_ptr to const Mesh1D 
  std::shared_ptr<const Mesh1D> time_mesh_ptr; 

  // shared_ptr to BcPtr_t or BCXDPtr_t
  ANYBC_SHAREDPTR_T bcs;
  
  // List of underlying Eigen::VectorXd from Discretization1D (XD)  
  CONTAINER_T ICs;

  // boolean flag. 
  bool time_dep_flag=true; 
}; 

template<typename LHS_EXPR, typename RHS_EXPR>
class GenSolver
{
  // Member Data -------------------------------------
  LHS_EXPR& m_lhs; // expression of time derivatives 
  RHS_EXPR& m_rhs; // expression of spatial derivatives 

  // Constructors + Destructor ===========================
  GenSolver()=delete; 
  GenSolver(LHS_EXPR& l_init, RHS_EXPR& r_init)
    : m_lhs(l_init), m_rhs(r_init)
  {}
  GenSolver(const GenSolver& other)=delete;
  // destructor 
  ~GenSolver()=default;  
}; 

#endif

