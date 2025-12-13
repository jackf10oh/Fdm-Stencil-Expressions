// FdmPlugin.hpp
//
//
//
// JAF 

#ifndef BVPPLUGIN_H
#define BVPPLUGIN_H 

#define LINOP_PLUGIN FdmPlugin
#define CUSTOM_IDENTITY_MATRIX_STORAGE Eigen::SparseMatrix<double,Eigen::ColMajor>

#include<iostream>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/LU>
#include<eigen3/Eigen/SparseCore>
#include<eigen3/Eigen/SparseLU>
#include "../LinOps/Discretization.hpp"
#include "../LinOps/LinOpTraits.hpp"
#include "BoundaryCond.hpp"

// also declared in BoundaryCond.hpp 
using MatrixStorage_t = Eigen::SparseMatrix<double, Eigen::ColMajor>; 

template<typename BaseDerived>
class FdmPlugin
{
  public:
    // member data 
    BcPtr_t lbc_ptr; // ptr to left bc
    BcPtr_t rbc_ptr; // ptr to right bc
  protected:
    double m_current_time; // current time
    MatrixStorage_t m_stencil; // stencil matrix of FDM scheme // don't need to store? 

  public:
    // constructors ---------------------------------------------------------
    FdmPlugin()
      : lbc_ptr(std::make_shared<BoundaryCond>()), 
      rbc_ptr(std::make_shared<BoundaryCond>()), 
      m_current_time(0.0)
    {}
    // destructors ---------------------------------------------------------
    ~FdmPlugin()=default; 

    // Member Funcs -------------------------------------------------------
    // set the time of L/R BCs to t 
    void SetTime_impl(double t){}; 
    void SetTime(double t)
    {
      // if we are an expression 
      if constexpr (is_expr_crtp<typename BaseDerived::Derived_t>::value){
        // std::cout << "Branch hit: "; 
        // conditionally call SetTime() on L/R nodes
        auto& expr = static_cast<typename BaseDerived::Derived_t&>(*this);
        if constexpr(is_linop_crtp<typename BaseDerived::Derived_t::LStorage_t>::value) 
        {
          // std::cout << "L "; 
          expr.Lhs().SetTime(t);
        }
        if constexpr(is_linop_crtp<typename BaseDerived::Derived_t::RStorage_t>::value)
        {
          // std::cout << "R"; 
          expr.Rhs().SetTime(t);
        } 
        // std::cout << std::endl; 
      }
      // set current time, update time in boundary conditions as well 
      m_current_time = t; 
      lbc_ptr->SetTime(t); 
      rbc_ptr->SetTime(t); 

      // we aren't an expression call the actual implementor 
      static_cast<typename BaseDerived::Derived_t*>(this)->SetTime_impl(t);
    }
    // get the time this operator is set to
    auto Time() const {return m_current_time;}; 
    
    // Get Underlying Matrix of Derived     
    decltype(auto) GetMat()
    {
      return static_cast<typename BaseDerived::Derived_t*>(this)->GetMat(); 
    };
    decltype(auto) GetMat() const
    {
      return static_cast<const typename BaseDerived::Derived_t*>(this)->GetMat(); 
    };

    // Use .apply() from LinOpBase and set BCs afterward 
    Discretization1D explicit_step(const Discretization1D& d) const
    {
      // temp store of result 
      Discretization1D result = static_cast<const typename BaseDerived::Derived_t*>(this)->apply(d);

      // fix the boundary conditions 
      lbc_ptr->SetSolL(result);
      lbc_ptr->SetSolR(result);

      return result; 
    };

    Discretization1D apply(const Discretization1D& d) const
    {
      // temp store of result 
      Discretization1D result(d.mesh()); 
      result = GetMat() * d.values();
      return result; 
    };

    // New functionaly. useful for implicit schemes 
    Discretization1D solve_implicit(const Discretization1D& d) 
    {
      // applies m_stencil=m_stencil for most derived classes. 
      m_stencil = GetMat(); 
      lbc_ptr->SetStencilL(m_stencil, d.mesh());
      rbc_ptr->SetStencilR(m_stencil, d.mesh());  
      // copy construct a new Discretization
      Discretization1D imp_sol = d; 
      lbc_ptr->SetImpSolL(imp_sol);
      lbc_ptr->SetImpSolR(imp_sol);
      
      // now we have the expression A*x = b
      // where x is the solution at timestep n+1 
      // and the rhs b is the given discretization 

      Discretization1D result(d.mesh());
      Eigen::SparseLU<MatrixStorage_t, Eigen::COLAMDOrdering<int> > solver(m_stencil);
      result = solver.solve(imp_sol.values());
      return result; 
    }
    // Compose stencil(linop) with another stencil(linop) 
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      // typical result of composition 
      auto result = static_cast<BaseDerived*>(this)->compose(std::forward<DerivedInner>(InnerOp));
      // if I have a left BC give it to result 
      if(lbc_ptr) result.lbc_ptr = lbc_ptr;
      else result.lbc_ptr = InnerOp.lbc_ptr;
      // if I have a right BC give it to result 
      if(rbc_ptr) result.rbc_ptr = rbc_ptr; 
      else result.rbc_ptr = InnerOp.rbc_ptr;
      // give result 
      return result;
    };

    // test function
    void print(){std::cout << "Hello world, from plugin!" << std::endl; }
};

#endif // BvpPlugin.hpp