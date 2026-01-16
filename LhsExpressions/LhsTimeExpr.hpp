// LhsTimeExpr.hpp
//
//
//
// JAF 1/15/2026 

#ifndef LHSTIMEEXPR_H
#define LHSTIMEEXPR_H

#include<Eigen/Core>
#include "LinOps/Discretization.hpp"
#include "Utilities/FornbergCalc.hpp"

// cout << c.GetMat().cwiseInverse() << endl; 
// std::tuple_cat 

template<typename Derived>
struct TimeDerivBase
{
  protected:
    // Member Data -----------------------
    std::size_t m_order; 
    std::size_t m_n_nodes; 
  public:
    TimeDerivBase()=default; 
    TimeDerivBase(std::size_t order_init,std::size_t  n_nodes_init)
      : m_order(order_init), m_n_nodes(n_nodes_init)
    {}
}; 

struct TimeDeriv : public TimeDerivBase<TimeDeriv>
{
  auto CoeffFromWeights(){}
}; 

struct CoeffTimeDeriv : public TimeDerivBase<CoeffTimeDeriv>
{

};

struct TimeDerivSum{
  
}; 

template<std::size_t Order = 1, std::size_t N=Order+1>
struct LhsTimeDeriv
{
  // Member Data ---------------------------------
  std::array<double,N> m_stored_times; 
  std::array<LinOps::Discretization1D,N-1> m_stored_sols; 
  FornCalc m_weights_calc; 

  // Constructors =======================================
  LhsTimeDeriv()
    : m_weights_calc(N,Order)
  {}
  // Member Funcs =======================================
  // from a time. set weights
  void SetWeightsFromTime(double t)
  {
    // set last entry in m_stored_times
    m_stored_times[N-1] = t; 
    // recalculate m_weights_calc 
    m_weights_calc.Calculate(t, m_stored_times.cbegin(), m_stored_times.cend(), Order);     
  }
  
  // consume a time. push back all previous
  void ConsumeTime(double t)
  {
    // iterate from m_stored_times[0] ... [n-2]
    auto end = std::prev(m_stored_times.end(),2); 
    auto it=m_stored_times.begin();  
    for(; it!=end; it++){
      // updated according to m_stored_times[i] = [i+1]
      *it = *std::next(it); 
    }
    *it = t; // m_stored_times[n-1] = t; 
    // all values in m_stored_time have been left shifted by 1
    // first value dropped.
    // second from right most value == t 
    // right most value unassigned 
  } // end ConsumeTime 

  // consume a solution. push back all previous 
  void ConsumeSolution(LinOps::Discretization1D sol)
  { 
    /* same idea as ConsumeTime(). with move semantics*/
    // iterate from 1st to 2nd to last of m_stored_sols
    std::size_t end = N-2; 
    for(std::size_t i=0; i<N-2; i++){
      m_stored_sols[i] = std::move(m_stored_sols[i+1]); 
    }
    // last m_stored_sol get sol input
    m_stored_sols[N-1] = std::move(sol);   
  }

  // // Build RHS starting point from m_stored_sols[0, ..., N-2]
  // Eigen::VectorXd BuildRhs()
  // {
  //   // start with oldest solution.  
  // }

}; 

#endif