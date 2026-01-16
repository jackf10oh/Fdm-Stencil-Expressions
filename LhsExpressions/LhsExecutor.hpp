// LhsExecutor.hpp
//
//
//
// JAF 1/16/2025 

#ifndef LHSEXECUTOR_H
#define LHSEXECUTOR_H 

template<typename LHS_EXPR, typename TUP_T = decltype(std::declval<LHS_EXPR>().toTuple())>
struct LhsExecutor
{
  // Member Data ---------------------------------
  // LHS_EXPR m_expr; 
  TUP_T m_sum_tup; 
  std::size_t m_order; 
  std::size_t m_num_nodes; 
  std::vector<double> m_stored_times; 
  std::vector<LinOps::Discretization1D> m_stored_sols; 
  FornCalc m_weights_calc; 

  // Constructors + Destructor =======================================
  LhsExecutor()=delete; 
  LhsExecutor(LHS_EXPR& expr_init) 
    :m_order(expr_init.Order()), 
    m_num_nodes(expr_init.Order()+1),  
    m_weights_calc(m_order, m_num_nodes), 
    m_stored_times(m_num_nodes), 
    m_stored_sols(m_order), 
    m_sum_tup(expr_init.toTuple()) 
  {}
  LhsExecutor(const LhsExecutor& other)=delete; 
  // destructor 
  ~LhsExecutor()=default; 

  // Member Funcs =======================================
  // from a time. set weights
  void SetWeightsFromTime(double t)
  {
    // set last entry in m_stored_times
    m_stored_times[m_num_nodes-1] = t; 
    // recalculate m_weights_calc 
    m_weights_calc.Calculate(t, m_stored_times.cbegin(), m_stored_times.cend(), m_order);     
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
    std::size_t end = m_num_nodes-2; 
    for(std::size_t i=0; i<end; i++){
      m_stored_sols[i] = std::move(m_stored_sols[i+1]); 
    }
    // last m_stored_sol get sol input
    m_stored_sols[m_num_nodes-2] = std::move(sol);   
  }

  // // Build RHS starting point from m_stored_sols[0, ..., N-2]
  // Eigen::VectorXd BuildRhs()
  // {
  //   // start with oldest solution.  
  // }

}; 

#endif // LhsExecutor.hpp