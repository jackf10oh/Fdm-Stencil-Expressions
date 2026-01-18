// NthTimeDeriv.hpp
//
//
//
// JAF 1/15/2026 

#ifndef NTHTIMEDERIV_H
#define NTHTIMEDERIV_H 

#include "TimeDerivBase.hpp"

// ===============================================================
class NthTimeDeriv : public TimeDerivBase<NthTimeDeriv> 
{
  private:
    // Member Data ---------------

  public:
    // Constructors + Destructor ====================================
    // NthTimeDeriv()=delete; // necessary?  
    NthTimeDeriv(std::size_t order=1)
      : TimeDerivBase<NthTimeDeriv>(order)
    {}
    NthTimeDeriv(const NthTimeDeriv& other)
      : TimeDerivBase<NthTimeDeriv>(other.m_order)
    {} 
    // destructor 
    ~NthTimeDeriv()=default; 

    // Member Funcs =================================================== 
    template<typename Cont>
    decltype(auto) CoeffAt(const Cont& v, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      std::size_t offset = this->m_order * n_nodes_per_row; 
      return v[ith_node + offset]; 
    } 

    // using LhsBase<NthTimeDeriv>::toTuple; 
    std::string toString() const {return "hi from NthTimeDeriv"; }; 

}; 

#endif // NthTimeDeriv.hpp