// NthPartialDeriv.hpp
//
//
//
// JAF 2/10/2026 

#ifndef NTHPARTIALDERIV_H
#define NTHPARTIALDERIV_H 

#include "DirectionalDiffOpBase.hpp"

namespace DiffOps{

template<std::size_t ORDER_, std::size_t DIRECTION_>
class NthPartialDeriv : public DirectionalDiffOpBase<NthPartialDeriv<ORDER_,DIRECTION_>, ORDER_, DIRECTION_> 
{
  public:
    // Constructors + Destructor ===================================================
    //default 
    NthPartialDeriv()=default; 
    // from mesh 
    NthPartialDeriv(const LinOps::MeshXD_SPtr_t& m)
    {this->set_mesh(m);};
    // copy 
    NthPartialDeriv(const NthPartialDeriv& other)=default; 

    // destructor
    ~NthPartialDeriv()=default; 

    // Member Funcs ---------------------------------
    template<typename Cont, std::size_t N = 0>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node, const std::array<double, N> = {}) const 
    {
      return weights[ORDER_ * n_nodes_per_row + ith_node]; 
    }
}; 

} // end namespace DiffOps


#endif // NthPartialDeriv.hpp 