// experimental_NthDerivOp.hpp
//
// uses the same #ifndef macro as NthDerivOp.hpp 
// this will overwrite any occurences of NthDerivOp class
//
// JAF 2/4/2026 

#ifndef DIFFOP_NTHDERIVOP_H
#define DIFFOP_NTHDERIVOP_H

#include "DiffOpBase.hpp"

namespace DiffOps{

template<std::size_t ORDER>
class NthDerivOp : public DiffOpBase<NthDerivOp<ORDER>, ORDER> 
{
  public:
    // Constructors + Destructor ===================================================
    //default 
    NthDerivOp()=default; 
    // from mesh 
    NthDerivOp(const LinOps::Mesh1D_SPtr_t& m)
    {set_mesh(m);};
    // copy 
    NthDerivOp(const NthDerivOp& other)=default; 

    // destructor
    ~NthDerivOp()=default; 

    // Member Funcs ---------------------------------
    template<typename Cont>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return weights[ORDER * n_nodes_per_row + ith_node]; 
    } 
    using DiffOpBase<NthDerivOp, ORDER>::get_weak_mesh1d; 
    using DiffOpBase<NthDerivOp, ORDER>::get_mesh1d; 
    using DiffOpBase<NthDerivOp, ORDER>::GetMat; 
    using DiffOpBase<NthDerivOp, ORDER>::set_mesh; 
}; 

} // end namespace LinOps

#endif // experimental_NthDerivOp.hpp