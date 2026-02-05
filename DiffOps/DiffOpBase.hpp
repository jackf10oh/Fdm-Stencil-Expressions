// DiffOpBase.hpp
//
// CRTP base that differential operator inherit from. overrides exression templating in LinOpBase for other DiffOps
//
// JAF 2/4/2026 

#ifndef DIFFOPBASE_H
#define DIFFOPBASE_H

#include "../LinOps/include/LinOps/LinearOpBase.hpp"

template<typename Derived, std::size_t ORDER>
class DiffOpBase : public LinOps::LinOpMixIn<DiffOpBase>, public LinOps::LinOpBase1D<DiffOpBase>
{
  protected:
    Mesh1D_WPtr_t m_mesh_ptr; 
    LinOps::MatrixStorage_t m_Mat; 

  public:
    // Type Defs ------------------------- 
    using DERIVED_T = Derived; // LinOpMixin / Base can still access grand children 

    // Member Funcs =====================================================

    // Derived classes must implement ... 
    template<typename Cont>
    decltype(auto) CoeffAt(std::size_t mat_row_i, const Mesh1D_SPtr_t& mesh, const Cont& weights, std::size_t n_nodes_per_row, std::size_t jth_node) const 
    {
      return static_cast<const Derived*>(this)->CoeffAt(mat_row_i, mesh, weights, n_nodes_per_row, ith_node); 
    } 

    // Defaults Implemented ----------------------------------------------------- 
    // current maximum order or Derivative expression  
    std::size_t Order() const {return ORDER; };

    // Matrix Getters 
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };  
    
    // mesh getters 
    Mesh1D_WPtr_t get_weak_mesh1d() const { return m_mesh_ptr; }
    Mesh1D_SPtr_t get_mesh1d() const { return m_mesh_ptr.lock(); } 

    // set the mesh domain the derivative operator works on 
    void set_mesh(const Mesh1D_SPtr_t& m){/* need to implement*/}
}; 

#endif
