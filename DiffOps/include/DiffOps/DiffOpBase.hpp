// DiffOpBase.hpp
//
// CRTP base that differential operator inherit from. overrides exression templating in LinOpBase for other DiffOps
//
// JAF 2/4/2026 

#ifndef DIFFOPBASE_H
#define DIFFOPBASE_H

#include <LinOps/LinearOpBase.hpp>
#include "DiffOpTraits.hpp"

namespace DiffOps{

using LinOps::MatrixStorage_t; 

namespace internal{

// leafs of an expression hold weak_ptr + Sparse Matrix
template<typename BASE_T, typename = void>
struct DiffOpBaseData
{
  LinOps::Mesh1D_WPtr_t m_mesh_ptr; 
  MatrixStorage_t m_Mat; 
}; 

// template have no member data ... 
template<typename EXPR_T>
struct DiffOpBaseData<EXPR_T, std::enable_if_t<DiffOps::traits::is_diffop_expr_crtp<EXPR_T>::value>>
{/* empty struct */}; 

} // end namespace DiffOps::internal 

template<typename DERIVED, std::size_t ORDER>
class DiffOpBase : public LinOps::LinOpMixIn<DiffOpBase<DERIVED,ORDER>>, public LinOps::LinOpBase1D<DiffOpBase<DERIVED,ORDER>>, private internal::DiffOpBaseData<DiffOpBase<DERIVED,ORDER>>
{
  public:
    // Type Defs ------------------------- 
    using DERIVED_T = DERIVED; // LinOpMixin / Base can still access grand children 
    using is_diffop_tag = void; 

    // Member Funcs =====================================================

    // Derived classes must implement ... 
    template<typename Cont>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return static_cast<const DERIVED*>(this)->CoeffAt(weights, n_nodes_per_row, ith_node); 
    } 

    // Defaults Implemented ----------------------------------------------------- 
    // current maximum order or Derivative expression  
    constexpr std::size_t Order() const {return ORDER; };

    // Matrix Getters 
    MatrixStorage_t& GetMat()
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_diffop_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<DERIVED&>(*this);  
        if constexpr(traits::is_diffop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_mesh1d(); 
        else if constexpr(traits::is_diffop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_mesh1d(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_Mat;
      } 
    }; 
    const MatrixStorage_t& GetMat() const 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_diffop_expr_crtp<DERIVED>::value){
        const auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_diffop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_mesh1d(); 
        else if constexpr(traits::is_diffop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_mesh1d(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_Mat;
      } 
    };  
    
    // mesh getters 
    LinOps::Mesh1D_WPtr_t get_weak_mesh1d() const
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_diffop_expr_crtp<DERIVED>::value){
        const auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_diffop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_mesh1d(); 
        else if constexpr(traits::is_diffop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_mesh1d(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_mesh_ptr;
      } 
    }

    // return Mesh1D pointed to
    LinOps::Mesh1D_SPtr_t get_mesh1d() const 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_diffop_expr_crtp<DERIVED>::value){
        const auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_diffop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_mesh1d(); 
        else if constexpr(traits::is_diffop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_mesh1d(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return this->m_mesh_ptr.lock();
      }
    } 

    void assign_to_weak_ptr(const LinOps::Mesh1D_SPtr_t& m)
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_diffop_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<DERIVED&>(*this);  
        if constexpr(traits::is_diffop_crtp<typename DERIVED::LStorage_t>::value){
          expr.Lhs().assign_to_weak_ptr(m); 
          return; 
        }
        else if constexpr(traits::is_diffop_crtp<typename DERIVED::RStorage_t>::value){
          expr.Rhs().assign_to_weak_ptr(m); 
          return; 
        }
        else static_assert(false, "cannot call assign_to_weak_ptr() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        this->m_mesh_ptr = m;
        return; 
      } 
    }; 

    // set the mesh domain the derivative operator works on 
    void set_mesh(const LinOps::Mesh1D_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      auto my_m = this->get_weak_mesh1d(); 
      if(!my_m.owner_before(m) && !m.owner_before(my_m)) return;
      // throw error on nullptr 
      if(!m) throw std::runtime_error("NthDerivOp .set_mesh(Mesh1D_SPtr_t) error: Mesh1D_SPtr_t is expired!"); 
      // point to new mesh 
      assign_to_weak_ptr(m); 

      // calculate new matrix entries 
      const std::size_t mesh_size = m->size();
      // resize Matrix
      GetMat().resize(mesh_size,mesh_size); 

      // skirts are # of nodes to side of x_bar that are used to calculate finite difference weights 
      constexpr int one_sided_skirt = ORDER;  
      constexpr int centered_skirt = (ORDER+1)/2;  

      // number of non zeros in the sparse matrix 
      std::size_t nnz = 2*centered_skirt*(1+one_sided_skirt) + (mesh_size-2*centered_skirt)*(1+2*centered_skirt);
      // resize NNZ data in Matrix 
      GetMat().resizeNonZeros(nnz);

      auto* outers_data = GetMat().outerIndexPtr(); 
      auto* vals_data = GetMat().valuePtr();
      auto* inners_data = GetMat().innerIndexPtr();

      // begin OpenMP parallel section. 
      #pragma omp parallel 
      {
        // instantiate stateful fornberg calculator. one per thread 
        FornCalc weight_calc(1+2*centered_skirt,ORDER);
        // first rows with forward stencil
        #pragma omp for nowait 
        for(std::size_t i=0; i<centered_skirt; i++)
        {
          auto left = m->cbegin()+i;
          auto right = left+one_sided_skirt+1; 
          weight_calc.Calculate(*left, left,right, ORDER); 
          std::size_t offset=i; 
          outers_data[i] = i*(1+one_sided_skirt); 
          for(int j=0; j<1+one_sided_skirt; ++j){
            inners_data[i*(1+one_sided_skirt)
                                      + offset] = offset;
            vals_data[i*(1+one_sided_skirt)
                                      + offset] = CoeffAt(weight_calc.m_arr, 1 + one_sided_skirt, j); 
            offset += 1; 
          }
        }
        // middle rows with centered centered stencil  
        #pragma omp for nowait 
        for(std::size_t i=centered_skirt;i<mesh_size-centered_skirt; i++)
        {
          auto left = m->cbegin()-centered_skirt+i;  
          auto right = left+(centered_skirt+1+centered_skirt);
          weight_calc.Calculate(m->at(i), left,right, ORDER);
          int offset = -centered_skirt;
          outers_data[i] = centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset); 
          for(int j=0; j< 1+2*centered_skirt; ++j){
            inners_data[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = i+offset;
            vals_data[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = CoeffAt(weight_calc.m_arr, 1+2*centered_skirt, j);
            offset += 1; 
          };
        }
        // last rows 
        #pragma omp for nowait 
        for(std::size_t i = mesh_size-centered_skirt; i<mesh_size; i++)
        {
          auto right = m->cbegin() + i; 
          auto left = right-(one_sided_skirt+1); 
          weight_calc.Calculate(*right, left,right, ORDER); 
          int offset= -one_sided_skirt;
          outers_data[i] = centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset); 
          for(int j=0; j < 1+one_sided_skirt; ++j){
            inners_data[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = i+offset;  
            vals_data[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = CoeffAt(weight_calc.m_arr, 1+one_sided_skirt, j);  
            offset += 1;
          }
        }
        // #pragma omp single
        // {
        //   std::cout << "Threads in this parallel region: " << omp_get_num_threads() << "\n";
        // }
      } 
      // End OpenMP parallel section. implicit barrier 
      outers_data[mesh_size] = nnz; 
    } 
}; 

} // end namespace DiffOps 

#endif
