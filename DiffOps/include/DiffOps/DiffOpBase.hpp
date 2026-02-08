// DiffOpBase.hpp
//
// CRTP base that differential operator inherit from. overrides exression templating in LinOpBase for other DiffOps
//
// JAF 2/4/2026 

#ifndef DIFFOPBASE_H
#define DIFFOPBASE_H

#include <LinOps/LinearOpBase.hpp>
#include "DiffOpTraits.hpp"
#include "DiffOpExpr.hpp"

namespace DiffOps{

using LinOps::MatrixStorage_t; 

template<typename DERIVED, std::size_t ORDER>
class DiffOpBase : public LinOps::LinOpMixIn<DiffOpBase<DERIVED,ORDER>>, public LinOps::LinOpBase1D<DiffOpBase<DERIVED,ORDER>>
{
  public:
    // Type Defs ------------------------- 
    using DERIVED_T = DERIVED; // LinOpMixin / Base can still access grand children 
    using is_diffop_tag = void; 
    constexpr static std::size_t ORDER_N = ORDER; 
  
  private:
    // Member Data ------------------------------------ 
    LinOps::Mesh1D_WPtr_t m_mesh_ptr; 
    MatrixStorage_t m_Mat; 

  public:
    // Member Funcs =====================================================
    // Derived classes must implement ... 
    template<typename Cont>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node) const 
    {
      return static_cast<const DERIVED*>(this)->CoeffAt(weights, n_nodes_per_row, ith_node); 
    } 

    // Defaults Implemented ----------------------------------------------------- 

    // Matrix Getters 
    MatrixStorage_t& GetMat(){ return m_Mat; }; 
    
    const MatrixStorage_t& GetMat() const { return m_Mat; };  
    
    // mesh getters 
    LinOps::Mesh1D_WPtr_t get_weak_mesh1d() const {return m_mesh_ptr; }

    LinOps::Mesh1D_SPtr_t get_mesh1d() const { return m_mesh_ptr.lock(); }

    // set the mesh domain the derivative operator works on 
    void set_mesh(const LinOps::Mesh1D_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)){
        /* if any leaf of the expression is time dependent this needs to be disabled...*/
        return;
      }
      // throw error on nullptr 
      if(!m) throw std::runtime_error("NthDerivOp .set_mesh(Mesh1D_SPtr_t) error: Mesh1D_SPtr_t is expired!"); 
      // point to new mesh 
      m_mesh_ptr = m; 

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
          std::size_t row_start = i*(1+one_sided_skirt); 
          outers_data[i] = row_start; 
          for(int offset=0; offset<1+one_sided_skirt; ++offset){
            inners_data[row_start+offset] = i+offset;
            vals_data[row_start+offset] = CoeffAt(weight_calc.m_arr, 1 + one_sided_skirt, offset); 
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

    
    // Operators ===================================================
    // L1 + L2 (lval) ---------------------------------------------- 
    template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    auto operator+(DIFFOP_T&& rhs) & {
        std::cout << "operator+() & called by DiffOpBase. " << std::endl; 
        constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
        return DiffOpExpr<typename DERIVED::DERIVED_T&, DIFFOP_T, internal::diffop_bin_add_op, MAX>(
        static_cast<typename DERIVED::DERIVED_T&>(*this),  // lhs
        std::forward<DIFFOP_T>(rhs) // rhs 
      );
    }
    
    // L1 + L2 (rval)  
    template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    auto operator+(DIFFOP_T&& rhs) && {
        std::cout << "operator+() && called by DiffOpBase. " << std::endl; 
        constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
        return DiffOpExpr<typename DERIVED::DERIVED_T&&, DIFFOP_T, internal::diffop_bin_add_op, MAX>(
        static_cast<typename DERIVED::DERIVED_T&&>(*this), // lhs 
        std::forward<DIFFOP_T>(rhs) // rhs
      );
    }

    // L1 - L2 (lval) ---------------------------------------------- 
    template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    auto operator-(DIFFOP_T&& rhs) & {
        std::cout << "operator-() & called by DiffOpBase. " << std::endl; 
        constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
        return DiffOpExpr<typename DERIVED::DERIVED_T&, DIFFOP_T, internal::diffop_bin_subtract_op, MAX>(
        static_cast<typename DERIVED::DERIVED_T&>(*this),  // lhs
        std::forward<DIFFOP_T>(rhs) // rhs 
      );
    }
    
    // L1 - L2 (rval)  
    template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    auto operator-(DIFFOP_T&& rhs) && {
        std::cout << "operator-() && called by DiffOpBase. " << std::endl; 
        constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
        return DiffOpExpr<typename DERIVED::DERIVED_T&&, DIFFOP_T, internal::diffop_bin_subtract_op, MAX>(
        static_cast<typename DERIVED::DERIVED_T&&>(*this), // lhs 
        std::forward<DIFFOP_T>(rhs) // rhs
      );
    }

    // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    template<typename SCALAR_T>
    auto left_scalar_mult_impl(SCALAR_T&& c) & {
      std::cout << "left_scalar_mult(c) & called by DiffOpBase." << std::endl; 
      return DiffOpExpr<SCALAR_T, typename DERIVED::DERIVED_T&, internal::diffop_left_mult_op, ORDER>(
        std::forward<SCALAR_T>(c), // lhs scalar
        static_cast<typename DERIVED::DERIVED_T&>(*this) // rhs 
      );
    }
    
    // Left multiply by a scalar: i.e. c*L (rval)
    template<typename SCALAR_T>
    auto left_scalar_mult_impl(SCALAR_T&& c) && {
      std::cout << "left_scalar_mult(c) && called by DiffOpBase." << std::endl; 
      return DiffOpExpr<SCALAR_T, typename DERIVED::DERIVED_T&&, internal::diffop_left_mult_op, ORDER>(
        std::forward<SCALAR_T>(c), // lhs scalar
        static_cast<typename DERIVED::DERIVED_T&>(*this) // rhs 
      );
    }

    // Unary negation (Lval) ---------------------------------------

    // Unary negation (rval) 
        
}; 

} // end namespace DiffOps 

#endif
