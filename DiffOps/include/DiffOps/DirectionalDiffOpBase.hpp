// DiffOpBase.hpp
//
// CRTP base that directional differential operator inherit from
//
// JAF 2/4/2026 

#ifndef DIFFOPBASE_H
#define DIFFOPBASE_H

#include <LinOps/LinearOpBase.hpp>
#include<Utilities/HighDimExpr.hpp> 
#include<Utilities/BlockDiagExpr.hpp> 
// #include "DiffOpTraits.hpp"
// #include "DiffOpExpr.hpp"

namespace DiffOps{

using LinOps::MatrixStorage_t; 

template<typename DERIVED, std::size_t ORDER_, std::size_t DIRECTION_>
class DirectionalDiffOpBase : public LinOps::LinOpBaseXD<DirectionalDiffOpBase<DERIVED,ORDER_,DIRECTION_>>
{
  public:
    // Type Defs ------------------------- 
    using DERIVED_T = DERIVED; // LinOpMixin / Base can still access grand children 
    using is_directional_diffop_tag = void; 
    constexpr static std::size_t ORDER = ORDER_; 
    constexpr static std::size_t DIRECTION = DIRECTION_; 
  
  private:
    // Member Data ------------------------------------ 
    LinOps::MeshXD_WPtr_t m_mesh_ptr; 
    MatrixStorage_t m_Mat; 
    std::size_t m_prod_before; 
    std::size_t m_prod_after; 

  public:
    // Member Funcs =====================================================
    // Derived classes must implement ... 
    template<typename Cont, std::size_t N = 0>
    double CoeffAt(const Cont& weights, std::size_t n_nodes_per_row, std::size_t ith_node, const std::array<double, N>& coords = {}) const 
    {
      return static_cast<const DERIVED*>(this)->CoeffAt(weights, n_nodes_per_row, ith_node, coords); 
    } 

    // Defaults Implemented ----------------------------------------------------- 

    // Matrix Getters 
    auto GetMat(){ return  make_BlockDiag(make_HighDim(m_Mat,m_prod_before),m_prod_after); }; 
    auto GetMat() const { return make_BlockDiag(make_HighDim(m_Mat,m_prod_before),m_prod_after); };  
    
    // mesh getters 
    LinOps::MeshXD_WPtr_t get_weak_meshxd() const {return m_mesh_ptr; }
    LinOps::MeshXD_SPtr_t get_meshxd() const { return m_mesh_ptr.lock(); }

    // set the mesh domain the derivative operator works on 
    void set_mesh(const LinOps::MeshXD_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)){
        /* if any leaf of the expression is time dependent this needs to be disabled...*/
        return;
      }
      // throw error on nullptr 
      if(!m) throw std::runtime_error("DirectionalDiffOpBase .set_mesh(MeshXD_SPtr_t) error: MeshXD_SPtr_t is expired!"); 
      // point to new mesh 
      m_mesh_ptr = m; 

      // calculate new matrix entries 
      /* if there's no coefficient operators in expression, just store the 1D matrix and return kronecker formula...*/
      this->m_prod_before = 1; 
      this->m_prod_after = m->sizes_middle_product(DIRECTION+1, m->dims()); 
      const std::size_t temp_prod_before = m->sizes_middle_product(0,DIRECTION); 

      const std::size_t mesh_size = m->dim_size(DIRECTION);

      // resize Matrix
      m_Mat.resize(temp_prod_before*mesh_size,temp_prod_before*mesh_size); 

      // skirts are # of nodes to side of x_bar that are used to calculate finite difference weights 
      constexpr int one_sided_skirt = ORDER;  
      constexpr int one_sided_size = 1 + one_sided_skirt; 
      constexpr int centered_skirt = (ORDER+1)/2;  
      constexpr int centered_size = 1 + 2*centered_skirt; 

      // number of non zeros in the sparse matrix 
      std::size_t nnz = 2*centered_skirt*(one_sided_size) + (mesh_size-2*centered_skirt)*(centered_size);
      // resize NNZ data in Matrix 
      m_Mat.resizeNonZeros(temp_prod_before*nnz);

      auto* outers_data = m_Mat.outerIndexPtr(); 
      auto* vals_data = m_Mat.valuePtr();
      auto* inners_data = m_Mat.innerIndexPtr();

      // begin OpenMP parallel section. 
      #pragma omp parallel 
      {
        // instantiate stateful fornberg calculator. one per thread 
        FornCalc weight_calc(centered_size,ORDER);
        /* each thread gets a std::array<double, N> for coords ... */

        // first rows with forward stencil
        #pragma omp for nowait 
        for(std::size_t i=0; i < temp_prod_before*centered_skirt; ++i)
        {
          auto left = m->GetMesh(DIRECTION)->cbegin() + (i/temp_prod_before);

          auto right = left + one_sided_size; 

          weight_calc.Calculate(*left, left,right, ORDER); 

          std::size_t row_start = i*(one_sided_size); 
          
          outers_data[i] = row_start; 
          for(int offset=0; offset < one_sided_size; ++offset){
            inners_data[row_start+offset] = (i%temp_prod_before) + (offset*temp_prod_before);
            vals_data[row_start+offset] = CoeffAt(weight_calc.m_arr, one_sided_size, offset/*, coords*/); 
          }
        }
        // middle rows with centered centered stencil  
        #pragma omp for nowait 
        for(std::size_t i=temp_prod_before*centered_skirt; i < (temp_prod_before*mesh_size)-(temp_prod_before*centered_skirt); ++i)
        {
          auto left = m->GetMesh(DIRECTION)->cbegin() - centered_skirt + (i/temp_prod_before);  
          auto right = left + centered_size;
          weight_calc.Calculate(left[centered_skirt], left,right, ORDER);

          int row_start = (temp_prod_before*centered_skirt)*(one_sided_size)
                        + (i-temp_prod_before*centered_skirt)*(centered_size);

          outers_data[i] = row_start; 
          int offset = -(centered_skirt);
          for(int j=0; j<centered_size; ++j){
            inners_data[row_start + j] = temp_prod_before*((i-temp_prod_before*centered_skirt) / temp_prod_before) + (i%temp_prod_before) + (centered_skirt+offset)*temp_prod_before;
            vals_data[row_start + j] = CoeffAt(weight_calc.m_arr, centered_size, j/*,  coords*/);
            offset += 1; 
          };
        }
        // last rows 
        #pragma omp for nowait 
        for(std::size_t i = (temp_prod_before*mesh_size)-(temp_prod_before*centered_skirt); i < (temp_prod_before*mesh_size); ++i)
        {
          auto right = m->GetMesh(DIRECTION)->cbegin() + (i/temp_prod_before) + 1; 
          auto left = std::prev(right,one_sided_size); 
          weight_calc.Calculate(*std::prev(right), left,right, ORDER); 

          int row_start = (temp_prod_before*centered_skirt)*(one_sided_size)
                        + (temp_prod_before*mesh_size - 2*temp_prod_before*centered_skirt)*(centered_size) 
                        + (i - ((temp_prod_before*mesh_size)-(temp_prod_before*centered_skirt)))*(one_sided_size);
          
          outers_data[i] = row_start;
          int offset= -(one_sided_skirt);
          for(int j=0; j < one_sided_size; ++j){
            inners_data[row_start + j] = temp_prod_before*(i/temp_prod_before) + (i%temp_prod_before) + (offset*temp_prod_before); //
            vals_data[row_start + j] = CoeffAt(weight_calc.m_arr, one_sided_size, j /*, coords...*/);  
            offset += 1;
          }
        }
        // #pragma omp single
        // {
        //   std::cout << "Threads in this parallel region: " << omp_get_num_threads() << "\n";
        // }
      } 
      // End OpenMP parallel section. implicit barrier 
      outers_data[temp_prod_before*mesh_size] = temp_prod_before*nnz; 
    } 

    
    // // Operators ===================================================
    // // L1 + L2 (lval) ---------------------------------------------- 
    // template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    // auto operator+(DIFFOP_T&& rhs) & {
    //     std::cout << "operator+() & called by DiffOpBase. " << std::endl; 
    //     constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
    //     return DiffOpExpr<typename DERIVED::DERIVED_T&, DIFFOP_T, internal::diffop_bin_add_op, MAX>(
    //     static_cast<typename DERIVED::DERIVED_T&>(*this),  // lhs
    //     std::forward<DIFFOP_T>(rhs) // rhs 
    //   );
    // }
    
    // // L1 + L2 (rval)  
    // template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    // auto operator+(DIFFOP_T&& rhs) && {
    //     std::cout << "operator+() && called by DiffOpBase. " << std::endl; 
    //     constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
    //     return DiffOpExpr<typename DERIVED::DERIVED_T&&, DIFFOP_T, internal::diffop_bin_add_op, MAX>(
    //     static_cast<typename DERIVED::DERIVED_T&&>(*this), // lhs 
    //     std::forward<DIFFOP_T>(rhs) // rhs
    //   );
    // }

    // // L1 - L2 (lval) ---------------------------------------------- 
    // template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    // auto operator-(DIFFOP_T&& rhs) & {
    //     std::cout << "operator-() & called by DiffOpBase. " << std::endl; 
    //     constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
    //     return DiffOpExpr<typename DERIVED::DERIVED_T&, DIFFOP_T, internal::diffop_bin_subtract_op, MAX>(
    //     static_cast<typename DERIVED::DERIVED_T&>(*this),  // lhs
    //     std::forward<DIFFOP_T>(rhs) // rhs 
    //   );
    // }
    
    // // L1 - L2 (rval)  
    // template<typename DIFFOP_T, typename = std::enable_if_t<traits::is_diffop_crtp<DIFFOP_T>::value>>
    // auto operator-(DIFFOP_T&& rhs) && {
    //     std::cout << "operator-() && called by DiffOpBase. " << std::endl; 
    //     constexpr std::size_t MAX = std::max(ORDER, DIFFOP_T::ORDER_N); // constexpr since C++14 
    //     return DiffOpExpr<typename DERIVED::DERIVED_T&&, DIFFOP_T, internal::diffop_bin_subtract_op, MAX>(
    //     static_cast<typename DERIVED::DERIVED_T&&>(*this), // lhs 
    //     std::forward<DIFFOP_T>(rhs) // rhs
    //   );
    // }

    // // Left multiply by a scalar: i.e. c*L (lval)------------------------------------------------------------------------- 
    // template<typename SCALAR_T>
    // auto left_scalar_mult_impl(SCALAR_T&& c) & {
    //   std::cout << "left_scalar_mult(c) & called by DiffOpBase." << std::endl; 
    //   return DiffOpExpr<SCALAR_T, typename DERIVED::DERIVED_T&, internal::diffop_left_mult_op, ORDER>(
    //     std::forward<SCALAR_T>(c), // lhs scalar
    //     static_cast<typename DERIVED::DERIVED_T&>(*this) // rhs 
    //   );
    // }
    
    // // Left multiply by a scalar: i.e. c*L (rval)
    // template<typename SCALAR_T>
    // auto left_scalar_mult_impl(SCALAR_T&& c) && {
    //   std::cout << "left_scalar_mult(c) && called by DiffOpBase." << std::endl; 
    //   return DiffOpExpr<SCALAR_T, typename DERIVED::DERIVED_T&&, internal::diffop_left_mult_op, ORDER>(
    //     std::forward<SCALAR_T>(c), // lhs scalar
    //     static_cast<typename DERIVED::DERIVED_T&>(*this) // rhs 
    //   );
    // }

    // // Unary negation (Lval) ---------------------------------------

    // // Unary negation (rval) 
        
}; 

} // end namespace DiffOps 

#endif
