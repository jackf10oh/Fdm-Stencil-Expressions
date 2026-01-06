// NthDerivOp.hpp
//
//
//
// JAF 12/5/2025

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<cstdint>
#include<Eigen/Core>
#include<omp.h>

#include "../../Utilities/FornbergCalc.hpp"
#include "../FdmPlugin.hpp"
#include "../../LinOps/LinearOpBase.hpp"

class NthDerivOp : public LinOpBase<NthDerivOp>
{
  public:
    // member data 
    std::size_t m_order; 
    std::vector<MatrixStorage_t::StorageIndex> m_inners; 
    std::vector<MatrixStorage_t::StorageIndex> m_outers; 
    std::vector<MatrixStorage_t::Scalar> m_vals; 
  public:
    // constructors ---------------------------------
    NthDerivOp(const MeshPtr_t& m, std::size_t order=1)
      : m_order(order)
    {set_mesh(m);};
    NthDerivOp(std::size_t order=1, const MeshPtr_t& m=nullptr)
      : m_order(order)
    {set_mesh(m);};
    
    // Destructors -----------------------------------
    ~NthDerivOp()=default; 
    
    // Member Funcs ---------------------------------
    std::size_t Order() const {return m_order; };
    auto GetMat(){ 
      return Eigen::Map<MatrixStorage_t>(
        m_mesh_ptr->size(), m_mesh_ptr->size(), m_vals.size(), 
        m_outers.data(), m_inners.data(), m_vals.data(), 
        0 // 0 flags the sparse matrix as compressed 
      ); 
    }; 
    auto GetMat() const { 
      return Eigen::Map<const MatrixStorage_t>(
        m_mesh_ptr->size(), m_mesh_ptr->size(), m_vals.size(), 
        m_outers.data(), m_inners.data(), m_vals.data(), 
        0 // 0 flags the sparse matrix as compressed 
      );     
    };     
    // set the mesh the derivative operator works on 
    void set_mesh(MeshPtr_t m)
    {
      //  do nothing if meshes are == or on nullptr 
      if(m==nullptr || m==m_mesh_ptr) return;  
      
      // new pointer set. proceed to recalculate Matrix 
      m_mesh_ptr = m; 
      const std::size_t mesh_size = m_mesh_ptr->size();

      // resize matrix to fit
      // m_stencil.resize(mesh_size,mesh_size);
      // set all entries to zero
      // m_stencil.setZero(); // no longer needed since setFromTriples 
      std::size_t one_sided_skirt = m_order;  
      std::size_t centered_skirt = (m_order+1)/2;  

      // number of non zeros in the sparse matrix 
      std::size_t nnz = 2*centered_skirt*(1+m_order) + (mesh_size-2*centered_skirt)*(1+2*centered_skirt);
      // resize member data 
      m_outers.resize(mesh_size+1);
      m_vals.resize(nnz);
      m_inners.resize(nnz);
      // // allocate full list of coeff triples 
      // typedef Eigen::Triplet<double> T;
      // std::vector<T> tripletList;
      // tripletList.resize(nnz); // not sure if this is correct size of lise...

      // begin OpenMP parallel section. **** COMMENTED OUT: std::vector writing is not thread safe ******
      #pragma omp parallel 
      {
        // instantiate stateful fornberg calculator. one per thread 
        FornCalc weight_calc(1+2*centered_skirt,m_order);
        // first rows with forward stencil
        #pragma omp for nowait 
        for(std::size_t i=0; i<centered_skirt; i++)
        {
          auto left = m->cbegin()+i;
          auto right = left+one_sided_skirt+1; 
          auto weights = weight_calc.GetWeights(*left, left,right, m_order); 
          std::size_t offset=i; 
          m_outers[i] = i*(1+one_sided_skirt); 
          for(auto& w : weights){
            m_inners[i*(1+one_sided_skirt)
                                      + offset] = offset;
            m_vals[i*(1+one_sided_skirt)
                                      + offset] = w; 
            offset++; 
          }
        }
        // middle rows with centered centered stencil  
        #pragma omp for nowait 
        for(std::size_t i=centered_skirt;i<mesh_size-centered_skirt; i++)
        {
          auto left = m->cbegin()-centered_skirt+i;  
          auto right = left+(centered_skirt+1+centered_skirt);
          auto weights = weight_calc.GetWeights(m->at(i), left,right, m_order);
          int offset = -centered_skirt;
          m_outers[i] = centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset); 
          for(auto& w : weights){
            m_inners[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = i+offset;
            m_vals[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = w;
            offset++; 
          };
        }
        // last rows 
        #pragma omp for nowait 
        for(std::size_t i = mesh_size-centered_skirt; i<mesh_size; i++)
        {
          auto right = m->cbegin() + i; 
          auto left = right-(one_sided_skirt+1); 
          auto weights = weight_calc.GetWeights(*right, left,right, m_order); 
          int offset= -one_sided_skirt;
          m_outers[i] = centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset); 
          for(auto& w : weights){
            m_inners[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = i+offset;  
            m_vals[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = w;  
            offset++;
          }
        }
        // #pragma omp single
        // {
        //   std::cout << "Threads in this parallel region: " << omp_get_num_threads() << "\n";
        // }
      } 
      // End OpenMP parallel section. implicit barrier 
      m_outers[m_outers.size()-1] = m_inners.size(); 
    }

    // Operators (overrides LinOpBase)
    // composition of linear of L1(L2( . ))
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      // derivative of a scalar multiple 
      if constexpr(is_scalar_multiply_expr<std::remove_cv_t<std::remove_reference_t<DerivedInner>>>::value){
        double c = InnerOp.Lhs(); 
        return c * compose(InnerOp.Rhs()); 
      }
      if constexpr(is_add_expr<std::remove_cv_t<std::remove_reference_t<DerivedInner>>>::value){
        return compose(InnerOp.Lhs()) + compose(InnerOp.Rhs()); 
      }
      // if taking derivative of another derivative 
      else if constexpr(std::is_same<std::remove_cv_t<std::remove_reference_t<DerivedInner>>,NthDerivOp>::value){
        // specialize composing to Derivatives to add order
        auto result =  NthDerivOp(m_order+InnerOp.Order(), m_mesh_ptr); 
        // if i have left a BC give it to result 
        if(lbc_ptr) result.lbc_ptr=lbc_ptr; 
        else result.lbc_ptr = InnerOp.lbc_ptr;
        // if i have a right BC give it to result 
        if(rbc_ptr) result.rbc_ptr=rbc_ptr; 
        else result.rbc_ptr = InnerOp.rbc_ptr;  
        // NthDerivOp result fully constructed
        return result; 
      }
      // all other cases use default LinOpBase implementation 
      else{
        static_assert(false, "Derivative only meant to be composed with other derivatives & scalar products"); 
        // return LinOpBase::compose(std::forward<DerivedInner>(InnerOp)); 
      }
    }; 
}; 

#endif // NthDerivOp.hpp
