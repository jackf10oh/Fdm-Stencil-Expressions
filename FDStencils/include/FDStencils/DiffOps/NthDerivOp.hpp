// NthDerivOp.hpp
//
//
//
// JAF 12/5/2025

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<cstdint>
#include<Eigen/Core>

#include<Utilities/FornbergCalc.hpp>
#include "../FdmPlugin.hpp"
#include<LinOps/LinearOpBase.hpp>

namespace Fds{

class NthDerivOp : public LinOps::LinOpBase<NthDerivOp>
{
  private:
    // member data 
    MatrixStorage_t m_stencil; 
    std::size_t m_order; 
  public:
    // constructors ---------------------------------
    NthDerivOp(const LinOps::Mesh1D_SPtr_t& m, std::size_t order=1)
      : m_order(order)
    {set_mesh(m);};
    NthDerivOp(std::size_t order=1)
      : m_order(order)
    {};
    
    // Destructors -----------------------------------
    ~NthDerivOp()=default; 
    
    // Member Funcs ---------------------------------
    std::size_t Order() const {return m_order; };
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };     
    // set the mesh the derivative operator works on 
    void set_mesh(const LinOps::Mesh1D_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;

      // store the mesh
      m_mesh_ptr = m;  

      // perform work on m
      const std::size_t mesh_size = m->size();

      // resize matrix to fit
      m_stencil.resize(mesh_size,mesh_size);
      // set all entries to zero
      // m_stencil.setZero(); // no longer needed since setFromTriples 
      std::size_t one_sided_skirt = m_order;  
      std::size_t centered_skirt = (m_order+1)/2;  

      // allocate full list of coeff triples 
      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.resize(2*centered_skirt*(1+m_order) + (mesh_size-2*centered_skirt)*(1+2*centered_skirt)); // not sure if this is correct size of lise...

      // begin OpenMP parallel section. **** COMMENTED OUT: std::vector writing is not thread safe ******
      // #pragma omp parallel 
      {
        // instantiate stateful fornberg calculator. one per thread 
        FornCalc weight_calc(1+2*((m_order+1)/2),m_order);
        // first rows with forward stencil
        // #pragma omp for nowait 
        for(std::size_t i=0; i<centered_skirt; i++)
        {
          auto left = m->cbegin()+i;
          auto right = left+one_sided_skirt+1; 
          auto weights = weight_calc.GetWeights(*left, left,right, m_order); 
          std::size_t offset=i; 
          for(auto& w : weights){
            tripletList[i*(1+one_sided_skirt) 
                        + offset] = T(i,offset,w);
            offset++; 
          }
        }
        // middle rows with centered centered stencil  
        // #pragma omp for nowait 
        for(std::size_t i=centered_skirt;i<mesh_size-centered_skirt; i++)
        {
          auto left = m->cbegin()-centered_skirt+i;  
          auto right = left+(centered_skirt+1+centered_skirt);
          auto weights = weight_calc.GetWeights(m->at(i), left,right, m_order);
          int offset = -centered_skirt;
          for(auto& w : weights){
            tripletList[centered_skirt*(1+one_sided_skirt) 
                        + (i-centered_skirt)*(1+2*centered_skirt)
                        +(centered_skirt+offset)] = T(i,i+offset,w);
            offset++; 
          };
        }
        // last rows 
        // #pragma omp for nowait 
        for(std::size_t i = mesh_size-centered_skirt; i<mesh_size; i++)
        {
          auto right = m->cbegin() + i; 
          auto left = right-(one_sided_skirt+1); 
          auto weights = weight_calc.GetWeights(*right, left,right, m_order); 
          int offset= -one_sided_skirt;
          for(auto it=weights.cbegin(); it!=weights.cend(); it++){
            tripletList[centered_skirt*(1+one_sided_skirt)
                        + (mesh_size-2*centered_skirt)*(1+2*centered_skirt) 
                        + (i-(mesh_size-centered_skirt))*(1+one_sided_skirt)
                        + (one_sided_skirt+offset)] = T(i,i+offset,*it);  
            offset++;
          }
        }
      } 
      // End OpenMP parallel section. implicit barrier 
      
      // Eigen 5.0 :(. might be faster if sorting time > time saved from parallelism 
      // m_stencil.insertFromSortedTriplets()(tripletList.begin(), tripletList.end());   
      m_stencil.setFromTriplets(tripletList.begin(), tripletList.end()); 
    }

    // "Operators" (overrides LinOpBase)
    // composition of linear of L1(L2( . ))
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      // derivative of a scalar multiple 
      using cleaned_rhs_t = std::remove_cv_t<std::remove_reference_t<DerivedInner>>; 
      if constexpr(LinOps::traits::is_scalar_multiply_expr<cleaned_rhs_t>::value){
        double c = InnerOp.Lhs(); 
        return c * compose(InnerOp.Rhs()); 
      }
      else if constexpr(LinOps::traits::is_add_expr<cleaned_rhs_t>::value){
        return compose(InnerOp.Lhs()) + compose(InnerOp.Rhs()); 
      }
      else if constexpr(LinOps::traits::is_subtraction_expr<cleaned_rhs_t>::value){
        return compose(InnerOp.Lhs()) - compose(InnerOp.Rhs()); 
      }
      else if constexpr(LinOps::traits::is_negation_expr<cleaned_rhs_t>::value){
        return - compose(InnerOp.Lhs()); 
      }
      // if taking derivative of another derivative 
      else if constexpr(std::is_same<cleaned_rhs_t,NthDerivOp>::value){
        // specialize composing to Derivatives to add order
        return NthDerivOp(m_order+InnerOp.Order()); 
      }
      // all other cases use default LinOpBase implementation 
      else{
        static_assert(false, "Derivative only meant to be composed with other derivatives & scalar products"); 
        // return LinOpBase::compose(std::forward<DerivedInner>(InnerOp)); 
      }
    }; 

}; 

} // end namespace Fds 

#endif // NthDerivOp.hpp
