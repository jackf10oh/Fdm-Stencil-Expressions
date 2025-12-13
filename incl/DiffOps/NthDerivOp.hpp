// NthDerivOp.hpp
//
//
//
// JAF 12/5/2025

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<eigen3/Eigen/Core>

#include "../Utilities/Fornberg.hpp"
#include "../Utilities/FornbergCalc.hpp"
#include "../../LinOps/Discretization.hpp"
#include "../FdmPlugin.hpp"
#include "../../LinOps/LinearOpBase.hpp"

class NthDerivOp : public LinOpBase<NthDerivOp>
{
  private:
    // member data 
    std::size_t m_order; 
  public:
    // constructors ---------------------------------
    NthDerivOp(MeshPtr_t m, std::size_t order=1)
      : m_order(order)
    {set_mesh(m);};
    NthDerivOp(std::size_t order=1, MeshPtr_t m=nullptr)
      : m_order(order)
    {set_mesh(m);};
    
    // Destructors -----------------------------------
    ~NthDerivOp()=default; 
    
    // Member Funcs ---------------------------------
    std::size_t Order() const {return m_order; };
    MatrixStorage_t& GetMat(){ return m_stencil; }; 
    const MatrixStorage_t& GetMat() const { return m_stencil; };     
    // set the mesh the derivative operator works on 
    void set_mesh(MeshPtr_t m)
    {
      //  do nothing if meshes are ==
      if(m==m_mesh_ptr) return;  
      // resize 0 on nullptr
      if(m==nullptr){m_stencil.resize(0,0); return;};
      
      // new pointer set. proceed to recalculate Matrix 
      m_mesh_ptr = m; 

      // resize matrix to fit
      m_stencil.resize(m->size(),m->size());
      // set all entries to zero
      // m_stencil.setZero(); // no longer needed since setFromTriples 

      typedef Eigen::Triplet<double> T;
      std::vector<T> tripletList;
      tripletList.reserve(2*m_order*m_stencil.rows()); // not sure if this is correct size of lise...

      // instantiate stateful fornberg calculator 
      FornCalc weight_calc(1+2*((m_order+1)/2),m_order);

      // first rows with forward stencil 
      std::size_t skirt = m_order; 
      std::size_t i=0; 
      auto left = m->cbegin();
      auto right = left+skirt+1; 
      for(std::size_t end=(m_order+1)/2; i<end; i++, left++, right++)
      {
        auto weights = weight_calc.GetWeights(m->at(i),left,right,m_order); 
        std::size_t offset=0; 
        for(auto& w : weights){
          tripletList.push_back(T(i,offset++,w));
        }
      }
      // middle rows with centered centered stencil  
      skirt = (m_order+1)/2; // num of points to left/right of center for stencil 
      i=skirt; 
      left = m->cbegin(); 
      right = left+(skirt+1+skirt); 
      for(;i<m->size()-skirt; i++, left++, right++)
      {
        auto weights = weight_calc.GetWeights(m->at(i),left,right,m_order);
        int offset = -skirt;
        for(auto& w : weights){
          tripletList.push_back(T(i,i+offset++,w));
        };
      }
      // last rows 
      skirt = m_order; 
      i = m->size()-1; 
      right = m->cend(); 
      left = right-skirt-1; 
      for(; i>m->size()-(m_order+1)/2-1; i--,left--,right--)
      {
        auto weights = weight_calc.GetWeights(m->at(i),left,right,m_order); 
        int offset=0;
        for(auto it=weights.rbegin(); it!=weights.rend(); it++){
          tripletList.push_back(T(i,i-offset++,*it)); 
        }
      }
      m_stencil.setFromTriplets(tripletList.begin(), tripletList.end()); 
    }

    // Operators (overrides LinOpBase)
    // composition of linear of L1(L2( . ))
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      if constexpr(is_scalar_multiply_expr<std::remove_cv_t<std::remove_reference_t<DerivedInner>>>::value){
        double c = InnerOp.Lhs(); 
        return c * compose(InnerOp.Rhs()); 
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