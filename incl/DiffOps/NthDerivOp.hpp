// NthDerivOp.hpp
//
//
//
// JAF 12/5/2025

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<eigen3/Eigen/Core>

#include "../Utilities/Fornberg.hpp"
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
    Eigen::MatrixXd& GetMat(){ return m_stencil; }; 
    const Eigen::MatrixXd& GetMat() const { return m_stencil; }; 
    Discretization1D apply(const Discretization1D& d) const
    {
      Discretization1D result(d.mesh()); 
      result = m_stencil*d.values(); 
      return result; 
    }
    
    // set the mesh the derivative operator works on 
    void set_mesh(MeshPtr_t m)
    {
      //  do nothing if meshes are ==
      if(m==m_mesh_ptr) return;  
      // resize 0 on nullptr
      if(m==nullptr){m_stencil.resize(0,0); return;};
      
      m_mesh_ptr = m; 

      // resize matrix to fit
      m_stencil.resize(m->size(),m->size());
      // set all entries to zero
      m_stencil.setZero(); 
      // first rows with forward stencil 
      int skirt = m_order; 
      int i=0; 
      auto left = m->cbegin();
      auto right = left+skirt+1; 
      for(int end=(m_order+1)/2; i<end; i++, left++, right++)
      {
        auto weights = FornWeights(m->at(i), m_order, std::vector(left,right)); 
        std::copy(weights.begin(), 
                  weights.end(), 
                  m_stencil.row(i).begin()+i
        );
      }
      // middle rows with centered centered stencil  
      skirt = (m_order+1)/2; // num of points to left/right of center for stencil 
      i=skirt; 
      left = m->cbegin(); 
      right = left+(skirt+1+skirt); 
      for(;i<m->size()-skirt; i++, left++, right++)
      {
        auto weights = FornWeights(m->at(i),m_order,std::vector(left,right)); 
        std::copy(weights.begin(),
                  weights.end(), 
                  m_stencil.row(i).begin()+(i-skirt)
        );
      }
      // last rows 
      skirt = m_order; 
      i = m->size()-1; 
      right = m->cend(); 
      left = right-skirt-1; 
      for(; i>=m->size()-(m_order+1)/2; i--,left--,right--)
      {
        auto weights = FornWeights(m->at(i), m_order, std::vector(left,right)); 
        std::copy(weights.rbegin(),
                  weights.rend(),
                  m_stencil.row(i).reverse().begin()+std::distance(right,m->cend())
        );
      }
    }

    // Operators (overrides LinOpBase)
    // composition of linear of L1(L2( . ))
    template<typename DerivedInner> 
    auto compose(DerivedInner&& InnerOp)
    {
      // if taking derivative of product u*v
      if constexpr(is_compose_expr<std::remove_cv_t<std::remove_reference_t<DerivedInner>>>::value){
        decltype(auto) u = InnerOp.Lhs(); 
        decltype(auto) v = InnerOp.Rhs(); 
        return compose(u).compose(v) + u.compose(compose(v));
      }
      // if taking derivative of scalar multiple c*u 
      else if constexpr(is_scalar_multiply_expr<std::remove_cv_t<std::remove_reference_t<DerivedInner>>>::value){
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
        return LinOpBase::compose(std::forward<DerivedInner>(InnerOp)); 
      }
    }; 
}; 

#endif // NthDerivOp.hpp