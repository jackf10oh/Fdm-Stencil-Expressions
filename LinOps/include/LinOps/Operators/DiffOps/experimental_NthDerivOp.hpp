// experimental_NthDerivOp.hpp
//
// uses the same #ifndef macro as NthDerivOp.hpp 
// this will overwrite any occurences of NthDerivOp class
//
// JAF 2/4/2026 

#ifndef NTHDERIVOP_H
#define NTHDERIVOP_H

#include<cstdint>
#include<Eigen/Core>
#include<omp.h>

#include<Utilities/FornbergCalc.hpp> 
#include "../../LinearOpBase.hpp"

namespace LinOps{

class NthDerivOp : public LinOpMixIn<NthDerivOp>, public LinOpBase1D<NthDerivOp>
{
  public:
    // Member Data --------------------------- 
    Mesh1D_WPtr_t m_mesh_ptr; 
    std::size_t m_order; 
    std::vector<MatrixStorage_t::StorageIndex> m_inners; 
    std::vector<MatrixStorage_t::StorageIndex> m_outers; 
    std::vector<MatrixStorage_t::Scalar> m_vals; 

  public:
    // Constructors + Destructor ===================================================
    // from mesh + order 
    NthDerivOp(const LinOps::Mesh1D_SPtr_t& m, std::size_t order=1)
      : m_order(order)
    {set_mesh(m);};
    // from order 
    NthDerivOp(std::size_t order=1)
      : m_order(order)
    {};
    // destructor
    ~NthDerivOp()=default; 

    // Member Funcs ---------------------------------
    // current order or Derivative 
    std::size_t Order() const {return m_order; };

    // Matrix Getters 
    auto GetMat(){ 
      return Eigen::Map<MatrixStorage_t>(
        m_outers.size()-1, m_outers.size()-1, m_vals.size(), 
        m_outers.data(), m_inners.data(), m_vals.data(), 
        0 // 0 flags the sparse matrix as compressed 
      ); 
    }; 
    auto GetMat() const { 
      return Eigen::Map<const MatrixStorage_t>(
        m_outers.size()-1, m_outers.size()-1, m_vals.size(), 
        m_outers.data(), m_inners.data(), m_vals.data(), 
        0 // 0 flags the sparse matrix as compressed 
      );     
    };     

    // mesh getters 
    Mesh1D_WPtr_t get_weak_mesh1d() const { return m_mesh_ptr; }
    Mesh1D_SPtr_t get_mesh1d() const { return m_mesh_ptr.lock(); } 

    // set the mesh domain the derivative operator works on 
    void set_mesh(const Mesh1D_SPtr_t& m)
    {
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // throw error on nullptr 
      if(!m) throw std::runtime_error("NthDerivOp .set_mesh(Mesh1D_SPtr_t) error: Mesh1D_SPtr_t is expired!"); 
      // point to new mesh 
      m_mesh_ptr = m;  

      // calculate new matrix entries 
      const std::size_t mesh_size = m->size();

      // skirts are # of nodes to side of x_bar that are used to calculate finite difference weights 
      std::size_t one_sided_skirt = m_order;  
      std::size_t centered_skirt = (m_order+1)/2;  

      // number of non zeros in the sparse matrix 
      std::size_t nnz = 2*centered_skirt*(1+m_order) + (mesh_size-2*centered_skirt)*(1+2*centered_skirt);
      // resize member data 
      m_outers.resize(mesh_size+1);
      m_vals.resize(nnz);
      m_inners.resize(nnz);

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
}; 

} // end namespace LinOps

#endif // experimental_NthDerivOp.hpp
