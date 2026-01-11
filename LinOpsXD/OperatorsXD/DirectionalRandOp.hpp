// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H

#include<iostream>
#include<Eigen/Sparse> 
#include<unsupported/Eigen/KroneckerProduct> 

#include "../LinearOpBaseXD.hpp" 

namespace LinOps{

class DirectionalRandOp: public LinOpBaseXD<DirectionalRandOp> 
{
  private:
    // private type defs
    typedef typename CUSTOM_LINOPSXD_SPARSE_MATRIX_STORAGE Matrix_t; 
    // member data 
    Matrix_t m_mat; 
    std::size_t m_direction; // which Mesh1D the operator acts on. 
  public:
    // constructors
    // direction only 
    DirectionalRandOp(std::size_t dir_init): m_direction(dir_init){ m_mesh_ptr=MeshXDPtr_t{}; }  
    // mesh + direction 
    DirectionalRandOp(MeshXDPtr_t m=MeshXDPtr_t{}, std::size_t dir_init=0) 
      : m_direction(dir_init) 
    {set_mesh(m);}; 


    // destructors 
    ~DirectionalRandOp()=default; 
    // member functions 
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 
    DiscretizationXD apply(const DiscretizationXD& d){ 

      // result copies mesh and dims from d 
      DiscretizationXD result;
      result.mesh() = d.mesh(); 
      result.dims_list() = d.dims_list();  

      // use move operator from result Eigen::VectorXd
      result = m_mat * d.values(); 

      return result; 
    }; 
    
    void set_mesh(MeshXDPtr_t m){
      // ensure we aren't resetting the mesh again
      if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // do nothing on nullptr. or throw an error 
      auto locked = m.lock(); 
      // do nothing on nullptr
      if(!locked) return; 
      // store the mesh
      m_mesh_ptr = m;   

      // perform work on locked ...
      
      // check this->direction < mesh-> # dims 
      if(m_direction >= locked->dims()) throw std::runtime_error("DirectionalNthDerivOp.set_mesh() error: direction >= MeshXD.dims()"); 

      // get a random matrix R that matches specific direction 
      std::size_t s = locked->dim_size(m_direction); 
      Matrix_t R = Eigen::MatrixXd::Random(s,s).sparseView();      
      // identity matrix for calculations later 
      Matrix_t I;  

      std::size_t prod_before = locked->sizes_middle_product(0,m_direction); 
      std::size_t prod_after = (locked->dims() > m_direction) ? locked->sizes_middle_product(m_direction+1, locked->dims()) : 1; 

      if(prod_before>1){
        I.resize(prod_before, prod_before); 
        I.setIdentity(); 
        m_mat = Eigen::KroneckerProductSparse(R, I); 
      }
      else{
        m_mat = std::move(R); 
      }; 
      if(prod_after>1){
        I.resize(prod_after, prod_after); 
        I.setIdentity(); 
        Matrix_t temp = Eigen::KroneckerProductSparse(I, m_mat);  
        m_mat = std::move(temp); 
      }; 
    } // end set_mesh()  
}; 

} // end namespace LinOps 

#endif 