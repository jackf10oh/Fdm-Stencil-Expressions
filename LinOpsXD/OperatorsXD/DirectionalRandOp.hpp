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
      if(!locked) return; 
      m_mesh_ptr = m; // store the mesh  
      // perform work on locked 

      // identity matrix for calculations later 
      Matrix_t I; 

      // get a random matrix R that matches specific direction 
      std::size_t s = locked->dim_size(m_direction); 
      auto temp = Eigen::MatrixXd::Random(s,s); 
      Matrix_t R = temp.sparseView(); 

      if(locked->dims()>m_direction)
      {
        // product of last meshes > direction 
        std::size_t prod = locked->sizes_middle_product(m_direction+1, locked->dims()); 
        
        // resize identity matrix for last meshes after direction 
        I.resize(prod,prod); I.setIdentity(); 

        // update R according to kronecker 
        Matrix_t intermediate = Eigen::KroneckerProductSparse(I,R); 
        R = std::move(intermediate); 
      }

      if(m_direction>0)
      {
        // product of first meshes < direction 
        std::size_t prod = locked->sizes_middle_product(0, m_direction); 
          
        // resize identity matrix for first meshes before direction 
        I.resize(prod,prod); I.setIdentity(); 

        // update R according to kronecker 
        Matrix_t intermediate = Eigen::KroneckerProductSparse(R,I); 
        R = std::move(intermediate); 
      }
      m_mat = std::move(R); // take ownership of high dimensional R with m_mat member data. 
    } // end set_mesh()  
}; 

} // end namespace LinOps 

#endif 