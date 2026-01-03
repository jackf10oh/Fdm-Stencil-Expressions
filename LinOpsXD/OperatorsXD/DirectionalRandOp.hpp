// DirectionalRandOp.hpp
//
//
//
// JAF 1/3/2025 

#ifndef DIRECTIONALRANDOP_H
#define DIRECTIONALRANDOP_H

#include<iostream>
#include "../LinearOpXDBase.hpp" 

class DirectionalRandOp: public LinOpXDBase<DirectionalRandOp> 
{
  private:
    // private type defs
    typedef typename Eigen::SparseMatrix<double, Eigen::RowMajor> Matrix_t; 
    // member data 
    Matrix_t m_mat; 
    std::size_t m_direction; // which Mesh1D the operator acts on. 
  public:
    // constructors
    // direction only 
    DirectionalRandOp(std::size_t dir_init): m_direction(dir_init){m_mesh_ptr=nullptr;}; 
    // mesh + direction 
    DirectionalRandOp(MeshXDPtr_t m=nullptr, std::size_t dir_init=0) 
      : m_direction(dir_init) 
    {set_mesh(m);}; 


    // destructors 
    ~DirectionalRandOp()=default; 
    // member functions 
    auto& GetMat(){ return m_mat; }; 
    const auto& GetMat() const { return m_mat; }; 
    DiscretizationXD apply(const DiscretizationXD& d){ return d; }; 
    void set_mesh(MeshXDPtr_t m){
      // do nothing on nullptr or same ptr 
      if(m==nullptr || m==m_mesh_ptr) return; 
      if(m->dims() < m_direction) throw std::runtime_error("MeshXD must have dims >= m_direction of DirectionalRandOp!"); 
      m_mesh_ptr = m; 

      std::size_t s = m_mesh_ptr->dim_size(m_direction); 
      auto temp = Eigen::MatrixXd::Random(s,s); 
      Matrix_t R = temp.sparseView(); 

      if(m_mesh_ptr->dims()>m_direction)
      {
        // product of last meshes > direction 
        std::size_t prod = 1;
        for(std::size_t i=m_direction+1; i<m_mesh_ptr->dims(); i++) prod *= m_mesh_ptr->dim_size(i); 
        
        // make identity matrix for last meshes after direction 
        Matrix_t I(prod,prod); I.setIdentity(); 

        // update R according to kronecker 
        Matrix_t intermediate = Eigen::KroneckerProductSparse(I,R); 
        R = std::move(intermediate); 
      }

      if(m_direction>0)
      {
        // product of first meshes < direction 
        std::size_t prod = 1;
        for(std::size_t i=0; i<m_direction; i++) prod *= m_mesh_ptr->dim_size(i); 
          
        // make identity matrix for first meshes before direction 
        Matrix_t I(prod,prod); I.setIdentity(); 

        // update R according to kronecker 
        Matrix_t intermediate = Eigen::KroneckerProductSparse(R,I); 
        R = std::move(intermediate); 
      }
      m_mat = std::move(R);
    } // end set_mesh()  
}; 

#endif 