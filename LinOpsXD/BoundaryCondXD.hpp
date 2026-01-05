// BoundaryCondXD.hpp
//
//
//
// JAF 1/2/2025 

#ifndef BOUNDARYCONDXD_H
#define BOUNDARYCONDXD_H 

#include<iostream>
#include<vector>
#include<tuple>
#include<unsupported/Eigen/KroneckerProduct>
#include "../DiffOps/BoundaryCond.hpp"
#include "../Utilities/FillStencil.hpp"
#include "../Utilities/SparseDiagExpr.hpp"
#include "MeshXD.hpp"

class BoundaryCondXD
{
  public:
    // member data. list of boundary conditions. 1 per Dimension 
    std::vector<std::pair<BcPtr_t,BcPtr_t>> m_bc_list; 

    // Constructors
    BoundaryCondXD()=default; 
    BoundaryCondXD(const BoundaryCondXD& other)=default;
    // Destructors 
    ~BoundaryCondXD()=default; 

    // set a DiscretizationXD to an explicit solution 
    void SetSol(DiscretizationXD& Sol, const MeshXDPtr_t& mesh){
      // check args are compaitble ---------------------------------
      if(m_bc_list.empty()) throw std::runtime_error("BoundaryCondXD: m_bc_list must be non empty!"); 
      if(m_bc_list.size()!=mesh->dims()) throw std::invalid_argument("BoundaryCondXD:  length of boundary condition list must be == to # of dims of MeshXDPtr");
      if(m_bc_list.size()!=Sol.dims()) throw std::invalid_argument("BoundaryCondXD: length of boundary condition list must be == to # of dims of DiscretizationXD");
      for(std::size_t i=0; i<Sol.dims(); i++){
        if(Sol.dim_size(i)!=mesh->dim_size(i)) throw std::invalid_argument("BoundaryCondXD: Dimension sizes of Discretization1D and MeshXD must match!");
      };

      // this lambda takes 1 BcPtr_t and applies it to Sol
      // without double assigning to corners/edges of XDim space 
      auto set_dim_boundaries = [&](std::size_t dim, const std::pair<BcPtr_t,BcPtr_t>& bc_pair){
        auto mesh_1dim = mesh->GetMesh(dim); 
        auto views = Sol.OneDim_views(dim); 
        // iterate through views that look like Mesh1D
        for(std::size_t i=0; i<views.size(); i++){
          // determine if it has been set by lower dimension 
          bool set_by_low_dim = false; 
          std::size_t s1 = 1; 
          for(int dim_i=0; dim_i<dim; dim_i++){
            std::size_t s2 = Sol.dim_size(dim_i); 
            // in first group of values
            if((i/s1)%s2 == 0) set_by_low_dim = true;  
            // in last group of values
            if((i/s1)%s2 == s2-1) set_by_low_dim = true;  
            // next check uses a large bucket 
            s1 *= s2; 
          }
          // if it hasn't, use BC on it. 
          if(!set_by_low_dim){
            bc_pair.first->SetSolL(views[i], mesh_1dim); 
            bc_pair.second->SetSolR(views[i], mesh_1dim); 
          }
        }
      }; // end set_dim_boundaries lambda 

      // apply lambda + corresponding bc from m_bc_list to Sol 
      for(std::size_t j=0; j<m_bc_list.size(); j++) set_dim_boundaries(j,m_bc_list[j]); 

      // void return type       
    }

    // set a DiscretizationXD to an implicit solution 
    void SetImpSol(DiscretizationXD& Sol, const MeshXDPtr_t& mesh){
      // check args are compaitble ---------------------------------
      if(m_bc_list.empty()) throw std::runtime_error("BoundaryCondXD: m_bc_list must be non empty!"); 
      if(m_bc_list.size()!=mesh->dims()) throw std::invalid_argument("BoundaryCondXD:  length of boundary condition list must be == to # of dims of MeshXDPtr");
      if(m_bc_list.size()!=Sol.dims()) throw std::invalid_argument("BoundaryCondXD: length of boundary condition list must be == to # of dims of DiscretizationXD");
      for(std::size_t i=0; i<Sol.dims(); i++){
        if(Sol.dim_size(i)!=mesh->dim_size(i)) throw std::invalid_argument("BoundaryCondXD: Dimension sizes of Discretization1D and MeshXD must match!"); 
      }; 

      // this lambda takes 1 BcPtr_t and applies it to Sol
      // without double assigning to corners/edges of XDim space 
      auto set_dim_boundaries_imp = [&](std::size_t dim, const std::pair<BcPtr_t,BcPtr_t>& bc_pair){
        auto mesh_1dim = mesh->GetMesh(dim); 
        auto views = Sol.OneDim_views(dim); 
        // iterate through views that look like Mesh1D
        for(std::size_t i=0; i<views.size(); i++){
          // determine if it has been set by lower dimension 
          bool set_by_low_dim = false; 
          std::size_t s1 = 1; 
          for(int dim_i=0; dim_i<dim; dim_i++){
            std::size_t s2 = Sol.dim_size(dim_i); 
            // in first group of values
            if((i/s1)%s2 == 0) set_by_low_dim = true;  
            // in last group of values
            if((i/s1)%s2 == s2-1) set_by_low_dim = true;  
            // next check uses a large bucket 
            s1 *= s2; 
          }
          // if it hasn't, use BC on it. 
          if(!set_by_low_dim){
            bc_pair.first->SetImpSolL(views[i], mesh_1dim); 
            bc_pair.second->SetImpSolR(views[i], mesh_1dim); 
          }
        }
      }; // end set_dim_boundaries lambda 
      // apply lambda + corresponding bc from m_bc_list to Sol 
      std::size_t dim_i=0; 
      for(auto& bc : m_bc_list) set_dim_boundaries_imp(dim_i++, bc); 
      // void return type       
    }

    // set a Matrixs' row according to m_bc_list. making it an implicit stencil  
    void SetStencilImp(MatrixStorage_t& Mat, const MeshXDPtr_t& mesh)
    {
      // check args are compaitble ---------------------------------
      if(m_bc_list.empty()) throw std::runtime_error("m_bc_list must be non empty!"); 
      if(m_bc_list.size()!=mesh->dims()) throw std::invalid_argument("length of boundary condition list must be == to # of dims of MeshXDPtr"); 
      // if(m_bc_list.size() > 2) throw std::runtime_error("BoundaryCondXD doesn't support dims >= 3 yet!"); 

      // initializations ---------------------------------------------
      // s1 = ith dims size, s2 = cumulative product of all dims < ith dim
      std::size_t s1=mesh->dim_size(0), s2=1; 
      // lambda to give an NxN identity matrix 
      MatrixStorage_t cap_mat; 
      auto I = [&cap_mat](std::size_t N) -> const MatrixStorage_t& {cap_mat.resize(N,N); cap_mat.setIdentity(); return cap_mat;}; 
      // resulting mask for replacing rows in Mat input ---------------
      MatrixStorage_t mask(s1,s1);  
      
      // base case: 
      // set stencil's first/last rows maually with BcPtr_t  
      m_bc_list[0].first->SetStencilL(mask,mesh->GetMesh(0));   
      m_bc_list[0].second->SetStencilR(mask,mesh->GetMesh(0)); 

      // recursive case: 
      for(std::size_t ith_dim=1; ith_dim<mesh->dims(); ith_dim++)
      {
        // take previous mesh to higher dimension 
        s1 = mesh->dim_size(ith_dim); 
        MatrixStorage_t mask_temp = Eigen::KroneckerProductSparse(I(s1), mask); 
        
        // fill new empty rows with ith BcPtr_t pair
        MatrixStorage_t sp_diag = make_SparseDiag( flat_stencil(m_bc_list[ith_dim],mesh->GetMesh(ith_dim)), s2 ); 
        MatrixStorage_t fill_rows = Eigen::KroneckerProductSparse(sp_diag, I(s1));
        fill_stencil(mask_temp, fill_rows); 

        // take ownership of higher dim temp with mask 
        mask = std::move(mask_temp); 

        // // increment the cumulative product 
        s2 *= s1; 
      }

      // overwrite rows in Mat with nonzero rows from mask 
      overwrite_stencil(Mat, mask); 

      // void return type
    }
  private:
    MatrixStorage_t flat_stencil(const std::pair<BcPtr_t,BcPtr_t>& p, const MeshPtr_t& mesh){

    // empty singuglar row of size == mesh size 
    MatrixStorage_t result(1, mesh->size()); 

    // set the row according to left boundary condition. 
    p.first->SetStencilL(result, mesh); 

    // store it into a temp vector 
    std::vector<double> temp(result.valuePtr(), result.valuePtr()+result.nonZeros()); 

    // set the row according to right boundary condition 
    p.second->SetStencilR(result, mesh); 

    // copy left side back in to result 
    for(auto i=0; i<temp.size(); i++) result.coeffRef(0,i) = temp[i]; 

    return result; 
  }; 
};

#endif // BOundaryCondXD.hpp 


/* 
// order of priority for m_bc_list to be applied 
corners / edges of XDim space are given priority to whichever BC would cover it first 
0 is unaffected by BCs. 
1 is first BC in m_bc_list,
2 is 2nd,
...  

// e.g. in 2D --------------------------
1 1 1 1 1 1 1
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
1 1 1 1 1 1 1


// e.g. in 3D -----------------------

// slice of bottom 
1 1 1 1 1 1 1
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
2 3 3 3 3 3 2
1 1 1 1 1 1 1

// slice directly off bottom. (3 no longer applies)
1 1 1 1 1 1 1
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
2 0 0 0 0 0 2
1 1 1 1 1 1 1

*/