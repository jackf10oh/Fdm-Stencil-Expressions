// BoundaryCondXD.hpp
//
//
//
// JAF 1/2/2025 

#ifndef BOUNDARYCONDXD_H
#define BOUNDARYCONDXD_H 

#include<vector>
#include<tuple>
#include "../DiffOps/BoundaryCond.hpp"
#include "../DiffOps/Utilities/FillStencil.hpp"
#include "../DiffOps/Utilities/SparseDiagExpr.hpp"
#include "MeshXD.hpp"

auto flat_stencil = [](const std::pair<BcPtr_t,BcPtr_t>& p, const MeshPtr_t& mesh){

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

    // set a DiscretizationXD to an implicit solution 

    // get a Matrix that to apply as a mask over XD fdm stencils 
    void SetStencilImp(MatrixStorage_t& Mat, const MeshXDPtr_t& mesh)
    {
      if(m_bc_list.empty()) throw std::runtime_error("m_bc_list must be non empty!"); 
      if(m_bc_list.size()!=mesh->dims()) throw std::invalid_argument("length of boundary condition list must be == to # of dims of MeshXDPtr"); 
      if(m_bc_list.size() > 2) throw std::runtime_error("BoundaryCondXD doesn't support dims >= 3 yet!"); 
      // resulting mask for the operation of replacing rows in Mat input
      // MatrixStorage_t mask(mesh->sizes_product());  

      MatrixStorage_t I; // will set to identity inside loop 
      
      // build up rows we are going to set in Mat -----------------------

      // first build the 1 dim case with rows set by BcPtr_t's 
      std::size_t s1=mesh->dim_size(0);  
      MatrixStorage_t mask(s1,s1); 
      m_bc_list[0].first->SetStencilL(mask,mesh->GetMesh(0));   
      m_bc_list[0].second->SetStencilR(mask,mesh->GetMesh(0)); 

      if(mesh->dims()==2)
      {
        std::size_t s2 = mesh->dim_size(1); 
        I.resize(s2,s2); I.setIdentity(); 
        // take mask into higher dimension 
        MatrixStorage_t mask_temp = Eigen::KroneckerProductSparse(I,mask);

        // fill new empty rows 
        I.resize(s1,s1); I.setIdentity(); 
        MatrixStorage_t fill_rows =  Eigen::KroneckerProductSparse(make_SparseDiag(flat_stencil(m_bc_list[1], mesh->GetMesh(1))),I); 
        fill_stencil(mask_temp,fill_rows); 
        mask = std::move(mask_temp); 
      }
      
      // overwrite rows in Mat with nonzero rows from mask 
      overwrite_stencil(Mat, mask); 
      // void return type
    }
};

#endif // BOundaryCondXD.hpp