// MeshXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef MESHXD_H
#define MESHXD_H 

#include<cstdint>
#include<vector>
#include<memory>
#include<algorithm>
#include<numeric>
#include "../LinOps/Mesh.hpp"

// forward declaration -> aliases
class MeshXD; 
using MeshXDPtr_t = std::shared_ptr<MeshXD>;

class MeshXD
{
  private:
    // member data --------------------------------------------------------------------------
    std::vector<MeshPtr_t> m_mesh_vec; // dynamic array of meshes. 
  public:
    // constructors ------------------------------------------------------------------------- 
    // uniformly on [0,1] with 5 steps per dim, with n_dims
    MeshXD(std::size_t n_dims_init) 
      : m_mesh_vec(n_dims_init)
    {
      for(std::size_t dim_i=0; dim_i<n_dims_init; dim_i++){
        m_mesh_vec[dim_i] = std::make_shared<Mesh1D>(0.0, 1.0, 5); 
      }; 
    };
    // uniformly on [l,r] with n_steps per dim, with n_dims
    MeshXD(double left=0.0, double right=1.0, std::size_t n_steps=11,std::size_t n_dims=1)
      : m_mesh_vec(n_dims)
    {
      MeshPtr_t original_ptr = std::make_shared<Mesh1D>(left,right,n_steps); 
      for(auto& ptr : m_mesh_vec) ptr=original_ptr;  
    }
    // uniformly on [ L(i), R(i) ] with NSteps(i) per dim with NSteps.size() dims 
    MeshXD(const std::vector<std::pair<double,double>>& axes_vec, const std::vector<std::size_t>& nsteps_vec)
      : m_mesh_vec(axes_vec.size())
    {
      if(axes_vec.size()!=nsteps_vec.size()) throw std::invalid_argument("# of axes != # of steps given");
      auto nstep_it = nsteps_vec.begin();
      auto axes_it=axes_vec.begin();
      for(auto write=m_mesh_vec.begin(); write!=m_mesh_vec.end(); write++){
        *write = std::make_shared<Mesh1D>(axes_it->first, axes_it->second, *nstep_it); 
        axes_it++; 
        nstep_it++; 
      }
    }
    // Copy 
    MeshXD(const MeshXD& other)=default; 
    // destructors --------------------------------------------------------------------------
    virtual ~MeshXD()=default;
    // member functions ---------------------------------------------------------------------
    // full size of XD mesh. i.e. axis1.size() * ... * axisn.size()
    std::size_t sizes_product(){return std::accumulate(m_mesh_vec.begin(), m_mesh_vec.end(),1ul, [](std::size_t rolling, const MeshPtr_t& axis){return rolling*(axis->size());});} 
    // product of axes up to dim exclusively [first, dim)
    std::size_t sizes_middle_product(std::size_t start, std::size_t end){
      if(start > end) throw std::invalid_argument("start index must be <= end index for middle product"); 
      if(end > m_mesh_vec.size()) throw std::invalid_argument("end index must be <= # of dims in MeshXD"); 
      std::size_t prod = 1; 
      for(std::size_t i=start; i<end; i++) prod *= m_mesh_vec[i]->size(); 
      return prod; 
    } 
    // size of a specific axis 
    std::size_t dim_size(std::size_t i){return m_mesh_vec.at(i)->size();} 
    // number of dimensions 
    std::size_t dims(){return m_mesh_vec.size(); } 
    // get a specific mesh 
    MeshPtr_t& GetMesh(std::size_t i){return m_mesh_vec[i];} 
    const MeshPtr_t& GetMesh(std::size_t i) const {return m_mesh_vec[i];} 
    MeshPtr_t& GetMeshAt(std::size_t i){return m_mesh_vec.at(i);}
    const MeshPtr_t& GetMeshAt(std::size_t i) const {return m_mesh_vec.at(i);}

};

#endif // MeshXD.hpp