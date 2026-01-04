// DiscretizationXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef DISCRETIZATIONXD_H
#define DISCRETIZATIONXD_H

#include<vector>
#include<Eigen/Core>
#include "MeshXD.hpp"

struct DiscretizationXD
{
  private:
    // private typedefs ==================================================== 
    typedef Eigen::Stride<0,Eigen::Dynamic> Stride_t; 
    typedef Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t> StrideView_t; 

    // member data ==========================================================
    Eigen::VectorXd m_vals; // flattened array of values 
    std::vector<std::size_t> m_dims; 
    MeshXDPtr_t m_mesh_ptr; 

  public:
    // constructors ==========================================================
    // Default ----------------------------------------------- 
    DiscretizationXD()=default; 
    // from size of m_vals. assume just 1 dimension ----------
    DiscretizationXD(std::size_t size_init): m_mesh_ptr(), m_vals(size_init), m_dims(1,size_init){}; 
    // from a MeshXDPtr --------------------------------------
    DiscretizationXD(const MeshXDPtr_t& mesh_init) 
      : m_mesh_ptr(mesh_init), m_vals(mesh_init ? mesh_init->sizes_product() : 0), m_dims(mesh_init ? mesh_init->dims() : 0)
    {
      if(m_mesh_ptr){
        for(std::size_t i=0; i<m_mesh_ptr->dims(); ++i) m_dims[i] = m_mesh_ptr->dim_size(i); 
      }
    }
    // copy constructor --------------------------------------
    DiscretizationXD(const DiscretizationXD& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals), m_dims(other.m_dims)
    {}; 
    // destructors ==========================================================
    ~DiscretizationXD()=default; 
    // member functions ==========================================================
    // get underlying values 
    Eigen::VectorXd& values(){return m_vals; }
    const Eigen::VectorXd& values() const {return m_vals; } 

    // --------------------------------------------------- ????????????? 
    std::vector<StrideView_t> dim_values_view(std::size_t ith_dim=0)
    {
      // i has to be one of the dimensions of DiscretizationXD 
      if(ith_dim >= dims()) throw std::invalid_argument("Discretization1D::dim_values_view(i) i must be < Discretization1D.dims().");

      std::size_t ith_dim_size = m_dims[ith_dim]; 
      std::size_t num_copies = sizes_product() / ith_dim_size; 
      std::size_t mod = sizes_middle_product(0, ith_dim); 
      std::size_t scale = mod * ith_dim_size; 

      std::vector<StrideView_t> result; 
      result.reserve(m_dims[ith_dim]); 

      Stride_t stride(0,mod); 

      // iterate through the copies 
      for(std::size_t n=0; n<num_copies; n++)
      {
        // offset from start of current copy 
        std::size_t offset = (mod ? n % mod : n) + (scale * (n/mod));  
        // begin data ptr of copy  
        auto begin = m_vals.data()+offset;  

        // MemView of current copy 
        result.emplace_back(begin, ith_dim_size, stride);
      }
      return result; 
    }

    // get underlying list of dim sizes
    std::vector<std::size_t>& dims_list(){ return m_dims; }
    const std::vector<std::size_t>& dims_list() const { return m_dims; } 

    // get underlying MeshXDPtr_t 
    MeshXDPtr_t& mesh(){return m_mesh_ptr;}; 
    const MeshXDPtr_t& mesh() const {return m_mesh_ptr;}; 

    // number of dimensions 
    std::size_t dims() const{ return m_dims.size(); };
    // size of ith dimension  
    std::size_t dim_size(std::size_t i) const {return m_dims[i]; }
    
    // product of all dimensions' sizes 
    std::size_t sizes_product() const { return std::accumulate(m_dims.begin(), m_dims.end(), 1, std::multiplies{}); } 
    // product of dimensions in [start,end)
    std::size_t sizes_middle_product(std::size_t start, std::size_t end){
      if(start > end) throw std::invalid_argument("start index must be <= end index for middle product"); 
      if(end > m_dims.size()) throw std::invalid_argument("end index must be <= # of dims in MeshXD"); 
      return std::accumulate(m_dims.begin()+start, m_dims.begin()+end, 1, std::multiplies{}); 
    }

    // set discretization to same size as meshxd's sizes_product
    void match_mesh(MeshXDPtr_t m) { 
      if(!m) return; // do nothing on nullptr 
      m_mesh_ptr=m; 
      m_vals.resize(m->sizes_product()); 
      m_dims.resize(m->dims()); 
      for(std::size_t i=0; i<m->dims(); ++i) m_dims[i] = m->dim_size(i); 
    }
    // set discretization to a constant
    void set_init(double val){ m_vals.setConstant(val);}
    // set vector to match a mesh size and set it constant 
    void match_mesh(MeshXDPtr_t m, double val){ match_mesh(m), m_vals.setConstant(val); } 

    // Operators ----------------------------------------------------
    DiscretizationXD& operator=(const DiscretizationXD& other) = default;
    DiscretizationXD& operator=(DiscretizationXD&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      other.m_mesh_ptr = nullptr; 
      m_vals = std::move(other.m_vals); 
      m_dims = std::move(other.m_dims);
      return *this;
    }; 
    DiscretizationXD& operator=(const Eigen::VectorXd& other){
      // only care about the size of the VectorXd if we already have a MeshXDPtr_t 
      if(m_mesh_ptr != nullptr)
      {
        if(other.size() != m_mesh_ptr->sizes_product()) throw std::invalid_argument("size of VectorXd must be == to m_mesh_ptr->sizes_product()"); 
      } 
      m_vals=other;
      // assume that m_dims is set to match m_mesh_ptr 
      return *this;
    }
    DiscretizationXD& operator=(Eigen::VectorXd&& other){ 
      // only care about the size of the VectorXd if we already have a MeshXDPtr_t 
      if(m_mesh_ptr != nullptr)
      {
        if(other.size() != m_mesh_ptr->sizes_product()) throw std::invalid_argument("size of VectorXd must be == to m_mesh_ptr->sizes_product()"); 
      } 

      m_vals=std::move(other);       
      // assume that m_dims is set to match m_mesh_ptr 
      return *this;
    }
}; 

#endif // DiscretizationXd.hpp 