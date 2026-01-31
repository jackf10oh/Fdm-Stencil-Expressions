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
#include<LinOps/LinOpTraits.hpp> // callable_traits<T> 

namespace LinOps{

struct DiscretizationXD
{
  private:
    // private typedefs ==================================================== 
    typedef Eigen::Stride<0,Eigen::Dynamic> Stride_t; 
    typedef Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Stride_t> StrideView_t; 

    // member data ==========================================================
    Eigen::VectorXd m_vals; // flattened array of values 
    MeshXD_WPtr_t m_mesh_ptr; 

  public:
    // Constructors / Destructors  ==========================================================
    // Default ----------------------------------------------- 
    DiscretizationXD()=default; 
    // from size of m_vals. assume just 1 dimension ----------
    DiscretizationXD(std::size_t size_init): m_mesh_ptr(), m_vals(size_init){}; 
    // from a MeshXDPtr --------------------------------------
    DiscretizationXD(const MeshXD_SPtr_t& mesh_init) 
      : m_mesh_ptr(mesh_init)
    {
      resize(mesh_init); 
    }
    // copy --------------------------------------
    DiscretizationXD(const DiscretizationXD& other)
      : m_mesh_ptr(other.m_mesh_ptr), m_vals(other.m_vals)
    {}; 

    // copy from Eigen::VectorXd
    DiscretizationXD(const Eigen::VectorXd& other, MeshXD_WPtr_t mesh_init = MeshXD_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(other)
    {}; 

    // move ---------------------------------------- 
    DiscretizationXD(DiscretizationXD&& other)
    : m_mesh_ptr(std::move(other.m_mesh_ptr)), m_vals(std::move(other.m_vals))
    {}; 
    
    // move from Eigen::VectorXD
    DiscretizationXD(Eigen::VectorXd&& other, MeshXD_WPtr_t mesh_init = MeshXD_WPtr_t{})
      : m_mesh_ptr(mesh_init), m_vals(std::move(other))
    {}; 
    // destructor 
    ~DiscretizationXD()=default; 

    // member functions ==========================================================
    // get underlying values 
    Eigen::VectorXd& values(){return m_vals; }
    const Eigen::VectorXd& values() const {return m_vals; } 

    // Give a list of Eigen::Map<>. each Map looks like a Discretization1D on a Mesh1d  
    std::vector<StrideView_t> OneDim_views(std::size_t ith_dim=0)
    {
      return m_mesh_ptr.lock()->OneDim_views(m_vals, ith_dim); 
    }

    // get underlying MeshXD_SPtr_t 
    MeshXD_SPtr_t mesh() const {return m_mesh_ptr.lock();}; 

    // number of dimensions 
    std::size_t dims() const{ return m_mesh_ptr.lock()->dims(); };
    // size of ith dimension  
    std::size_t dim_size(std::size_t i) const {return m_mesh_ptr.lock()->dim_size(i); }
    
    // product of all dimensions' sizes 
    std::size_t sizes_product() const { return m_vals.size(); } 
    // product of dimensions in [start,end)
    std::size_t sizes_middle_product(std::size_t start, std::size_t end){
      return m_mesh_ptr.lock()->sizes_middle_product(start,end); 
    }

    // store a new MeshXD_WPtr_t
    auto& set_mesh(MeshXD_WPtr_t m){ m_mesh_ptr = m; return *this; } 
    // set discretization to same size as meshxd's sizes_product
    DiscretizationXD& resize(const MeshXD_SPtr_t& m) { 
      m_mesh_ptr=m; 
      m_vals.conservativeResize(m->sizes_product()); 
      return *this; 
    }
    // set discretization to a constant
    auto& set_init(double val){ m_vals.setConstant(val); return *this; }
    // set vector to match a mesh size and set it constant 
    auto& set_init(const MeshXD_SPtr_t& m, double val){ resize(m), m_vals.setConstant(val); return *this;} 

    // set discretizations values according to callable type F
    template<
    typename F,
    typename = std::enable_if_t<
      std::is_same_v<double, typename traits::callable_traits<F>::result_type>
      >
    >
    auto& set_init(const MeshXD_SPtr_t& m, F func)
    {
      // stores mesh, resizes dims + vals 
      resize(m); 

      // set all values according to func(x0, x1, ... , xn) 

      // check there are enough dimensions to use callable type F
      static constexpr std::size_t num_args = traits::callable_traits<F>::num_args; 
      if(dims() < num_args) throw std::invalid_argument("# dims of MeshXD_SPtr_t must be >= # args in callable F"); 

      // some initializations before looping ------------------------------------------
      // how many values are in flattened array 
      std::size_t flat_end = sizes_middle_product(0,num_args);   
      // # of times to copy/paste first [0,flat_end) vales 
      std::size_t n_layers = m_vals.size() / flat_end; 
      
      // stores n args for std::apply(func,args_arr) later
      std::array<double, num_args> args_arr;  

      // holds .begin() iterators for each Mesh1D inside m_mesh_vec 
      std::array<std::vector<double>::const_iterator, num_args> domain_arr;
      for(std::size_t dim=0; dim<num_args; dim++) domain_arr[dim] = m->GetMesh(dim)->cbegin(); 

      // holds sizes_middle_product(i) for i=0,...,num_args // calculate all of cumulative_prod_arr
      std::array<std::size_t, num_args> cumulative_prod_arr; 
      for(std::size_t dim=0; dim<num_args; dim++) cumulative_prod_arr[dim] = sizes_middle_product(0, dim); 

      // iterate through flattened first layer 
      for(std::size_t flat_i=0; flat_i<flat_end; flat_i++){
        // fill args_arr 
        for(std::size_t dim=0; dim < num_args; dim++){
          std::size_t dim_i = flat_i;
          dim_i /= cumulative_prod_arr[dim]; 
          dim_i %= m->dim_size(dim); 
          args_arr[dim] = domain_arr[dim][dim_i]; 
        }
        // write func(arg_arr) to m_vals 
        m_vals[flat_i] = std::apply(func, args_arr);  
      } // end first layer 

      // if we need to copy into more layers 
      if(flat_end != m_vals.size()){
        // ither through views
        for(std::size_t ith_view=0; ith_view<flat_end; ith_view++){
          // fill in all layers with first layers value 
          for(std::size_t layer=1; layer<n_layers; layer++){
            m_vals[ith_view + flat_end*layer] = m_vals[ith_view]; 
          } // end for loop through values of ith layer 
        } // end for loop through layers  
      } // end if 

      // end of function
      return *this;
    }

    // set discretization based on current stored mesh + callable F 
    template<
    typename F,
    typename = std::enable_if_t<
      std::is_same_v<double, typename traits::callable_traits<F>::result_type>
      >
    >
    auto& set_init(F func)
    {
      auto locked = m_mesh_ptr.lock(); 
      set_init(locked, std::move(func));    
      // end of function
      return *this;   
    }
       
    // Operators ----------------------------------------------------
    DiscretizationXD& operator=(const DiscretizationXD& other) = default;
    DiscretizationXD& operator=(DiscretizationXD&& other){
      m_mesh_ptr = std::move(other.m_mesh_ptr); 
      // other.m_mesh_ptr = nullptr; 
      m_vals = std::move(other.m_vals); 
      return *this;
    }; 
}; 

} // end namespace LinOps 

#endif // DiscretizationXd.hpp