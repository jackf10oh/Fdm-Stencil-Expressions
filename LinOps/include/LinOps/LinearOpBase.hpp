// LinearOpBase.hpp
//
// CTRP base for 1D differential operator L
// where L works on function discretizations across x0 , ..., xN 
//
// JAF 12/7/2025

#ifndef LINEAROPBASE_H
#define LINEAROPBASE_H

#include<cstdint>
#include<Eigen/Core>
#include "Mesh.hpp"
#include "MeshXD.hpp"
#include "Discretization.hpp" 
#include "DiscretizationXD.hpp" 
#include "LinOpMixIn.hpp" 

namespace LinOps{

using MatrixStorage_t = Eigen::SparseMatrix<double,Eigen::RowMajor>; 

template<typename DERIVED>
class LinOpBase1D : public internal::LinOpMixIn<LinOpBase1D<DERIVED>> 
{
  public:
    // Type Defs -----------------------------------------------------
    using DERIVED_T = DERIVED; // so LinOpMixIn can access grand child class

    // Member Funcs  ======================================================

    // defaults given. can be overriden -------------------------------------

    // multiply the underlying expression with Discretization's underlying vecXd
    Discretization1D apply(const Discretization1D& d) const 
    {
      Eigen::VectorXd v = internal::LinOpMixIn<LinOpBase1D<DERIVED_T>>::GetMat() * d.values();  // calculate A*b
      Discretization1D result(std::move(v), get_weak_mesh1d()); // move A*b into result's values
      return result;
    };

    // return weak_ptr of Mesh1D pointed to. default: just convert the shared_ptr 
    Mesh1D_WPtr_t get_weak_mesh1d() const
    {
      return static_cast<const DERIVED_T*>(this)->get_mesh1d();
    }

    // Must implement! -------------------------
    // fit operator to a mesh of rectangular domain 
    void set_mesh(const Mesh1D_SPtr_t& m) 
    {
      // // ensure we aren't resetting the mesh again
      // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
      // // on nullptr throw an error  
      // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
      // m_mesh_ptr = m; // store the mesh  
      // // perform work on m 
      static_cast<DERIVED_T*>(this)->set_mesh(m);
    };
    
    // return Mesh1D pointed to 
    Mesh1D_SPtr_t get_mesh1d() const 
    {
      return static_cast<const DERIVED_T*>(this)->get_mesh1d();
    } 

}; // end LinOpBase1D<>  

// template<typename DERIVED>
// class LinOpBaseXD : public internal::LinOpMixIn<LinOpBaseXD<DERIVED>> 
// {
//   public:
//     // Type Defs -----------------------------------------------------
//     using DERIVED_T = DERIVED; // so LinOpMixIn can access grand child class

//     // Member Funcs  ======================================================

//     // defaults given. can be overriden -------------------------------------

//     // multiply the underlying expression with Discretization's underlying vecXd
//     DiscretizationXD apply(const DiscretizationXD& d) const 
//     {
//       Eigen::VectorXd v = internal::LinOpMixIn<LinOpBaseXD<DERIVED_T>>::GetMat() * d.values();  // calculate A*b
//       DiscretizationXD result(std::move(v), get_weak_meshxd()); // move A*b into result's values
//       return result;
//     };

//     // return weak_ptr of Mesh1D pointed to. default: just convert the shared_ptr 
//     MeshXD_WPtr_t get_weak_meshxd() const
//     {
//       return static_cast<const DERIVED_T*>(this)->get_mesh1d();
//     }

//     // Must implement! -------------------------
//     // fit operator to a mesh of rectangular domain 
//     void set_mesh(const MeshXD_SPtr_t& m) 
//     {
//       // // ensure we aren't resetting the mesh again
//       // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
//       // // on nullptr throw an error  
//       // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
//       // m_mesh_ptr = m; // store the mesh  
//       // // perform work on m 
//       static_cast<DERIVED_T*>(this)->set_mesh(m);
//     };
    
//     // return Mesh1D pointed to 
//     MeshXD_SPtr_t get_meshxd() const 
//     {
//       return static_cast<const DERIVED_T*>(this)->get_mesh1d();
//     } 

// }; // end LinOpBaseXD<>  

} // end namespace LinOps 

#endif // LinearOpBase.hpp