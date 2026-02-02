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
#include<Eigen/SparseCore>
#include<Eigen/Sparse>
#include "Mesh.hpp"
#include "MeshXD.hpp"
#include "Discretization.hpp" 
#include "DiscretizationXD.hpp" 
#include "LinOpMixIn.hpp" 

namespace LinOps{

using MatrixStorage_t = Eigen::SparseMatrix<double,Eigen::RowMajor>; 

template<typename DERIVED>
class LinOpBase1D 
{
  public:
    // Type Defs -----------------------------------------------------
    typedef struct{} is_1dim_linop_tag; // to tell if a class derived from LinOpBase1D<>

    // Member Funcs  ======================================================

    // defaults given. can be overriden -------------------------------------

    // multiply the underlying expression with Discretization's underlying vecXd
    Discretization1D apply(const Discretization1D& d) const 
    {
      Eigen::VectorXd v = static_cast<const DERIVED*>(this)->GetMat() * d.values();  // calculate A*b
      Discretization1D result(std::move(v), get_weak_mesh1d()); // move A*b into result's values
      return result;
    };

    // return weak_ptr of Mesh1D pointed to. default: just convert the shared_ptr 
    Mesh1D_WPtr_t get_weak_mesh1d() const
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_weak_mesh1d(); 
        else if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_weak_mesh1d(); 
        else static_assert(false, "cannot call get_weak_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return static_cast<const DERIVED*>(this)->get_mesh1d();
      }
    }

    // Must implement! -------------------------
    // fit operator to a mesh of rectangular domain 
    void set_mesh(const Mesh1D_SPtr_t& m) 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<DERIVED&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::LStorage_t>::value) expr.Lhs().set_mesh(m); 
        if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::RStorage_t>::value) expr.Rhs().set_mesh(m); 
      }
      // Non Expression case ... 
      else{
        // // ensure we aren't resetting the mesh again
        // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
        // // on nullptr throw an error  
        // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
        // m_mesh_ptr = m; // store the mesh  
        // // perform work on m 
        static_cast<DERIVED*>(this)->set_mesh(m);
      }
    };
    
    // return Mesh1D pointed to
    Mesh1D_SPtr_t get_mesh1d() const 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_mesh1d(); 
        else if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_mesh1d(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return static_cast<const DERIVED*>(this)->get_mesh1d();
      }
    } 

}; // end LinOpBase1D<>  

template<typename DERIVED>
class LinOpBaseXD 
{
  public:
    // Type Defs -----------------------------------------------------
    typedef struct{} is_xdim_linop_tag; // to tell if a class derived from LinOpBase1D<>

    // Member Funcs  ======================================================

    // defaults given. can be overriden -------------------------------------

    // multiply the underlying expression with Discretization's underlying vecXd
    DiscretizationXD apply(const DiscretizationXD& d) const 
    {
      Eigen::VectorXd v = static_cast<const DERIVED*>(this)->GetMat() * d.values();  // calculate A*b
      DiscretizationXD result(std::move(v), get_weak_meshxd()); // move A*b into result's values
      return result;
    };

    // return weak_ptr of Mesh1D pointed to. default: just convert the shared_ptr 
    MeshXD_WPtr_t get_weak_meshxd() const
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_xdim_linop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_weak_mesh1d(); 
        else if constexpr(traits::is_xdim_linop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_weak_mesh1d(); 
        else static_assert(false, "cannot call get_weak_meshxd() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return static_cast<const DERIVED*>(this)->get_meshxd();
      }
    }

    // Must implement! -------------------------
    // fit operator to a mesh of rectangular domain 
    void set_mesh(const MeshXD_SPtr_t& m) 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<DERIVED&>(*this);  
        if constexpr(traits::is_xdim_linop_crtp<typename DERIVED::LStorage_t>::value) expr.Lhs().set_mesh(m); 
        if constexpr(traits::is_xdim_linop_crtp<typename DERIVED::RStorage_t>::value) expr.Rhs().set_mesh(m); 
      }
      // Non Expression case ... 
      else{
        // // ensure we aren't resetting the mesh again
        // if(!m_mesh_ptr.owner_before(m) && !m.owner_before(m_mesh_ptr)) return;
        // // on nullptr throw an error  
        // if(!m) throw std::runtime_error("set_mesh error: std::shared_ptr<const Mesh1D> is expried"); 
        // m_mesh_ptr = m; // store the mesh  
        // // perform work on m 
        static_cast<DERIVED*>(this)->set_mesh(m);
      }
    };
    
    // return Mesh1D pointed to
    MeshXD_SPtr_t get_meshxd() const 
    {
      // if DERIVED is an expression 
      if constexpr(traits::is_expr_crtp<DERIVED>::value){
        auto& expr = static_cast<const DERIVED&>(*this);  
        if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::LStorage_t>::value) return expr.Lhs().get_meshxd(); 
        else if constexpr(traits::is_1dim_linop_crtp<typename DERIVED::RStorage_t>::value) return expr.Rhs().get_meshxd(); 
        else static_assert(false, "cannot call get_mesh1d() on expr with no LinOpBase1D's in it!"); 
      }
      // Non Expression case ... 
      else{
        return static_cast<const DERIVED*>(this)->get_meshxd();
      }
    } 

}; // end LinOpBaseXD<>  

} // end namespace LinOps 

#endif // LinearOpBase.hpp