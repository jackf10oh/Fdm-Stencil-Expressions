// DiscretizationXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef DISCRETIZATIONXD_H
#define DISCRETIZATIONXD_H

#include<Eigen/Core>
#include "MeshXD.hpp"

struct DiscretizationXD
{
  private:
    // member data 
    Eigen::VectorXd m_vals; // flattened array of values 
    MeshXDPtr_t m_mesh_ptr; 
  public:
    // constructors
    DiscretizationXD()=default; 
    // destructors -----------------------------------------------------------------
    ~DiscretizationXD()=default; 
    // member functions ------------------------------------------------------------ 
    MeshXDPtr_t& mesh(){return m_mesh_ptr;}; 
    const MeshXDPtr_t& mesh() const {return m_mesh_ptr;}; 
}; 

#endif // DiscretizationXd.hpp 