// DiscretizationXD.hpp
//
//
//
// JAF 12/26/2025

#ifndef DISCRETIZATIONXD_H
#define DISCRETIZATIONXD_H

#include<Eigen/Core>

struct DiscretizationXD
{
  private:
    // member data 
    Eigen::VectorXd m_vals; // flattened array of values 
  public:
    // constructors
    DiscretizationXD()=default; 
    // destructors -----------------------------------------------------------------
    ~DiscretizationXD()=default; 
    // member functions ------------------------------------------------------------ 
}; 

#endif // DiscretizationXd.hpp 