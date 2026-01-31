// IOp.hpp 
//
// the identity operator
// 
// JAF 12/7/2025 

#ifndef IDENTITYOP_H
#define IDENTITYOP_H

#include<variant> 

#include "../LinearOpBase.hpp"

namespace LinOps{

class IOp : public LinOpMixIn<IOp>, public LinOpBase1D<IOp>, public LinOpBaseXD<IOp>
{
  private:
    std::variant<Mesh1D_WPtr_t, MeshXD_WPtr_t> m_mesh_ptr; 
    MatrixStorage_t m_Mat; 
  public: 
    // Constructors --------------------------
    IOp(std::size_t s_init=0)
    {
      resize(s_init); 
    } 
    IOp(const Mesh1D_SPtr_t& m)
    {
      set_mesh(m);
    } 
    IOp(const MeshXD_SPtr_t& m)
    {
      set_mesh(m);
    } 
    // Destructors ----------------------------
    ~IOp()=default;

    // Member Funcs  ======================================================

    // matrix getters 
    MatrixStorage_t& GetMat() { return m_Mat; };
    const MatrixStorage_t& GetMat() const { return m_Mat; };

    // Identity just returns inputs as outputs
    Discretization1D apply(const Discretization1D& d_arr) const { return d_arr; } 
    DiscretizationXD apply(const DiscretizationXD& d_arr) const { return d_arr; } 

    // return weak_ptr of Mesh1D pointed to
    Mesh1D_WPtr_t get_weak_mesh1d() const { 
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr); 
      return Mesh1D_WPtr_t{}; 
    }

    // return weak_ptr of MeshXD pointed to
    MeshXD_WPtr_t get_weak_meshxd() const { 
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr); 
      return MeshXD_WPtr_t{}; 
    }

    // return Mesh1D pointed to 
    Mesh1D_SPtr_t get_mesh1d() const {
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr)) return std::get<Mesh1D_WPtr_t>(m_mesh_ptr).lock(); 
      return Mesh1D_SPtr_t{};
    } 

    // return MeshXD pointed to 
    MeshXD_SPtr_t get_meshxd() const {
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr)) return std::get<MeshXD_WPtr_t>(m_mesh_ptr).lock(); 
      return MeshXD_SPtr_t{};
    } 

    // fit operator to a domain mesh 
    void set_mesh(const Mesh1D_SPtr_t& m) 
    {
      // ensure we aren't resetting the mesh again
      if(std::holds_alternative<Mesh1D_WPtr_t>(m_mesh_ptr))
      {
        auto my_m = std::get<Mesh1D_WPtr_t>(m_mesh_ptr); 
        if(!my_m.owner_before(m) && !m.owner_before(my_m)) return;
      }
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("IOp.set_mesh(m) error: std::shared_ptr<const Mesh1D> is expried"); 
      m_mesh_ptr.emplace<Mesh1D_WPtr_t>(m); // store the mesh  
      // perform work on m 
      // set m_Mat to I(s,s) 
      std::size_t s = m->size(); 
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };

    // fit operator to a domain mesh 
    void set_mesh(const MeshXD_SPtr_t& m) 
    {
      // ensure we aren't resetting the mesh again
      if(std::holds_alternative<MeshXD_WPtr_t>(m_mesh_ptr))
      {
        auto my_m = std::get<MeshXD_WPtr_t>(m_mesh_ptr); 
        if(!my_m.owner_before(m) && !m.owner_before(my_m)) return;
      }
      // on nullptr throw an error  
      if(!m) throw std::runtime_error("IOp.set_mesh(m) error: std::shared_ptr<const Mesh1D> is expried"); 
      m_mesh_ptr.emplace<MeshXD_WPtr_t>(m); // store the mesh  
      // perform work on m 
      // set m_Mat to I(s,s) 
      std::size_t s = m->sizes_product(); 
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };
    
    // resize Identity to s; 
    void resize(std::size_t s=0)
    {
      m_Mat.resize(s, s); 
      m_Mat.setIdentity(); 
    };

}; // End IOp 

} // end namespace LinOps 

#endif