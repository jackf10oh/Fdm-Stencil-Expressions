// LinOpExpr.hpp
//
//
//
// JAF 12/7/2025

#ifndef LINOPEXPR_H
#define LINOPEXPR_H

#include<cstdint>
#include<type_traits>
#include<thread>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinearOpBase.hpp"
#include "LinOpTraits.hpp"

// TRAITS ======================================================================
// number of LinOpBase's stored in an expression tree 
// terminating case 
template<typename T, typename = void>
struct expr_traits
{
  static constexpr std::size_t num_leaves =
      is_linop_crtp<T>::value ? 1 : 0;
};

// recursive case 
template<typename EXPR_T>
struct expr_traits<EXPR_T, std::enable_if_t<is_expr_crtp<EXPR_T>::value>>
{
  static constexpr std::size_t num_leaves = expr_traits<typename EXPR_T::LStorage_t>::num_leaves \
                                            + expr_traits<typename EXPR_T::RStorage_t>::num_leaves; 
}; 

// Expression of L,R, & BinOp --------------------------------------------
template<typename Lhs_t, typename Rhs_t, typename BinaryOp_t>
class LinOpExpr : public LinOpBase<LinOpExpr<Lhs_t, Rhs_t, BinaryOp_t>>
{
  public:
    // types and constexpr flags
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    using RStorage_t = typename Storage_t<Rhs_t>::type;
  private:
    // member data   
    LStorage_t m_Lhs;
    RStorage_t m_Rhs;
    BinaryOp_t m_BinOp; 
    
  public:
    // Constructors 
    LinOpExpr(LStorage_t A, RStorage_t B, BinaryOp_t bin_op)
      : m_Lhs(A), m_Rhs(B), m_BinOp(bin_op)
    {};

    // destructors
    ~LinOpExpr()=default;

    // member funcs
    // getters
    LStorage_t& Lhs(){ return m_Lhs; }; 
    RStorage_t& Rhs(){ return m_Rhs; }; 
    // returns combination bin_op(A,B) of 2 stored LinOps----------------
    decltype(auto) GetMat()
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    decltype(auto) GetMat() const
    {
      return m_BinOp(m_Lhs, m_Rhs); 
    };
    // apply bin_op(A,B) * x --------------------
    // Discretization1D apply(const Discretization1D& d_arr) const 
    // {
    //   // we might constexpr branch this to L.apply(R.apply(d.values())) later ...
    //   // potentially introduces bug if 
    //   // R(.) maps boundaries -> L(.) uses boundaries 
    // };
    // sets both stored diffops to work on a mesh ------------------
    template<typename ARRAY_ITER_T>
    void set_mesh_bumper(ARRAY_ITER_T& it, MeshPtr_t m){
      constexpr bool left = is_linop_crtp<Lhs_t>::value; 
      constexpr bool right = is_linop_crtp<Rhs_t>::value;
      // if LHS is a LINOP or an EXPR 
      if constexpr(left){
        // place thread if LHS is strictly a LINOP
        using left_cleaned = std::remove_const_t<std::remove_reference_t<Lhs_t>>; 
        if constexpr(!is_expr_crtp<left_cleaned>::value){
          std::cout <<"Linop hit" <<std::endl; 
          std::cout << typeid(left_cleaned).name() << std::endl; 
          if(m!=m_Lhs.mesh()){
            std::cout <<"New mesh! creating thread..." <<std::endl; 
            std::thread t1([&](MeshPtr_t m){m_Lhs.set_mesh(m);}, m);  
            *it = std::move(t1); 
          }
          it++; 
        }
        else{
          // recursive set_mesh_bumper if LHS is strictly a EXPR
          m_Lhs.template set_mesh_bumper(it, m); 
        }
      } 
      // if RHS is a LINOP or an EXPR 
      if constexpr(right){
        // place thread if RHS is strictly a LINOP
        using right_cleaned = std::remove_const_t<std::remove_reference_t<Rhs_t>>; 
        if constexpr(!is_expr_crtp<right_cleaned>::value){
          std::cout <<"Linop hit" <<std::endl; 
          std::cout << typeid(right_cleaned).name() << std::endl; 
          if(m!=m_Rhs.mesh()){
            std::cout <<"New mesh! creating thread..." <<std::endl; 
            std::thread t2([&](MeshPtr_t m){m_Rhs.set_mesh(m);}, m);
            *it = std::move(t2); 
          }
          it++; 
        }
        else{ 
          // recursive set_mesh_bumper if RHS is strictly a EXPR 
          m_Rhs.template set_mesh_bumper(it, m); 
        }
      }
      // neither LHS or RHS is a linop or expr 
      return; 
    }
    void set_mesh(MeshPtr_t m)
    {
      constexpr std::size_t N = expr_traits<LinOpExpr<Lhs_t,Rhs_t, BinaryOp_t>>::num_leaves; 
      std::cout << N << " leaves" << std::endl; 

      std::array<std::thread, N> thread_list;
      auto it = thread_list.begin();  
      set_mesh_bumper(it, m); 
      for(std::size_t i=0; i<N; i++) {
        thread_list[i].join(); 
        // std::cout<<"joining: "<<i<<std::endl;
      } 
    }; 
    const MeshPtr_t& mesh() const 
    {
      // if LHS is from linop base. given priority over RHS 
      if constexpr(is_linop_crtp<LStorage_t>::value)
      {
        // and it has a mesh 
        const MeshPtr_t& result = m_Lhs.mesh(); 
        if(result) return result; 
      }
      // if RHS is from linop base. 
      if constexpr(is_linop_crtp<RStorage_t>::value)
      {
        // and it has a mesh 
        const MeshPtr_t& result = m_Rhs.mesh(); 
        if(result) return result; 
      }
      // any other other cases give the stored mesh in the expression. presumably a nullptr. 
      return this->m_mesh_ptr; 
    }
};

// Specialization for Unary operators --------------------------------------------
template<typename Lhs_t, typename UnaryOp_t>
class LinOpExpr<Lhs_t, void, UnaryOp_t> : public LinOpBase<LinOpExpr<Lhs_t, void, UnaryOp_t>>
{
  public:
    // types and constexpr flags
    using LStorage_t = typename Storage_t<Lhs_t>::type;
    // using RStorage_t = typename Storage_t<Rhs_t>::type; // not storing a second argument anymore 
  private:
    // member data   
    LStorage_t m_Lhs;
    // RStorage_t m_Rhs; // not storing a second argument anymore 
    UnaryOp_t m_UnarOp; 
    
  public:
    // Constructors 
    LinOpExpr(LStorage_t A, UnaryOp_t unar_op)
      : m_Lhs(A), m_UnarOp(unar_op)
    {};

    // destructors
    ~LinOpExpr()=default;

    // member funcs
    // getters
    LStorage_t& Lhs(){ return m_Lhs; }; 
    // RStorage_t& Rhs(){ return m_Rhs; }; 
    void Rhs(){}; // since nothing is stored. but still declared so that decltype(Rhs()) is still usable   
    // returns combination bin_op(A,B) of 2 stored LinOps----------------
    decltype(auto) GetMat()
    {
      // return m_BinOp(m_Lhs, m_Rhs); 
      return m_UnarOp(m_Lhs); 
    };
    decltype(auto) GetMat() const
    {
      // return m_BinOp(m_Lhs, m_Rhs); 
      return m_UnarOp(m_Lhs); 
    };
    // apply bin_op(A,B) * x --------------------
    // Discretization1D apply(const Discretization1D& d_arr) const 
    // {
    //   // we might constexpr branch this to L.apply(R.apply(d.values())) later ...
    //   // potentially introduces bug if 
    //   // R(.) maps boundaries -> L(.) uses boundaries 
    // };
    // sets both stored diffops to work on a mesh ------------------
    // sets both stored diffops to work on a mesh ------------------
    template<typename ARRAY_ITER_T>
    void set_mesh_bumper(ARRAY_ITER_T& it, MeshPtr_t m){
      constexpr bool left = is_linop_crtp<Lhs_t>::value; 
      // constexpr bool right = is_linop_crtp<Rhs_t>::value;
      // if LHS is a LINOP or an EXPR 
      if constexpr(left){
        // place thread if LHS is strictly a LINOP
        using left_cleaned = std::remove_const_t<std::remove_reference_t<Lhs_t>>; 
        if constexpr(!is_expr_crtp<left_cleaned>::value){
          std::cout <<"Linop hit" <<std::endl; 
          std::cout << typeid(left_cleaned).name() << std::endl; 
          if(m!=m_Lhs.mesh()){
            std::cout <<"New mesh! creating thread..." <<std::endl; 
            std::thread t1([&](MeshPtr_t m){m_Lhs.set_mesh(m);}, m);  
            *it = std::move(t1); 
          }
          it++; 
        }
        else{
          // recursive set_mesh_bumper if LHS is strictly a EXPR
          m_Lhs.template set_mesh_bumper(it, m); 
        }
      } 
      // if RHS is a LINOP or an EXPR 
      // if constexpr(right){
      //   // place thread if RHS is strictly a LINOP
      //   using right_cleaned = std::remove_const_t<std::remove_reference_t<Rhs_t>>; 
      //   if constexpr(!is_expr_crtp<right_cleaned>::value){
      //     std::cout <<"Linop hit: creating thread" <<std::endl; 
      //     std::cout << typeid(right_cleaned).name() << std::endl; 
      //     std::thread t2([&](MeshPtr_t m){m_Rhs.set_mesh(m);}, m);
      //     *it = std::move(t2); 
      //     it++; 
      //   }
      //   else{ 
      //     // recursive set_mesh_bumper if RHS is strictly a EXPR 
      //     m_Rhs.template set_mesh_bumper(it, m); 
      //   }
      // }
      // neither LHS or RHS is a linop or expr 
      return; 
    }
    void set_mesh(MeshPtr_t m)
    {
      constexpr std::size_t N = expr_traits<LinOpExpr<Lhs_t,void, UnaryOp_t>>::num_leaves; 
      std::cout << N << " leaves" << std::endl; 

      std::array<std::thread, N> thread_list;
      auto it = thread_list.begin();  
      set_mesh_bumper(it, m); 
      for(std::size_t i=0; i<N; i++) {
        thread_list[i].join(); 
        // std::cout<<"joining: "<<i<<std::endl;
      } 
    }; 
    // return a const ref to stored mesh of LHS or RHS. 
    const MeshPtr_t& mesh() const 
    {
      // if LHS is from linop base. given priority over RHS 
      if constexpr(is_linop_crtp<LStorage_t>::value)
      {
        // and it has a mesh 
        const MeshPtr_t& result = m_Lhs.mesh(); 
        if(result) return result; 
      }
      // // if RHS is from linop base. dont need to check RHS anymore  
      // if constexpr(is_linop_crtp<RStorage_t>::value)
      // {
      //   // and it has a mesh 
      //   const MeshPtr_t& result = m_Rhs.mesh(); 
      //   if(result) return result; 
      // }
      // any other other cases give the stored mesh in the expression. presumably a nullptr. 
      return this->m_mesh_ptr; 
    }
};

#endif // LinOpExpr.hpp