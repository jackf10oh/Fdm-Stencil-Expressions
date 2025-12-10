// LinOpExpr.hpp
//
//
//
// JAF 12/7/2025

#ifndef LINOPEXPR_H
#define LINOPEXPR_H

#include<cstdint>
#include<eigen3/Eigen/Core>
#include<type_traits>
#include "Mesh.hpp"
#include "Discretization.hpp" 
#include "LinearOpBase.hpp"
#include "LinOpTraits.hpp"

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
    Discretization1D apply(const Discretization1D& d_arr) const 
    {
      // we might constexpr branch this to L.apply(R.apply(d.values())) later ...
      // potentially introduces bug if 
      // R(.) maps boundaries -> L(.) uses boundaries 
      Discretization1D result(d_arr.mesh()); 
      result = GetMat()*d_arr.values(); 
      return result; 
    };
    // sets both stored diffops to work on a mesh ------------------
    void set_mesh(MeshPtr_t m)
    {
      if constexpr(is_linop_crtp<Lhs_t>::value) m_Lhs.set_mesh(m);
      if constexpr(is_linop_crtp<Rhs_t>::value) m_Rhs.set_mesh(m);
    }; 
    const MeshPtr_t& mesh() const 
    {
      // if we are a scalar multiply expression
      if constexpr(is_scalar_multiply_expr<LinOpExpr>::value)
      {
        return m_Rhs().mesh();
      }
      else
      {
        return m_Lhs.mesh() ? m_Lhs.mesh() : m_Rhs.mesh();
      }
    }
};

// operator to add linops L1+L2 
template<
  typename DerivedL, 
  typename DerivedR,
  std::enable_if_t<
    std::conjunction_v<
      is_linop_crtp<DerivedL>,
      is_linop_crtp<DerivedR>
    >,
    int 
  > = 0
>
auto operator+(DerivedL&& Lhs, DerivedR&& Rhs)
{
  using LStorage_t = typename Storage_t<DerivedL>::type;
  using RStorage_t = typename Storage_t<DerivedR>::type;

  auto bin_op = [](const LStorage_t& A, const RStorage_t& B){ return A.GetMat()+B.GetMat(); };
  using Op_t = make_flagged_t<decltype(bin_op), OperatorAddition_t>;

  return LinOpExpr<DerivedL, DerivedR, Op_t>(
    std::forward<DerivedL>(Lhs),
    std::forward<DerivedR>(Rhs),
    static_cast<Op_t>(bin_op)
  );
};

// operator to difference linops L1-L2 
template<
  typename DerivedL, 
  typename DerivedR,
  std::enable_if_t<
    std::conjunction_v<
      is_linop_crtp<DerivedL>,
      is_linop_crtp<DerivedR>
    >,
    int 
  > = 0
>
auto operator-(DerivedL&& Lhs, DerivedR&& Rhs)
{
  using LStorage_t = typename Storage_t<DerivedL>::type;
  using RStorage_t = typename Storage_t<DerivedR>::type;

  auto bin_op = [](const LStorage_t& A, const RStorage_t& B){ return A.GetMat()-B.GetMat(); };
  using Op_t = make_flagged_t<decltype(bin_op), OperatorAddition_t>;

  return LinOpExpr<DerivedL, DerivedR, Op_t>(
    std::forward<DerivedL>(Lhs),
    std::forward<DerivedR>(Rhs),
    static_cast<Op_t>(bin_op)
  );
};

// operator for scalar multiplication c*L
template<
  typename Scalar_t, 
  typename Derived,
  std::enable_if_t<
    std::conjunction_v<
      std::is_convertible<Scalar_t,double>,
      is_linop_crtp<Derived>
    >,
    int 
  > = 0
>
auto operator*(Scalar_t c, Derived&& L)
{
  using LStorage_t = double;
  using RStorage_t = typename Storage_t<Derived>::type;

  auto bin_op = [](const LStorage_t& scalar, const RStorage_t& Mat){ return scalar*Mat.GetMat(); };
  using Op_t = make_flagged_t<decltype(bin_op), ScalarMultiply_t>;

  return LinOpExpr<double, Derived, Op_t>(
    static_cast<double>(c),
    std::forward<Derived>(L),
    static_cast<Op_t>(bin_op)
  );
};

template<
  typename Derived, 
  std::enable_if_t<
    is_linop_crtp<Derived>::value,
    int 
  > = 0
>
auto operator-(Derived&& L)
{
  return (-1.0) * std::forward<Derived>(L); 
}

#endif // LinOpExpr.hpp