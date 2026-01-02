// make_diagonal.hpp
//
// takes a Eigen::SparseMatrix<double,Eigen::RowMajor> matrix of rows==1 and makes it into a diagonal matrix 
// tightly based on eigen docs for make circulant
// https://libeigen.gitlab.io/eigen/docs-5.0/TopicNewExpressionType.html
// 
// JAF 1/2/2026 

#include<Eigen/Core>
#include<Eigen/Sparse>

// forward declaration -----------------------------------------------
template <class ArgType>
class Circulant;

// type traits -----------------------------------------------------
namespace Eigen {
namespace internal {
template <class ArgType>
struct traits<Circulant<ArgType> > {
  typedef Eigen::Dense StorageKind;
  typedef Eigen::MatrixXpr XprKind;
  typedef typename ArgType::StorageIndex StorageIndex;
  typedef typename ArgType::Scalar Scalar;
  enum {
    Flags = Eigen::ColMajor,
    RowsAtCompileTime = ArgType::RowsAtCompileTime,
    ColsAtCompileTime = ArgType::RowsAtCompileTime,
    MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = ArgType::MaxRowsAtCompileTime
  };
};
}  // namespace internal
}  // namespace Eigen

// expression class --------------------------------------------------------------- 
template <class ArgType>
class Circulant : public Eigen::MatrixBase<Circulant<ArgType> > {
 public:
  Circulant(const ArgType& arg) : m_arg(arg) {
    EIGEN_STATIC_ASSERT(ArgType::ColsAtCompileTime == 1, YOU_TRIED_CALLING_A_VECTOR_METHOD_ON_A_MATRIX);
  }
 
  typedef typename Eigen::internal::ref_selector<Circulant>::type Nested;
 
  typedef Eigen::Index Index;
  Index rows() const { return m_arg.rows(); }
  Index cols() const { return m_arg.rows(); }
 
  typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
  ArgTypeNested m_arg;
};

// the evaluator ----------------------------------------------------------------------- 
namespace Eigen {
namespace internal {
template <typename ArgType>
struct evaluator<Circulant<ArgType> > : evaluator_base<Circulant<ArgType> > {
  typedef Circulant<ArgType> XprType;
  typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
  typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
 
  enum { CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost, Flags = Eigen::ColMajor };
 
  evaluator(const XprType& xpr) : m_argImpl(xpr.m_arg), m_rows(xpr.rows()) {}
 
  CoeffReturnType coeff(Index row, Index col) const {
    Index index = row - col;
    if (index < 0) index += m_rows;
    return m_argImpl.coeff(index);
  }
 
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  const Index m_rows;
};
}  // namespace internal
}  // namespace Eigen

// the entry point --------------------------------------------------------------------------- 
template <class ArgType>
Circulant<ArgType> makeCirculant(const Eigen::MatrixBase<ArgType>& arg) {
  return Circulant<ArgType>(arg.derived());
}