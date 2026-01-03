// make_diagonal.hpp
//
// takes a Eigen::SparseMatrix<double,Eigen::RowMajor> matrix of rows==1 and makes it into a diagonal matrix 
// tightly based on eigen docs for make circulant
// https://libeigen.gitlab.io/eigen/docs-5.0/TopicNewExpressionType.html
// 
// JAF 1/2/2026 

#ifndef SPARSEDIAGEXPR_H
#define SPARSEDIAGEXPR_H

#include<Eigen/Core>
#include<Eigen/Sparse>

// forward declaration -----------------------------------------------
template <class ArgType>
class SparseDiag;

// type traits =======================================================================
namespace Eigen {
namespace internal {
template <class ArgType>
struct traits<SparseDiag<ArgType> > {
  typedef Eigen::Sparse StorageKind;
  typedef Eigen::MatrixXpr XprKind;
  typedef typename ArgType::StorageIndex StorageIndex;
  typedef typename ArgType::Scalar Scalar;
  enum {
    Flags = Eigen::RowMajor,
    RowsAtCompileTime = ArgType::RowsAtCompileTime,
    ColsAtCompileTime = ArgType::ColsAtCompileTime,
    MaxRowsAtCompileTime = ArgType::MaxRowsAtCompileTime,
    MaxColsAtCompileTime = ArgType::MaxColsAtCompileTime
  };
};
}  // namespace internal
}  // namespace Eigen

// expression class ======================================================================= 
template <class ArgType>
class SparseDiag : public Eigen::SparseMatrixBase<SparseDiag<ArgType> > {
  public:
    // typedefs 
    typedef typename Eigen::internal::ref_selector<SparseDiag>::type Nested;
    typedef Eigen::Index Index;
    typedef typename Eigen::internal::ref_selector<ArgType>::type ArgTypeNested;
    
    // constructors 
    SparseDiag(const ArgType& arg_init, std::size_t repeats_init=1) : m_arg(arg_init), m_num_repeats(repeats_init) {
      static_assert( std::is_same<ArgType, MatrixStorage_t>::value, "Invalid arg type. SparseDiag must be constructed from MatrixStorage_t."); 
      if(arg_init.rows()!=1) throw std::invalid_argument("SparseDiag must be constructed from MatrixStorage_t."); 
      if(repeats_init<1) throw std::invalid_argument("SparseDiag must be constructed with num_repeats >= 1."); 
    }
    
    // member functions 
    Index rows() const { return m_num_repeats * m_arg.cols(); }
    Index cols() const { return m_num_repeats * m_arg.cols(); }

    // member data 
    ArgTypeNested m_arg;
    std::size_t m_num_repeats; 
  
};

// the evaluator =======================================================================
namespace Eigen {
namespace internal {
template <typename ArgType>
struct evaluator<SparseDiag<ArgType> > : evaluator_base<SparseDiag<ArgType> > {

  // typedefs -------------------------------------------------- 
  typedef SparseDiag<ArgType> XprType;
  typedef typename nested_eval<ArgType, XprType::ColsAtCompileTime>::type ArgTypeNested;
  typedef typename remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  typedef typename XprType::Index Index; 
  typedef typename XprType::Scalar Scalar; 

  // custom InnerIterator ----------------------------------
  struct InnerIterator{
    // constructor 
    InnerIterator(const evaluator& eval, Index row_idx)
      : m_eval(eval), m_row(row_idx)
    {
      // iterator is only valid if row_idx is in bounds and underlying row has nonzero entry 
      m_valid = (m_row < m_eval.rows()) && \
                      (m_eval.m_argImpl.coeff(0,m_row / m_eval.m_num_repeats) != Scalar(0.0)); 
    };
    // member data  
    operator bool() const { return m_valid; }; 
    void operator++(){ m_valid=false; }; 
    Index row() const { return m_row; }; 
    Index col() const { return m_row; };
    Index index() const { return m_row; };
    Scalar value() const { return m_eval.m_argImpl.coeff(0,m_row / m_eval.m_num_repeats); };
    // member data
    const evaluator& m_eval; 
    Index m_row;    
    bool m_valid; 
  }; 

  // flags ------------------------------------------------------
  enum { CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost, Flags = Eigen::RowMajor };
 
  // constructors ------------------------------------------------------ 
  evaluator(const XprType& xpr) 
    : m_argImpl(xpr.m_arg), m_rows(xpr.rows()), m_cols(xpr.cols()), m_num_repeats(xpr.m_num_repeats)
  {};
 
  // member functions ------------------------------------------------------
  Index rows() const {return m_rows; }; 
  Index cols() const {return m_cols; }; 
 
  // member data ------------------------------------------------------
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  Index m_rows; 
  Index m_cols; 
  Index m_num_repeats; 
};
}  // namespace internal
}  // namespace Eigen

// the entry point ======================================================================= 
template <class ArgType>
SparseDiag<ArgType> make_SparseDiag(const Eigen::SparseMatrixBase<ArgType>& arg, std::size_t num_repeats=1) {
  return SparseDiag<ArgType>(arg.derived(), num_repeats);
}

#endif // SparseDiagExpr.hpp