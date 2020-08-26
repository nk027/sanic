
#ifndef DTNS_PRECONDITIONER_H
#define DTNS_PRECONDITIONER_H

using namespace Eigen;

template <typename _Scalar>
class DTNSPreconditioner {
  typedef _Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, 1> Vect;
public:
    typedef typename Vect::StorageIndex StorageIndex;
    enum {
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };
    
    DTNSPreconditioner() : m_isInitialized(false) {}
    
    template<typename MatType>
    explicit DTNSPreconditioner(const MatType& mat)
      : m_tns(mat.cols()) {
      compute(mat);
    }
    
    Index rows() const {return m_tns.size();}
    Index cols() const {return m_tns.size();}
    
    template<typename MatType>
    DTNSPreconditioner& analyzePattern(const MatType& ) {return *this;}
    
    template<typename MatType>
    DTNSPreconditioner& factorize(const MatType& mat) {
      // 1st order diagonal TNS - A_ii + 1
      m_tns.resize(mat.cols());
      for(int j = 0; j < mat.outerSize(); ++j) {
        typename MatType::InnerIterator it(mat, j);
        while(it && it.index() != j) {++it;}
        if(it && it.index() == j) {m_tns(j) = Scalar(1) + it.value();}
      }
      m_isInitialized = true;
      return *this;
    }
    
    template<typename MatType>
    DTNSPreconditioner& compute(const MatType& mat) {
      return factorize(mat);
    }
    
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const { 
      x = m_tns.array() * b.array();
    }
    
    template<typename Rhs> inline const Solve<DTNSPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const {
      eigen_assert(m_isInitialized && "DTNSPreconditioner is not initialized.");
      eigen_assert(m_tns.size() == b.rows()
        && "DTNSPreconditioner::solve(): invalid number of rows of the rhs b");
      return Solve<DTNSPreconditioner, Rhs>(*this, b.derived());
    }
    
    ComputationInfo info() {return Success;}
    
  protected:
    Vect m_tns;
    bool m_isInitialized;
};

#endif
