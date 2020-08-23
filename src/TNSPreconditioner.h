
#ifndef TNS_PRECONDITIONER_H
#define TNS_PRECONDITIONER_H

using namespace Eigen;

template <typename _Scalar>
class TNSPreconditioner {
  typedef _Scalar Scalar;
  typedef Matrix<Scalar, Dynamic, Dynamic> Matr;
  typedef Matrix<Scalar, Dynamic, 1> Vect;
public:
    typedef typename Vect::StorageIndex StorageIndex;
    enum {
      RowsAtCompiletime = Dynamic,
      ColsAtCompileTime = Dynamic,
      MaxColsAtCompileTime = Dynamic
    };
    
    TNSPreconditioner() : m_isInitialized(false) {}
    
    template<typename MatType>
    explicit TNSPreconditioner(const MatType& mat)
      : m_tns(mat.cols()) { // Adapt this if going for a matrix
      compute(mat);
    }
    
    Index rows() const {return m_tns.size();}
    Index cols() const {return m_tns.size();}
    
    template<typename MatType>
    TNSPreconditioner& analyzePattern(const MatType& ) {return *this;}
    
    template<typename MatType>
    TNSPreconditioner& factorize(const MatType& mat) {
      // TNS - Requires m_tns as matrix
      // 1st order TNS - I + A
      // m_tns = mat;
      // m_tns.diagonal().array() += 1;
      // 2nd order TNS - I + A + A^2
      // m_tns += m_tns * m_tns;
      
      // Diagonal TNS - Requires m_tns as vector
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
    TNSPreconditioner& compute(const MatType& mat) {
      return factorize(mat);
    }
    
    template<typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const { 
      // x = m_tns * b; // matrix
      x = m_tns.array() * b.array(); // vector
    }
    
    template<typename Rhs> inline const Solve<TNSPreconditioner, Rhs>
    solve(const MatrixBase<Rhs>& b) const {
      eigen_assert(m_isInitialized && "TNSPreconditioner is not initialized.");
      eigen_assert(m_tns.size() == b.rows()
        && "TNSPreconditioner::solve(): invalid number of rows of the rhs b");
      return Solve<TNSPreconditioner, Rhs>(*this, b.derived());
    }
    
    ComputationInfo info() {return Success;}
    
  protected:
    Vect m_tns;
    // Matr m_tns;
    bool m_isInitialized;
};

#endif
