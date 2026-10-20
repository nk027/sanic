#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

namespace {

typedef Eigen::Ref<const Eigen::MatrixXd> DenseRhs;

// Dense LU does not expose info(); full pivoting additionally reveals rank.
template <typename Solver>
void check_lu(const Solver& solver) {
  if(solver.info() != Eigen::Success) Rcpp::stop("LU decomposition or solve failed.");
}

void check_lu(const Eigen::PartialPivLU<Eigen::MatrixXd>&) {}

void check_lu(const Eigen::FullPivLU<Eigen::MatrixXd>& solver) {
  if(solver.rank() < solver.cols()) {
    Rcpp::stop("Matrix is singular to working precision.");
  }
}

struct LUFactor {
  const Eigen::Index rows;
  explicit LUFactor(Eigen::Index n) : rows(n) {}
  virtual ~LUFactor() {}
  virtual Eigen::MatrixXd solve(const DenseRhs& b) = 0;
};

template <typename Solver>
struct TypedLUFactor : LUFactor {
  Solver solver;

  template <typename MatrixA>
  explicit TypedLUFactor(const MatrixA& a) : LUFactor(a.rows()) {
    solver.compute(a);
    check_lu(solver);
  }

  Eigen::MatrixXd solve(const DenseRhs& b) {
    if(b.rows() != rows) {
      Rcpp::stop("The factorization and 'b' must have the same number of rows.");
    }
    if(!b.allFinite()) Rcpp::stop("Matrix 'b' must contain only finite values.");
    // Eigen's sparse triangular solve can access column zero for an empty RHS.
    if(b.cols() == 0) return Eigen::MatrixXd(rows, 0);
    Eigen::MatrixXd x = solver.solve(b);
    check_lu(solver);
    if(!x.allFinite()) Rcpp::stop("Solving failed.");
    return x;
  }
};

SEXP lu_tag() {
  return Rf_install("sanic_lu_factor");
}

// Only retained factors need heap allocation and an external pointer.
struct RetainLU {
  typedef SEXP Result;

  template <typename Solver, typename MatrixA>
  Result run(const MatrixA& a) const {
    Rcpp::XPtr<LUFactor> ptr(new TypedLUFactor<Solver>(a), true, lu_tag());
    ptr.attr("class") = "sanic_lu";
    return ptr;
  }
};

struct SolveLU {
  typedef Eigen::MatrixXd Result;
  const DenseRhs b;
  explicit SolveLU(const DenseRhs& rhs) : b(rhs) {}

  template <typename Solver, typename MatrixA>
  Result run(const MatrixA& a) const {
    TypedLUFactor<Solver> factor(a);
    return factor.solve(b);
  }
};

template <typename Action>
typename Action::Result dispatch_lu(const Eigen::Map<Eigen::MatrixXd>& a,
                                  unsigned int pivot, const Action& action) {
  if(a.rows() == 0 || a.rows() != a.cols()) {
    Rcpp::stop("Please provide a nonempty square matrix 'a'.");
  }
  if(!a.allFinite()) Rcpp::stop("Matrix 'a' must contain only finite values.");
  if(pivot == 0) {
    Rcpp::warning("No unpivoted scheme available -- using default.");
  } else if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }
  if(pivot == 2) {
    return action.template run<Eigen::FullPivLU<Eigen::MatrixXd> >(a);
  }
  return action.template run<Eigen::PartialPivLU<Eigen::MatrixXd> >(a);
}

template <typename Action>
typename Action::Result dispatch_lu(const Eigen::Map<Eigen::SparseMatrix<double> >& a,
                                  unsigned int ordering, const Action& action) {
  if(a.rows() == 0 || a.rows() != a.cols()) {
    Rcpp::stop("Please provide a nonempty square matrix 'a'.");
  }
  for(int i = 0; i < a.nonZeros(); ++i) {
    if(!std::isfinite(a.valuePtr()[i])) {
      Rcpp::stop("Matrix 'a' must contain only finite values.");
    }
  }
  // SparseLU owns its working matrix/factors even with a mapped input type.
  // Accept the R mapping directly rather than copying it into a temporary.
  typedef Eigen::Map<Eigen::SparseMatrix<double> > Sparse;
  if(ordering == 0) {
    return action.template run<Eigen::SparseLU<Sparse, Eigen::AMDOrdering<int> > >(a);
  }
  if(ordering == 2) {
    return action.template run<Eigen::SparseLU<Sparse, Eigen::NaturalOrdering<int> > >(a);
  }
  if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }
  return action.template run<Eigen::SparseLU<Sparse, Eigen::COLAMDOrdering<int> > >(a);
}

} // namespace

// [[Rcpp::export]]
Eigen::MatrixXd solve_LU(const Eigen::Map<Eigen::MatrixXd> a,
                         const Eigen::Map<Eigen::MatrixXd> b,
                         unsigned int pivot = 1) {
  return dispatch_lu(a, pivot, SolveLU(b));
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_SLU(const Eigen::Map<Eigen::SparseMatrix<double> > a,
                          const Eigen::Map<Eigen::MatrixXd> b,
                          unsigned int ordering = 1) {
  return dispatch_lu(a, ordering, SolveLU(b));
}

// [[Rcpp::export]]
SEXP factorize_LU(const Eigen::Map<Eigen::MatrixXd> a, unsigned int pivot) {
  return dispatch_lu(a, pivot, RetainLU());
}

// [[Rcpp::export]]
SEXP factorize_SLU(const Eigen::Map<Eigen::SparseMatrix<double> > a,
                   unsigned int ordering) {
  return dispatch_lu(a, ordering, RetainLU());
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_LU_factor(SEXP factor, Rcpp::Nullable<Rcpp::NumericMatrix> b) {
  if(TYPEOF(factor) != EXTPTRSXP || R_ExternalPtrTag(factor) != lu_tag()) {
    Rcpp::stop("Invalid LU factorization.");
  }
  Rcpp::XPtr<LUFactor> ptr(factor);
  if(!ptr) {
    Rcpp::stop("LU factorization is no longer valid; recreate it with factorize_lu().");
  }
  if(b.isNull()) {
    return ptr->solve(Eigen::MatrixXd::Identity(ptr->rows, ptr->rows));
  }
  Rcpp::NumericMatrix rhs(b);
  Eigen::Map<Eigen::MatrixXd> mapped(rhs.begin(), rhs.nrow(), rhs.ncol());
  return ptr->solve(mapped);
}
