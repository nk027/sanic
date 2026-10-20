#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template <typename Solver, typename MatrixA, typename MatrixB>
Eigen::MatrixXd run_lu(Solver& solver, const MatrixA& a, const MatrixB& b) {
  solver.compute(a);
  Eigen::MatrixXd x = solver.solve(b);
  if(!x.allFinite()) {
    Rcpp::stop("Solving failed.");
  }
  return x;
}

template <typename Solver, typename MatrixA, typename MatrixB>
Eigen::MatrixXd run_sparse_lu(Solver& solver, const MatrixA& a, const MatrixB& b) {
  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed.");
  }
  Eigen::MatrixXd x = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solving failed.");
  }
  return x;
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_LU(
  const Eigen::Map <Eigen::MatrixXd> a,
  const Eigen::Map <Eigen::MatrixXd> b,
  unsigned int pivot = 1) {

  if(pivot == 0) {
    Rcpp::warning("No unpivoted scheme available -- using default.");
  } else if(pivot == 2) {
    Eigen::FullPivLU<Eigen::MatrixXd> solver;
    solver.compute(a);
    if(solver.rank() < a.cols()) {
      Rcpp::stop("Matrix is singular to working precision.");
    }
    Eigen::MatrixXd x = solver.solve(b);
    if(!x.allFinite()) Rcpp::stop("Solving failed.");
    return x;
  } else if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }
  Eigen::PartialPivLU<Eigen::MatrixXd> solver;
  return run_lu(solver, a, b);
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_SLU(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int ordering = 1) {

  if(ordering == 0) {
    Eigen::SparseLU<Eigen::Map<Eigen::SparseMatrix<double>>, Eigen::AMDOrdering<int>> solver;
    return run_sparse_lu(solver, a, b);
  }
  if(ordering == 2) {
    Eigen::SparseLU<Eigen::Map<Eigen::SparseMatrix<double>>, Eigen::NaturalOrdering<int>> solver;
    return run_sparse_lu(solver, a, b);
  }
  if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }
  Eigen::SparseLU<Eigen::Map<Eigen::SparseMatrix<double>>, Eigen::COLAMDOrdering<int>> solver;
  return run_sparse_lu(solver, a, b);
}
