#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template <typename Solver, typename MatrixA, typename MatrixB>
Eigen::MatrixXd run_cholesky(Solver& solver, const MatrixA& a, const MatrixB& b) {
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
Eigen::MatrixXd solve_LL(
  const Eigen::Map <Eigen::MatrixXd> a,
  const Eigen::Map <Eigen::MatrixXd> b,
  unsigned int pivot = 1) {

  if(pivot == 0) {
    Eigen::LLT<Eigen::MatrixXd> solver;
    return run_cholesky(solver, a, b);
  }
  if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }
  Eigen::LDLT<Eigen::MatrixXd> solver;
  return run_cholesky(solver, a, b);
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_SLL(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int pivot = 1, unsigned int ordering = 0) {

  if(ordering == 1) {
    Rcpp::warning("No COLAMD ordering available -- using default.");
  } else if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }
  const bool natural = ordering == 2;

  if(pivot == 0) {
    if(natural) {
      Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> solver;
      return run_cholesky(solver, a, b);
    }
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> solver;
    return run_cholesky(solver, a, b);
  }
  if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }
  if(natural) {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int>> solver;
    return run_cholesky(solver, a, b);
  }
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int>> solver;
  return run_cholesky(solver, a, b);
}
