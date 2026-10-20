#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template <typename Solver, typename MatrixA, typename MatrixB>
Eigen::MatrixXd run_qr(Solver& solver, const MatrixA& a, const MatrixB& b) {
  solver.compute(a);
  return solver.solve(b);
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_HQR(
  const Eigen::Map<Eigen::MatrixXd> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int pivot = 1) {

  if(pivot == 0) {
    Eigen::HouseholderQR<Eigen::MatrixXd> solver;
    return run_qr(solver, a, b);
  }
  if(pivot == 2) {
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> solver;
    return run_qr(solver, a, b);
  }
  if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver;
  return run_qr(solver, a, b);
}

// [[Rcpp::export]]
Eigen::MatrixXd solve_SQR(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int ordering = 1) {

  if(ordering == 0) {
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::AMDOrdering<int>> solver;
    return run_qr(solver, a, b);
  }
  if(ordering == 2) {
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int>> solver;
    return run_qr(solver, a, b);
  }
  if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }
  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
  return run_qr(solver, a, b);
}
