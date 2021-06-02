
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_HQR(
  const Eigen::Map<Eigen::MatrixXd> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int pivot = 1) {

  Eigen::ColPivHouseholderQR <Eigen::MatrixXd> qr;
  if(pivot == 0) {
    Eigen::HouseholderQR <Eigen::MatrixXd> qr;
  } else if(pivot == 2) {
    Eigen::FullPivHouseholderQR <Eigen::MatrixXd> qr;
  } else if(pivot > 2) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  qr.compute(a);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SQR(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int ordering = 1) {

  Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
  if(ordering == 0) {
    Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
  } else if(ordering == 2) {
    Eigen::SparseQR < Eigen::MappedSparseMatrix<double>, Eigen::NaturalOrdering<int> > solver;
  } else if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }

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
