
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_CPHQR(
  Eigen::Map<Eigen::MatrixXd> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(a);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SQR(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::SparseLU<Eigen::SparseMatrix<double> > qr(a);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}
