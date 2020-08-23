
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_HQR(
  Eigen::Map<Eigen::MatrixXd> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_CPHQR(
  Eigen::Map<Eigen::MatrixXd> A,
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(A);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SQR(
  Eigen::MappedSparseMatrix<double> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::SparseLU<Eigen::SparseMatrix<double> > qr(A);
  Eigen::MatrixXd x = qr.solve(b);

  return x;
}
