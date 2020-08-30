
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_PPLU(
  Eigen::Map<Eigen::MatrixXd> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::PartialPivLU<Eigen::MatrixXd> lu(a);
  Eigen::MatrixXd x = lu.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLU(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::SparseLU<Eigen::SparseMatrix<double> > lu(a);
  Eigen::MatrixXd x = lu.solve(b);

  return x;
}
