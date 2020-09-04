
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_LDLT(
  Eigen::Map<Eigen::MatrixXd> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::LDLT<Eigen::MatrixXd> chol(a);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLDLT(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b) {

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol(a);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}
