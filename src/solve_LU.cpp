
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_PPLU(
  Eigen::Map<Eigen::MatrixXd> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
    
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
  Eigen::MatrixXd x = lu.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLU(
  Eigen::MappedSparseMatrix<double> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
    
  Eigen::SparseLU<Eigen::SparseMatrix<double> > lu(A);
  Eigen::MatrixXd x = lu.solve(b);

  return x;
}
