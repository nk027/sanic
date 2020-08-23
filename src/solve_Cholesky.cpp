
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_LLT(
  Eigen::Map<Eigen::MatrixXd> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::LLT<Eigen::MatrixXd> chol(A);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_LDLT(
  Eigen::Map<Eigen::MatrixXd> A, 
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::LDLT<Eigen::MatrixXd> chol(A);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLLT(
  Eigen::MappedSparseMatrix<double> A,
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > chol(A);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLDLT(
  Eigen::MappedSparseMatrix<double> A,
  Eigen::Map<Eigen::MatrixXd> b) {
  
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol(A);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}
