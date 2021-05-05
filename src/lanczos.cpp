
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd lanczos(
  const Eigen::MappedSparseMatrix<double> a,
  int iter = 0) {

  if(!iter) {
    iter = a.rows();
  }

  Eigen::VectorXd v = Eigen::VectorXd::Random(a.rows()).normalized();
  Eigen::VectorXd w_p = a * v;
  double alpha = w_p.adjoint() * v;
  Eigen::VectorXd w = w_p - alpha * v;
  double beta;

  for(int i = 1; i <= iter; i++) {
    beta = w.norm();
    if(beta != 0) {
      v = w.normalized();
    } else {
      v = Eigen::VectorXd::Random(a.rows()).normalized();
    }
    w_p = a * v;
    alpha = w_p.adjoint() * v;
    w = w_p - alpha * v - beta * v2;
  }

  return v;
}
