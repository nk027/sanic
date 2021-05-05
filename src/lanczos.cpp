
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

  // This is just a skeleton
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


// [[Rcpp::export]]
Eigen::MatrixXd arnoldi(
  const Eigen::MappedSparseMatrix<double> a,
  double tol = 0, int iter = 0) {

  if(!iter) {
    iter = a.rows();
  }
  if(!tol) {
    tol = Eigen::NumTraits<double>::dummy_precision();
  }

  Eigen::MatrixXd q;
  Eigen::MatrixXd h;
  Eigen::MatrixXd v;
  q.col(0) = Eigen::VectorXd::Random(a.rows()).normalized();

  for(int i = 0; i <= iter; i++) {
    v = a * q.col(i);
    for(int j = 0; j <= i, j++) {
      h(i, j) = q.col(i).adjoint() * z;
      v = z - h(i, j) * q.col(i);
    }
    h(i + 1, i) = v.norm();
    if(h(i + 1, i) <= tol) {
      return q;
    } else {
      q.col(i + 1) = v / h(i + 1, i);
    }
  }

  return q;
}

