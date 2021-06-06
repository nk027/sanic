
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
bool is_symmetric_D(const Eigen::Map<Eigen::MatrixXd> x,
  double tol = 0) {

  if(!tol) {
    tol = Eigen::NumTraits<double>::epsilon();
  }

  return x.isApprox(x.transpose(), tol);
  // for(int i = 1; i <= x.rows() - 1; i++) {
  //   for(int j = 0; j < i; j++) {
  //     if(fabs(x(i, j) - x(j, i)) >= tol) {
  //       return false;
  //     }
  //   }
  // }

  // return true;
}


// [[Rcpp::export]]
bool is_symmetric_S(const Eigen::MappedSparseMatrix<double> x,
  double tol = 0) {

  if(!tol) {
    tol = Eigen::NumTraits<double>::epsilon();
  }

  return x.isApprox(x.transpose(), tol);
}
