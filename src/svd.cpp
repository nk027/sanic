
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List svd_J(
  const Eigen::Map<Eigen::MatrixXd> a,
  unsigned int type = 0) {

  Eigen::JacobiSVD < Eigen::MatrixXd, Eigen::ColPivHouseholderQRPreconditioner> svd;

  if(type == 0) {
    svd.compute(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
  } else if(type == 1) {
    svd.compute(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
  } else if(type == 2) {
    svd.compute(a);
  } else if(type > 2) {
    Rcpp::stop("No valid type requested.");
  }

  if(type == 2) {
    return Rcpp::List::create(
      Rcpp::Named("values") = svd.singularValues());
  } else {
    return Rcpp::List::create(
      Rcpp::Named("d") = svd.singularValues(),
      Rcpp::Named("u") = svd.matrixU(),
      Rcpp::Named("v") = svd.matrixV());
  }
}


// [[Rcpp::export]]
Rcpp::List svd_BDC(
  const Eigen::Map<Eigen::MatrixXd> a,
  unsigned int type = 0) {

  Eigen::BDCSVD <Eigen::MatrixXd> svd;

  if(type == 0) {
    svd.compute(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
  } else if(type == 1) {
    svd.compute(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
  } else if(type == 2) {
    svd.compute(a);
  } else if(type > 2) {
    Rcpp::stop("No valid type requested.");
  }

  if(type == 2) {
    return Rcpp::List::create(
      Rcpp::Named("values") = svd.singularValues());
  } else {
    return Rcpp::List::create(
      Rcpp::Named("d") = svd.singularValues(),
      Rcpp::Named("u") = svd.matrixU(),
      Rcpp::Named("v") = svd.matrixV());
  }
}
