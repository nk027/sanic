
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List svd_J(
  const Eigen::Map<Eigen::MatrixXd> a,
  int type = 0,
  int precond = 0) {

  Eigen::JacobiSVD < Eigen::MatrixXd, Eigen::ColPivHouseholderQRPreconditioner> svd;
  if(precond == 1) {
    Eigen::JacobiSVD < Eigen::MatrixXd, Eigen::FullPivHouseholderQRPreconditioner > svd;
  } else if(precond == 2) {
    Eigen::JacobiSVD < Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner > svd;
  } else if(precond == 3) {
    Eigen::JacobiSVD < Eigen::MatrixXd, Eigen::NoQRPreconditioner > svd;
  } else if(precond > 3) {
    Rcpp::warning("No valid preconditioner requested - using default.");
  }

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
      Rcpp::Named("values") = svd.singularValues(),
      Rcpp::Named("U") = svd.matrixU(),
      Rcpp::Named("V") = svd.matrixV());
  }
}


// [[Rcpp::export]]
Rcpp::List svd_BDC(
  const Eigen::Map<Eigen::MatrixXd> a,
  int type = 0) {

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
      Rcpp::Named("values") = svd.singularValues(),
      Rcpp::Named("U") = svd.matrixU(),
      Rcpp::Named("V") = svd.matrixV());
  }
}
