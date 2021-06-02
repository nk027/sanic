
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List eigen_SA(
  const Eigen::Map<Eigen::MatrixXd> a,
  bool vectors = true) {

  Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> eig;

  if(vectors) {
    eig.compute(a, Eigen::ComputeEigenvectors);
  } else {
    eig.compute(a, Eigen::EigenvaluesOnly);
  }

  if(vectors) {
    return Rcpp::List::create(
      Rcpp::Named("values") = eig.eigenvalues().real(),
      Rcpp::Named("vectors") = eig.eigenvectors());
  } else {
    return Rcpp::List::create(
      Rcpp::Named("values") = eig.eigenvalues().real());
  }
}


// [[Rcpp::export]]
Rcpp::List eigen_SQ(
  const Eigen::Map<Eigen::MatrixXd> a,
  bool vectors = true) {

  Eigen::EigenSolver < Eigen::MatrixXd > eig;

  eig.compute(a, vectors);

  if(vectors) {
    return Rcpp::List::create(
      Rcpp::Named("values") = eig.eigenvalues(),
      Rcpp::Named("vectors") = eig.eigenvectors());
  } else {
    return Rcpp::List::create(
      Rcpp::Named("values") = eig.eigenvalues());
  }
}
