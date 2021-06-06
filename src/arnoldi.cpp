
#include <RcppEigen.h>
#include "UpperHessenbergEigen.h"

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List lanczos_E(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map <Eigen::VectorXd> b,
  double tol = 0, unsigned int iter = 0,
  bool eigen = true,
  bool orthogonalise = true) {

  if(!iter || iter > a.rows()) {
    iter = a.rows();
  }
  if(!tol) {
    tol = Eigen::NumTraits<double>::epsilon();
  }

  // Orthonormal matrix
  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(a.rows(), iter);
  q.col(0) = b.normalized(); // Eigen::VectorXd::Random(a.rows()).normalized();
  // Diagonal and subdiagonal vectors
  Eigen::VectorXd d = Eigen::VectorXd(iter);
  Eigen::VectorXd s = Eigen::VectorXd(iter - 1);
  // Temporary vector
  Eigen::VectorXd v = Eigen::VectorXd(a.rows());

  // Initial iteration
  v = a * q.col(0);
  d(0) = q.col(0).dot(v);
  v = v - d(0) * q.col(0);

  // Iterate
  for(int i = 1; i < iter; i++) {

    if(orthogonalise) {
      for(int j = 0; j <= i; j++) { // Orthogonalise using modified Gram-Schmidt
        v = v - q.col(j).dot(v) * q.col(j);
      }
    }

    if(v.norm() < tol) { // Zero vector -- hopefully converged
      break;
    }

    // Subdiagonal coefficient
    s(i - 1) = v.norm();
    q.col(i) = v / s(i - 1);

    v = a * q.col(i);

    // Diagonal coefficient
    d(i) = q.col(i).dot(v);
    v = v - d(i) * q.col(i) - s(i - 1) * q.col(i - 1);
  }


  if(!eigen) {
    return Rcpp::List::create(
    Rcpp::Named("diagonal") = d,
    Rcpp::Named("subdiagonal") = s,
    Rcpp::Named("Q") = q
    );
  }

  // Solve eigenproblem
  Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> tri;
  tri.computeFromTridiagonal(d, s, false);

  return Rcpp::List::create(
    Rcpp::Named("diagonal") = d,
    Rcpp::Named("subdiagonal") = s,
    Rcpp::Named("Q") = q,
    Rcpp::Named("values") = tri.eigenvalues()
  );
}


// [[Rcpp::export]]
Rcpp::List arnoldi_E(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map <Eigen::VectorXd> b,
  double tol = 0, unsigned int iter = 0,
  bool eigen = true) {

   if(!iter || iter > a.rows()) {
    iter = a.rows();
  }
  if(!tol) {
    tol = Eigen::NumTraits<double>::epsilon();
  }

  // Arnoldi matrix
  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(a.rows(), iter);
  q.col(0) = b.normalized(); // Eigen::VectorXd::Random(a.rows()).normalized();
  // Hessenberg matrix
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(iter, iter);
  // Temporary vector
  Eigen::VectorXd v = Eigen::VectorXd(a.rows());

  // Iterate
  for(int i = 0; i < iter; i++) {

    v = a * q.col(i);

    for(int j = 0; j <= i; j++) { // Orthogonalise using modified Gram-Schmidt
      h(j, i) = q.col(j).dot(v);
      v = v - h(j, i) * q.col(j);
    }

    if(v.norm() < tol) { // Zero vector -- hopefully converged
      break;
    }
    if(i < iter - 1) {
      h(i + 1, i) = v.norm();
      q.col(i + 1) = v / h(i + 1, i);
    }
  }

  if(!eigen) {
    return Rcpp::List::create(
      Rcpp::Named("H") = h,
      Rcpp::Named("Q") = q
    );
  }

  // Solve eigenproblem
  Spectra::UpperHessenbergEigen <> eig;
  eig.compute(h.topLeftCorner(iter, iter), false);

  return Rcpp::List::create(
    Rcpp::Named("H") = h,
    Rcpp::Named("Q") = q,
    Rcpp::Named("values") = eig.eigenvalues()
  );
}


// // [[Rcpp::export]]
// Rcpp::List hessenberg_E(
//   const Eigen::Map<Eigen::MatrixXd> a) {

//   // Compute Hessenberg form
//   Eigen::HessenbergDecomposition <Eigen::MatrixXd> hess;
//   hess.compute(a);
//   Eigen::MatrixXd h = hess.matrixH();
//   Eigen::MatrixXd q = hess.matrixQ();

//   return Rcpp::List::create(
//     Rcpp::Named("H") = h,
//     Rcpp::Named("Q") = q
//   );
// }


// // [[Rcpp::export]]
// Rcpp::List tridiagonal_E(
//   const Eigen::Map<Eigen::MatrixXd> a) {

//   // Compute tridiagonal form
//   Eigen::Tridiagonalization <Eigen::MatrixXd> tri;
//   tri.compute(a);
//   Eigen::VectorXd d = tri.diagonal();
//   Eigen::VectorXd s = tri.subDiagonal();

//   return Rcpp::List::create(
//     Rcpp::Named("diagonal") = d,
//     Rcpp::Named("subdiagonal") = s
//   );
// }
