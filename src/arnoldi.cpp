
#include <RcppEigen.h>
#include "UpperHessenbergEigen.h"

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Rcpp::List lanczos_E(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
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

  if(b.size() != a.rows() || !b.allFinite() || b.norm() == 0) {
    Rcpp::stop("Starting vector 'b' must be finite, non-zero, and match the matrix dimensions.");
  }

  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(a.rows(), iter);
  q.col(0) = b.normalized();
  Eigen::VectorXd d = Eigen::VectorXd::Zero(iter);
  Eigen::VectorXd s = Eigen::VectorXd::Zero(iter - 1);
  Eigen::VectorXd v = Eigen::VectorXd(a.rows());
  unsigned int m = 1;

  v = a * q.col(0);
  d(0) = q.col(0).dot(v);
  v = v - d(0) * q.col(0);

  for(unsigned int i = 1; i < iter; i++) {
    if(orthogonalise) {
      for(unsigned int j = 0; j <= i; j++) {
        v = v - q.col(j).dot(v) * q.col(j);
      }
    }

    const double beta = v.norm();
    if(beta < tol) break;

    s(i - 1) = beta;
    q.col(i) = v / beta;
    m = i + 1;

    v = a * q.col(i);
    d(i) = q.col(i).dot(v);
    v = v - d(i) * q.col(i) - beta * q.col(i - 1);
  }

  const Eigen::VectorXd d_out = d.head(m);
  const Eigen::VectorXd s_out = s.head(m - 1);
  const Eigen::MatrixXd q_out = q.leftCols(m);


  if(!eigen) {
    return Rcpp::List::create(
    Rcpp::Named("diagonal") = d_out,
    Rcpp::Named("subdiagonal") = s_out,
    Rcpp::Named("Q") = q_out
    );
  }

  // Solve eigenproblem
  Eigen::SelfAdjointEigenSolver <Eigen::MatrixXd> tri;
  tri.computeFromTridiagonal(d_out, s_out, false);

  return Rcpp::List::create(
    Rcpp::Named("diagonal") = d_out,
    Rcpp::Named("subdiagonal") = s_out,
    Rcpp::Named("Q") = q_out,
    Rcpp::Named("values") = tri.eigenvalues()
  );
}


// [[Rcpp::export]]
Rcpp::List arnoldi_E(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map <Eigen::VectorXd> b,
  double tol = 0, unsigned int iter = 0,
  bool eigen = true) {

   if(!iter || iter > a.rows()) {
    iter = a.rows();
  }
  if(!tol) {
    tol = Eigen::NumTraits<double>::epsilon();
  }

  if(b.size() != a.rows() || !b.allFinite() || b.norm() == 0) {
    Rcpp::stop("Starting vector 'b' must be finite, non-zero, and match the matrix dimensions.");
  }

  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(a.rows(), iter);
  q.col(0) = b.normalized();
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(iter, iter);
  Eigen::VectorXd v = Eigen::VectorXd(a.rows());
  unsigned int m = 0;

  for(unsigned int i = 0; i < iter; i++) {
    v = a * q.col(i);
    for(unsigned int j = 0; j <= i; j++) {
      h(j, i) = q.col(j).dot(v);
      v = v - h(j, i) * q.col(j);
    }

    m = i + 1;
    const double beta = v.norm();
    if(beta < tol || i == iter - 1) break;
    h(i + 1, i) = beta;
    q.col(i + 1) = v / beta;
  }
  const Eigen::MatrixXd h_out = h.topLeftCorner(m, m);
  const Eigen::MatrixXd q_out = q.leftCols(m);

  if(!eigen) {
    return Rcpp::List::create(
      Rcpp::Named("H") = h_out,
      Rcpp::Named("Q") = q_out
    );
  }

  // Solve eigenproblem
  Spectra::UpperHessenbergEigen <> eig;
  eig.compute(h_out, false);

  return Rcpp::List::create(
    Rcpp::Named("H") = h_out,
    Rcpp::Named("Q") = q_out,
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
