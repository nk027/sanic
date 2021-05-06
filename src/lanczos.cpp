
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// // [[Rcpp::export]]
// Eigen::MatrixXd lanczos(
//   const Eigen::MappedSparseMatrix<double> a,
//   int iter = 0) {

//   if(!iter) {
//     iter = a.rows();
//   }

//   Eigen::VectorXd v = Eigen::VectorXd::Random(a.rows()).normalized();
//   Eigen::VectorXd w_p = a * v;
//   double alpha = w_p.adjoint() * v;
//   Eigen::VectorXd w = w_p - alpha * v;
//   double beta;

//   // This is just a skeleton
//   for(int i = 1; i <= iter; i++) {
//     beta = w.norm();
//     if(beta != 0) {
//       v = w.normalized();
//     } else {
//       v = Eigen::VectorXd::Random(a.rows()).normalized();
//     }
//     w_p = a * v;
//     alpha = w_p.adjoint() * v;
//     w = w_p - alpha * v - beta * v2;
//   }

//   return v;
// }


// [[Rcpp::export]]
Rcpp::List arnoldi(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map <Eigen::VectorXd> b,
  double tol = 0, int iter = 0) {

  if(!iter) {
    iter = a.rows();
  }
  if(!tol) {
    tol = Eigen::NumTraits<double>::dummy_precision();
  }

  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(a.rows(), iter);
  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(iter, a.rows());
  Eigen::VectorXd v(a.rows());
  q.col(0) = b.normalized();
  // Eigen::VectorXd::Random(a.rows()).normalized();

  for(int i = 0; i <= iter; i++) {
    // Rcpp::Rcout << "i" << i << "\t";
    v = a * q.col(i);
    for(int j = 0; j <= i; j++) {
      // Rcpp::Rcout << "j" << j << "\t";
      h(j, i) = q.col(j).dot(v);
      v = v - h(j, i) * q.col(j);
    }
    if(i < iter) {
      h(i + 1, i) = v.norm();
      if(h(i + 1, i) >= tol) {
        q.col(i + 1) = v / h(i + 1, i);
      } else {
        goto exit;
      }
    }
  }

  exit:
  Eigen::RealSchur <Eigen::MatrixXd> schur;
  schur.computeFromHessenberg(h, q, false);
  Eigen::MatrixXd t = schur.matrixT();

  return Rcpp::List::create(Rcpp::Named("Q") = q, Rcpp::Named("H") = h, Rcpp::Named("T") = t.diagonal());
}


// [[Rcpp::export]]
Rcpp::List hessenberg(
  const Eigen::MappedSparseMatrix<double> a) {

  Eigen::HessenbergDecomposition <Eigen::MatrixXd> hess;
  hess.compute(a);
  Eigen::MatrixXd h = hess.matrixH();
  Eigen::MatrixXd q = hess.matrixQ();

  Eigen::RealSchur <Eigen::MatrixXd> schur;
  schur.computeFromHessenberg(h, q, false);
  Eigen::MatrixXd t = schur.matrixT();

  return Rcpp::List::create(Rcpp::Named("Q") = q, Rcpp::Named("H") = h, Rcpp::Named("T") = t.diagonal());
}
