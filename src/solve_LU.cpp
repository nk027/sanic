
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_LU(
  const Eigen::Map <Eigen::MatrixXd> a,
  const Eigen::Map <Eigen::MatrixXd> b,
  int pivot = 0) {

  Eigen::PartialPivLU <Eigen::MatrixXd> lu;
  if(pivot == 1) {
    Eigen::FullPivLU <Eigen::MatrixXd> lu;
  } else if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested - using default.");
  }

  lu.compute(a);
  Eigen::MatrixXd x = lu.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLU(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  int ordering = 0) {

  Eigen::SparseLU < Eigen::MappedSparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
  if(ordering == 1) {
    Eigen::SparseLU < Eigen::MappedSparseMatrix<double>, Eigen::AMDOrdering<int> > solver;
  } else if(ordering == 2) {
    Eigen::SparseLU < Eigen::MappedSparseMatrix<double>, Eigen::NaturalOrdering<int> > solver;
  } else if(ordering > 2) {
    Rcpp::warning("No valid ordering requested - using default.");
  }

  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Decomposition failed.");
  }
  Eigen::MatrixXd x = solver.solve(b);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solving failed.");
  }

  return x;
}
