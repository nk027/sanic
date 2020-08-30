
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_BiCGSTAB(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b, Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0, bool verbose = false) {

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

  if(tol) {
    solver.setTolerance(tol);
  }
  if(iter) {
    solver.setMaxIterations(iter);
  }

  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Setup failed.");
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solving failed.");
  }

  if(verbose) {
    Rcpp::Rcout << "Iterations: " << solver.iterations() << "\n";
    Rcpp::Rcout << "Estimated error: " << solver.error() << "\n";
  }

  return x;
}


// [[Rcpp::export]]
Eigen::VectorXd solve_LSCG(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b, Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0, bool verbose = false) {

  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
  if(tol) {
    solver.setTolerance(tol);
  }
  if(iter) {
    solver.setMaxIterations(iter);
  }

  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Setup failed.");
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solving failed.");
  }

  if(verbose) {
    Rcpp::Rcout << "Iterations: " << solver.iterations() << "\n";
    Rcpp::Rcout << "Estimated error: " << solver.error() << "\n";
  }

  return x;
}


// [[Rcpp::export]]
Eigen::VectorXd solve_CG(
  Eigen::MappedSparseMatrix<double> a,
  Eigen::Map<Eigen::MatrixXd> b, Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0, int verbose = false) {

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver;
  if(tol) {
    solver.setTolerance(tol);
  }
  if(iter) {
    solver.setMaxIterations(iter);
  }

  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Setup failed.");
  }

  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Solving failed.");
  }

  if(verbose) {
    Rcpp::Rcout << "Iterations: " << solver.iterations() << "\n";
    Rcpp::Rcout << "Estimated error: " << solver.error() << "\n";
  }

  return x;
}
