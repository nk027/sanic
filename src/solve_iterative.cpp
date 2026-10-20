#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template <typename Solver, typename MatrixA, typename MatrixB>
Rcpp::NumericMatrix run_iterative(
  Solver& solver, const MatrixA& a, const MatrixB& b, const MatrixB& x0,
  double tol, unsigned int iter, bool verbose) {

  if(tol) solver.setTolerance(tol);
  if(iter) solver.setMaxIterations(iter);
  solver.compute(a);
  if(solver.info() != Eigen::Success) {
    Rcpp::stop("Setup failed.");
  }
  Eigen::MatrixXd x = solver.solveWithGuess(b, x0);
  const bool converged = solver.info() == Eigen::Success;
  const int iterations = solver.iterations();
  const double error = solver.error();
  if(!converged) {
    Rcpp::warning("Iterative solver did not converge successfully.");
  }
  if(verbose) {
    Rcpp::Rcout << "Iterations: " << iterations << "\n";
    Rcpp::Rcout << "Estimated error: " << error << "\n";
  }
  Rcpp::NumericMatrix result = Rcpp::wrap(x);
  result.attr("converged") = converged;
  result.attr("iterations") = iterations;
  result.attr("error") = error;
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_BiCGSTAB(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, unsigned int iter = 0,
  unsigned int precond = 1,
  bool verbose = false) {

  if(precond == 0) {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;
    return run_iterative(solver, a, b, x0, tol, iter, verbose);
  }
  if(precond == 2) {
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
    return run_iterative(solver, a, b, x0, tol, iter, verbose);
  }
  if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;
  return run_iterative(solver, a, b, x0, tol, iter, verbose);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_LSCG(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, unsigned int iter = 0,
  unsigned int precond = 1,
  bool verbose = false) {

  if(precond == 0) {
    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;
    return run_iterative(solver, a, b, x0, tol, iter, verbose);
  }
  if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
  return run_iterative(solver, a, b, x0, tol, iter, verbose);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix solve_CG(
  const Eigen::Map<Eigen::SparseMatrix<double>> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, unsigned int iter = 0,
  unsigned int precond = 1,
  bool verbose = false) {

  if(precond == 0) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner> solver;
    return run_iterative(solver, a, b, x0, tol, iter, verbose);
  }
  if(precond == 2) {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;
    return run_iterative(solver, a, b, x0, tol, iter, verbose);
  }
  if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested -- using default.");
  }
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
  return run_iterative(solver, a, b, x0, tol, iter, verbose);
}
