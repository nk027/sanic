
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_BiCGSTAB(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0,
  int precond = 0,
  bool verbose = false) {

  Eigen::BiCGSTAB < Eigen::SparseMatrix<double> > solver;
  if(precond == 1) {
    Eigen::BiCGSTAB < Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
  } else if(precond == 2) {
    Eigen::BiCGSTAB < Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > solver;
  } else if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested - using default.");
  }

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
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0,
  int precond = 0,
  bool verbose = false) {

  Eigen::LeastSquaresConjugateGradient < Eigen::SparseMatrix<double> > solver;
  if(precond == 1) {
    Eigen::LeastSquaresConjugateGradient < Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > solver;
  } else if(precond > 1) {
    Rcpp::warning("No valid preconditioner requested - using default.");
  }

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
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  const Eigen::Map<Eigen::MatrixXd> x0,
  double tol = 0, int iter = 0,
  int precond = 0,
  bool verbose = false) {

  Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper > solver;
  if(precond == 1) {
    Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IncompleteCholesky<double> > solver;
  } else if(precond == 2) {
    Eigen::ConjugateGradient < Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner > solver;
  } else if(precond > 2) {
    Rcpp::warning("No valid preconditioner requested - using default.");
  }

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
