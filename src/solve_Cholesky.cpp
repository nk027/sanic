
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
Eigen::MatrixXd solve_LL(
  const Eigen::Map <Eigen::MatrixXd> a,
  const Eigen::Map <Eigen::MatrixXd> b,
  unsigned int pivot = 1) {

  Eigen::LDLT <Eigen::MatrixXd> chol;
  if(pivot == 0) {
    Eigen::LLT <Eigen::MatrixXd> chol;
  } else if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
  }

  chol.compute(a);
  Eigen::MatrixXd x = chol.solve(b);

  return x;
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_SLL(
  const Eigen::MappedSparseMatrix<double> a,
  const Eigen::Map<Eigen::MatrixXd> b,
  unsigned int pivot = 1, unsigned int ordering = 0) {

  Eigen::SimplicialLDLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > solver;
  if(ordering == 1) {
    Rcpp::warning("No COLAMD ordering available -- using default.");
  } else if(ordering == 2) {
    Eigen::SimplicialLDLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > solver;
  } else if(ordering > 2) {
    Rcpp::warning("No valid ordering requested -- using default.");
  }

  if(pivot == 0) {
    Eigen::SimplicialLLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::AMDOrdering<int> > solver;
    if(ordering == 1) {
        Rcpp::warning("No COLAMD ordering available -- using default.");
    } else if(ordering == 2) {
      Eigen::SimplicialLLT < Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::NaturalOrdering<int> > solver;
    } else if(ordering > 2) {
      Rcpp::warning("No valid ordering requested -- using default.");
    }
  } else if(pivot > 1) {
    Rcpp::warning("No valid pivoting scheme requested -- using default.");
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
