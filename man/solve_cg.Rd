% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_iterative.R
\name{solve_cg}
\alias{solve_cg}
\title{Solve Systems of Equations Using Iterative Methods}
\usage{
solve_cg(
  a,
  b,
  x0,
  type = c("BiCGSTAB", "LSCG", "CG"),
  iter,
  tol,
  precond = 1L,
  verbose = FALSE
)
}
\arguments{
\item{a}{Square numeric matrix with the coefficients of the linear system.
Dense and sparse matrices are supported, but the format must be sparse (see
\code{\link{sparsify}}). Dense matrices are coerced automatically.}

\item{b}{Numeric vector or matrix at the right-hand side of the linear
system. If missing, 'b' is set to an identity matrix and 'a' is
inverted.}

\item{x0}{Numeric vector or matrix with an initial guess. Must be of the
same dimension as 'b'.}

\item{type}{Character scalar. Whether to use the BiCGSTAB, least squares
CG or classic CG method.}

\item{iter}{Integer scalar with the maximum number of iterations. Defaults
to the theoretical maximum, i.e. the number of columns in 'a'.}

\item{tol}{Numeric scalar with the desired tolerance. Defaults to the
machine precision.}

\item{precond}{Integer scalar indicating the type of preconditioner to be
used. Defaults to diagonal preconditioning. See the Details for further
information.}

\item{verbose}{Logical scalar. Whether to print iterations and tolerance.}
}
\value{
Solves for \eqn{x} and returns a numeric matrix with the results.
}
\description{
Iterative solvers using the Conjugate Gradient method for sparse systems of
equations \eqn{Ax = b}{Ax = b}. Three different types are available: (1)
stabilized bi-conjugate gradient (BiCGSTAB) for square matrices, (2)
conjugate gradient for rectangular least-squares (LSCG), and (3) classic
conjugate gradient (CG) for symmetric positive definite matrices.
}
\details{
Preconditioners can be set to \code{0} for no / identity preconditioning,
\code{1} (default) for Jacobi / diagonal preconditioning, or \code{2} for
incomplete factorisation. Not all schemes are available for every type:

* \code{type = "BiCGSTAB"} The default is \code{precond = 1} for diagonal
preconditioning. Set \code{precond = 0} for no preconditioning, or
\code{precond = 2} for an incomplete LUT preconditioner.
* \code{type = "LSCG"} The default is \code{precond = 1} for diagonal least
squares preconditioning. Set \code{precond = 0} for no preconditioning.
* \code{type = "CG"} The default is \code{precond = 1} for diagonal
preconditioning. Set \code{precond = 0} for no preconditioning, or
\code{precond = 2} for an incomplete Cholesky preconditioner.
}
\examples{
set.seed(42)
x <- rnorm(3)

# Solve via BiCGSTAB for square matrices
A <- matrix(rnorm(9), nrow = 3, ncol = 3)
b <- A \%*\% x
norm(solve_cg(A, b, type = "B") - x)

# Solve via LSCG for rectangular matrices
A <- matrix(rnorm(12), nrow = 4, ncol = 3)
b <- A \%*\% x
norm(solve_cg(A, b, type = "LS") - x)

# Solve via classic CG for symmetric matrices
A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
b <- A \%*\% x
norm(solve_cg(A, b, type = "CG") - x)

# The input matrix A should always be in sparse format
A <- sparsify(crossprod(matrix(rnorm(9), nrow = 3, ncol = 3)))
# The right-hand side should be a dense matrix
b <- as.matrix(A \%*\% x)

# We can check the speed of convergence and quality directly
solve_cg(A, b, verbose = TRUE)
# And provide guesses as starting value
solve_cg(A, b, x0 = x, verbose = TRUE)

}
