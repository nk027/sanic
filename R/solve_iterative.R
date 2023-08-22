
#' Solve Systems of Equations Using Iterative Methods
#'
#' Iterative solvers using the Conjugate Gradient method for sparse systems of
#' equations \eqn{Ax = b}{Ax = b}. Three different types are available: (1)
#' stabilized bi-conjugate gradient (BiCGSTAB) for square matrices, (2)
#' conjugate gradient for rectangular least-squares (LSCG), and (3) classic
#' conjugate gradient (CG) for symmetric positive definite matrices.
#'
#' @details
#' Preconditioners can be set to \code{0} for no / identity preconditioning,
#' \code{1} (default) for Jacobi / diagonal preconditioning, or \code{2} for
#' incomplete factorisation. Not all schemes are available for every type:
#'
#' * \code{type = "BiCGSTAB"} The default is \code{precond = 1} for diagonal
#' preconditioning. Set \code{precond = 0} for no preconditioning, or
#' \code{precond = 2} for an incomplete LUT preconditioner.
#' * \code{type = "LSCG"} The default is \code{precond = 1} for diagonal least
#' squares preconditioning. Set \code{precond = 0} for no preconditioning.
#' * \code{type = "CG"} The default is \code{precond = 1} for diagonal
#' preconditioning. Set \code{precond = 0} for no preconditioning, or
#' \code{precond = 2} for an incomplete Cholesky preconditioner.
#'
#' @inheritParams solve_chol
#' @param a Square numeric matrix with the coefficients of the linear system.
#' Dense and sparse matrices are supported, but the format must be sparse (see
#' \code{\link{sparsify}}). Dense matrices are coerced automatically.
#' @param type Character scalar. Whether to use the BiCGSTAB, least squares
#' CG or classic CG method.
#' @param x0 Numeric vector or matrix with an initial guess. Must be of the
#' same dimension as 'b'.
#' @param iter Integer scalar with the maximum number of iterations. Defaults
#' to the theoretical maximum, i.e. the number of columns in 'a'.
#' @param tol Numeric scalar with the desired tolerance. Defaults to the
#' machine precision.
#' @param precond Integer scalar indicating the type of preconditioner to be
#' used. Defaults to diagonal preconditioning. See the Details for further
#' information.
#' @param verbose Logical scalar. Whether to print iterations and tolerance.
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#'
#' @export
#' @examples
#' set.seed(42)
#' x <- rnorm(3)
#'
#' # Solve via BiCGSTAB for square matrices
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' b <- A %*% x
#' norm(solve_cg(A, b, type = "B") - x)
#'
#' # Solve via LSCG for rectangular matrices
#' A <- matrix(rnorm(12), nrow = 4, ncol = 3)
#' b <- A %*% x
#' norm(solve_cg(A, b, type = "LS") - x)
#'
#' # Solve via classic CG for symmetric matrices
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' b <- A %*% x
#' norm(solve_cg(A, b, type = "CG") - x)
#'
#' # The input matrix A should always be in sparse format
#' A <- sparsify(crossprod(matrix(rnorm(9), nrow = 3, ncol = 3)))
#' # The right-hand side should be a dense matrix
#' b <- as.matrix(A %*% x)
#'
#' # We can check the speed of convergence and quality directly
#' solve_cg(A, b, verbose = TRUE)
#' # And provide guesses as starting value
#' solve_cg(A, b, x0 = x, verbose = TRUE)
#'
solve_cg <- function(a, b, x0,
  type = c("BiCGSTAB", "LSCG", "CG"),
  iter, tol, precond = 1L, verbose = FALSE) {

  # Checks -----
  type <- match.arg(type)

  if(is.matrix(a)) {a <- sparsify(a)} # Has to be sparse
  if(!is_sparse(a)) {stop("Please provide 'a' as sparse 'dgCMatrix'.")}

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(type != "LSCG" && !is_square(a)) {
    stop("Please provide a square matrix 'a' or use `type = 'LSCG'`.")
  }
  if(dim(a)[1L] != length(b)) {
    stop("The dimensions of the matrix 'a' must match the length of 'b'.")
  }

  precond <- int_check(precond, 0L, 2L)
  verbose <- isTRUE(verbose)

  if(!missing(x0)) { # Gotta check it
    if(!is.numeric(x0) || (!is.vector(x0) && !is.matrix(x0))) {
      stop("Please provide the guess 'x0' as a numeric vector or matrix.")
    }
  } else { # Construct it
    x0 <- if(is.vector(b)) {
      rep(0, min(dim(a)))
    } else {
      matrix(0, dim(a)[2L], dim(b)[2L])
    }
  }
  if(!missing(tol)) {
    tol <- num_check(tol, min = 0, max = Inf,
      msg = "Please provide a valid tolerance 'tol'.")
  } else {tol <- 0} # Uses Eigen's default (i.e. machine precision)
  if(!missing(iter)) {
    iter <- int_check(tol, min = 0, max = max(dim(a)),
      msg = "Please provide a valid number of iterations 'iter'.")
  } else {iter <- 0L} # Uses Eigen's default (i.e. the number of columns in a)

  # Execute -----
  if(type == "BiCGSTAB") {
    return(solve_BiCGSTAB(a, b, x0 = x0,
      tol = tol, iter = iter, precond = precond, verbose = verbose))
  } else if(type == "LSCG") {
    return(solve_LSCG(a, b, x0 = x0,
      tol = tol, iter = iter, precond = precond, verbose = verbose))
  } else { # CG
    return(solve_CG(a, b, x0 = x0,
      tol = tol, iter = iter, precond = precond, verbose = verbose))
  }
}
