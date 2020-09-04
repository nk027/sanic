
#' Solve a System of Equations using Iterative Methods
#'
#' Function to use Conjugate Gradient (CG) methods to solver systems of
#' equations.
#'
#' @inheritParams solve_chol
#' @param type Character scalar. Whether to use the BiCGSTAB, least squares
#' CG or classic CG method.
#' @param x0 Numeric vector or matrix with an initial guess. Must be of the
#' same dimension as 'b'.
#' @param iter Integer scalar with the maximum number of iterations. Defaults
#' to the theoretical maximum, i.e. the number of columns in 'a'.
#' @param tol Numeric scalar with the desired tolerance. Defaults to the
#' machine precision.
#' @param verbose Logical scalar. Whether to print iterations and tolerance.
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#'
#' @export
#' @examples
#' # Solve via least squares or bi-conjugate gradient methods
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' # The matrix A should be of class 'dgCMatrix' (otherwise it is converted)
#' A <- sparsify(A)
#' x <- rnorm(3)
#' b  <- A %*% x
#'
#' x_bi <- solve_cg(A, b)
#' x_ls <- solve_cg(A, b, type = "LS")
#'
#' # Solve via conjugate gradient for symmetric matrices
#' AA <- A %*% A
#' b <- AA %*% x
#' x_cg <- solve_cg(AA, b, type = "CG")
#'
solve_cg <- function(a, b, x0,
  type = c("BiCGSTAB", "LSCG", "CG"),
  tol, iter, verbose = FALSE) {

  # Checks -----
  type <- match.arg(type)

  if(is.matrix(a)) {a <- sparsify(a)} # Has to be sparse
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_sparse) {stop("Please provide 'a' as sparse 'dgCMatrix'.")}

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  verbose <- isTRUE(verbose)

  if(!missing(x0)) { # Gotta check it
    if(!is.numeric(x0) || (!is.vector(x0) && !is.matrix(x0)) ||
      !all(class(b) == class(x0))) {
      stop("Please provide the guess 'x0' as 'a' numeric vector or matrix.")
    }
  } else { # Construct it
    x0 <- if(is.vector(b)) {
      rep(0, length(b))
    } else {
      matrix(0, dim(b)[1L], dim(b)[2L])
    }
  }
  if(!missing(tol)) {
    if(!is.numeric(tol) || tol <= 0) {stop("Please provide a valid 'tol'.")}
  } else {tol <- 0} # Uses Eigen's default (i.e. machine precision)
  if(!missing(iter)) {
    if(!is.numeric(iter) || iter < 1) {stop("Please provide a valid 'iter'.")}
    iter <- as.integer(iter)
  } else {iter <- 0L} # Uses Eigen's default (i.e. the number of columns in a)

  # Execute -----
  if(type == "BiCGSTAB") {
    return(solve_BiCGSTAB(a, b, x0 = x0,
      tol = tol, iter = iter, verbose = verbose))
  } else if(type == "LSCG") {
    return(solve_LSCG(a, b, x0 = x0,
      tol = tol, iter = iter, verbose = verbose))
  } else { # CG
    return(solve_CG(a, b, x0 = x0,
      tol = tol, iter = iter, verbose = verbose))
  }
}
