#' Solve Systems of Equations
#'
#' Solve systems of equations \eqn{Ax = b}{Ax = b} using an automatically
#' chosen direct method (see \code{\link{solve_chol}}). Methods are chosen for
#' speed at reasonable accuracy. Please choose a suitable method manually if
#' numerical stability is the main consideration.
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#' @param a Numeric matrix with the coefficients of the linear system, dense or
#'   sparse (see \code{\link{sparsify}}), or a \code{sanic_lu} factorization
#'   from \code{\link{factorize_lu}}.
#' @param b Numeric vector or matrix at the right-hand side of the linear
#'   system. If missing, 'b' is set to an identity matrix and 'a' is
#'   inverted.
#' @param ... Dispatched to methods in the solvers.
#'
#' @export
#' @examples
#' set.seed(42)
#' x <- rnorm(3)
#'
#' # Solve using a general matrix
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' b <- A %*% x
#' norm(solve2(A, b) - x)
#'
#' # Solve using a symmetric matrix
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' b <- A %*% x
#' norm(solve2(A, b) - x)
#'
#' # Solve using a square matrix
#' A <- matrix(rnorm(12), nrow = 4, ncol = 3)
#' b <- A %*% x
#' norm(solve2(A, b) - x)
#'
#' # Solve reusing an LU factorization
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' f <- factorize_lu(A)
#' norm(solve2(f, A %*% x) - x)
#'
solve2 <- function(a, b, ...) {
  UseMethod("solve2")
}

#' @rdname solve2
#' @export
solve2.default <- function(a, b, ...) {
  # Check whether a and b fit together
  if (!is_square(a)) {
    # Use QR for rectangular problems
    return(solve_qr(a, b, ...))
  }
  if (maybe_symmetric(a)) {
    # Cholesky for symmetric ones
    if (is_symmetric(a, tol = 0)) {
      # if(all(diag(a) > 0) || all(diag(a) < 0)) { # What was the point of this?
      # return(solve_chol(a, b, ...))
      # } else {
      ans <- tryCatch(solve_chol(a, b, ...), error = function(e) NULL)
      if (
        !is.null(ans) &&
          all(is.finite(ans)) &&
          isTRUE(all.equal(
            as.matrix(a) %*% ans,
            if (missing(b)) diag(dim(a)[1L]) else as.matrix(b),
            tolerance = sqrt(.Machine$double.eps)
          ))
      ) {
        return(ans)
      }
      return(solve_lu(a, b, ...))
    }
  }
  return(solve_lu(a, b, ...)) # LU for the rest
}

#' @rdname solve2
#' @export
solve2.sanic_lu <- function(a, b, ...) {
  solve(a, b, ...)
}
