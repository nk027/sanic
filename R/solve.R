#' Solve Systems of Equations
#'
#' Solve systems of equations \eqn{Ax = b}{Ax = b} using an automatically
#' chosen direct method (see \code{\link{solve_chol}}). Methods are chosen for
#' speed at reasonable accuracy. Please choose a suitable method manually if
#' numerical stability is the main consideration.
#'
#' @inheritParams solve_chol
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
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
solve2 <- function(a, b, ...) {

  # Check whether a and b fit together
  if(!is_square(a)) { # Use QR for rectangular problems
    return(solve_qr(a, b, ...))
  }
  if(maybe_symmetric(a)) { # Cholesky for symmetric ones
    if(is_symmetric(a, tol = 0)) {
      # if(all(diag(a) > 0) || all(diag(a) < 0)) { # What was the point of this?
        # return(solve_chol(a, b, ...))
      # } else {
      return(solve_chol(a, b, ...))
      # }
    }
  }
  return(solve_lu(a, b, ...)) # LU for the rest
}
