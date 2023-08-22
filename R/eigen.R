
#' Spectral Decomposition
#'
#' Solvers for eigenproblems around the matrix \eqn{A}. Compute eigenvalues
#' \eqn{\lambda}{l} and eigenvectors \eqn{v}{v} of \eqn{A}{A}, such that
#' \eqn{Av = \lambda v}{Av = lv}.
#'
#' @param a Square numeric matrix.
#' @param symmetric Logical scalar indicating whether 'a' is symmetric. By
#' default symmetry is checked up to machine precision, which may take a long
#' time for symmetric matrices.
#' @param vectors Logical scalar indicating whether eigenvectors should be
#' computed and returned.
#'
#' @return Solves the eigenproblem and returns a list with eigenvalues in
#' the \code{"values"} slot and, if requested, eigenvectors in the slot
#' \code{"vectors"}.
#'
#' @export
#' @examples
#' set.seed(42)
#' # Compute eigenvalues and eigenvectors for a square matrix
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' ev <- eigen2(A, symmetric = FALSE)
#'
#' # Compute eigenvalues and eigenvectors for a symmetric matrix
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' ev <- eigen2(A, symmetric = TRUE)
#' # Check reconstruction
#' norm(A %*% ev$vectors - ev$vectors %*% diag(ev$values))
#'
eigen2 <- function(a, symmetric, vectors = TRUE) {

  # Checks -----
  if(!is.matrix(a)) {stop("Please provide 'a' as matrix.")}
  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}
  if(missing(symmetric)) {symmetric <- is_symmetric(a, tol = 0)}
  symmetric <- isTRUE(symmetric)
  vectors <- isTRUE(vectors)

  # Execute -----
  if(symmetric) { # Self-adjoint solver
    return(eigen_SA(a, vectors = vectors))
  } else { # Square solver
    return(eigen_SQ(a, vectors = vectors))
  }
}
