
#' Solve Systems of Equations Using Direct Methods
#'
#' Direct solvers using Cholesky, LU, or QR decompositions for systems of
#' equations \eqn{Ax = b}{Ax = b}. Dense or sparse methods are used depending
#' on the format of the input matrix (see \code{\link{sparsify}}).
#'
#' @details
#' Pivoting schemes for dense matrices can be set to \code{0} for no pivoting,
#' \code{1} (default) for partial pivoting, or \code{2} for full pivoting. Not
#' all schemes are available for every decomposition:
#'
#' * \code{solve_chol()} The default is \code{pivot = 1} for the robust LDLT
#' decomposition of \eqn{A}{A}, such that \eqn{A = P'LDL^*P}{A = P'LDL^*P}. For
#' the LDLT \eqn{A}{A} needs to be positive or negative semidefinite. Set
#' \code{pivot = 0} for the plain LLT decomposition of \eqn{A}{A}, such that
#' \eqn{A = LL^* = U^*U}{A = LL^* = U^*U}. For the LLT \eqn{A}{A} needs to be
#' positive definite and preferably numerically stable.
#' * \code{solve_lu()} The default is \code{pivot = 1} for the partial pivoting
#' LU decomposition of \eqn{A}, such that \eqn{A = PLU}{A = PLU}. For this
#' scheme \eqn{A} needs to be invertible and preferably numerically stable. Set
#' \code{pivot = 2} for the complete pivoting LU decomposition of \eqn{A}{A},
#' such that \eqn{A = P^{-1}LUQ^{-1}}{A = PLUQ}. This scheme is applicable to
#' square matrices, rank-revealing, and stable.
#' \code{solve_qr()} The default is \code{pivot = 1} for the column pivoting
#' Householder QR decomposition of \eqn{A}{A}, such that \eqn{AP = QR}{AP = QR}.
#' This scheme is generally applicable, rank-revealing, and stable. Set
#' \code{pivot = 2} for the full pivoting Householder QR decomposition of
#' \eqn{A}{A}, such that \eqn{PAP' = QR}{PAP' = QR}. This scheme is generally
#' applicable, rank-revealing, and optimally stable. Set \code{pivot = 0} for
#' an unpivoted Householder QR decomposition of \eqn{A}{A}, such that
#' \eqn{A = QR}{A = QR}. This scheme is generally applicable, but not as stable
#' as pivoted variants.
#'
#' Ordering schemes for sparse matrices can be set to \code{0} for approximate
#' minimum degree (AMD) ordering, \code{1} for column approximate minimum degree
#' (COLAMD) ordering, or \code{2} for natural ordering. Not all orderings are
#' available for every decomposition:
#'
#' * \code{solve_chol()} The default is \code{ordering = 0} for AMD ordering.
#' Set \code{ordering = 2} for natural ordering.
#' * \code{solve_lu()} The default is \code{ordering = 1} for COLAMD ordering.
#' Set \code{ordering = 0} for AMD or \code{ordering = 2} for natural ordering.
#' * \code{solve_qr()} The default is \code{ordering = 1} for COLAMD ordering.
#' Set \code{ordering = 0} for AMD or \code{ordering = 2} for natural ordering.
#'
#' @param a Square numeric matrix with the coefficients of the linear system.
#' Both dense and sparse matrices are supported (see \code{\link{sparsify}}).
#' @param b Numeric vector or matrix at the right-hand side of the linear
#' system. If missing, 'b' is set to an identity matrix and 'a' is
#' inverted.
#' @param pivot Integer scalar indicating the pivoting scheme to be used.
#' Defaults to partial pivoting. See the Details for further information.
#' @param ordering Integer scalar indicating the ordering scheme to be used.
#' See the Details for further information.
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#'
#' @export
#' @examples
#' set.seed(42)
#' x <- rnorm(3)
#'
#' # Solve via QR for general matrices
#' A <- matrix(rnorm(12), nrow = 4, ncol = 3)
#' b <- A %*% x
#' norm(solve_qr(A, b) - x)
#'
#' # Solve via LU for square matrices
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' b <- A %*% x
#' norm(solve_lu(A, b) - x)
#'
#' # Solve via Cholesky for symmetric matrices
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' b <- A %*% x
#' norm(solve_chol(A, b) - x)
#'
#' # Sparse methods are available for the 'dgCMatrix' class from Matrix
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' b <- A %*% x
#' norm(solve_qr(sparsify(A), b))
#' norm(solve_lu(sparsify(A), b))
#' norm(solve_chol(sparsify(A), b))
#'
solve_chol <- function(a, b, pivot = 1L, ordering = 0L) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- is_sparse(a)
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix or sparse 'dgCMatrix'.")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  pivot <- int_check(pivot, 0L, 2L)
  ordering <- int_check(ordering, 0L, 2L)

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_LL(a, b, pivot = pivot))
  } else { # Sparse
    return(solve_SLL(a, b, pivot = pivot, ordering = ordering))
  }
}


#' @rdname solve_chol
#' @export
solve_lu <- function(a, b, pivot = 1L, ordering = 1L) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- is_sparse(a)
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix or sparse 'dgCMatrix'.")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  pivot <- int_check(pivot, 0L, 2L)
  ordering <- int_check(ordering, 0L, 2L)

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_LU(a, b, pivot = pivot))
  } else { # Sparse
    return(solve_SLU(a, b, ordering = ordering))
  }
}


#' @rdname solve_chol
#' @export
solve_qr <- function(a, b, pivot = 1L, ordering = 1L) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- is_sparse(a)
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix or sparse 'dgCMatrix'.")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  b_rows <- if(is.null(dim(b))) {length(b)} else {dim(b)[1L]}
  is_fitting <- dim(a)[1L] - b_rows == 0L
  if(!is_fitting) {stop("Both 'a' and 'b' must have the same number of rows.")}

  pivot <- int_check(pivot, 0L, 2L)
  ordering <- int_check(ordering, 0L, 2L)

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_HQR(a, b, pivot = pivot))
  } else { # Sparse
    return(solve_SQR(a, b, ordering = ordering))
  }
}
