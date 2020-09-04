
#' Solve a System of Equations Using Direct Methods
#'
#' Functions to access specific direct solvers for systems of equations.
#'
#' @param a Square numeric matrix with the coefficients of the linear system.
#' Both dense and sparse matrices are supported (see \link{sparsify}).
#' @param b Numeric vector or matrix at the right-hand side of the linear
#' system. If missing, 'b' is set to an identity matrix and 'a' is
#' inverted.
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#'
#' @export
#' @examples
#' # Solve via LU and QR for general matrices
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' x <- rnorm(3)
#' b  <- A %*% x
#'
#' x_lu <- solve_lu(A, b)
#' x_qr <- solve_qr(A, b)
#'
#' # Solve via Cholesky for symmetric matrices
#' AA <- crossprod(A)
#' b <- AA %*% x
#'
#' x_chol <- solve_chol(AA, b)
#'
#' # Sparse methods are available for the 'dgCMatrix' class from Matrix
#' x_slu <- solve_lu(sparsify(A), b)
#'
solve_chol <- function(a, b) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_LDLT(a, b))
  } else { # Sparse
    return(solve_SLDLT(a, b))
  }
}


#' @rdname solve_chol
#' @export
solve_lu <- function(a, b) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(inherits(b, "Matrix")) {b <- as.matrix(b)}
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_PPLU(a, b))
  } else { # Sparse
    return(solve_SLU(a, b))
  }
}


#' @rdname solve_chol
#' @export
solve_qr <- function(a, b) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
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

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_CPHQR(a, b))
  } else { # Sparse
      return(solve_SQR(a, b))
  }
}
