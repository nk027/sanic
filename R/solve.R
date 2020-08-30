
#' @title Solve a System of Equations using Cholesky factorization
#'
#'
#'
#' @param a
#' @param b rhs
#' @param type Whether to decompose to LDLT or LLT
#'
#' @export
#' @examples
#' solve_chol(crossprod(matrix(rnorm(9), 3)), rnorm(3))
solve_chol <- function(a, b, type = c("LDLT", "LLT")) {

  # Checks -----
  type <- match.arg(type)

  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide a as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else {
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide b as a numeric vector or matrix.")
    }

  is_square <- diff(dim(a)) == 0L
  if(!is_square) {stop("Please provide a square matrix a.")}

  # Execute -----
  if(is_matrix) { # Dense
    if(type == "LDLT") {
      return(solve_LDLT(a, b))
    } else {
      return(solve_LLT(a, b))
    }
  } else { # Sparse
    if(type == "LDLT") {
      return(solve_SLDLT(a, b))
    } else {
      return(solve_SLLT(a, b))
    }
  }
}


#' @title Solve a System of Equations using LU factorization
#'
#' @param a lhs
#' @param b rhs
#'
#' @export
#' @examples
#' solve_lu(matrix(rnorm(9), 3), rnorm(3))
solve_lu <- function(a, b) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide a as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else {
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide b as a numeric vector or matrix.")
    }

  is_square <- diff(dim(a)) == 0L
  if(!is_square) {stop("Please provide a square matrix a.")}

  # Execute -----
  if(is_matrix) { # Dense
    return(solve_PPLU(a, b))
  } else { # Sparse
    return(solve_SLU(a, b))
  }
}


#' @title Solve using QR
#'
#' @param a lhs
#' @param b rhs
#' @param pivot Whether to pivot
#'
#' @export
#' @examples
#' solve_qr(matrix(rnorm(9), 3), rnorm(3))
solve_qr <- function(a, b, pivot = TRUE) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide a as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else {
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide b as a numeric vector or matrix.")
    }
  }

  b_rows <- if(is.null(dim(b))) {length(b)} else {dim(b)[1L]}
  is_fitting <- diff(dim(a)[1L], b_rows) == 0L
  if(!is_fitting) {stop("Both a and b must have the same number of rows.")}

  pivot <- isTRUE(pivot)

  # Execute -----
  if(is_matrix) { # Dense
    if(pivot) {
      return(solve_CPHQR(a, b))
    } else {
      return(solve_HQR(a, b))
    }
  } else { # Sparse
    if(pivot) {
      return(solve_SQR(a, b))
    } else {
      stop("No sparse solver without pivoting is available.")
    }
  }
}
