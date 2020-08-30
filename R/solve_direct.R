
#' @title Solve a System of Equations Using Direct Methods
#'
#'
#'
#' @param a Square numeric matrix with the coefficients of the linear system.
#' Both dense and sparse matrices are supported (see \link{sparsify}).
#' @param b Numeric vector or matrix at the right-hand side of the linear
#' system. If missing, 'b' is set to an identity matrix and 'a' is
#' inverted.
#' @param type Character scalar. Whether to decompose 'a' as LDLT or LLT.
#' @param pivot Logical scalar. Whether to perform column pivoting.
#'
#' @returns Solves for \eqn{x} and returns a numeric matrix with the results.
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
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

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


#' @rdname solve_chol
solve_lu <- function(a, b) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
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
solve_qr <- function(a, b, pivot = TRUE) {

  # Checks -----
  is_matrix <- is.matrix(a)
  is_sparse <- inherits(a, "dgCMatrix")
  if(!is_matrix && !is_sparse) {
    stop("Please provide 'a' as matrix (or sparse 'dgCMatrix').")
  }

  if(missing(b)) {b <- diag(dim(a)[1L])} else { # Invert without b
    if(!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
      stop("Please provide 'b' as a numeric vector or matrix.")
    }
  }

  b_rows <- if(is.null(dim(b))) {length(b)} else {dim(b)[1L]}
  is_fitting <- diff(dim(a)[1L], b_rows) == 0L
  if(!is_fitting) {stop("Both 'a' and 'b' must have the same number of rows.")}

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
