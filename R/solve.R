
#' @title Solve using Cholesky
#'
#' @param A lhs
#' @param b rhs
#' @param type Whether to decompose to LDLT or LLT
#'
#' @export
#' @examples
#' solve_chol(crossprod(matrix(rnorm(9), 3)), rnorm(3))
solve_chol <- function(A, b, type = c("LDLT", "LLT")) {

  # Checks -----
  is_matrix <- inherits(A, "matrix")
  is_sparse <- inherits(A, "dgCMatrix")
  if(!is_matrix && !is_sparse) stop("Please provide a matrix.")

  is_square <- diff(dim(A)) == 0
  if(!is_square) stop("Please provide a square matrix (or use SVD).")

  type <- match.arg(type)

  if(missing(b)) b <- diag(nrow(A))

  # Execute -----
  if(is_matrix)
    if(type == "LDLT") return(solve_LDLT(A, b)) else return(solve_LLT(A, b))
  if(type == "LDLT") return(solve_SLDLT(A, b)) else return(solve_SLLT(A, b))
}


#' @title Solve using LU
#'
#' @param A lhs
#' @param b rhs
#'
#' @export
#' @examples
#' solve_lu(matrix(rnorm(9), 3), rnorm(3))
solve_lu <- function(A, b) {

  # Checks -----
  is_matrix <- inherits(A, "matrix")
  is_sparse <- inherits(A, "dgCMatrix")
  if(!is_matrix && !is_sparse) stop("Please provide a matrix.")

  is_square <- diff(dim(A)) == 0
  if(!is_square) stop("Please provide a square matrix (or use SVD).")

  if(missing(b)) b <- diag(nrow(A))

  # Execute -----
  if(is_matrix) return(solve_PPLU(A, b))
  return(solve_SLU(A, b))
}


#' @title Solve using QR
#'
#' @param A lhs
#' @param b rhs
#' @param pivot Whether to pivot
#'
#' @export
#' @examples
#' solve_qr(matrix(rnorm(9), 3), rnorm(3))
solve_qr <- function(A, b, pivot = FALSE) {

  # Checks -----
  is_matrix <- inherits(A, "matrix")
  is_sparse <- inherits(A, "dgCMatrix")
  if(!is_matrix && !is_sparse) stop("Please provide a matrix.")

  if(missing(b)) b <- diag(nrow(A))

  is_square <- diff(dim(A)) == 0
  if(!is_square) stop("Please provide a square matrix (or use SVD).")

  # Execute -----
  if(is_matrix)
    if(pivot) return(solve_CPHQR(A, b)) else return(solve_HQR(A, b))
  return(solve_SQR(A, b))
}


#' @title Solve using conjugate gradient
#'
#' @param A lhs
#' @param b rhs
#' @param type Whether to use the BiCGSTAB, least squares CG or CG solver
#' @param x0 Initial guess
#' @param iter Iterations
#' @param tol Tolerance
#' @param verbose Whether to print iterations and tolerance
#'
#' @importFrom Matrix Matrix
#' @importFrom methods as
#'
#' @export
#' @examples
#' solve_cg(matrix(rnorm(9), 3), rnorm(3))
solve_cg <- function(A, b,
  type = c("BiCGSTAB", "TNBiCGSTAB", "LSCG", "CG"),
  x0, tol, iter, verbose = FALSE) {

  # Checks -----
  is_sparse <- inherits(A, "dgCMatrix")
  is_matrix <- inherits(A, "matrix")
  if(!is_sparse)
    if(is_matrix) A <- as(A, "sparseMatrix") else stop("Please provide a matrix.")

  is_square <- diff(dim(A)) == 0
  if(!is_square) stop("Please provide a square matrix (or use SVD).")

  type <- match.arg(type)
  verbose <- isTRUE(verbose)

  if(!missing(x0)) stopifnot(is.numeric(x0)) else {
    x0 <- if(is.vector(b)) rep(0, length(b)) else matrix(0, nrow(b), ncol(b))
  }
  if(!missing(tol)) stopifnot(is.numeric(tol) && iter > 0) else tol <- 0L
  if(!missing(iter)) stopifnot(is.integer(iter) && iter > 0) else iter <- 0L

  if(missing(b)) b <- diag(nrow(A))

  # Execute -----
  if(type == "BiCGSTAB")
    return(solve_BCGST(A, b, x0 = x0, tol = tol, iter = iter, verbose = verbose))
  if(type == "DTNSBiCGSTAB")
    return(solve_DTNSBCGST(A, b, x0 = x0, tol = tol, iter = iter, verbose = verbose))
  if(type == "LSCG")
    return(solve_CGLS(A, b, x0 = x0, tol = tol, iter = iter, verbose = verbose))
  # if(type == "CG")
  return(solve_CGLS(A, b, x0 = x0, tol = tol, iter = iter, verbose = verbose))
}
