
#' Krylov Subspace Spectral Decomposition
#'
#' Arnoldi iteration and Lanczos method to iteratively approximate the
#' Hessenberg or tridiagonal form of a matrix \eqn{A}{A} and find its
#' eigenvalues.
#'
#' @inheritParams eigen2
#' @inheritParams solve_cg
#' @param b Arbitrary numeric non-zero vector used to construct the basis.
#' @param orthogonalise Logical scalar indicating whether to use plain Lanczos
#' or full reorthogonalisation. Defaults to reorthogonalisation.
#'
#' @return Returns a list with slots \code{"H"} for the Hessenberg form of 'a'
#' or slots \code{"diagonal"} and \code{"subdiagonal"} for its triangular form,
#' slot \code{"Q"} with the orthonormal basis, and, if requested, eigenvalues
#' in the slot \code{"values"}.
#'
#' @export
#' @examples
#' \dontshow{set.seed(42)}
#' # Compute Hessenberg of a square matrix
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' ks <- arnoldi(A, symmetric = FALSE)
#'
#' # Compute tridiagonal of a symmetric matrix
#' A <- crossprod(matrix(rnorm(9), nrow = 3, ncol = 3))
#' ks <- lanczos(A)
#' ks <- arnoldi(A, symmetric = TRUE) # Short-hand
#'
arnoldi <- function(a, b, symmetric,
  iter = nrow(a), tol = .Machine$double.eps,
  vectors = TRUE, orthogonalise = TRUE) {

  # Checks -----
  if(is.matrix(a)) {a <- sparsify(a)} # Has to be sparse
  if(!is_sparse(a)) {stop("Please provide 'a' as sparse 'dgCMatrix'.")}
  if(!is_square(a)) {stop("Please provide a square matrix 'a'.")}

  tol <- num_check(tol, min = 0, max = Inf,
    msg = "Please provide a valid tolerance 'tol'.")
  iter <- int_check(iter, min = 1L, max = nrow(a),
    msg = "Please provide a valid number of iterations 'iter'.")

  if(missing(symmetric)) {symmetric <- is_symmetric(a, tol = 0, checks = FALSE)}
  symmetric <- isTRUE(symmetric)
  if(missing(b)) {
    b <- rnorm(nrow(a))
  } else {stopifnot(length(b) == nrow(a))}
  vectors <- isTRUE(vectors)
  orthogonalise <- isTRUE(orthogonalise)

  # Execute -----
  if(symmetric) { # Lanczos algorithm
    return(lanczos_E(a, b = b,
      tol = tol, iter = iter, orthogonalise = orthogonalise))
  } else { # Arnoldi iteration
    return(arnoldi_E(a, b = b,
      tol = tol, iter = iter))
  }
}


#' @rdname arnoldi
#' @export
lanczos = function(a, b,
  iter = nrow(a), tol = .Machine$double.eps,
  vectors = TRUE, orthogonalise = TRUE) {
  arnoldi(a, symmetric = TRUE, ...)
}
