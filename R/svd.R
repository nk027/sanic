
#' Singular Value Decomposition
#'
#' Solvers for generalized eigenproblems around the matrix \eqn{A}. Compute
#' singular values \eqn{\Sigma}{s}, left singular vectors \eqn{U}{U} and
#' right singular vectors \eqn{V}{V} of \eqn{A}{A}, such that
#' \eqn{A = U \Sigma V^*}{A = USV^*}. Two different types are available: (1)
#' bidiagonal divide and conquer strategy (BDC) SVD, and (2) two-sided Jacobi
#' SVD for small matrices (<16) and high accuracy.
#'
#' @param a Numeric matrix.
#' @param type Character scalar. Whether to use BDC or Jacobi SVD.
#' @param vectors Logical scalar indicating whether singular vectors should be
#' computed and returned.
#' @param thin Logical scalar indicating whether singular vectors should be
#' returned in thinned or full format.
#'
#' @return Solves the generalised eigenproblem and returns a list with
#' singular values in the \code{"d"} component and, if requested, singular
#' vectors in the components \code{"u"} and \code{"v"}.
#'
#' @export
#' @examples
#' set.seed(42)
#' # Compute singular values and vectors using BDC
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' sv <- svd2(A)
#'
#' # Compute singular values using Jacobi
#' A <- matrix(rnorm(9), nrow = 3, ncol = 3)
#' sv <- svd2(A, type = "J", vectors = FALSE)
#'
#' # Compute singular values and full vectors using BDC
#' A <- matrix(rnorm(12), nrow = 4, ncol = 3)
#' sv <- svd2(A, type = "B", thin = FALSE)
#' A <- matrix(rnorm(12), nrow = 3, ncol = 4)
#' sv <- svd2(A, type = "B", thin = FALSE)
#'
svd2 <- function(a, type = c("BDC", "Jacobi"), vectors = TRUE, thin = TRUE) {

  # Checks -----
  type <- match.arg(type)

  if(!is.matrix(a)) {stop("Please provide 'a' as matrix or sparse 'dgCMatrix'.")}
  vectors <- isTRUE(vectors)
  thin <- isTRUE(thin)

  thin_full <- if(vectors) {if(thin) {0} else {1}} else {2}

  # Execute -----
  if(type == "Jacobi") {
    return(svd_J(a, type = thin_full))
  } else {
    return(svd_BDC(a, type = thin_full))
  }
}
