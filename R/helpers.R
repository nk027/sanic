
#' @title The Way
#'
#' @param x Matrix or something.
#'
#' @return
the_way <- function(x) {

  is_square(x)
  is.complex(x)
  is.finite(x)
  isSymmetric.matrix(x)

  # Consider Sparsity (`coop::sparsity()`)
  # Direct | Iterative
  # Cholesky > LU > QR
  # Diagonal dominance?

  # More: Spectral and SV decompositions, more pivoting, stationary solvers
}

#' @title Transform a Matrix to Be Sparse.
#'
#' Concise function to transform dense to sparse matrices of class
#' \code{dgCMatrix} (see \link[Matrix]{sparseMatrix}).
#'
#' @param x Numeric matrix to transform to a sparse 'dgCMatrix'.
#'
#' @return Returns 'x' as \code{dgCMatrix}.
#'
#' @export
#' @examples
#' sparsify(matrix(rnorm(9L), 3L))
sparsify <- function(x) {
  return(as(x, "dgCMatrix"))
}

#' @title Check Squarity
#'
#' @param x Matrix or something else with dimensions.
#'
#' @return Returns a logical scalar indicating squarity.
is_square <- function(x) {
  isTRUE(diff(dim(x)) == 0L)
}