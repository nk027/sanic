
#' Transform a Matrix to Be Sparse.
#'
#' Concise function to transform dense to sparse matrices of class
#' \code{dgCMatrix} (see \link[Matrix]{sparseMatrix}).
#'
#' @param x Numeric matrix to transform to a sparse 'dgCMatrix'.
#'
#' @return Returns 'x' as \code{dgCMatrix}.
#'
#' @importFrom Matrix Matrix
#' @importFrom methods as
#'
#' @export
#' @examples
#' sparsify(matrix(rnorm(9L), 3L))
sparsify <- function(x) {
  return(as(x, "dgCMatrix"))
}

#' Check Squarity
#'
#' @param x Matrix or something else with dimensions.
#'
#' @return Returns a logical scalar indicating squarity.
#'
#' @noRd
is_square <- function(x) {
  isTRUE(diff(dim(x)) == 0L)
}
