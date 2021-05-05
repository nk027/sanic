
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

#' Check numeric scalar
#'
#' Check whether an object is bounded and coercible to a numeric value.
#'
#' @param x Numeric scalar.
#' @param min Numeric scalar. Minimum value of \emph{x}.
#' @param max Numeric scalar. Maximum value of \emph{x}.
#' @param fun Function to apply to \emph{x} before returning.
#' @param msg String fed to \code{\link[base]{stop}} if an error occurs.
#'
#' @return Returns \code{fun(x)}.
#'
#' @noRd
num_check <- function(
  x, min = 0, max = Inf,
  msg = "Please check the numeric parameters.",
  fun = as.numeric) {

  if(!is.numeric(x) || length(x) != 1 || x < min || x > max) {stop(msg)}

  return(fun(x))
}


#' @noRd
int_check <- function(
  x, min = 0L, max = Inf,
  msg = "Please check the integer parameters.") {

  num_check(x, min, max, msg, fun = as.integer)
}
