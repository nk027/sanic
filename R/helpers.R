
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

#' Check Sparsity
#'
#' @param x Matrix.
#'
#' @return Returns a logical scalar indicating sparsity.
#'
#' @noRd
is_sparse <- function(x) {
  isTRUE(inherits(x, "dgCMatrix"))
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


#' Check Symmetry
#'
#' @param x A numeric matrix.
#' @param tol A numeric scalar with the desired tolerance. The default value of
#' 0 is coerced to the machine precision.
#' @param checks A logical scalar indicating whether 'x' should be checked
#' for squarity and type before checking symmetry.
#'
#' @return Returns a logical scalar indicating symmetry.
#'
#' @noRd
is_symmetric <- function(x, tol = 0) {
  is_sparse <- is_sparse(x)
  if(!is.matrix(x) && !is_sparse) {stop("Please provide a matrix")}
  if(!is_square(x)) {return(FALSE)}
  if(is_sparse) {
    isTRUE(is_symmetric_S(x, tol = tol))
  } else {
    isTRUE(is_symmetric_D(x, tol = tol))
  }
}


#' Check Non-Symmetry Lazily
#'
#' @param x A numeric matrix.
#' @param tol A numeric scalar with the desired tolerance. The default value of
#' 0 is coerced to the machine precision.
#'
#' @return Returns a logical scalar possibly indicating symmetry.
#'
#' @noRd
maybe_symmetric <- function(x, tol = .Machine$double.eps, checks = TRUE) {
  if(isTRUE(checks)) {
    if(!is_square(x)) {return(FALSE)}
    if(!is.matrix(x)) {stop("Please provide a matrix")}
  }
  pos <- c(1, nrow(x))
  isTRUE(all.equal.numeric(x[pos, ], x[, pos], tolerance = tol))
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
