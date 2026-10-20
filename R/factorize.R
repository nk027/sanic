#' Reusable LU Factorization
#'
#' Compute an LU factorization once and reuse it for different right-hand sides.
#'
#' @param a For \code{factorize_lu()}, a nonempty square numeric matrix, dense
#'   or sparse (\code{dgCMatrix}). For \code{solve()}, a \code{sanic_lu}
#'   factorization.
#' @param pivot Pivoting scheme for dense matrices, as in
#'   \code{\link{solve_lu}}.
#' @param ordering Ordering scheme for sparse matrices, as in
#'   \code{\link{solve_lu}}.
#' @details
#' The factorization owns its data: changing or removing \code{a} does not
#' affect subsequent solves. Dense LU uses partial pivoting by default;
#' \code{pivot = 2} uses full pivoting and checks for singularity. As with
#' \code{solve_lu()}, partial pivoting assumes an invertible matrix and does not
#' perform a rank check. Sparse LU ignores \code{pivot}.
#'
#' Factorizations are session-local external pointers. They cannot be reused
#' after serialization (for example, with \code{saveRDS()}) or in another R
#' process; recreate them from the original matrix instead.
#' @return
#' \code{factorize_lu()} returns an object of class \code{sanic_lu}. Its
#' \code{solve()} method returns a numeric matrix, including for vector
#' right-hand sides. Omitting \code{b} computes the inverse.
#' @export
#' @examples
#' A <- matrix(c(4, 1, 2, 3), 2)
#' f <- factorize_lu(A)
#' solve(f, c(1, 2))
#' solve(f, cbind(c(1, 2), c(3, 4)))
#' solve(f)
factorize_lu <- function(a, pivot = 1L, ordering = 1L) {
  if ((!is.matrix(a) || !is.numeric(a)) && !is_sparse(a)) {
    stop("Please provide 'a' as a numeric matrix or sparse 'dgCMatrix'.")
  }
  if (!is_square(a) || nrow(a) == 0L) {
    stop("Please provide a nonempty square matrix 'a'.")
  }
  pivot <- int_check(pivot, 0L, 2L)
  ordering <- int_check(ordering, 0L, 2L)
  if (is_sparse(a)) {
    factorize_SLU(a, ordering)
  } else {
    storage.mode(a) <- "double"
    factorize_LU(a, pivot)
  }
}

#' @rdname factorize_lu
#' @param b Numeric vector or matrix of right-hand sides. If missing, compute
#'   the inverse. Sparse right-hand sides are converted to dense matrices.
#' @param ... Unused; additional arguments are rejected.
#' @export
solve.sanic_lu <- function(a, b, ...) {
  if (length(list(...))) {
    stop("Additional arguments are not supported for LU factorizations.")
  }
  if (missing(b)) {
    return(solve_LU_factor(a, NULL))
  }
  if (inherits(b, "Matrix")) {
    b <- as.matrix(b)
  }
  if (!is.numeric(b) || (!is.vector(b) && !is.matrix(b))) {
    stop("Please provide 'b' as a numeric vector or matrix.")
  }
  b <- as.matrix(b)
  storage.mode(b) <- "double"
  solve_LU_factor(a, b)
}
