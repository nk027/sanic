#' Solve Systems of Equations
#'
#' @inheritParams solve_chol
#'
#' @return Solves for \eqn{x} and returns a numeric matrix with the results.
#'
#' @noRd
solve2 <- function(a, b,
  symmetric = FALSE, sparse = FALSE) {

  if(isTRUE(sparse) && !is_sparse(a)) {a <- sparsify(a)}
  symmetric <- isTRUE(symmetric)

  if(!is_square(a)) {
    return(solve_qr(a, b, pivot = 1L))
  }
  if(symmetric || maybe_symmetric(a)) {
    if(symmetric || is_symmetric(a, tol = 0, checks = FALSE)) {
      if(all(diag(a) > 0) || all(diag(a) < 0)) {
        return(solve_chol(a, b, pivot = 0L))
      } else {
        return(solve_chol(a, b, pivot = 1L))
      }
    }
  }
  return(solve_lu(a, b, pivot = 1L))
}
