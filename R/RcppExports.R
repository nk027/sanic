# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

lanczos_E <- function(a, b, tol = 0, iter = 0L, eigen = TRUE, orthogonalise = TRUE) {
    .Call('_sanic_lanczos_E', PACKAGE = 'sanic', a, b, tol, iter, eigen, orthogonalise)
}

arnoldi_E <- function(a, b, tol = 0, iter = 0L, eigen = TRUE) {
    .Call('_sanic_arnoldi_E', PACKAGE = 'sanic', a, b, tol, iter, eigen)
}

eigen_SA <- function(a, vectors = TRUE) {
    .Call('_sanic_eigen_SA', PACKAGE = 'sanic', a, vectors)
}

eigen_SQ <- function(a, vectors = TRUE) {
    .Call('_sanic_eigen_SQ', PACKAGE = 'sanic', a, vectors)
}

solve_LL <- function(a, b, pivot = 1L) {
    .Call('_sanic_solve_LL', PACKAGE = 'sanic', a, b, pivot)
}

solve_SLL <- function(a, b, pivot = 1L, ordering = 0L) {
    .Call('_sanic_solve_SLL', PACKAGE = 'sanic', a, b, pivot, ordering)
}

solve_LU <- function(a, b, pivot = 1L) {
    .Call('_sanic_solve_LU', PACKAGE = 'sanic', a, b, pivot)
}

solve_SLU <- function(a, b, ordering = 1L) {
    .Call('_sanic_solve_SLU', PACKAGE = 'sanic', a, b, ordering)
}

solve_HQR <- function(a, b, pivot = 1L) {
    .Call('_sanic_solve_HQR', PACKAGE = 'sanic', a, b, pivot)
}

solve_SQR <- function(a, b, ordering = 1L) {
    .Call('_sanic_solve_SQR', PACKAGE = 'sanic', a, b, ordering)
}

solve_BiCGSTAB <- function(a, b, x0, tol = 0, iter = 0L, precond = 1L, verbose = FALSE) {
    .Call('_sanic_solve_BiCGSTAB', PACKAGE = 'sanic', a, b, x0, tol, iter, precond, verbose)
}

solve_LSCG <- function(a, b, x0, tol = 0, iter = 0L, precond = 1L, verbose = FALSE) {
    .Call('_sanic_solve_LSCG', PACKAGE = 'sanic', a, b, x0, tol, iter, precond, verbose)
}

solve_CG <- function(a, b, x0, tol = 0, iter = 0L, precond = 1L, verbose = FALSE) {
    .Call('_sanic_solve_CG', PACKAGE = 'sanic', a, b, x0, tol, iter, precond, verbose)
}

svd_J <- function(a, type = 0L, precond = 0L) {
    .Call('_sanic_svd_J', PACKAGE = 'sanic', a, type, precond)
}

svd_BDC <- function(a, type = 0L) {
    .Call('_sanic_svd_BDC', PACKAGE = 'sanic', a, type)
}

is_symmetric_D <- function(x, tol = 0) {
    .Call('_sanic_is_symmetric_D', PACKAGE = 'sanic', x, tol)
}

is_symmetric_S <- function(x, tol = 0) {
    .Call('_sanic_is_symmetric_S', PACKAGE = 'sanic', x, tol)
}

