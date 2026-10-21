library(tinytest)
library(Matrix)

expect_same <- function(current, target, tolerance = 1e-10) {
  expect_equal(unname(current), unname(target), tolerance = tolerance)
}

set.seed(42)
A <- crossprod(matrix(rnorm(16), 4, 4)) + diag(4)
b <- rnorm(4)
B <- cbind(b, 2 * b, -b)
A_sparse <- sparsify(A)
expected <- solve(A, B)
rect <- matrix(rnorm(12), 4, 3)
rect_b <- rnorm(4)
rect_expected <- matrix(qr.solve(rect, rect_b), ncol = 1)

# Dense Cholesky ------------------------------------------------------------
for (pivot in 0:1) {
  expect_same(solve_chol(A, b, pivot = pivot), expected[, 1, drop = FALSE],
    tolerance = 1e-10)
  expect_same(solve_chol(A, B, pivot = pivot), expected,
    tolerance = 1e-10)
}
expect_same(solve_chol(A, Matrix(B, sparse = TRUE)), expected,
  tolerance = 1e-10)
expect_same(solve_chol(A), solve(A), tolerance = 1e-10)

# Sparse Cholesky -----------------------------------------------------------
for (ordering in c(0, 2)) {
  expect_same(solve_chol(A_sparse, B, ordering = ordering), expected,
    tolerance = 1e-10)
}
expect_same(solve_chol(A_sparse), solve(A), tolerance = 1e-10)

# Dense LU ------------------------------------------------------------------
for (pivot in 0:2) {
  expect_same(suppressWarnings(solve_lu(A, B, pivot = pivot)), expected,
    tolerance = 1e-10)
  expect_same(suppressWarnings(solve_lu(A, b, pivot = pivot)),
    expected[, 1, drop = FALSE], tolerance = 1e-10)
}
expect_warning(solve_lu(A, b, pivot = 0), "No unpivoted")
expect_same(solve_lu(A, Matrix(B, sparse = TRUE)), expected,
  tolerance = 1e-10)
expect_same(solve_lu(A), solve(A), tolerance = 1e-10)

# Sparse LU -----------------------------------------------------------------
for (ordering in 0:2) {
  expect_same(solve_lu(A_sparse, B, ordering = ordering), expected,
    tolerance = 1e-10)
}
expect_same(solve_lu(A_sparse), solve(A), tolerance = 1e-10)

# Dense QR ------------------------------------------------------------------
for (pivot in 0:2) {
  expect_same(solve_qr(rect, rect_b, pivot = pivot), rect_expected,
    tolerance = 1e-10)
}
x_rect <- solve_qr(rect, rect_b)
expect_same(nrow(x_rect), ncol(rect))
expect_same(ncol(x_rect), 1L)
expect_same(solve_qr(A, B), expected, tolerance = 1e-10)
expect_same(solve_qr(A), solve(A), tolerance = 1e-10)
expect_same(solve_qr(A, Matrix(B, sparse = TRUE)), expected,
  tolerance = 1e-10)

# Sparse QR -----------------------------------------------------------------
for (ordering in 0:2) {
  expect_same(suppressWarnings(solve_qr(sparsify(rect), rect_b, ordering = ordering)),
    rect_expected, tolerance = 1e-10)
}
expect_warning(solve_qr(sparsify(rect), rect_b, ordering = 0), "AMD")
expect_same(solve_qr(A_sparse, B, ordering = 0), expected, tolerance = 1e-10)
expect_same(solve_qr(A_sparse, B), expected, tolerance = 1e-10)
expect_same(solve_qr(A_sparse), solve(A), tolerance = 1e-10)

# Input validation ----------------------------------------------------------
expect_error(solve_chol(list(1)), "matrix")
expect_error(solve_chol(A, letters[1:4]), "numeric")
expect_error(solve_chol(rect), "square")
expect_error(solve_chol(A, pivot = 3), "parameters")
expect_error(solve_chol(A, ordering = -1), "parameters")
expect_error(solve_lu(list(1)), "matrix")
expect_error(solve_lu(A, "x"), "numeric")
expect_error(solve_lu(rect), "square")
expect_error(solve_lu(A, pivot = 3), "parameters")
expect_error(solve_lu(A, ordering = 3), "parameters")
expect_error(solve_qr(list(1)), "matrix")
expect_error(solve_qr(A, "x"), "numeric")
expect_error(solve_qr(rect, b[-1]), "same number of rows")
expect_error(solve_qr(A, pivot = 3), "parameters")
expect_error(solve_qr(A, ordering = 3), "parameters")

# Automatic direct solver selection -----------------------------------------
expect_same(solve2(rect, rect_b), solve_qr(rect, rect_b),
  tolerance = 1e-10)
expect_same(solve2(A, b), solve_chol(A, b), tolerance = 1e-10)

# Sparse inputs are routed through the same selection (regression).
expect_same(solve2(A_sparse, B), expected, tolerance = 1e-10)
expect_same(solve2(A_sparse, b), expected[, 1, drop = FALSE], tolerance = 1e-10)
expect_same(solve2(A_sparse), solve(A), tolerance = 1e-10)
expect_same(solve2(sparsify(rect), rect_b), rect_expected, tolerance = 1e-10)

# A symmetric matrix that is not positive definite falls back to LU.
indefinite <- matrix(c(0, 1, 1, 0), 2)
indefinite_b <- c(1, 2)
expect_same(solve2(indefinite, indefinite_b),
  solve_lu(indefinite, indefinite_b), tolerance = 1e-10)

# A general square matrix uses LU.
general <- matrix(c(0, 2, 1, 3, 4, 1, 2, 0, 5), 3)
general_b <- c(1, 2, 3)
expect_same(solve2(general, general_b), solve_lu(general, general_b),
  tolerance = 1e-10)

# Less common Cholesky paths ------------------------------------------------
expect_warning(solve_chol(A, b, pivot = 2), "No valid pivoting")
expect_same(suppressWarnings(solve_chol(A, b, pivot = 2)), expected[, 1, drop = FALSE],
  tolerance = 1e-10)
expect_warning(solve_chol(A_sparse, b, ordering = 1), "No COLAMD")
expect_same(suppressWarnings(solve_chol(A_sparse, b, ordering = 1)),
  expected[, 1, drop = FALSE], tolerance = 1e-10)
expect_warning(solve_chol(A_sparse, b, pivot = 2), "No valid pivoting")
expect_same(solve_chol(A_sparse, b, pivot = 0, ordering = 0),
  expected[, 1, drop = FALSE], tolerance = 1e-10)
expect_same(solve_chol(A_sparse, b, pivot = 0, ordering = 2),
  expected[, 1, drop = FALSE], tolerance = 1e-10)
