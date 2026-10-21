library(tinytest)
library(Matrix)

expect_same <- function(current, target, tolerance = 1e-8) {
  expect_equal(as.vector(unname(current)), as.vector(unname(target)),
    tolerance = tolerance)
}

set.seed(7)
A <- crossprod(matrix(rnorm(25), 5, 5)) + diag(5)
A_sparse <- sparsify(A)
b <- rnorm(5)
B <- cbind(b, -b)
rect <- matrix(rnorm(20), 5, 4)
rect_sparse <- sparsify(rect)
rect_b <- as.numeric(rect %*% c(1, -2, 3, 0.5))
rect_x <- matrix(qr.solve(rect, rect_b), ncol = 1)

# BiCGSTAB ------------------------------------------------------------------
x <- solve_cg(A_sparse, b)
expect_true(is.matrix(x))
expect_true(isTRUE(attr(x, "converged")))
expect_true(is.integer(attr(x, "iterations")))
expect_true(is.numeric(attr(x, "error")))
expect_same(as.numeric(A %*% x), b, tolerance = 1e-8)

expect_same(A %*% solve_cg(A_sparse, B), B, tolerance = 1e-8)
expect_same(A %*% solve_cg(A, b), b, tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, Matrix(B, sparse = TRUE)), B,
  tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse), diag(5), tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, b, type = "B", precond = 0),
  b, tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, b, x0 = b, tol = 1e-10, iter = 5),
  b, tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, B, x0 = B), B, tolerance = 1e-8)

# LSCG ----------------------------------------------------------------------
expect_same(rect %*% solve_cg(rect_sparse, rect_b, type = "LSCG"),
  rect_b, tolerance = 1e-8)
expect_same(rect %*% solve_cg(rect_sparse, rect_b, type = "LS",
  precond = 0), rect_b, tolerance = 1e-8)
expect_same(rect %*% solve_cg(rect_sparse, rect_b, type = "LSCG", x0 = rect_x),
  rect_b, tolerance = 1e-8)

# Classic CG ----------------------------------------------------------------
expect_same(A %*% solve_cg(A_sparse, b, type = "CG"), b,
  tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, b, type = "C", precond = 0),
  b, tolerance = 1e-8)
expect_same(A %*% solve_cg(A_sparse, b, type = "CG", precond = 2),
  b, tolerance = 1e-8)

# Input validation ----------------------------------------------------------
expect_error(solve_cg(list(1), b), "sparse")
expect_error(solve_cg(A_sparse, letters[1:5]), "numeric")
expect_error(solve_cg(rect_sparse, rect_b, type = "CG"), "square")
expect_error(solve_cg(A_sparse, b[-1]), "dimensions")
expect_error(solve_cg(A_sparse, B[-1, ]), "dimensions")
expect_error(solve_cg(A_sparse, b, x0 = "x"), "guess")
expect_error(solve_cg(A_sparse, b, precond = 3), "parameters")
expect_error(solve_cg(A_sparse, b, tol = -1), "tolerance")
expect_error(solve_cg(A_sparse, b, iter = -1), "iterations")
expect_error(solve_cg(A_sparse, b, iter = 6), "iterations")
expect_error(solve_cg(A_sparse, b, type = "unknown"))

# Verbose output and non-convergence are reported.
expect_stdout(solve_cg(A_sparse, b, verbose = TRUE), "Iterations:")
expect_warning(solve_cg(A_sparse, b, iter = 1, tol = 1e-16),
  "did not converge")
expect_warning(solve_cg(rect_sparse, rect_b, type = "LSCG", precond = 2),
  "No valid preconditioner")
