library(tinytest)
library(Matrix)

A <- matrix(c(0, 2, 1, 3, 4, 1, 2, 0, 5), 3)
B <- matrix(as.numeric(seq_len(6)), 3)
b <- B[, 1]

check_factor <- function(a, ...) {
  f <- factorize_lu(a, ...)
  expect_true(inherits(f, "sanic_lu"))
  for (rhs in list(b, B, B * 2, Matrix(B, sparse = TRUE))) {
    x <- solve(f, rhs)
    expect_true(is.matrix(x))
    expect_equal(A %*% x, as.matrix(rhs), tolerance = 1e-12)
    expect_equal(x, solve_lu(a, rhs, ...), tolerance = 1e-12)
  }
  expect_equal(solve(f), solve(A), tolerance = 1e-12)
  expect_equal(dim(solve(f, matrix(numeric(), 3, 0))), c(3L, 0L))
  expect_equal(solve_lu(a, ...), solve(A), tolerance = 1e-12)
  expect_equal(dim(solve_lu(a, matrix(numeric(), 3, 0), ...)), c(3L, 0L))
  f
}

for (pivot in 1:2) {
  check_factor(A, pivot = pivot)
}
for (ordering in 0:2) {
  check_factor(sparsify(A), ordering = ordering)
}
f <- suppressWarnings(factorize_lu(A, pivot = 0))
expect_equal(solve(f, B), solve(A, B), tolerance = 1e-12)

# Integer matrices and right-hand sides are accepted.
f <- factorize_lu(matrix(c(4L, 1L, 2L, 3L), 2))
expect_true(is.matrix(solve(f, 1:2)))

# Factors survive mutation/removal of the input, and garbage collection.
for (sparse in c(FALSE, TRUE)) {
  for (scheme in if (sparse) 0:2 else 1:2) {
    original <- if (sparse) sparsify(A) else A
    f <- if (sparse) {
      factorize_lu(original, ordering = scheme)
    } else {
      factorize_lu(original, pivot = scheme)
    }
    if (sparse) {
      original@x[] <- 0
    } else {
      original[,] <- 0
    }
    rm(original)
    invisible(gc())
    expect_equal(solve(f, B), solve(A, B), tolerance = 1e-12)
  }
}

f <- factorize_lu(A)
expect_error(solve(f, 1:2), "same number of rows")
expect_error(solve(f, matrix(1, 2, 3)), "same number of rows")
expect_error(solve(f, c(1, NA, 2)), "finite")
expect_error(solve(f, letters[1:3]), "numeric")
expect_error(solve(f, NULL), "numeric")
expect_error(solve(f, B, tol = 1e-6), "Additional arguments")
expect_error(factorize_lu(matrix(1, 2, 3)), "square")
expect_error(factorize_lu(matrix(numeric(), 0, 0)), "nonempty")
expect_error(factorize_lu(matrix("a", 2, 2)), "numeric")
expect_error(factorize_lu(A, pivot = 3), "parameters")
expect_error(factorize_lu(A, ordering = 3), "parameters")
for (sparse in c(FALSE, TRUE)) {
  bad <- A
  bad[1, 1] <- Inf
  if (sparse) {
    bad <- sparsify(bad)
  }
  expect_error(factorize_lu(bad), "finite")
  expect_error(solve_lu(bad, B), "finite")
  a <- if (sparse) sparsify(A) else A
  expect_error(solve_lu(a, c(1, NA, 2)), "finite")
  expect_error(solve_lu(a, c(1, 2)), "same number of rows")
  expect_error(solve_lu(a, matrix(1, 2, 3)), "same number of rows")
}
expect_error(factorize_lu(matrix(1, 3, 3), pivot = 2), "singular")
expect_error(factorize_lu(sparsify(matrix(1, 3, 3))), "failed")
expect_error(solve_lu(matrix(1, 3, 3), b, pivot = 2), "singular")
expect_error(solve_lu(sparsify(matrix(1, 3, 3)), b), "failed")

# Serializing does not preserve native state, and must fail safely on reuse.
expired <- unserialize(serialize(f, NULL))
expect_error(solve(expired, B), "no longer valid")
expect_error(solve(expired), "no longer valid")
forged <- structure(list(), class = "sanic_lu")
expect_error(solve(forged, B), "Invalid LU")
forged <- structure(new("externalptr"), class = "sanic_lu")
expect_error(solve(forged, B), "Invalid LU")

# Existing one-shot entry points retain their behavior.
expect_equal(solve_lu(A, B), solve(A, B), tolerance = 1e-12)
expect_equal(solve_lu(A), solve(A), tolerance = 1e-12)
expect_equal(solve2(A, B), solve(A, B), tolerance = 1e-12)

# Reusable factorizations are accepted by the solve2() generic.
f <- factorize_lu(A)
expect_equal(solve2(f, B), solve(A, B), tolerance = 1e-12)
expect_equal(solve2(f, b), solve(f, b), tolerance = 1e-12)
expect_equal(solve2(f), solve(A), tolerance = 1e-12)
expect_error(solve2(f, B, tol = 1e-6), "Additional arguments")
expect_error(solve2(f, 1:2), "same number of rows")
for (ordering in 0:2) {
  fs <- factorize_lu(sparsify(A), ordering = ordering)
  expect_equal(solve2(fs, B), solve(A, B), tolerance = 1e-12)
}
