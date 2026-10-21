library(tinytest)
library(Matrix)

set.seed(3)
A <- as(Diagonal(4, x = rep(2, 4)), "dgCMatrix")
b <- c(1, 2, 3, 4)

# A scalar multiple of the identity has a one-dimensional Krylov space.
la <- lanczos(A, b = b, tol = 1e-12)
expect_equal(length(la$diagonal), 1L)
expect_equal(length(la$subdiagonal), 0L)
expect_equal(ncol(la$Q), 1L)
expect_equal(length(la$values), 1L)
expect_equal(unname(la$values), 2)

la_plain <- lanczos(A, b = b, eigen = FALSE, orthogonalise = FALSE,
  tol = 1e-12)
expect_equal(names(la_plain), c("diagonal", "subdiagonal", "Q"))

ar <- arnoldi(A, b = b, symmetric = FALSE, tol = 1e-12)
expect_equal(dim(ar$H), c(1L, 1L))
expect_equal(ncol(ar$Q), 1L)
expect_equal(length(ar$values), 1L)
expect_equal(unname(Mod(ar$values)), 2)

ar_plain <- arnoldi(A, b = b, symmetric = FALSE, eigen = FALSE,
  tol = 1e-12)
expect_equal(names(ar_plain), c("H", "Q"))

# Dense input is coerced to sparse, and symmetry is detected by default.
expect_equal(arnoldi(as.matrix(A), b = b, tol = 1e-12)$values,
  lanczos(as.matrix(A), b = b, tol = 1e-12)$values)
expect_equal(ncol(arnoldi(A, b = b, tol = 1e-12)$Q), 1L)
expect_equal(ncol(arnoldi(A, b = rnorm(4), tol = 1e-12)$Q), 1L)

# Full requested Krylov dimension remains available when there is no early
# breakdown (a diagonal matrix with distinct eigenvalues and a dense start).
B <- as(Diagonal(4, x = 1:4), "dgCMatrix")
full <- arnoldi(B, b = b, symmetric = FALSE, iter = 3)
expect_equal(dim(full$H), c(3L, 3L))
expect_equal(ncol(full$Q), 3L)

full_lanczos <- lanczos(as(B, "dgCMatrix"), b = b, iter = 2)
expect_equal(length(full_lanczos$diagonal), 2L)
expect_equal(ncol(full_lanczos$Q), 2L)

# Input validation ----------------------------------------------------------
expect_error(arnoldi(list(1), b = b, symmetric = FALSE), "sparse")
expect_error(arnoldi(matrix(1, 2, 3), b = c(1, 2), symmetric = FALSE),
  "square")
expect_error(arnoldi(A, b = b, symmetric = FALSE, tol = -1), "tolerance")
expect_error(arnoldi(A, b = b, symmetric = FALSE, iter = 0), "iterations")
expect_error(arnoldi(A, b = b, symmetric = FALSE, iter = 5), "iterations")
expect_error(arnoldi(A, b = rep(0, 4), symmetric = FALSE), "Starting vector")
expect_error(lanczos(A, b = rep(0, 4)), "Starting vector")
expect_error(lanczos(A, b = c(1, 2)), "length")

# A random start vector is generated when 'b' is missing.
random_start <- arnoldi(A, symmetric = TRUE, tol = 1e-12)
expect_equal(length(random_start$values), 1L)
