library(tinytest)

set.seed(11)
S <- crossprod(matrix(rnorm(16), 4, 4))
G <- matrix(rnorm(16), 4, 4)
wide <- matrix(rnorm(12), 3, 4)
tall <- matrix(rnorm(12), 4, 3)

# Eigenvalue decomposition --------------------------------------------------
ev <- eigen2(S, symmetric = TRUE)
expect_equal(sort(ev$values), sort(eigen(S, symmetric = TRUE)$values),
  tolerance = 1e-10)
expect_equal(S %*% ev$vectors, ev$vectors %*% diag(ev$values),
  tolerance = 1e-8)

ev_auto <- eigen2(S)
expect_equal(sort(ev_auto$values), sort(ev$values), tolerance = 1e-10)

ev_general <- eigen2(G, symmetric = FALSE)
expect_equal(sort(Mod(ev_general$values)),
  sort(Mod(eigen(G, only.values = TRUE)$values)), tolerance = 1e-8)

values_only <- eigen2(G, symmetric = FALSE, vectors = FALSE)
expect_equal(names(values_only), "values")

values_sym <- eigen2(S, symmetric = TRUE, vectors = FALSE)
expect_equal(names(values_sym), "values")
expect_equal(sort(values_sym$values), sort(ev$values), tolerance = 1e-10)
expect_equal(sort(Mod(values_only$values)),
  sort(Mod(eigen(G, only.values = TRUE)$values)), tolerance = 1e-8)

expect_error(eigen2(as.vector(S)), "matrix")
expect_error(eigen2(tall), "square")

# Singular value decomposition ----------------------------------------------
for (a in list(G, tall, wide)) {
  sv <- svd2(a)
  k <- length(sv$d)
  expect_equal(sv$d, svd(a)$d, tolerance = 1e-8)
  expect_equal(a, sv$u %*% diag(sv$d) %*% t(sv$v),
    tolerance = 1e-8)

  sv_full <- svd2(a, type = "B", thin = FALSE)
  expect_equal(sv_full$d, svd(a)$d, tolerance = 1e-8)
  expect_equal(dim(sv_full$u), c(nrow(a), nrow(a)))
  expect_equal(dim(sv_full$v), c(ncol(a), ncol(a)))
  expect_equal(a, sv_full$u[, seq_len(k)] %*% diag(sv_full$d) %*%
    t(sv_full$v[, seq_len(k)]), tolerance = 1e-8)

  sv_jacobi <- svd2(a, type = "Jacobi")
  expect_equal(sv_jacobi$d, svd(a)$d, tolerance = 1e-8)

  sv_jacobi_full <- svd2(a, type = "Jacobi", thin = FALSE)
  expect_equal(sv_jacobi_full$d, svd(a)$d, tolerance = 1e-8)

  sv_bdc_values <- svd2(a, vectors = FALSE)
  expect_equal(names(sv_bdc_values), "d")
  expect_equal(sv_bdc_values$d, svd(a)$d, tolerance = 1e-8)

  sv_values <- svd2(a, type = "J", vectors = FALSE)
  expect_equal(names(sv_values), "d")
  expect_equal(sv_values$d, svd(a)$d, tolerance = 1e-8)
}

expect_error(svd2(as.vector(G)), "numeric matrix")
expect_error(svd2(matrix(letters[1:4], 2)), "numeric matrix")
expect_error(svd2(G + 1i), "numeric matrix")
expect_error(svd2(G, type = "unknown"))
