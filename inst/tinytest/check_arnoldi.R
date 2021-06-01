
Rcpp::sourceCpp("src/lanczos.cpp")

N <- 500
z <- rnorm(N) # Random vector to start with
x1 <- sanic::sparsify(y1 <- matrix(rnorm(N^2), N)) # Non-symmetric
x2 <- sanic::sparsify(y2 <- crossprod(matrix(rnorm(N^2), N))) # Symmetric

mb(
  "eigen" = eigen(x1, FALSE, TRUE),
  "pracma" = pracma::arnoldi(y1, z, m = N),
  "arnoldi" = arnoldi(x1, z, iter = N, tol = 1e-12),
  "hessenberg" = hessenberg(x1),
  "eigen_symm" = eigen(x2, TRUE, TRUE),
  "mgcv" = mgcv::slanczos(x2, k = N, tol = 1e-12),
  "lanczos" = lanczos(x2, z, iter = N, tol = 1e-12),
  "tridiagonal" = tridiagonal(x2),
  times = 10
)

out <- list(
  "eigen" = eigen(x1, FALSE, TRUE),
  "pracma" = pracma::arnoldi(y1, z, m = N),
  "arnoldi" = arnoldi(x1, z, iter = N, tol = 1e-12),
  "hessenberg" = hessenberg(x1),
  "eigen_symm" = eigen(x2, TRUE, TRUE),
  "mgcv" = mgcv::slanczos(x2, k = N, tol = 1e-12),
  "arnoldi_symm" = arnoldi(x2, z, iter = N, tol = 1e-12),
  "lanczos" = lanczos(x2, z, iter = N, tol = 1e-12),
  "tridiagonal" = tridiagonal(x2)
)

ev <- sort(Re(out$eigen$values))
summary(ev - sort(out$arnoldi$values))
summary(ev - sort(out$hessenberg$values))

ev <- sort(out$eigen_symm$values)
summary(ev - sort(out$mgcv$values))
summary(ev - sort(out$arnoldi_symm$values))
summary(ev - sort(out$lanczos$values)) # Duplicated Ritz values
summary(ev - sort(out$tridiagonal$values))

plot(out$lanczos$values)


mgcv::slanczos(x2)$values
sort(eigen(x2, symmetric = TRUE)$values)
lanczos(x2, z, tol = 1e-12)$T
tridiagonal(x2)$T

arnoldi(x2, z)
plot(lanczos(x2, z)$T)
points(lanczos2(x2, z)$T, col = "red")
tridiagonal(x2)
eigen(x2, symmetric = TRUE)$values

sort(hessenberg(x)$T)

sort(arnoldi(x, z)$T)

arnoldi(x, z, iter = 10)$T

a1 <- arnoldi(x)
a2 <- arnoldi(x, iter = N / 2)
h1 <- hessenberg(x)
t1 <- tridiagonal(x)

plot(sort(a1$T, TRUE))
points(sort(a2$T, TRUE), col = "#800000")
points(sort(h1$T, TRUE), col = "#008080")
points(sort(t1$E, TRUE), col = "#000080")

mb(arnoldi(x), arnoldi(x, iter = N / 2), hessenberg(x), tridiagonal(x),
  times = 2)
