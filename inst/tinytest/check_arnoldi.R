
Rcpp::sourceCpp("src/lanczos.cpp")

N <- 1000
z <- rnorm(N)
x <- sanic::sparsify(y <- matrix(rnorm(N^2), N))
x <- sanic::sparsify(y <- crossprod(matrix(rnorm(N^2), N)))

mb(pracma::arnoldi(y, z, m = N), arnoldi(x, z, iter = N))

a <- arnoldi(x, z, iter = N)
b <- pracma::arnoldi(y, z, m = N)
all.equal(Re(eigen(a$H)$values), Re(eigen(b$H)$values))

mb(
  eigen(y, only.values = TRUE),
  eigen(arnoldi(x, z, iter = N)$H, only.values = TRUE),
  eigen(pracma::arnoldi(y, z, m = N)$H, only.values = TRUE),
  times = 10
)

mb(
  hessenberg(x),
  arnoldi(x, z, tol = 1e-12)
)

mgcv::slanczos(x)

sort(eigen(x, symmetric = TRUE)$values)

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
