
Rcpp::sourceCpp("src/lanczos.cpp")

N <- 10
x <- sanic::sparsify(y <- matrix(rnorm(N^2), N))
x <- sanic::sparsify(y <- crossprod(matrix(rnorm(N^2), N)))
z <- rnorm(N)

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
