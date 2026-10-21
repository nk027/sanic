# sanic: Solving Ax = b Nimbly in C++

[![CRAN](https://www.r-pkg.org/badges/version/sanic)](https://CRAN.R-project.org/package=sanic)
[![month](https://cranlogs.r-pkg.org/badges/sanic)](https://www.r-pkg.org/pkg/sanic)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/sanic)](https://www.r-pkg.org/pkg/sanic)

**sanic** provides direct and iterative solvers for linear systems in **R**,
powered by the [*Eigen*](https://eigen.tuxfamily.org) **C++** library. It
supports dense and sparse systems, reusable LU factorizations, eigenvalue and
singular value decompositions, and Arnoldi and Lanczos projections.

## Installation

Install the released version from [CRAN](https://CRAN.R-project.org/package=sanic) or the development version from GitHub:

```r
install.packages("sanic")
remotes::install_github("nk027/sanic")
```

## Testing

The test suite uses [tinytest](https://github.com/markvanderloo/tinytest) and
runs when the `RUN_TESTS` environment variable is true. Set it for a single
check:

```sh
RUN_TESTS=true R CMD check --no-manual sanic_0.0.2.tar.gz
```

or add `RUN_TESTS=true` to a project-local `.Renviron` and run the suite
interactively after installing with tests:

```sh
R CMD INSTALL --install-tests .
cd tests && Rscript -e 'tinytest::test_package("sanic")'
```

Coverage can be measured with `covr`:

```r
RUN_TESTS=true Rscript -e 'Sys.setenv(RUN_TESTS = "true"); covr::package_coverage()'
```

## Usage

Use `solve2()` to select a direct solver automatically: QR for rectangular
input, a Cholesky attempt for symmetric input (with LU fallback), and LU
otherwise.

```r
library(sanic)

A <- matrix(c(4, 1, 1, 3), nrow = 2)
b1 <- c(1, 2)
b2 <- c(3, 4)

x <- solve2(A, b1)
x
#            [,1]
# [1,] 0.09090909
# [2,] 0.63636364
```

For repeated LU solves with the same matrix, compute the factorization once:

```r
f <- factorize_lu(A)
x1 <- solve(f, b1)
x2 <- solve(f, b2)
cbind(x1, x2)
```

Factors are session-local and must be recreated from the original matrix after
saving/loading or when moving to another R process.

Call `solve_lu()`, `solve_qr()`, or `solve_chol()` directly to choose a method
explicitly. For iterative solves, use `solve_cg()` and select the method explicitly. The default is BiCGSTAB; for classic conjugate gradient use:

```r
x_cg <- solve_cg(sparsify(A), b1, type = "CG", tol = 1e-10)
max(abs(A %*% x_cg - b1))
```

Classic CG requires a symmetric positive definite matrix. Iterative methods
offer preconditioners and convergence controls; heed convergence warnings and
check the residual.

Use `eigen2()` for eigenvalue decompositions and `svd2()` for singular value
decompositions of dense matrices. For large sparse matrices, `arnoldi()` and
`lanczos()` construct smaller Krylov projections. Set `iter` below the matrix dimension to build a smaller projection.

For a visual guide to the methods and their assumptions, see the [usage
vignette](vignettes/basic-usage.Rmd), or open it after installation:

```r
vignette("basic-usage", package = "sanic")
```

For timings against base R and the `Matrix` package, see the [benchmarking
vignette](vignettes/benchmarks.Rmd):

```r
vignette("benchmarks", package = "sanic")
```

## Linear system solvers

The pivoting options below refer to dense methods; sparse methods use their own
pivoting and ordering schemes. For Cholesky, positive definite input is a safe
starting point for both LLT and LDLT.

| Solver | Function | Notes | Sparse support | Reference |
|---|---|---|---|---|
| LU decomposition | `solve_lu()` | Partial or full pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html) |
| Householder QR | `solve_qr()` | Column, full, or no pivoting; rectangular least squares | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivHouseholderQR.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR.html) |
| Cholesky decomposition | `solve_chol()` | LDLT by default; LLT with `pivot = 0` | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLDLT.html), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLLT.html) |
| Iterative Krylov methods | `solve_cg()` | BiCGSTAB (default); LSCG for least squares; CG for symmetric positive definite systems; preconditioners available | Converted¹ | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LeastSquaresConjugateGradient.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html) |

¹ Uses sparse matrices internally; dense inputs are converted.

## Spectral methods

| Method | Function | Notes | Sparse support | Reference |
|---|---|---|---|---|
| Eigenvalue decomposition | `eigen2()` | General square and symmetric matrices | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html) |
| Singular value decomposition | `svd2()` | Divide-and-conquer SVD for large problems; Jacobi SVD for small problems | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html) |
| Arnoldi iteration | `arnoldi()` | Hessenberg projection; uses Lanczos for symmetric matrices | Converted¹ | [1](https://en.wikipedia.org/wiki/Arnoldi_iteration) |
| Lanczos algorithm | `lanczos()` | Tridiagonal projection for symmetric matrices | Converted¹ | [1](https://en.wikipedia.org/wiki/Lanczos_algorithm) |

¹ Uses sparse matrices internally; dense inputs are converted.

## Citation

To cite **sanic**, please use:

> Kuschnig, N. (2026). *sanic: Solving Ax = b Nimbly in C++* (R package version 0.0.2). https://CRAN.R-project.org/package=sanic

BibTeX:

```bibtex
@Manual{sanic,
  title  = {{sanic}: Solving {A}x = b Nimbly in {C++}},
  author = {Nikolas Kuschnig},
  year   = {2026},
  note   = {R package version 0.0.2},
  url    = {https://CRAN.R-project.org/package=sanic}
}
```

## License

See [LICENSE](LICENSE) for package licensing and third-party attribution,
including the MPL-2.0 header adapted from Eigen and Spectra.
