
sanic: Solving Ax = b Nimbly in C++
=======

[![CRAN](http://www.r-pkg.org/badges/version/sanic)](http://cran.r-project.org/package=sanic)
[![codecov](https://codecov.io/gh/nk027/sanic/branch/master/graph/badge.svg)](https://codecov.io/gh/nk027/sanic)
[![month](http://cranlogs.r-pkg.org/badges/sanic)](http://www.r-pkg.org/pkg/sanic)
[![total](http://cranlogs.r-pkg.org/badges/grand-total/sanic)](http://www.r-pkg.org/pkg/sanic)

Routines for solving large systems of linear equations in **R**. Direct and iterative solvers from the [*Eigen*](https://eigen.tuxfamily.org) **C++** library are made available. Solvers include Cholesky, LU, QR, and Krylov subspace methods (Conjugate Gradient, BiCGSTAB). Both dense and sparse problems are supported.

Installation
-------

**sanic** is available on [CRAN](https://CRAN.R-project.org/package=sanic). The development version can be installed from GitHub.
``` r
install.packages("sanic")
devtools::install_github("nk027/sanic")
```

Solvers
-------

Solver | Function | Notes | Sparse | Reference
--- | --- | --- | --- | ---
LU decomposition | `solve_lu()` | Partial pivoting, full pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU) [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU)
Householder QR decomposition | `solve_qr()` | Column pivoting, full pivoting, no pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR)
Cholesky decomposition | `solve_chol()` | LDLT for semidefinite problems, LLT for positive definite problems | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LLT) [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLDLT), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLLT)
Conjugate Gradient (CG) | `solve_cg()` | Biconjugate gradient stabilised (BiCGTAB) for square problems, least squares (LSCG) for rectangular problems, classic CG for symmetric positive definite problems, preconditioners | Always | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LeastSquaresConjugateGradient), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient)

Eigenproblems
--------

Solver | Function | Notes | Sparse | Reference
--- | --- | --- | --- | ---
Spectral decomposition | `eigen2()` | Square and symmetric problems | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver)
Singular value decomposition | `svd2()` | Bidiagonal Divide and Conquer SVD for large and Jacobi SVD for small problems | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD)
Arnoldi iteration | `arnoldi()` | Square problems using an iteratively constructed Hessenberg matrix (cf. `hessenberg()`) | Yes | [1](https://en.wikipedia.org/wiki/Arnoldi_iteration)
Lanczos algorithm | `lanczos()` | Symmetric problems using an iteratively constructed tridiagonal matrix(cf. `tridiagonal()`) | Yes | [1](https://en.wikipedia.org/wiki/Lanczos_algorithm)
