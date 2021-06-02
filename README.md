
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

Solver | Notes | Sparse method | Reference
--- | --- | --- | ---
LU decomposition | Partial pivoting, blocking | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU)
Householder QR decomposition | Column pivoting, good reliability | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR)
Cholesky decomposition | Semidefinite symmetric problems, pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLDLT)
Conjugate Gradient (CG) | Symmetric problems, Jacobi preconditioner | Always | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient)
Least Squares (LS) CG | Rectangular LS problems, LS Jacobi preconditioner | Always | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LeastSquaresConjugateGradient)
Biconjugate gradient stabilised | Square problems, Jacobi preconditioner | Always | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB)
