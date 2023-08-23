
sanic: Solving Ax = b Nimbly in C++
=======

[![CRAN](https://www.r-pkg.org/badges/version/sanic)](https://CRAN.R-project.org/package=sanic)
[![month](https://cranlogs.r-pkg.org/badges/sanic)](https://www.r-pkg.org/pkg/sanic)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/sanic)](https://www.r-pkg.org/pkg/sanic)

Routines for solving large systems of linear equations in **R**. Direct and iterative solvers from the [*Eigen*](https://eigen.tuxfamily.org) **C++** library are made available. Solvers include Cholesky, LU, QR, and Krylov subspace methods (Conjugate Gradient, BiCGSTAB). Both dense and sparse problems are supported.

Installation
-------

**sanic** is available on [CRAN](https://CRAN.R-project.org/package=sanic). The development version can be installed from GitHub.
```r
install.packages("sanic")
devtools::install_github("nk027/sanic")
```

Usage
-------

To solve a linear system of equations, use the `solve2()` function for automatic dispatch to a specific solver or access the LU, QR, Cholesky or Conjugate Gradient solvers directly.

To solve an eigenproblem, use the `eigen2()` or `svd2()` functions, or the `arnoldi()` function.

Solvers
-------

Solver | Function | Notes | Sparse | Reference
--- | --- | --- | --- | ---
LU decomposition | `solve_lu()` | Partial pivoting, full pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU.html)
Householder QR decomposition | `solve_qr()` | Column pivoting, full pivoting, no pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1FullPivLU.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR.html), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR.html)
Cholesky decomposition | `solve_chol()` | LDLT for semidefinite problems, LLT for positive definite problems | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html) [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLDLT.html), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLLT.html)
Conjugate Gradient (CG) | `solve_cg()` | Biconjugate gradient stabilised (BiCGTAB) for square problems, least squares (LSCG) for rectangular problems, classic CG for symmetric positive definite problems, preconditioners | Always | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LeastSquaresConjugateGradient.html), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html)

Eigenproblems
-------

Solver | Function | Notes | Sparse | Reference
--- | --- | --- | --- | ---
Spectral decomposition | `eigen2()` | Square and symmetric problems | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html)
Singular value decomposition | `svd2()` | Bidiagonal Divide and Conquer SVD for large and Jacobi SVD for small problems | No | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BDCSVD.html), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1JacobiSVD.html)
Arnoldi iteration | `arnoldi()` | Square problems using an iteratively constructed Hessenberg matrix | Always | [1](https://en.wikipedia.org/wiki/Arnoldi_iteration)
Lanczos algorithm | `lanczos()` | Symmetric problems using an iteratively constructed tridiagonal matrix | Always | [1](https://en.wikipedia.org/wiki/Lanczos_algorithm)
