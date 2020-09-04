
# sanic: Solving Ax = b Nimbly in C++

Routines for solving large systems of linear equations in **R**. Direct and iterative solvers from the [*Eigen*](https://eigen.tuxfamily.org) **C++** library are made available. Solvers include Cholesky, LU, QR, and Krylov subspace methods (Conjugate Gradient, BiCGSTAB). Both dense and sparse problems are supported.

# Solvers

Solver | Notes | Sparse method | Reference
--- | --- | --- | ---
LU decomposition | Partial pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1PartialPivLU), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseLU)
Householder QR decomposition | Column pivoting | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ColPivHouseholderQR), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1HouseholderQR), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SparseQR)
Cholesky decomposition | LDLT and LLT | Yes | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT), [2](https://eigen.tuxfamily.org/dox/classEigen_1_1LLT), [3](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLDLT), [4](https://eigen.tuxfamily.org/dox/classEigen_1_1SimplicialLLT)
Conjugate Gradient (CG) | Initial guess, Preconditioner | Default | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient)
Least Squares CG | Initial guess, Preconditioner | Default | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1LeastSquaresConjugateGradient)
Biconjugate gradient stabilised | Initial guess, Preconditioner | Default | [1](https://eigen.tuxfamily.org/dox/classEigen_1_1BiCGSTAB)
