
# v0.1.0, CRAN update

- Add `tinytest` for exported solvers
  - Around 80% total coverage at the time
- Add `factorize_lu()` and `solve()` method for reusable factorizations
  - Share the LU core with one-shot solves
  - Check nonfinite inputs and RHS dimensions consistently
- Merge @eddelbuettel's Eigen 5.0.* PR
- Add attributes to iterative solver results
- Add an illustrative vignette on solvers
- Several bug fixes and minor tweaks

# v0.0.2, CRAN update

- Provide solvers for eigenvalue problems
  - Arnoldi and Lanczos methods for sparse matrices
    - Solver for the Upper Hessenberg form taken from the Spectra library
  - Interface to Eigen solvers for dense matrices
- Fix small bugs and tweak documentation
  - Fix symmetry check, default value for ordering, and return for CG solver
  - Remove dimension constraints for the QR solvers
  - Thanks to helpful comments by Bettina Grün

# v0.0.1, CRAN submission

- Provide interface to Eigen solvers
  - Direct: LU, QR, Cholesky
  - Iterative: Conjugate Gradient, BiCGSTAB, Least Squares CG
