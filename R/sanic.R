#' @title Solving Ax = b Nimbly in C++
#'
#' @description Routines for solving large systems of linear equations in R.
#' Direct and iterative solvers from the Eigen C++ library are made available.
#' Solvers include Cholesky, LU, QR, and Krylov subspace methods (Conjugate
#' Gradient, BiCGSTAB). Both dense and sparse problems are supported.
#'
#' @docType package
#' @name sanic
#' @importFrom Rcpp sourceCpp
#' @useDynLib sanic
NULL
