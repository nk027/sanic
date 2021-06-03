// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// lanczos_E
Rcpp::List lanczos_E(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map <Eigen::VectorXd> b, double tol, unsigned int iter, bool eigen, bool orthogonalise);
RcppExport SEXP _sanic_lanczos_E(SEXP aSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP iterSEXP, SEXP eigenSEXP, SEXP orthogonaliseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::VectorXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type eigen(eigenSEXP);
    Rcpp::traits::input_parameter< bool >::type orthogonalise(orthogonaliseSEXP);
    rcpp_result_gen = Rcpp::wrap(lanczos_E(a, b, tol, iter, eigen, orthogonalise));
    return rcpp_result_gen;
END_RCPP
}
// arnoldi_E
Rcpp::List arnoldi_E(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map <Eigen::VectorXd> b, double tol, unsigned int iter, bool eigen);
RcppExport SEXP _sanic_arnoldi_E(SEXP aSEXP, SEXP bSEXP, SEXP tolSEXP, SEXP iterSEXP, SEXP eigenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::VectorXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type eigen(eigenSEXP);
    rcpp_result_gen = Rcpp::wrap(arnoldi_E(a, b, tol, iter, eigen));
    return rcpp_result_gen;
END_RCPP
}
// hessenberg_E
Rcpp::List hessenberg_E(const Eigen::Map<Eigen::MatrixXd> a);
RcppExport SEXP _sanic_hessenberg_E(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(hessenberg_E(a));
    return rcpp_result_gen;
END_RCPP
}
// tridiagonal_E
Rcpp::List tridiagonal_E(const Eigen::Map<Eigen::MatrixXd> a);
RcppExport SEXP _sanic_tridiagonal_E(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(tridiagonal_E(a));
    return rcpp_result_gen;
END_RCPP
}
// eigen_SA
Rcpp::List eigen_SA(const Eigen::Map<Eigen::MatrixXd> a, bool vectors);
RcppExport SEXP _sanic_eigen_SA(SEXP aSEXP, SEXP vectorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type vectors(vectorsSEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_SA(a, vectors));
    return rcpp_result_gen;
END_RCPP
}
// eigen_SQ
Rcpp::List eigen_SQ(const Eigen::Map<Eigen::MatrixXd> a, bool vectors);
RcppExport SEXP _sanic_eigen_SQ(SEXP aSEXP, SEXP vectorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< bool >::type vectors(vectorsSEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_SQ(a, vectors));
    return rcpp_result_gen;
END_RCPP
}
// solve_LL
Eigen::MatrixXd solve_LL(const Eigen::Map <Eigen::MatrixXd> a, const Eigen::Map <Eigen::MatrixXd> b, unsigned int pivot);
RcppExport SEXP _sanic_solve_LL(SEXP aSEXP, SEXP bSEXP, SEXP pivotSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type pivot(pivotSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_LL(a, b, pivot));
    return rcpp_result_gen;
END_RCPP
}
// solve_SLL
Eigen::MatrixXd solve_SLL(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, unsigned int pivot, unsigned int ordering);
RcppExport SEXP _sanic_solve_SLL(SEXP aSEXP, SEXP bSEXP, SEXP pivotSEXP, SEXP orderingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type pivot(pivotSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ordering(orderingSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_SLL(a, b, pivot, ordering));
    return rcpp_result_gen;
END_RCPP
}
// solve_LU
Eigen::MatrixXd solve_LU(const Eigen::Map <Eigen::MatrixXd> a, const Eigen::Map <Eigen::MatrixXd> b, unsigned int pivot);
RcppExport SEXP _sanic_solve_LU(SEXP aSEXP, SEXP bSEXP, SEXP pivotSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map <Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type pivot(pivotSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_LU(a, b, pivot));
    return rcpp_result_gen;
END_RCPP
}
// solve_SLU
Eigen::MatrixXd solve_SLU(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, unsigned int ordering);
RcppExport SEXP _sanic_solve_SLU(SEXP aSEXP, SEXP bSEXP, SEXP orderingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ordering(orderingSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_SLU(a, b, ordering));
    return rcpp_result_gen;
END_RCPP
}
// solve_HQR
Eigen::MatrixXd solve_HQR(const Eigen::Map<Eigen::MatrixXd> a, const Eigen::Map<Eigen::MatrixXd> b, unsigned int pivot);
RcppExport SEXP _sanic_solve_HQR(SEXP aSEXP, SEXP bSEXP, SEXP pivotSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type pivot(pivotSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_HQR(a, b, pivot));
    return rcpp_result_gen;
END_RCPP
}
// solve_SQR
Eigen::MatrixXd solve_SQR(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, unsigned int ordering);
RcppExport SEXP _sanic_solve_SQR(SEXP aSEXP, SEXP bSEXP, SEXP orderingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ordering(orderingSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_SQR(a, b, ordering));
    return rcpp_result_gen;
END_RCPP
}
// solve_BiCGSTAB
Eigen::MatrixXd solve_BiCGSTAB(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, const Eigen::Map<Eigen::MatrixXd> x0, double tol, unsigned int iter, unsigned int precond, bool verbose);
RcppExport SEXP _sanic_solve_BiCGSTAB(SEXP aSEXP, SEXP bSEXP, SEXP x0SEXP, SEXP tolSEXP, SEXP iterSEXP, SEXP precondSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type precond(precondSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_BiCGSTAB(a, b, x0, tol, iter, precond, verbose));
    return rcpp_result_gen;
END_RCPP
}
// solve_LSCG
Eigen::MatrixXd solve_LSCG(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, const Eigen::Map<Eigen::MatrixXd> x0, double tol, unsigned int iter, unsigned int precond, bool verbose);
RcppExport SEXP _sanic_solve_LSCG(SEXP aSEXP, SEXP bSEXP, SEXP x0SEXP, SEXP tolSEXP, SEXP iterSEXP, SEXP precondSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type precond(precondSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_LSCG(a, b, x0, tol, iter, precond, verbose));
    return rcpp_result_gen;
END_RCPP
}
// solve_CG
Eigen::MatrixXd solve_CG(const Eigen::MappedSparseMatrix<double> a, const Eigen::Map<Eigen::MatrixXd> b, const Eigen::Map<Eigen::MatrixXd> x0, double tol, unsigned int iter, unsigned int precond, bool verbose);
RcppExport SEXP _sanic_solve_CG(SEXP aSEXP, SEXP bSEXP, SEXP x0SEXP, SEXP tolSEXP, SEXP iterSEXP, SEXP precondSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MappedSparseMatrix<double> >::type a(aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type precond(precondSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(solve_CG(a, b, x0, tol, iter, precond, verbose));
    return rcpp_result_gen;
END_RCPP
}
// svd_J
Rcpp::List svd_J(const Eigen::Map<Eigen::MatrixXd> a, unsigned int type, unsigned int precond);
RcppExport SEXP _sanic_svd_J(SEXP aSEXP, SEXP typeSEXP, SEXP precondSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type precond(precondSEXP);
    rcpp_result_gen = Rcpp::wrap(svd_J(a, type, precond));
    return rcpp_result_gen;
END_RCPP
}
// svd_BDC
Rcpp::List svd_BDC(const Eigen::Map<Eigen::MatrixXd> a, unsigned int type);
RcppExport SEXP _sanic_svd_BDC(SEXP aSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type a(aSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(svd_BDC(a, type));
    return rcpp_result_gen;
END_RCPP
}
// is_symmetric_E
bool is_symmetric_E(const Eigen::Map<Eigen::MatrixXd> x, double tol);
RcppExport SEXP _sanic_is_symmetric_E(SEXP xSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(is_symmetric_E(x, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sanic_lanczos_E", (DL_FUNC) &_sanic_lanczos_E, 6},
    {"_sanic_arnoldi_E", (DL_FUNC) &_sanic_arnoldi_E, 5},
    {"_sanic_hessenberg_E", (DL_FUNC) &_sanic_hessenberg_E, 1},
    {"_sanic_tridiagonal_E", (DL_FUNC) &_sanic_tridiagonal_E, 1},
    {"_sanic_eigen_SA", (DL_FUNC) &_sanic_eigen_SA, 2},
    {"_sanic_eigen_SQ", (DL_FUNC) &_sanic_eigen_SQ, 2},
    {"_sanic_solve_LL", (DL_FUNC) &_sanic_solve_LL, 3},
    {"_sanic_solve_SLL", (DL_FUNC) &_sanic_solve_SLL, 4},
    {"_sanic_solve_LU", (DL_FUNC) &_sanic_solve_LU, 3},
    {"_sanic_solve_SLU", (DL_FUNC) &_sanic_solve_SLU, 3},
    {"_sanic_solve_HQR", (DL_FUNC) &_sanic_solve_HQR, 3},
    {"_sanic_solve_SQR", (DL_FUNC) &_sanic_solve_SQR, 3},
    {"_sanic_solve_BiCGSTAB", (DL_FUNC) &_sanic_solve_BiCGSTAB, 7},
    {"_sanic_solve_LSCG", (DL_FUNC) &_sanic_solve_LSCG, 7},
    {"_sanic_solve_CG", (DL_FUNC) &_sanic_solve_CG, 7},
    {"_sanic_svd_J", (DL_FUNC) &_sanic_svd_J, 3},
    {"_sanic_svd_BDC", (DL_FUNC) &_sanic_svd_BDC, 2},
    {"_sanic_is_symmetric_E", (DL_FUNC) &_sanic_is_symmetric_E, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_sanic(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
