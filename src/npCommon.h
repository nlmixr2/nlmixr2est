#ifndef nlmixr2est_npCommon_h
#define nlmixr2est_npCommon_h
// Shared numerical routines for the nonparametric estimation engines (npag,
// npb): Burke's interior-point weight solver, and (later milestones) the Sobol
// initial grid, QR rank-revealing condensation, and adaptive grid expansion.
// These are plain-arma routines with no dependency on the private inner.cpp
// state, so they are unit-testable in isolation.
#include <RcppArmadillo.h>

// Burke's primal-dual interior point method (Yamada 2021, Appendix A).
//
// Given the non-negative likelihood matrix psi (n_sub x n_point) with
// psi(i,k) = p(y_i | support point k), solve the convex program
//   maximize   sum_i log( sum_k lambda_k psi(i,k) )
//   subject to lambda_k >= 0, sum_k lambda_k = 1
// for the support-point weights lambda.
//
// Returns the weight vector (length n_point, non-negative, sums to 1) and writes
// the objective value (the maximized log-likelihood) to *obj.  Non-finite psi
// entries are an error; negative entries are coerced to their absolute value
// (matching the reference implementation).  Throws (Rcpp stop) if the Newton
// system is not positive definite.
arma::vec npBurke(const arma::mat& psi, double* obj);

#endif // nlmixr2est_npCommon_h
