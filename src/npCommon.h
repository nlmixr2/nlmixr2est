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

// Sobol low-discrepancy initial grid: n support points over the box
// [lower, upper] (one row per point, ncol = length(lower) = neta).  Uses the
// same Boost Sobol engine + half-step scaling as the importance sampler.
arma::mat npSobolGrid(int n, const arma::vec& lower, const arma::vec& upper);

// Weight-threshold condensation (Yamada Alg 3): 0-based indices of the support
// points whose weight exceeds max(lambda) * ratio (ratio defaults to 1e-3).
arma::uvec npCondenseWeights(const arma::vec& lambda, double ratio);

// QR rank-revealing condensation: 0-based indices of a maximal set of
// linearly-independent support points (columns of psi), found by a
// column-pivoted QR of the row-normalized psi (drop columns whose relative
// R-diagonal is below tol, default 1e-8).  Mirrors pmcore qr.rs + the npag
// condensation keep-rule.
arma::uvec npCondenseQR(const arma::mat& psi, double tol);

// Adaptive-grid expansion (Yamada Alg 2): for each support point (row of theta)
// add up to 2*d daughters at +/- eps*(upper-lower) per dimension; discard
// daughters outside the box or closer than minDist (scaled L1) to any original
// point.  Returns the original points followed by the accepted daughters.
arma::mat npExpandGrid(const arma::mat& theta, double eps, const arma::vec& lower,
                       const arma::vec& upper, double minDist);

#endif // nlmixr2est_npCommon_h
