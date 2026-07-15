#ifndef nlmixr2est_np_h
#define nlmixr2est_np_h
// Thin numeric interface between the FOCEI inner machinery (inner.cpp) and the
// nonparametric estimation kernels: NPAG (npag.cpp) and NPB (npb.cpp).  Like
// imp.h, the kernels work only with plain arma/double types and never see the
// internal focei_options / focei_ind structs, which stay private to inner.cpp.
#include <RcppArmadillo.h>

// ---- implemented in inner.cpp (shared with the imp interface) ----

// Conditional data log-likelihood log p(y_i | eta) for subject id, WITHOUT the
// Omega prior (the sum of the per-observation conditional log-densities from
// likInner0).  Returns -Inf on a non-finite (bad-solve) observation.  This is
// the single new numeric primitive the nonparametric engines need; a support
// point k evaluated for subject i gives psi(i,k) = exp(npEvalCondLik(eta_k, i)).
double npEvalCondLik(double *eta, int id);

// Build the Psi matrix (nSub x nPoint) on the already set-up FOCEi inner solve:
// psi(i,k) = p(y_i | support point k) where etaPoints is nPoint x neta.  Parallel
// over base subjects.  Requires vaeInnerSetup_ (or foceiSetup_) already run.
void npBuildPsiCore(const arma::mat& etaPoints, int cores, arma::mat& psi);

// Nonparametric adaptive-grid EM driver; called from foceiFitCpp_ when
// est=="npag" (in place of foceiOuter).
void npagOuter(Rcpp::Environment e);

// Nonparametric Bayes (stick-breaking DP) driver; called from foceiFitCpp_ when
// est=="npb" (in place of foceiOuter).
void npbOuter(Rcpp::Environment e);

#endif // nlmixr2est_np_h
