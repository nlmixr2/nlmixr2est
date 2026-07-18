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

// Extended-least-squares residual objective at fixed per-subject etas (postEta,
// nsub x neta): sum over observations of (f-dv)^2/r + log(r), on the transform-both-
// sides scale, re-solving each subject at its eta.  The log(r) term keeps the residual
// scale from drifting to zero on a flexible support.  R_PosInf on a bad solve.
double npResidELS(const arma::mat& postEta);

// Pin the current solve for a residual-only optimization: cache each subject's
// (per-component) ODE states at the posterior etas so npResidELS can recompute r
// with the current residual params WITHOUT re-integrating (set op_focei.freezeOde
// around the resid optimizer).  Clear releases the cache.
void npResidFreezeBuild(const arma::mat& postEta);
void npResidFreezeClear();

// Per-endpoint moment residual sums at fixed per-subject etas: for each 0-based
// endpoint (obsEndpoint gives the endpoint of each observation in subject-major
// getIndIx order), the additive sum(err^2), the proportional sum((err/f)^2), and the
// observation count (nEnd x 3).  The saem-style estimate for warm-starting the
// residual optimization.
arma::mat npResidMoments(const arma::mat& postEta, const arma::ivec& obsEndpoint, int nEnd);

// Per-observation 0-based endpoint index (subject-major getIndIx order) from the cached
// CMT covariate (rxode2 getIndCmt): matches each observation's cmt to endpointCmt (the
// per-endpoint cmt values, predDf order).  All-zeros for a single-endpoint model.
arma::ivec npBuildObsEndpoint(const std::vector<int>& endpointCmt);

// Build the Psi matrix (nSub x nPoint) on the already set-up FOCEi inner solve:
// psi(i,k) = p(y_i | support point k) where etaPoints is nPoint x neta.  Parallel
// over base subjects.  Requires vaeInnerSetup_ (or foceiSetup_) already run.
void npBuildPsiCore(const arma::mat& etaPoints, int cores, arma::mat& psi);

// Unnormalized Psi at a residual-error multiplier gamma (single absolute scale,
// no per-row offset).  For the D(F) certificate's cross-matrix likelihood ratios.
void npBuildPsiCoreGamma(const arma::mat& etaPoints, int cores, double gamma,
                         arma::mat& psi);

// ODE-freeze cache for residual-parameter optimization: npFreezeBuild solves and
// caches each (support point, subject) once; npFreezePsiScaled rebuilds the
// (normalized) Psi from the cached states without re-integrating; npFreezeClear
// releases the cache.
void npFreezeBuild(const arma::mat& etaPoints, int cores);
void npFreezePsiScaled(const arma::mat& etaPoints, int cores, arma::mat& psi, double* offset);
void npFreezeClear();

// Build Psi at a residual-error multiplier gamma, reusing likInner0's full
// conditional likelihood so censoring and transform-both-sides are handled at the
// scaled error.  Row-normalized (log-sum-exp) for numerical stability; the true
// log-likelihood is burke_objf + *offset.
// rowMax (optional): the per-subject pre-normalization max log conditional
// likelihood over the support points.  A subject whose value is non-finite
// (-Inf) had zero density at every support point (degenerate); lets the caller
// detect that without a separate raw Psi build.
void npBuildPsiCoreScaled(const arma::mat& etaPoints, int cores, double gamma,
                          arma::mat& psi, double* offset, arma::vec* rowMax = nullptr);

// Subject-parallel conditional-likelihood contributions for the npb
// support-location MH step: for each physical subject, k = z[subject]; when
// occupied, computes the current/proposed conditional log-likelihoods at the
// cluster's current/proposed support eta (curLoc[k] / propLoc[k]) into
// curContrib/propContrib (indexed by physical subject id).  No RNG inside; the
// caller keeps the (serial) proposal + accept/reject draws.  Implemented in
// inner.cpp, used by npb.cpp.
void npbSupportMHContrib(const std::vector<int>& z, const std::vector<char>& occ,
                         std::vector<std::vector<double> >& curLoc,
                         std::vector<std::vector<double> >& propLoc,
                         std::vector<double>& curContrib,
                         std::vector<double>& propContrib);

// Shared fit finalization for the nonparametric engines: given the discrete
// mixing distribution (support points in eta space + weights) and the per-subject
// posterior-mean etas, summarize into the population theta shift + Omega, push
// into the FOCEi state (imp M-step helpers), build the fit env (impMapPass), and
// set e$objective = -2*objf.  Returns the (full) support-point covariance Omega
// so the caller can attach it.  Implemented in npag.cpp.
// gamma (the fitted assay-error multiplier, 1 if not optimized) is folded into
// the variance-scale residual thetas named by residScaleIdx (0-based fullTheta
// indices) so the reported residual reflects the estimate.
arma::mat npFinalizeFit(Rcpp::Environment e, arma::mat& support,
                        const arma::vec& weights, arma::mat postEta,
                        double objf, const arma::mat& omModel,
                        double gamma, const std::vector<int>& residScaleIdx,
                        const std::vector<int>& injEtaIdx,
                        const std::vector<int>& injThetaIdx);

// Residual/regressor theta optimization at the posterior-mean etas given the mixing
// distribution (support, weights): bounded bobyqa over the fullTheta indices in idx
// (kinds in kind) on the extended-least-squares objective sum((f-dv)^2/r + log(r)),
// warm-started from the saem-style per-endpoint moment.  optEnd[j]/optProp[j] tag
// idx[j] (0-based endpoint or -1; 1 if proportional); obsEndpoint gives each
// observation's endpoint in subject-major getIndIx order.  `freeze` is vestigial
// (ELS re-solves).  Leaves fullTheta at the optimum; returns the ELS value (R_PosInf
// on failure, thetas restored).  Implemented in npag.cpp, shared with npb.cpp.
double npOptimizeResid(const arma::mat& support, const arma::vec& weights,
                       const std::vector<int>& idx, const std::vector<int>& kind,
                       int cores, const std::vector<double>& lower,
                       const std::vector<double>& upper, bool freeze,
                       const arma::ivec& obsEndpoint = arma::ivec(),
                       const std::vector<int>& optEnd = std::vector<int>(),
                       const std::vector<int>& optProp = std::vector<int>(),
                       bool reDerive = false);

// Nonparametric adaptive-grid EM driver; called from foceiFitCpp_ when
// est=="npag" (in place of foceiOuter).
void npagOuter(Rcpp::Environment e);

// Nonparametric Bayes (stick-breaking DP) driver; called from foceiFitCpp_ when
// est=="npb" (in place of foceiOuter).
void npbOuter(Rcpp::Environment e);

#endif // nlmixr2est_np_h
