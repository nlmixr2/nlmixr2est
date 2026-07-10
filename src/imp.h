#ifndef nlmixr2est_imp_h
#define nlmixr2est_imp_h
// Thin numeric interface between the FOCEI inner machinery (inner.cpp) and the
// importance-sampling EM kernel (imp.cpp).  imp.cpp holds the IS/EM algorithm
// and works only with plain arma/double types -- it never sees the internal
// focei_options / focei_ind structs, which stay private to inner.cpp.
#include <RcppArmadillo.h>

// ---- implemented in inner.cpp (see the live op_focei / inds_focei state) ----

// Number of (base) subjects and number of etas in the current setup.
int impNsub();
int impNeta();

// Importance-sampling controls carried on op_focei (set in foceiSetup_ from the
// impmap control): samples per subject, proposal-variance inflation gamma, and
// the solve's OpenMP core count.
int impNsample();
double impGammaProp();
int impCores();
int impSetSolveCores(int cores);   // override solve cores (returns previous); forces mixture EM serial

// 0.5 * log|Omega^-1| = -0.5 * log|Omega| (importance-sampling objective normalizer).
double impLogDetOmegaInv5();

// Maximum EM iterations (from the impmap control).
int impNiter();

// Omega diagonal parameterization ("sqrt"/"log"/"identity") for the EM Omega update.
std::string impDiagXform();

// Convergence controller / proposal-scale adaptation controls (from impmapControl):
double impIaccept();      // target effective-sample fraction that gamma adapts toward
double impIscaleMin();    // lower bound for the adapted gamma
double impIscaleMax();    // upper bound for the adapted gamma
int impNconvWindow();     // trailing-iteration window for the convergence check
double impCtol();         // relative windowed-convergence tolerance (derived from sigdig if unset)

int impMuGroupN();                    // number of mu-referenced covariate groups (diagnostic)

// M-step helpers (the EM loop is in impOuter):
void impSetEta(int id, const arma::vec& eta);      // overwrite subject id's eta
void impGetEta(int id, arma::vec& eta);            // read subject id's eta
void impGetOmega(arma::mat& Om);                   // current Omega (for its zero pattern)
double impUpdateMuThetas();                        // mu-referenced covariate regression (updateMuGroups)
void impMuInterceptStep();                         // simple mu intercept EM update (no covariates)
void impReMap();                                   // re-optimize all conditional modes (innerOpt)
void impSetOmega(const arma::mat& Omega, const std::string& diagXform); // install new Omega
void impSyncInitParToFullTheta();                  // sync optimizer reference to converged fullTheta
void impGetEstPar(arma::vec& par);                 // current estimated free-parameter vector (EM convergence)

// Run a single MAP pass over all subjects at the initial parameters (reuses the
// FOCEI posthoc path, foceiOuterFinal) and populate the fit environment `e`.
void impMapPass(Rcpp::Environment e);

// Copy subject `id`'s current MAP mode (eta, length neta) into `mode`.
void impGetMode(int id, arma::vec& mode);

// Subject `id`'s individual objective contribution at its current mode.
double impGetIndLik(int id);

// Eta Hessian of the negative inner joint objective at subject `id`'s mode
// (the FOCEI Laplace information matrix; positive-definite at the mode).
// Returns false if the solve/Hessian could not be formed.
bool impGetHessian(int id, arma::mat& H);

// Joint log-density log(pi_i(eta)) = log(l(y_i|phi,theta) * h(eta|Omega)) at an
// arbitrary eta for subject `id` (the importance-sampling weight numerator).
double impEvalJointLik(const arma::vec& eta, int id);

// Accumulate subject `id`'s IS-weighted score (into `g`, length nSens) and
// Gauss-Newton Hessian (into `H`, nSens x nSens) for the non-mu structural
// thetas, from its samples `S` (nsamp x neta) and normalized weights `zk`.
void impThetaScore(int id, const arma::mat& S, const arma::vec& zk,
                   arma::vec& g, arma::mat& H);

// Number of non-mu structural thetas (the length of impThetaSensIdx).
int impThetaSensN();

// 0-based eta indices whose Omega diagonal is fixed (held across the EM update).
void impGetOmegaFixedEta(std::vector<int>& idx);

// ---- mixture (sub-population) support ----
int impNmix();                                     // number of components (1 if none)
double impMixProb(int j);                          // population proportion of component j (0-based)
void impSetMixThetas(const arma::vec& theta);      // install absolute $MIX thetas (mean-posterior EM) + recompute proportions

// ---- Monte-Carlo covariance support (implemented in inner.cpp) ----
int impNtheta();                                   // number of thetas
bool impCovEnabled();                              // whether impCov=TRUE was requested
void impGetEstThetaIdx(std::vector<int>& idx);     // fullTheta indices of the estimated thetas
void impGetCovParList(std::vector<int>& idx);      // fullTheta index of every free param (fixedTrans order)
double impGetFullThetaVal(int idx);                // current value of fullTheta[idx]
void impSetThetaAll(int idx, double val);          // set fullTheta[idx] on every subject (FD perturb)
void impForceResolve(int id);                      // force likInner0 to re-solve subject id
int impOmegaN();                                   // number of parameterized Omega free parameters
double impGetOmegaThetaVal(int m);                 // current value of Omega free parameter m
void impSetOmegaThetaAll(int m, double val);       // set Omega free param m + rebuild omegaInv/logdet

// Clear the persistent inner neqOverride (multi-endpoint pool) at fit end so it
// does not leak into a subsequent fit sharing the global solve context.
void impClearInnerNeqOverride();
// Apply a Newton step: add `step` (length nSens) to the non-mu structural thetas
// and propagate to the full parameter vector.
void impUpdateStructThetas(const arma::vec& step);

// ---- implemented in imp.cpp ----

// Monte-Carlo observed-information covariance for the estimated thetas: FD
// Hessian of the importance-sampling -2LL over fixed (common-random-number)
// samples.  Stashes impCovTheta / impSeTheta / impCovThetaIdx on `e`.
void impComputeCov(Rcpp::Environment e);

// Importance-sampling EM driver; called from foceiFitCpp_ when est=="impmap"
// (in place of foceiOuter).  Module M1: a single MAP pass plus per-subject
// mode/Hessian/likelihood collection into `e`.
void impOuter(Rcpp::Environment e);

#endif // nlmixr2est_imp_h
