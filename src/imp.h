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

// 0.5 * log|Omega^-1| = -0.5 * log|Omega| (importance-sampling objective normalizer).
double impLogDetOmegaInv5();

// Maximum EM iterations (from the impmap control).
int impNiter();

// Omega diagonal parameterization ("sqrt"/"log"/"identity") for the EM Omega update.
std::string impDiagXform();

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

// Solve the theta-sensitivity model for subject `id` at `eta`, filling `dfdth`
// (nobs x nSens) with d(prediction)/d(theta) at each observation.  Column s maps
// to the non-mu structural theta impThetaSensIdx[s].  Returns false if the
// theta-sensitivity model is not available.
bool impThetaSensDfDtheta(int id, const arma::vec& eta, arma::mat& dfdth);

// Accumulate subject `id`'s IS-weighted score (into `g`, length nSens) and
// Gauss-Newton Hessian (into `H`, nSens x nSens) for the non-mu structural
// thetas, from its samples `S` (nsamp x neta) and normalized weights `zk`.
void impThetaScore(int id, const arma::mat& S, const arma::vec& zk,
                   arma::vec& g, arma::mat& H);

// Number of non-mu structural thetas (the length of impThetaSensIdx).
int impThetaSensN();
// Apply a Newton step: add `step` (length nSens) to the non-mu structural thetas
// and propagate to the full parameter vector.
void impUpdateStructThetas(const arma::vec& step);

// ---- implemented in imp.cpp ----

// Importance-sampling EM driver; called from foceiFitCpp_ when est=="impmap"
// (in place of foceiOuter).  Module M1: a single MAP pass plus per-subject
// mode/Hessian/likelihood collection into `e`.
void impOuter(Rcpp::Environment e);

#endif // nlmixr2est_imp_h
