// NPAG -- Nonparametric Adaptive Grid (Yamada 2021).  The estimation loop lives
// entirely here and reuses the FOCEI inner likelihood machinery in inner.cpp
// (via the np.h / imp.h numeric interfaces) to fill the Psi matrix of per-subject
// conditional likelihoods at each support point.  Called from foceiFitCpp_ in
// place of foceiOuter when est=="npag".
#include <RcppArmadillo.h>
#include "np.h"
#include "npCommon.h"
#include "imp.h"   // M-step helpers (impSetEta/impSetOmega/impMapPass/...) reused for finalization

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Residual-error theta optimization (per user: gamma is only a warm start).
// With the support points and weights FIXED, optimize the residual thetas
// (add/prop/lnorm/lambda/ar -- every endpoint) against the nonparametric -2LL
// with the bounded bobyqa optimizer.  This recovers what a single global gamma
// cannot: the add/prop ratio and per-endpoint magnitudes.
// Posterior-mean support eta of each subject given the current mixing distribution
// (support nspp x neta, weights nspp): E[eta | y_i] = sum_k w_k psi_ik phi_k /
// sum_k w_k psi_ik -- the individual random-effect estimate the residual step uses.
static arma::mat npPosteriorEta(const arma::mat& support, const arma::vec& weights,
                                int cores) {
  arma::mat psi;
  npBuildPsiCore(support, cores, psi);          // nsub x nspp
  arma::mat wpsi = psi;
  wpsi.each_row() %= weights.t();               // w_k psi_ik
  arma::vec denom = arma::sum(wpsi, 1);
  denom = arma::clamp(denom, 1e-300, arma::datum::inf);
  arma::mat num = wpsi * support;               // nsub x neta
  num.each_col() /= denom;
  return num;
}

// Context is passed to the R-callable objective through file-static pointers.  The
// residual (err) thetas and any structural "regressor" theta are optimized together
// against ONE extended-least-squares objective sum_obs((f-dv)^2/r + log(r)) at the
// posterior-mean etas.  The log(r) term keeps the residual scale from drifting to zero
// on a flexible support (the marginal likelihood rewards a vanishing residual because
// each support point can then fit its subjects arbitrarily well).  When the set
// includes a regressor (which shifts the ODE predictions), gNpReDerive re-derives the
// posterior-mean etas at each candidate so the eta grid cannot stale-absorb the
// structural shift -- that is what identifies the regressor.
namespace {
  std::vector<int> gNpOptIdx;              // fullTheta indices being optimized
  std::vector<int> gNpOptKind;            // 0 = free, 1 = positive (SD), 2 = corr
  arma::mat gNpPostEta;                    // per-subject etas for the ELS objective
  const arma::mat* gNpSupport = nullptr;   // support grid (for re-deriving the etas)
  const arma::vec* gNpWeights = nullptr;   // weights (for re-deriving the etas)
  int gNpCores = 1;
  bool gNpReDerive = false;                // re-derive the etas each candidate (regressor)

  // clamp a candidate to the parameter's natural range (belt-and-suspenders with
  // the bobyqa bounds): SD >= 0, continuous-AR correlation in (0,1).
  double npResidVal(double v, int kind) {
    if (kind == 1) return std::fabs(v);                                 // SD >= 0
    if (kind == 2) return std::max(0.0, std::min(0.999, v));            // AR corr in (0,1)
    return v;                                                           // lambda etc.
  }

  double npResidObjVal(const double *p) {
    for (size_t j = 0; j < gNpOptIdx.size(); ++j)
      impSetThetaAll(gNpOptIdx[j], npResidVal(p[j], gNpOptKind[j]));
    if (gNpReDerive)
      gNpPostEta = npPosteriorEta(*gNpSupport, *gNpWeights, gNpCores);
    return npResidELS(gNpPostEta);
  }
  double npResidObjR(Rcpp::NumericVector p) { return npResidObjVal(p.begin()); }
}

// Optimize the residual thetas in idx (kinds in kind) at the posterior-mean etas
// given (support, weights), using bobyqa BOUNDED by the ini-block lower/upper on the
// extended-least-squares objective.  Warm-starts every variance-scale theta from the
// saem-style per-endpoint moment (additive: sqrt(mean err^2); proportional:
// sqrt(mean (err/f)^2)); when the optimized set is exactly one such scale per
// endpoint, the moment IS the ELS optimum, so it is used directly (no bobyqa).
// optEnd[j]/optProp[j] describe idx[j]: optEnd 0-based endpoint (or -1 if not a
// variance scale), optProp 1 for proportional else 0.  obsEndpoint gives each
// observation's endpoint (subject-major getIndIx order) for the moment estimate.
// `freeze` is now vestigial (ELS always re-solves).  Leaves fullTheta at the
// optimum; returns the ELS value (R_PosInf if bobyqa failed, thetas left at start).
double npOptimizeResid(const arma::mat& support, const arma::vec& weights,
                       const std::vector<int>& idx,
                       const std::vector<int>& kind, int cores,
                       const std::vector<double>& lower,
                       const std::vector<double>& upper,
                       bool freeze,
                       const arma::ivec& obsEndpoint,
                       const std::vector<int>& optEnd,
                       const std::vector<int>& optProp,
                       bool reDerive) {
  (void) freeze;
  int n = (int)idx.size();
  if (n == 0) return R_PosInf;
  gNpOptIdx = idx; gNpOptKind = kind;
  gNpSupport = &support; gNpWeights = &weights; gNpCores = cores;
  gNpReDerive = reDerive;
  gNpPostEta = npPosteriorEta(support, weights, cores);

  // saem-style moment warm start: per-endpoint additive / proportional SD at the
  // individual predictions.  Set each variance-scale theta to its moment.
  bool haveWarm = !optEnd.empty() && (int)optEnd.size() == n;
  bool allSimpleScale = haveWarm && !reDerive;   // a lone scale per endpoint, no regressor
  if (haveWarm) {
    int nEnd = 0;
    for (int j = 0; j < n; ++j) nEnd = std::max(nEnd, optEnd[j] + 1);
    for (int j = 0; j < (int)obsEndpoint.n_elem; ++j) nEnd = std::max(nEnd, (int)obsEndpoint[j] + 1);
    arma::mat mom = npResidMoments(gNpPostEta, obsEndpoint, nEnd);
    std::vector<int> endCount(std::max(nEnd, 1), 0);
    for (int j = 0; j < n; ++j) if (optEnd[j] >= 0) endCount[optEnd[j]]++;
    for (int j = 0; j < n; ++j) {
      int e = optEnd[j];
      if (e < 0 || kind[j] != 1 || endCount[e] > 1) { allSimpleScale = false; continue; }
      double nn = mom(e, 2);
      if (nn <= 0.0) { allSimpleScale = false; continue; }
      double v = std::sqrt(mom(e, (optProp[j] == 1) ? 1 : 0) / nn);
      if (j < (int)lower.size() && std::isfinite(lower[j])) v = std::max(v, lower[j]);
      if (j < (int)upper.size() && std::isfinite(upper[j])) v = std::min(v, upper[j]);
      if (std::isfinite(v) && v > 0.0) impSetThetaAll(idx[j], v);
    }
    // simple path (one scale per endpoint, no regressor): the moment IS the estimate.
    if (allSimpleScale) return npResidELS(gNpPostEta);
  }

  std::vector<double> start(n);
  Rcpp::NumericVector par0(n), lo(n), hi(n);
  for (int j = 0; j < n; ++j) {
    double v = impGetFullThetaVal(idx[j]);
    if (kind[j] == 1) v = std::fabs(v);
    start[j] = v; par0[j] = v;
    lo[j] = (j < (int)lower.size()) ? lower[j] : R_NegInf;
    hi[j] = (j < (int)upper.size()) ? upper[j] : R_PosInf;
    // an unbounded regressor (kind 0, no ini bound) must not let bobyqa wander to an
    // extreme value where the profiled fit is spuriously flat; bracket it generously
    // around the start (in the parameter's own -- often log -- scale).
    if (kind[j] == 0) {
      if (!std::isfinite(lo[j])) lo[j] = v - 6.0;
      if (!std::isfinite(hi[j])) hi[j] = v + 6.0;
    }
  }
  // .boundedResidOpt: bobyqa (>=2 params) / bounded optimize (1 param), honoring
  // the ini-block bounds so the optimizer stays in the valid region.
  Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
  Rcpp::Function boundedOpt = nlmixr2[".boundedResidOpt"];
  Rcpp::InternalFunction fn(&npResidObjR);
  Rcpp::List ret = boundedOpt(Rcpp::_["par"] = par0, Rcpp::_["fn"] = fn,
                              Rcpp::_["lower"] = lo, Rcpp::_["upper"] = hi,
                              Rcpp::_["control"] = Rcpp::List::create(Rcpp::_["maxfun"] = 200 * n * n));
  double f = as<double>(ret["value"]);
  gNpReDerive = false;
  if (!ISNA(f)) {
    Rcpp::NumericVector x = ret["x"];
    for (int j = 0; j < n; ++j) impSetThetaAll(idx[j], npResidVal(x[j], kind[j]));
    return f;
  }
  for (int j = 0; j < n; ++j) impSetThetaAll(idx[j], npResidVal(start[j], kind[j]));
  return R_PosInf;
}

// NPAG tuning constants (Yamada 2021): eps convergence floor / objf tolerance /
// F-candidate tolerance / minimum scaled support-point distance.
struct npagCtl {
  int points = 2028;    // initial Sobol grid size
  int cycles = 100;     // maximum cycles
  int cores = 1;
  double ratio = 1e-3;  // weight-threshold condensation (max*ratio)
  double qrTol = 1e-8;  // QR rank-revealing tolerance
  double epsInit = 0.2; // initial adaptive-grid fraction
  double thetaE = 1e-4; // eps floor
  double thetaG = 1e-4; // objf-change tolerance that halves eps
  double thetaF = 1e-2; // successive-F tolerance for the exit test
  double thetaD = 1e-4; // minimum scaled distance for daughter points
  bool gammaOptimize = false; // optimize the residual-error magnitude (gamma)
  double gammaInit = 1.0;     // initial gamma multiplier
  double gammaDelta = 0.1;    // initial gamma step fraction
  bool trace = false;         // record the per-cycle parameter history (scale.h)
  arma::mat omModel;          // model Omega sparsity (for the trace push mask)
  int residMode = 1;          // residual theta opt: 0 none, 1 alternate, 2 final
  std::vector<int> residOptIdx;  // fullTheta indices of the non-fixed residual (err) thetas
  std::vector<int> residOptKind; // per-idx: 0 free (lambda), 1 SD (>0), 2 corr
  std::vector<double> residOptLower; // ini-block lower bounds (bobyqa)
  std::vector<double> residOptUpper; // ini-block upper bounds (bobyqa)
  std::vector<int> residScaleIdx; // variance-scale subset (gamma warm-start fold)
  std::vector<int> residOptEnd;   // per-idx 0-based endpoint (-1 if not a variance scale)
  std::vector<int> residOptProp;  // per-idx 1 if proportional, else 0 (moment warm start)
  arma::ivec obsEndpoint;         // per-observation endpoint (subject-major getIndIx order)
  std::vector<int> regressIdx;    // fullTheta indices of structural regressor thetas
  std::vector<double> regressLower;
  std::vector<double> regressUpper;
  bool mixOptimize = false;   // estimate mix() proportions via the in-cycle EM update
  bool residFreeze = true;    // freeze the ODE during resid opt (safe: err params only)
  std::vector<int> muExpandEtaIdx;   // 0-based eta indices of mu-expanded (injected) etas
  std::vector<int> muExpandThetaIdx; // paired 0-based theta (fullTheta) indices
};

// Support covariance installed into op_focei, masked by the model's Omega
// sparsity (diagonal always + modeled off-diagonals), with fixed-Omega diagonals
// restored to the model (fixed) value -- those etas are pinned to 0 in the grid
// so their support-point variance is 0.
static arma::mat npMaskedOmega(const arma::mat& Omega, const arma::mat& omModel) {
  int neta = (int)Omega.n_rows;
  arma::mat out(neta, neta, arma::fill::zeros);
  bool haveModel = ((int)omModel.n_rows == neta && (int)omModel.n_cols == neta);
  for (int i = 0; i < neta; ++i)
    for (int j = 0; j < neta; ++j)
      if (i == j || (haveModel && omModel(i, j) != 0.0)) out(i, j) = Omega(i, j);
  std::vector<int> fixedEta;
  impGetOmegaFixedEta(fixedEta);
  for (size_t f = 0; f < fixedEta.size(); ++f) {
    int fi = fixedEta[f];
    if (fi >= 0 && fi < neta && haveModel) out(fi, fi) = omModel(fi, fi);
  }
  return out;
}

// The NPAG adaptive-grid cycle (Yamada Alg 1).  Requires the FOCEi inner solve
// already set up (vaeInnerSetup_ / foceiSetup_), so npBuildPsiCore can fill Psi.
// Returns the final support points (eta space), their weights, the log-
// likelihood, cycle count, and a converged flag.
//
// objf is the nonparametric marginal log-likelihood
//   sum_i log(sum_k lambda_k p(y_i | support point k))
// using the inner conditional-density constant convention (and, in the gamma
// path, a per-subject log-sum-exp offset).  This is NOT the NONMEM/FOCEI OFV
// convention, so the resulting -2LL is not directly comparable to nlmixr2's
// FOCEI/SAEM/FOCE -2LL -- compare NPAG to NPAG (or to Pmetrics NPAG) instead.
struct npagResult {
  arma::mat theta;   // nspp x neta support points (eta space)
  arma::vec lambda;  // nspp weights (sum 1)
  arma::mat psi;     // nsub x nspp final conditional-likelihood matrix
  double objf;       // maximized log-likelihood
  double gamma;      // final residual-error multiplier (1 if not optimized)
  int cycle;
  bool converged;
};

static npagResult npagRunCycle(const arma::vec& lower, const arma::vec& upper,
                               const npagCtl& ctl) {
  arma::mat theta = npSobolGrid(ctl.points, lower, upper);
  double eps = ctl.epsInit, f0 = -1e30, f1 = 0.0, lastObj = -1e30, objf = R_NegInf;
  // gamma is held at 1: the per-cycle warm-start folds its multiplier straight
  // into the variance-scale thetas, so the main/condensation Psi builds and the
  // finalization all see gamma == 1.
  double gamma = 1.0, gDelta = ctl.gammaDelta;
  int cycle = 0;
  bool converged = false;
  arma::vec lam;
  arma::mat psi;
  // The residual step uses extended least squares (which does not collapse) at the
  // posterior-mean etas.  A structural regressor is appended to the SAME step (kind 0,
  // no endpoint) and identified by re-deriving those etas per candidate so the eta
  // grid cannot stale-absorb the structural shift.
  bool useResidOpt = !ctl.residOptIdx.empty();
  bool useRegress = !ctl.regressIdx.empty();
  bool doResidOpt = useResidOpt || useRegress;
  // combined optimized set: residual (err) thetas, then the regressor thetas.
  std::vector<int> optIdx = ctl.residOptIdx, optKind = ctl.residOptKind;
  std::vector<int> optEnd = ctl.residOptEnd, optProp = ctl.residOptProp;
  std::vector<double> optLo = ctl.residOptLower, optHi = ctl.residOptUpper;
  for (size_t g = 0; g < ctl.regressIdx.size(); ++g) {
    optIdx.push_back(ctl.regressIdx[g]);
    optKind.push_back(0);                       // free
    optEnd.push_back(-1);                       // not a variance scale
    optProp.push_back(0);
    optLo.push_back(g < ctl.regressLower.size() ? ctl.regressLower[g] : R_NegInf);
    optHi.push_back(g < ctl.regressUpper.size() ? ctl.regressUpper[g] : R_PosInf);
  }
  for (;;) {
    cycle++;
    if (cycle > 1) {
      theta = npExpandGrid(theta, eps, lower, upper, ctl.thetaD);
    }
    // estimation: Psi at the current gamma, then Burke weights.  The gamma path
    // builds Psi through likInner0 (npBuildPsiCoreScaled) so censoring and
    // transform-both-sides are handled at the scaled residual error.
    double obj0 = 0.0, offset = 0.0;
    if (ctl.gammaOptimize) {
      npBuildPsiCoreScaled(theta, ctl.cores, gamma, psi, &offset);
    } else {
      // per-row log-sum-exp normalized (the removed maxima are returned in offset)
      // so a hard subject's conditional density cannot underflow a whole Psi row to
      // zero -- which would abort condensation (npCondenseQR).  Burke's weights are
      // invariant to per-row scaling; the objective adds offset back.
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &offset);
    }
    // a subject whose conditional density is 0 at EVERY support point makes the
    // Burke IPM (and condensation) degenerate -- report it clearly.  The usual cause
    // is a transform-both-sides link evaluated where the structural prediction is
    // non-positive (e.g. lnorm/log at an observation whose model prediction is 0).
    // Check on the RAW (unnormalized) build: the per-row log-sum-exp normalization of
    // the working psi masks a zero row (exp(-Inf) becomes a NaN the offset absorbs),
    // so a raw exp(condLik) is what reliably exposes the degeneracy.  One extra Psi
    // build, only at cycle 1.
    if (cycle == 1) {
      arma::mat psiRaw;
      npBuildPsiCore(theta, ctl.cores, psiRaw);
      arma::vec rowSum = arma::sum(psiRaw, 1);
      arma::uvec bad = arma::find_nonfinite(rowSum);
      if (bad.is_empty()) bad = arma::find(rowSum <= 0.0);
      if (!bad.is_empty()) {
        throw std::runtime_error(
          "npag/npb: " + std::to_string(bad.n_elem) + " subject(s) (first at index " +
          std::to_string((int)bad[0] + 1) + ") have zero conditional density at every "
          "support point.  This usually means a transform-both-sides error model was "
          "evaluated where the model prediction is non-positive (e.g. an 'lnorm' or log "
          "link at an observation whose prediction is 0).  Check the residual error "
          "model and the observations near a zero prediction.");
      }
    }
    lam = npBurke(psi, &obj0);
    // condensation: weight threshold, then QR rank-revealing
    arma::uvec wk = npCondenseWeights(lam, ctl.ratio);
    theta = theta.rows(wk); psi = psi.cols(wk);
    arma::uvec qk = npCondenseQR(psi, ctl.qrTol);
    theta = theta.rows(qk); psi = psi.cols(qk);
    // re-solve the weights on the condensed grid
    if (ctl.gammaOptimize) {
      npBuildPsiCoreScaled(theta, ctl.cores, gamma, psi, &offset);
      double bObj = 0.0; lam = npBurke(psi, &bObj); objf = bObj + offset;
    } else {
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &offset);
      double bObj = 0.0; lam = npBurke(psi, &bObj); objf = bObj + offset;
    }
    // residual-error magnitude warm-start (gamma): a fast global up/down search
    // from 1 on the condensed grid, then FOLD the winning multiplier into the
    // variance-scale thetas so the thetas carry the estimate (finalization needs
    // no fold; gamma stays effectively 1).  This only warm-starts the overall
    // magnitude -- the per-endpoint values and add/prop ratio come from the
    // Nelder-Mead step below.
    if (ctl.gammaOptimize) {
      double gUp = 1.0 * (1.0 + gDelta), gDn = 1.0 / (1.0 + gDelta);
      arma::mat psiU, psiD;
      double offU = 0.0, offD = 0.0, bU = 0.0, bD = 0.0;
      npBuildPsiCoreScaled(theta, ctl.cores, gUp, psiU, &offU);
      arma::vec lamU = npBurke(psiU, &bU); double objU = bU + offU;
      npBuildPsiCoreScaled(theta, ctl.cores, gDn, psiD, &offD);
      arma::vec lamD = npBurke(psiD, &bD); double objD = bD + offD;
      double gsel = 1.0;
      if (objU > objf) { gsel = gUp; objf = objU; lam = lamU; psi = psiU; gDelta *= 4.0; }
      if (objD > objf) { gsel = gDn; objf = objD; lam = lamD; psi = psiD; gDelta *= 4.0; }
      gDelta *= 0.5;
      if (gDelta <= 0.01) gDelta = 0.1;
      if (gDelta > 1.0) gDelta = 1.0;
      for (size_t s = 0; s < ctl.residScaleIdx.size(); ++s) {
        int idx = ctl.residScaleIdx[s];
        impSetThetaAll(idx, impGetFullThetaVal(idx) * gsel);
      }
    }
    // residual + regressor optimization (extended least squares at the posterior-mean
    // etas; the etas are re-derived per candidate when a regressor is present so the
    // structural shift is identified).  "alternate" runs it every cycle (block-
    // coordinate ascent), then re-solves the weights at the new thetas.
    if (ctl.residMode == 1 && doResidOpt) {
      npOptimizeResid(theta, lam, optIdx, optKind, ctl.cores, optLo, optHi,
                      ctl.residFreeze, ctl.obsEndpoint, optEnd, optProp, useRegress);
      double off = 0.0, b = 0.0;
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
      lam = npBurke(psi, &b); objf = b + off;
    }
    // mixture-proportion EM update (support points + weights fixed): reweights the
    // components, then rebuild the marginalized Psi + Burke weights so the objf
    // reflects the new proportions -- block-coordinate ascent with the resid step.
    if (ctl.mixOptimize && impNmix() > 1) {
      npMixEMUpdate(theta, lam, ctl.cores);
      if (ctl.gammaOptimize) {
        double off = 0.0, b = 0.0;
        npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
        lam = npBurke(psi, &b); objf = b + off;
      } else {
        npBuildPsiCore(theta, ctl.cores, psi);
        lam = npBurke(psi, &objf);
      }
    }
    // per-cycle parameter-history row (shared scale.h printer).  Pushing Omega to
    // op_focei is safe here -- npag sums llikObs, which does not use omegaInv --
    // so impGetEstPar reflects the current support-point variance while the
    // objective is the nonparametric -2LL.
    if (ctl.trace) {
      arma::rowvec meanEta = lam.t() * theta;
      arma::mat Om(theta.n_cols, theta.n_cols, arma::fill::zeros);
      for (int k = 0; k < (int)theta.n_rows; ++k) {
        arma::rowvec d = theta.row(k) - meanEta; Om += lam[k] * (d.t() * d);
      }
      impSetOmega(npMaskedOmega(Om, ctl.omModel), impDiagXform());
      arma::vec par; impGetEstPar(par);
      impIterPrintRow(par, -2.0 * objf);
    }
    // evaluation (Yamada Alg 1 convergence controller); objf is the true
    // (offset-corrected) log-likelihood in both paths.
    if (std::fabs(lastObj - objf) <= ctl.thetaG && eps > ctl.thetaE) {
      eps /= 2.0;
      if (eps <= ctl.thetaE) {
        f1 = objf;
        if (std::fabs(f1 - f0) <= ctl.thetaF) { converged = true; break; }
        f0 = f1; eps = ctl.epsInit;
      }
    }
    if (cycle >= ctl.cycles) break;
    lastObj = objf;
  }
  // "final" mode: optimize the residual + regressor thetas once at the converged support.
  if (ctl.residMode == 2 && doResidOpt) {
    npOptimizeResid(theta, lam, optIdx, optKind, ctl.cores, optLo, optHi,
                    ctl.residFreeze, ctl.obsEndpoint, optEnd, optProp, useRegress);
    double off = 0.0, b = 0.0;
    npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
    lam = npBurke(psi, &b); objf = b + off;
  }
  // final support refinement: with the residual + regressor thetas now held CONSTANT,
  // re-optimize the support (adaptive grid + Burke) against the marginal likelihood.
  // The saem-style ELS residual is larger than the marginal-likelihood-optimal residual,
  // which leaves the (sparse) support marginal-LL-suboptimal -- refining it at the fixed
  // residual restores the nonparametric MLE of the mixing distribution and the D(F)
  // certificate (~0) without changing the reported residual or any structural parameter.
  if (doResidOpt) {
    double eps2 = ctl.epsInit, g0 = -1e30, g1 = 0.0, last2 = -1e30;
    for (int rc = 0; rc < ctl.cycles; ++rc) {
      theta = npExpandGrid(theta, eps2, lower, upper, ctl.thetaD);
      double off = 0.0;
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
      double b0 = 0.0; lam = npBurke(psi, &b0);
      arma::uvec wk = npCondenseWeights(lam, ctl.ratio);
      theta = theta.rows(wk); psi = psi.cols(wk);
      arma::uvec qk = npCondenseQR(psi, ctl.qrTol);
      theta = theta.rows(qk); psi = psi.cols(qk);
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
      double b = 0.0; lam = npBurke(psi, &b); objf = b + off;
      if (std::fabs(last2 - objf) <= ctl.thetaG && eps2 > ctl.thetaE) {
        eps2 /= 2.0;
        if (eps2 <= ctl.thetaE) {
          g1 = objf;
          if (std::fabs(g1 - g0) <= ctl.thetaF) break;
          g0 = g1; eps2 = ctl.epsInit;
        }
      }
      last2 = objf;
    }
  }
  npagResult r;
  r.theta = theta; r.lambda = lam; r.psi = psi; r.objf = objf; r.gamma = gamma;
  r.cycle = cycle; r.converged = converged;
  return r;
}

// est="npag" driver.  Runs the adaptive-grid cycle, summarizes the discrete
// mixing distribution into a population mean (mu-referenced theta shift) +
// covariance (Omega) and per-subject posterior-mean etas, pushes them into the
// FOCEi state via the shared imp M-step helpers, builds the fit environment
// (impMapPass -> foceiOuterFinal), then overrides the objective with the
// nonparametric -2LL and attaches the support-point distribution.  Called from
// foceiFitCpp_ in place of foceiOuter.
void npagOuter(Environment e) {
  List control = e["control"];
  arma::vec lower = as<arma::vec>(control["npBoxLower"]);
  arma::vec upper = as<arma::vec>(control["npBoxUpper"]);
  npagCtl ctl;
  ctl.points = as<int>(control["npPoints"]);
  ctl.cycles = as<int>(control["npCycles"]);
  ctl.cores = impCores();
  ctl.gammaOptimize = as<bool>(control["npGammaOptimize"]);
  ctl.trace = true;
  if (control.containsElementNamed("npResidScaleIdx")) {
    IntegerVector ri = control["npResidScaleIdx"];
    ctl.residScaleIdx.assign(ri.begin(), ri.end());
  }
  if (control.containsElementNamed("npResidOptIdx")) {
    IntegerVector ro = control["npResidOptIdx"];
    ctl.residOptIdx.assign(ro.begin(), ro.end());
  }
  if (control.containsElementNamed("npResidOptKind")) {
    IntegerVector rk = control["npResidOptKind"];
    ctl.residOptKind.assign(rk.begin(), rk.end());
  }
  if (control.containsElementNamed("npResidMode"))
    ctl.residMode = as<int>(control["npResidMode"]);
  if (control.containsElementNamed("npMixOptimize"))
    ctl.mixOptimize = as<bool>(control["npMixOptimize"]);
  if (control.containsElementNamed("npResidFreeze"))
    ctl.residFreeze = as<bool>(control["npResidFreeze"]);
  if (control.containsElementNamed("npMuExpandEtaIdx")) {
    IntegerVector me = control["npMuExpandEtaIdx"];
    ctl.muExpandEtaIdx.assign(me.begin(), me.end());
  }
  if (control.containsElementNamed("npMuExpandThetaIdx")) {
    IntegerVector mt = control["npMuExpandThetaIdx"];
    ctl.muExpandThetaIdx.assign(mt.begin(), mt.end());
  }
  if (control.containsElementNamed("npResidOptLower")) {
    NumericVector rl = control["npResidOptLower"];
    ctl.residOptLower.assign(rl.begin(), rl.end());
  }
  if (control.containsElementNamed("npResidOptUpper")) {
    NumericVector ru = control["npResidOptUpper"];
    ctl.residOptUpper.assign(ru.begin(), ru.end());
  }
  if (control.containsElementNamed("npResidOptEnd")) {
    IntegerVector re = control["npResidOptEnd"];
    ctl.residOptEnd.assign(re.begin(), re.end());
  }
  if (control.containsElementNamed("npResidOptProp")) {
    IntegerVector rp = control["npResidOptProp"];
    ctl.residOptProp.assign(rp.begin(), rp.end());
  }
  if (control.containsElementNamed("npObsEndpoint")) {
    IntegerVector oe = control["npObsEndpoint"];
    ctl.obsEndpoint.set_size(oe.size());
    for (int j = 0; j < oe.size(); ++j) ctl.obsEndpoint[j] = oe[j];
  }
  if (control.containsElementNamed("npRegressIdx")) {
    IntegerVector gi = control["npRegressIdx"];
    ctl.regressIdx.assign(gi.begin(), gi.end());
  }
  if (control.containsElementNamed("npRegressLower")) {
    NumericVector gl = control["npRegressLower"];
    ctl.regressLower.assign(gl.begin(), gl.end());
  }
  if (control.containsElementNamed("npRegressUpper")) {
    NumericVector gu = control["npRegressUpper"];
    ctl.regressUpper.assign(gu.begin(), gu.end());
  }
  std::vector<int> residScaleIdx = ctl.residScaleIdx;  // for the finalization fold (no-op at gamma=1)
  // foceiSetup_ defers building op_focei.omegaInv on the maxOuterIterations=0
  // path; likInner0 still evaluates the Omega-prior term (omegaInv * eta), so
  // rebuild the inverse from the initial Omega before the cycle calls likInner0.
  // (npag ignores the prior -- it sums llikObs -- but likInner0 always forms it.)
  arma::mat omModel;            // model's initial Omega (for the sparsity mask)
  impGetOmega(omModel);
  if ((int)omModel.n_rows == impNeta() && omModel.n_rows > 0) {
    impSetOmega(omModel, impDiagXform());
  }
  ctl.omModel = omModel;
  impUpdateMixProbs();          // populate mixProb from the mix() proportion thetas (npag skips updateTheta)
  impIterPrintStart();          // shared scale.h iteration printer + parHist
  npagResult r;
  try {
    r = npagRunCycle(lower, upper, ctl);
  } catch (const std::exception& ex) {
    Rcpp::stop(std::string("npag cycle failed: ") + ex.what());
  }

  int nsub = impNsub();
  int neta = impNeta();
  int nspp = (int)r.theta.n_rows;

  // per-subject posterior-mean eta: eta_i = sum_k w_ik phi_k,
  // w_ik = lambda_k psi(i,k) / sum_j lambda_j psi(i,j)  (row scaling of psi cancels)
  arma::mat postEta(nsub, neta, arma::fill::zeros);
  for (int i = 0; i < nsub; ++i) {
    arma::rowvec num = r.lambda.t() % r.psi.row(i);
    double den = arma::accu(num);
    arma::rowvec w = (den > 0.0) ? (num / den)
      : arma::rowvec(nspp, arma::fill::value(1.0 / (double)nspp));
    postEta.row(i) = w * r.theta;   // (1 x nspp)(nspp x neta)
  }
  // global-optimality certificate (Yamada Sec 2.9): D(F) = max_theta D(theta,F),
  // D(theta,F) = sum_i p(y_i|theta)/p(y_i|F) - N.  Evaluated at the FINAL gamma
  // (so p(y_i|theta) and p(y_i|F) share the fitted residual scale) over a fresh
  // Sobol scan of the box; at the NPML max D(theta,F) ~ 0 (a large positive value
  // means the grid missed a mode -- not the global optimum).
  double npagDF = NA_REAL;
  {
    arma::mat psiSup;
    npBuildPsiCoreGamma(r.theta, ctl.cores, r.gamma, psiSup);  // nsub x nspp (unnormalized)
    arma::vec pyf = psiSup * r.lambda;                 // p(y_i|F)
    pyf = arma::clamp(pyf, 1e-300, arma::datum::inf);
    int nScan = std::max(2048, 4 * ctl.points);
    arma::mat scan = npSobolGrid(nScan, lower, upper);
    arma::mat psiScan;
    npBuildPsiCoreGamma(scan, ctl.cores, r.gamma, psiScan);    // nsub x nScan
    double dmax = -arma::datum::inf;
    for (int k = 0; k < nScan; ++k) {
      double d = arma::accu(psiScan.col(k) / pyf) - (double)nsub;
      if (d > dmax) dmax = d;
    }
    // the support points themselves should give ~0; include them for the bound
    for (int k = 0; k < (int)psiSup.n_cols; ++k) {
      double d = arma::accu(psiSup.col(k) / pyf) - (double)nsub;
      if (d > dmax) dmax = d;
    }
    npagDF = dmax;
  }

  arma::mat finalSupport = r.theta;   // collapsed in-place for mu-expanded etas
  arma::mat Omega = npFinalizeFit(e, finalSupport, r.lambda, postEta, r.objf, omModel,
                                  r.gamma, residScaleIdx,
                                  ctl.muExpandEtaIdx, ctl.muExpandThetaIdx);
  impIterPrintGet(e);          // closing rule + stash e$parHistData
  e["npagDF"] = npagDF;        // global-optimality certificate
  // nonparametric outputs (support-point distribution + trace)
  e["npagSupport"] = wrap(finalSupport);   // nspp x neta (eta space)
  e["npagWeights"] = wrap(r.lambda);
  e["npagPosteriorEta"] = wrap(postEta);   // nsub x neta
  e["npagOmega"] = wrap(Omega);
  e["npagGamma"] = r.gamma;
  e["npagNspp"] = nspp;
  e["npagCycles"] = r.cycle;
  e["npagConverged"] = r.converged;
  e["npagLogLik"] = r.objf;
}

// Shared finalization (see np.h): summarize the discrete mixing distribution into
// the population theta shift + Omega, push into the FOCEi state, build the fit
// env, and set the nonparametric objective.  Returns the full support covariance.
arma::mat npFinalizeFit(Environment e, arma::mat& support,
                        const arma::vec& weights, arma::mat postEta,
                        double objf, const arma::mat& omModel,
                        double gamma, const std::vector<int>& residScaleIdx,
                        const std::vector<int>& injEtaIdx,
                        const std::vector<int>& injThetaIdx) {
  int nsub = impNsub();
  int neta = impNeta();
  arma::rowvec meanEta = weights.t() * support;               // 1 x neta
  // mu-expanded (injected) etas -> recover the parameter as a FIXED effect: fold
  // the injected eta's support-weighted mean into its paired theta, then collapse
  // that eta's random effect (zero its support column + per-subject eta, so it
  // contributes no BSV).  This is the saem "fix eta to 0 (with the mean in theta)".
  for (size_t x = 0; x < injEtaIdx.size(); ++x) {
    int j = injEtaIdx[x];
    int t = (x < injThetaIdx.size()) ? injThetaIdx[x] : -1;
    if (j < 0 || j >= neta || t < 0) continue;
    impSetThetaAll(t, impGetFullThetaVal(t) + meanEta[j]);
    support.col(j).zeros();
    postEta.col(j).zeros();
    meanEta[j] = 0.0;
  }
  arma::mat Omega(neta, neta, arma::fill::zeros);
  for (int k = 0; k < (int)support.n_rows; ++k) {
    arma::rowvec d = support.row(k) - meanEta;
    Omega += weights[k] * (d.t() * d);
  }
  // push the population mean into the mu-referenced thetas (mean-shift), then
  // re-center the per-subject etas around it so individual params are unchanged.
  for (int i = 0; i < nsub; ++i) { arma::vec ev = postEta.row(i).t(); impSetEta(i, ev); }
  impMuInterceptStep();
  impUpdateMuThetas();
  for (int i = 0; i < nsub; ++i) {
    arma::vec ev = (postEta.row(i) - meanEta).t();
    impSetEta(i, ev);
  }
  // Install the support-point covariance masked by the MODEL's Omega sparsity:
  // the diagonal always, plus off-diagonals only where the model has correlated
  // etas (a full covariance on a diagonal model would create free parameters that
  // do not exist).  The full covariance is returned for reporting.
  // npMaskedOmega restores fixed-Omega diagonals to the model value; mirror that
  // in the reported covariance too.
  arma::mat omInstall = npMaskedOmega(Omega, omModel);
  std::vector<int> fixedEta;
  impGetOmegaFixedEta(fixedEta);
  bool haveModel = ((int)omModel.n_rows == neta && (int)omModel.n_cols == neta);
  for (size_t f = 0; f < fixedEta.size(); ++f) {
    int fi = fixedEta[f];
    if (fi >= 0 && fi < neta && haveModel) Omega(fi, fi) = omModel(fi, fi);
  }
  // positive-definite floor: a support dimension can collapse to ~0 variance (a
  // near-point-mass -- common for a mu-expanded fixed-effect eta), which makes the
  // installed Omega singular and the posthoc MAP's Omega^-1 fail.  Floor the
  // diagonal so the inverse is well-defined; the tiny variance still reads as "no
  // BSV" for that parameter.
  for (int d = 0; d < neta; ++d) if (omInstall(d, d) < 1e-6) omInstall(d, d) = 1e-6;
  // mu-expanded etas carry a FIXED omega only so they are excluded from the free
  // omega objective (and, in mixtures, do not trip the multi-free-eta omega setup).
  // At finalization they are collapsed to a fixed effect: their support-mean is
  // already folded into the theta, so force their installed AND reported BSV to ~0,
  // overriding the fixed-Omega restore above.
  for (size_t x = 0; x < injEtaIdx.size(); ++x) {
    int j = injEtaIdx[x];
    if (j >= 0 && j < neta) { omInstall(j, j) = 1e-6; Omega(j, j) = 0.0; }
  }
  impSetOmega(omInstall, impDiagXform());
  // fold the fitted assay-error multiplier (gamma) into the variance-scale
  // residual coefficients: r was fit as gamma^2 * r(theta), which equals
  // r(gamma*theta) because r is homogeneous degree 2 in those coefficients.  The
  // post-cycle npResidScale is already back to 1, so scaling the thetas makes the
  // reported residual, the tables, and the objective all consistent.
  if (gamma != 1.0) {
    for (size_t s = 0; s < residScaleIdx.size(); ++s) {
      int idx = residScaleIdx[s];
      impSetThetaAll(idx, impGetFullThetaVal(idx) * gamma);
    }
  }
  impSyncInitParToFullTheta();
  impMapPass(e);
  e["objective"] = -2.0 * objf;
  return Omega;
}

//' Diagnostic: NPAG objective at a fixed grid and residual multiplier gamma
//' @param etaPoints support points, one per row
//' @param cores threads
//' @param gamma residual-error multiplier
//' @return offset-corrected marginal log-likelihood
//' @keywords internal
//' @export
// [[Rcpp::export]]
double npObjAtGamma_(arma::mat etaPoints, int cores, double gamma) {
  arma::mat psi; double offset = 0.0;
  npBuildPsiCoreScaled(etaPoints, cores, gamma, psi, &offset);
  double b = 0.0; arma::vec lam = npBurke(psi, &b);
  return b + offset;
}

//' Run the NPAG adaptive-grid cycle on a set-up inner problem
//'
//' Requires the FOCEi inner problem to be set up (\code{.npInnerSetup}).  Runs
//' the full Yamada adaptive-grid cycle (Sobol grid -> Psi -> Burke IPM ->
//' condensation -> expansion -> convergence) and returns the discrete mixing
//' distribution.  Exposed for testing ahead of the full fit-object wiring.
//'
//' @param lower,upper Numeric vectors, the per-eta support-point box.
//' @param points Initial Sobol grid size.
//' @param cycles Maximum cycles.
//' @param cores OpenMP threads.
//' @param gammaOptimize Optimize the residual-error magnitude (gamma) each cycle
//'   (only valid for uncensored normal endpoints).
//' @return A list with \code{support} (support points, eta space; one per row),
//'   \code{weights}, \code{objf} (log-likelihood), \code{gamma}, \code{cycles},
//'   and \code{converged}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npagCycle_(arma::vec lower, arma::vec upper, int points = 2028,
                      int cycles = 100, int cores = 1, bool gammaOptimize = false) {
  npagCtl ctl;
  ctl.points = points; ctl.cycles = cycles; ctl.cores = cores;
  ctl.gammaOptimize = gammaOptimize;
  npagResult r = npagRunCycle(lower, upper, ctl);
  Rcpp::NumericMatrix support(r.theta.n_rows, r.theta.n_cols);
  std::copy(r.theta.begin(), r.theta.end(), support.begin());
  Rcpp::NumericVector weights(r.lambda.begin(), r.lambda.end());
  return Rcpp::List::create(
    Rcpp::Named("support") = support,
    Rcpp::Named("weights") = weights,
    Rcpp::Named("objf") = r.objf,
    Rcpp::Named("gamma") = r.gamma,
    Rcpp::Named("cycles") = r.cycle,
    Rcpp::Named("converged") = r.converged);
}

//' Burke interior-point weight solver (nonparametric maximum likelihood)
//'
//' Solves the convex nonparametric-maximum-likelihood weight problem for a fixed
//' set of support points: given the likelihood matrix \code{psi} (subjects in
//' rows, support points in columns) it returns the maximum-likelihood mixing
//' weights and the objective (log-likelihood).  Exposed for testing the C++
//' interior-point routine against golden fixtures.
//'
//' @param psi Numeric matrix, \code{psi[i, k] = p(y_i | support point k)}, with
//'   subjects in rows and support points in columns.
//' @return A list with \code{weights} (length \code{ncol(psi)}, non-negative,
//'   summing to 1) and \code{objective} (the maximized log-likelihood).
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npIpmBurke(arma::mat psi) {
  double obj = 0.0;
  arma::vec lam = npBurke(psi, &obj);
  Rcpp::NumericVector weights(lam.begin(), lam.end());  // plain vector, not Nx1
  return Rcpp::List::create(Rcpp::Named("weights") = weights,
                            Rcpp::Named("objective") = obj);
}

//' Sobol initial grid over a box (nonparametric engines)
//'
//' @param n Number of support points.
//' @param lower,upper Numeric vectors giving the per-dimension box bounds.
//' @return Numeric matrix, one support point per row.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix npSobolGrid_(int n, arma::vec lower, arma::vec upper) {
  arma::mat g = npSobolGrid(n, lower, upper);
  Rcpp::NumericMatrix out(g.n_rows, g.n_cols);
  std::copy(g.begin(), g.end(), out.begin());
  return out;
}

//' Condense support points (nonparametric engines)
//'
//' @param lambda Support-point weights.
//' @param psi Conditional-likelihood matrix (subjects x support points).
//' @param ratio Weight-threshold ratio (keep weight > max*ratio).
//' @param tol QR rank-revealing tolerance.
//' @return List with 1-based kept indices from the weight threshold
//'   (\code{weightKeep}) and from the subsequent QR pass (\code{qrKeep}).
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npCondense_(arma::vec lambda, arma::mat psi, double ratio = 1e-3,
                       double tol = 1e-8) {
  arma::uvec wk = npCondenseWeights(lambda, ratio);
  // apply the weight keep to psi columns, then QR-condense the survivors
  arma::mat psiW = psi.cols(wk);
  arma::uvec qk = npCondenseQR(psiW, tol);
  arma::uvec qkOrig = wk.elem(qk);   // map QR keep back to original indices
  return Rcpp::List::create(
    Rcpp::Named("weightKeep") = Rcpp::IntegerVector(wk.begin(), wk.end()) + 1,
    Rcpp::Named("qrKeep") = Rcpp::IntegerVector(qkOrig.begin(), qkOrig.end()) + 1);
}
