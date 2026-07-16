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

// SAEM's Nelder-Mead (src/neldermead.cpp): func(p, fval) writes the objective at
// p into *fval.  Reused to optimize the residual-error thetas.
typedef void (*np_fn_ptr)(double *, double *);
extern "C" void nelder_fn(np_fn_ptr func, int n, double *start, double *step,
                          int itmax, double ftol_rel, double rcoef, double ecoef,
                          double ccoef, int *iconv, int *it, int *nfcall,
                          double *ynewlo, double *xmin, int *iprint);

// ---------------------------------------------------------------------------
// Residual-error theta optimization (per user: gamma is only a warm start).
// With the support points and weights FIXED, optimize the residual thetas
// (add/prop/lnorm/lambda/ar -- every endpoint) against the nonparametric -2LL
// with the same Nelder-Mead SAEM uses for its residuals.  This recovers what a
// single global gamma cannot: the add/prop ratio and per-endpoint magnitudes.
// Context is passed to the C-callback through file-static pointers (as saem does).
namespace {
  const arma::mat* gNpSupport = nullptr;   // fixed support grid (nspp x neta)
  const arma::vec* gNpWeights = nullptr;   // fixed weights (nspp)
  std::vector<int> gNpOptIdx;              // fullTheta indices being optimized
  std::vector<int> gNpOptKind;            // 0 = free, 1 = positive (SD), 2 = corr
  int gNpCores = 1;

  // map an unconstrained Nelder-Mead coordinate to the parameter's valid range
  double npResidVal(double v, int kind) {
    if (kind == 1) return std::fabs(v);                                 // SD > 0
    if (kind == 2) return std::max(-0.999, std::min(0.999, v));         // corr in (-1,1)
    return v;                                                           // lambda etc.
  }

  // Core objective: install the candidate residual thetas, rebuild Psi on the
  // fixed support, and return -2 * sum_i log(sum_k lambda_k Psi_ik).  Uses the
  // per-row log-sum-exp normalized builder (offset carries the removed row maxima)
  // so a small residual variance cannot underflow a whole subject row to zero.
  double npResidObjVal(const double *p) {
    for (size_t j = 0; j < gNpOptIdx.size(); ++j)
      impSetThetaAll(gNpOptIdx[j], npResidVal(p[j], gNpOptKind[j]));
    // frozen Psi: the support/structural solve is cached (npFreezeBuild), so this
    // only recomputes f/r for the candidate residual thetas -- no re-integration.
    arma::mat psi; double offset = 0.0;
    npFreezePsiScaled(*gNpSupport, gNpCores, psi, &offset);
    arma::vec pyf = psi * (*gNpWeights);
    pyf = arma::clamp(pyf, 1e-300, arma::datum::inf);
    return -2.0 * (arma::accu(arma::log(pyf)) + offset);
  }
  // Nelder-Mead (C, src/neldermead.cpp) callback signature.
  void npResidObjFn(double *p, double *fval) { *fval = npResidObjVal(p); }
  // newuoa (minqa via R .newuoa) callback: R-callable through Rcpp::InternalFunction.
  double npResidObjR(Rcpp::NumericVector p) { return npResidObjVal(p.begin()); }
}

// Optimize the residual thetas in idx (kinds in kind) with support/weights fixed.
// type: 0 = Nelder-Mead (C), 1 = newuoa (minqa via R .newuoa; falls back to
// Nelder-Mead if it fails).  Leaves fullTheta at the optimum, returns the -2LL.
static double npOptimizeResid(const arma::mat& support, const arma::vec& weights,
                              const std::vector<int>& idx,
                              const std::vector<int>& kind, int cores, int type) {
  int n = (int)idx.size();
  if (n == 0) return R_NegInf;
  gNpSupport = &support; gNpWeights = &weights;
  gNpOptIdx = idx; gNpOptKind = kind; gNpCores = cores;
  // solve + cache the ODE states once; every objective evaluation below reuses
  // them and only recomputes f/r for the candidate residual thetas.
  npFreezeBuild(support, cores);
  std::vector<double> start(n), step(n), xmin(n);
  for (int j = 0; j < n; ++j) {
    double v = impGetFullThetaVal(idx[j]);
    if (kind[j] == 1) v = std::fabs(v);
    start[j] = v;
    step[j] = std::max(0.1 * std::fabs(v), 0.05);
  }
  if (type == 1) {
    Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
    Rcpp::Function newuoa = nlmixr2[".newuoa"];
    Rcpp::InternalFunction fn(&npResidObjR);
    Rcpp::NumericVector par0(n);
    for (int j = 0; j < n; ++j) par0[j] = start[j];
    Rcpp::List ret = newuoa(Rcpp::_["par"] = par0, Rcpp::_["fn"] = fn,
                            Rcpp::_["control"] = Rcpp::List::create(
                              Rcpp::_["rhobeg"] = 0.2, Rcpp::_["rhoend"] = 1e-5,
                              Rcpp::_["maxfun"] = 200 * n * n));
    double f = as<double>(ret["value"]);
    if (!ISNA(f)) {
      Rcpp::NumericVector x = ret["x"];
      for (int j = 0; j < n; ++j) impSetThetaAll(idx[j], npResidVal(x[j], kind[j]));
      npFreezeClear();
      return f;
    }
    // newuoa failed -> fall through to Nelder-Mead
  }
  int iconv = 0, it = 0, nfcall = 0, iprint = -1; double ynewlo = 0.0;
  nelder_fn(npResidObjFn, n, start.data(), step.data(), 300, 1e-6, 1.0, 2.0, 0.5,
            &iconv, &it, &nfcall, &ynewlo, xmin.data(), &iprint);
  for (int j = 0; j < n; ++j) impSetThetaAll(idx[j], npResidVal(xmin[j], kind[j]));
  npFreezeClear();
  return ynewlo;
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
  int residType = 1;          // optimizer: 0 Nelder-Mead, 1 newuoa (default)
  std::vector<int> residOptIdx;  // fullTheta indices of ALL non-fixed residual thetas
  std::vector<int> residOptKind; // per-idx: 0 free (lambda), 1 SD (>0), 2 corr
  std::vector<int> residScaleIdx; // variance-scale subset (gamma warm-start fold)
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
  // A single variance-scale residual parameter (e.g. pure additive or pure
  // proportional error) is exactly what the gamma up/down search estimates -- it
  // folds gamma into that one theta -- so skip the (redundant) Nelder-Mead/newuoa
  // step there.  The residual optimizer is only needed for >= 2 residual params
  // (combined error, multiple endpoints) or a lone transform/correlation param.
  bool singleScale = (ctl.residOptIdx.size() == 1 && ctl.residScaleIdx.size() == 1);
  bool useResidOpt = !ctl.residOptIdx.empty() && !(singleScale && ctl.gammaOptimize);
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
      npBuildPsiCore(theta, ctl.cores, psi);
    }
    // a subject whose conditional density is 0 at EVERY support point makes the
    // Burke IPM (and condensation) degenerate to an empty problem -- report it
    // clearly instead of letting Armadillo throw "Mat::max(): object has no
    // elements".  The usual cause is a transform-both-sides link evaluated where
    // the structural prediction is non-positive (e.g. lnorm/log at an observation
    // whose model prediction is 0).
    if (cycle == 1) {
      arma::vec rowSum = arma::sum(psi, 1);
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
      lam = npBurke(psi, &objf);
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
    // residual theta optimization (per user: gamma is only a warm start).  With
    // the support points and weights FIXED, optimize ALL non-fixed residual thetas
    // (add/prop/lnorm/lambda/ar, every endpoint) with Nelder-Mead, then re-solve
    // the weights at the new thetas -- block-coordinate ascent.  "alternate" runs
    // this every cycle; "final" (residMode 2) runs it once after the loop.
    if (ctl.residMode == 1 && useResidOpt) {
      npOptimizeResid(theta, lam, ctl.residOptIdx, ctl.residOptKind, ctl.cores, ctl.residType);
      double off = 0.0, b = 0.0;
      npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
      lam = npBurke(psi, &b); objf = b + off;
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
  // "final" mode: optimize the residual thetas once at the converged support.
  if (ctl.residMode == 2 && useResidOpt) {
    npOptimizeResid(theta, lam, ctl.residOptIdx, ctl.residOptKind, ctl.cores, ctl.residType);
    double off = 0.0, b = 0.0;
    npBuildPsiCoreScaled(theta, ctl.cores, 1.0, psi, &off);
    lam = npBurke(psi, &b); objf = b + off;
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
  if (control.containsElementNamed("npResidType"))
    ctl.residType = as<int>(control["npResidType"]);
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

  arma::mat Omega = npFinalizeFit(e, r.theta, r.lambda, postEta, r.objf, omModel,
                                  r.gamma, residScaleIdx);
  impIterPrintGet(e);          // closing rule + stash e$parHistData
  e["npagDF"] = npagDF;        // global-optimality certificate
  // nonparametric outputs (support-point distribution + trace)
  e["npagSupport"] = wrap(r.theta);        // nspp x neta (eta space)
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
arma::mat npFinalizeFit(Environment e, const arma::mat& support,
                        const arma::vec& weights, const arma::mat& postEta,
                        double objf, const arma::mat& omModel,
                        double gamma, const std::vector<int>& residScaleIdx) {
  int nsub = impNsub();
  int neta = impNeta();
  arma::rowvec meanEta = weights.t() * support;               // 1 x neta
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
