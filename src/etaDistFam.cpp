// Maximum-likelihood M-step for a declared random-effect distribution.
//
// In the (y, eta) augmentation the complete-data likelihood factors as
//   log p(y | eta) + log p(eta | theta_dist)
// and theta_dist appears ONLY in the second term, so its M-step is a pure
// distribution fit to the etas -- no data term, no ODE solve -- exactly as the
// residual step fits accumulated residuals rather than re-solving.
// rxEtaDistExpand() breaks that by rewriting eta = Q(phi(z)) with z ~ N(0,1),
// which moves theta_dist into the DATA likelihood and lands it in a
// derivative-free search over ODE solves.  EM lets the augmentation be chosen
// freely: sample in z-space, take this M-step in eta-space.
//
// Not circular: prior draws of z would make Q(phi(z)) exactly family(theta_old)
// and return theta_old, but these are POSTERIOR draws -- which is why their
// measured spread is below 1 rather than equal to it -- so they carry data
// information.
//
// Estimated in the family's NATIVE parameters.  The objective is C++ and stays
// there: newuoa is reached through R, but the function it is handed is an
// Rcpp::InternalFunction pointing back at C++, so the per-evaluation work never
// crosses the boundary -- only the single optimizer call per M-step does.
// Positive parameters are optimized on the log scale because neither newuoa nor
// nelder_fn is bounded.
#include <RcppArmadillo.h>
#include <rxode2llPtrs.h>
#include <n1qn1c.h>
#include "etaDistFam.h"
#include "etaDistExpr.h"

// Defines the rxode2ll function pointers and the .Call entry that fills them
// from rxode2ll::.rxode2llPtr() at load (R/zzz.R).  External pointers rather
// than R_GetCCallable: the latter leaves this package compiled against a cached
// address and a typedef'd signature, so a rxode2ll update would need a rebuild
// here and a reload would leave the pointers dangling.
//
// The macro names its entry iniRxode2llPtrs; init.c registers the
// package-prefixed symbol, so rename it the way inner.cpp does for n1qn1.
extern "C" {
#define iniRxode2llPtrs _nlmixr2est_iniRxode2llPtrs
iniRxode2ll
}
#include <algorithm>

typedef void (*fn_ptr) (double *, double *);
extern "C" void nelder_fn(fn_ptr func, int n, double *start, double *step,
                          int itmax, double ftol_rel, double rcoef, double ecoef,
                          double ccoef, int *iconv, int *it, int *nfcall,
                          double *ynewlo, double *xmin, int *iprint);

// nelder_fn takes a plain function pointer, so the objective's data has to be
// reachable at file scope.  Both estimators call the M-step from serial code
// (saem between MCMC sweeps, imp between E-steps), so this is not shared across
// threads.
static std::vector<double> gEtaDistVals;   // etas being fit
static std::vector<double> gEtaDistW;      // matching weights (empty => unit)
static int gEtaDistFam = 0;
static int gEtaDistNa = 0;
static int gEtaDistPos = 0;

static inline void gEtaDistUnpack(const double *p, double *a) {
  for (int i = 0; i < gEtaDistNa; ++i) {
    a[i] = (gEtaDistPos & (1 << i)) ? std::exp(p[i]) : p[i];
  }
}

static double gEtaDistObj(const double *p) {
  double a[4];
  gEtaDistUnpack(p, a);
  for (int i = 0; i < gEtaDistNa; ++i) if (!std::isfinite(a[i])) return 1e300;
  double nll = 0.0;
  const size_t n = gEtaDistVals.size();
  const bool wtd = !gEtaDistW.empty();
  for (size_t i = 0; i < n; ++i) {
    double wi = wtd ? gEtaDistW[i] : 1.0;
    if (wi == 0.0) continue;
    double l = rxEtaDistLogD(gEtaDistFam, gEtaDistVals[i], a);
    if (!std::isfinite(l)) return 1e300;
    nll -= wi*l;
  }
  return std::isfinite(nll) ? nll : 1e300;
}
static void gEtaDistNmFn(double *p, double *fx) { *fx = gEtaDistObj(p); }

// ---- exact gradients, from rxode2ll ---------------------------------------
//
// Every family this dispatch fits has a Stan-backed log-density in rxode2ll
// with analytic derivatives, reached through the external pointers installed at
// load (see the iniRxode2ll above).  With them the family MLE is a smooth,
// low-dimensional problem WITH a gradient, so it can be handed to n1qn1 instead
// of a derivative-free search.
//
// Writes the log-density to *ll and d(logD)/d(native parameter) into g[0..na).
// Returns false when the family has no pointer available -- the caller then
// falls back to the derivative-free path, which is always correct.
//
// ret[] is the caller-allocated cache rxode2ll's <Fam>Full() uses: 3 + 2*na
// doubles, so 9 covers the widest family here (3 parameters).
static bool rxEtaDistGradD1(rxLlik1_t f, rxLlik1_t d0,
                            double x, const double *a, double *ll, double *g) {
  if (f == NULL || d0 == NULL) return false;
  double ret[9]; std::fill_n(ret, 9, 0.0);
  *ll = f(ret, x, a[0]);
  g[0] = d0(ret, x, a[0]);
  return true;
}
static bool rxEtaDistGradD2(rxLlik2_t f, rxLlik2_t d0, rxLlik2_t d1,
                            double x, const double *a, double *ll, double *g) {
  if (f == NULL || d0 == NULL || d1 == NULL) return false;
  double ret[9]; std::fill_n(ret, 9, 0.0);
  *ll = f(ret, x, a[0], a[1]);
  g[0] = d0(ret, x, a[0], a[1]);
  g[1] = d1(ret, x, a[0], a[1]);
  return true;
}
static bool rxEtaDistGradD3(rxLlik3_t f, rxLlik3_t d0, rxLlik3_t d1, rxLlik3_t d2,
                            double x, const double *a, double *ll, double *g) {
  if (f == NULL || d0 == NULL || d1 == NULL || d2 == NULL) return false;
  double ret[9]; std::fill_n(ret, 9, 0.0);
  *ll = f(ret, x, a[0], a[1], a[2]);
  g[0] = d0(ret, x, a[0], a[1], a[2]);
  g[1] = d1(ret, x, a[0], a[1], a[2]);
  g[2] = d2(ret, x, a[0], a[1], a[2]);
  return true;
}

bool rxEtaDistGradD(int fam, double x, const double *a, double *ll, double *g) {
  switch (fam) {
  case RXETADIST_NORM:
    return rxEtaDistGradD2(_p_rxLlikNorm, _p_rxLlikNormDmean, _p_rxLlikNormDsd, x, a, ll, g);
  case RXETADIST_STUDENTT:
    return rxEtaDistGradD3(_p_rxLlikT, _p_rxLlikTDdf, _p_rxLlikTDmean, _p_rxLlikTDsd, x, a, ll, g);
  case RXETADIST_CAUCHY:
    return rxEtaDistGradD2(_p_rxLlikCauchy, _p_rxLlikCauchyDlocation, _p_rxLlikCauchyDscale, x, a, ll, g);
  case RXETADIST_DBLEXP:
    return rxEtaDistGradD2(_p_rxLlikDblExp, _p_rxLlikDblExpDMu, _p_rxLlikDblExpDSigma, x, a, ll, g);
  case RXETADIST_LOGIS:
    return rxEtaDistGradD2(_p_rxLlikLogis, _p_rxLlikLogisDLocation, _p_rxLlikLogisDScale, x, a, ll, g);
  case RXETADIST_GUMBEL:
    return rxEtaDistGradD2(_p_rxLlikGumbel, _p_rxLlikGumbelDMu, _p_rxLlikGumbelDBeta, x, a, ll, g);
  case RXETADIST_LNORM:
    return rxEtaDistGradD2(_p_rxLlikLnorm, _p_rxLlikLnormDMeanlog, _p_rxLlikLnormDSdlog, x, a, ll, g);
  case RXETADIST_CHISQ:
    return rxEtaDistGradD1(_p_rxLlikChisq, _p_rxLlikChisqDdf, x, a, ll, g);
  case RXETADIST_INVCHISQ:
    return rxEtaDistGradD1(_p_rxLlikInvChisq, _p_rxLlikInvChisqDNu, x, a, ll, g);
  case RXETADIST_SCINVCHISQ:
    return rxEtaDistGradD2(_p_rxLlikScaledInvChisq, _p_rxLlikScaledInvChisqDNu, _p_rxLlikScaledInvChisqDSigma, x, a, ll, g);
  case RXETADIST_EXP:
    return rxEtaDistGradD1(_p_rxLlikExp, _p_rxLlikExpDrate, x, a, ll, g);
  case RXETADIST_GAMMA:
    return rxEtaDistGradD2(_p_rxLlikGamma, _p_rxLlikGammaDshape, _p_rxLlikGammaDrate, x, a, ll, g);
  case RXETADIST_INVGAMMA:
    return rxEtaDistGradD2(_p_rxLlikInvGamma, _p_rxLlikInvGammaDAlpha, _p_rxLlikInvGammaDBeta, x, a, ll, g);
  case RXETADIST_WEIBULL:
    return rxEtaDistGradD2(_p_rxLlikWeibull, _p_rxLlikWeibullDshape, _p_rxLlikWeibullDscale, x, a, ll, g);
  case RXETADIST_FRECHET:
    return rxEtaDistGradD2(_p_rxLlikFrechet, _p_rxLlikFrechetDAlpha, _p_rxLlikFrechetDSigma, x, a, ll, g);
  case RXETADIST_RAYLEIGH:
    return rxEtaDistGradD1(_p_rxLlikRayleigh, _p_rxLlikRayleighDSigma, x, a, ll, g);
  case RXETADIST_PARETO:
    return rxEtaDistGradD2(_p_rxLlikPareto, _p_rxLlikParetoDYMin, _p_rxLlikParetoDAlpha, x, a, ll, g);
  case RXETADIST_PARETO2:
    return rxEtaDistGradD3(_p_rxLlikParetoType2, _p_rxLlikParetoType2DMu, _p_rxLlikParetoType2DLambda, _p_rxLlikParetoType2DAlpha, x, a, ll, g);
  case RXETADIST_BETA:
    return rxEtaDistGradD2(_p_rxLlikBeta, _p_rxLlikBetaDshape1, _p_rxLlikBetaDshape2, x, a, ll, g);
  case RXETADIST_BETAPROP:
    return rxEtaDistGradD2(_p_rxLlikBetaProportion, _p_rxLlikBetaProportionDMu, _p_rxLlikBetaProportionDKappa, x, a, ll, g);
  case RXETADIST_UNIF:
    return rxEtaDistGradD2(_p_rxLlikUnif, _p_rxLlikUnifDalpha, _p_rxLlikUnifDbeta, x, a, ll, g);
  default: return false;
  }
}

// Objective AND gradient of the family MLE at the OPTIMIZER's coordinates.
// Positive parameters are carried on the log scale, so the chain rule is
// d/d(log p) = p * d/dp.  Returns false if any observation has no usable
// gradient, which sends the caller back to the derivative-free path.
static bool gEtaDistObjGrad(const double *p, double *fx, double *gr) {
  double a[4];
  gEtaDistUnpack(p, a);
  for (int i = 0; i < gEtaDistNa; ++i) if (!std::isfinite(a[i])) return false;
  double nll = 0.0;
  std::fill_n(gr, gEtaDistNa, 0.0);
  const size_t n = gEtaDistVals.size();
  const bool wtd = !gEtaDistW.empty();
  double g1[4], ll;
  for (size_t i = 0; i < n; ++i) {
    double wi = wtd ? gEtaDistW[i] : 1.0;
    if (wi == 0.0) continue;
    if (!rxEtaDistGradD(gEtaDistFam, gEtaDistVals[i], a, &ll, g1)) return false;
    if (!std::isfinite(ll)) return false;
    nll -= wi*ll;
    for (int k = 0; k < gEtaDistNa; ++k) {
      if (!std::isfinite(g1[k])) return false;
      gr[k] -= wi*g1[k];
    }
  }
  if (!std::isfinite(nll)) return false;
  // chain rule onto the optimizer's scale
  for (int k = 0; k < gEtaDistNa; ++k) {
    if (gEtaDistPos & (1 << k)) gr[k] *= a[k];
    if (!std::isfinite(gr[k])) return false;
  }
  *fx = nll;
  return true;
}

static int gEtaDistN1Bad = 0;
static void gEtaDistN1Cost(int *ind, int *nn, double *x, double *f, double *g,
                           int *ti, float *tr, double *td, int *id) {
  (void)ti; (void)tr; (void)td; (void)id;
  if (gEtaDistN1Bad) return;
  double fv = 0.0, gg[4];
  if (!gEtaDistObjGrad(x, &fv, gg)) { gEtaDistN1Bad = 1; return; }
  if (*ind == 2 || *ind == 4) *f = fv;
  if (*ind == 3 || *ind == 4) for (int i = 0; i < *nn; ++i) g[i] = gg[i];
}


// newuoa reaches the same objective through R (nlmixr2est:::.newuoa), the way
// SAEM's own _saemType==2 residual step does.  The objective itself stays in
// C++ -- Rcpp::InternalFunction hands newuoa a pointer, so the per-evaluation
// work never crosses into R; only the one optimizer call per M-step does.
double gEtaDistNewuoaFn(Rcpp::NumericVector p) {
  std::vector<double> pv(p.begin(), p.end());
  double v = gEtaDistObj(pv.data());
  if (!std::isfinite(v)) return 1e300;
  return v;
}

// ---- native arguments -> the user's thetas, in C++ -------------------------
//
// The M-step fits a family's NATIVE parameters; the declaration writes those as
// expressions over ini() thetas, so the fitted values have to be mapped back.
// That map was the last part of this M-step that had to call R: an R
// Nelder-Mead over an objective that eval()'d the argument expressions, invoked
// from inside the C++ loop.
//
// Same objective as the R version it replaces: relative, on the log scale, so
// arguments on very different scales weigh comparably, with a sign penalty.
static std::vector<std::vector<etaDistTok> > gEtaDistRpn;
static std::vector<double> gEtaDistTarget;
static int gEtaDistNth = 0;

static double gEtaDistMapObj(const double *p) {
  double v = 0.0;
  for (size_t k = 0; k < gEtaDistRpn.size(); ++k) {
    double a = etaDistExprEval(gEtaDistRpn[k], p, gEtaDistNth);
    if (!std::isfinite(a)) return 1e10;
    double t = gEtaDistTarget[k];
    double la = std::log(std::max(std::fabs(a), 1e-300));
    double lt = std::log(std::max(std::fabs(t), 1e-300));
    v += (la - lt)*(la - lt);
    if ((a < 0) != (t < 0)) v += 1e3;
  }
  return std::isfinite(v) ? v : 1e10;
}
static void gEtaDistMapNmFn(double *p, double *fx) { *fx = gEtaDistMapObj(p); }

// Returns false when any expression is outside the C++ grammar, or the solve
// does not converge -- the caller then keeps the R route, which handles the
// general case.
bool rxEtaDistArgsToThetas(const std::vector<std::string> &exprs,
                           const std::vector<std::string> &thetaNames,
                           const double *start, const double *target,
                           double *out) {
  int nth = (int)thetaNames.size();
  if (nth <= 0 || nth > 32 || exprs.size() != (size_t)0 + exprs.size()) return false;
  gEtaDistRpn.assign(exprs.size(), std::vector<etaDistTok>());
  for (size_t k = 0; k < exprs.size(); ++k) {
    if (!etaDistExprParse(exprs[k], thetaNames, gEtaDistRpn[k])) return false;
  }
  gEtaDistTarget.assign(target, target + exprs.size());
  gEtaDistNth = nth;
  std::vector<double> st(start, start + nth), stp((size_t)nth), xm(start, start + nth);
  for (int i = 0; i < nth; ++i)
    stp[(size_t)i] = (std::fabs(st[(size_t)i]) > 1e-8) ? 0.1*std::fabs(st[(size_t)i]) : 0.1;
  int iconv, it, nfcall, iprint = 0;
  double ynewlo = R_PosInf;
  // SAEM's own simplex (src/neldermead.cpp), not an R one
  nelder_fn(gEtaDistMapNmFn, nth, st.data(), stp.data(), 200*nth, 1e-10,
            1.0, 2.0, 0.5, &iconv, &it, &nfcall, &ynewlo, xm.data(), &iprint);
  if (!std::isfinite(ynewlo) || ynewlo > 1e-6) return false;
  for (int i = 0; i < nth; ++i) {
    if (!std::isfinite(xm[(size_t)i])) return false;
    out[i] = xm[(size_t)i];
  }
  return true;
}

//[[Rcpp::export]]
Rcpp::NumericVector rxEtaDistArgsToThetasTest_(Rcpp::CharacterVector exprs,
                                               Rcpp::CharacterVector thetaNames,
                                               Rcpp::NumericVector start,
                                               Rcpp::NumericVector target) {
  std::vector<std::string> e, tn;
  for (int i = 0; i < exprs.size(); ++i) e.push_back(Rcpp::as<std::string>(exprs[i]));
  for (int i = 0; i < thetaNames.size(); ++i) tn.push_back(Rcpp::as<std::string>(thetaNames[i]));
  std::vector<double> out((size_t)tn.size(), 0.0);
  if (!rxEtaDistArgsToThetas(e, tn, start.begin(), target.begin(), out.data())) {
    return Rcpp::NumericVector(0);
  }
  Rcpp::NumericVector r(tn.size());
  for (size_t i = 0; i < out.size(); ++i) r[i] = out[i];
  r.names() = thetaNames;
  return r;
}

bool rxEtaDistMleW(int fam, const std::vector<double> &vals,
                   const std::vector<double> *w, double *a0) {
  int na = rxEtaDistNarg(fam);
  if (na <= 0 || vals.size() < 2) return false;
  gEtaDistFam = fam; gEtaDistNa = na; gEtaDistPos = rxEtaDistPosMask(fam);
  gEtaDistVals = vals;
  gEtaDistW.clear();
  if (w != nullptr) {
    if (w->size() != vals.size()) return false;
    // A weight vector that sums to nothing carries no information; refusing is
    // the honest answer, and leaves the caller's parameters where they were.
    double sw = 0.0;
    for (size_t i = 0; i < w->size(); ++i) {
      double wi = (*w)[i];
      if (!std::isfinite(wi) || wi < 0.0) return false;
      sw += wi;
    }
    if (!(sw > 0.0)) return false;
    gEtaDistW = *w;
  }
  std::vector<double> st(na), stp(na), xm(na);
  for (int i = 0; i < na; ++i) {
    double v = (gEtaDistPos & (1 << i)) ? std::log(a0[i]) : a0[i];
    if (!std::isfinite(v)) return false;
    st[i] = v; xm[i] = v;
    // nelder_fn derives nothing from the start, so give every coordinate a
    // usable step even when it starts at zero
    stp[i] = (std::fabs(v) > 1e-8) ? 0.1*std::fabs(v) : 0.1;
  }
  // newuoa first.  It builds a quadratic model from its interpolation points,
  // so on a smooth low-dimensional MLE like this it needs far fewer objective
  // evaluations than a simplex walking downhill -- and this objective is
  // evaluated over every sampled eta, every M-step call, for every declared
  // distribution.  Nelder-Mead is kept as the fallback for when newuoa returns
  // nothing usable, mirroring what SAEM's own _saemType==2 residual step does.
  double ynewlo = R_PosInf;
  bool haveMin = false;
  // n1qn1 first, when rxode2ll can supply exact derivatives for this family.
  // A quasi-Newton with the true gradient converges in far fewer objective
  // evaluations than a derivative-free method, and every evaluation here walks
  // all the sampled etas.  Falls through to newuoa when the family has no
  // pointer, when n1qn1 is unavailable, or when any observation yields a
  // non-finite gradient.
  {
    double llProbe, gProbe[4];
    double aProbe[4];
    gEtaDistUnpack(st.data(), aProbe);
    bool haveGrad = (n1qn1_ != NULL) && !gEtaDistVals.empty() &&
      rxEtaDistGradD(fam, gEtaDistVals[0], aProbe, &llProbe, gProbe);
    if (haveGrad) {
      std::vector<double> x(st.begin(), st.end()), gg((size_t)na, 0.0);
      std::vector<double> zm((size_t)(na*(na+13)/2 + 1), 0.0);
      std::vector<double> var((size_t)na, 0.1);
      double f = 0.0, eps = 1e-8;
      int nn = na, mode = 1, niter = 100*na, nsim = 100*na, impr = 0, izs = 0;
      float rzs = 0; double dzs = 0; int idz = 0;
      gEtaDistN1Bad = 0;
      n1qn1_(gEtaDistN1Cost, &nn, x.data(), &f, gg.data(), var.data(), &eps,
             &mode, &niter, &nsim, &impr, zm.data(), &izs, &rzs, &dzs, &idz);
      if (getenv("NLMIXR2_ETADIST_OPT") != NULL)
        Rprintf("etaDistMle fam=%d n1qn1 bad=%d f=%.6g\n", fam, gEtaDistN1Bad, f);
      if (!gEtaDistN1Bad && std::isfinite(f) && f < 1e300) {
        bool ok = true;
        for (int i = 0; i < na; ++i) if (!std::isfinite(x[(size_t)i])) ok = false;
        if (ok) {
          for (int i = 0; i < na; ++i) xm[(size_t)i] = x[(size_t)i];
          ynewlo = f;
          haveMin = true;
        }
      }
    }
  }
  if (!haveMin && getenv("NLMIXR2_ETADIST_OPT") != NULL)
    Rprintf("etaDistMle fam=%d falling back to newuoa\n", fam);
  if (!haveMin) {
    int npt = 2*na + 1;
    Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
    Rcpp::Function newuoa = nlmixr2[".newuoa"];
    Rcpp::InternalFunction fnRef(&gEtaDistNewuoaFn);
    Rcpp::NumericVector par0(na);
    double rhobeg = 0.0;
    for (int i = 0; i < na; ++i) {
      par0[i] = st[(size_t)i];
      double sp = std::fabs(stp[(size_t)i]);
      if (sp > rhobeg) rhobeg = sp;
    }
    if (!(rhobeg > 0.0)) rhobeg = 0.1;
    Rcpp::List ret;
    bool ok = true;
    try {
      ret = newuoa(Rcpp::_["par"] = par0, Rcpp::_["fn"] = fnRef,
                   Rcpp::_["control"] = Rcpp::List::create(
                     Rcpp::_["rhobeg"] = rhobeg,
                     Rcpp::_["rhoend"] = 1e-8,
                     Rcpp::_["npt"] = npt,
                     Rcpp::_["maxfun"] = 200*na));
    } catch (...) {
      ok = false;
    }
    if (ok && ret.containsElementNamed("value") && ret.containsElementNamed("par")) {
      double f = Rcpp::as<double>(ret["value"]);
      Rcpp::NumericVector xx = ret["par"];
      if (std::isfinite(f) && f < 1e300 && (int)xx.size() == na) {
        for (int i = 0; i < na; ++i) xm[(size_t)i] = xx[i];
        ynewlo = f;
        haveMin = true;
      }
    }
  }

  if (!haveMin) {
    int iconv, it, nfcall, iprint = 0;
    nelder_fn(gEtaDistNmFn, na, st.data(), stp.data(), 200*na, 1e-8,
              1.0, 2.0, 0.5, &iconv, &it, &nfcall, &ynewlo, xm.data(), &iprint);
  }
  if (!std::isfinite(ynewlo) || ynewlo >= 1e300) return false;
  double a[4];
  gEtaDistUnpack(xm.data(), a);
  for (int i = 0; i < na; ++i) {
    if (!std::isfinite(a[i])) return false;
    a0[i] = a[i];
  }
  return true;
}

bool rxEtaDistMle(int fam, const std::vector<double> &vals, double *a0) {
  return rxEtaDistMleW(fam, vals, nullptr, a0);
}

bool rxEtaDistSpreadOk(const std::vector<double> &w, double lo, double hi,
                       double *sdOut, const std::vector<double> *wt) {
  if (sdOut != nullptr) *sdOut = NA_REAL;
  if (wt != nullptr && wt->size() < w.size()) return false;
  size_t n = 0;
  double sw = 0.0, m = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    m += wi*w[i]; sw += wi; n++;
  }
  if (n < 2 || !(sw > 0.0)) return false;
  m /= sw;
  double v = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    double d = w[i] - m; v += wi*d*d;
  }
  // reliability weights: the unbiased normalizer is sum(w) - sum(w^2)/sum(w),
  // which reduces to (n - 1) when every weight is 1
  double sw2 = 0.0;
  for (size_t i = 0; i < w.size(); ++i) {
    if (!std::isfinite(w[i])) continue;
    double wi = (wt == nullptr) ? 1.0 : (*wt)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    sw2 += wi*wi;
  }
  double den = sw - sw2/sw;
  if (!(den > 0.0)) den = sw;
  v /= den;
  double sd = std::sqrt(v > 0.0 ? v : 0.0);
  if (sdOut != nullptr) *sdOut = sd;
  return std::isfinite(sd) && sd >= lo && sd <= hi;
}

// Closed-form M-step for a Gaussian copula's correlation.
//
// The SAMPLE CORRELATION, not the raw product-moment mean(z1*z2).  The latent
// pair is bivariate normal with unit variances, so the constrained MLE is the
// product-moment -- but ONLY when the draws actually have unit variance, and
// they do not: an unmixed chain gives a pooled spread of 2.15, whence
// mean(z^2) ~ 4.6 and a true correlation of 0.3 comes out as 0.3*4.6 = 1.38,
// which clamps to 0.999.  It then STAYS there, because 0.999*mean(z^2) is
// itself ~0.999 -- and a pinned correlation makes the copula partner's latent
// numerically equal to its partner's (w_k = rho*z_j + sqrt(1-rho^2)*z_k with
// rho ~ 1), so BOTH declared families end up fitted to the same draws.
//
// Measured on Bauer's gamma model: the correlation was pinned at 0.999 by the
// first M-step of every diverging run and never moved again.
//
// Normalizing by the observed spreads costs nothing, is bounded in [-1, 1] by
// construction, and agrees with the product-moment exactly when the draws do
// have unit variance -- so it is strictly the safer estimator here.
double rxEtaDistCorMleW(const std::vector<double> &z1,
                        const std::vector<double> &z2,
                        const std::vector<double> *w) {
  size_t n = std::min(z1.size(), z2.size());
  if (n < 2) return NA_REAL;
  if (w != nullptr && w->size() < n) return NA_REAL;
  double sw = 0.0, m1 = 0.0, m2 = 0.0; size_t m = 0;
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z1[i]) || !std::isfinite(z2[i])) continue;
    double wi = (w == nullptr) ? 1.0 : (*w)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    m1 += wi*z1[i]; m2 += wi*z2[i]; sw += wi; m++;
  }
  if (m < 2 || !(sw > 0.0)) return NA_REAL;
  m1 /= sw; m2 /= sw;
  double s11 = 0.0, s22 = 0.0, s12 = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(z1[i]) || !std::isfinite(z2[i])) continue;
    double wi = (w == nullptr) ? 1.0 : (*w)[i];
    if (!std::isfinite(wi) || wi <= 0.0) continue;
    double d1 = z1[i] - m1, d2 = z2[i] - m2;
    s11 += wi*d1*d1; s22 += wi*d2*d2; s12 += wi*d1*d2;
  }
  if (!(s11 > 0.0) || !(s22 > 0.0)) return NA_REAL;
  double r = s12 / std::sqrt(s11*s22);
  if (!std::isfinite(r)) return NA_REAL;
  return std::max(std::min(r, 0.999), -0.999);
}

double rxEtaDistCorMle(const std::vector<double> &z1,
                       const std::vector<double> &z2) {
  return rxEtaDistCorMleW(z1, z2, nullptr);
}
