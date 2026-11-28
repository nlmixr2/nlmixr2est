// Importance-sampling EM kernel for est="impmap".
//
// Peer of saem.cpp / inner.cpp holding the NONMEM-style Monte Carlo importance-
// sampling EM.  It reuses the FOCEI inner machinery (mode + Hessian + joint-
// density evaluation) through the thin numeric interface in imp.h, so all
// model-evaluation heavy lifting stays in inner.cpp.
//
// Module M1: single MAP pass collecting each subject's mode/Hessian/likelihood.
// Module M2: per-subject multivariate-normal proposal + thread-safe threefry
// sampler.
// Module M3: importance weights, the individual objective contribution, and the
// E-step conditional mean and variance.
// Module M4: the EM iteration -- {re-MAP, E-step, M-step} for nIter steps, where
// the M-step updates the mu-referenced thetas and Omega from the weighted
// conditional moments, then finalizes the fit at the converged estimates.
#include <RcppArmadillo.h>
#include <rxode2ptr.h>
#include <boost/random/sobol.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <ctime>
#include <atomic>
#include "nmMcmcRng.h"
#include "nmParallel.h"
#include "impQrng.h"
#include "imp.h"
#include "logSumExp.h"
#include "odeSwap.h" // odeSwapAnyNdiffSet()
#include "utilc.h"   // RSprintf (covariance-step progress header)
#include "nmProgress.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// Largest per-component structural-theta Newton step the M-step will take in one
// EM iteration, relative to that theta's own magnitude (floored at 1 so a theta
// sitting near zero still gets a usable absolute allowance).
#define IMP_MSTEP_TRUST 0.5

// Loose sanity bound on an UNDAMPED step, relative to each theta's own
// magnitude.  A well-conditioned Hessian is allowed a large move -- that is
// ordinary optimizer progress -- so this sits far above anything a healthy
// M-step proposes and only catches a step that is plainly not a step.
#define IMP_MSTEP_SANE 20.0

// Counts of M-step Newton steps that needed Levenberg-Marquardt damping, and of
// iterations whose structural thetas could not be updated at all.  Reported once
// at the end of the fit: a run that damps constantly has an ill-conditioned
// theta-sensitivity Hessian and its structural estimates deserve a second look.
//
// NOT inside the _OPENMP guard: the M-step reads them unconditionally, and a
// toolchain without OpenMP (macOS clang) then fails to compile with "use of
// undeclared identifier 'nMStepDamped'".  Only <omp.h> belongs in the guard.
static std::atomic<int> nMStepDamped(0);
static std::atomic<int> nMStepSkipped(0);

using namespace Rcpp;

// ---- quasi-random (QRPEM) sampling elements --------------------------------
// Base Sobol point set on (0,1)^neta: isample rows, one point per row.  The
// boost engine emits one coordinate per call, cycling through the dimensions
// of consecutive points, and skips the trivial zero point; the top 53 bits are
// scaled into (0,1) with a half-step offset so no coordinate is ever 0 or 1.
// Raw engine output, row-major (isample x neta).  Held as the 64-bit integers
// so a scramble can act on the digits; the unscrambled path converts with the
// exact expression it always used.
static std::vector<uint64_t> impSobolBits(int isample, int neta) {
  boost::random::sobol eng((std::size_t)neta);
  std::vector<uint64_t> V((size_t)isample * (size_t)neta);
  for (int k = 0; k < isample; ++k) {
    for (int j = 0; j < neta; ++j) {
      V[(size_t)k * (size_t)neta + (size_t)j] = (uint64_t)eng();
    }
  }
  return V;
}

static arma::mat impSobolU0(int isample, int neta) {
  std::vector<uint64_t> V = impSobolBits(isample, neta);
  arma::mat U0(isample, neta);
  for (int k = 0; k < isample; ++k) {
    for (int j = 0; j < neta; ++j) {
      uint64_t v = V[(size_t)k * (size_t)neta + (size_t)j];
      U0(k, j) = std::ldexp((double)(v >> 11) + 0.5, -53);
    }
  }
  return U0;
}

// Scrambled Sobol points -> N(0,1).  The scramble seed is derived
// arithmetically from (base, id, iter, dimension), so this draws nothing from
// the threefry engine and replaces the Cranley-Patterson shift rather than
// composing with it.
static arma::mat impQrZScrambled(const std::vector<uint64_t>& V, int nrow, int ncol,
                                 int type, uint32_t base, int id, int iter) {
  arma::mat Z(nrow, ncol);
  for (int j = 0; j < ncol; ++j) {
    uint32_t sj = impQrngSeed(base, id, iter, j);
    for (int k = 0; k < nrow; ++k) {
      double u = impQrScrambleU(V[(size_t)k * (size_t)ncol + (size_t)j], type, sj);
      if (u < 1e-12) u = 1e-12;
      else if (u > 1.0 - 1e-12) u = 1.0 - 1e-12;
      Z(k, j) = R::qnorm(u, 0.0, 1.0, 1, 0);
    }
  }
  return Z;
}

// Uniform Sobol points -> N(0,1): optional Cranley-Patterson shift (mod 1),
// clamp away from 0/1, then the inverse normal CDF per coordinate.
static arma::mat impQrZ(const arma::mat& U0, const arma::vec* shift) {
  arma::mat Z(U0.n_rows, U0.n_cols);
  for (arma::uword j = 0; j < U0.n_cols; ++j) {
    double s = (shift != nullptr) ? (*shift)[j] : 0.0;
    for (arma::uword k = 0; k < U0.n_rows; ++k) {
      double u = U0(k, j) + s;
      u -= std::floor(u);
      if (u < 1e-12) u = 1e-12;
      else if (u > 1.0 - 1e-12) u = 1.0 - 1e-12;
      Z(k, j) = R::qnorm(u, 0.0, 1.0, 1, 0);
    }
  }
  return Z;
}

// ---- multivariate-t proposal (NONMEM DF) -----------------------------------
// A Gaussian proposal scaled by gamma can only be made WIDER, never
// heavier-tailed, and importance weights blow up (infinite variance) when the
// proposal's tails are lighter than the target's -- a failure Pareto k-hat
// detects and xi/Kish ESS cannot.  Switching the proposal to a multivariate t
// with `df` degrees of freedom fixes the SHAPE: t tails are polynomial, so for
// small df they dominate any Gaussian target tail.  df = 0 means "Gaussian",
// the historical path.
//
// Only rxNormEng/rxUnifEng are available on the thread-safe engine, so the
// chi-square is drawn by INVERSE CDF from a single uniform, using boost's
// quantile function.  Inverse CDF rather than a rejection sampler
// (Marsaglia-Tsang) on purpose: a rejection loop consumes a VARIABLE number of
// draws, so the stream position would depend on the accept/reject history and
// the fit would stop being reproducible under any change that shifts it.  One
// uniform in, one chi-square out, and it works for any real df > 0.  It also
// makes the quasi-random path natural -- a Sobol coordinate can feed it
// directly.
static inline double impChisqQuantile(double u, double df) {
  if (u <= 0.0) u = 1e-12;
  if (u >= 1.0) u = 1.0 - 1e-12;
  return boost::math::quantile(boost::math::chi_squared_distribution<double>(df), u);
}

// Scale factor turning a standard normal vector into a multivariate-t one:
// t = z * sqrt(df / chisq_df).  df <= 0 -> 1 (Gaussian: no scaling and NO
// uniform drawn, so the Gaussian path consumes exactly the stream it always
// did and stays bit-identical).
static inline double impTScale(double df) {
  if (df <= 0.0) return 1.0;
  double w = impChisqQuantile(rxUnifEng(0.0, 1.0), df);
  if (!(w > 0.0) || !R_finite(w)) return 1.0;
  return std::sqrt(df / w);
}

// Log of the peak-normalized proposal kernel's RECIPROCAL, i.e. the term that
// enters the log importance weight q_k.  Gaussian: d'Hd/(2 gamma).
// t: ((df+p)/2) * log1p(d'Hd/(gamma*df)).  Both are 0 at d = 0.
static inline double impLogKernelRecip(double quad, double gammaId, double df, int p) {
  // NB the Gaussian branch reproduces the original operation order exactly --
  // (1/(2*gamma)) * quad, not quad/(2*gamma).  They differ in the last bits,
  // and est="imp" (whose quadratic form has a different magnitude) shows the
  // reassociation as ~1e-10 drift, which would break the bit-identity contract
  // this branch has held to at every step.
  if (df <= 0.0) return (1.0 / (2.0 * gammaId)) * quad;
  return 0.5 * (df + (double)p) * std::log1p(quad / (gammaId * df));
}

// Additive correction to the Gaussian normalizer Ci when a t proposal is used:
//   (p/2) log(df/2) + lgamma(df/2) - lgamma((df+p)/2)
// -> 0 as df -> Inf (Gamma(x+a)/Gamma(x) ~ x^a), recovering the Gaussian.
static inline double impCiTCorrection(double df, int p) {
  if (df <= 0.0) return 0.0;
  return 0.5 * (double)p * std::log(0.5 * df) +
    std::lgamma(0.5 * df) - std::lgamma(0.5 * (df + (double)p));
}

// ---- generalized proposal families -----------------------------------------
// The kernel already factors the proposal as (peak-normalized kernel) x (peak
// density), with the peak's gamma and |H| dependence folded into Ci and only a
// gamma-free, family-specific peak RATIO left over -- that ratio is exactly
// what impCiTCorrection is for the t.  Everything below generalizes those two
// pieces plus the radial scale, so a family is three numbers rather than a new
// code path.
//
// Two invariants every family must hold, because the rest of the kernel is
// built on them:
//   * the kernel is EXACTLY 0 at d = 0.  qCenter (below) omits the kernel term
//     for that reason, and xi is measured against qCenter.
//   * the density depends on d only through quad = d' H d (elliptical).  A
//     skewed or shifted-component family would break qCenter, xi, iaccept and
//     Ci simultaneously.  Do not add one without revisiting all four.
//
// MVN and MVT are the historical df<=0/df>0 pair and delegate to the untouched
// functions above, so they stay bit-identical.
enum impPropType { impPropMvn = 0, impPropMvt = 1, impPropMvl = 2, impPropMix = 3 };

struct impProp {
  int type = impPropMvn;
  double df = 0.0;          // mvt
  double kappa = 1.0;       // mvl: S = kappa * Sigma
  double rScale = 1.0;      // mvl: sqrt(kappa), the Gamma scale
  int nmix = 0;             // mix: 2 or 3 components
  double c[3] = {1.0, 1.0, 1.0};   // mix: variance multipliers
  double lw[3] = {0.0, 0.0, 0.0};  // mix: log w_k - (p/2) log c_k
  double cw[3] = {0.0, 0.0, 0.0};  // mix: cumulative weights (component pick)
  double logPeak = 0.0;     // mix: logSumExp_k(lw_k)
  double corr = 0.0;        // Ci correction: log(gaussPeak) - log(familyPeak)
};

// Gamma(shape, 1) by inverse CDF, the same one-uniform-in convention as
// impChisqQuantile (a rejection sampler would consume a variable number of
// draws and break reproducibility).
static inline double impGammaQuantile(double u, double shape) {
  if (u <= 0.0) u = 1e-12;
  if (u >= 1.0) u = 1.0 - 1e-12;
  return boost::math::quantile(boost::math::gamma_distribution<double>(shape, 1.0), u);
}

// unweighted log-sum-exp over the proposal's log-weights; the max-shifted sum
// itself is rxLogSumExp() (logSumExp.h), shared with the mixture marginals
static inline double impLogSumExp3(const double* v, int n) {
  return rxLogSumExp(v, n);
}

// Build a subject's proposal from the fit-constant spec plus its own df.
// kappa = 1/(p+1) makes the MVL COVARIANCE-MATCHED: the spherical Laplace with
// S = Sigma has Cov = (p+1) Sigma, so without this gamma, iscaleMin and
// iscaleMax would silently mean something different per neta, and est="imp"'s
// Sigma <- gamma * condVar recursion would inflate by (p+1) each iteration.
static impProp impPropBuild(int type, double df, int nmix, const double* cIn,
                            const double* wIn, int p) {
  impProp pr;
  pr.type = type;
  double dp = (double)p;
  if (type == impPropMvt) {
    pr.df = df;
    pr.corr = impCiTCorrection(df, p);
  } else if (type == impPropMvl) {
    pr.kappa = 1.0 / (dp + 1.0);
    pr.rScale = std::sqrt(pr.kappa);
    // log(gaussPeak) - log(mvlPeak) for f(x) ~ exp(-sqrt(x' S^-1 x)), whose
    // normalizer is Gamma(p/2) / (2 pi^(p/2) Gamma(p) |S|^(1/2)).
    pr.corr = std::log(2.0) - 0.5 * dp * std::log(2.0) +
      std::lgamma(dp) - std::lgamma(0.5 * dp) + 0.5 * dp * std::log(pr.kappa);
  } else if (type == impPropMix) {
    pr.nmix = nmix;
    double acc = 0.0;
    for (int k = 0; k < nmix; ++k) {
      pr.c[k] = cIn[k];
      pr.lw[k] = std::log(wIn[k]) - 0.5 * dp * std::log(cIn[k]);
      acc += wIn[k];
      pr.cw[k] = acc;
    }
    pr.cw[nmix - 1] = 1.0;   // guard against fp drift in the last bucket
    pr.logPeak = impLogSumExp3(pr.lw, nmix);
    pr.corr = -pr.logPeak;
  }
  return pr;
}

// Fit-constant proposal spec, read once from the control.
struct impPropSpec { int type = impPropMvn; int nmix = 0; double c[3]; double w[3]; };

static impPropSpec impPropSpecFromControl(double dfCtl) {
  impPropSpec sp;
  sp.c[0] = sp.c[1] = sp.c[2] = 1.0;
  sp.w[0] = sp.w[1] = sp.w[2] = 0.0;
  int t = impProposalType();          // 0 auto, 1 normal, 2 t, 3 laplace, 4 mixture
  if (t == 0) sp.type = (dfCtl > 0.0) ? impPropMvt : impPropMvn;
  else if (t == 1) sp.type = impPropMvn;
  else if (t == 2) sp.type = impPropMvt;
  else if (t == 3) sp.type = impPropMvl;
  else {
    sp.type = impPropMix;
    std::vector<double> cv, wv;
    impPropMixGet(cv, wv);
    sp.nmix = (int)cv.size();
    if (sp.nmix < 2 || sp.nmix > 3 || (int)wv.size() != sp.nmix) {
      sp.type = impPropMvn; sp.nmix = 0; return sp;
    }
    double sw = 0.0;
    for (int k = 0; k < sp.nmix; ++k) sw += wv[k];
    // The scales are used AS GIVEN.  An earlier version rescaled them so
    // sum_k w_k c_k == 1 ("so gamma means for the mixture what it means for
    // MVN"), which quietly defeated the method: with the documented
    // c = (1, 9), w = (0.9, 0.1) it made the DOMINANT component 0.56x the
    // Laplace covariance, so 90% of draws came from a proposal far too narrow
    // and the weight tail got WORSE than a plain normal's (measured on the
    // 3-eta fixture: max k-hat 0.974 with 0.63 failing subjects, against
    // -0.221 for normal).  It also contradicted the control's own validation,
    // which requires propMixScale[1] == 1 precisely so component 1 IS the
    // Laplace-approximation covariance.
    //
    // A defensive mixture is deliberately over-dispersed -- that is the
    // mechanism, not a defect.  Its covariance is (sum_k w_k c_k) * gamma *
    // Sigma; gamma still scales the whole mixture, so iscaleMin/iscaleMax
    // still bound it, they just bound a proposal that starts wider than 1.
    for (int k = 0; k < sp.nmix; ++k) {
      sp.c[k] = cv[k];
      sp.w[k] = wv[k] / sw;
    }
  }
  return sp;
}

// Radial/scale multiplier on z.  MVN draws NOTHING (stream unchanged); every
// other family draws exactly one uniform, so the draw count per sample is a
// deterministic function of the family and the stream stays reproducible.
static inline double impPropScale(const impProp& pr, const arma::vec& z) {
  switch (pr.type) {
  case impPropMvn: return 1.0;
  case impPropMvt: return impTScale(pr.df);
  case impPropMvl: {
    double nz = arma::norm(z, 2);
    if (!(nz > 0.0)) return 1.0;
    double r = pr.rScale * impGammaQuantile(rxUnifEng(0.0, 1.0), (double)z.n_elem);
    if (!(r > 0.0) || !R_finite(r)) return 1.0;
    return r / nz;
  }
  default: {
    double u = rxUnifEng(0.0, 1.0);
    int k = 0;
    while (k + 1 < pr.nmix && u > pr.cw[k]) ++k;
    return std::sqrt(pr.c[k]);
  }
  }
}

// -log of the peak-normalized kernel; exactly 0 at quad == 0 for every family.
static inline double impPropLogKernelRecip(const impProp& pr, double quad,
                                           double gammaId, int p) {
  switch (pr.type) {
  case impPropMvn:
  case impPropMvt: return impLogKernelRecip(quad, gammaId, pr.df, p);
  case impPropMvl: return std::sqrt(quad / (gammaId * pr.kappa));
  default: {
    // MIXTURE DENSITY, never the drawn component's -- using the component's
    // would give plausible-looking but wrong weights.  Max-shifted: the naive
    // form underflows to log(0) at large quad and would return +Inf, silently
    // zeroing the very sample the wide component exists to cover.
    double t[3];
    for (int k = 0; k < pr.nmix; ++k) {
      t[k] = pr.lw[k] - quad / (2.0 * pr.c[k] * gammaId);
    }
    return pr.logPeak - impLogSumExp3(t, pr.nmix);
  }
  }
}


// ---- Pareto k-hat (PSIS) ----------------------------------------------------
// In-kernel copy of the Zhang & Stephens (2009) empirical-Bayes GPD fit that
// R/impPsis.R computes post-fit (validated there against loo::psis).  Needed
// HERE because AUTO must react to tail failure while the EM is still running,
// and the cheap in-sample statistics (xi, Kish ESS) are structurally blind to
// it.  Pure arithmetic on the weight vector; no RNG, no allocation beyond the
// tail copy, so it is safe inside the parallel E-step.
static double impPsisKhat(const arma::vec& w) {
  // Mirrors R/impPsis.R line for line (that version is validated against
  // loo::psis); the earlier attempt to be clever with Armadillo expressions
  // silently returned ~0 for weights whose true k-hat was 2.0.
  std::vector<double> x;
  x.reserve(w.n_elem);
  for (arma::uword i = 0; i < w.n_elem; ++i) {
    double v = w[i];
    if (R_finite(v) && v > 0.0) x.push_back(v);
  }
  int S = (int)x.size();
  if (S < 25) return NA_REAL;
  int tail = (int)std::floor(0.2 * S);
  int t2 = (int)std::floor(3.0 * std::sqrt((double)S));
  if (t2 < tail) tail = t2;
  if (tail < 5) return NA_REAL;
  std::sort(x.begin(), x.end());
  // threshold = largest weight NOT in the tail (R: .srt[.s - .tail], 1-based)
  double u = x[S - tail - 1];
  std::vector<double> exc;
  exc.reserve(tail);
  for (int i = S - tail; i < S; ++i) {
    double e = x[i] - u;
    if (e > 0.0) exc.push_back(e);
  }
  int n = (int)exc.size();
  if (n < 5) return NA_REAL;
  std::sort(exc.begin(), exc.end());
  int m = 30 + (int)std::floor(std::sqrt((double)n));
  // R: .xstar <- x[floor(.n/4 + 0.5)]  -- 1-based, so subtract 1 here
  int qi = (int)std::floor((double)n / 4.0 + 0.5) - 1;
  if (qi < 0) qi = 0;
  if (qi >= n) qi = n - 1;
  double xstar = exc[qi];
  if (!(xstar > 0.0) || !R_finite(xstar)) return NA_REAL;
  std::vector<double> theta(m), l(m);
  double lmax = R_NegInf;
  for (int j = 0; j < m; ++j) {
    theta[j] = 1.0 / exc[n - 1] +
      (1.0 - std::sqrt((double)m / ((double)(j + 1) - 0.5))) / 3.0 / xstar;
    double k = 0.0;
    for (int i = 0; i < n; ++i) k += std::log1p(-theta[j] * exc[i]);
    k /= (double)n;
    double lj = R_NegInf;
    double a = -theta[j] / k;
    if (R_finite(k) && k != 0.0 && a > 0.0) lj = (double)n * (std::log(a) - k - 1.0);
    l[j] = lj;
    if (lj > lmax) lmax = lj;
  }
  if (!R_finite(lmax)) return NA_REAL;
  double sw = 0.0, tw = 0.0;
  for (int j = 0; j < m; ++j) {
    double wj = std::exp(l[j] - lmax);
    if (!R_finite(wj)) continue;
    sw += wj; tw += theta[j] * wj;
  }
  if (!(sw > 0.0) || !R_finite(sw)) return NA_REAL;
  double thetaHat = tw / sw;
  double kh = 0.0;
  for (int i = 0; i < n; ++i) kh += std::log1p(-thetaHat * exc[i]);
  return kh / (double)n;
}

// ---- SIR (sampling-importance-resampling) element ---------------------------
// Systematic (stratified) resampling: sirN indices drawn with copy counts
// proportional to the normalized weights (each count is within 1 of
// sirN * zk_norm[i], the systematic-resampling guarantee).  u0 in [0,1) is the
// single stratified offset; no RNG inside.
static arma::uvec impSirIndex(const arma::vec& zk, int sirN, double u0) {
  arma::uvec idx(sirN);
  arma::uword n = zk.n_elem;
  double tot = arma::accu(zk);
  if (!(tot > 0.0) || !R_finite(tot)) {
    // degenerate weights: uniform strided coverage
    for (int r = 0; r < sirN; ++r)
      idx[r] = (arma::uword)((double)r * n / sirN) % n;
    return idx;
  }
  arma::uword i = 0;
  double c = zk[0] / tot;
  for (int r = 0; r < sirN; ++r) {
    double p = (u0 + r) / (double)sirN;
    while (p > c && i + 1 < n) { ++i; c += zk[i] / tot; }
    idx[r] = i;
  }
  return idx;
}

// Test hook: 1-based systematic-resampling indices.
//[[Rcpp::export]]
IntegerVector impSirIndex_(NumericVector zk, int sirN, double u0) {
  if (zk.size() < 1) stop("'zk' must be non-empty");
  if (sirN < 1) stop("'sirN' must be positive");
  if (u0 < 0.0 || u0 >= 1.0) stop("'u0' must be in [0, 1)");
  arma::vec z(zk.begin(), zk.size());
  arma::uvec idx = impSirIndex(z, sirN, u0);
  IntegerVector out(sirN);
  for (int r = 0; r < sirN; ++r) out[r] = (int)idx[r] + 1;
  return out;
}

// Test hook: a proposal family's peak-normalized log-kernel and its Ci
// correction, so the closed forms can be pinned against numeric quadrature
// rather than against themselves.
//[[Rcpp::export]]
NumericVector impPropKernel_(std::string type, double df, NumericVector mixScale,
                             NumericVector mixWeight, double quad, double gamma,
                             int p) {
  if (p < 1) stop("'p' must be positive");
  if (!(gamma > 0.0)) stop("'gamma' must be positive");
  int t = (type == "normal") ? impPropMvn : ((type == "t") ? impPropMvt :
    ((type == "laplace") ? impPropMvl : ((type == "mixture") ? impPropMix : -1)));
  if (t < 0) stop("'type' must be \"normal\", \"t\", \"laplace\" or \"mixture\"");
  int nmix = 0;
  double cc[3] = {1.0, 1.0, 1.0}, ww[3] = {0.0, 0.0, 0.0};
  if (t == impPropMix) {
    nmix = (int)mixScale.size();
    if (nmix < 2 || nmix > 3 || (int)mixWeight.size() != nmix) {
      stop("'mixScale'/'mixWeight' must have matching length 2 or 3");
    }
    double sw = 0.0;
    for (int k = 0; k < nmix; ++k) sw += mixWeight[k];
    for (int k = 0; k < nmix; ++k) { cc[k] = mixScale[k]; ww[k] = mixWeight[k] / sw; }
  }
  impProp pr = impPropBuild(t, df, nmix, cc, ww, p);
  NumericVector out = NumericVector::create(
    _["logKernelRecip"] = impPropLogKernelRecip(pr, quad, gamma, p),
    _["ciCorrection"] = pr.corr,
    _["kappa"] = pr.kappa);
  return out;
}

// Test hook: the (optionally shifted or scrambled) quasi-random N(0,1) point
// set.  scramble is "none"/"owen"/"lms"; a scramble REPLACES the shift, the
// same way the kernel treats them, so supplying both is an error.
//[[Rcpp::export]]
NumericMatrix impQrPoints_(int isample, int neta, Nullable<NumericVector> shift,
                           std::string scramble = "none", int seed = 42) {
  if (isample < 1 || neta < 1) stop("'isample' and 'neta' must be positive");
  int type = (scramble == "none") ? impQrScrambleNone :
    ((scramble == "owen") ? impQrScrambleOwen :
     ((scramble == "lms") ? impQrScrambleLms : -1));
  if (type < 0) stop("'scramble' must be \"none\", \"owen\" or \"lms\"");
  if (type != impQrScrambleNone) {
    if (shift.isNotNull()) stop("'shift' does not apply when scrambling");
    std::vector<uint64_t> V = impSobolBits(isample, neta);
    return wrap(impQrZScrambled(V, isample, neta, type, (uint32_t)seed, 0, 0));
  }
  arma::mat U0 = impSobolU0(isample, neta);
  if (shift.isNotNull()) {
    NumericVector s(shift);
    if ((int)s.size() != neta) stop("'shift' must have length 'neta'");
    arma::vec sv(s.begin(), neta);
    return wrap(impQrZ(U0, &sv));
  }
  return wrap(impQrZ(U0, nullptr));
}

// impGetHessian() solves through odeSwapSolveInd(), which writes rx->ndiff -- a
// field on the single shared rx_solve struct -- right before every solve of a
// slot whose OWN need is nonzero (see odeSwap.cpp).  A subject that falls back
// to the doFD/pred path (fInd->doFD, set per-subject) solves a DIFFERENT slot
// than a subject still on the plain inner path, concurrently in the E-step's
// parallel loop -- but odeSwapSolveInd() already makes that pairing benign
// (pred needs no derivative, so it never writes).  This lock is
// defense-in-depth for the case that alone does not cover: two DIFFERENT slots
// that both need a NONZERO but DIFFERING ndiff, concurrently.  Only a
// linCmt()-with-non-mu-structural-theta model ever has a slot with
// ndiffSet=true at all, so gate the lock on that and pay nothing otherwise.
static bool impGetHessianNdiffSafe(int id, arma::mat& H) {
  if (!odeSwapAnyNdiffSet()) return impGetHessian(id, H);
  bool gotH;
#ifdef _OPENMP
#pragma omp critical (odeSwapNdiffHessian)
#endif
  { gotH = impGetHessian(id, H); }
  return gotH;
}

// Per-expanded-subject proposal (mode, information H_i, log|H_i|, lower Cholesky
// of gamma * H_i^-1) for one E-step.
struct ImpEProp {
  std::vector<arma::vec> modes;
  std::vector<arma::mat> Hs;
  std::vector<double> logDetH;
  std::vector<arma::mat> cholL;
  std::vector<char> haveL;
  explicit ImpEProp(int n) : modes(n), Hs(n), logDetH(n, 0.0), cholL(n), haveL(n, 0) {}
};

// Proposal from a covariance Sig (information Sig^-1); true when id's slot is filled
static bool impPropFromSigma(ImpEProp& P, int id, const arma::mat& Sig, double gammaId) {
  arma::mat Hi;
  double ldv, lds;
  if (!(arma::inv_sympd(Hi, Sig) && arma::log_det(ldv, lds, Hi) && lds > 0)) return false;
  arma::mat L;
  if (!arma::chol(L, gammaId * Sig, "lower")) return false;
  P.Hs[id] = Hi; P.logDetH[id] = ldv; P.cholL[id] = L; P.haveL[id] = 1;
  return true;
}

// Proposal from the MAP information H (the curvature is kept even if chol fails)
static void impPropFromHessian(ImpEProp& P, int id, const arma::mat& H, double gammaId) {
  arma::mat Sigma;
  double ldv, lds;
  if (!(arma::inv_sympd(Sigma, H) && arma::log_det(ldv, lds, H) && lds > 0)) return;
  P.Hs[id] = H;
  P.logDetH[id] = ldv;
  Sigma *= gammaId;
  arma::mat L;
  if (arma::chol(L, Sigma, "lower")) { P.cholL[id] = L; P.haveL[id] = 1; }
}

// est="imp": no MAP search.  The proposal is centered at each subject's running
// conditional mean (the current eta) with covariance gamma * V_i, where V_i is
// that subject's conditional variance from the PREVIOUS E-step (passed in via
// condVar, before it is overwritten).  On the first iteration -- or a
// pseudo-subject / degenerate V_i -- fall back to the population Omega, which is
// over-dispersed and always available.
static void impPropImp(ImpEProp& P, int id, int neta, const std::vector<arma::mat>& condVar,
                       const arma::mat& impOmega, double gammaId) {
  arma::vec eta(neta);
  impGetEta(id, eta);
  P.modes[id] = eta;
  const bool useCv = (size_t)id < condVar.size() && condVar[id].n_rows == (arma::uword)neta &&
    arma::any(arma::vectorise(condVar[id]) != 0.0);
  impPropFromSigma(P, id, useCv ? condVar[id] : impOmega, gammaId);
}

// est="impmap": the proposal is the MAP mode and Hessian.
static void impPropMap(ImpEProp& P, int id, int neta, const arma::mat& impOmega,
                       double gammaId, std::atomic<int>& nOmegaFallback) {
  arma::vec mode(neta);
  arma::mat H(neta, neta, arma::fill::zeros);
  impGetMode(id, mode);
  P.modes[id] = mode;
  if (impGetHessianNdiffSafe(id, H)) impPropFromHessian(P, id, H, gammaId);
  if (P.haveL[id]) return;
  // The MAP Hessian did not yield a proposal for this subject.  The
  // Hessian is the curvature AT the mode, so it is only guaranteed
  // positive definite at a converged one -- and the inner solve is not
  // required to have converged (a pass that exhausted maxInnerIterations
  // still returns a finite objective and can be the eta reported).  At
  // such a point inv_sympd()/chol() legitimately fail.
  //
  // Dropping the subject here is the wrong answer to that: it silently
  // removes an individual from the E-step, and if it happens widely the
  // whole fit collapses to the failure objective with nothing said.
  // est="imp" already meets its own version of this by centering on the
  // population Omega, which is over-dispersed and always available; do
  // the same.  An over-dispersed proposal costs effective sample size --
  // the importance weights correct for it -- where a missing subject
  // costs correctness.
  if (impPropFromSigma(P, id, impOmega, gammaId)) {
    nOmegaFallback.fetch_add(1, std::memory_order_relaxed);
  }
}

// Every expanded subject's proposal from the per-component MAP modes.  Each
// subject fills only its own slot (no reduction, no RNG), so this is parallelized
// over subjects; the per-subject inner Hessian (impGetHessian) is the impmap
// serial cost.  cores already carries the mixture / pool-sizing serial guard.
static void impEStepProposals(ImpEProp& P, int nExp, int neta, const arma::vec& gammaVec,
                              int cores, bool isImp, const std::vector<arma::mat>& condVar) {
  // Subjects whose MAP Hessian gave no proposal and fell back to the population
  // Omega.  Atomic because the fill loop below is parallel over subjects.
  std::atomic<int> nOmegaFallback(0);
  // Needed by BOTH branches: est="imp" centers on it when a subject has no
  // usable conditional variance yet, and est="impmap" falls back to it when the
  // MAP Hessian will not give a proposal.
  arma::mat impOmega;
  impGetOmega(impOmega);
  bool doParProp = (cores > 1);
  nmForEachSubject(rx, nExp, cores, doParProp, [&](int id) {
    if (isImp) {
      impPropImp(P, id, neta, condVar, impOmega, gammaVec[id]);
    } else {
      impPropMap(P, id, neta, impOmega, gammaVec[id], nOmegaFallback);
    }
  });
  // Say so.  A proposal built on the population Omega rather than on the
  // subject's own curvature is a weaker proposal, and the reason it happened
  // (an unconverged MAP) is worth knowing about; silently substituting it is
  // how a fit ends up quietly worse than it looks.
  if (nOmegaFallback.load(std::memory_order_relaxed) > 0) {
    RSprintf("impmap: %d of %d subject proposals fell back to the population omega (MAP Hessian not usable)\n",
             nOmegaFallback.load(std::memory_order_relaxed), nExp);
  }
}

// QRPEM: one Sobol base point set per E-step (read-only in the parallel loop);
// with qrShift=FALSE the N(0,1) points are also fixed and shared.
struct ImpQrSet {
  bool on = false;
  int scramble = impQrScrambleNone;
  bool shift = false;
  bool refresh = true;
  int nrow = 0;
  arma::mat U0, Z0;
  std::vector<uint64_t> bits;
};

static ImpQrSet impQrSetup(int nrow, int neta) {
  ImpQrSet q;
  q.on = impQrEnabled();
  q.scramble = impQrScramble();
  // Scrambling IS the randomization, so it replaces the Cranley-Patterson
  // shift rather than composing with it; qrRefresh still says whether the
  // randomization is redrawn each iteration.
  q.shift = (q.scramble == impQrScrambleNone) && impQrShiftEnabled();
  q.refresh = impQrRefreshEnabled();
  if (!q.on) return q;
  q.nrow = nrow;
  if (q.scramble == impQrScrambleNone) {
    q.U0 = impSobolU0(nrow, neta);
    if (!q.shift) q.Z0 = impQrZ(q.U0, nullptr);
  } else {
    q.bits = impSobolBits(nrow, neta);
  }
  return q;
}

// Sobol points -> N(0,I) for one subject's E-step.  With qrShift the shift comes
// from this subject's seeded stream; qrRefresh=FALSE instead uses the
// fit-constant qrPinSeed (captured once in impOuter -- seed0 advances every
// E-step) so the whole fit reuses one shift per subject.
static arma::mat impQrZSubject(const ImpQrSet& q, int nsId, int neta, uint32_t qrPinSeed,
                               int id, int iter) {
  if (q.scramble != impQrScrambleNone) {
    // qrRefresh=FALSE pins the scramble to the subject, so every
    // iteration sees the same scrambled set (a deterministic EM map).
    arma::mat Z = impQrZScrambled(q.bits, q.nrow, neta, q.scramble, qrPinSeed, id,
                                  q.refresh ? iter : 0);
    if (nsId < q.nrow) Z = Z.rows(0, nsId - 1);
    return Z;
  }
  if (!q.shift) return q.Z0;
  if (!q.refresh) nmSetSeedEng1(qrPinSeed + (uint32_t)(id * 2));
  arma::vec sh(neta);
  for (int jj = 0; jj < neta; ++jj) sh[jj] = rxUnifEng(0.0, 1.0);
  return impQrZ(q.U0, &sh);
}

// est="its": 2*neta equally weighted sigma points, mode +/- sqrt(neta) L e_j
static void impSigmaPoints(arma::mat& S, const arma::vec& mode, const arma::mat& L, int neta) {
  const double sq = std::sqrt((double)neta);
  for (int jj = 0; jj < neta; ++jj) {
    S.row(2 * jj) = (mode + sq * L.col(jj)).t();
    S.row(2 * jj + 1) = (mode - sq * L.col(jj)).t();
  }
}

// Quasi-random samples from the N(0,I) points Z
static void impQrSamples(arma::mat& S, const arma::mat& Z, const arma::vec& mode,
                         const arma::mat& L, const impProp& pr, int nsId) {
  // The MVN/MVT branch is kept TEXTUALLY as it was: the two forms should
  // compile to the same gemv, but "should" is not the bit-identity
  // contract this path holds to.
  if (pr.type <= impPropMvt) {
    for (int k = 0; k < nsId; ++k) {
      // Quasi-random t: the Sobol point supplies the Gaussian direction and
      // a fresh uniform supplies the chi-square scale.  (A fully
      // quasi-random t would need an extra Sobol coordinate; this keeps the
      // low-discrepancy structure where it does the most good.)
      double ts = impTScale(pr.df);
      S.row(k) = (mode + L * (Z.row(k).t() * ts)).t();
    }
  } else {
    for (int k = 0; k < nsId; ++k) {
      arma::vec z = Z.row(k).t();
      double ts = impPropScale(pr, z);
      S.row(k) = (mode + L * (z * ts)).t();
    }
  }
}

// Pseudo-random samples from this subject's seeded stream
static void impRngSamples(arma::mat& S, const arma::vec& mode, const arma::mat& L,
                          const impProp& pr, int nsId, int neta) {
  for (int k = 0; k < nsId; ++k) {
    arma::vec z(neta);
    for (int jj = 0; jj < neta; ++jj) z[jj] = rxNormEng(0.0, 1.0);
    // evaluated AFTER the z fill, so the draw order (and therefore the
    // stream position) is unchanged for MVN/MVT
    double ts = impPropScale(pr, z);
    S.row(k) = (mode + L * (z * ts)).t();
  }
}

// #958 combined build: impEvalJointLik() just solved this subject's
// INNER model at this exact (theta, eta) -- when combSens is loaded,
// that solve already carries the M-step's d(f)/d(theta)/d(V)/d(theta)
// columns, so harvest them here (impThetaSensCollect's reuseSolve path,
// a pure read -- no ODE solve, no R API call, safe under the E-step's
// own parallel region) instead of letting the M-step re-solve every
// sample from scratch.  Skipped whenever the caller disabled harvesting
// (SIR resampling only needs a subset of the E-step's samples, so eagerly
// harvesting all of them there would be wasted work -- see impOuter's
// harvestSens gate).  A sample whose solve fell back to the pred-model FD
// path (impLastInnerSolveUsable == false) leaves ind->solve holding the
// PRED model's trajectory, not the inner model's, so it is left
// unharvested (sampleOk=0) rather than triggering impThetaSensCollect's
// own solve-with-retry path -- that path is R-API-unsafe from a bare
// OpenMP region (it is only ever run under the M-step's explicit
// impInnerParallelOn() guard), which this loop does not set.  doFD alone
// is not enough: it is a STICKY per-subject flag (only flips permanently,
// inside the inner Newton optimizer's own convergence handling) and does
// not track whether THIS ARBITRARY sampled eta's solve produced a usable
// result, so a bad solve at doFD==0 (e.g. an extreme early-EM proposal
// draw) would otherwise still reach impThetaSensCollect, which falls back
// to impThetaSensFD()'s finite-difference re-solves -- the same
// R-API-unsafe path.  qk is -impEvalJointLik(eta, id) + <finite kernel
// term>, so its own finiteness is direct evidence this exact eta's inner
// solve succeeded (impEvalJointLik/likInner0 return a non-finite value on
// a bad solve); require both.
static void impHarvestSample(impThetaSensData& acc, int id, const arma::vec& eta, double qk,
                             int neta) {
  // Every sample appends exactly one slot; an unharvested one appends an
  // empty placeholder with sampleOk=0 so the per-sample indexing lines up.
  impThetaSensData c;
  bool ok = false;
  if (impLastInnerSolveUsable(id) && R_finite(qk)) {
    arma::mat srow(1, neta);
    srow.row(0) = eta.t();
    impThetaSensCollect(id, srow, c, /*reuseSolve=*/true);
    if (c.nobs > 0) {
      acc.nobs = c.nobs;
      acc.dvv = c.dvv; acc.limv = c.limv;
      acc.censv = c.censv; acc.distv = c.distv;
      acc.ddvmat = c.ddvmat; acc.dlimmat = c.dlimmat;
    }
    ok = c.nobs > 0 && !c.sampleOk.empty() && c.sampleOk[0] != 0;
  }
  acc.fvec.push_back(ok ? c.fvec[0] : arma::vec());
  acc.Vvec.push_back(ok ? c.Vvec[0] : arma::vec());
  acc.dfmat.push_back(ok ? c.dfmat[0] : arma::mat());
  acc.dVmat.push_back(ok ? c.dVmat[0] : arma::mat());
  acc.sampleOk.push_back(ok ? 1 : 0);
  if (ok) impThetaSensHarvestTick();
}

// Per-expanded-subject E-step results (combined per base subject afterwards).
struct ImpEOut {
  arma::mat cm;
  std::vector<arma::mat> cv;
  arma::vec Li, Neff, Xi, Khat;
  std::vector<char> ok;
  ImpEOut(int nExp, int neta)
    : cm(nExp, neta, arma::fill::zeros),
      cv(nExp, arma::mat(neta, neta, arma::fill::zeros)),
      Li(nExp), Neff(nExp), Xi(nExp), Khat(nExp), ok(nExp, 0) {
    Li.fill(NA_REAL); Neff.fill(NA_REAL); Xi.fill(NA_REAL); Khat.fill(NA_REAL);
  }
};

// Everything one subject's E-step reads or writes.
struct ImpECtx {
  int neta, nExp, iter;
  const arma::ivec& isampleVec;
  const arma::vec& gammaVec;
  const std::vector<impProp>& props;
  double negHalfLogDetOmega;
  bool its, harvestSens;
  uint32_t seed0, qrPinSeed;
  const ImpEProp& P;
  const ImpQrSet& qr;
  ImpEOut& O;
  std::vector<arma::mat>& outS;
  std::vector<arma::vec>& outZk;
  std::vector<impThetaSensData>& outSens;
};

// Fill the subject's sample matrix S
static void impEStepDraw(const ImpECtx& c, int id, arma::mat& S, int nsId) {
  const impProp& pr = c.props[id];
  if (c.its) {
    impSigmaPoints(S, c.P.modes[id], c.P.cholL[id], c.neta);
  } else if (c.qr.on) {
    arma::mat Z = impQrZSubject(c.qr, nsId, c.neta, c.qrPinSeed, id, c.iter);
    impQrSamples(S, Z, c.P.modes[id], c.P.cholL[id], pr, nsId);
  } else {
    impRngSamples(S, c.P.modes[id], c.P.cholL[id], pr, nsId, c.neta);
  }
}

// log importance weights; returns the number of usable samples.
// impEvalJointLik returns the NEGATIVE log joint density (minimized by the
// inner optimizer), so it enters with a minus; the proposal density is
// subtracted, contributing +(1/(2 gamma)) d'H d.
static int impEStepLogWeights(const ImpECtx& c, int id, const arma::mat& S, arma::vec& q) {
  const impProp& pr = c.props[id];
  const int nsId = (int)S.n_rows;
  int nGood = 0;
  for (int k = 0; k < nsId; ++k) {
    arma::vec eta = S.row(k).t();
    arma::vec d = eta - c.P.modes[id];
    double qk = -impEvalJointLik(eta, id) +
      impPropLogKernelRecip(pr, arma::as_scalar(d.t() * c.P.Hs[id] * d), c.gammaVec[id], c.neta);
    if (c.harvestSens) impHarvestSample(c.outSens[id], id, eta, qk, c.neta);
    // A proposal sample whose inner solve fails (NA/NaN, even after the
    // FOCEI tolerance-relaxation retry) is dropped from the weighted mean by
    // giving it -Inf log-weight (weight 0).
    if (R_finite(qk)) { q[k] = qk; ++nGood; } else { q[k] = R_NegInf; }
  }
  return nGood;
}

// Weights, objective contribution, diagnostics and conditional moments of one
// subject from its log weights q.
static void impEStepMoments(const ImpECtx& c, int id, const arma::mat& S, const arma::vec& q,
                            double qCenter) {
  const impProp& pr = c.props[id];
  const int nsId = (int)S.n_rows;
  const int neta = c.neta;
  double gammaId = c.gammaVec[id];
  double qmax = q.max();
  arma::vec w = arma::exp(q - qmax);
  double sumw = arma::accu(w);
  arma::vec zk = w / sumw;
  c.outS[id] = S;
  c.outZk[id] = zk;
  c.O.ok[id] = 1;
  double logMeanExp = c.its ? qCenter : qmax + std::log(sumw / (double)nsId);
  double Ci = c.negHalfLogDetOmega + 0.5 * neta * std::log(gammaId) -
    0.5 * c.P.logDetH[id] + pr.corr;
  c.O.Li[id] = -(logMeanExp + Ci);          // per-component negative log-likelihood
  c.O.Neff[id] = 1.0 / arma::accu(arma::square(zk));
  // NONMEM xi: mean weight normalized to the proposal center (see impEStepSubject).
  // xi is NOT family-invariant: against a Gaussian target it comes out as
  // gamma^(-p/2) * exp(-corr), so a family with corr != 0 would drive the
  // controller (which inverts xi = gamma^(-neta/2) analytically) to
  // a DIFFERENT gamma -- for covariance-matched MVL at p=3 that is ~3x too
  // wide, or pinned at iscaleMax, with gammaStable still reporting
  // converged.  Adding corr back puts every family on gamma^(-p/2), so
  // iaccept / iscaleMin / iscaleMax / pExp keep their tuning verbatim.
  //
  // Gated off for MVN and MVT.  MVN's corr is exactly 0.0 so it is a no-op
  // there either way; MVT's is not (2.5% at df=30), and moving xi moves
  // gamma and therefore the whole fit -- the AUTO ladder's df 30/20
  // constants were tuned with this term absent.  By the same token MVT is
  // not covariance-matched either (Cov = df/(df-2) Sigma); both are
  // deliberate, not oversights.
  double xiCorr = (pr.type <= impPropMvt) ? 0.0 : pr.corr;
  if (R_finite(qCenter)) c.O.Xi[id] = std::exp(logMeanExp - qCenter + xiCorr);
  c.O.Khat[id] = impPsisKhat(zk);   // tail index of this subject's weights
  arma::rowvec pbar = zk.t() * S;
  c.O.cm.row(id) = pbar;
  arma::mat B(neta, neta, arma::fill::zeros);
  for (int k = 0; k < nsId; ++k) {
    arma::rowvec dr = S.row(k) - pbar;
    B += zk[k] * (dr.t() * dr);
  }
  if (c.its) {
    // Hs is H itself (gamma = 1) or, after a fallback, Omega^-1
    c.O.cm.row(id) = c.P.modes[id].t();
    arma::mat V;
    if (arma::inv_sympd(V, c.P.Hs[id])) B = V;
  }
  c.O.cv[id] = B;
}

// One expanded subject's E-step.
static void impEStepSubject(const ImpECtx& c, int id) {
  // fresh per-(iter, expanded-subject) stream, independent of thread count
  nmSetSeedEng1(c.seed0 + (uint32_t)((c.iter * c.nExp + id) * 2));
  c.O.cm.row(id) = c.P.modes[id].t();
  if (!c.P.haveL[id]) return;
  // est="its": no sampling.  The moments are the mode and H^-1, the objective is
  // the Laplace term, and the M-step sees 2*neta equally weighted sigma points
  // (mode +/- sqrt(neta) L e_j) so the residual update carries H^-1 as well.
  int nsId = c.its ? 2 * c.neta : (int)c.isampleVec[id];   // this subject's sample count
  arma::mat S(nsId, c.neta);
  impEStepDraw(c, id, S, nsId);
  // NONMEM-style xi (NM7 eq. 1.73/1.75): the mean importance weight
  // normalized so the proposal kernel is 1 at its center.  q at the center
  // is just the log joint density there (the kernel term vanishes at d=0),
  // so xi = exp(logMeanExp - qCenter) once logMeanExp is formed.
  //
  // xi == 1 exactly when the proposal matches the posterior shape, < 1 when
  // the proposal is over-dispersed, and > 1 when the posterior has heavier
  // tails than the Gaussian proposal -- xi is NOT bounded above.  That
  // heavy-tail signal is what the per-subject gamma controller
  // targets, and it is a different statistic from the Kish
  // effective sample size (they are not comparable).
  //
  // Evaluated HERE -- after the draws are complete but before the weight
  // loop -- deliberately: the sample matrix S is already filled, so this
  // cannot perturb the sampling stream, and the weight loop still leaves
  // the subject's eta at S.row(isample-1) exactly as before, preserving
  // bit-identity with the pre-xi kernel.
  double qCenter = -impEvalJointLik(c.P.modes[id], id);
  arma::vec q(nsId);
  int nGood = impEStepLogWeights(c, id, S, q);
  // the inner solves above re-seed the engine per subject
  // (setSeedEng1(getRxSeed1()+id)); restore this block's sampling seed so
  // that leak never carries into a subsequent draw.
  nmRestoreMcmcSeed();
  if (nGood == 0) return;
  if (c.its) {
    // equal weights (failed points dropped); Laplace objective at the mode
    if (!R_finite(qCenter)) return;
    for (int k = 0; k < nsId; ++k) q[k] = R_finite(q[k]) ? qCenter : R_NegInf;
  }
  impEStepMoments(c, id, S, q, qCenter);
}

// Combined (per base subject) E-step results.
struct ImpEDest {
  arma::mat& condMean;
  std::vector<arma::mat>& condVar;
  arma::vec& Li;
  arma::vec& Neff;
  arma::vec& Xi;
  arma::mat& aMat;
};

// Smallest usable component objective of base subject i (+Inf when none)
static double impMixMinL(const ImpEOut& O, int i, int nsub, int Nm) {
  double minL = R_PosInf;
  for (int j = 0; j < Nm; ++j) {
    int eid = i + j * nsub;
    if (O.ok[eid] && R_finite(O.Li[eid]) && O.Li[eid] < minL) minL = O.Li[eid];
  }
  return minL;
}

// Unnormalized posterior mixture weights a_j exp(minL - L_ij); returns their sum
static double impMixWeights(const ImpEOut& O, int i, int nsub, int Nm, double minL,
                            arma::vec& wj) {
  double sumW = 0.0;
  for (int j = 0; j < Nm; ++j) {
    int eid = i + j * nsub;
    if (O.ok[eid] && R_finite(O.Li[eid])) {
      wj[j] = impMixProb(j) * std::exp(minL - O.Li[eid]);
      sumW += wj[j];
    }
  }
  return sumW;
}

// Responsibility-weighted Neff and xi of base subject i
static void impMixDiag(const ImpEOut& O, int i, int nsub, int Nm, ImpEDest& D) {
  double nf = 0.0;
  for (int j = 0; j < Nm; ++j) {
    double a = D.aMat(i, j);
    if (a > 0 && R_finite(O.Neff[i + j * nsub])) nf += a * O.Neff[i + j * nsub];
  }
  if (nf > 0) D.Neff[i] = nf;
  // Responsibility-weighted xi, matching how Neff is combined above.
  double xf = 0.0; bool anyXi = false;
  for (int j = 0; j < Nm; ++j) {
    double a = D.aMat(i, j);
    if (a > 0 && R_finite(O.Xi[i + j * nsub])) { xf += a * O.Xi[i + j * nsub]; anyXi = true; }
  }
  if (anyXi) D.Xi[i] = xf;
}

// Responsibility-combined moments of base subject i from its components
static void impMixMoments(const ImpEOut& O, int i, int nsub, int Nm, int neta,
                          const arma::vec& wj, double sumW, ImpEDest& D) {
  arma::rowvec cm(neta, arma::fill::zeros);
  for (int j = 0; j < Nm; ++j) {
    double a = wj[j] / sumW; D.aMat(i, j) = a;
    if (a > 0) cm += a * O.cm.row(i + j * nsub);
  }
  D.condMean.row(i) = cm;
  // Component-weighted conditional variance: sum_j a_ij (B_ij + m_ij m_ij') - m_i m_i'.
  arma::mat cv(neta, neta, arma::fill::zeros);
  for (int j = 0; j < Nm; ++j) {
    double a = D.aMat(i, j); if (a <= 0) continue;
    arma::rowvec mj = O.cm.row(i + j * nsub);
    cv += a * (O.cv[i + j * nsub] + mj.t() * mj);
  }
  cv -= cm.t() * cm;
  D.condVar[i] = cv;
}

// Combine base subject i over its mixture components.  Posterior mixture weights
// a_ij = a_j exp(-L_ij) / sum_k a_k exp(-L_ik), with the min-L shift for
// numerical stability (as saem does).
static void impMixCombine(const ImpEOut& O, int i, int nsub, int Nm, int neta,
                          std::vector<arma::vec>& outZk, ImpEDest& D) {
  double minL = impMixMinL(O, i, nsub, Nm);
  if (!R_finite(minL)) { D.condMean.row(i) = O.cm.row(i); D.aMat(i, 0) = 1.0; return; }
  arma::vec wj(Nm, arma::fill::zeros);
  double sumW = impMixWeights(O, i, nsub, Nm, minL, wj);
  if (sumW <= 0.0) { D.condMean.row(i) = O.cm.row(i); D.aMat(i, 0) = 1.0; return; }
  impMixMoments(O, i, nsub, Nm, neta, wj, sumW, D);
  D.Li[i] = minL - std::log(sumW);           // marginal (mixture) negative log-likelihood
  impMixDiag(O, i, nsub, Nm, D);
  // Fold the responsibility a_ij into each component's importance weights so the
  // M-step, iterating the expanded subjects, forms the component-weighted score.
  for (int j = 0; j < Nm; ++j) {
    int eid = i + j * nsub;
    if (outZk[eid].n_elem > 0) outZk[eid] *= D.aMat(i, j);
  }
}

// Combine over mixture components per base subject.  For a single component
// (Nm == 1) this is bit-identical to the non-mixture E-step.
static void impEStepCombine(const ImpEOut& O, int nsub, int Nm, int neta,
                            std::vector<arma::vec>& outZk, ImpEDest& D) {
  D.condMean.set_size(nsub, neta);
  D.condVar.assign(nsub, arma::mat(neta, neta, arma::fill::zeros));
  D.Li.set_size(nsub); D.Li.fill(NA_REAL);
  D.Neff.set_size(nsub); D.Neff.fill(NA_REAL);
  D.Xi.set_size(nsub); D.Xi.fill(NA_REAL);
  D.aMat.set_size(nsub, Nm); D.aMat.zeros();
  for (int i = 0; i < nsub; ++i) {
    if (Nm != 1) {
      impMixCombine(O, i, nsub, Nm, neta, outZk, D);
      continue;
    }
    D.condMean.row(i) = O.cm.row(i);
    D.condVar[i] = O.cv[i];
    D.Li[i] = O.Li[i];
    D.Neff[i] = O.Neff[i];
    D.Xi[i] = O.Xi[i];
    D.aMat(i, 0) = 1.0;
  }
}

// On the last iteration, stash the per-subject sampler diagnostics (first
// component for a mixture).
static void impEStepStash(Environment& eStash, const ImpEProp& P, int nsub, int neta,
                          const std::vector<arma::mat>& outS,
                          const std::vector<arma::vec>& outZk) {
  arma::mat modeMat(nsub, neta);
  List hessList(nsub), samplesList(nsub), weightList(nsub);
  for (int id = 0; id < nsub; ++id) {
    modeMat.row(id) = P.modes[id].t();
    hessList[id] = P.haveL[id] ? wrap(P.Hs[id]) : R_NilValue;
    samplesList[id] = wrap(outS[id]);
    // Normalized importance weights, kept so a TAIL-sensitive diagnostic
    // (Pareto k-hat) can be computed post-fit.  xi and the Kish ESS are both
    // means over samples drawn FROM the proposal, so both are structurally
    // blind to tails the proposal never visits -- which is the failure mode a
    // heavier-tailed (t) proposal exists to prevent.  k-hat estimates the
    // tail index of the weight distribution and does see it.  Pareto k-hat is
    // scale-invariant, so the mixture responsibility folded into these
    // weights does not affect it.
    weightList[id] = wrap(outZk[id]);
  }
  eStash["impEtaMode"] = wrap(modeMat);
  eStash["impEtaHess"] = hessList;
  eStash["impSamples"] = samplesList;
  eStash["impWeights"] = weightList;
}

// One importance-sampling E-step at the current conditional modes: draw
// proposal samples, form importance weights, and return each subject's
// conditional mean, variance, objective contribution, and effective sample
// size.  `iter` shifts the per-subject RNG stream so successive iterations use
// fresh (still thread-count-independent) samples.
// gammaVec carries ONE proposal scale per EXPANDED subject (id = i + j*nsub), so
// a mixture's components each get their own scale.  Under gammaMethod="global"
// every element holds the same value and this is arithmetically identical to the
// previous scalar gamma; the per-subject controller is what makes them diverge.
// isampleVec carries one sample count per EXPANDED subject.  Under the classic
// (scalar) setting every entry is the same and this is arithmetically identical
// to the previous fixed-count E-step; per-subject counts let a badly covered
// subject buy more samples without charging every other subject for them, which
// is the NM7 Technical Guide's own mechanism (its derivation is Gaussian
// throughout and never mentions a t proposal).
static void impEStep(int nsub, int neta, const arma::ivec& isampleVec,
                     const arma::vec& gammaVec,
                     const std::vector<impProp>& props, int cores,
                     int iter, double negHalfLogDetOmega, bool isImp,
                     arma::mat& condMean, std::vector<arma::mat>& condVar,
                     arma::vec& Li, arma::vec& Neff, arma::vec& Xi,
                     arma::vec& NeffExpOut, arma::vec& XiExpOut, arma::vec& KhatExpOut,
                     std::vector<arma::mat>& outS, std::vector<arma::vec>& outZk,
                     arma::mat& aMat, Environment* eStash,
                     uint32_t qrPinSeed, bool harvestSens,
                     std::vector<impThetaSensData>& outSens) {
  int Nm = impNmix();
  int nExp = nsub * Nm; // expanded pseudo-subjects: component j (1-based) is i + (j-1)*nsub
  ImpEProp P(nExp);
  impEStepProposals(P, nExp, neta, gammaVec, cores, isImp, condVar);

  ImpEOut O(nExp, neta);
  outS.assign(nExp, arma::mat());
  outZk.assign(nExp, arma::vec());
  if (harvestSens) outSens.assign(nExp, impThetaSensData());
  // The engine array is allocated + seeded by rxWithSeed() at the fit start
  // (impmap.R), so we don't call seedEng() here.
  // Base the per-subject stream on the control's impSeed, NOT the ambient
  // rxSeed: the FOCEI solve setup resets rxSeed to its default before the
  // E-step runs, so getRxSeed1() would ignore impSeed entirely.  The base is
  // salted per consumer (E-step / QR shift / SIR / covariance) and offset per
  // (iteration, expanded-subject) so every step draws an independent but
  // reproducible, thread-count-independent stream.
  uint32_t seed0 = (uint32_t)impBaseSeed();
  bool doPar = (cores > 1);
  ImpQrSet qr = impQrSetup((int)isampleVec.max(), neta);
  ImpECtx c{neta, nExp, iter, isampleVec, gammaVec, props, negHalfLogDetOmega,
            impIsIts(), harvestSens, seed0, qrPinSeed, P, qr, O, outS, outZk, outSens};
  // Parallelize over base subjects, iterating a subject's mixture components
  // serially within one thread.  A base subject's expanded pseudo-subjects (id =
  // i + j*nsub) share the same underlying rxode2 solve rows, so evaluating two of
  // its components concurrently would race on that solve buffer -- keeping them on
  // one thread avoids the race while preserving nsub-way parallelism.  For Nm == 1
  // this is the same nsub-iteration loop as before (bit-identical, same seeds).
  nmForEachSubject(rx, nsub, cores, doPar, [&](int i) {
    for (int j = 0; j < Nm; ++j) impEStepSubject(c, i + j * nsub);
  });

  ImpEDest D{condMean, condVar, Li, Neff, Xi, aMat};
  impEStepCombine(O, nsub, Nm, neta, outZk, D);

  // Per-EXPANDED-subject xi, which is what the individual-gamma controller
  // drives (one gamma_i per expanded subject).  `Xi` above is the
  // responsibility-combined per-base-subject value used for reporting.  The
  // same holds for the effective sample size: AUTO reallocates samples per
  // expanded subject, so it needs NeffExp, not the per-base-subject `Neff`.
  NeffExpOut = O.Neff;
  XiExpOut = O.Xi;
  KhatExpOut = O.Khat;
  if (eStash != nullptr) impEStepStash(*eStash, P, nsub, neta, outS, outZk);
}

// Monte-Carlo observed-information covariance for the estimated thetas and Omega.
//
// At the converged estimate this recomputes each subject's proposal (mode +
// Hessian), draws one fixed set of importance samples, and takes the finite-
// difference Hessian of the importance-sampling -2LL with respect to the thetas
// and the parameterized Omega elements (an Omega perturbation rebuilds omegaInv
// and the -0.5 log|Omega| normalizer).  Because the SAME samples are reused for
// every perturbation (common random numbers) the reweighted objective is a
// deterministic smooth function of the parameters, so the FD Hessian is well
// behaved.  The observed information is 0.5 * d2(-2LL)/dpar2 and the covariance
// is its inverse.  (mu-referenced thetas still need mode tracking -- a later
// increment -- since the fixed samples do not follow the mode shift.)
// gammaVec carries the CONVERGED proposal scales from the EM (one per expanded
// subject; the first nsub entries are component 0, which is the indexing this
// function uses).  It previously read impGammaProp(), the control's *initial*
// gamma -- wrong whenever the scale adapted during the fit, and badly wrong
// under gammaMethod="individual" where the scales are per subject.  The
// covariance must be evaluated with the same proposal the fit converged on.
// Central finite-difference Hessian of `f` at `par0`, with the same relative
// step the MC covariance has always used.  Split out of impComputeCov() so the
// FD bookkeeping is readable on its own; `f` reuses the fixed importance
// samples (common random numbers), which is what makes this well behaved.
static arma::mat impFdHessian(const arma::vec& par0,
                              const std::function<double(const arma::vec&)>& f) {
  const int np = (int)par0.n_elem;
  arma::vec hstep(np);
  for (int j = 0; j < np; ++j) {
    double a = std::fabs(par0[j]);
    hstep[j] = 1e-3 * (a > 1e-3 ? a : 1.0);
  }
  double f0 = f(par0);
  arma::mat hess(np, np, arma::fill::zeros);
  for (int j = 0; j < np; ++j) {
    arma::vec p = par0; p[j] = par0[j] + hstep[j]; double fp = f(p);
    p = par0; p[j] = par0[j] - hstep[j]; double fm = f(p);
    hess(j, j) = (fp - 2.0 * f0 + fm) / (hstep[j] * hstep[j]);
  }
  for (int a = 0; a < np; ++a) {
    for (int b = a + 1; b < np; ++b) {
      arma::vec p = par0; p[a] += hstep[a]; p[b] += hstep[b]; double fpp = f(p);
      p = par0; p[a] += hstep[a]; p[b] -= hstep[b]; double fpm = f(p);
      p = par0; p[a] -= hstep[a]; p[b] += hstep[b]; double fmp = f(p);
      p = par0; p[a] -= hstep[a]; p[b] -= hstep[b]; double fmm = f(p);
      double v = (fpp - fpm - fmp + fmm) / (4.0 * hstep[a] * hstep[b]);
      hess(a, b) = v; hess(b, a) = v;
    }
  }
  return hess;
}

// src/inner.cpp: the covariance probe tolerances on the live solve, and back
void covProbeSolveTolPush(double *sav);
void covProbeSolveTolPop(const double *sav);
struct ImpCovSolveTolScope {
  double sav[2];
  ImpCovSolveTolScope() { covProbeSolveTolPush(sav); }
  ~ImpCovSolveTolScope() { covProbeSolveTolPop(sav); }
};

// Per-pseudo-subject proposal at the converged estimate: mode, information H,
// lower Cholesky of gamma*H^-1.  For a mixture these run over the EXPANDED
// subjects (component j is i + j*nsub), exactly as the E-step does: the objective
// differenced is the mixture MARGINAL, so every component needs its own proposal.
// Component 0 alone left the proportions' directions exactly flat and their SEs at 0.
static void impCovProposals(ImpEProp& P, int nExp, int neta, const arma::vec& gammaVec) {
  arma::vec mode(neta);
  arma::mat H(neta, neta, arma::fill::zeros);
  for (int id = 0; id < nExp; ++id) {
    impGetMode(id, mode);
    P.modes[id] = mode;
    if (!impGetHessian(id, H)) continue;
    H = 0.5 * (H + H.t()); // guard against tiny numerical asymmetry
    arma::mat Sigma; double ldv, lds;
    if (!(arma::inv_sympd(Sigma, H) && arma::log_det(ldv, lds, H) && lds > 0)) continue;
    Sigma *= gammaVec[id];
    arma::mat L;
    if (arma::chol(L, Sigma, "lower")) {
      P.Hs[id] = H; P.cholL[id] = L; P.logDetH[id] = ldv; P.haveL[id] = 1;
    }
  }
}

// The covariance step's fixed N(0,I) points for one subject
static arma::mat impCovQrZ(const ImpQrSet& q, int isample, int neta, int id, int covIter) {
  if (q.scramble != impQrScrambleNone) {
    // A fixed sample set with common random numbers across the FD
    // perturbations is what this step needs, and any valid randomization
    // gives one.  Use the LAST iteration the EM actually ran (covIter) so
    // it is the converged fit's randomization rather than its first
    // iteration's; under qrRefresh=FALSE the key is iteration-independent
    // anyway and this is the same set the whole fit used.
    return impQrZScrambled(q.bits, isample, neta, q.scramble,
                           (uint32_t)impBaseSeed() + 0x9E3779B1u, id, covIter);
  }
  if (!q.shift) return q.Z0;
  arma::vec sh(neta);
  for (int j = 0; j < neta; ++j) sh[j] = rxUnifEng(0.0, 1.0);
  return impQrZ(q.U0, &sh);
}

// Draw one fixed sample set per subject (serial, deterministic seed).  The
// engine is already allocated + seeded by rxWithSeed() at the fit start.
// With qr=TRUE the fixed set is the Sobol point set (per-subject shifted
// when qrShift; qrRefresh is irrelevant here -- the set never refreshes).
static void impCovDraw(std::vector<arma::mat>& Ss, const ImpEProp& P,
                       const std::vector<impProp>& props, int nExp, int isample, int neta,
                       int covIter) {
  ImpQrSet q = impQrSetup(isample, neta);
  uint32_t seed0 = (uint32_t)impBaseSeed() + 0x2545F491u;
  setRxThreadId(0);
  for (int id = 0; id < nExp; ++id) {
    if (!P.haveL[id]) continue;
    setSeedEng1(seed0 + (uint32_t)(id * 2 + 1));
    arma::mat S(isample, neta);
    if (q.on) {
      arma::mat Z = impCovQrZ(q, isample, neta, id, covIter);
      for (int k = 0; k < isample; ++k) {
        arma::vec z = Z.row(k).t();
        S.row(k) = (P.modes[id] + P.cholL[id] * (z * impPropScale(props[id], z))).t();
      }
    } else {
      impRngSamples(S, P.modes[id], P.cholL[id], props[id], isample, neta);
    }
    Ss[id] = S;
  }
  setRxThreadId(-1);
}

// Importance-sampling -2LL at a parameter vector, reusing the fixed samples.
// Each subject's contribution is an independent scalar; accumulate them into a
// per-subject buffer under the parallel loop, then reduce in id order so the sum
// is bit-identical to the serial `obj += ...` accumulation.  objBuf is allocated
// once and refilled per call (the FD Hessian calls it O(np^2) times).
struct ImpCovObj {
  const ImpEProp& P;
  const std::vector<arma::mat>& Ss;
  const std::vector<impProp>& props;
  const arma::vec& gammaVec;
  const std::vector<int>& pl;
  int nsub, neta, isample, Nm, nExp, ntheta, np, cores;
  bool doParCov;
  arma::vec par0;
  std::vector<double> objBuf;
  std::vector<char> objGood;
  std::vector<double> mixLl;
  // Each fixed sample's data part of the -log joint (impEvalJointLik() less its eta
  // prior), kept from the evaluation at the estimates: a point that moves only Omega
  // parameters changes only the prior, so it needs no solve.
  bool omReuse = false;
  std::vector<arma::vec> dataPart;
  bool haveData = false;
  int nReused = 0;
  // Progress bar over the finite-difference covariance evaluations, like the
  // focei covariance step.  The objective is called f0 (1) + 2*np (diagonal) +
  // 2*np*(np-1) (off-diagonal) = 1 + 2*np*np times; tick once per call.
  bool covProg = false;
  clock_t covT0 = 0;
  int covTot = 0, covCur = 0, covTick = 0;

  // Perturb parameter j: a theta on the subjects' parameter pointers, or an
  // Omega free parameter (which also rebuilds omegaInv / logDetOmegaInv5).
  void setPar(int j, double val) const {
    if (pl[j] < ntheta) impSetThetaAll(pl[j], val);
    else impSetOmegaThetaAll(pl[j] - ntheta, val);
  }

  bool thetasAtCentre(const arma::vec& par) const {
    for (int j = 0; j < np; ++j) {
      if (pl[j] < ntheta && par[j] != par0[j]) return false;
    }
    return true;
  }

  // -log joint at sample k of subject id
  double negLogJoint(int id, int k, const arma::vec& eta, bool reuse, bool record) {
    if (reuse) return dataPart[id][k] + impEtaPriorHalf(eta);
    double v = impEvalJointLik(eta, id);
    if (record) dataPart[id][k] = v - impEtaPriorHalf(eta);
    return v;
  }

  void subject(int id, double negHalfLogDetOmega, bool reuse, bool record) {
    if (!P.haveL[id]) return;
    if (!reuse) impForceResolve(id);
    if (record) dataPart[id].set_size(isample);
    // This subject's converged proposal scale (all equal under "global").
    // Must use the SAME proposal FAMILY and parameters as the E-step or the
    // reweighted objective is not the one the fit converged on.
    const impProp& pr = props[id];
    arma::vec q(isample); int nGood = 0;
    for (int k = 0; k < isample; ++k) {
      arma::vec eta = Ss[id].row(k).t();
      arma::vec d = eta - P.modes[id];
      double qk = -negLogJoint(id, k, eta, reuse, record) +
        impPropLogKernelRecip(pr, arma::as_scalar(d.t() * P.Hs[id] * d), gammaVec[id], neta);
      if (R_finite(qk)) { q[k] = qk; ++nGood; } else q[k] = R_NegInf;
    }
    if (nGood == 0) return;
    double qmax = q.max();
    double sumw = arma::accu(arma::exp(q - qmax));
    double logMeanExp = qmax + std::log(sumw / (double)isample);
    double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gammaVec[id]) - 0.5 * P.logDetH[id]
      + pr.corr;
    objBuf[id] = logMeanExp + Ci;   // log-likelihood; combined below
    objGood[id] = 1;
  }

  // Reduce in id order so the sum is bit-identical to a serial accumulation.
  // For a mixture each physical subject's contribution is the marginal
  // log(sum_m p_m L_im) over its components, which is what makes the mixture
  // proportions enter the objective at all.  A subject with no usable
  // component contributes 0, as the single-component path always did.
  double reduce() {
    double obj = 0.0;
    if (Nm == 1) {
      for (int id = 0; id < nsub; ++id) obj += objGood[id] ? -2.0 * objBuf[id] : 0.0;
      return obj;
    }
    for (int i = 0; i < nsub; ++i) {
      for (int m = 0; m < Nm; ++m) {
        mixLl[(size_t)m] = objGood[i + m*nsub] ? objBuf[i + m*nsub] : NA_REAL;
      }
      double lse = impMixLogSumExp(mixLl);
      if (R_FINITE(lse)) obj += -2.0 * lse;
    }
    return obj;
  }

  void init() {
    objBuf.assign(nExp, 0.0);
    objGood.assign(nExp, 0);
    mixLl.resize((size_t)Nm);
    par0.set_size(np);
    for (int j = 0; j < np; ++j)
      par0[j] = (pl[j] < ntheta) ? impGetFullThetaVal(pl[j])
                                 : impGetOmegaThetaVal(pl[j] - ntheta);
    omReuse = impCovReuseOn() && Nm == 1 && !impIsFo();
    dataPart.resize(nExp);
    covProg = impCovProgress();
    covT0 = clock();
    covTot = 1 + 2 * np * np;
  }

  double operator()(const arma::vec& par) {
    for (int j = 0; j < np; ++j) setPar(j, par[j]);
    const bool reuse = omReuse && haveData && thetasAtCentre(par);
    const bool record = omReuse && !haveData;
    if (reuse) ++nReused;
    // Re-read after setting: an Omega perturbation changes -0.5 log|Omega|.
    double negHalfLogDetOmega = impLogDetOmegaInv5();
    std::fill(objBuf.begin(), objBuf.end(), 0.0);
    std::fill(objGood.begin(), objGood.end(), 0);
    nmForEachSubject(rx, nExp, cores, doParCov, [&](int id) {
      subject(id, negHalfLogDetOmega, reuse, record);
    });
    double obj = reduce();
    if (covProg) covTick = par_progress(covCur++, covTot, covTick, 1, covT0, 0);
    if (record) haveData = true;
    return obj;
  }
};

// Publish the covariance: the theta block on e, then the full matrix through
// .impCovInstall().
static void impCovPublish(Environment& e, const arma::mat& cov, const arma::mat& info,
                          const std::vector<int>& pl, int ntheta, const arma::vec& par0) {
  int np = (int)pl.size();
  int nTh = 0;
  for (int j = 0; j < np; ++j) if (pl[j] < ntheta) ++nTh;
  arma::vec se(np);
  for (int j = 0; j < np; ++j) se[j] = (cov(j, j) > 0) ? std::sqrt(cov(j, j)) : NA_REAL;
  // Full covariance in the free-parameter (estimation) order: the thetas, then the
  // Omega parameters, which are entries of chol(Omega^-1).
  e["impCovInternal"] = wrap(cov);
  e["impCovThetaN"] = nTh;
  IntegerVector thIdxR(nTh);
  { int t = 0; for (int j = 0; j < np; ++j) if (pl[j] < ntheta) thIdxR[t++] = pl[j] + 1; }
  e["impCovThetaIdx"] = thIdxR;
  if (nTh > 0) {
    e["impCovTheta"] = wrap(arma::mat(cov.submat(0, 0, nTh - 1, nTh - 1)));
    e["impSeTheta"] = wrap(arma::vec(se.subvec(0, nTh - 1)));
  }
  // .impCovInstall() (R/impmap.R) maps the Omega rows to the variances and
  // covariances they are reported as, checks the result and publishes it as the
  // fit's covariance, before foceiFinalizeTables builds the SE / CI / correlation
  // tables from it.  It needs d(Omega)/d(p) of each Omega parameter in the
  // covariance, the parameter values (to say which parameterization that is),
  // and the information, to repair one that is not positive definite.
  List dOm = impOmegaParDeriv();
  List dOmCov(np - nTh);
  NumericVector omPar(np - nTh);
  for (int j = nTh; j < np; ++j) {
    dOmCov[j - nTh] = dOm[pl[j] - ntheta];
    omPar[j - nTh] = par0[j];
  }
  arma::mat Om;
  impGetOmega(Om);
  Environment nlmixr2 = Environment::namespace_env("nlmixr2est");
  Function covInstall = nlmixr2[".impCovInstall"];
  // held first: Function's argument list allocates while it is built
  RObject covR = wrap(cov), omR = wrap(Om), infoR = wrap(info);
  covInstall(e, covR, thIdxR, dOmCov, omR, omPar, infoR);
}

static void impComputeCov(Environment e, const arma::vec& gammaVec,
                          const std::vector<impProp>& props, int covIter) {
  // The proposal Hessians below and every objective evaluation re-score each
  // subject away from its mode, at perturbed parameters; the fit reports the
  // log-likelihoods the final MAP pass left at the estimates.
  CovLlikObsGuard _llikObsGuard;
  int nsub = impNsub();
  int neta = impNeta();
  int isample = impNsample();
  int Nm = impNmix();
  int nExp = nsub * Nm;   // expanded pseudo-subjects: component j is i + j*nsub
  // Parallelize the reweighted-objective subject loop (each subject contributes an
  // independent scalar to the -2LL) when the inner solves are thread-safe: no
  // mixture (expanded pseudo-subjects share solve rows) and no multi-endpoint
  // pool-sizing (neqOverride against the larger theta-sens pool) -- the same
  // contract as the E-step.  The one-time fixed-sample draw stays serial.
  int cores = std::max(1, impCores());
  bool doParCov = (cores > 1) && (impNmix() == 1) && !impPoolSizing();

  // Free parameters in the optimizer's (fixedTrans) order -- the same order the
  // fit's covariance uses.  pl[j] is a theta fullTheta index when < ntheta, else
  // the Omega parameter (pl[j] - ntheta).
  std::vector<int> pl;
  impGetCovParList(pl);
  int np = (int)pl.size();
  if (np == 0) return;
  int ntheta = impNtheta();

  ImpEProp P(nExp);
  impCovProposals(P, nExp, neta, gammaVec);
  std::vector<arma::mat> Ss(nExp);
  impCovDraw(Ss, P, props, nExp, isample, neta, covIter);

  ImpCovObj f{P, Ss, props, gammaVec, pl, nsub, neta, isample, Nm, nExp, ntheta, np, cores,
              doParCov};
  f.init();
  bool covInPlace = f.covProg && nmProgressInPlace();
  if (f.covProg) RSprintf("calculating covariance matrix\n");

  arma::mat Hess = impFdHessian(f.par0, std::ref(f));
  // objective evaluations that needed no solve
  e["impCovReused"] = f.nReused;
  if (f.covProg) nmProgressEnd(f.covTot, f.covTick, f.covT0, covInPlace);
  // Restore the converged estimates.
  for (int j = 0; j < np; ++j) f.setPar(j, f.par0[j]);
  // every pseudo-subject, not just component 0: the perturbed solves above ran
  // over nExp, so stopping at nsub leaves the other components' cached solves
  // sitting at the LAST perturbation rather than the converged estimates
  for (int id = 0; id < nExp; ++id) impForceResolve(id);

  // Observed information = 0.5 * Hess(-2LL) (symmetrized); covariance = inverse.
  arma::mat info = 0.25 * (Hess + Hess.t());
  arma::mat cov;
  if (!arma::inv_sympd(cov, info) && !arma::pinv(cov, info)) {
    cov = arma::mat(np, np, arma::fill::zeros); cov.fill(NA_REAL);
  }
  impCovPublish(e, cov, info, pl, ntheta, f.par0);
}

// State of one importance-sampling EM run (impOuter).
struct ImpOuterState {
  explicit ImpOuterState(Environment& env) : e(env) {}
  Environment& e;
  int nIter = 0, nsub = 0, neta = 0, isample = 0, cores = 1;
  std::string diagXform;
  double gamma = 1.0;
  // Convergence controller + proposal-scale adaptation (NONMEM ISCALE/IACCEPT/CTYPE).
  double iaccept = 0.0, iscaleMin = 0.0, iscaleMax = 0.0, ctol = 0.0;
  int nConvWindow = 0, mapIter = 0, nBurn = 0, nIterTotal = 0;
  bool burnFreezeOmega = false, eOnly = false, frozen = false;

  arma::mat condMean;
  std::vector<arma::mat> condVar;
  std::vector<arma::mat> sampS;
  std::vector<arma::vec> sampZk;
  arma::vec Li, Neff, Xi, NeffExp, XiExp, KhatExp;
  arma::mat aMat;                 // posterior mixture responsibilities (nsub x Nmix)
  int Nmix = 1;
  int nExp = 0;                   // expanded pseudo-subjects for the mixture E/M-step
  double obj = R_PosInf;
  int nSens = 0;
  // est="imp": skip the per-iteration MAP search; the E-step proposal is centered
  // at the running conditional mean with covariance gamma*Omega.
  bool isImp = false;
  // Per-expanded-subject proposal scale.  Under the current (global) rule every
  // element tracks the scalar `gamma`, so the E-step arithmetic is unchanged;
  // holding it as a vector is what lets a per-subject controller drive the
  // elements apart without touching the E-step again.
  arma::vec gammaVec;
  // Per-expanded-subject sample count.  Uniform unless AUTO raises it for a
  // badly covered subject; holding it as a vector is what lets that happen
  // without charging every other subject for the extra samples.
  arma::ivec isampleVec;
  // Counts actually USED by the last E-step.  isampleVec is reallocated after
  // the E-step, so reporting it directly pairs next iteration's counts with
  // this iteration's Neff and can show an effective-sample FRACTION above 1,
  // which is impossible.
  arma::ivec isampleUsed;
  arma::vec dfUsed;               // df actually used last E-step
  // Same for gamma.  The controllers all mutate their vectors AFTER the E-step,
  // for the next iteration, so anything consuming "the converged proposal" must
  // read the snapshot rather than the live vector -- including impComputeCov.
  arma::vec gammaUsed;
  double gammaScalarUsed = 1.0;   // scalar counterpart of gammaUsed
  // Per-expanded-subject proposal df and acceptance target.  Uniform unless
  // AUTO differentiates them.
  arma::vec dfVec;
  // Proposal FAMILY, per expanded subject, mirroring dfVec/dfFloor.  Kept as a
  // separate vector rather than folding df into a struct: AUTO mutates dfVec in
  // place at five sites and its tuned comparisons are arithmetic on it.
  impPropSpec propSpec;
  arma::ivec typeVec, typeFloor;
  std::vector<impProp> props, propsUsed;
  // Improvability state for the AUTO df escalation (see AUTO step 3).  noImp
  // counts consecutive iterations in which the current rung has failed to
  // improve on what the subject manages without any escalation.
  // Lowest df a subject may be returned to.  Withdrawal undoes the k-hat
  // ESCALATION, never the df that step 1 assigned: a non-normal endpoint needs a
  // heavy-tailed proposal by construction, so dropping it to a Gaussian because
  // k-hat did not improve would break the model, not repair it.
  arma::vec dfFloor;
  // What this subject's k-hat looks like WITHOUT any escalation -- the
  // counterfactual withdrawal has to judge a rung against.  Refreshed on every
  // iteration the subject spends at its floor, and frozen while it is escalated.
  arma::vec kAtFloor;
  arma::ivec noImp;
  arma::ivec escDead;             // 1 = escalation withdrawn, do not retry
  arma::vec iacceptVec;
  bool autoOn = false;
  int isampleBudget = 0;          // AUTO reallocates within this total
  bool nonmemSparse = false;
  int dfPatience = 0;             // <= 0 disables withdrawal entirely
  // gammaMethod="individual": each subject's gamma_i is adapted two-sided so
  // that its own xi_i approaches iaccept.  Otherwise the legacy global rule
  // drives one shared scale off the MEAN Kish effective-sample fraction.
  //
  // NONMEM specifies the OBJECTIVE of this loop but not its update law.  The
  // NM7 Technical Guide gives the target ("gamma is continually adjusted so
  // that xi_i approximates IACCEPT", note after eq. 1.76), the per-subject
  // granularity (gamma_i in eq. 1.90; "adjusted for each subject" at 1.474),
  // and the [ISCALE_MIN, ISCALE_MAX] bounds -- and nothing further.  The
  // functional form below (multiplicative, exponent 2/neta, symmetric cap) is
  // ours; the guide publishes no formula to copy.
  bool gammaInd = false;
  // Largest relative per-step change the controller applied last iteration;
  // starts large so convergence cannot trip before the controller has run.
  double gammaStepPrev = R_PosInf;
  int savedCores = -1;

  uint32_t qrPinSeed = 0;
  bool sir = false;
  int sirN = 0;
  uint32_t sirSeed = 0;
  bool harvestSens = false;
  std::vector<impThetaSensData> outSens;

  arma::mat Om0, omMask;
  std::vector<int> omFloored;     // etas the installed omega has at the 1e-6 floor
  // Eta indices whose Omega diagonal is fix()ed: their rows/columns are held at
  // the starting Omega through every EM update.
  std::vector<int> omFixedEta;

  std::vector<double> objTrace, gammaTrace, neffTrace, xiTrace;
  std::vector<arma::vec> parHist;
  bool converged = false;
  int iterRun = 0;
  // this iteration's mean effective-sample fraction and mean xi
  double accFrac = 0.0, xiMean = NA_REAL;
};

// Drift tolerance on the gamma window mean for gammaRule="target".  Looser than
// the floor rule's 1e-3 because that rule's gamma is exactly constant once
// coverage is healthy, while a target-tracking gamma carries the Monte-Carlo
// noise of xi.  PROVISIONAL -- a constant tuned for this rule, to be re-measured
// with the rest of the target-rule set.
static const double impGammaTargetTol = 1e-2;

static void impOuterReadControls(ImpOuterState& st) {
  st.nIter = impNiter();
  st.diagXform = impDiagXform();
  st.nsub = impNsub();
  st.neta = impNeta();
  st.isample = impNsample();
  st.gamma = impGammaProp();
  st.cores = std::max(1, impCores());
  st.iaccept = impIaccept();
  st.iscaleMin = impIscaleMin();
  st.iscaleMax = impIscaleMax();
  st.nConvWindow = impNconvWindow();
  st.ctol = impCtol();
  st.mapIter = std::max(0, impMapIter());
  // Burn-in iterations run BEFORE the nIter budget rather than out of it, so
  // nBurn never silently shortens the fit.  They exist to let the gamma / AUTO
  // controllers settle, so those run normally; only the Omega M-step (when
  // burnFreezeOmega) and the convergence test are held back.
  st.nBurn = std::max(0, impNburn());
  st.burnFreezeOmega = impBurnFreezeOmega();
  // nIter = 0 is an E-step-only evaluation at the supplied parameters (NONMEM
  // EONLY=1): one E-step, no M-step, and no burn-in (it would move the parameters).
  st.eOnly = (st.nIter <= 0);
  // A frozen run (the post-fit "imp" covariance at a fit's estimates) keeps the
  // parameters where they were supplied for all nIter iterations: E-steps only,
  // with the proposal controllers adapting between them, no M-step, no burn-in
  // and no convergence test.
  st.frozen = impFrozen();
  if (st.eOnly || st.frozen) st.nBurn = 0;
  st.nIterTotal = st.eOnly ? 1 : st.nBurn + st.nIter;
  st.Nmix = impNmix();
  st.nExp = st.nsub * st.Nmix;
  st.nSens = impThetaSensN();
  st.isImp = impIsImp();
}

static void impOuterInitVectors(ImpOuterState& st) {
  const int nExp = st.nExp;
  st.gammaVec.set_size(nExp); st.gammaVec.fill(st.gamma);
  st.isampleVec.set_size(nExp); st.isampleVec.fill(st.isample);
  st.isampleUsed.set_size(nExp); st.isampleUsed.fill(st.isample);
  st.dfUsed.set_size(nExp); st.dfUsed.fill(impDf());
  st.gammaUsed.set_size(nExp); st.gammaUsed.fill(st.gamma);
  st.gammaScalarUsed = st.gamma;
  st.dfVec.set_size(nExp); st.dfVec.fill(impDf());
  st.propSpec = impPropSpecFromControl(impDf());
  st.typeVec.set_size(nExp); st.typeVec.fill(st.propSpec.type);
  st.typeFloor.set_size(nExp); st.typeFloor.fill(st.propSpec.type);
  st.props.resize(nExp);
  st.propsUsed.resize(nExp);
  st.dfFloor.set_size(nExp); st.dfFloor.fill(impDf());
  st.kAtFloor.set_size(nExp); st.kAtFloor.fill(NA_REAL);
  st.noImp.zeros(nExp);
  st.escDead.zeros(nExp);
  st.iacceptVec.set_size(nExp); st.iacceptVec.fill(st.iaccept);
  st.autoOn = impAutoEnabled();
  st.isampleBudget = st.isample * nExp;
}

// A per-subject isample= vector is given per BASE subject; tile it across
// the mixture's expanded pseudo-subjects so component j of subject i gets
// subject i's count.
static void impOuterTileIsample(ImpOuterState& st) {
  std::vector<int> isv;
  impNsampleVecGet(isv);
  if ((int)isv.size() > 1 && (int)isv.size() != st.nsub) {
    Rf_error("'isample' has length %d but there are %d subjects; a per-subject "
             "isample vector must have one entry per subject",
             (int)isv.size(), st.nsub);
  }
  if ((int)isv.size() != st.nsub) return;
  for (int i = 0; i < st.nsub; ++i) {
    int v = isv[i] < 1 ? 1 : isv[i];
    for (int j = 0; j < st.Nmix; ++j) st.isampleVec[i + j * st.nsub] = v;
  }
}

// ---- AUTO step 1: proposal df from each subject's data sparsity ----------
// Decided once, up front, from nobs/neta: a subject with few observations
// relative to the number of random effects has a poorly-determined, often
// heavy-tailed conditional density that a normal proposal cannot cover, and
// no sample count repairs that (measured: boosting a failing subject 10x
// moved its k-hat 0.76 -> 3.28, because extra draws reach further into the
// tail the proposal misses).  Heavier tails do repair it, cheaply.
//
// The trigger is the tutorial's, verbatim: "When there are fewer data points
// than there are ETAs to be estimated (sparse data) or data are categorical,
// then DF should be set to a nonzero number to use the t-distribution sampler
// and/or decrease the acceptance rate IACCEPT from the default value of 0.4
// to about 0.2" (Bauer, NONMEM Tutorial Part II, 2019).  So the test is
// nobs < neta -- NOT a wider sparsity band -- and categorical/non-normal data
// qualifies regardless of how much of it there is.
//
// On the VALUES: IACCEPT ~ 0.2 is the tutorial's own number.  DF = 4 is NOT
// -- the tutorial says only "DF should be set to a nonzero number".  NONMEM
// does not publish what AUTO=1 actually picks, so 4 is nlmixr2's choice,
// informed by the measurement that even df = 20 cleared every failing
// subject on theophylline (max k-hat 1.61 -> -0.58) for 0.25% of the
// effective sample size, so a heavier default is safe and cheap.  Any
// comparison against NONMEM should treat these numbers as nlmixr2's, not as
// a reproduction of NONMEM's internals.
// The nobs<neta half of that trigger is GATED by default, and this is a
// deliberate divergence from the tutorial.  Measured on a fixture built to the
// tutorial's own definition of sparse (2 observations, 3 etas), 8 seeds vs an
// isample=8000 reference, applying it makes every number worse:
//
//                     objRMSE   OmegaRMSE   max k-hat   subjects k>0.7
//   no adaptation      0.1517     0.02379      1.065         3.88
//   global df 30       0.1666     0.02695      1.075         5.00
//   AUTO (this rule)   0.2059     0.03163      1.088         4.88
//
// The reason is visible in the reference, which itself reads max k-hat 0.794:
// with fewer observations than random effects the individual posterior is not
// identified, so the heavy tail is STRUCTURAL and no proposal shape repairs
// it.  Escalating there spends accuracy on a problem the proposal does not
// own.  So sparsity alone no longer assigns a t proposal; k-hat has to ask for
// it, and the improvability check below withdraws it if it does not help.
//
// The categorical/non-normal half is NOT gated -- measured on a sparse Poisson
// fixture AUTO is the best of the three methods (objRMSE 0.00341 vs 0.00525
// for no adaptation), even though its k-hat never fails.
//
// impmapControl(autoNonmemSparse=TRUE) restores the unconditional tutorial
// behavior.  It is a real option, not a courtesy: NONMEM's own testing is not
// published, the models it was tuned on are not ours, and someone who measures
// the reverse on their own problem should be able to have the documented rule.
static void impAutoStep1(ImpOuterState& st) {
  bool nonNormal = impAutoNonNormal();
  // The df ladder is defined only for the MVN/MVT axis -- it encodes
  // "Gaussian" as df <= 0, and its constants (30/20, the k-hat 0.7/1.0
  // thresholds, autoDfPatience) were all swept against a Gaussian
  // baseline.  A subject already on an exponential-tailed MVL or on a
  // mixture has a different baseline and no measurement behind it, so
  // AUTO leaves its family alone.  The sample-budget reallocation below
  // is family-agnostic and still applies to every subject.
  const bool dfAxis = (st.propSpec.type == impPropMvn || st.propSpec.type == impPropMvt);
  for (int i = 0; i < st.nsub; ++i) {
    bool sparse = st.nonmemSparse && (impNobs(i) < st.neta);
    // TUNED: start at a MILD t (df 20) rather than the heaviest.  df 20 was
    // measured to clear every failing subject on theophylline for 0.25% of
    // the effective sample size, whereas df 4 costs far more Monte-Carlo
    // noise for no extra tail benefit.  The k-hat escalation below drops it
    // further only on evidence.
    // TUNED (df sweep, theophylline, 12 seeds, RMSE vs an isample=20000
    // reference): df 30 clears the tail failure completely (max k-hat
    // 2.58 -> -0.02, bad subjects 2.08 -> 0.08) for a 7% RMSE cost
    // (0.0240 -> 0.0257).  Entering at 20 costs 33% and at 8 costs 3x for no
    // extra tail benefit, so the earlier ladder was simply mistuned.
    // impDf() rather than 0 for the untriggered case: a global df= is an
    // explicit request and must survive auto=TRUE.  At the default df = 0
    // this is identical to the old expression.
    double dfI = (dfAxis && (sparse || nonNormal)) ? 30.0 : impDf();
    // TUNED: iaccept is left alone here.  Lowering it to 0.2 forces gamma
    // wide, and widening a Gaussian is the lever measured NOT to fix tails
    // while costing a lot of ESS -- on a Poisson fixture whose k-hat was
    // already -1.33 (no failure at all) it cut ESS 0.549 -> 0.412 and
    // doubled the objective noise for nothing.  It is now applied only where
    // k-hat says the proposal is genuinely struggling (below).
    double iaI = st.iaccept;
    int typeI = st.propSpec.type;
    if (dfAxis) typeI = dfI > 0.0 ? impPropMvt : impPropMvn;
    for (int j = 0; j < st.Nmix; ++j) {
      st.dfVec[i + j * st.nsub] = dfI;
      st.dfFloor[i + j * st.nsub] = dfI;   // withdrawal may not go below this
      st.typeVec[i + j * st.nsub] = typeI;
      st.typeFloor[i + j * st.nsub] = typeI;
      st.iacceptVec[i + j * st.nsub] = iaI;
    }
  }
}

// Seeds, SIR and the #958 sensitivity harvest
static void impOuterSetupSampling(ImpOuterState& st) {
  // Fit-constant salted base for the pinned (qrRefresh=FALSE) Cranley-Patterson
  // shift streams; only consumed in that mode so every other path's draw
  // stream is unchanged.
  st.qrPinSeed = (uint32_t)impBaseSeed() + 0x9E3779B1u;
  // Salted base for the SIR stratified offsets; SIR only engages where it
  // actually shrinks that subject's sample (sirN < isample).
  st.sir = impSirEnabled();
  st.sirN = impSirN();
  st.sirSeed = (uint32_t)impBaseSeed() + 0x7F4A7C15u;
  // #958 combined build: harvest the M-step's theta-sensitivity columns
  // straight off the E-step's own per-sample inner solve (impEStep) instead
  // of the M-step re-solving every sample -- see impHarvestSample.
  // Off under SIR: the M-step there only needs a resampled SUBSET of the
  // E-step's samples (sirN < isample), so eagerly harvesting every sample in
  // the E-step would do MORE solving, not less; that case keeps the existing
  // (unfused) M-step path unchanged.
  //
  // The E-step loop is parallel over subjects, and a subject
  // whose inner solve fell back to doFD/pred solves a DIFFERENT registered
  // slot than a subject still on the plain inner path -- concurrently, in
  // the same loop.  That used to be able to race on rx->ndiff (a field on
  // the single shared rx_solve struct); it no longer can, because
  // odeSwapSolveInd() now skips writing rx->ndiff entirely for a slot whose
  // own declared need is 0 (odeSlotPred never needs a Jacobian derivative --
  // "rxPred implements no sensitivities" -- so it never touches the shared
  // field), leaving only same-value writes from concurrent inner solves,
  // which is benign.  See odeSwap.cpp's odeSwapSolveInd for the fix and its
  // KNOWN GAP note.
  st.harvestSens = impCombSensEnabled() && (st.nSens > 0) && !st.sir;
}

// Omega structure mask: only the elements estimated in the model (nonzero in
// the starting Omega) are updated; the rest stay zero so the parameterization
// keeps the same number of Omega thetas.
static void impOuterSetupOmega(ImpOuterState& st) {
  impGetOmega(st.Om0);
  st.omMask = arma::conv_to<arma::mat>::from(st.Om0 != 0.0);
  impGetOmegaFixedEta(st.omFixedEta);
}

// Before the E-step: re-MAP, flatten the global gamma and build the proposals.
static void impOuterPreEStep(ImpOuterState& st, int iter) {
  // MAP-assist period: impReMap() is the mu-referenced FOCEI inner problem,
  // so re-centering every iteration is the dominant cost of est="impmap".
  // mapIter = 0 keeps the startup MAP and never re-centers.
  if (iter > 0 && !st.isImp && st.mapIter > 0 && (iter % st.mapIter) == 0) impReMap();
  // Global rule: every expanded subject shares the scalar scale.  Under
  // "individual" gammaVec already carries the per-subject scales the
  // controller set at the end of the previous iteration (and the constructor
  // fill for iteration 0), so it must NOT be flattened here.
  if (!st.gammaInd) st.gammaVec.fill(st.gamma);
  st.isampleUsed = st.isampleVec;
  st.dfUsed = st.dfVec;
  // Materialize this iteration's proposals from (family, df, spec).  Built
  // BEFORE the E-step so the parallel loop only ever reads them, and
  // snapshotted in the same breath as dfUsed because the controllers below
  // advance the live vectors for an iteration that never ran.
  for (int id = 0; id < st.nExp; ++id) {
    st.props[id] = impPropBuild(st.typeVec[id], st.dfVec[id], st.propSpec.nmix,
                                st.propSpec.c, st.propSpec.w, st.neta);
  }
  st.propsUsed = st.props;
  st.gammaUsed = st.gammaVec;
  st.gammaScalarUsed = st.gamma;
}

// The E-step's objective, mean effective-sample fraction and mean xi
static void impOuterSummarize(ImpOuterState& st) {
  st.obj = 0.0;
  for (int id = 0; id < st.nsub; ++id) if (R_finite(st.Li[id])) st.obj += 2.0 * st.Li[id];
  // Mean effective-sample fraction (importance-sampling "acceptance ratio").
  double accFrac = 0.0; int nAcc = 0;
  for (int id = 0; id < st.nsub; ++id)
    if (R_finite(st.Neff[id])) { accFrac += st.Neff[id] / (double)st.isampleUsed[id]; ++nAcc; }
  if (nAcc > 0) accFrac /= (double)nAcc;
  st.accFrac = accFrac;
  // Mean NONMEM-style xi (a DIFFERENT statistic from accFrac -- see impEStep).
  double xiMean = 0.0; int nXi = 0;
  for (int id = 0; id < st.nsub; ++id)
    if (R_finite(st.Xi[id])) { xiMean += st.Xi[id]; ++nXi; }
  if (nXi > 0) xiMean /= (double)nXi; else xiMean = NA_REAL;
  st.xiMean = xiMean;
  // `gamma` is the scalar summary of the scale actually used this iteration:
  // itself under "global", the mean of the per-subject scales under
  // "individual".  It is what gets traced and reported as impGammaUsed.
  if (st.gammaInd) st.gamma = arma::mean(st.gammaVec);
}

// Withdrawal.  Only for a subject actually moved off its floor -- testing
// dfVec > 0 instead stranded subjects that were never escalated at all --
// and only once the ladder is EXHAUSTED (!wantEsc), so "this is not
// helping" means the heaviest applicable rung is not helping.  Note
// "escalated" is dfVec != dfFloor, not an inequality: escalating a
// Gaussian floor RAISES df (0 -> 30) while escalating a t floor LOWERS it
// (30 -> 20).  Returns true when the escalation was withdrawn.
static bool impAutoWithdraw(ImpOuterState& st, int id, double kh, bool wantEsc) {
  if (st.nonmemSparse || st.dfPatience <= 0 || wantEsc ||
      st.dfVec[id] == st.dfFloor[id] || !R_finite(st.kAtFloor[id])) {
    return false;
  }
  // Non-strict on the margin: a rung that buys exactly the 0.1 it is
  // asked for has earned its place.  Strict `<` withdrew a proposal that
  // took k-hat 1.25 -> 0.75 against a floor of 0.85.
  bool improved = (kh <= 0.7) || (kh <= st.kAtFloor[id] - 0.1);
  if (improved) {
    st.noImp[id] = 0;
    return false;
  }
  if (++st.noImp[id] < st.dfPatience) return false;
  // back to the step-1 floor, NOT to 0: for a non-normal endpoint that
  // floor is the t proposal the model requires, and withdrawal is
  // permanent (escDead), so dropping below it could never be undone.
  st.dfVec[id] = st.dfFloor[id];
  st.typeVec[id] = st.typeFloor[id];
  st.escDead[id] = 1;
  return true;
}

// AUTO step 3 for one expanded subject
static void impAutoEscalateOne(ImpOuterState& st, int id) {
  // the ladder walks the MVN/MVT axis only -- see AUTO step 1 above
  if (st.typeVec[id] != impPropMvn && st.typeVec[id] != impPropMvt) return;
  double kh = st.KhatExp[id];
  if (!R_finite(kh)) return;              // no usable k-hat: leave alone
  if (st.escDead[id]) return;             // measured not to help here
  // IMPROVABILITY.  A t proposal repairs a tail the proposal shape is
  // missing; it cannot repair one the DATA creates.  With fewer
  // observations than random effects the individual posterior is not
  // identified and k-hat stays high however heavy the proposal gets
  // (measured: 1.065 under a Gaussian, 1.075 at df 30, 1.088 under AUTO --
  // and the isample=8000 reference itself reads 0.794).  Escalating there
  // costs accuracy for nothing, so withdraw it once the evidence is in.
  //
  // This is NOT the circular relaxation the note below rejects.  That one
  // backs off because k-hat became GOOD, which is evidence the remedy is
  // working.  This backs off because k-hat stayed BAD after escalating,
  // which is evidence the remedy does not apply.  Withdrawal is final so
  // the pair cannot oscillate.
  //
  // Both halves are load bearing, and gating the trigger alone is NOT
  // enough.  On the sparse fixture, objective RMSE by method: 0.206 with
  // the tutorial rule applied, 0.181 gating it but never withdrawing,
  // 0.152 with AUTO off entirely, 0.102 gating AND withdrawing.  Gating on
  // its own is still worse than not adapting at all; the withdrawal is
  // what turns it into a win.
  //
  // autoDfPatience TUNED to 2 (sweep over 0,1,2,3,5; 8 seeds):
  //
  //   patience   sparse obj   theo3 obj   sparse k>0.7   theo3 k>0.7
  //     0 (off)     0.18065      0.01654        5.12          0.25
  //     1           0.18740      0.01418        4.50          0.50
  //     2           0.10198      0.01588        4.62          0.38
  //     3           0.12873      0.01654        5.00          0.25
  //     5           0.18267      0.01654        5.12          0.25
  //
  // 2 is the sparse optimum, and the curve rises either side of it.  It is
  // not the best objective on theo3 (1 is), but 1 has the worst tail on
  // every fixture, and tail behaviour is weighted higher: infinite-variance
  // weights are a correctness problem with unbounded error while
  // Monte-Carlo noise is bounded.  Same reasoning that puts auto on at all.
  //
  // Measured BEFORE the monotone-baseline fix, 1 looked best on sparse
  // (0.097 against 0.142).  That was the easy-baseline bug withdrawing at
  // roughly the right moment by accident; with the bar correct the ordering
  // inverts.  A reminder that a constant tuned against unverified code is
  // only as good as the code under it.
  //
  // Patience also interacts with nIter: at 3 and 5 the counter cannot
  // accumulate inside a 12-iteration fit on a fixture whose k-hat does
  // improve, which is why those rows match switching withdrawal off.
  //
  // Decide the rung this k-hat asks for FIRST, because withdrawal must not
  // pre-empt an escalation that is still available.  Running the withdrawal
  // check first retired a subject at df 30 carrying its last strike the
  // moment its k-hat crossed 1.0 -- it went to the floor and escDead
  // without ever trying df 20, which is precisely the rung that k-hat was
  // asking for.
  double want = 0.0;
  bool wantEsc = false;
  if (kh > 0.7) {
    // lightest tail plausibly heavy enough for this severity
    want = (kh > 1.0) ? 20.0 : 30.0;
    // one-way: only ever go heavier, so the shape cannot oscillate
    wantEsc = (st.typeVec[id] == impPropMvn || want < st.dfVec[id]);
  }
  // The counterfactual, kept FRESH.  Three separate review findings all
  // reduce to getting this value wrong:
  //   * frozen at a pre-deterioration reading, a rung that genuinely
  //     repairs a since-degraded subject scores as failing and is
  //     withdrawn;
  //   * taken from the k-hat at the moment of escalation, a noise spike
  //     sets an easy bar and a subject stalls at a proposal achieving
  //     nothing;
  //   * taken as the running minimum, an old good value is pinned and a
  //     working rung is withdrawn again.
  // Recording it while the subject sits at its floor solves all three: it
  // tracks deterioration that happens BEFORE any escalation, and it stops
  // tracking the moment an escalation could be responsible for the value.
  if (st.dfVec[id] == st.dfFloor[id]) st.kAtFloor[id] = kh;
  if (impAutoWithdraw(st, id, kh, wantEsc)) return;
  if (wantEsc) {
    // Each rung gets a full patience window; inheriting strikes withdrew a
    // heavier rung after a single iteration for want of an improvement it
    // had not had time to deliver.
    st.noImp[id] = 0;
    st.dfVec[id] = want;
    st.typeVec[id] = impPropMvt;
  }
  // NOTE deliberately one-way on SUCCESS.  Relaxing on a low k-hat is
  // circular: once a t proposal is in place k-hat drops precisely BECAUSE
  // it is working (measured -1.3 to -1.9 on subjects that read 2.03 under a
  // Gaussian), so "safe now" is evidence the fix is needed, not that it is
  // not.  Relaxing on it walks straight back into the failure and
  // oscillates.  Efficiency comes from picking the right rung immediately
  // -- the severity mapping above -- not from backing off afterwards.
}

// ---- AUTO step 3: escalate df when k-hat reports tail failure ----------
// The tutorial's df trigger (nobs < neta, or categorical) is necessary but
// NOT sufficient: measured on plain theophylline -- 11 observations against
// 1 eta, transformably normal, so neither condition fires -- four of twelve
// subjects still have infinite-variance weights (max k-hat 2.36), which a
// t proposal clears completely.  NONMEM has no tail diagnostic to notice
// this; we do, so AUTO uses it.
//
// This is nlmixr2 going BEYOND the documented NONMEM rule, not reproducing
// it.  Escalation is one-way and stepwise (Gaussian -> 8 -> 4) so a noisy
// k-hat cannot oscillate the proposal shape from iteration to iteration.
//
// Severity-driven, and two-way.  The earlier version walked a fixed
// ladder one rung per iteration regardless of HOW bad k-hat was, and
// never came back -- so a subject with k-hat 0.75 took the same path as
// one at 3.0, and any subject that did not need a t proposal kept paying
// for it.  Both waste accuracy, which is what the gate was failing on.
//
// Mapping re-derived after the imp/impmap/focei pooling fixes, which removed
// the tail failures the ORIGINAL sweep was calibrated on -- see
// plans/imp-auto-reinstrument.md and the appendix in the nlmixr2 imp article.
// Sweep: theophylline + 3 etas, 8 seeds, RMSE vs an isample=8000 reference.
//
//   df   objRMSE   omegaRMSE   max k-hat   subjects k>0.7
//    0    0.0113     0.00406      0.941        2.38
//   30    0.0095     0.00281      0.593        0.25
//   20    0.0097     0.00255      0.484        0.12
//   12    0.0105     0.00224      0.405        0.00
//    8    0.0116     0.00223      0.273        0.00
//
// The objective RMSE is MINIMIZED at df 30 and degrades past it -- df 8 is
// worse than no t proposal at all -- while the tail and Omega keep improving.
// 30/20 is that trade-off's sweet spot; heavier rungs buy Omega at a real
// objective cost.
//
// The old k>2 -> df 12 rung is gone.  Post-fix nothing in the sweep exceeds
// k-hat ~1.1 (theophylline 3-eta peaks at 0.94, a deliberately sparse
// nobs<neta fixture at 1.07), so it was unreachable, and where df 12 could be
// measured it was worse than 20 on the objective.  An unreachable rung
// calibrated on a regime that no longer occurs is not a safety margin.
static void impAutoEscalate(ImpOuterState& st) {
  for (int id = 0; id < st.nExp; ++id) impAutoEscalateOne(st, id);
}

// ---- AUTO step 2: reallocate the sample budget across subjects --------
// The tutorial's trigger is "many ETAs or the objective function has large
// stochastic fluctuations".  Both show up as a low effective-sample
// FRACTION for that subject, so the budget is shifted toward the subjects
// whose weights are least effective and away from the ones that are
// already cheap to resolve.  The total is held at isample*nExp, so this is
// load-balancing rather than a cost increase.
//
// Deliberately NOT driven by k-hat: a bad tail is not fixable by more
// samples (measured: 10x on a failing subject moved k-hat 0.76 -> 3.28),
// that is what the df assignment above is for.  Sample count is the right
// lever only for ordinary Monte-Carlo noise.
static void impAutoRealloc(ImpOuterState& st) {
  const int nExp = st.nExp;
  arma::vec need(nExp, arma::fill::ones);
  for (int id = 0; id < nExp; ++id) {
    // NeffExp, not Neff: id runs over the EXPANDED subjects, and Neff has
    // one (responsibility-combined) entry per base subject.
    double fr = R_finite(st.NeffExp[id]) ? st.NeffExp[id] / (double)st.isampleUsed[id] : 1.0;
    if (!(fr > 0.0) || !R_finite(fr)) fr = 1.0;
    // Only deviate from uniform where the effective-sample fraction is
    // genuinely deficient.  Reallocating on small differences just adds
    // per-subject sampling noise for no gain; a subject already above the
    // acceptance target does not need more draws.
    // NONMEM's stated behaviour is two-sided: more budget to subjects with
    // noisy likelihoods or many active random effects, and LESS to
    // data-rich subjects whose conditional mode is easy to resolve.  A
    // threshold rule ("boost only below iaccept") captures only the first
    // half and is completely inert when every subject is above the target,
    // which is the common case -- so allocate PROPORTIONALLY to inverse
    // effective-sample fraction instead, normalized to the fixed budget.
    // Easy subjects then automatically fund the hard ones.
    double w = 1.0 / std::max(0.02, std::min(1.0, fr));
    need[id] = std::min(6.0, w);
  }
  double tot = arma::accu(need);
  if (!(tot > 0.0 && R_finite(tot))) return;
  for (int id = 0; id < nExp; ++id) {
    int v = (int)std::lround((double)st.isampleBudget * need[id] / tot);
    if (v < 25) v = 25;                             // floor: keep PSIS usable
    if (v > 20 * st.isample) v = 20 * st.isample;   // ceiling: bound the cost
    st.isampleVec[id] = v;
  }
}

// #958 fused build: outSens[eid] was already filled sample-by-sample
// during THIS iteration's E-step (impEStep), off the same inner solve
// impEvalJointLik used for the importance weights -- no M-step solve
// at all for a subject that harvested cleanly.  A subject every one of
// whose samples fell back to the pred-model FD path in the E-step
// (outSens[eid].nobs == 0, impEStep leaves it unharvested rather than
// re-solving from an OpenMP region -- see impHarvestSample) falls back here to
// the ordinary (serial-safe) solve-based read, so no subject's M-step
// contribution is silently dropped.
static void impMStepThetaHarvested(ImpOuterState& st, arma::vec& g, arma::mat& H) {
  for (int eid = 0; eid < st.nExp; ++eid) {
    if (st.sampS[eid].n_rows == 0) continue;
    if (st.outSens[eid].nobs > 0) {
      impThetaAccumOne(st.outSens[eid], st.sampZk[eid], g, H);
    } else {
      impThetaSensData c;
      impThetaSensCollect(eid, st.sampS[eid], c);
      impThetaAccumOne(c, st.sampZk[eid], g, H);
    }
  }
}

// One subject's theta-sensitivity solves (SIR-resampled if enabled); true when used
static bool impMStepCollectOne(ImpOuterState& st, int iter, int eid, impThetaSensData& coll,
                               arma::vec& zksub) {
  if (st.sampS[eid].n_rows == 0) return false;
  if (!(st.sir && st.sirN < (int)st.sampS[eid].n_rows)) {
    zksub = st.sampZk[eid];
    impThetaSensCollect(eid, st.sampS[eid], coll);
    return true;
  }
  // SIR acceleration: an equal-weight systematic resample stands in
  // for the full weighted sample, cutting the theta-sens solves from
  // isample to sirN per subject.  For a mixture, zk sums to the
  // component responsibility a_ij (folded in by impEStep), so the
  // equal weights carry a_ij/sirN to preserve the component mass.
  double aij = arma::accu(st.sampZk[eid]);
  if (!(aij > 0.0 && R_finite(aij))) return false;
  // Seed this subject's stream on the current thread; the drawn u0 is
  // a pure function of the (iter, eid) seed, so it is thread-count
  // invariant.
  nmSetSeedEng1(st.sirSeed + (uint32_t)((iter * st.nExp + eid) * 2));
  double u0 = rxUnifEng(0.0, 1.0);
  arma::uvec ridx = impSirIndex(st.sampZk[eid], st.sirN, u0);
  arma::mat Ssir = st.sampS[eid].rows(ridx);
  zksub.set_size(st.sirN); zksub.fill(aij / (double)st.sirN);
  impThetaSensCollect(eid, Ssir, coll);
  return true;
}

static void impMStepThetaCollect(ImpOuterState& st, int iter, arma::vec& g, arma::mat& H) {
  const int nExp = st.nExp;
  std::vector<impThetaSensData> coll(nExp);
  // Only zksub is retained per subject (the serial accumulation needs it); the
  // sample matrix is passed straight into impThetaSensCollect (sampS[eid] as-is,
  // or a local SIR resample), so there is no per-subject copy of the full sample.
  std::vector<arma::vec> zksub(nExp);
  std::vector<char> useSub(nExp, 0);
  // Parallelize the theta-sensitivity solves over subjects; each subject writes
  // only its own slot (coll[eid]) and the score/Hessian accumulation below stays
  // serial in eid order, so the result is bit-identical to the serial M-step.
  // Only engage when the solve method is thread-safe (liblsoda); the
  // inner-parallel flag (impInnerParallelOn/Off) defers worker-thread R-API
  // warnings.  impThetaSensCollect reads calc_lhs into a thread-local buffer and
  // loosens tolerance per-subject, so no shared solve state is touched.  cores
  // already carries the mixture / pool-sizing serial guard (forced to 1).
  bool doParM = (st.cores > 1) && impMStepParallelOk();
  if (doParM) impInnerParallelOn();
  nmForEachSubject(rx, nExp, st.cores, doParM, nmStatic, [&](int eid) {
    if (impMStepCollectOne(st, iter, eid, coll[eid], zksub[eid])) useSub[eid] = 1;
  });
  if (doParM) impInnerParallelOff();
  // Serial accumulation in eid order -- bit-identical to the serial theta score.
  for (int eid = 0; eid < nExp; ++eid) {
    if (useSub[eid]) impThetaAccumOne(coll[eid], zksub[eid], g, H);
  }
}

// Levenberg-Marquardt, escalated until the step lands inside a trust region
// relative to each theta's own magnitude.  H is Gauss-Newton and hence positive
// semi-definite, so the damping only improves its conditioning and rotates the
// step toward a scaled gradient step.
static bool impMStepDampedStep(const arma::vec& g, const arma::mat& H, arma::vec& step) {
  double hscale = H.is_finite() ? arma::abs(H.diag()).max() : 0.0;
  if (!R_finite(hscale) || hscale <= 0.0) hscale = 1.0;
  // Start AT lambda = 1e-8, not at 0.  Reaching here means H already
  // failed the conditioning test, so an undamped try is not a cheap first
  // guess: arma::solve() answers a singular system with an APPROXIMATE
  // solution (and prints "system is singular; rcond: ...; attempting
  // approx solution" to the user's console), and that approximate step is
  // then accepted whenever it happens to land inside the trust region.
  // The trust region bounds a step's SIZE; it says nothing about whether
  // the Hessian behind it determined its direction.
  double lambda = 1e-8;
  for (int tryK = 0; tryK < 12; ++tryK) {
    arma::mat Hd = H;
    Hd.diag() += lambda * hscale;
    arma::vec cand;
    // no_approx: a singular system must FAIL so the ridge escalates,
    // rather than silently returning a least-squares answer.
    if (arma::solve(cand, Hd, g, arma::solve_opts::no_approx) &&
        cand.is_finite() &&
        impStructStepRel(cand) <= IMP_MSTEP_TRUST) {
      step = cand;
      return true;
    }
    lambda *= 100.0;
  }
  return false;
}

// The Newton step on the non-mu structural thetas is only as trustworthy
// as the IS-weighted Gauss-Newton Hessian behind it.  When the importance
// weights degenerate -- a badly covered subject, or a MAP that moved
// because mceta= changed which mode the inner search starts from -- H goes
// near-singular.  arma::solve still SUCCEEDS on such an H and still returns
// a FINITE step; it is just astronomically large, so a bare
// `solve() && is_finite()` guard lets one EM iteration throw the whole
// parameter vector off the map.  Measured on Bauer's gamma model at
// mceta=10: a single iteration moved log(relative variance) from -2.46 to
// +60.3 (a relative variance of 1.6e26) and the run reported a garbage
// objective, while the same fit at mceta=0 was fine.
//
// So the step is damped rather than merely finite-checked:
//   * Levenberg-Marquardt on the diagonal, escalated until the step lands
//     inside the trust region.  H is Gauss-Newton, hence positive
//     semi-definite, so adding lambda*diag-scale only ever improves its
//     conditioning and rotates the step toward a scaled gradient step.
//   * The trust region is relative to each theta's own magnitude (floored
//     at 1) so it means the same thing for a theta of 0.1 and one of 1e5.
// The undamped Newton step is taken whenever the Hessian behind it is
// actually trustworthy, so a healthy iteration is bit-identical to the
// unguarded code.  "Trustworthy" is measured directly, by conditioning --
// NOT by how big the step is, since a large step off a well-determined
// Hessian is a legitimate move and damping it would slow every ordinary
// fit.  The sanity bound alongside it is deliberately loose: it exists to
// catch a step no honest M-step would ever propose, not to shape normal
// ones.
static void impMStepThetaStep(const arma::vec& g, const arma::mat& H) {
  arma::vec step;
  double rc = 0.0;
  if (H.is_finite()) rc = arma::rcond(H);
  bool wellCond = R_finite(rc) && rc > 1e-10;
  if (wellCond && arma::solve(step, H, g) && step.is_finite() &&
      impStructStepRel(step) <= IMP_MSTEP_SANE) {
    impUpdateStructThetas(step);
  } else if (impMStepDampedStep(g, H, step)) {
    nMStepDamped++;
    impUpdateStructThetas(step);
  } else {
    // Even heavily damped the step will not come inside the trust region
    // (or will not solve at all).  Leaving the thetas alone is the safe
    // move -- this iteration's E-step simply does not inform them -- and
    // it is reported rather than silent, since a run that spends every
    // iteration here is not estimating those thetas at all.
    nMStepSkipped++;
  }
}

// M-step.  First a Newton step on the non-mu structural thetas from the
// IS-weighted score and Gauss-Newton Hessian accumulated over subjects/samples
// -- done before the mu updates (which shift thetas/etas) so it sees the
// E-step parameters, and before impSetEta since impThetaScore's re-solves
// overwrite the etas.  Skipped if the Hessian is not usable (thetas unchanged).
static void impMStepTheta(ImpOuterState& st, int iter) {
  OdeSwapEsBatch thetaEvents(impCombSensEnabled() ? odeSlotInner : odeSlotThetaSens);
  arma::vec g(st.nSens, arma::fill::zeros);
  arma::mat H(st.nSens, st.nSens, arma::fill::zeros);
  // Batched gradient pass over the theta-sensitivity model (impThetaScore
  // switches ind->neqOverride to thetaSensNeq); single shared pool, no swap.
  // For a mixture the loop runs over the expanded pseudo-subjects and the
  // component weight a_ij is already folded into sampZk, so this forms the
  // component-weighted score directly.
  //
  // The score/Hessian accumulation is one long sequential floating-point sum
  // threaded across subjects, so it is kept serial (in eid order) for bit-
  // identity; only the solve-heavy half (impThetaSensCollect) is parallelized.
  // Each subject's samples (SIR-resampled if enabled) and collected sensitivity
  // outputs are stashed per subject, then accumulated in eid order.  cores
  // already carries the mixture / pool-sizing serial guard (forced to 1).
  if (st.harvestSens) {
    impMStepThetaHarvested(st, g, H);
  } else {
    impMStepThetaCollect(st, iter, g, H);
  }
  g /= (double)st.nsub;
  H /= (double)st.nsub;
  // MAP correction: fold in the ini({}) prior's score/curvature for these
  // thetas, if any (no-op otherwise) -- see impPriorStructThetaCorrect().
  impPriorStructThetaCorrect(st.nsub, g, H);
  impMStepThetaStep(g, H);
}

// Mixture: mean-posterior EM update of the $MIX proportions (the stable M-step
// SAEM uses).  For constant proportions the fixed point of the NONMEM
// Gauss-Newton score d_ij = a_{Nm,i}/a_Nm - a_ji/a_j is exactly a_j = mean_i
// a_ij, and this update stays in the simplex (a full Newton step overshoots to
// the boundary and collapses a component).  Install the target proportions as
// absolute $MIX thetas via the multinomial logit theta_m = log(a_m / a_Nm),
// floored away from 0 so a transiently-empty component can recover.  Uses
// impmap's own responsibilities a_ij, separate from FOCEI's Laplace mixProb.
static void impMStepMix(ImpOuterState& st) {
  const int Nmix = st.Nmix;
  arma::vec aStar(Nmix, arma::fill::zeros);
  for (int i = 0; i < st.nsub; ++i)
    for (int j = 0; j < Nmix; ++j) aStar[j] += st.aMat(i, j);
  aStar /= (double)st.nsub;
  double aFloor = 1e-3;
  for (int j = 0; j < Nmix; ++j) if (aStar[j] < aFloor) aStar[j] = aFloor;
  aStar /= arma::accu(aStar);
  arma::vec thetaM(Nmix - 1);
  for (int m = 0; m < Nmix - 1; ++m) thetaM[m] = std::log(aStar[m] / aStar[Nmix - 1]);
  if (thetaM.is_finite()) impSetMixThetas(thetaM);
}

// Seed each subject's eta with its conditional mean, then update the mu-
// referenced population parameters -- covariate groups by regression
// (updateMuGroups) and simple intercepts by the mean-shift (impMuInterceptStep)
// -- both of which recenter the etas to mean-zero residuals.  Omega is then the
// average recentered conditional moment, masked to the estimated structure.
static void impMStepMuOmega(ImpOuterState& st, bool burnIter) {
  const int nsub = st.nsub, neta = st.neta;
  for (int id = 0; id < nsub; ++id) {
    arma::vec cm = st.condMean.row(id).t();
    impSetEta(id, cm);
  }
  impUpdateMuThetas();
  impMuInterceptStep();

  arma::vec r(neta);
  arma::mat Omega(neta, neta, arma::fill::zeros);
  for (int id = 0; id < nsub; ++id) {
    impGetEta(id, r);
    Omega += r * r.t() + st.condVar[id];
  }
  Omega /= (double)nsub;
  // MAP correction: fold in the ini({}) prior's omega term(s), if any (no-op
  // otherwise) -- see impPriorOmegaCorrect().  Runs BEFORE the structure
  // mask/fixed-restore below, since a generic (non-invWishart) term's
  // correction can otherwise leak weight outside the declared structure or
  // a fix()ed row.
  if (!impPenOmegaCorrect(Omega, nsub, st.omMask, st.omFixedEta, st.Om0)) impPriorOmegaCorrect(Omega, nsub);
  Omega %= st.omMask;
  // Restore fix()ed Omega rows/columns to their starting values.
  for (size_t k = 0; k < st.omFixedEta.size(); ++k) {
    int fi = st.omFixedEta[k];
    if (fi >= 0 && fi < neta) { Omega.row(fi) = st.Om0.row(fi); Omega.col(fi) = st.Om0.col(fi); }
  }
  // burnFreezeOmega: install nothing, so Omega (and every quantity
  // impSetOmega rebuilds from it -- omegaInv, cholOmegaInv, logDetOmegaInv5,
  // the Omega thetas in fullTheta) stays at its starting value.  The thetas
  // still move; the whole M-step above ran.
  if (!(burnIter && st.burnFreezeOmega)) st.omFloored = impSetOmega(Omega, st.diagXform);
}

// Mean of trace[lo, mid) and trace[mid, hi): the drift of the window mean
static double impWindowMeanDrift(const std::vector<double>& tr, int n, int nConvWindow,
                                 double& mOld) {
  int hw = nConvWindow / 2;              // >= 1, and nConvWindow - hw >= 1
  double mNew = 0.0;
  mOld = 0.0;
  for (int k = n - hw; k < n; ++k) mNew += tr[k];
  for (int k = n - nConvWindow; k < n - hw; ++k) mOld += tr[k];
  mNew /= (double)hw;
  mOld /= (double)(nConvWindow - hw);
  return std::fabs(mNew - mOld);
}

// Objective half of the convergence test
static double impConvObjMetric(const ImpOuterState& st, int n) {
  if (impGammaRuleTarget() && st.nConvWindow >= 2) {
    // TARGET rule: mean|delta obj| is a NOISE measure, not a drift measure --
    // it has a floor at the Monte-Carlo noise level and never reaches ctol no
    // matter how settled the fit is.  The target rule deliberately widens the
    // proposal to put xi on iaccept, which costs effective sample size and so
    // RAISES that floor: measured on the 3-eta theophylline fixture,
    // mean|delta| is 0.00055-0.00093 under "floor" (ESS 0.95-0.96, all below
    // a ctol of 1e-3) against 0.00181-0.00264 under "target" (ESS 0.69-0.71,
    // all above it), so the fit could never converge however long it ran.
    // Use the drift of the window MEAN instead -- the same noise-robust form
    // the gamma gate uses.  Scoped to this rule so the floor rule's
    // tuned behaviour is bit-identical.
    double mOld;
    return impWindowMeanDrift(st.objTrace, n, st.nConvWindow, mOld) /
      std::max(1.0, std::fabs(st.obj));
  }
  double s = 0.0;
  for (int k = n - st.nConvWindow; k < n; ++k) s += std::fabs(st.objTrace[k] - st.objTrace[k - 1]);
  return (s / (double)st.nConvWindow) / std::max(1.0, std::fabs(st.obj));
}

// Proposal-scale half of the convergence test.  Under "individual" the traced
// scalar is a MEAN, so offsetting moves across subjects could look settled while
// individual scales are still swinging.  Gate on the largest per-subject step the
// controller actually applied instead, which cannot be masked that way.
// gammaRule="target" tracks a MONTE-CARLO statistic, so gamma keeps jittering
// around its fixed point and never stops moving the way the one-sided floor
// rule's gamma does (which is literally constant on a healthy model, making
// the window test vacuous).  Gate it on the drift of the window MEAN
// instead: noise averages out over nConvWindow, a genuine trend does not.
static bool impConvGammaStable(const ImpOuterState& st, int n) {
  if (st.gammaInd) return st.gammaStepPrev <= 1e-3;
  if (impGammaRuleTarget() && st.nConvWindow >= 2) {
    // Split the window in half and compare the two means.  BOTH halves must be
    // non-empty: at nConvWindow=1 the older half has no elements, mOld stays 0,
    // and the test degenerates to gamma <= gammaTargetTol -- unsatisfiable,
    // since gamma is bounded below by iscaleMin (0.1 by default), so the fit
    // could never converge.  Hence the nConvWindow >= 2 guard here and the
    // fall-through to the instantaneous test below.
    double mOld;
    double d = impWindowMeanDrift(st.gammaTrace, n, st.nConvWindow, mOld);
    return d <= impGammaTargetTol * std::max(1.0, mOld);
  }
  // Instantaneous test.  The floor rule's gamma is exactly constant once
  // coverage is healthy so 1e-3 is ample; a target-tracking gamma carries
  // xi's Monte-Carlo noise, so it gets the same looser tolerance the
  // window-mean test uses (without it, nConvWindow=1 could never converge).
  double gWin0 = st.gammaTrace[n - st.nConvWindow - 1];
  double tol = impGammaRuleTarget() ? impGammaTargetTol : 1e-3;
  return std::fabs(st.gamma - gWin0) <= tol * std::max(1.0, gWin0);
}

// Windowed convergence check (NONMEM-style CTYPE).  Requires all of:
//  (a) the mean absolute objective change over the trailing nConvWindow
//      iterations, relative to the current objective, is below ctol
//      (window-averaging suppresses the per-iteration Monte-Carlo noise);
//  (b) the estimates have stopped drifting -- the max relative net change of
//      any parameter across the window is below ctol -- so a low-leverage
//      parameter (e.g. a small-endpoint sigma barely moving the objective)
//      still converging does not trip an objective-only test;
//  (c) the proposal scale gamma has settled, so objective drift while gamma
//      is still adapting is not mistaken for convergence.
//  (d) the whole trailing window lies AFTER the burn-in.  Not just "do not
//      converge while burning in": a window straddling the boundary would
//      average burn-in iterations, whose Omega is frozen under
//      burnFreezeOmega, so the parameter-drift half (b) would read a drift
//      that no M-step was allowed to produce.  objTrace gets exactly one
//      push per iteration, so at nBurn = 0 this is identical to the size
//      test it sits beside and the default path is unchanged.
static bool impOuterConverged(const ImpOuterState& st, int iter, const arma::vec& parNow) {
  if (st.frozen || st.nConvWindow <= 0 || !R_finite(st.obj) ||
      iter < st.nBurn + st.nConvWindow ||
      (int)st.objTrace.size() < st.nConvWindow + 1) {
    return false;
  }
  int n = (int)st.objTrace.size();
  double objMetric = impConvObjMetric(st, n);
  const arma::vec& parOld = st.parHist[n - st.nConvWindow - 1];
  double parMetric = arma::max(arma::abs(parNow - parOld) / (arma::abs(parOld) + 1e-8));
  // A settled parameter still has ~1-2% Monte-Carlo net drift across the
  // window, so this gate sits above that floor: it is a safety against
  // premature stopping on a systematically-drifting low-leverage parameter,
  // not a precision requirement (the objMetric gate carries the precision).
  double parTol = 0.02;
  if (impIsIts()) {
    // deterministic iteration: relative change (floored at 1) below ctol
    parMetric = arma::max(arma::abs(parNow - parOld) /
                          arma::clamp(arma::abs(parOld), 1.0, arma::datum::inf));
    parTol = st.ctol;
  }
  bool gammaStable = impConvGammaStable(st, n);
  return gammaStable && objMetric < st.ctol && parMetric < parTol;
}

// gammaMethod="individual": xi_i falls monotonically as the proposal
// widens, so xi_i ABOVE the target means subject i's proposal is too
// narrow for its posterior (the heavy-tail case) and gamma_i must grow;
// xi_i below target means it is needlessly over-dispersed and gamma_i
// shrinks.  Two-sided, so a bad early iteration is not paid for forever.
//
// The damping exponent MUST scale with the number of random effects.
// For a Gaussian posterior matched by H_i the statistic is
// xi = gamma^(-neta/2), so the update gamma' = gamma * (xi/iaccept)^p
// linearizes (in log gamma) to an error multiplier
//
//   lambda = 1 - p * neta / 2
//
// A FIXED p (e.g. the natural-looking sqrt, p = 1/2) is therefore only
// stable for neta < 8: at neta = 8 lambda = -1 and the controller enters
// a period-2 limit cycle that never settles, so the fit can never meet
// the gamma-stability gate.  That is not a hypothetical -- an 8-eta model
// oscillates gamma 1.32/1.24/1.30/1.23... and xi 0.36/0.45/0.37/0.46...
// indefinitely.  Models with 8+ random effects are ordinary here.
//
// p = 2/neta gives lambda = 0 at every dimension.  This is not tuned
// damping: xi = gamma^(-neta/2) inverts exactly, so
// gamma * (xi/target)^(2/neta) IS the analytic solve for the gamma that
// hits the target -- a FIXED exponent is the thing that needs justifying.
// The symmetric 1.25x cap still bounds a single noisy iteration.  Weights
// correct for gamma, so this moves variance, not the estimates.
static double impGammaIndFactor(double xiId, double target, double pExp) {
  const double capUp = 1.25, capDn = 1.0 / 1.25;
  if (!R_finite(xiId)) return capUp;      // xi = +Inf: as under-dispersed as it gets
  if (xiId <= 0.0) return capDn;          // xi = 0: proposal far too wide
  double fac = std::pow(xiId / target, pExp);
  if (fac > capUp) fac = capUp;
  else if (fac < capDn) fac = capDn;
  return fac;
}

static void impAdaptGammaInd(ImpOuterState& st) {
  const double pExp = 2.0 / std::max(1.0, (double)st.neta);
  double maxStep = 0.0;
  for (int id = 0; id < st.nExp; ++id) {
    double xiId = st.XiExp[id];
    if (ISNAN(xiId)) continue;      // no usable xi (dead subject): leave alone
    double fac = impGammaIndFactor(xiId, st.iacceptVec[id], pExp);
    double g = st.gammaVec[id] * fac;
    if (g < st.iscaleMin) g = st.iscaleMin;
    if (g > st.iscaleMax) g = st.iscaleMax;
    // Measure the step actually taken (after clamping) so a subject pinned
    // at a bound reads as settled rather than as perpetually moving.
    double applied = (st.gammaVec[id] > 0.0) ? std::fabs(g / st.gammaVec[id] - 1.0) : 0.0;
    if (applied > maxStep) maxStep = applied;
    st.gammaVec[id] = g;
  }
  st.gammaStepPrev = maxStep;
}

// gammaRule="target" -- NONMEM's rule for the SHARED scale.  The NM7
// Technical Guide (note after eq. 1.76) says gamma is "continually adjusted
// so that xi_i approximates IACCEPT", i.e. two-sided on xi, not one-sided on
// the Kish fraction.  Reuse the per-subject branch's analytic inversion:
// xi = gamma^(-neta/2) inverts exactly, so pExp = 2/neta gives an error
// multiplier of 0 at every dimension.  A FIXED exponent (the sqrt the
// "floor" branch uses) is the form measured to enter a period-2 limit
// cycle at neta >= 8, which one-sidedness is currently the only thing
// suppressing -- so the two-sided rule must NOT reuse it.
//
// xiMean is NA when no subject had a usable xi, and a single non-finite
// subject would otherwise poison gamma for the rest of the fit through
// pow(NaN, .), so guard before applying.
//
// LIMIT OF THE SHARED FORM.  The exponent 2/neta inverts xi = gamma^(-neta/2)
// exactly for ONE posterior.  Over heterogeneous subjects sharing one gamma
// the mean responds with an effective exponent sum(w_i d_i); if a
// heavy-tailed subject with d_i > neta/2 dominates the mean enough that
// sum(w_i d_i) > neta, the linearized map has eigenvalue < -1 and gamma
// oscillates between the caps instead of settling.  The 1.25x caps bound it
// (it cannot diverge) and the window-mean gate will simply not report
// convergence.  gammaMethod="individual" has no such exposure -- each
// subject inverts its own xi -- so that is the remedy when it bites.
static void impAdaptGammaTarget(ImpOuterState& st) {
  if (!(st.iaccept > 0 && R_finite(st.xiMean) && st.xiMean > 0.0)) return;
  const double pExp = 2.0 / std::max(1.0, (double)st.neta);
  double fac = std::pow(st.xiMean / st.iaccept, pExp);
  if (fac > 1.25) fac = 1.25;
  else if (fac < 1.0 / 1.25) fac = 1.0 / 1.25;
  double g = st.gamma * fac;
  if (g < st.iscaleMin) g = st.iscaleMin;
  if (g > st.iscaleMax) g = st.iscaleMax;
  // measure the step AFTER clamping, so a gamma pinned at a bound reads as
  // settled rather than as perpetually moving (mirrors gammaStepPrev in the
  // per-subject branch); the convergence gate compares gamma across the
  // window, which a clamped gamma satisfies for free.
  st.gamma = g;
}

// gammaRule="floor" (default) -- iaccept is a one-sided FLOOR on the mean
// Kish effective-sample fraction: gamma is left at its efficient starting
// value while coverage is healthy and inflated only when it drops below.
static void impAdaptGammaFloor(ImpOuterState& st) {
  if (!(st.iaccept > 0 && st.accFrac > 0 && st.accFrac < st.iaccept)) return;
  double fac = std::sqrt(st.iaccept / st.accFrac);
  if (fac > 1.25) fac = 1.25;
  st.gamma *= fac;
  if (st.gamma < st.iscaleMin) st.gamma = st.iscaleMin;
  if (st.gamma > st.iscaleMax) st.gamma = st.iscaleMax;
}

// Adapt the proposal scale gamma (NONMEM ISCALE) treating iaccept as an
// effective-sample-fraction FLOOR, not a rigid target: the effective sample
// size is maximized near gamma = 1 (the Laplace proposal matches a Gaussian
// posterior), so deliberately forcing accFrac down to iaccept would only add
// Monte-Carlo noise.  gamma is left at the efficient starting value while
// coverage is healthy (accFrac >= iaccept) and inflated only when coverage
// drops below the floor (heavy-tailed/skewed posterior), with a gentle
// per-step cap and clamped to [iscaleMin, iscaleMax].  The importance weights
// correct for gamma, so this changes only the variance, not the estimates.
//
// Iteration 0 does not adapt.  Its weights are normalized against the STARTING
// mode/Hessian at the initial estimates, which is not yet a meaningful
// reference, so both statistics the controller reads are transients: measured
// on theophylline, xi = 4.1e6 against a settled ~0.39 and accFrac = 0.257
// against a settled ~0.53.  Adapting off that inflated gamma to 1.248
// (sqrt(0.4/0.257), just under the 1.25 cap) on a problem whose proposal is
// healthy from iteration 1 onward -- and the "global" branch is
// one-sided, so it could never come back down.
static void impAdaptGamma(ImpOuterState& st, int iter) {
  if (iter == 0) return;   // nothing to adapt from yet
  if (st.gammaInd) {
    impAdaptGammaInd(st);
  } else if (impGammaRuleTarget()) {
    impAdaptGammaTarget(st);
  } else {
    impAdaptGammaFloor(st);
  }
}

// One EM iteration; true when the loop stops here.
static bool impOuterIter(ImpOuterState& st, int iter) {
  const bool burnIter = (iter < st.nBurn);
  impOuterPreEStep(st, iter);
  // Stash the E-step diagnostics on every iteration so the fit environment
  // reflects the last iteration actually run (the loop may stop early).
  impEStep(st.nsub, st.neta, st.isampleVec, st.gammaVec, st.props, st.cores, iter,
           impLogDetOmegaInv5(), st.isImp,
           st.condMean, st.condVar, st.Li, st.Neff, st.Xi, st.NeffExp, st.XiExp, st.KhatExp,
           st.sampS, st.sampZk, st.aMat, &st.e, st.qrPinSeed,
           st.harvestSens, st.outSens);
  impOuterSummarize(st);
  if (st.autoOn) impAutoEscalate(st);
  if (st.autoOn && st.nExp > 1) impAutoRealloc(st);
  st.objTrace.push_back(st.obj);
  st.gammaTrace.push_back(st.gamma);
  st.neffTrace.push_back(st.accFrac);
  st.xiTrace.push_back(st.xiMean);
  st.iterRun = iter + 1;

  if (st.eOnly) {
    arma::vec parNow; impGetEstPar(parNow);
    impIterPrintRow(parNow, st.obj);
    return true;
  }
  if (st.nSens > 0 && !st.frozen) impMStepTheta(st, iter);
  if (st.Nmix > 1 && !st.frozen) impMStepMix(st);
  if (!st.frozen) impMStepMuOmega(st, burnIter);

  // Record the current estimates for the parameter-stability half of the test.
  arma::vec parNow; impGetEstPar(parNow);
  st.parHist.push_back(parNow);
  // Print this iteration + record it into the parameter-history walk.
  impIterPrintRow(parNow, st.obj);
  if (impOuterConverged(st, iter, parNow)) {
    st.converged = true;
    return true;
  }
  impAdaptGamma(st, iter);
  return false;
}

// Finalize the fit at the converged estimates.
static void impOuterFinalize(ImpOuterState& st) {
  // Close the iteration print and stash the parameter-history walk (e$parHistData).
  impIterPrintGet(st.e);
  // Drop the objective stashed by the initial MAP pass first: the finalize path
  // keeps an existing e$objective (re-adjusting it as if it were an unadjusted
  // -2LL), which would publish the initial-parameter objective instead of the
  // converged FOCEi evaluation.
  if (st.e.exists("objective")) st.e.remove("objective");
  // the omega the last M-step installed is the one reported
  impWarnOmegaFloor(st.e, st.omFloored);
  impSyncInitParToFullTheta();
  impMapPass(st.e);

  // Monte-Carlo observed-information covariance (theta block) at the converged
  // estimate, before the inner neqOverride is cleared (the inner solves use it).
  // Selected by impmapControl(covMethod="imp") (the default).
  // Snapshots, not the live vectors: the adaptation steps have already moved
  // gammaVec/dfVec on for an iteration that never ran, so passing them would
  // evaluate the covariance at a proposal the fit never actually used (and
  // disagree with the reported impDfInd / impGammaInd).
  if (impCovEnabled()) {
    // the finite differences of the objective run at the covariance probe
    // tolerances, as FOCEi's covariance step does
    ImpCovSolveTolScope _tol;
    impComputeCov(st.e, st.gammaUsed, st.propsUsed, impQrRefreshEnabled() ? st.iterRun - 1 : 0);
  }

  // Clear the multi-endpoint inner neqOverride so it does not leak into a later fit.
  impClearInnerNeqOverride();
}

// Stash the last-iteration E-step diagnostics.
static void impOuterStashEStep(ImpOuterState& st) {
  Environment& e = st.e;
  List condVarList(st.nsub);
  for (int id = 0; id < st.nsub; ++id) condVarList[id] = wrap(st.condVar[id]);
  e["impCondMean"] = wrap(st.condMean);
  e["impCondVar"]  = condVarList;
  e["impLi"]       = wrap(st.Li);
  e["impNeff"]     = wrap(st.Neff);
  e["impXi"]       = wrap(st.Xi);
  // Per-expanded-subject proposal scales actually used, and which controller
  // produced them.  Under "global" every element equals impGammaUsed.
  e["impGammaInd"] = wrap(st.gammaUsed);
  e["impGammaMethod"] = st.gammaInd ? "individual" : "global";
  e["impObj"]      = st.obj;
  // Snapshot, so the scalar cannot disagree with the impGammaInd vector: under
  // gammaMethod="global" the scalar gamma is adapted at the END of the loop for
  // an iteration that never runs, which would report gamma_start*fac here
  // against gamma_start in the vector.
  e["impGammaUsed"] = st.gammaScalarUsed;
  e["impNsample"]  = st.isample;
  e["impNsampleInd"] = wrap(st.isampleUsed);  // counts USED by the last E-step
  e["impDfInd"]    = wrap(st.dfUsed);          // df USED by the last E-step
  e["impIacceptInd"] = wrap(st.iacceptVec);    // per-expanded-subject acceptance target
  e["impKhatIter"] = wrap(st.KhatExp);         // last-iteration in-kernel k-hat
  e["impAuto"]     = st.autoOn;
  e["impQr"]       = impQrEnabled();
}

// Resolved proposal family, and the per-expanded-subject families actually
// used by the last E-step (so an AUTO ladder that moved a subject
// normal -> t is visible).
static void impOuterStashProposal(ImpOuterState& st) {
  Environment& e = st.e;
  const char* pn[4] = {"normal", "t", "laplace", "mixture"};
  e["impProposal"] = std::string(pn[st.propSpec.type]);
  CharacterVector pind(st.nExp);
  for (int id = 0; id < st.nExp; ++id) pind[id] = pn[st.propsUsed[id].type];
  e["impPropInd"] = pind;
  if (st.propSpec.type == impPropMix) {
    NumericVector cs(st.propSpec.nmix), ws(st.propSpec.nmix);
    for (int k = 0; k < st.propSpec.nmix; ++k) {
      cs[k] = st.propSpec.c[k];
      ws[k] = st.propSpec.w[k];
    }
    e["impPropMixScale"] = cs;      // as given: see impPropSpecFromControl
    e["impPropMixWeight"] = ws;
  }
  int scr = impQrScramble();
  const char* scrName = "none";
  if (scr == impQrScrambleOwen) scrName = "owen";
  else if (scr == impQrScrambleLms) scrName = "lms";
  e["impQrScramble"] = scrName;
  e["impSir"]      = impSirEnabled();
  e["impSirSample"] = impSirN();
}

// An M-step that needed damping (or had to be skipped) on most iterations
// means the theta-sensitivity Hessian is ill-conditioned under the current
// importance weights, and the structural thetas are being driven by a
// regularized step rather than a Newton step.  That is far better than the
// blow-up it replaces, but it is not something to leave unsaid.
static void impOuterMStepNote(int iterRun) {
  int nd = nMStepDamped.load(std::memory_order_relaxed);
  int ns = nMStepSkipped.load(std::memory_order_relaxed);
  if (ns > 0) {
    RSprintf("imp: %d of %d M-step iterations could not update the structural thetas (ill-conditioned theta Hessian)\n",
             ns, iterRun);
  } else if (nd > 0 && iterRun > 0 && nd * 2 > iterRun) {
    RSprintf("imp: %d of %d M-step Newton steps needed damping (ill-conditioned theta Hessian)\n",
             nd, iterRun);
  }
}

static void impOuterStashRun(ImpOuterState& st) {
  Environment& e = st.e;
  const int nIterRep = st.eOnly ? 0 : st.nIterTotal;
  const int iterRep = st.eOnly ? 0 : st.iterRun;
  e["impNiter"]    = nIterRep;                    // includes nBurn: the budget impIter is measured against
  e["impEonly"]    = st.eOnly;
  e["impMapIter"]  = st.isImp ? 0 : st.mapIter;   // 0 under est="imp": no re-centering at all
  // Burn-in iterations are the FIRST nBurn rows of every trace and of $parHist.
  e["impNburn"]    = st.nBurn;
  e["impBurnFreezeOmega"] = st.burnFreezeOmega;
  e["impIter"]     = iterRep;
  e["impMStepDamped"]  = nMStepDamped.load(std::memory_order_relaxed);
  e["impMStepSkipped"] = nMStepSkipped.load(std::memory_order_relaxed);
  impOuterMStepNote(st.iterRun);
  e["impConverged"] = st.converged;
  e["impObjTrace"] = wrap(st.objTrace);
  e["impGammaTrace"] = wrap(st.gammaTrace);
  e["impNeffFrac"] = wrap(st.neffTrace);
  e["impXiTrace"]  = wrap(st.xiTrace);
  e["impMuGroupN"] = impMuGroupN();
}

void impOuter(Environment e) {
  // Model registration may leave the theta-sensitivity event callbacks installed.
  OdeSwapEsBatch innerEvents(odeSlotInner);
  ImpOuterState st(e);
  impOuterReadControls(st);
  impOuterInitVectors(st);
  impOuterTileIsample(st);
  st.gammaInd = impGammaIndividual();

  // Force the EM serial when the per-subject inner solves are not thread-safe:
  //  (a) mixtures -- the expanded pseudo-subjects' per-component solves race and
  //      the proportions/Omega drift and collapse;
  //  (b) multi-endpoint pool-sizing -- the inner MAP runs with ind->neqOverride
  //      against a pool sized for the larger theta-sensitivity model, whose
  //      per-thread work arrays are not safe under the override, so the parallel
  //      E-step non-deterministically rejects a subject's samples (neff collapse).
  // The plain single-endpoint path keeps full parallelism; restored before return.
  if ((st.Nmix > 1 || impPoolSizing()) && st.cores > 1) {
    st.savedCores = impSetSolveCores(1);
    st.cores = 1;
  }

  // Initial MAP at the starting parameters.
  impMapPass(e);

  // NB placed AFTER the initial MAP pass on purpose: fInd->nObs is only filled
  // in during a solve, so reading it before impMapPass() returns 0 for every
  // subject and makes the whole dataset look sparse (which silently gave every
  // subject df = 4).
  st.nonmemSparse = impAutoNonmemSparse();
  st.dfPatience = impAutoDfPatience();
  if (st.autoOn) impAutoStep1(st);

  impOuterSetupSampling(st);
  impOuterSetupOmega(st);
  nMStepDamped = 0;
  nMStepSkipped = 0;

  // Iteration print + parameter-history capture (shared scale.h machinery).
  impIterPrintStart();
  for (int iter = 0; iter < st.nIterTotal; ++iter) {
    if (impOuterIter(st, iter)) break;
  }

  impOuterFinalize(st);
  impOuterStashEStep(st);
  impOuterStashProposal(st);
  impOuterStashRun(st);

  // Restore the solve core count if it was forced serial for the mixture EM.
  if (st.savedCores > 0) impSetSolveCores(st.savedCores);
}
