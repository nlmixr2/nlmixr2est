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
#include <ctime>
#include "nmMcmcRng.h"
#include "imp.h"
#include "utilc.h"   // RSprintf (covariance-step progress header)
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// ---- quasi-random (QRPEM) sampling elements --------------------------------
// Base Sobol point set on (0,1)^neta: isample rows, one point per row.  The
// boost engine emits one coordinate per call, cycling through the dimensions
// of consecutive points, and skips the trivial zero point; the top 53 bits are
// scaled into (0,1) with a half-step offset so no coordinate is ever 0 or 1.
static arma::mat impSobolU0(int isample, int neta) {
  boost::random::sobol eng((std::size_t)neta);
  arma::mat U0(isample, neta);
  for (int k = 0; k < isample; ++k) {
    for (int j = 0; j < neta; ++j) {
      uint64_t v = (uint64_t)eng();
      U0(k, j) = std::ldexp((double)(v >> 11) + 0.5, -53);
    }
  }
  return U0;
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

// Same, but from a supplied uniform (quasi-random path).
static inline double impTScaleU(double u, double df) {
  if (df <= 0.0) return 1.0;
  double w = impChisqQuantile(u, df);
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

// Test hook: the (optionally shifted) quasi-random N(0,1) point set.
//[[Rcpp::export]]
NumericMatrix impQrPoints_(int isample, int neta, Nullable<NumericVector> shift) {
  if (isample < 1 || neta < 1) stop("'isample' and 'neta' must be positive");
  arma::mat U0 = impSobolU0(isample, neta);
  if (shift.isNotNull()) {
    NumericVector s(shift);
    if ((int)s.size() != neta) stop("'shift' must have length 'neta'");
    arma::vec sv(s.begin(), neta);
    return wrap(impQrZ(U0, &sv));
  }
  return wrap(impQrZ(U0, nullptr));
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
                     const arma::vec& dfVec, int cores,
                     int iter, double negHalfLogDetOmega, bool isImp,
                     arma::mat& condMean, std::vector<arma::mat>& condVar,
                     arma::vec& Li, arma::vec& Neff, arma::vec& Xi,
                     arma::vec& XiExpOut, arma::vec& KhatExpOut,
                     std::vector<arma::mat>& outS, std::vector<arma::vec>& outZk,
                     arma::mat& aMat, Environment* eStash,
                     uint32_t qrPinSeed) {
  int Nm = impNmix();
  int nExp = nsub * Nm; // expanded pseudo-subjects: component j (1-based) is i + (j-1)*nsub

  // Per-expanded-subject proposal (mode, information H_i, log|H_i|, lower Cholesky
  // of gamma * H_i^-1) -- from the per-component MAP modes.  Each subject fills only
  // its own slot (no reduction, no RNG), so this is parallelized over subjects; the
  // per-subject inner Hessian (impGetHessian) is the impmap serial cost.  cores
  // already carries the mixture / pool-sizing serial guard.
  std::vector<arma::vec> modes(nExp);
  std::vector<arma::mat> Hs(nExp);
  std::vector<double> logDetH(nExp, 0.0);
  std::vector<arma::mat> cholL(nExp);
  std::vector<char> haveL(nExp, 0);
  // est="imp": no MAP search.  The proposal is centered at each subject's running
  // conditional mean (the current eta) with covariance gamma * V_i, where V_i is
  // that subject's conditional variance from the PREVIOUS E-step (passed in via
  // condVar, before it is overwritten below).  On the first iteration -- or a
  // pseudo-subject / degenerate V_i -- fall back to the population Omega, which is
  // over-dispersed and always available.
  arma::mat impOmega;
  if (isImp) impGetOmega(impOmega);
  bool doParProp = (cores > 1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doParProp)
#endif
  for (int id = 0; id < nExp; ++id) {
#ifdef _OPENMP
    if (doParProp) setRxThreadId(omp_get_thread_num());
#endif
    double gammaId = gammaVec[id];
    if (isImp) {
      arma::vec eta(neta);
      impGetEta(id, eta);
      modes[id] = eta;
      arma::mat Sig;
      if ((size_t)id < condVar.size() && condVar[id].n_rows == (arma::uword)neta &&
          arma::any(arma::vectorise(condVar[id]) != 0.0)) {
        Sig = condVar[id];
      } else {
        Sig = impOmega;
      }
      arma::mat Hi;
      double ldv, lds;
      if (arma::inv_sympd(Hi, Sig) && arma::log_det(ldv, lds, Hi) && lds > 0) {
        arma::mat L;
        if (arma::chol(L, gammaId * Sig, "lower")) {
          Hs[id] = Hi; logDetH[id] = ldv; cholL[id] = L; haveL[id] = 1;
        }
      }
    } else {
      arma::vec mode(neta);
      arma::mat H(neta, neta, arma::fill::zeros);
      impGetMode(id, mode);
      modes[id] = mode;
      if (impGetHessian(id, H)) {
        arma::mat Sigma;
        double ldv, lds;
        if (arma::inv_sympd(Sigma, H) && arma::log_det(ldv, lds, H) && lds > 0) {
          Hs[id] = H;
          logDetH[id] = ldv;
          Sigma *= gammaId;
          arma::mat L;
          if (arma::chol(L, Sigma, "lower")) { cholL[id] = L; haveL[id] = 1; }
        }
      }
    }
#ifdef _OPENMP
    if (doParProp) setRxThreadId(-1);
#endif
  }

  // Per-expanded-subject E-step results (combined per base subject below).
  arma::mat cmExp(nExp, neta, arma::fill::zeros);
  std::vector<arma::mat> cvExp(nExp, arma::mat(neta, neta, arma::fill::zeros));
  arma::vec LiExp(nExp); LiExp.fill(NA_REAL);
  arma::vec NeffExp(nExp); NeffExp.fill(NA_REAL);
  arma::vec XiExp(nExp); XiExp.fill(NA_REAL);
  arma::vec KhatExp(nExp); KhatExp.fill(NA_REAL);
  std::vector<char> okExp(nExp, 0);
  outS.assign(nExp, arma::mat());
  outZk.assign(nExp, arma::vec());
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
  // QRPEM: one Sobol base point set per E-step (read-only in the parallel
  // loop); with qrShift=FALSE the N(0,1) points are also fixed and shared.
  bool qr = impQrEnabled();
  bool qrShift = impQrShiftEnabled();
  bool qrRefresh = impQrRefreshEnabled();
  arma::mat qrU0, qrZ0;
  if (qr) {
    qrU0 = impSobolU0((int)isampleVec.max(), neta);
    if (!qrShift) qrZ0 = impQrZ(qrU0, nullptr);
  }
  // Parallelize over base subjects, iterating a subject's mixture components
  // serially within one thread.  A base subject's expanded pseudo-subjects (id =
  // i + j*nsub) share the same underlying rxode2 solve rows, so evaluating two of
  // its components concurrently would race on that solve buffer -- keeping them on
  // one thread avoids the race while preserving nsub-way parallelism.  For Nm == 1
  // this is the same nsub-iteration loop as before (bit-identical, same seeds).
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doPar)
#endif
  for (int i = 0; i < nsub; ++i) {
#ifdef _OPENMP
    if (doPar) setRxThreadId(omp_get_thread_num());
#endif
    for (int j = 0; j < Nm; ++j) {
      int id = i + j * nsub;
      // fresh per-(iter, expanded-subject) stream, independent of thread count
      nmSetSeedEng1(seed0 + (uint32_t)((iter * nExp + id) * 2));
      cmExp.row(id) = modes[id].t();
      if (!haveL[id]) continue;
      // This subject's own proposal scale (all equal under gammaMethod="global").
      double gammaId = gammaVec[id];
      double dfId = dfVec[id];   // this subject's proposal df; 0 = Gaussian
      int nsId = (int)isampleVec[id];   // this subject's sample count
      arma::mat S(nsId, neta);
      if (qr) {
        // Sobol points -> N(0,I) -> proposal.  With qrShift the shift comes
        // from this subject's seeded stream; qrRefresh=FALSE instead uses the
        // fit-constant qrPinSeed (captured once in impOuter -- seed0 advances
        // every E-step) so the whole fit reuses one shift per subject.
        arma::mat Z;
        if (qrShift) {
          if (!qrRefresh) nmSetSeedEng1(qrPinSeed + (uint32_t)(id * 2));
          arma::vec sh(neta);
          for (int jj = 0; jj < neta; ++jj) sh[jj] = rxUnifEng(0.0, 1.0);
          Z = impQrZ(qrU0, &sh);
        } else {
          Z = qrZ0;
        }
        for (int k = 0; k < nsId; ++k) {
          // Quasi-random t: the Sobol point supplies the Gaussian direction and
          // a fresh uniform supplies the chi-square scale.  (A fully
          // quasi-random t would need an extra Sobol coordinate; this keeps the
          // low-discrepancy structure where it does the most good.)
          double ts = impTScale(dfId);
          S.row(k) = (modes[id] + cholL[id] * (Z.row(k).t() * ts)).t();
        }
      } else {
        for (int k = 0; k < nsId; ++k) {
          arma::vec z(neta);
          for (int jj = 0; jj < neta; ++jj) z[jj] = rxNormEng(0.0, 1.0);
          double ts = impTScale(dfId);
          S.row(k) = (modes[id] + cholL[id] * (z * ts)).t();
        }
      }
      // NONMEM-style xi (NM7 eq. 1.73/1.75): the mean importance weight
      // normalized so the proposal kernel is 1 at its center.  q at the center
      // is just the log joint density there (the kernel term vanishes at d=0),
      // so xi = exp(logMeanExp - qCenter) once logMeanExp is formed below.
      //
      // xi == 1 exactly when the proposal matches the posterior shape, < 1 when
      // the proposal is over-dispersed, and > 1 when the posterior has heavier
      // tails than the Gaussian proposal -- xi is NOT bounded above.  That
      // heavy-tail signal is what the per-subject gamma controller
      // targets, and it is a different statistic from the Kish
      // effective sample size below (they are not comparable).
      //
      // Evaluated HERE -- after the draws are complete but before the weight
      // loop -- deliberately: the sample matrix S is already filled, so this
      // cannot perturb the sampling stream, and the weight loop still leaves
      // the subject's eta at S.row(isample-1) exactly as before, preserving
      // bit-identity with the pre-xi kernel.
      double qCenter = -impEvalJointLik(modes[id], id);
      // log importance weights.  impEvalJointLik returns the NEGATIVE log joint
      // density (minimized by the inner optimizer), so it enters with a minus;
      // the proposal density is subtracted, contributing +(1/(2 gamma)) d'H d.
      arma::vec q(nsId);
      int nGood = 0;
      for (int k = 0; k < nsId; ++k) {
        arma::vec eta = S.row(k).t();
        arma::vec d = eta - modes[id];
        double qk = -impEvalJointLik(eta, id) +
          impLogKernelRecip(arma::as_scalar(d.t() * Hs[id] * d), gammaId, dfId, neta);
        // A proposal sample whose inner solve fails (NA/NaN, even after the
        // FOCEI tolerance-relaxation retry) is dropped from the weighted mean by
        // giving it -Inf log-weight (weight 0).
        if (R_finite(qk)) { q[k] = qk; ++nGood; } else { q[k] = R_NegInf; }
      }
      // the inner solves above re-seed the engine per subject
      // (setSeedEng1(getRxSeed1()+id)); restore this block's sampling seed so
      // that leak never carries into a subsequent draw.
      nmRestoreMcmcSeed();
      if (nGood == 0) continue;
      double qmax = q.max();
      arma::vec w = arma::exp(q - qmax);
      double sumw = arma::accu(w);
      arma::vec zk = w / sumw;
      outS[id] = S;
      outZk[id] = zk;
      okExp[id] = 1;
      double logMeanExp = qmax + std::log(sumw / (double)nsId);
      double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gammaId) -
        0.5 * logDetH[id] + impCiTCorrection(dfId, neta);
      LiExp[id] = -(logMeanExp + Ci);          // per-component negative log-likelihood
      NeffExp[id] = 1.0 / arma::accu(arma::square(zk));
      // NONMEM xi: mean weight normalized to the proposal center (see above).
      if (R_finite(qCenter)) XiExp[id] = std::exp(logMeanExp - qCenter);
      KhatExp[id] = impPsisKhat(zk);   // tail index of this subject's weights
      arma::rowvec pbar = zk.t() * S;
      cmExp.row(id) = pbar;
      arma::mat B(neta, neta, arma::fill::zeros);
      for (int k = 0; k < nsId; ++k) {
        arma::rowvec dr = S.row(k) - pbar;
        B += zk[k] * (dr.t() * dr);
      }
      cvExp[id] = B;
    }
#ifdef _OPENMP
    if (doPar) setRxThreadId(-1);
#endif
  }

  // Combine over mixture components per base subject.  For a single component
  // (Nm == 1) this is bit-identical to the non-mixture E-step.
  condMean.set_size(nsub, neta);
  condVar.assign(nsub, arma::mat(neta, neta, arma::fill::zeros));
  Li.set_size(nsub); Li.fill(NA_REAL);
  Neff.set_size(nsub); Neff.fill(NA_REAL);
  Xi.set_size(nsub); Xi.fill(NA_REAL);
  aMat.set_size(nsub, Nm); aMat.zeros();
  for (int i = 0; i < nsub; ++i) {
    if (Nm == 1) {
      condMean.row(i) = cmExp.row(i);
      condVar[i] = cvExp[i];
      Li[i] = LiExp[i];
      Neff[i] = NeffExp[i];
      Xi[i] = XiExp[i];
      aMat(i, 0) = 1.0;
      continue;
    }
    // Posterior mixture weights a_ij = a_j exp(-L_ij) / sum_k a_k exp(-L_ik),
    // with the min-L shift for numerical stability (as saem does).
    double minL = R_PosInf;
    for (int j = 0; j < Nm; ++j) {
      int eid = i + j * nsub;
      if (okExp[eid] && R_finite(LiExp[eid]) && LiExp[eid] < minL) minL = LiExp[eid];
    }
    if (!R_finite(minL)) { condMean.row(i) = cmExp.row(i); aMat(i, 0) = 1.0; continue; }
    arma::vec wj(Nm, arma::fill::zeros); double sumW = 0.0;
    for (int j = 0; j < Nm; ++j) {
      int eid = i + j * nsub;
      if (okExp[eid] && R_finite(LiExp[eid])) {
        wj[j] = impMixProb(j) * std::exp(minL - LiExp[eid]);
        sumW += wj[j];
      }
    }
    if (sumW <= 0.0) { condMean.row(i) = cmExp.row(i); aMat(i, 0) = 1.0; continue; }
    arma::rowvec cm(neta, arma::fill::zeros);
    for (int j = 0; j < Nm; ++j) {
      double a = wj[j] / sumW; aMat(i, j) = a;
      if (a > 0) cm += a * cmExp.row(i + j * nsub);
    }
    condMean.row(i) = cm;
    // Component-weighted conditional variance: sum_j a_ij (B_ij + m_ij m_ij') - m_i m_i'.
    arma::mat cv(neta, neta, arma::fill::zeros);
    for (int j = 0; j < Nm; ++j) {
      double a = aMat(i, j); if (a <= 0) continue;
      arma::rowvec mj = cmExp.row(i + j * nsub);
      cv += a * (cvExp[i + j * nsub] + mj.t() * mj);
    }
    cv -= cm.t() * cm;
    condVar[i] = cv;
    Li[i] = minL - std::log(sumW);           // marginal (mixture) negative log-likelihood
    double nf = 0.0;
    for (int j = 0; j < Nm; ++j) {
      double a = aMat(i, j);
      if (a > 0 && R_finite(NeffExp[i + j * nsub])) nf += a * NeffExp[i + j * nsub];
    }
    if (nf > 0) Neff[i] = nf;
    // Responsibility-weighted xi, matching how Neff is combined above.
    double xf = 0.0; bool anyXi = false;
    for (int j = 0; j < Nm; ++j) {
      double a = aMat(i, j);
      if (a > 0 && R_finite(XiExp[i + j * nsub])) { xf += a * XiExp[i + j * nsub]; anyXi = true; }
    }
    if (anyXi) Xi[i] = xf;
    // Fold the responsibility a_ij into each component's importance weights so the
    // M-step, iterating the expanded subjects, forms the component-weighted score.
    for (int j = 0; j < Nm; ++j) {
      int eid = i + j * nsub;
      if (outZk[eid].n_elem > 0) outZk[eid] *= aMat(i, j);
    }
  }

  // Per-EXPANDED-subject xi, which is what the individual-gamma controller
  // drives (one gamma_i per expanded subject).  `Xi` above is the
  // responsibility-combined per-base-subject value used for reporting.
  XiExpOut = XiExp;
  KhatExpOut = KhatExp;

  // On the last iteration, stash the per-subject sampler diagnostics (first
  // component for a mixture).
  if (eStash != nullptr) {
    arma::mat modeMat(nsub, neta);
    List hessList(nsub), samplesList(nsub), weightList(nsub);
    for (int id = 0; id < nsub; ++id) {
      modeMat.row(id) = modes[id].t();
      hessList[id] = haveL[id] ? wrap(Hs[id]) : R_NilValue;
      samplesList[id] = wrap(outS[id]);
      // Normalized importance weights, kept so a TAIL-sensitive diagnostic
      // (Pareto k-hat) can be computed post-fit.  xi and the Kish ESS are both
      // means over samples drawn FROM the proposal, so both are structurally
      // blind to tails the proposal never visits -- which is the failure mode a
      // heavier-tailed (t) proposal exists to prevent.  k-hat estimates the
      // tail index of the weight distribution and does see it.  Pareto k-hat is
      // scale-invariant, so the mixture responsibility folded into these
      // weights above does not affect it.
      weightList[id] = wrap(outZk[id]);
    }
    (*eStash)["impEtaMode"] = wrap(modeMat);
    (*eStash)["impEtaHess"] = hessList;
    (*eStash)["impSamples"] = samplesList;
    (*eStash)["impWeights"] = weightList;
  }
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
void impComputeCov(Environment e, const arma::vec& gammaVec,
                   const arma::vec& dfVec) {
  int nsub = impNsub();
  int neta = impNeta();
  int isample = impNsample();
  // Parallelize the reweighted-objective subject loop (each subject contributes an
  // independent scalar to the -2LL) when the inner solves are thread-safe: no
  // mixture (expanded pseudo-subjects share solve rows) and no multi-endpoint
  // pool-sizing (neqOverride against the larger theta-sens pool) -- the same
  // contract as the E-step.  The one-time fixed-sample draw below stays serial.
  int cores = impCores();
  if (cores < 1) cores = 1;
  bool doParCov = (cores > 1) && (impNmix() == 1) && !impPoolSizing();

  // Free parameters in the optimizer's (fixedTrans) order -- the same order the
  // fit's covariance uses.  pl[j] is a theta fullTheta index when < ntheta, else
  // the Omega parameter (pl[j] - ntheta).
  std::vector<int> pl;
  impGetCovParList(pl);
  int np = (int)pl.size();
  if (np == 0) return;
  int ntheta = impNtheta();
  int nTh = 0;
  for (int j = 0; j < np; ++j) if (pl[j] < ntheta) ++nTh;

  // Per-subject proposal: mode, information H, lower Cholesky of gamma*H^-1, and
  // one fixed sample matrix.
  std::vector<arma::vec> modes(nsub);
  std::vector<arma::mat> Hs(nsub), Ls(nsub), Ss(nsub);
  std::vector<double> logDetH(nsub, 0.0);
  std::vector<char> ok(nsub, 0);
  arma::vec mode(neta);
  arma::mat H(neta, neta, arma::fill::zeros);
  for (int id = 0; id < nsub; ++id) {
    impGetMode(id, mode);
    modes[id] = mode;
    if (impGetHessian(id, H)) {
      H = 0.5 * (H + H.t()); // guard against tiny numerical asymmetry
      arma::mat Sigma; double ldv, lds;
      if (arma::inv_sympd(Sigma, H) && arma::log_det(ldv, lds, H) && lds > 0) {
        Sigma *= gammaVec[id];
        arma::mat L;
        if (arma::chol(L, Sigma, "lower")) {
          Hs[id] = H; Ls[id] = L; logDetH[id] = ldv; ok[id] = 1;
        }
      }
    }
  }
  // Draw one fixed sample set per subject (serial, deterministic seed).  The
  // engine is already allocated + seeded by rxWithSeed() at the fit start.
  // With qr=TRUE the fixed set is the Sobol point set (per-subject shifted
  // when qrShift; qrRefresh is irrelevant here -- the set never refreshes).
  bool qr = impQrEnabled();
  bool qrShift = impQrShiftEnabled();
  arma::mat qrU0, qrZ0;
  if (qr) {
    qrU0 = impSobolU0(isample, neta);
    if (!qrShift) qrZ0 = impQrZ(qrU0, nullptr);
  }
  uint32_t seed0 = (uint32_t)impBaseSeed() + 0x2545F491u;
  setRxThreadId(0);
  for (int id = 0; id < nsub; ++id) {
    if (!ok[id]) continue;
    setSeedEng1(seed0 + (uint32_t)(id * 2 + 1));
    arma::mat S(isample, neta);
    if (qr) {
      arma::mat Z;
      if (qrShift) {
        arma::vec sh(neta);
        for (int j = 0; j < neta; ++j) sh[j] = rxUnifEng(0.0, 1.0);
        Z = impQrZ(qrU0, &sh);
      } else {
        Z = qrZ0;
      }
      for (int k = 0; k < isample; ++k) {
        S.row(k) = (modes[id] + Ls[id] * (Z.row(k).t() * impTScale(dfVec[id]))).t();
      }
    } else {
      for (int k = 0; k < isample; ++k) {
        arma::vec z(neta);
        for (int j = 0; j < neta; ++j) z[j] = rxNormEng(0.0, 1.0);
        S.row(k) = (modes[id] + Ls[id] * (z * impTScale(dfVec[id]))).t();
      }
    }
    Ss[id] = S;
  }
  setRxThreadId(-1);

  // Perturb parameter j: a theta on the subjects' parameter pointers, or an
  // Omega free parameter (which also rebuilds omegaInv / logDetOmegaInv5).
  auto setPar = [&](int j, double val) {
    if (pl[j] < ntheta) impSetThetaAll(pl[j], val);
    else impSetOmegaThetaAll(pl[j] - ntheta, val);
  };

  // Importance-sampling -2LL at a parameter vector, reusing the fixed samples.
  // Each subject's contribution is an independent scalar; accumulate them into a
  // per-subject buffer under the parallel loop, then reduce in id order so the sum
  // is bit-identical to the serial `obj += ...` accumulation.  objBuf is allocated
  // once here and refilled per call (the FD Hessian calls evalObj O(np^2) times).
  std::vector<double> objBuf(nsub, 0.0);
  // Progress bar over the finite-difference covariance evaluations, like the
  // focei covariance step.  evalObj is called f0 (1) + 2*np (diagonal) +
  // 2*np*(np-1) (off-diagonal) = 1 + 2*np*np times; tick once per call.
  bool covProg = impCovProgress();
  clock_t covT0 = clock();
  int covTot = 1 + 2 * np * np, covCur = 0, covTick = 0;
  if (covProg) RSprintf("calculating covariance matrix\n");
  auto evalObj = [&](const arma::vec& par) -> double {
    for (int j = 0; j < np; ++j) setPar(j, par[j]);
    // Re-read after setting: an Omega perturbation changes -0.5 log|Omega|.
    double negHalfLogDetOmega = impLogDetOmegaInv5();
    std::fill(objBuf.begin(), objBuf.end(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doParCov)
#endif
    for (int id = 0; id < nsub; ++id) {
#ifdef _OPENMP
      if (doParCov) setRxThreadId(omp_get_thread_num());
#endif
      if (ok[id]) {
        impForceResolve(id);
        // This subject's converged proposal scale (all equal under "global").
        // Must use the SAME proposal shape (df) as the E-step or the reweighted
        // objective is not the one the fit converged on.
        double dfCov = dfVec[id];
        arma::vec q(isample); int nGood = 0;
        for (int k = 0; k < isample; ++k) {
          arma::vec eta = Ss[id].row(k).t();
          arma::vec d = eta - modes[id];
          double qk = -impEvalJointLik(eta, id) +
            impLogKernelRecip(arma::as_scalar(d.t() * Hs[id] * d), gammaVec[id], dfCov, neta);
          if (R_finite(qk)) { q[k] = qk; ++nGood; } else q[k] = R_NegInf;
        }
        if (nGood != 0) {
          double qmax = q.max();
          double sumw = arma::accu(arma::exp(q - qmax));
          double logMeanExp = qmax + std::log(sumw / (double)isample);
          double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gammaVec[id]) - 0.5 * logDetH[id]
            + impCiTCorrection(dfCov, neta);
          objBuf[id] = -2.0 * (logMeanExp + Ci);
        }
      }
#ifdef _OPENMP
      if (doParCov) setRxThreadId(-1);
#endif
    }
    double obj = 0.0;
    for (int id = 0; id < nsub; ++id) obj += objBuf[id];
    if (covProg) covTick = par_progress(covCur++, covTot, covTick, 1, covT0, 0);
    return obj;
  };

  arma::vec par0(np);
  for (int j = 0; j < np; ++j)
    par0[j] = (pl[j] < ntheta) ? impGetFullThetaVal(pl[j])
                               : impGetOmegaThetaVal(pl[j] - ntheta);
  double f0 = evalObj(par0);
  arma::vec hstep(np);
  for (int j = 0; j < np; ++j) {
    double a = std::fabs(par0[j]);
    hstep[j] = 1e-3 * (a > 1e-3 ? a : 1.0);
  }
  arma::vec fp(np), fm(np);
  for (int j = 0; j < np; ++j) {
    arma::vec p = par0; p[j] = par0[j] + hstep[j]; fp[j] = evalObj(p);
    p = par0; p[j] = par0[j] - hstep[j]; fm[j] = evalObj(p);
  }
  arma::mat Hess(np, np, arma::fill::zeros);
  for (int j = 0; j < np; ++j)
    Hess(j, j) = (fp[j] - 2.0 * f0 + fm[j]) / (hstep[j] * hstep[j]);
  for (int a = 0; a < np; ++a) {
    for (int b = a + 1; b < np; ++b) {
      arma::vec p = par0; p[a] += hstep[a]; p[b] += hstep[b]; double fpp = evalObj(p);
      p = par0; p[a] += hstep[a]; p[b] -= hstep[b]; double fpm = evalObj(p);
      p = par0; p[a] -= hstep[a]; p[b] += hstep[b]; double fmp = evalObj(p);
      p = par0; p[a] -= hstep[a]; p[b] -= hstep[b]; double fmm = evalObj(p);
      double v = (fpp - fpm - fmp + fmm) / (4.0 * hstep[a] * hstep[b]);
      Hess(a, b) = v; Hess(b, a) = v;
    }
  }
  if (covProg) par_progress(covTot, covTot, covTick, 1, covT0, 1);   // close the bar
  // Restore the converged estimates.
  for (int j = 0; j < np; ++j) setPar(j, par0[j]);
  for (int id = 0; id < nsub; ++id) impForceResolve(id);

  // Observed information = 0.5 * Hess(-2LL) (symmetrized); covariance = inverse.
  arma::mat info = 0.25 * (Hess + Hess.t());
  arma::mat cov;
  if (!arma::inv_sympd(cov, info)) {
    if (!arma::pinv(cov, info)) { cov = arma::mat(np, np, arma::fill::zeros); cov.fill(NA_REAL); }
  }
  arma::vec se(np);
  for (int j = 0; j < np; ++j) se[j] = (cov(j, j) > 0) ? std::sqrt(cov(j, j)) : NA_REAL;
  // Full covariance in free-parameter order (matches the fit's covariance layout).
  e["impCov"] = wrap(cov);
  e["impSe"] = wrap(se);
  e["impCovThetaN"] = nTh;
  IntegerVector thIdxR(nTh);
  { int t = 0; for (int j = 0; j < np; ++j) if (pl[j] < ntheta) thIdxR[t++] = pl[j] + 1; }
  e["impCovThetaIdx"] = thIdxR;
  if (nTh > 0) {
    e["impCovTheta"] = wrap(arma::mat(cov.submat(0, 0, nTh - 1, nTh - 1)));
    e["impSeTheta"] = wrap(arma::vec(se.subvec(0, nTh - 1)));
  }
  // Publish as the fit's covariance so the standard SE / CI / correlation table
  // machinery (foceiFinalizeTables) picks it up.
  if (cov.is_finite()) {
    e["cov"] = wrap(cov);
    e["covMethod"] = CharacterVector::create("imp");
  }
}

void impOuter(Environment e) {
  int nIter = impNiter();
  std::string diagXform = impDiagXform();
  int nsub = impNsub();
  int neta = impNeta();
  int isample = impNsample();
  double gamma = impGammaProp();
  int cores = impCores();
  if (cores < 1) cores = 1;

  // Convergence controller + proposal-scale adaptation (NONMEM ISCALE/IACCEPT/CTYPE).
  double iaccept = impIaccept();
  // Drift tolerance on the gamma window mean for gammaRule="target".  Looser than
  // the floor rule's 1e-3 because that rule's gamma is exactly constant once
  // coverage is healthy, while a target-tracking gamma carries the Monte-Carlo
  // noise of xi.  PROVISIONAL -- a constant tuned for this rule, to be re-measured
  // with the rest of the target-rule set.
  const double gammaTargetTol = 1e-2;
  double iscaleMin = impIscaleMin();
  double iscaleMax = impIscaleMax();
  int nConvWindow = impNconvWindow();
  double ctol = impCtol();

  arma::mat condMean;
  std::vector<arma::mat> condVar;
  std::vector<arma::mat> sampS;
  std::vector<arma::vec> sampZk;
  arma::vec Li, Neff, Xi, XiExp, KhatExp;
  arma::mat aMat;                 // posterior mixture responsibilities (nsub x Nmix)
  int Nmix = impNmix();
  int nExp = nsub * Nmix;         // expanded pseudo-subjects for the mixture E/M-step
  double obj = R_PosInf;
  int nSens = impThetaSensN();
  // est="imp": skip the per-iteration MAP search; the E-step proposal is centered
  // at the running conditional mean with covariance gamma*Omega.
  bool isImp = impIsImp();
  // Per-expanded-subject proposal scale.  Under the current (global) rule every
  // element tracks the scalar `gamma`, so the E-step arithmetic is unchanged;
  // holding it as a vector is what lets a per-subject controller drive the
  // elements apart without touching the E-step again.
  arma::vec gammaVec(nExp); gammaVec.fill(gamma);
  // Per-expanded-subject sample count.  Uniform unless AUTO raises it for a
  // badly covered subject; holding it as a vector is what lets that happen
  // without charging every other subject for the extra samples.
  arma::ivec isampleVec(nExp); isampleVec.fill(isample);
  // Counts actually USED by the last E-step.  isampleVec is reallocated after
  // the E-step, so reporting it directly pairs next iteration's counts with
  // this iteration's Neff and can show an effective-sample FRACTION above 1,
  // which is impossible.
  arma::ivec isampleUsed(nExp); isampleUsed.fill(isample);
  arma::vec dfUsed(nExp); dfUsed.fill(impDf());   // df actually used last E-step
  // Same for gamma.  The controllers all mutate their vectors AFTER the E-step,
  // for the next iteration, so anything consuming "the converged proposal" must
  // read the snapshot rather than the live vector -- including impComputeCov.
  arma::vec gammaUsed(nExp); gammaUsed.fill(gamma);
  double gammaScalarUsed = gamma;   // scalar counterpart of gammaUsed
  // Per-expanded-subject proposal df and acceptance target.  Uniform unless
  // AUTO differentiates them.
  arma::vec dfVec(nExp); dfVec.fill(impDf());
  // Improvability state for the AUTO df escalation (see AUTO step 3).  noImp
  // counts consecutive iterations in which the current rung has failed to
  // improve on what the subject manages without any escalation.
  // Lowest df a subject may be returned to.  Withdrawal undoes the k-hat
  // ESCALATION, never the df that step 1 assigned: a non-normal endpoint needs a
  // heavy-tailed proposal by construction, so dropping it to a Gaussian because
  // k-hat did not improve would break the model, not repair it.
  arma::vec dfFloor(nExp); dfFloor.fill(impDf());
  // What this subject's k-hat looks like WITHOUT any escalation -- the
  // counterfactual withdrawal has to judge a rung against.  Refreshed on every
  // iteration the subject spends at its floor, and frozen while it is escalated.
  arma::vec kAtFloor(nExp); kAtFloor.fill(NA_REAL);
  arma::ivec noImp(nExp, arma::fill::zeros);
  arma::ivec escDead(nExp, arma::fill::zeros);   // 1 = escalation withdrawn, do not retry
  arma::vec iacceptVec(nExp); iacceptVec.fill(iaccept);
  bool autoOn = impAutoEnabled();
  int isampleBudget = isample * nExp;   // AUTO reallocates within this total
  {
    // A per-subject isample= vector is given per BASE subject; tile it across
    // the mixture's expanded pseudo-subjects so component j of subject i gets
    // subject i's count.
    std::vector<int> isv;
    impNsampleVecGet(isv);
    if ((int)isv.size() > 1 && (int)isv.size() != nsub) {
      Rf_error("'isample' has length %d but there are %d subjects; a per-subject "
               "isample vector must have one entry per subject",
               (int)isv.size(), nsub);
    }
    if ((int)isv.size() == nsub) {
      for (int i = 0; i < nsub; ++i) {
        int v = isv[i] < 1 ? 1 : isv[i];
        for (int j = 0; j < Nmix; ++j) isampleVec[i + j * nsub] = v;
      }
    }
  }
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
  bool gammaInd = impGammaIndividual();
  // Largest relative per-step change the controller applied last iteration;
  // starts large so convergence cannot trip before the controller has run.
  double gammaStepPrev = R_PosInf;

  // Force the EM serial when the per-subject inner solves are not thread-safe:
  //  (a) mixtures -- the expanded pseudo-subjects' per-component solves race and
  //      the proportions/Omega drift and collapse;
  //  (b) multi-endpoint pool-sizing -- the inner MAP runs with ind->neqOverride
  //      against a pool sized for the larger theta-sensitivity model, whose
  //      per-thread work arrays are not safe under the override, so the parallel
  //      E-step non-deterministically rejects a subject's samples (neff collapse).
  // The plain single-endpoint path keeps full parallelism; restored before return.
  int savedCores = -1;
  if ((Nmix > 1 || impPoolSizing()) && cores > 1) { savedCores = impSetSolveCores(1); cores = 1; }

  // Initial MAP at the starting parameters.
  impMapPass(e);

  // NB placed AFTER the initial MAP pass on purpose: fInd->nObs is only filled
  // in during a solve, so reading it before impMapPass() returns 0 for every
  // subject and makes the whole dataset look sparse (which silently gave every
  // subject df = 4).
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
  const bool nonmemSparse = impAutoNonmemSparse();
  const int dfPatience = impAutoDfPatience();  // <= 0 disables withdrawal entirely
  if (autoOn) {
    bool nonNormal = impAutoNonNormal();
    for (int i = 0; i < nsub; ++i) {
      bool sparse = nonmemSparse && (impNobs(i) < neta);
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
      double dfI = (sparse || nonNormal) ? 30.0 : impDf();
      // TUNED: iaccept is left alone here.  Lowering it to 0.2 forces gamma
      // wide, and widening a Gaussian is the lever measured NOT to fix tails
      // while costing a lot of ESS -- on a Poisson fixture whose k-hat was
      // already -1.33 (no failure at all) it cut ESS 0.549 -> 0.412 and
      // doubled the objective noise for nothing.  It is now applied only where
      // k-hat says the proposal is genuinely struggling (below).
      double iaI = iaccept;
      for (int j = 0; j < Nmix; ++j) {
        dfVec[i + j * nsub] = dfI;
        dfFloor[i + j * nsub] = dfI;   // withdrawal may not go below this
        iacceptVec[i + j * nsub] = iaI;
      }
    }
  }


  // Fit-constant base seed for the pinned (qrRefresh=FALSE) Cranley-Patterson
  // shift streams; only consumed in that mode so every other path's draw
  // stream is unchanged.
  // Fit-constant salted base for the pinned (qrRefresh=FALSE) shift streams.
  uint32_t qrPinSeed = (uint32_t)impBaseSeed() + 0x9E3779B1u;
  // Salted base for the SIR stratified offsets; SIR only engages when it can
  // shrink the sample (sirN < isample) so every other path is unchanged.
  // SIR only engages where it actually shrinks that subject's sample.
  bool sir = impSirEnabled();
  int sirN = impSirN();
  uint32_t sirSeed = (uint32_t)impBaseSeed() + 0x7F4A7C15u;

  // Omega structure mask: only the elements estimated in the model (nonzero in
  // the starting Omega) are updated; the rest stay zero so the parameterization
  // keeps the same number of Omega thetas.
  arma::mat Om0;
  impGetOmega(Om0);
  arma::mat omMask = arma::conv_to<arma::mat>::from(Om0 != 0.0);
  // Eta indices whose Omega diagonal is fix()ed: their rows/columns are held at
  // the starting Omega through every EM update.
  std::vector<int> omFixedEta;
  impGetOmegaFixedEta(omFixedEta);

  std::vector<double> objTrace, gammaTrace, neffTrace, xiTrace;
  std::vector<arma::vec> parHist;
  bool converged = false;
  int iterRun = 0;

  // Iteration print + parameter-history capture (shared scale.h machinery).
  impIterPrintStart();

  arma::vec r(neta);
  for (int iter = 0; iter < nIter; ++iter) {
    if (iter > 0 && !isImp) impReMap();
    // Stash the E-step diagnostics on every iteration so the fit environment
    // reflects the last iteration actually run (the loop may stop early).
    // Global rule: every expanded subject shares the scalar scale.  Under
    // "individual" gammaVec already carries the per-subject scales the
    // controller set at the end of the previous iteration (and the constructor
    // fill for iteration 0), so it must NOT be flattened here.
    if (!gammaInd) gammaVec.fill(gamma);
    isampleUsed = isampleVec;
    dfUsed = dfVec;
    gammaUsed = gammaVec;
    gammaScalarUsed = gamma;
    impEStep(nsub, neta, isampleVec, gammaVec, dfVec, cores, iter, impLogDetOmegaInv5(), isImp,
             condMean, condVar, Li, Neff, Xi, XiExp, KhatExp, sampS, sampZk, aMat, &e, qrPinSeed);
    obj = 0.0;
    for (int id = 0; id < nsub; ++id) if (R_finite(Li[id])) obj += 2.0 * Li[id];
    // Mean effective-sample fraction (importance-sampling "acceptance ratio").
    double accFrac = 0.0; int nAcc = 0;
    for (int id = 0; id < nsub; ++id)
      if (R_finite(Neff[id])) { accFrac += Neff[id] / (double)isampleUsed[id]; ++nAcc; }
    if (nAcc > 0) accFrac /= (double)nAcc;
    // Mean NONMEM-style xi (a DIFFERENT statistic from accFrac -- see impEStep).
    double xiMean = 0.0; int nXi = 0;
    for (int id = 0; id < nsub; ++id)
      if (R_finite(Xi[id])) { xiMean += Xi[id]; ++nXi; }
    if (nXi > 0) xiMean /= (double)nXi; else xiMean = NA_REAL;
    // `gamma` is the scalar summary of the scale actually used this iteration:
    // itself under "global", the mean of the per-subject scales under
    // "individual".  It is what gets traced and reported as impGammaUsed.
    if (gammaInd) gamma = arma::mean(gammaVec);
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
    if (autoOn) {
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
      for (int id = 0; id < nExp; ++id) {
        double kh = KhatExp[id];
        if (!R_finite(kh)) continue;              // no usable k-hat: leave alone
        if (escDead[id]) continue;                // measured not to help here
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
        //   patience   sparse objRMSE   theo3 objRMSE   theo3 OmegaRMSE
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
        // Establish the baseline for a df that step 1 pre-assigned (the
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
          wantEsc = (dfVec[id] <= 0.0 || want < dfVec[id]);
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
        if (dfVec[id] == dfFloor[id]) kAtFloor[id] = kh;

        // Withdrawal.  Only for a subject actually moved off its floor -- testing
        // dfVec > 0 instead stranded subjects that were never escalated at all --
        // and only once the ladder is EXHAUSTED (!wantEsc), so "this is not
        // helping" means the heaviest applicable rung is not helping.  Note
        // "escalated" is dfVec != dfFloor, not an inequality: escalating a
        // Gaussian floor RAISES df (0 -> 30) while escalating a t floor LOWERS it
        // (30 -> 20).
        if (!nonmemSparse && dfPatience > 0 && !wantEsc &&
            dfVec[id] != dfFloor[id] && R_finite(kAtFloor[id])) {
          // Non-strict on the margin: a rung that buys exactly the 0.1 it is
          // asked for has earned its place.  Strict `<` withdrew a proposal that
          // took k-hat 1.25 -> 0.75 against a floor of 0.85.
          bool improved = (kh <= 0.7) || (kh <= kAtFloor[id] - 0.1);
          if (improved) {
            noImp[id] = 0;
          } else if (++noImp[id] >= dfPatience) {
            // back to the step-1 floor, NOT to 0: for a non-normal endpoint that
            // floor is the t proposal the model requires, and withdrawal is
            // permanent (escDead), so dropping below it could never be undone.
            dfVec[id] = dfFloor[id];
            escDead[id] = 1;
            continue;
          }
        }

        if (wantEsc) {
          // Each rung gets a full patience window; inheriting strikes withdrew a
          // heavier rung after a single iteration for want of an improvement it
          // had not had time to deliver.
          noImp[id] = 0;
          dfVec[id] = want;
        }
        // NOTE deliberately one-way on SUCCESS.  Relaxing on a low k-hat is
        // circular: once a t proposal is in place k-hat drops precisely BECAUSE
        // it is working (measured -1.3 to -1.9 on subjects that read 2.03 under a
        // Gaussian), so "safe now" is evidence the fix is needed, not that it is
        // not.  Relaxing on it walks straight back into the failure and
        // oscillates.  Efficiency comes from picking the right rung immediately
        // -- the severity mapping above -- not from backing off afterwards.
      }
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
    if (autoOn && nExp > 1) {
      arma::vec need(nExp, arma::fill::ones);
      for (int id = 0; id < nExp; ++id) {
        double fr = R_finite(Neff[id]) ? Neff[id] / (double)isampleUsed[id] : 1.0;
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
      if (tot > 0.0 && R_finite(tot)) {
        for (int id = 0; id < nExp; ++id) {
          int v = (int)std::lround((double)isampleBudget * need[id] / tot);
          if (v < 25) v = 25;                       // floor: keep PSIS usable
          if (v > 20 * isample) v = 20 * isample;   // ceiling: bound the cost
          isampleVec[id] = v;
        }
      }
    }
    objTrace.push_back(obj);
    gammaTrace.push_back(gamma);
    neffTrace.push_back(accFrac);
    xiTrace.push_back(xiMean);
    iterRun = iter + 1;

    // M-step.  First a Newton step on the non-mu structural thetas from the
    // IS-weighted score and Gauss-Newton Hessian accumulated over subjects/samples
    // -- done before the mu updates (which shift thetas/etas) so it sees the
    // E-step parameters, and before impSetEta since impThetaScore's re-solves
    // overwrite the etas.  Skipped if the Hessian is not usable (thetas unchanged).
    if (nSens > 0) {
      arma::vec g(nSens, arma::fill::zeros);
      arma::mat H(nSens, nSens, arma::fill::zeros);
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
      // outputs are stashed per subject, then accumulated in eid order below.  cores
      // already carries the mixture / pool-sizing serial guard (forced to 1 above).
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
      bool doParM = (cores > 1) && impMStepParallelOk();
      if (doParM) impInnerParallelOn();
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static) if(doParM)
#endif
      for (int eid = 0; eid < nExp; ++eid) {
#ifdef _OPENMP
        if (doParM) setRxThreadId(omp_get_thread_num());
#endif
        if (sampS[eid].n_rows != 0) {
          if (sir && sirN < (int)sampS[eid].n_rows) {
            // SIR acceleration: an equal-weight systematic resample stands in
            // for the full weighted sample, cutting the theta-sens solves from
            // isample to sirN per subject.  For a mixture, zk sums to the
            // component responsibility a_ij (folded in by impEStep), so the
            // equal weights carry a_ij/sirN to preserve the component mass.
            double aij = arma::accu(sampZk[eid]);
            if (aij > 0.0 && R_finite(aij)) {
              // Seed this subject's stream on the current thread; the drawn u0 is
              // a pure function of the (iter, eid) seed, so it is thread-count
              // invariant.
              nmSetSeedEng1(sirSeed + (uint32_t)((iter * nExp + eid) * 2));
              double u0 = rxUnifEng(0.0, 1.0);
              arma::uvec ridx = impSirIndex(sampZk[eid], sirN, u0);
              arma::mat Ssir = sampS[eid].rows(ridx);
              zksub[eid].set_size(sirN); zksub[eid].fill(aij / (double)sirN);
              impThetaSensCollect(eid, Ssir, coll[eid]);
              useSub[eid] = 1;
            }
          } else {
            zksub[eid] = sampZk[eid];
            impThetaSensCollect(eid, sampS[eid], coll[eid]);
            useSub[eid] = 1;
          }
        }
#ifdef _OPENMP
        if (doParM) setRxThreadId(-1);
#endif
      }
      if (doParM) impInnerParallelOff();
      // Serial accumulation in eid order -- bit-identical to the serial theta score.
      for (int eid = 0; eid < nExp; ++eid) {
        if (useSub[eid]) impThetaAccumOne(coll[eid], zksub[eid], g, H);
      }
      g /= (double)nsub;
      H /= (double)nsub;
      arma::vec step;
      if (arma::solve(step, H, g) && step.is_finite()) impUpdateStructThetas(step);
    }

    // Mixture: mean-posterior EM update of the $MIX proportions (the stable M-step
    // SAEM uses).  For constant proportions the fixed point of the NONMEM
    // Gauss-Newton score d_ij = a_{Nm,i}/a_Nm - a_ji/a_j is exactly a_j = mean_i
    // a_ij, and this update stays in the simplex (a full Newton step overshoots to
    // the boundary and collapses a component).  Install the target proportions as
    // absolute $MIX thetas via the multinomial logit theta_m = log(a_m / a_Nm),
    // floored away from 0 so a transiently-empty component can recover.  Uses
    // impmap's own responsibilities a_ij, separate from FOCEI's Laplace mixProb.
    if (Nmix > 1) {
      arma::vec aStar(Nmix, arma::fill::zeros);
      for (int i = 0; i < nsub; ++i)
        for (int j = 0; j < Nmix; ++j) aStar[j] += aMat(i, j);
      aStar /= (double)nsub;
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
    for (int id = 0; id < nsub; ++id) {
      arma::vec cm = condMean.row(id).t();
      impSetEta(id, cm);
    }
    impUpdateMuThetas();
    impMuInterceptStep();

    arma::mat Omega(neta, neta, arma::fill::zeros);
    for (int id = 0; id < nsub; ++id) {
      impGetEta(id, r);
      Omega += r * r.t() + condVar[id];
    }
    Omega /= (double)nsub;
    Omega %= omMask;
    // Restore fix()ed Omega rows/columns to their starting values.
    for (size_t k = 0; k < omFixedEta.size(); ++k) {
      int fi = omFixedEta[k];
      if (fi >= 0 && fi < neta) { Omega.row(fi) = Om0.row(fi); Omega.col(fi) = Om0.col(fi); }
    }
    impSetOmega(Omega, diagXform);

    // Record the current estimates for the parameter-stability half of the test.
    arma::vec parNow; impGetEstPar(parNow);
    parHist.push_back(parNow);
    // Print this iteration + record it into the parameter-history walk.
    impIterPrintRow(parNow, obj);

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
    if (nConvWindow > 0 && R_finite(obj) &&
        (int)objTrace.size() >= nConvWindow + 1) {
      int n = (int)objTrace.size();
      double objMetric;
      if (impGammaRuleTarget() && nConvWindow >= 2) {
        // TARGET rule: mean|delta obj| is a NOISE measure, not a drift measure --
        // it has a floor at the Monte-Carlo noise level and never reaches ctol no
        // matter how settled the fit is.  The target rule deliberately widens the
        // proposal to put xi on iaccept, which costs effective sample size and so
        // RAISES that floor: measured on the 3-eta theophylline fixture,
        // mean|delta| is 0.00055-0.00093 under "floor" (ESS 0.95-0.96, all below
        // a ctol of 1e-3) against 0.00181-0.00264 under "target" (ESS 0.69-0.71,
        // all above it), so the fit could never converge however long it ran.
        // Use the drift of the window MEAN instead -- the same noise-robust form
        // the gamma gate below uses.  Scoped to this rule so the floor rule's
        // tuned behaviour is bit-identical.
        int hw = nConvWindow / 2;
        double mNew = 0.0, mOld = 0.0;
        for (int k = n - hw; k < n; ++k) mNew += objTrace[k];
        for (int k = n - nConvWindow; k < n - hw; ++k) mOld += objTrace[k];
        mNew /= (double)hw;
        mOld /= (double)(nConvWindow - hw);
        objMetric = std::fabs(mNew - mOld) / std::max(1.0, std::fabs(obj));
      } else {
        double s = 0.0;
        for (int k = n - nConvWindow; k < n; ++k) s += std::fabs(objTrace[k] - objTrace[k - 1]);
        objMetric = (s / (double)nConvWindow) / std::max(1.0, std::fabs(obj));
      }
      const arma::vec& parOld = parHist[n - nConvWindow - 1];
      double parMetric = arma::max(arma::abs(parNow - parOld) / (arma::abs(parOld) + 1e-8));
      // A settled parameter still has ~1-2% Monte-Carlo net drift across the
      // window, so this gate sits above that floor: it is a safety against
      // premature stopping on a systematically-drifting low-leverage parameter,
      // not a precision requirement (the objMetric gate carries the precision).
      double parTol = 0.02;
      double gWin0 = gammaTrace[n - nConvWindow - 1];
      // Under "individual" the traced scalar is a MEAN, so offsetting moves
      // across subjects could look settled while individual scales are still
      // swinging.  Gate on the largest per-subject step the controller actually
      // applied instead, which cannot be masked that way.
      // gammaRule="target" tracks a MONTE-CARLO statistic, so gamma keeps jittering
      // around its fixed point and never stops moving the way the one-sided floor
      // rule's gamma does (which is literally constant on a healthy model, making
      // the window test below vacuous).  Gate it on the drift of the window MEAN
      // instead: noise averages out over nConvWindow, a genuine trend does not.
      bool gammaStable;
      if (gammaInd) {
        gammaStable = (gammaStepPrev <= 1e-3);
      } else if (impGammaRuleTarget() && nConvWindow >= 2) {
        // Split the window in half and compare the two means.  BOTH halves must be
        // non-empty: at nConvWindow=1 the older half has no elements, mOld stays 0,
        // and the test degenerates to gamma <= gammaTargetTol -- unsatisfiable,
        // since gamma is bounded below by iscaleMin (0.1 by default), so the fit
        // could never converge.  Hence the nConvWindow >= 2 guard here and the
        // fall-through to the instantaneous test below.
        int hw = nConvWindow / 2;              // >= 1, and nConvWindow - hw >= 1
        double mNew = 0.0, mOld = 0.0;
        for (int k = n - hw; k < n; ++k) mNew += gammaTrace[k];
        for (int k = n - nConvWindow; k < n - hw; ++k) mOld += gammaTrace[k];
        mNew /= (double)hw;
        mOld /= (double)(nConvWindow - hw);
        gammaStable = (std::fabs(mNew - mOld) <= gammaTargetTol * std::max(1.0, mOld));
      } else {
        // Instantaneous test.  The floor rule's gamma is exactly constant once
        // coverage is healthy so 1e-3 is ample; a target-tracking gamma carries
        // xi's Monte-Carlo noise, so it gets the same looser tolerance the
        // window-mean test uses (without it, nConvWindow=1 could never converge).
        double tol = impGammaRuleTarget() ? gammaTargetTol : 1e-3;
        gammaStable = (std::fabs(gamma - gWin0) <= tol * std::max(1.0, gWin0));
      }
      if (gammaStable && objMetric < ctol && parMetric < parTol) { converged = true; break; }
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
    // healthy from iteration 1 onward -- and the "global" branch below is
    // one-sided, so it could never come back down.
    if (iter == 0) {
      // nothing to adapt from yet
    } else if (gammaInd) {
      // xi_i falls monotonically as the proposal
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
      const double pExp = 2.0 / std::max(1.0, (double)neta);
      const double capUp = 1.25, capDn = 1.0 / 1.25;
      double maxStep = 0.0;
      for (int id = 0; id < nExp; ++id) {
        double xiId = XiExp[id];
        if (ISNAN(xiId)) continue;      // no usable xi (dead subject): leave alone
        double fac;
        if (!R_finite(xiId)) {
          fac = capUp;                  // xi = +Inf: as under-dispersed as it gets
        } else if (xiId <= 0.0) {
          fac = capDn;                  // xi = 0: proposal far too wide
        } else {
          fac = std::pow(xiId / iacceptVec[id], pExp);
          if (fac > capUp) fac = capUp;
          else if (fac < capDn) fac = capDn;
        }
        double g = gammaVec[id] * fac;
        if (g < iscaleMin) g = iscaleMin;
        if (g > iscaleMax) g = iscaleMax;
        // Measure the step actually taken (after clamping) so a subject pinned
        // at a bound reads as settled rather than as perpetually moving.
        double applied = (gammaVec[id] > 0.0) ? std::fabs(g / gammaVec[id] - 1.0) : 0.0;
        if (applied > maxStep) maxStep = applied;
        gammaVec[id] = g;
      }
      gammaStepPrev = maxStep;
    } else if (impGammaRuleTarget()) {
      // gammaRule="target" -- NONMEM's rule for the SHARED scale.  The NM7
      // Technical Guide (note after eq. 1.76) says gamma is "continually adjusted
      // so that xi_i approximates IACCEPT", i.e. two-sided on xi, not one-sided on
      // the Kish fraction.  Reuse the per-subject branch's analytic inversion:
      // xi = gamma^(-neta/2) inverts exactly, so pExp = 2/neta gives an error
      // multiplier of 0 at every dimension.  A FIXED exponent (the sqrt the
      // "floor" branch below uses) is the form measured to enter a period-2 limit
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
      // (it cannot diverge) and the window-mean gate below will simply not report
      // convergence.  gammaMethod="individual" has no such exposure -- each
      // subject inverts its own xi -- so that is the remedy when it bites.
      if (iaccept > 0 && R_finite(xiMean) && xiMean > 0.0) {
        const double pExp = 2.0 / std::max(1.0, (double)neta);
        double fac = std::pow(xiMean / iaccept, pExp);
        if (fac > 1.25) fac = 1.25;
        else if (fac < 1.0 / 1.25) fac = 1.0 / 1.25;
        double g = gamma * fac;
        if (g < iscaleMin) g = iscaleMin;
        if (g > iscaleMax) g = iscaleMax;
        // measure the step AFTER clamping, so a gamma pinned at a bound reads as
        // settled rather than as perpetually moving (mirrors gammaStepPrev in the
        // per-subject branch); the convergence gate compares gamma across the
        // window, which a clamped gamma satisfies for free.
        gamma = g;
      }
    } else if (iaccept > 0 && accFrac > 0 && accFrac < iaccept) {
      // gammaRule="floor" (default) -- iaccept is a one-sided FLOOR on the mean
      // Kish effective-sample fraction: gamma is left at its efficient starting
      // value while coverage is healthy and inflated only when it drops below.
      double fac = std::sqrt(iaccept / accFrac);
      if (fac > 1.25) fac = 1.25;
      gamma *= fac;
      if (gamma < iscaleMin) gamma = iscaleMin;
      if (gamma > iscaleMax) gamma = iscaleMax;
    }
  }

  // Close the iteration print and stash the parameter-history walk (e$parHistData).
  impIterPrintGet(e);

  // Finalize the fit at the converged estimates.  Drop the objective stashed by
  // the initial MAP pass first: the finalize path keeps an existing e$objective
  // (re-adjusting it as if it were an unadjusted -2LL), which would publish the
  // initial-parameter objective instead of the converged FOCEi evaluation.
  if (e.exists("objective")) e.remove("objective");
  impSyncInitParToFullTheta();
  impMapPass(e);

  // Monte-Carlo observed-information covariance (theta block) at the converged
  // estimate, before the inner neqOverride is cleared (the inner solves use it).
  // Selected by impmapControl(covMethod="imp") (the default).
  // Snapshots, not the live vectors: the adaptation steps have already moved
  // gammaVec/dfVec on for an iteration that never ran, so passing them would
  // evaluate the covariance at a proposal the fit never actually used (and
  // disagree with the reported impDfInd / impGammaInd).
  if (impCovEnabled()) impComputeCov(e, gammaUsed, dfUsed);

  // Clear the multi-endpoint inner neqOverride so it does not leak into a later fit.
  impClearInnerNeqOverride();

  // Stash the last-iteration E-step diagnostics.
  List condVarList(nsub);
  for (int id = 0; id < nsub; ++id) condVarList[id] = wrap(condVar[id]);
  e["impCondMean"] = wrap(condMean);
  e["impCondVar"]  = condVarList;
  e["impLi"]       = wrap(Li);
  e["impNeff"]     = wrap(Neff);
  e["impXi"]       = wrap(Xi);
  // Per-expanded-subject proposal scales actually used, and which controller
  // produced them.  Under "global" every element equals impGammaUsed.
  e["impGammaInd"] = wrap(gammaUsed);
  e["impGammaMethod"] = gammaInd ? "individual" : "global";
  e["impObj"]      = obj;
  // Snapshot, so the scalar cannot disagree with the impGammaInd vector: under
  // gammaMethod="global" the scalar gamma is adapted at the END of the loop for
  // an iteration that never runs, which would report gamma_start*fac here
  // against gamma_start in the vector.
  e["impGammaUsed"] = gammaScalarUsed;
  e["impNsample"]  = isample;
  e["impNsampleInd"] = wrap(isampleUsed);  // counts USED by the last E-step
  e["impDfInd"]    = wrap(dfUsed);          // df USED by the last E-step
  e["impIacceptInd"] = wrap(iacceptVec);    // per-expanded-subject acceptance target
  e["impKhatIter"] = wrap(KhatExp);         // last-iteration in-kernel k-hat
  e["impAuto"]     = autoOn;
  e["impQr"]       = impQrEnabled();
  e["impSir"]      = impSirEnabled();
  e["impSirSample"] = impSirN();
  e["impNiter"]    = nIter;
  e["impIter"]     = iterRun;
  e["impConverged"] = converged;
  e["impObjTrace"] = wrap(objTrace);
  e["impGammaTrace"] = wrap(gammaTrace);
  e["impNeffFrac"] = wrap(neffTrace);
  e["impXiTrace"]  = wrap(xiTrace);
  e["impMuGroupN"] = impMuGroupN();

  // Restore the solve core count if it was forced serial for the mixture EM.
  if (savedCores > 0) impSetSolveCores(savedCores);
}
