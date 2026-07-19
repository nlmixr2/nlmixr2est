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
#include "nmMcmcRng.h"
#include "imp.h"
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
static void impEStep(int nsub, int neta, int isample, double gamma, int cores,
                     int iter, double negHalfLogDetOmega, bool isImp,
                     arma::mat& condMean, std::vector<arma::mat>& condVar,
                     arma::vec& Li, arma::vec& Neff,
                     std::vector<arma::mat>& outS, std::vector<arma::vec>& outZk,
                     arma::mat& aMat, Environment* eStash,
                     uint32_t qrPinSeed) {
  double invGamma2 = 1.0 / (2.0 * gamma);
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
        if (arma::chol(L, gamma * Sig, "lower")) {
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
          Sigma *= gamma;
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
    qrU0 = impSobolU0(isample, neta);
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
      arma::mat S(isample, neta);
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
        for (int k = 0; k < isample; ++k) {
          S.row(k) = (modes[id] + cholL[id] * Z.row(k).t()).t();
        }
      } else {
        for (int k = 0; k < isample; ++k) {
          arma::vec z(neta);
          for (int jj = 0; jj < neta; ++jj) z[jj] = rxNormEng(0.0, 1.0);
          S.row(k) = (modes[id] + cholL[id] * z).t();
        }
      }
      // log importance weights.  impEvalJointLik returns the NEGATIVE log joint
      // density (minimized by the inner optimizer), so it enters with a minus;
      // the proposal density is subtracted, contributing +(1/(2 gamma)) d'H d.
      arma::vec q(isample);
      int nGood = 0;
      for (int k = 0; k < isample; ++k) {
        arma::vec eta = S.row(k).t();
        arma::vec d = eta - modes[id];
        double qk = -impEvalJointLik(eta, id) +
          invGamma2 * arma::as_scalar(d.t() * Hs[id] * d);
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
      double logMeanExp = qmax + std::log(sumw / (double)isample);
      double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gamma) -
        0.5 * logDetH[id];
      LiExp[id] = -(logMeanExp + Ci);          // per-component negative log-likelihood
      NeffExp[id] = 1.0 / arma::accu(arma::square(zk));
      arma::rowvec pbar = zk.t() * S;
      cmExp.row(id) = pbar;
      arma::mat B(neta, neta, arma::fill::zeros);
      for (int k = 0; k < isample; ++k) {
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
  aMat.set_size(nsub, Nm); aMat.zeros();
  for (int i = 0; i < nsub; ++i) {
    if (Nm == 1) {
      condMean.row(i) = cmExp.row(i);
      condVar[i] = cvExp[i];
      Li[i] = LiExp[i];
      Neff[i] = NeffExp[i];
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
    // Fold the responsibility a_ij into each component's importance weights so the
    // M-step, iterating the expanded subjects, forms the component-weighted score.
    for (int j = 0; j < Nm; ++j) {
      int eid = i + j * nsub;
      if (outZk[eid].n_elem > 0) outZk[eid] *= aMat(i, j);
    }
  }

  // On the last iteration, stash the per-subject sampler diagnostics (first
  // component for a mixture).
  if (eStash != nullptr) {
    arma::mat modeMat(nsub, neta);
    List hessList(nsub), samplesList(nsub);
    for (int id = 0; id < nsub; ++id) {
      modeMat.row(id) = modes[id].t();
      hessList[id] = haveL[id] ? wrap(Hs[id]) : R_NilValue;
      samplesList[id] = wrap(outS[id]);
    }
    (*eStash)["impEtaMode"] = wrap(modeMat);
    (*eStash)["impEtaHess"] = hessList;
    (*eStash)["impSamples"] = samplesList;
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
void impComputeCov(Environment e) {
  int nsub = impNsub();
  int neta = impNeta();
  int isample = impNsample();
  double gamma = impGammaProp();
  double invGamma2 = 1.0 / (2.0 * gamma);
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
        Sigma *= gamma;
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
        S.row(k) = (modes[id] + Ls[id] * Z.row(k).t()).t();
      }
    } else {
      for (int k = 0; k < isample; ++k) {
        arma::vec z(neta);
        for (int j = 0; j < neta; ++j) z[j] = rxNormEng(0.0, 1.0);
        S.row(k) = (modes[id] + Ls[id] * z).t();
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
        arma::vec q(isample); int nGood = 0;
        for (int k = 0; k < isample; ++k) {
          arma::vec eta = Ss[id].row(k).t();
          arma::vec d = eta - modes[id];
          double qk = -impEvalJointLik(eta, id) +
            invGamma2 * arma::as_scalar(d.t() * Hs[id] * d);
          if (R_finite(qk)) { q[k] = qk; ++nGood; } else q[k] = R_NegInf;
        }
        if (nGood != 0) {
          double qmax = q.max();
          double sumw = arma::accu(arma::exp(q - qmax));
          double logMeanExp = qmax + std::log(sumw / (double)isample);
          double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gamma) - 0.5 * logDetH[id];
          objBuf[id] = -2.0 * (logMeanExp + Ci);
        }
      }
#ifdef _OPENMP
      if (doParCov) setRxThreadId(-1);
#endif
    }
    double obj = 0.0;
    for (int id = 0; id < nsub; ++id) obj += objBuf[id];
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
  double iscaleMin = impIscaleMin();
  double iscaleMax = impIscaleMax();
  int nConvWindow = impNconvWindow();
  double ctol = impCtol();

  arma::mat condMean;
  std::vector<arma::mat> condVar;
  std::vector<arma::mat> sampS;
  std::vector<arma::vec> sampZk;
  arma::vec Li, Neff;
  arma::mat aMat;                 // posterior mixture responsibilities (nsub x Nmix)
  int Nmix = impNmix();
  int nExp = nsub * Nmix;         // expanded pseudo-subjects for the mixture E/M-step
  double obj = R_PosInf;
  int nSens = impThetaSensN();
  // est="imp": skip the per-iteration MAP search; the E-step proposal is centered
  // at the running conditional mean with covariance gamma*Omega.
  bool isImp = impIsImp();

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

  // Fit-constant base seed for the pinned (qrRefresh=FALSE) Cranley-Patterson
  // shift streams; only consumed in that mode so every other path's draw
  // stream is unchanged.
  // Fit-constant salted base for the pinned (qrRefresh=FALSE) shift streams.
  uint32_t qrPinSeed = (uint32_t)impBaseSeed() + 0x9E3779B1u;
  // Salted base for the SIR stratified offsets; SIR only engages when it can
  // shrink the sample (sirN < isample) so every other path is unchanged.
  bool sir = impSirEnabled() && impSirN() < isample;
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

  std::vector<double> objTrace, gammaTrace, neffTrace;
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
    impEStep(nsub, neta, isample, gamma, cores, iter, impLogDetOmegaInv5(), isImp,
             condMean, condVar, Li, Neff, sampS, sampZk, aMat, &e, qrPinSeed);
    obj = 0.0;
    for (int id = 0; id < nsub; ++id) if (R_finite(Li[id])) obj += 2.0 * Li[id];
    // Mean effective-sample fraction (importance-sampling "acceptance ratio").
    double accFrac = 0.0; int nAcc = 0;
    for (int id = 0; id < nsub; ++id)
      if (R_finite(Neff[id])) { accFrac += Neff[id] / (double)isample; ++nAcc; }
    if (nAcc > 0) accFrac /= (double)nAcc;
    objTrace.push_back(obj);
    gammaTrace.push_back(gamma);
    neffTrace.push_back(accFrac);
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
      double s = 0.0;
      for (int k = n - nConvWindow; k < n; ++k) s += std::fabs(objTrace[k] - objTrace[k - 1]);
      double objMetric = (s / (double)nConvWindow) / std::max(1.0, std::fabs(obj));
      const arma::vec& parOld = parHist[n - nConvWindow - 1];
      double parMetric = arma::max(arma::abs(parNow - parOld) / (arma::abs(parOld) + 1e-8));
      // A settled parameter still has ~1-2% Monte-Carlo net drift across the
      // window, so this gate sits above that floor: it is a safety against
      // premature stopping on a systematically-drifting low-leverage parameter,
      // not a precision requirement (the objMetric gate carries the precision).
      double parTol = 0.02;
      double gWin0 = gammaTrace[n - nConvWindow - 1];
      bool gammaStable = std::fabs(gamma - gWin0) <= 1e-3 * std::max(1.0, gWin0);
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
    if (iaccept > 0 && accFrac > 0 && accFrac < iaccept) {
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
  if (impCovEnabled()) impComputeCov(e);

  // Clear the multi-endpoint inner neqOverride so it does not leak into a later fit.
  impClearInnerNeqOverride();

  // Stash the last-iteration E-step diagnostics.
  List condVarList(nsub);
  for (int id = 0; id < nsub; ++id) condVarList[id] = wrap(condVar[id]);
  e["impCondMean"] = wrap(condMean);
  e["impCondVar"]  = condVarList;
  e["impLi"]       = wrap(Li);
  e["impNeff"]     = wrap(Neff);
  e["impObj"]      = obj;
  e["impGammaUsed"] = gamma;
  e["impNsample"]  = isample;
  e["impQr"]       = impQrEnabled();
  e["impSir"]      = impSirEnabled();
  e["impSirSample"] = impSirN();
  e["impNiter"]    = nIter;
  e["impIter"]     = iterRun;
  e["impConverged"] = converged;
  e["impObjTrace"] = wrap(objTrace);
  e["impGammaTrace"] = wrap(gammaTrace);
  e["impNeffFrac"] = wrap(neffTrace);
  e["impMuGroupN"] = impMuGroupN();

  // Restore the solve core count if it was forced serial for the mixture EM.
  if (savedCores > 0) impSetSolveCores(savedCores);
}
