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
// E-step conditional mean and variance.  The M-step (Omega / theta updates) is
// added later.
#include <RcppArmadillo.h>
#include <rxode2ptr.h>
#include "imp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

void impOuter(Environment e) {
  // Single MAP pass (reuses the FOCEI posthoc path); leaves each subject's
  // conditional mode in the live inner state.
  impMapPass(e);

  int nsub = impNsub();
  int neta = impNeta();
  int isample = impNsample();
  double gamma = impGammaProp();
  int cores = impCores();
  if (cores < 1) cores = 1;

  double negHalfLogDetOmega = impLogDetOmegaInv5(); // = -0.5 * log|Omega|

  arma::mat modeMat(nsub, neta, arma::fill::zeros);
  arma::vec indLik(nsub, arma::fill::zeros);
  List hessList(nsub);

  // Per-subject mode, inner information matrix H_i, log|H_i|, and the lower
  // Cholesky factor of the proposal covariance Sigma_i = gamma * H_i^-1
  // (Sigma_i = L_i L_i').  The linear algebra is done serially, before the RNG.
  std::vector<arma::vec> modes(nsub);
  std::vector<arma::mat> Hs(nsub);
  std::vector<double> logDetH(nsub, 0.0);
  std::vector<arma::mat> cholL(nsub);
  std::vector<char> haveL(nsub, 0);
  arma::vec mode(neta);
  arma::mat H(neta, neta, arma::fill::zeros);
  for (int id = 0; id < nsub; ++id) {
    impGetMode(id, mode);
    modes[id] = mode;
    modeMat.row(id) = mode.t();
    indLik[id] = impGetIndLik(id);
    if (impGetHessian(id, H)) {
      hessList[id] = wrap(H);
      arma::mat Sigma;
      double ldv, lds;
      if (arma::inv_sympd(Sigma, H) && arma::log_det(ldv, lds, H) && lds > 0) {
        Hs[id] = H;
        logDetH[id] = ldv;
        Sigma *= gamma;
        arma::mat L;
        if (arma::chol(L, Sigma, "lower")) {
          cholL[id] = L;
          haveL[id] = 1;
        }
      }
    } else {
      hessList[id] = R_NilValue;
    }
  }

  // For each subject: draw isample proposal samples phi_k = mode_i + L_i * z,
  // z ~ N(0,I); then form importance weights w_k = exp(q_k) with
  //   q_k = log(l(y_i|phi_k,theta) * h(phi_k|Omega)) [joint density kernel]
  //         + (1/(2 gamma)) (phi_k - mode_i)' H_i (phi_k - mode_i) [minus log proposal].
  // Normalized weights z_k = w_k / sum_k w_k give the conditional mean phi_bar_i
  // and variance B_i; L_i = -log( mean_k exp(q_k + C_i) ) is the individual
  // objective contribution, C_i = -0.5 log|Omega| + 0.5 n log gamma - 0.5 log|H_i|.
  //
  // Each subject reseeds its thread's threefry stream to seed0 + id*2 (offset by
  // 2 to decorrelate from the ODE solver's per-subject stream seed0 + id) AFTER
  // setRxThreadId sets the OpenMP thread id, so a subject's draws depend only on
  // (seed0, id) and are identical regardless of the thread count.
  seedEng(cores);
  uint32_t seed0 = getRxSeed1(cores);
  double invGamma2 = 1.0 / (2.0 * gamma);
  std::vector<arma::mat> samples(nsub);
  std::vector<arma::mat> condVar(nsub);
  arma::mat condMean(nsub, neta, arma::fill::zeros);
  arma::vec Li(nsub); Li.fill(NA_REAL);
  arma::vec Neff(nsub); Neff.fill(NA_REAL);
  bool doPar = (cores > 1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doPar)
#endif
  for (int id = 0; id < nsub; ++id) {
#ifdef _OPENMP
    if (doPar) setRxThreadId(omp_get_thread_num());
#endif
    setSeedEng1(seed0 + (uint32_t)id * 2u);
    arma::mat S(isample, neta, arma::fill::zeros);
    condVar[id] = arma::mat(neta, neta, arma::fill::zeros);
    condMean.row(id) = modes[id].t();
    if (haveL[id]) {
      for (int k = 0; k < isample; ++k) {
        arma::vec z(neta);
        for (int j = 0; j < neta; ++j) z[j] = rxNormEng(0.0, 1.0);
        S.row(k) = (modes[id] + cholL[id] * z).t();
      }
      // importance weights on the log scale.  impEvalJointLik returns the
      // NEGATIVE log joint density (the inner optimizer minimizes it), so the
      // joint-density term enters with a minus sign; the proposal density is
      // subtracted, contributing +(1/(2 gamma)) d' H d.
      arma::vec q(isample);
      for (int k = 0; k < isample; ++k) {
        arma::vec eta = S.row(k).t();
        arma::vec d = eta - modes[id];
        q[k] = -impEvalJointLik(eta, id) +
          invGamma2 * arma::as_scalar(d.t() * Hs[id] * d);
      }
      double qmax = q.max();
      arma::vec w = arma::exp(q - qmax);
      double sumw = arma::accu(w);
      arma::vec zk = w / sumw;
      double logMeanExp = qmax + std::log(sumw / (double)isample);
      double Ci = negHalfLogDetOmega + 0.5 * neta * std::log(gamma) -
        0.5 * logDetH[id];
      Li[id] = -(logMeanExp + Ci);
      Neff[id] = 1.0 / arma::accu(arma::square(zk));
      arma::rowvec pbar = zk.t() * S; // weighted conditional mean (1 x neta)
      condMean.row(id) = pbar;
      arma::mat B(neta, neta, arma::fill::zeros);
      for (int k = 0; k < isample; ++k) {
        arma::rowvec dr = S.row(k) - pbar;
        B += zk[k] * (dr.t() * dr);
      }
      condVar[id] = B;
    }
    samples[id] = S;
#ifdef _OPENMP
    if (doPar) setRxThreadId(-1);
#endif
  }

  List samplesList(nsub), condVarList(nsub);
  for (int id = 0; id < nsub; ++id) {
    samplesList[id] = wrap(samples[id]);
    condVarList[id] = wrap(condVar[id]);
  }
  // objective (-2 log likelihood) summed over subjects with a valid proposal
  double obj = 0.0;
  for (int id = 0; id < nsub; ++id) {
    if (R_finite(Li[id])) obj += 2.0 * Li[id];
  }

  e["impEtaMode"]   = wrap(modeMat);
  e["impIndLik"]    = wrap(indLik);
  e["impEtaHess"]   = hessList;
  e["impSamples"]   = samplesList;
  e["impCondMean"]  = wrap(condMean);
  e["impCondVar"]   = condVarList;
  e["impLi"]        = wrap(Li);
  e["impNeff"]      = wrap(Neff);
  e["impObj"]       = obj;
  e["impGammaUsed"] = gamma;
  e["impNsample"]   = isample;
}
