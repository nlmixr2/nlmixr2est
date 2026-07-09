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
#include "imp.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// One importance-sampling E-step at the current conditional modes: draw
// proposal samples, form importance weights, and return each subject's
// conditional mean, variance, objective contribution, and effective sample
// size.  `iter` shifts the per-subject RNG stream so successive iterations use
// fresh (still thread-count-independent) samples.
static void impEStep(int nsub, int neta, int isample, double gamma, int cores,
                     int iter, double negHalfLogDetOmega,
                     arma::mat& condMean, std::vector<arma::mat>& condVar,
                     arma::vec& Li, arma::vec& Neff, Environment* eStash) {
  double invGamma2 = 1.0 / (2.0 * gamma);

  // Per-subject mode, information matrix H_i, log|H_i|, and the lower Cholesky
  // factor of the proposal covariance gamma * H_i^-1 -- computed serially.
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

  condMean.set_size(nsub, neta);
  condVar.assign(nsub, arma::mat(neta, neta, arma::fill::zeros));
  Li.set_size(nsub); Li.fill(NA_REAL);
  Neff.set_size(nsub); Neff.fill(NA_REAL);

  std::vector<arma::mat> samples(nsub);
  seedEng(cores);
  uint32_t seed0 = getRxSeed1(cores);
  bool doPar = (cores > 1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doPar)
#endif
  for (int id = 0; id < nsub; ++id) {
#ifdef _OPENMP
    if (doPar) setRxThreadId(omp_get_thread_num());
#endif
    // fresh per-(iter,subject) stream, independent of thread count
    setSeedEng1(seed0 + (uint32_t)((iter * nsub + id) * 2));
    condMean.row(id) = modes[id].t();
    if (haveL[id]) {
      arma::mat S(isample, neta);
      for (int k = 0; k < isample; ++k) {
        arma::vec z(neta);
        for (int j = 0; j < neta; ++j) z[j] = rxNormEng(0.0, 1.0);
        S.row(k) = (modes[id] + cholL[id] * z).t();
      }
      // log importance weights.  impEvalJointLik returns the NEGATIVE log joint
      // density (minimized by the inner optimizer), so it enters with a minus;
      // the proposal density is subtracted, contributing +(1/(2 gamma)) d'H d.
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
      arma::rowvec pbar = zk.t() * S;
      condMean.row(id) = pbar;
      arma::mat B(neta, neta, arma::fill::zeros);
      for (int k = 0; k < isample; ++k) {
        arma::rowvec dr = S.row(k) - pbar;
        B += zk[k] * (dr.t() * dr);
      }
      condVar[id] = B;
      if (eStash != nullptr) samples[id] = S;
    }
#ifdef _OPENMP
    if (doPar) setRxThreadId(-1);
#endif
  }

  // On the last iteration, stash the per-subject sampler diagnostics.
  if (eStash != nullptr) {
    arma::mat modeMat(nsub, neta);
    List hessList(nsub), samplesList(nsub);
    for (int id = 0; id < nsub; ++id) {
      modeMat.row(id) = modes[id].t();
      hessList[id] = haveL[id] ? wrap(Hs[id]) : R_NilValue;
      samplesList[id] = wrap(samples[id]);
    }
    (*eStash)["impEtaMode"] = wrap(modeMat);
    (*eStash)["impEtaHess"] = hessList;
    (*eStash)["impSamples"] = samplesList;
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

  arma::mat condMean;
  std::vector<arma::mat> condVar;
  arma::vec Li, Neff;
  double obj = R_PosInf;

  // Initial MAP at the starting parameters.
  impMapPass(e);

  // Omega structure mask: only the elements estimated in the model (nonzero in
  // the starting Omega) are updated; the rest stay zero so the parameterization
  // keeps the same number of Omega thetas.
  arma::mat Om0;
  impGetOmega(Om0);
  arma::mat omMask = arma::conv_to<arma::mat>::from(Om0 != 0.0);

  arma::vec r(neta);
  for (int iter = 0; iter < nIter; ++iter) {
    if (iter > 0) impReMap();
    Environment* stash = (iter == nIter - 1) ? &e : nullptr;
    impEStep(nsub, neta, isample, gamma, cores, iter, impLogDetOmegaInv5(),
             condMean, condVar, Li, Neff, stash);
    obj = 0.0;
    for (int id = 0; id < nsub; ++id) if (R_finite(Li[id])) obj += 2.0 * Li[id];

    // M-step.  Seed each subject's eta with its conditional mean, then update
    // the mu-referenced population parameters -- covariate groups by regression
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
    impSetOmega(Omega, diagXform);
  }

  // Finalize the fit at the converged estimates.
  impSyncInitParToFullTheta();
  impMapPass(e);

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
  e["impNiter"]    = nIter;
  e["impMuGroupN"] = impMuGroupN();
}
