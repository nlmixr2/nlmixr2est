// Importance-sampling EM (IMPMAP / IMP) kernel.
//
// Peer of saem.cpp / inner.cpp holding the NONMEM-style Monte Carlo importance-
// sampling EM.  It reuses the FOCEI inner machinery (mode + Hessian + joint-
// density evaluation) through the thin numeric interface in imp.h, so all
// model-evaluation heavy lifting stays in inner.cpp.
//
// Module M1: single MAP pass collecting each subject's mode/Hessian/likelihood.
// Module M2: per-subject multivariate-normal proposal + thread-safe threefry
// sampler.  The importance weights and EM updates are added in later modules.
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

  arma::mat modeMat(nsub, neta, arma::fill::zeros);
  arma::vec indLik(nsub, arma::fill::zeros);
  List hessList(nsub);

  // Per-subject mode and the lower Cholesky factor of the proposal covariance
  // Sigma_i = gamma * H_i^-1 (Sigma_i = L_i L_i'); the linear algebra is done
  // serially, before the RNG loop.
  std::vector<arma::vec> modes(nsub);
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
      if (arma::inv_sympd(Sigma, H)) {
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

  // Draw isample proposal samples per subject: phi_k = mode_i + L_i * z, z ~ N(0,I).
  // Each subject reseeds its thread's threefry stream to seed0 + id*2 -- offset by
  // 2 to decorrelate from the ODE solver's per-subject stream (seed0 + id) -- AFTER
  // setRxThreadId sets the OpenMP thread id, so a subject's draws depend only on
  // (seed0, id) and are identical regardless of the thread count.
  seedEng(cores);
  uint32_t seed0 = getRxSeed1(cores);
  std::vector<arma::mat> samples(nsub);
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
    if (haveL[id]) {
      for (int k = 0; k < isample; ++k) {
        arma::vec z(neta);
        for (int j = 0; j < neta; ++j) z[j] = rxNormEng(0.0, 1.0);
        S.row(k) = (modes[id] + cholL[id] * z).t();
      }
    }
    samples[id] = S;
#ifdef _OPENMP
    if (doPar) setRxThreadId(-1);
#endif
  }

  List samplesList(nsub);
  for (int id = 0; id < nsub; ++id) samplesList[id] = wrap(samples[id]);

  e["impEtaMode"]   = wrap(modeMat);
  e["impIndLik"]    = wrap(indLik);
  e["impEtaHess"]   = hessList;
  e["impSamples"]   = samplesList;
  e["impGammaUsed"] = gamma;
  e["impNsample"]   = isample;
}
