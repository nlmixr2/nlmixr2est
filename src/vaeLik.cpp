// vaeLik.cpp -- FOCE-linearized marginal -2LL for the VAE fit object, with
// M2/M3/M4 censoring handled exactly like the FOCEi inner problem (inner.cpp)
// via censEst.h, MIXTURE models combined like inner.cpp, parallelized over
// subjects with OpenMP.
//
// The decoder gives, at the encoder linearization point zLin_i (the individual
// transformed parameters), the prediction f0, the eta-Jacobian J0 = d f/d eta,
// and the residual variance R for every observation. The linearized model is
//   f(z) = f0 + J0 (z - zLin_i).
// For each subject (and each mixture component) a Newton step set drives z to
// the FOCE mode of the linearized conditional posterior (MAP-EBE) -- minimizing
//   -ll_data(z) + 0.5 (z - zPop_i)' Omega^-1 (z - zPop_i)
// where ll_data sums per-observation normal (or censored-normal) log-likelihoods
// (censored obs use doCensNormal1 for the value and censNormalPartials for the
// score/curvature). The Laplace marginal log-likelihood per component m is
//   ll_i^(m) = -0.5 [ -2 ll_data(z*) + (z*-zPop)' Omega^-1 (z*-zPop)
//                     + log|Omega| + log|H| ].
// Mixtures combine as in inner.cpp: L_i = sum_m mixProb[m] exp(ll_i^(m)), so the
// individual -2LL is -2 logsumexp_m( log mixProb[m] + ll_i^(m) ); the reported
// EBE (and mixest) is the argmax component's mode. nMix=1 (mixProb=1) reduces to
// the plain FOCE-linearized -2LL.
//
// Component-major layout: f0/J0/R and zLin/zPop are stacked over components
// (component m's block starts at m*totalObs obs and m*N subjects); y/cens/limit
// are shared across components (the data is the same, only the model differs).
#include "armahead.h"
#include "censEst.h"
#include "rxomp.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

static inline bool vaeIsObs(int cens, double lim) {
  return cens == 0 && !(R_FINITE(lim) && !ISNA(lim));
}

// [[Rcpp::export]]
List vaeFoceLik(arma::vec f0, arma::mat J0, arma::vec R, arma::vec y,
                arma::ivec cens, arma::vec limit, arma::mat zLin, arma::mat zPop,
                arma::vec omega, arma::ivec subjOffset, double Rscale,
                int adjLik, int maxit, double tol, int cores,
                int nMix, arma::vec mixProb) {
  const int N = zLin.n_rows / nMix;
  const int neta = zLin.n_cols;
  const int totalObs = subjOffset[subjOffset.n_elem - 1];
  const double ln2pi = std::log(2.0 * M_PI);
  const arma::vec oInv = 1.0 / omega;
  const double logDetOmega = arma::accu(arma::log(omega));
  const arma::mat OmegaInvD = arma::diagmat(oInv);
  const bool doParallel = (cores > 1);

  arma::vec obji(N, arma::fill::zeros);
  // one MAP eta per (component, subject): component-major, row m*N+i, matching
  // the fit etaMat (mix-1 etas for all subjects, then mix-2, ...).
  arma::mat zStar(N * nMix, neta, arma::fill::zeros);
  arma::ivec mixest(N, arma::fill::ones);

  // Laplace -2LL for one subject under one mixture component (its data block
  // starts at obs offset `base`); returns the objective and the mode in `zOut`.
  auto laplace = [&](int base, int o0, int o1, const arma::vec& zL,
                     const arma::vec& zP, arma::vec& zOut) -> double {
    arma::vec z = zL;
    arma::mat H(neta, neta);
    for (int it = 0; it < maxit; ++it) {
      arma::vec grad = OmegaInvD * (z - zP);
      H = OmegaInvD;
      for (int j = o0; j < o1; ++j) {
        const arma::rowvec Jj = J0.row(base + j);
        const double fj = f0[base + j] + arma::dot(Jj, (z - zL));
        const double Rj = R[base + j] * Rscale;
        double rf, rff;
        if (vaeIsObs(cens[j], limit[j])) {
          rf = -(y[j] - fj) / Rj; rff = 1.0 / Rj;
        } else {
          double cp[3];
          censNormalPartials((double)cens[j], y[j], limit[j], fj, Rj, 2, cp);
          rf = cp[0]; rff = cp[2];
        }
        grad += rf * Jj.t();
        H += rff * (Jj.t() * Jj);
      }
      arma::vec step;
      if (!arma::solve(step, H, grad)) break;
      z -= step;
      if (arma::norm(step) < tol) break;
    }
    double ll = 0.0;
    H = OmegaInvD;
    for (int j = o0; j < o1; ++j) {
      const arma::rowvec Jj = J0.row(base + j);
      const double fj = f0[base + j] + arma::dot(Jj, (z - zL));
      const double Rj = R[base + j] * Rscale;
      const double res = y[j] - fj;
      const double llg = -0.5 * (res * res / Rj + std::log(Rj) + ln2pi);
      ll += doCensNormal1((double)cens[j], y[j], limit[j], llg, fj, Rj, adjLik);
      double rff;
      if (vaeIsObs(cens[j], limit[j])) {
        rff = 1.0 / Rj;
      } else {
        double cp[3];
        censNormalPartials((double)cens[j], y[j], limit[j], fj, Rj, 2, cp);
        rff = cp[2];
      }
      H += rff * (Jj.t() * Jj);
    }
    const arma::vec ez = z - zP;
    double logDetH, sgn;
    arma::log_det(logDetH, sgn, H);
    zOut = z;
    return -2.0 * ll + arma::dot(ez, oInv % ez) + logDetOmega + logDetH;
  };

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int i = 0; i < N; ++i) {
#ifdef _OPENMP
    // Cross-DLL OpenMP thread-id fix (see inner.cpp): rxode2 and nlmixr2est each
    // statically link their own libgomp on Windows, so hand rxode2 our real
    // thread id inside the parallel region.
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    const int o0 = subjOffset[i], o1 = subjOffset[i + 1];
    // Per-subject work is wrapped in try/catch: an exception escaping the OpenMP
    // region calls std::terminate (see inner.cpp). A per-subject numerical
    // failure (singular Hessian, non-finite solve) falls back to a large
    // objective and the encoder mode.
    try {
      // per-component Laplace log-likelihood (mixture insertion): each component
      // gets its own MAP eta (as if nMix*N pseudo-subjects), stored per component.
      arma::vec compLL(nMix);
      int bestM = 0;
      double bestLL = -std::numeric_limits<double>::infinity();
      for (int m = 0; m < nMix; ++m) {
        arma::vec zMode;
        const double objM = laplace(m * totalObs, o0, o1,
                                    zLin.row(m * N + i).t(), zPop.row(m * N + i).t(), zMode);
        compLL[m] = std::log(mixProb[m]) - 0.5 * objM;
        zStar.row(m * N + i) = zMode.t();
        if (compLL[m] > bestLL) { bestLL = compLL[m]; bestM = m; }
      }
      // -2 logsumexp_m( log mixProb[m] + ll_i^(m) )
      double s = 0.0;
      for (int m = 0; m < nMix; ++m) s += std::exp(compLL[m] - bestLL);
      double val = -2.0 * (bestLL + std::log(s));
      obji[i] = R_FINITE(val) ? val : 1e30;
      mixest[i] = bestM + 1;
    } catch (...) {
      obji[i] = 1e30;
      for (int m = 0; m < nMix; ++m) zStar.row(m * N + i) = zLin.row(m * N + i);
      mixest[i] = 1;
    }
  }

  return List::create(_["objective"] = arma::accu(obji), _["obji"] = obji,
                      _["zStar"] = zStar, _["mixest"] = mixest);
}
