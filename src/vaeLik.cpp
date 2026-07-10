// vaeLik.cpp -- FOCE-linearized marginal -2LL for the VAE fit object, with
// M2/M3/M4 censoring handled exactly like the FOCEi inner problem (inner.cpp)
// via censEst.h, parallelized over subjects with OpenMP.
//
// The decoder gives, at the encoder linearization point zLin_i (the individual
// transformed parameters), the prediction f0, the eta-Jacobian J0 = d f/d eta,
// and the residual variance R for every observation. The linearized model is
//   f(z) = f0 + J0 (z - zLin_i).
// For each subject a Newton step set drives z to the FOCE mode of the
// linearized conditional posterior (MAP-EBE) -- minimizing
//   -ll_data(z) + 0.5 (z - zPop_i)' Omega^-1 (z - zPop_i)
// where ll_data sums per-observation normal (or censored-normal) log-likelihoods
// (censored obs use doCensNormal1 for the value and censNormalPartials for the
// score/curvature). The Laplace marginal -2LL is
//   -2 ll_data(z*) + (z*-zPop)' Omega^-1 (z*-zPop) + log|Omega| + log|H|,
// H = Omega^-1 + sum_j rho_ff_j J0_j J0_j'. No ODE re-solves inside the Hessian,
// so the covariance (numerical Hessian of this objective over the population
// parameters) is cheap. Additive/proportional scale enters through Rscale.
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
                int adjLik, int maxit, double tol, int cores) {
  const int N = zLin.n_rows;
  const int neta = zLin.n_cols;
  const double ln2pi = std::log(2.0 * M_PI);
  const arma::vec oInv = 1.0 / omega;
  const double logDetOmega = arma::accu(arma::log(omega));
  const arma::mat OmegaInvD = arma::diagmat(oInv);

  arma::vec obji(N, arma::fill::zeros);
  arma::mat zStar(N, neta, arma::fill::zeros);
  const bool doParallel = (cores > 1);

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
    const arma::vec zL = zLin.row(i).t();
    const arma::vec zP = zPop.row(i).t();
    arma::vec z = zL;
    arma::mat H(neta, neta);

    for (int it = 0; it < maxit; ++it) {
      arma::vec grad = OmegaInvD * (z - zP);
      H = OmegaInvD;
      for (int j = o0; j < o1; ++j) {
        const arma::rowvec Jj = J0.row(j);
        const double fj = f0[j] + arma::dot(Jj, (z - zL));
        const double Rj = R[j] * Rscale;
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

    // Laplace -2LL at the mode
    double ll = 0.0;
    H = OmegaInvD;
    for (int j = o0; j < o1; ++j) {
      const arma::rowvec Jj = J0.row(j);
      const double fj = f0[j] + arma::dot(Jj, (z - zL));
      const double Rj = R[j] * Rscale;
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
    obji[i] = -2.0 * ll + arma::dot(ez, oInv % ez) + logDetOmega + logDetH;
    zStar.row(i) = z.t();
  }

  return List::create(_["objective"] = arma::accu(obji), _["obji"] = obji,
                      _["zStar"] = zStar);
}
