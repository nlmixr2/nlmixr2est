// NPB -- Nonparametric Bayes: a truncated stick-breaking Dirichlet-process
// mixture sampled by a blocked Metropolis-within-Gibbs sampler (Tatarinova 2013,
// Ishwaran & James truncation).  Shares the conditional-likelihood primitive
// (npEvalCondLik / npBuildPsiCore) and the fit finalization (npFinalizeFit) with
// NPAG; it differs only in how the support points and weights are drawn.  Called
// from foceiFitCpp_ in place of foceiOuter when est=="npb".
//
// Model: F(theta) = sum_{k=1..K} w_k delta_{phi_k}, phi_k ~ G_0, v_k ~ Beta(1,alpha),
// w_1 = v_1, w_k = v_k prod_{j<k}(1-v_j), truncated at K (v_K = 1).  Base measure
// G_0 = N(0, diag(Omega_prior)) in eta space.  Each Gibbs sweep:
//   (a) cluster assignments  z_i ~ Categorical(w_k p(y_i | phi_k))
//   (b) stick weights        v_k | counts ~ Beta(1 + n_k, alpha + sum_{j>k} n_j)
//   (c) support locations    phi_k via a Gaussian random-walk MH against G_0 and
//                            the assigned subjects' conditional likelihoods
// Post-burn-in draws give the posterior mixing distribution, per-subject
// posterior-mean etas, and posterior draws of the population mean (credible
// intervals).  The reported objective is the nonparametric -2LL, NOT comparable
// to the NONMEM/FOCEI OFV.
#include <RcppArmadillo.h>
#include "np.h"
#include "npCommon.h"
#include "imp.h"

using namespace Rcpp;

// Gelman-Rubin potential-scale-reduction factor (R-hat) per parameter across
// chains.  Each chain matrix is (draws x params); returns one R-hat per param
// (~1 at convergence).  Returns all-ones for < 2 chains or < 2 draws.
static arma::vec npbGelmanRubin(const std::vector<arma::mat>& chains) {
  int m = (int)chains.size();
  int p = (m > 0) ? (int)chains[0].n_cols : 0;
  arma::vec rhat(p, arma::fill::ones);
  if (m < 2) return rhat;
  int n = (int)chains[0].n_rows;
  if (n < 2) return rhat;
  for (int j = 0; j < p; ++j) {
    arma::vec cm(m), cv(m);
    for (int c = 0; c < m; ++c) {
      arma::vec col = chains[c].col(j);
      cm[c] = arma::mean(col);
      cv[c] = arma::var(col);                 // (n-1) denominator
    }
    double grand = arma::mean(cm);
    double B = (double)n / (m - 1) * arma::accu(arma::square(cm - grand));
    double W = arma::mean(cv);
    if (W <= 0.0) { rhat[j] = 1.0; continue; }
    double Vhat = ((double)(n - 1) / n) * W + (1.0 / n) * B;
    rhat[j] = std::sqrt(Vhat / W);
  }
  return rhat;
}

void npbOuter(Environment e) {
  List control = e["control"];
  int K = as<int>(control["npPoints"]);          // truncation level
  int nBurn = as<int>(control["npBurnin"]);
  int nSamp = as<int>(control["npNsamp"]);
  double alpha = as<double>(control["npAlpha"]);
  double propSd = as<double>(control["npPropSd"]);
  int seed = as<int>(control["npSeed"]);
  int nchains = control.containsElementNamed("nchains") ? as<int>(control["nchains"]) : 1;
  if (nchains < 1) nchains = 1;
  int cores = impCores();

  int nsub = impNsub();
  int neta = impNeta();

  // likInner0 forms the Omega-prior term; rebuild omegaInv from the initial Omega
  // (see npagOuter) and use its diagonal as the base measure G_0 variance.
  arma::vec g0var(neta, arma::fill::ones);
  arma::mat omModel;
  impGetOmega(omModel);
  if ((int)omModel.n_rows == neta && neta > 0) {
    impSetOmega(omModel, impDiagXform());
    g0var = arma::clamp(arma::vec(omModel.diag()), 1e-6, arma::datum::inf);
  }
  arma::vec g0sd = arma::sqrt(g0var);

  // npb does not support mix() models (guarded in R); no mixture-proportion setup.
  Function setSeed("set.seed");

  // pooled across chains
  arma::mat postEtaAcc(nsub, neta, arma::fill::zeros);
  arma::mat meanDraws((size_t)nchains * nSamp, neta, arma::fill::zeros);
  std::vector<arma::rowvec> supAll;
  std::vector<double> wAll;
  std::vector<arma::mat> chainMeanDraws(nchains);   // per chain, for Gelman-Rubin
  arma::mat lastPsi;             // final-sweep Psi of the last chain (objective)
  arma::vec wLast;              // final weights of the last chain
  int total = nBurn + nSamp;

  // Independent chains (seed offset per chain) for Gelman-Rubin R-hat.
  for (int chain = 0; chain < nchains; ++chain) {
  // reproducible RNG seeded from the control (via R's RNG so set.seed also works)
  setSeed(seed + chain);
  GetRNGstate();
  // initialize support points from G_0, uniform weights
  arma::mat phi(K, neta);
  for (int k = 0; k < K; ++k)
    for (int j = 0; j < neta; ++j) phi(k, j) = R::rnorm(0.0, g0sd[j]);
  arma::vec w(K); w.fill(1.0 / K);
  std::vector<int> z(nsub, 0);
  arma::mat chainMd(nSamp, neta, arma::fill::zeros);
  int drawn = 0;

  for (int it = 0; it < total; ++it) {
    // (a) conditional likelihoods p(y_i | phi_k) and cluster assignments
    arma::mat psi;
    npBuildPsiCore(phi, cores, psi);   // nsub x K
    lastPsi = psi;
    std::vector<int> nk(K, 0);
    for (int i = 0; i < nsub; ++i) {
      arma::rowvec p = psi.row(i) % w.t();
      double s = arma::accu(p);
      int zi = 0;
      if (s > 0.0) {
        double u = R::unif_rand() * s, c = 0.0;
        for (int k = 0; k < K; ++k) { c += p[k]; zi = k; if (u <= c) break; }
      } else {
        zi = (int)(R::unif_rand() * K);
        if (zi >= K) zi = K - 1;
      }
      z[i] = zi; nk[zi]++;
    }
    // (b) stick-breaking weights  v_k ~ Beta(1 + n_k, alpha + sum_{j>k} n_j)
    std::vector<int> tail(K, 0);
    { int acc = 0; for (int k = K - 1; k >= 0; --k) { tail[k] = acc; acc += nk[k]; } }
    double cumLog1mV = 0.0;
    for (int k = 0; k < K; ++k) {
      double vv = (k == K - 1) ? 1.0 : R::rbeta(1.0 + nk[k], alpha + (double)tail[k]);
      w[k] = std::exp(cumLog1mV) * vv;
      cumLog1mV += std::log(std::max(1e-300, 1.0 - vv));
    }
    // (c) support locations: MH for occupied clusters, fresh G_0 draw otherwise
    for (int k = 0; k < K; ++k) {
      if (nk[k] == 0) {
        for (int j = 0; j < neta; ++j) phi(k, j) = R::rnorm(0.0, g0sd[j]);
        continue;
      }
      arma::rowvec cur = phi.row(k), prop = cur;
      for (int j = 0; j < neta; ++j) prop[j] += R::rnorm(0.0, propSd);
      double curLp = 0.0, propLp = 0.0;
      for (int j = 0; j < neta; ++j) {
        curLp += -0.5 * cur[j] * cur[j] / g0var[j];
        propLp += -0.5 * prop[j] * prop[j] / g0var[j];
      }
      std::vector<double> curEta(cur.begin(), cur.end()), propEta(prop.begin(), prop.end());
      for (int i = 0; i < nsub; ++i) {
        if (z[i] != k) continue;
        curLp += npEvalCondLik(&curEta[0], i);
        propLp += npEvalCondLik(&propEta[0], i);
      }
      if (std::log(R::unif_rand()) < propLp - curLp) phi.row(k) = prop;
    }
    // (d) accumulate post-burn-in draws
    if (it >= nBurn) {
      for (int i = 0; i < nsub; ++i) postEtaAcc.row(i) += phi.row(z[i]);
      arma::rowvec md = w.t() * phi;
      chainMd.row(drawn) = md;
      meanDraws.row((size_t)chain * nSamp + drawn) = md;
      for (int k = 0; k < K; ++k) if (w[k] > 1e-8) { supAll.push_back(phi.row(k)); wAll.push_back(w[k]); }
      drawn++;
    }
  }
  PutRNGstate();
  chainMeanDraws[chain] = chainMd;
  wLast = w;
  } // end chain loop

  arma::vec rhat = npbGelmanRubin(chainMeanDraws);
  arma::mat postEta = postEtaAcc / ((double)nchains * nSamp);

  // posterior mixing distribution E[F]: pooled post-burn-in support points with
  // their (renormalized) weights.
  int M = (int)supAll.size();
  arma::mat support(std::max(M, 1), neta, arma::fill::zeros);
  arma::vec weights(std::max(M, 1), arma::fill::value(1.0));
  if (M > 0) {
    double wsum = 0.0;
    for (int m = 0; m < M; ++m) { support.row(m) = supAll[m]; weights[m] = wAll[m]; wsum += wAll[m]; }
    if (wsum > 0.0) weights /= wsum;
  }

  // nonparametric marginal log-likelihood at the final draw (last chain)
  double objf = 0.0;
  for (int i = 0; i < nsub; ++i) {
    double s = arma::accu(lastPsi.row(i) % wLast.t());
    objf += std::log(std::max(1e-300, s));
  }

  // npb does not optimize the assay-error multiplier, so there is no gamma to
  // fold into the residual thetas (gamma = 1 makes the fold a no-op).
  // npb does not mu-expand (guarded in R), so no injected etas to collapse.
  arma::mat Omega = npFinalizeFit(e, support, weights, postEta, objf, omModel,
                                  1.0, std::vector<int>(),
                                  std::vector<int>(), std::vector<int>());

  e["npbSupport"] = wrap(support);          // pooled posterior support (E[F])
  e["npbWeights"] = wrap(weights);
  e["npbPosteriorEta"] = wrap(postEta);     // per-subject posterior mean eta
  e["npbMeanDraws"] = wrap(meanDraws);      // pooled posterior draws of the mean (CIs)
  e["npbRhat"] = wrap(rhat);                // Gelman-Rubin R-hat per eta (~1 converged)
  e["npbNchains"] = nchains;
  e["npbOmega"] = wrap(Omega);
  e["npbAlpha"] = alpha;
  e["npbK"] = K;
  e["npbBurnin"] = nBurn;
  e["npbNsamp"] = nSamp;
  e["npbLogLik"] = objf;
}
