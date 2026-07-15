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

void npbOuter(Environment e) {
  List control = e["control"];
  int K = as<int>(control["npPoints"]);          // truncation level
  int nBurn = as<int>(control["npBurnin"]);
  int nSamp = as<int>(control["npNsamp"]);
  double alpha = as<double>(control["npAlpha"]);
  double propSd = as<double>(control["npPropSd"]);
  int seed = as<int>(control["npSeed"]);
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

  // reproducible RNG seeded from the control (via R's RNG so set.seed also works)
  Function setSeed("set.seed");
  setSeed(seed);
  GetRNGstate();

  // initialize support points from G_0, uniform weights
  arma::mat phi(K, neta);
  for (int k = 0; k < K; ++k)
    for (int j = 0; j < neta; ++j) phi(k, j) = R::rnorm(0.0, g0sd[j]);
  arma::vec w(K); w.fill(1.0 / K);

  arma::mat postEtaAcc(nsub, neta, arma::fill::zeros);
  arma::mat meanDraws(nSamp, neta, arma::fill::zeros);
  std::vector<arma::rowvec> supAll;
  std::vector<double> wAll;
  std::vector<int> z(nsub, 0);
  arma::mat lastPsi;              // final-sweep Psi for the objective
  int drawn = 0;
  int total = nBurn + nSamp;

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
      meanDraws.row(drawn) = w.t() * phi;
      for (int k = 0; k < K; ++k) if (w[k] > 1e-8) { supAll.push_back(phi.row(k)); wAll.push_back(w[k]); }
      drawn++;
    }
  }
  PutRNGstate();

  arma::mat postEta = postEtaAcc / (double)nSamp;

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

  // nonparametric marginal log-likelihood at the final draw
  double objf = 0.0;
  for (int i = 0; i < nsub; ++i) {
    double s = arma::accu(lastPsi.row(i) % w.t());
    objf += std::log(std::max(1e-300, s));
  }

  // npb does not optimize the assay-error multiplier, so there is no gamma to
  // fold into the residual thetas (gamma = 1 makes the fold a no-op).
  arma::mat Omega = npFinalizeFit(e, support, weights, postEta, objf, omModel,
                                  1.0, std::vector<int>());

  e["npbSupport"] = wrap(support);          // pooled posterior support (E[F])
  e["npbWeights"] = wrap(weights);
  e["npbPosteriorEta"] = wrap(postEta);     // per-subject posterior mean eta
  e["npbMeanDraws"] = wrap(meanDraws);      // posterior draws of the population mean (CIs)
  e["npbOmega"] = wrap(Omega);
  e["npbAlpha"] = alpha;
  e["npbK"] = K;
  e["npbBurnin"] = nBurn;
  e["npbNsamp"] = nSamp;
  e["npbLogLik"] = objf;
}
