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

// Occupied support (weight > tiny) as a compact (support, weights) pair for the
// residual/regressor optimization, falling back to the max-weight cluster.  The
// truncated stick keeps K clusters but most carry ~0 weight; passing only the
// occupied ones keeps the frozen-solve cache small.
static void npbCompactSupport(const arma::mat& phi, const arma::vec& w,
                              arma::mat& sup, arma::vec& wt) {
  std::vector<arma::uword> occ;
  for (arma::uword k = 0; k < w.n_elem; ++k) if (w[k] > 1e-8) occ.push_back(k);
  if (occ.empty()) occ.push_back(w.index_max());
  arma::uvec idx(occ);
  sup = phi.rows(idx);
  wt = w(idx);
  double s = arma::accu(wt);
  if (s > 0.0) wt /= s;
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

  // mixture (mix()) support: populate mixProb from the proportion thetas so
  // npBuildPsiCore marginalizes over the components (npMixCondLik).  When the
  // proportions are estimated (npMixOptimize), sample them each sweep from
  // Dirichlet(alpha0 + counts); otherwise hold the ini proportions.
  impUpdateMixProbs();
  bool mixOpt = control.containsElementNamed("npMixOptimize") &&
    as<bool>(control["npMixOptimize"]);
  int nMix = impNmix();
  double mixAlpha0 = 1.0;   // symmetric Dirichlet prior on the proportions

  // residual/regressor optimization (shared with npag): with the mixing
  // distribution held fixed, bounded-bobyqa the residual-error thetas (and any
  // non-mu structural "regressor" thetas) against the nonparametric -2LL.  The
  // index/kind/bound vectors are built in R (.npFamilyFit), identical to npag.
  // npResidMode: 0 none, 1 alternate (re-fit during burn-in, then hold fixed for
  // the sampling phase so every collected draw shares the converged residuals),
  // 2 final (fit once at the converged draw).
  int residMode = control.containsElementNamed("npResidMode") ? as<int>(control["npResidMode"]) : 0;
  std::vector<int> residOptIdx, residOptKind, residOptEnd, residOptProp;
  std::vector<double> residOptLower, residOptUpper;
  std::vector<int> regressIdx;
  std::vector<double> regressLower, regressUpper;
  arma::ivec obsEndpoint;
  bool residFreeze = true;
  if (control.containsElementNamed("npResidOptIdx")) {
    IntegerVector ro = control["npResidOptIdx"];
    residOptIdx.assign(ro.begin(), ro.end());
  }
  if (control.containsElementNamed("npResidOptKind")) {
    IntegerVector rk = control["npResidOptKind"];
    residOptKind.assign(rk.begin(), rk.end());
  }
  if (control.containsElementNamed("npResidOptLower")) {
    NumericVector rl = control["npResidOptLower"];
    residOptLower.assign(rl.begin(), rl.end());
  }
  if (control.containsElementNamed("npResidOptUpper")) {
    NumericVector ru = control["npResidOptUpper"];
    residOptUpper.assign(ru.begin(), ru.end());
  }
  if (control.containsElementNamed("npResidOptEnd")) {
    IntegerVector re = control["npResidOptEnd"];
    residOptEnd.assign(re.begin(), re.end());
  }
  if (control.containsElementNamed("npResidOptProp")) {
    IntegerVector rp = control["npResidOptProp"];
    residOptProp.assign(rp.begin(), rp.end());
  }
  if (control.containsElementNamed("npEndpointCmt")) {
    IntegerVector ec = control["npEndpointCmt"];
    std::vector<int> endpointCmt(ec.begin(), ec.end());
    obsEndpoint = npBuildObsEndpoint(endpointCmt);
  }
  if (control.containsElementNamed("npRegressIdx")) {
    IntegerVector gi = control["npRegressIdx"];
    regressIdx.assign(gi.begin(), gi.end());
  }
  if (control.containsElementNamed("npRegressLower")) {
    NumericVector gl = control["npRegressLower"];
    regressLower.assign(gl.begin(), gl.end());
  }
  if (control.containsElementNamed("npRegressUpper")) {
    NumericVector gu = control["npRegressUpper"];
    regressUpper.assign(gu.begin(), gu.end());
  }
  if (control.containsElementNamed("npResidFreeze"))
    residFreeze = as<bool>(control["npResidFreeze"]);
  bool hasResidOpt = residMode != 0 && (!residOptIdx.empty() || !regressIdx.empty());
  bool useRegress = !regressIdx.empty();
  // combined optimized set: residual (err) thetas, then the regressor thetas (kind 0,
  // no endpoint).  A regressor triggers re-derivation of the posterior-mean etas.
  std::vector<int> optIdx = residOptIdx, optKind = residOptKind;
  std::vector<int> optEnd = residOptEnd, optProp = residOptProp;
  std::vector<double> optLo = residOptLower, optHi = residOptUpper;
  for (size_t g = 0; g < regressIdx.size(); ++g) {
    optIdx.push_back(regressIdx[g]);
    optKind.push_back(0);
    optEnd.push_back(-1);
    optProp.push_back(0);
    optLo.push_back(g < regressLower.size() ? regressLower[g] : R_NegInf);
    optHi.push_back(g < regressUpper.size() ? regressUpper[g] : R_PosInf);
  }
  // bound the in-burn-in optimization rounds (each is a bounded bobyqa with a Psi
  // rebuild) so a long burn-in does not multiply the cost.
  int residEvery = std::max(1, nBurn / 10);
  Function setSeed("set.seed");

  // pooled across chains
  std::vector<double> mixProbSum(std::max(nMix, 1), 0.0);   // post-burn-in proportion draws
  int mixProbN = 0;
  arma::mat postEtaAcc(nsub, neta, arma::fill::zeros);
  arma::mat meanDraws((size_t)nchains * nSamp, neta, arma::fill::zeros);
  std::vector<arma::rowvec> supAll;
  std::vector<double> wAll;
  std::vector<arma::mat> chainMeanDraws(nchains);   // per chain, for Gelman-Rubin
  arma::mat phiLast;            // final-sweep support of the last chain (objective)
  arma::vec wLast;              // final weights of the last chain
  int total = nBurn + nSamp;

  // shared scale.h iteration printer + parameter history (chain 0 sweeps).
  impIterPrintStart();

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
    // per-sweep parameter-history row (chain 0), reusing the just-built psi for
    // the marginal -2LL of the incoming (phi, w) draw.  Pushing the support
    // covariance into op_focei is safe (npBuildPsiCore sums llikObs and does not
    // use omegaInv) and is overwritten by npFinalizeFit.  impSetOmega calls back
    // into R (rxSymInvCholCreate), so -- exactly as the npOptimizeResid step does
    // -- bracket it with Put/GetRNGstate to preserve the sampler's RNG stream
    // (and thus reproducibility) across the R call.
    if (chain == 0) {
      double llit = 0.0;
      for (int i = 0; i < nsub; ++i)
        llit += std::log(std::max(1e-300, arma::accu(psi.row(i) % w.t())));
      arma::rowvec meanEta = w.t() * phi;
      arma::mat Om(neta, neta, arma::fill::zeros);
      for (int k = 0; k < K; ++k) {
        arma::rowvec d = phi.row(k) - meanEta; Om += w[k] * (d.t() * d);
      }
      PutRNGstate();
      impSetOmega(npMaskedOmega(Om, omModel), impDiagXform());
      GetRNGstate();
      arma::vec par; impGetEstPar(par);
      impIterPrintRow(par, -2.0 * llit);
    }
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
    // (a2) mixture proportions: Dirichlet Gibbs draw from each subject's assigned
    // support eta and its posterior component responsibility.
    if (mixOpt && nMix > 1) {
      arma::mat subEta(nsub, neta);
      for (int i = 0; i < nsub; ++i) subEta.row(i) = phi.row(z[i]);
      npbSampleMixProbs(subEta, mixAlpha0);
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
    // (c) support locations: MH for occupied clusters, fresh G_0 draw otherwise.
    // The proposal + accept/reject draws stay serial and in their original order
    // (so the RNG stream, and thus reproducibility, is unchanged); only the
    // per-subject conditional-likelihood solves are parallelized
    // (npbSupportMHContrib).  Serial pre-pass draws the proposals and priors:
    std::vector<char> occ(K, 0);
    std::vector<std::vector<double> > curLoc(K), propLoc(K);
    std::vector<double> curPrior(K, 0.0), propPrior(K, 0.0), acceptU(K, 0.0);
    for (int k = 0; k < K; ++k) {
      if (nk[k] == 0) {
        for (int j = 0; j < neta; ++j) phi(k, j) = R::rnorm(0.0, g0sd[j]);
        continue;
      }
      occ[k] = 1;
      arma::rowvec cur = phi.row(k), prop = cur;
      for (int j = 0; j < neta; ++j) prop[j] += R::rnorm(0.0, propSd);
      acceptU[k] = R::unif_rand();  // drawn now (before the solves) to keep RNG order
      double cP = 0.0, pP = 0.0;
      for (int j = 0; j < neta; ++j) {
        cP += -0.5 * cur[j] * cur[j] / g0var[j];
        pP += -0.5 * prop[j] * prop[j] / g0var[j];
      }
      curPrior[k] = cP; propPrior[k] = pP;
      curLoc[k].assign(cur.begin(), cur.end());
      propLoc[k].assign(prop.begin(), prop.end());
    }
    // parallel solves of each subject's cur/prop conditional log-likelihood
    std::vector<double> curContrib, propContrib;
    npbSupportMHContrib(z, occ, curLoc, propLoc, curContrib, propContrib);
    // serial reduction (in ascending subject order per cluster, matching the old
    // accumulation exactly) + accept/reject
    std::vector<double> curLp(K), propLp(K);
    for (int k = 0; k < K; ++k) { curLp[k] = curPrior[k]; propLp[k] = propPrior[k]; }
    for (int i = 0; i < nsub; ++i) {
      int k = z[i];
      if (k >= 0 && k < K && occ[k]) { curLp[k] += curContrib[i]; propLp[k] += propContrib[i]; }
    }
    for (int k = 0; k < K; ++k) {
      if (!occ[k]) continue;
      if (std::log(acceptU[k]) < propLp[k] - curLp[k])
        for (int j = 0; j < neta; ++j) phi(k, j) = propLoc[k][j];
    }
    // (c2) residual/regressor optimization (alternate): re-fit the residual thetas
    // against the current draw's mixing distribution during burn-in, then hold them
    // for the sampling phase so every collected draw shares the converged residual
    // scale.  Chain 0 only; later chains inherit the fitted residuals.  Bracketed by
    // Put/GetRNGstate because npOptimizeResid calls back into R (bobyqa).
    if (hasResidOpt && residMode == 1 && chain == 0 && it > 0 && it < nBurn &&
        (it % residEvery == 0)) {
      arma::mat rsup; arma::vec rwt;
      npbCompactSupport(phi, w, rsup, rwt);
      PutRNGstate();
      npOptimizeResid(rsup, rwt, optIdx, optKind, cores, optLo, optHi, residFreeze,
                      obsEndpoint, optEnd, optProp, useRegress);
      GetRNGstate();
    }
    // (d) accumulate post-burn-in draws
    if (it >= nBurn) {
      if (mixOpt && nMix > 1) {
        for (int m = 0; m < nMix; ++m) mixProbSum[m] += impMixProb(m);
        mixProbN++;
      }
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
  phiLast = phi;
  } // end chain loop

  // residual/regressor optimization (final): fit once at the converged draw.
  // (alternate mode has already converged the residuals during burn-in.)
  if (hasResidOpt && residMode == 2) {
    arma::mat rsup; arma::vec rwt;
    npbCompactSupport(phiLast, wLast, rsup, rwt);
    npOptimizeResid(rsup, rwt, optIdx, optKind, cores, optLo, optHi, residFreeze,
                    obsEndpoint, optEnd, optProp, useRegress);
  }

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

  // nonparametric marginal log-likelihood at the final draw (last chain), rebuilt
  // AFTER any residual optimization so it reflects the installed residual thetas.
  arma::mat finalPsi;
  npBuildPsiCore(phiLast, cores, finalPsi);
  double objf = 0.0;
  for (int i = 0; i < nsub; ++i) {
    double s = arma::accu(finalPsi.row(i) % wLast.t());
    objf += std::log(std::max(1e-300, s));
  }

  // install the posterior-MEAN mixture proportions (more stable than the last
  // draw) so the reported proportion thetas reflect the posterior, and stash the
  // mean for accessors.
  if (mixOpt && nMix > 1 && mixProbN > 0) {
    arma::vec pmean(nMix);
    for (int m = 0; m < nMix; ++m) pmean[m] = mixProbSum[m] / (double)mixProbN;
    for (int m = 0; m < nMix; ++m) pmean[m] = std::max(1e-4, pmean[m]);
    pmean /= arma::accu(pmean);
    Environment rxode2ns = Environment::namespace_env("rxode2");
    Function mlogit = as<Function>(rxode2ns["mlogit"]);
    NumericVector pin(nMix - 1);
    for (int m = 0; m < nMix - 1; ++m) pin[m] = pmean[m];
    NumericVector th = mlogit(pin);
    arma::vec thv(nMix - 1);
    for (int m = 0; m < nMix - 1; ++m) thv[m] = th[m];
    impSetMixThetas(thv);
    e["npbMixProb"] = wrap(pmean);
  }
  // npb does not optimize the assay-error multiplier, so there is no gamma to
  // fold into the residual thetas (gamma = 1 makes the fold a no-op).
  // npb does not mu-expand (guarded in R), so no injected etas to collapse.
  arma::mat Omega = npFinalizeFit(e, support, weights, postEta, objf, omModel,
                                  1.0, std::vector<int>(),
                                  std::vector<int>(), std::vector<int>());
  impIterPrintGet(e);          // closing rule + stash e$parHistData

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
  // A log/boxCox (lnorm) endpoint saw a non-positive prediction (floored by
  // rxode2); the fit still ran on the floored value -- surface it into $runInfo.
  if (impNpTbsDomainWarn()) {
    Rcpp::warning("lnorm/log endpoint has a <=0 prediction; drop 0-prediction obs");
  }
}
