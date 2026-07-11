// RPEM estimation engine (design/rpem/02, 04, 05).
//
// The whole E-step/M-step loop lives here in C++ (D17): drive rxode2's C-level
// solve in-process and read rx_pred_ straight from the solve buffer, with no
// per-iteration round-trip to R.
//
// This first increment provides the core E-step solve building block:
//   * rpemSetup()   -- bind the compiled rpem likelihood model + set up the
//                      rxode2 solve struct (mirrors nlm.cpp's nlmSetup).
//   * rpemSolvePop()-- solve every subject at a supplied per-subject parameter
//                      row (population THETA + that subject's ETA draw) and
//                      return per-subject -sum(rx_pred_) = -log p(Y_i|theta_i).
// The batched (subject x sample) E-step and the M-step build on top of this.

// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "inner.h"
#include "rxomp.h"
#include "censEst.h"   // doCensNormal1: reusable censored normal log-likelihood
#include <vector>
#include <functional>

#define _(String) (String)

// Solve one subject with the bound rpem predOnly model (rxPred function ptrs).
#define rpemPredOde(id) ind_solve(rx, id, rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)

struct rpemOptions {
  unsigned int ntheta = 0;   // number of constant params per subject (THETA + ETA); DV comes from data
  int *nobs = NULL;          // [nsub] observations per subject
  int *idS  = NULL;          // [nsub] start offset into flattened per-obs storage
  int *idF  = NULL;          // [nsub] final offset
  unsigned int nobsTot = 0;
  int nsub = 0;
  bool loaded = false;
  std::vector<int> buf;      // backing store for nobs/idS/idF
  // E-step sample store (populated by rpemEstepK1Draw) for M-step reuse
  int nGauss = 0;
  int nEta = 0;
  // Mixture components (mix()): 1 for a non-mixture fit.  The per-sample stores
  // below get a trailing component stride of nMix, so [i*nGauss+j] indexes the
  // (nMix==1) case unchanged and [(i*nGauss+j)*nMix + k] the mixture case.
  int nMix = 1;
  std::vector<double> logp;  // [nsub*nGauss] log p(Y_i | theta_ij)
  std::vector<double> etaS;  // [nsub*nGauss*nEta] drawn etas
  std::vector<double> ssv;   // [nsub*nGauss] additive SS = sum (DV-cp)^2
  std::vector<double> wssv;  // [nsub*nGauss] proportional WSS = sum ((DV-cp)/cp)^2
  // Per-observation raw structural prediction cp and observed DV, for the residual
  // M-steps that have no closed form and re-score per candidate residual parameter:
  // combined (add+prop), power (b*cp^c), and TBS lambda.  Block for sample (i,j)
  // starts at sampObsOff[i] + j*nobs_i, length nobs_i.
  std::vector<double> cpv;   // [nGauss*nobsTot] structural prediction cp
  std::vector<double> dvv;   // [nGauss*nobsTot] observed DV
  std::vector<long> sampObsOff; // [nsub] start offset of subject i's sample blocks
  // Per-observation transform code / bounds (constant across samples, indexed by
  // idS[i]+o), for the per-endpoint TBS lambda profile with boxCox/yeoJohnson.
  std::vector<int> yjv;      // [nobsTot] yj0 transform code
  std::vector<double> lowv;  // [nobsTot] transform lower bound
  std::vector<double> hiv;   // [nobsTot] transform upper bound
  // Per-observation censoring (constant across samples, indexed by idS[i]+o), for
  // BLQ / M2 / M3 / M4 handling via doCensNormal1.
  std::vector<int> censv;    // [nobsTot] CENS column (-1/0/1)
  std::vector<double> limv;  // [nobsTot] LIMIT column (or -Inf)
  bool anyCens = false;      // whether any observation is censored
};

rpemOptions rpemOp;

// golden-section maximizer (defined below); used by the power-residual profile in the
// mixture M-step as well as the non-mixture power / TBS M-steps.
static double rpemGolden(const std::function<double(double)> &f, double a, double b, int iter);

//[[Rcpp::export]]
RObject rpemFree() {
  rpemOp.buf.clear();
  rpemOp.nobs = rpemOp.idS = rpemOp.idF = NULL;
  rpemOp.nobsTot = 0;
  rpemOp.nsub = 0;
  rpemOp.loaded = false;
  rpemOp.cpv.clear();
  rpemOp.dvv.clear();
  rpemOp.sampObsOff.clear();
  rpemOp.yjv.clear();
  rpemOp.lowv.clear();
  rpemOp.hiv.clear();
  rpemOp.censv.clear();
  rpemOp.limv.clear();
  rpemOp.anyCens = false;
  return R_NilValue;
}

// Bind the rpem model and (re)establish the population solve struct `rx`, then
// recompute per-subject observation offsets.  Called at setup and again after any
// rxRmvn draw (which resets rxode2's global solve state and would otherwise leave
// our cached `rx` dangling to a 1-subject struct).
static void rpemDoSetup(RObject pred, List rxControl, NumericVector p, RObject data) {
  List mvp = rxode2::rxModelVars_(pred);
  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);
  rpemOp.ntheta = (unsigned int)p.size();
  rxode2::rxSolve_(pred, rxControl,
                   R_NilValue, R_NilValue,
                   p, data, R_NilValue,
                   1);         // setupOnly
  rx = getRxSolve_();
  rpemOp.anyCens = (hasRxCens(rx) != 0);
  rpemOp.nsub = getRxNsub(rx);
  rpemOp.buf.assign((size_t)rpemOp.nsub * 3u, 0);
  rpemOp.nobs = rpemOp.buf.data();
  rpemOp.idS  = rpemOp.nobs + rpemOp.nsub;
  rpemOp.idF  = rpemOp.idS + rpemOp.nsub;
  rpemOp.nobsTot = 0U;
  for (int id = 0; id < rpemOp.nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
    int no = 0;
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      if (getIndEvid(ind, j) == 0) { rpemOp.nobsTot++; no++; }
    }
    rpemOp.nobs[id] = no;
    if (id == 0) { rpemOp.idS[0] = 0; rpemOp.idF[0] = no - 1; }
    else { rpemOp.idS[id] = rpemOp.idF[id-1] + 1; rpemOp.idF[id] = rpemOp.idS[id] + no - 1; }
  }
  rpemOp.loaded = true;
}

//[[Rcpp::export]]
RObject rpemSetup(Environment e) {
  rpemFree();
  rpemDoSetup(e["predOnly"], as<List>(e["rxControl"]),
              as<NumericVector>(e["param"]), e["data"]);
  return R_NilValue;
}

// Solve subject `id` at parameter vector `par` (length ntheta); return the
// summed rx_pred_ over its observations (= -log p(Y_id | theta)).
// Returns the summed rx_pred_ (= -log p) for subject id.  If ssOut/wssOut are
// non-null, also accumulates the residual sum of squares SS = sum (DV - cp)^2
// (additive) and the proportional WSS = sum ((DV - cp)/cp)^2, from the cached
// structural prediction cp = rx_pred_f_ (lhs[1]) -- used by the residual M-step.
// If cpOut/dvOut are non-null they receive the per-observation raw structural
// prediction cp and observed DV (in observation order), for the numeric residual
// M-steps (combined, power, TBS lambda).
// Set subject `id`'s parameter pointers and (re)initialize its state, ready for a
// solve (serial ind_solve or a batched par_solve).
static inline void rpemSetSubject(const double *par, int id) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) setIndParPtr(ind, (int)i, par[i]);
  iniSubjectE(id, 1, ind, op, rx, rxPred.update_inis);
}

// Read subject `id`'s predictions AFTER its ODE has been solved (by ind_solve or
// par_solve): calc_lhs at each observation and accumulate the summed rx_pred_
// (= -log p), the additive/proportional sums of squares, and the raw per-obs
// cp / DV.  No solving happens here, so it is safe to call serially after a
// batched par_solve.
static inline double rpemReadSubject(int id, double *ssOut = nullptr, double *wssOut = nullptr,
                                     double *cpOut = nullptr, double *dvOut = nullptr,
                                     int *yjOut = nullptr, double *lowOut = nullptr, double *hiOut = nullptr,
                                     int *censOut = nullptr, double *limOut = nullptr) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  bool anyCens = rpemOp.anyCens;
  double s = 0.0, ss = 0.0, wss = 0.0, curT;
  int kk, oi = 0;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxPred.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxPred.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      double cp = lhs[1];                // rx_pred_f_ = structural prediction
      double dv = getIndDv(ind, kk);
      int cens = 0; double lim = R_NegInf;
      if (anyCens) {
        cens = getIndCens(ind, kk);
        if (hasRxLimit(rx)) { lim = getIndLimit(ind, kk); if (ISNA(lim)) lim = R_NegInf; }
      }
      double v;
      if (cens != 0 || (anyCens && R_FINITE(lim))) {
        // censored (M2/M3/M4): -loglik via the reusable censored-normal function,
        // using the per-obs variance rx_r_ = lhs[2].
        double ll = -lhs[0];             // uncensored (density) loglik
        v = -doCensNormal1((double)cens, dv, lim, ll, cp, lhs[2], 0);
      } else {
        v = lhs[0];                      // -density loglik (unchanged)
      }
      if (ISNA(v)) v = 0.0;
      s += v;
      double r = dv - cp;
      if (ssOut != nullptr) ss += r * r;
      if (wssOut != nullptr && cp != 0.0) { double rc = r / cp; wss += rc * rc; }
      if (cpOut != nullptr) cpOut[oi] = cp;
      if (dvOut != nullptr) dvOut[oi] = dv;
      // per-obs transform code / bounds (set by calc_lhs for this observation's
      // endpoint), for the per-endpoint TBS lambda profile.
      if (yjOut != nullptr) yjOut[oi] = getIndYj(ind);
      if (lowOut != nullptr) lowOut[oi] = getIndLogitLow(ind);
      if (hiOut != nullptr) hiOut[oi] = getIndLogitHi(ind);
      if (censOut != nullptr) censOut[oi] = cens;
      if (limOut != nullptr) limOut[oi] = lim;
      ++oi;
    }
  }
  if (ssOut != nullptr) *ssOut = ss;
  if (wssOut != nullptr) *wssOut = wss;
  return s;
}

// Serial solve of one subject (set params, solve its ODE, read predictions).
static inline double rpemSolveSubject(const double *par, int id,
                                      double *ssOut = nullptr, double *wssOut = nullptr,
                                      double *cpOut = nullptr, double *dvOut = nullptr,
                                      int *yjOut = nullptr, double *lowOut = nullptr, double *hiOut = nullptr,
                                      int *censOut = nullptr, double *limOut = nullptr) {
  rpemSetSubject(par, id);
  rpemPredOde(id);
  return rpemReadSubject(id, ssOut, wssOut, cpOut, dvOut, yjOut, lowOut, hiOut, censOut, limOut);
}

// parMat: nsub rows x ntheta cols (each subject's THETA+ETA). Returns per-subject
// summed rx_pred_ (= -log p(Y_i | theta_i)).
//[[Rcpp::export]]
NumericVector rpemSolvePop(NumericMatrix parMat) {
  if (!rpemOp.loaded) stop("'rpem' problem not loaded");
  int nsub = rpemOp.nsub;
  if (parMat.nrow() != nsub) stop("parMat must have one row per subject");
  if ((unsigned int)parMat.ncol() != rpemOp.ntheta) stop("parMat must have ntheta columns");
  NumericVector out(nsub);
  std::vector<double> row(rpemOp.ntheta);
  for (int id = 0; id < nsub; ++id) {
    for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = parMat(id, i);
    out[id] = rpemSolveSubject(row.data(), id);
  }
  return out;
}

// K=1 E-step Monte Carlo accumulation (design/rpem/04). parBig has
// nsub*nGauss rows x ntheta cols; rows [id*nGauss .. id*nGauss+nGauss-1] are the
// nGauss parameter draws (population THETA + sampled ETA) for subject id.
// Returns per-subject log n_i and lnL = sum_i log n_i, where
//   n_i = mean_j p(Y_i | theta_ij),  log p = -sum(rx_pred_)
// combined with a numerically stable log-sum-exp.
//[[Rcpp::export]]
List rpemEstepK1(NumericMatrix parBig, int nGauss) {
  if (!rpemOp.loaded) stop("'rpem' problem not loaded");
  int nsub = rpemOp.nsub;
  if (parBig.nrow() != nsub * nGauss) stop("parBig must have nsub*nGauss rows");
  if ((unsigned int)parBig.ncol() != rpemOp.ntheta) stop("parBig must have ntheta columns");
  NumericVector logn(nsub);
  std::vector<double> row(rpemOp.ntheta);
  std::vector<double> lp(nGauss);
  double lnL = 0.0;
  for (int id = 0; id < nsub; ++id) {
    double mx = R_NegInf;
    for (int j = 0; j < nGauss; ++j) {
      for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = parBig(id * nGauss + j, i);
      double logp = -rpemSolveSubject(row.data(), id);
      lp[j] = logp;
      if (logp > mx) mx = logp;
    }
    double s = 0.0;
    for (int j = 0; j < nGauss; ++j) s += exp(lp[j] - mx);
    double logni = mx + log(s) - log((double)nGauss);
    logn[id] = logni;
    lnL += logni;
  }
  return List::create(_["logn"] = logn, _["lnL"] = lnL);
}

// K=1 E-step.  The nGauss etas per subject ~ N(0, omega) are drawn IN R (before
// this call) and passed as etaMat, so the sampling RNG is fully decoupled from the
// solve -- the fit is reproducible regardless of the solve core count (D18, D21).
//   base:   length-ntheta template param row (ETA slots overwritten per draw)
//   etaIdx: 0-based positions of the ETA params in the row
//   etaMat: (nsub*nGauss) x nEta pre-drawn etas, layout [i*nGauss + j]
//   ncores: >1 solves each Monte Carlo sample's subjects in parallel via par_solve
// Solves, accumulates n_i via log-sum-exp, stores per-sample log p / etas / cp /
// DV for the M-step, and returns logn, lnL, and the etas.
// Optional mode-centered importance sampling (rpemControl(impInflate=)): etaMat is
// drawn per subject from N(ebe_i, cInflate*Omega) instead of the prior N(0, Omega),
// covering the posterior tails of high-variance random effects the prior proposal
// misses (which under-estimated the largest Omega in multi-eta models).  The stored
// rpemOp.logp becomes the IMPORTANCE LOG-WEIGHT logw = logLik + log N(eta;0,Omega) -
// log N(eta;ebe,cInflate*Omega), so every downstream weight / logn / MH acceptance is
// the correct self-normalized posterior; the raw residual SS stays in rpemOp.ssv.
// ebe=0, cInflate=1 (the default) reduces logw to logLik -- exactly the paper's prior
// sampling.  Returns the posterior mean (ebe) for the next iteration's proposal center.
//[[Rcpp::export]]
List rpemEstepK1Draw(Environment e, NumericVector base, IntegerVector etaIdx,
                     NumericMatrix etaMat, int nGauss, int ncores,
                     NumericMatrix ebeCenter, NumericVector omVec, double cInflate) {
  int nEta = etaIdx.size();

  RObject pred = e["predOnly"];
  List rxControl = as<List>(e["rxControl"]);
  NumericVector param = as<NumericVector>(e["param"]);
  RObject data = e["data"];

  // Establish the solve once (no in-C++ draw to clobber it).
  rpemDoSetup(pred, rxControl, param, data);
  int nsub = rpemOp.nsub;
  if ((unsigned int)base.size() != rpemOp.ntheta) stop("base must have ntheta entries");
  int nAll = nsub * nGauss;
  if (etaMat.nrow() != nAll || etaMat.ncol() != nEta) stop("etaMat must be (nsub*nGauss) x nEta");
  NumericMatrix etaAll = etaMat;

  rpemOp.nGauss = nGauss;
  rpemOp.nEta = nEta;
  rpemOp.nMix = 1;
  rpemOp.logp.assign((size_t)nAll, 0.0);
  rpemOp.etaS.assign((size_t)nAll * nEta, 0.0);
  rpemOp.ssv.assign((size_t)nAll, 0.0);
  rpemOp.wssv.assign((size_t)nAll, 0.0);
  // Per-obs raw cp / DV blocks for the numeric residual M-steps.
  rpemOp.sampObsOff.assign((size_t)nsub, 0);
  long acc = 0;
  for (int id = 0; id < nsub; ++id) { rpemOp.sampObsOff[id] = acc; acc += (long)nGauss * rpemOp.nobs[id]; }
  rpemOp.cpv.assign((size_t)acc, 0.0);
  rpemOp.dvv.assign((size_t)acc, 0.0);
  rpemOp.yjv.assign((size_t)rpemOp.nobsTot, 0);      // per subject-obs (idS[i]+o)
  rpemOp.lowv.assign((size_t)rpemOp.nobsTot, 0.0);
  rpemOp.hiv.assign((size_t)rpemOp.nobsTot, 1.0);
  rpemOp.censv.assign((size_t)rpemOp.nobsTot, 0);
  rpemOp.limv.assign((size_t)rpemOp.nobsTot, R_NegInf);

  // Copy all R objects into plain C++ buffers BEFORE the parallel region; the
  // OpenMP loop must not touch any Rcpp/R object (only these buffers + the
  // rxode2 C solve). etaS doubles as the shared eta buffer for the loop.
  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) rpemOp.etaS[k] = etaAll[k];
  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) baseBuf[i] = base[i];
  std::vector<int> etaIdxBuf(nEta);
  for (int a = 0; a < nEta; ++a) etaIdxBuf[a] = etaIdx[a];
  std::vector<double> lognV(nsub, 0.0);

  // Importance-weight buffers: proposal center ebe_i (nsub x nEta) and prior variances.
  if (ebeCenter.nrow() != nsub || ebeCenter.ncol() != nEta) stop("ebeCenter must be nsub x nEta");
  if ((int)omVec.size() != nEta) stop("omVec must have nEta entries");
  std::vector<double> ebeBuf((size_t)nsub * nEta), omBuf(nEta);
  for (size_t k = 0; k < (size_t)nsub * nEta; ++k) ebeBuf[k] = ebeCenter[k];
  for (int a = 0; a < nEta; ++a) omBuf[a] = omVec[a];
  double halfNetaLogC = 0.5 * nEta * log(cInflate);
  // logw - logLik: prior N(eta;0,Om) over proposal N(eta;ebe,cInflate*Om), diagonal Om.
  auto logRatio = [&](int id, int j) -> double {
    double lr = halfNetaLogC;
    for (int a = 0; a < nEta; ++a) {
      double eta = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
      double om = omBuf[a];
      double d = eta - ebeBuf[(size_t)id * nEta + a];
      lr += -0.5 * eta * eta / om + 0.5 * d * d / (cInflate * om);
    }
    return lr;
  };

  std::vector<double> row(rpemOp.ntheta), lp(nGauss);
  if (ncores > 1) {
    // Parallel E-step: for each Monte Carlo sample, set every subject's parameters
    // (population THETA + that subject's eta draw) and solve all subjects in one
    // par_solve call.  par_solve owns the per-thread ODE workspace, so this is the
    // thread-safe route (a manual OpenMP loop over ind_solve segfaults).  The cheap
    // prediction read (calc_lhs, no ODE work) stays serial.  Each subject's solve
    // is correct to solver tolerance for any thread count; the stochastic MH M-step
    // means cross-core fits can differ negligibly (Monte-Carlo level), but a fixed
    // core count is fully reproducible.
    for (int j = 0; j < nGauss; ++j) {
      for (int id = 0; id < nsub; ++id) {
        for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
        for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
        rpemSetSubject(row.data(), id);
        setIndSolve(getSolvingOptionsInd(rx, id), -1);
      }
      resetRxBadSolve(rx);
      par_solve(rx);
      for (int id = 0; id < nsub; ++id) {
        size_t r = (size_t)id * nGauss + j;
        double ssTmp = 0.0, wssTmp = 0.0;
        long ob = rpemOp.sampObsOff[id] + (long)j * rpemOp.nobs[id];
        long so = rpemOp.idS[id];
        double logLik = -rpemReadSubject(id, &ssTmp, &wssTmp, &rpemOp.cpv[ob], &rpemOp.dvv[ob],
                                         &rpemOp.yjv[so], &rpemOp.lowv[so], &rpemOp.hiv[so],
                                         &rpemOp.censv[so], &rpemOp.limv[so]);
        rpemOp.logp[r] = logLik + logRatio(id, j); rpemOp.ssv[r] = ssTmp; rpemOp.wssv[r] = wssTmp;
      }
    }
    for (int id = 0; id < nsub; ++id) {
      double mx = R_NegInf;
      for (int j = 0; j < nGauss; ++j) { double v = rpemOp.logp[(size_t)id * nGauss + j]; if (v > mx) mx = v; }
      double s = 0.0;
      for (int j = 0; j < nGauss; ++j) s += exp(rpemOp.logp[(size_t)id * nGauss + j] - mx);
      lognV[id] = mx + log(s) - log((double)nGauss);
    }
  } else {
    // Serial reference path (cores == 1): solve each subject with ind_solve.
    for (int id = 0; id < nsub; ++id) {
      double mx = R_NegInf;
      for (int j = 0; j < nGauss; ++j) {
        size_t r = (size_t)id * nGauss + j;
        for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
        for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[r * nEta + a];
        double ssTmp = 0.0, wssTmp = 0.0;
        long ob = rpemOp.sampObsOff[id] + (long)j * rpemOp.nobs[id];
        long so = rpemOp.idS[id];
        double logLik = -rpemSolveSubject(row.data(), id, &ssTmp, &wssTmp,
                                          &rpemOp.cpv[ob], &rpemOp.dvv[ob],
                                          &rpemOp.yjv[so], &rpemOp.lowv[so], &rpemOp.hiv[so],
                                          &rpemOp.censv[so], &rpemOp.limv[so]);
        double lw = logLik + logRatio(id, j);
        lp[j] = lw;
        rpemOp.logp[r] = lw;
        rpemOp.ssv[r] = ssTmp;
        rpemOp.wssv[r] = wssTmp;
        if (lw > mx) mx = lw;
      }
      double s = 0.0;
      for (int j = 0; j < nGauss; ++j) s += exp(lp[j] - mx);
      lognV[id] = mx + log(s) - log((double)nGauss);
    }
  }

  // Build the R return objects serially, after the parallel region.
  NumericVector logn(nsub);
  NumericMatrix etaOut(nAll, nEta);
  NumericVector logpOut(nAll);        // layout [i*nGauss + j]
  NumericMatrix ebeOut(nsub, nEta);   // posterior-mean eta (proposal center next iter)
  double lnL = 0.0;
  for (int id = 0; id < nsub; ++id) { logn[id] = lognV[id]; lnL += lognV[id]; }
  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) etaOut[k] = rpemOp.etaS[k];
  for (int k = 0; k < nAll; ++k) logpOut[k] = rpemOp.logp[k];
  // EBE_i = sum_j w_ij eta_ij with w_ij = softmax_j(logw) -- the self-normalized
  // importance weights (posterior mean), fed back as the next proposal center.
  for (int id = 0; id < nsub; ++id) {
    double mx = R_NegInf;
    for (int j = 0; j < nGauss; ++j) { double v = rpemOp.logp[(size_t)id * nGauss + j]; if (v > mx) mx = v; }
    double sw = 0.0;
    for (int j = 0; j < nGauss; ++j) sw += exp(rpemOp.logp[(size_t)id * nGauss + j] - mx);
    for (int a = 0; a < nEta; ++a) {
      double acc = 0.0;
      for (int j = 0; j < nGauss; ++j)
        acc += exp(rpemOp.logp[(size_t)id * nGauss + j] - mx) / sw *
               rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
      ebeOut(id, a) = acc;
    }
  }
  return List::create(_["logn"] = logn, _["lnL"] = lnL, _["eta"] = etaOut,
                      _["logp"] = logpOut, _["ebe"] = ebeOut);
}

// Mixture E-step (mix(), split-ETA).  The mixed parameter's typical value is
// selected by the rxode2 model via setIndMixest, so we solve every Monte Carlo
// sample once per component k (1..K), storing the per-component log p(Y_i|k,eta_ij)
// with a trailing stride of nMix=K.  The subject likelihood mixes the components,
//   n_i = sum_k w_k * mean_j p(Y_i | k, eta_ij),
// and the returned per-component log n_ik = logsumexp_j logp_ijk - log nGauss lets
// the M-step form the component posteriors tau_ik.  Shares one eta draw across
// components (etaMat), so the mixture differs only in the typical value.
//   w: length-K mixture weights (probabilities), used only for lnL/logn here.
//[[Rcpp::export]]
List rpemEstepMixDraw(Environment e, NumericVector base, IntegerVector etaIdx,
                      NumericMatrix etaMat, int nGauss, int ncores,
                      int K, NumericVector w) {
  int nEta = etaIdx.size();
  RObject pred = e["predOnly"];
  List rxControl = as<List>(e["rxControl"]);
  NumericVector param = as<NumericVector>(e["param"]);
  RObject data = e["data"];

  rpemDoSetup(pred, rxControl, param, data);
  int nsub = rpemOp.nsub;
  if ((unsigned int)base.size() != rpemOp.ntheta) stop("base must have ntheta entries");
  if (K < 1) stop("K must be >= 1");
  if (w.size() != K) stop("w must have K entries");
  int nAll = nsub * nGauss;
  if (etaMat.nrow() != nAll || etaMat.ncol() != nEta) stop("etaMat must be (nsub*nGauss) x nEta");

  rpemOp.nGauss = nGauss;
  rpemOp.nEta = nEta;
  rpemOp.nMix = K;
  rpemOp.logp.assign((size_t)nAll * K, 0.0);
  rpemOp.etaS.assign((size_t)nAll * nEta, 0.0);
  rpemOp.ssv.assign((size_t)nAll * K, 0.0);
  rpemOp.wssv.assign((size_t)nAll * K, 0.0);
  rpemOp.sampObsOff.assign((size_t)nsub, 0);
  long acc = 0;
  for (int id = 0; id < nsub; ++id) { rpemOp.sampObsOff[id] = acc; acc += (long)nGauss * rpemOp.nobs[id]; }
  rpemOp.cpv.assign((size_t)acc * K, 0.0);
  rpemOp.dvv.assign((size_t)acc * K, 0.0);
  rpemOp.yjv.assign((size_t)rpemOp.nobsTot, 0);
  rpemOp.lowv.assign((size_t)rpemOp.nobsTot, 0.0);
  rpemOp.hiv.assign((size_t)rpemOp.nobsTot, 1.0);
  rpemOp.censv.assign((size_t)rpemOp.nobsTot, 0);
  rpemOp.limv.assign((size_t)rpemOp.nobsTot, R_NegInf);

  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) rpemOp.etaS[k] = etaMat[k];
  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) baseBuf[i] = base[i];
  std::vector<int> etaIdxBuf(nEta);
  for (int a = 0; a < nEta; ++a) etaIdxBuf[a] = etaIdx[a];
  std::vector<double> row(rpemOp.ntheta);

  // Solve each Monte Carlo sample once per component; the yj/cens per-obs buffers
  // are constant across components/samples so only the last write is kept.
  for (int kc = 0; kc < K; ++kc) {
    for (int j = 0; j < nGauss; ++j) {
      if (ncores > 1) {
        for (int id = 0; id < nsub; ++id) {
          for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
          for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
          rpemSetSubject(row.data(), id);
          setIndMixest(getSolvingOptionsInd(rx, id), kc + 1);   // 1-based component
          setIndSolve(getSolvingOptionsInd(rx, id), -1);
        }
        resetRxBadSolve(rx);
        par_solve(rx);
        for (int id = 0; id < nsub; ++id) {
          size_t r = ((size_t)id * nGauss + j) * K + kc;
          double ssTmp = 0.0, wssTmp = 0.0;
          long ob = (rpemOp.sampObsOff[id] + (long)j * rpemOp.nobs[id]) * K + (long)kc * rpemOp.nobs[id];
          long so = rpemOp.idS[id];
          double logp = -rpemReadSubject(id, &ssTmp, &wssTmp, &rpemOp.cpv[ob], &rpemOp.dvv[ob],
                                         &rpemOp.yjv[so], &rpemOp.lowv[so], &rpemOp.hiv[so],
                                         &rpemOp.censv[so], &rpemOp.limv[so]);
          rpemOp.logp[r] = logp; rpemOp.ssv[r] = ssTmp; rpemOp.wssv[r] = wssTmp;
        }
      } else {
        for (int id = 0; id < nsub; ++id) {
          size_t r = ((size_t)id * nGauss + j) * K + kc;
          for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
          for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
          rpemSetSubject(row.data(), id);
          setIndMixest(getSolvingOptionsInd(rx, id), kc + 1);
          rpemPredOde(id);
          double ssTmp = 0.0, wssTmp = 0.0;
          long ob = (rpemOp.sampObsOff[id] + (long)j * rpemOp.nobs[id]) * K + (long)kc * rpemOp.nobs[id];
          long so = rpemOp.idS[id];
          double logp = -rpemReadSubject(id, &ssTmp, &wssTmp, &rpemOp.cpv[ob], &rpemOp.dvv[ob],
                                         &rpemOp.yjv[so], &rpemOp.lowv[so], &rpemOp.hiv[so],
                                         &rpemOp.censv[so], &rpemOp.limv[so]);
          rpemOp.logp[r] = logp; rpemOp.ssv[r] = ssTmp; rpemOp.wssv[r] = wssTmp;
        }
      }
    }
  }

  // Per-component log n_ik and mixture log n_i = logsumexp_k (log w_k + log n_ik).
  NumericVector logn(nsub);
  NumericMatrix lognik(nsub, K);       // log mean_j p(Y_i | k, eta_ij)
  NumericMatrix etaOut(nAll, nEta);
  double lnL = 0.0;
  std::vector<double> logw(K);
  for (int kc = 0; kc < K; ++kc) logw[kc] = log(w[kc]);
  for (int id = 0; id < nsub; ++id) {
    double mxk = R_NegInf;
    for (int kc = 0; kc < K; ++kc) {
      double mx = R_NegInf;
      for (int j = 0; j < nGauss; ++j) {
        double v = rpemOp.logp[((size_t)id * nGauss + j) * K + kc];
        if (v > mx) mx = v;
      }
      double s = 0.0;
      for (int j = 0; j < nGauss; ++j) s += exp(rpemOp.logp[((size_t)id * nGauss + j) * K + kc] - mx);
      double lnik = mx + log(s) - log((double)nGauss);
      lognik(id, kc) = lnik;
      double wk = logw[kc] + lnik;
      if (wk > mxk) mxk = wk;
    }
    double sk = 0.0;
    for (int kc = 0; kc < K; ++kc) sk += exp(logw[kc] + lognik(id, kc) - mxk);
    double lni = mxk + log(sk);
    logn[id] = lni; lnL += lni;
  }
  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) etaOut[k] = rpemOp.etaS[k];
  // Per-component per-sample log p (nAll x K, row = i*nGauss+j) for EBE weights.
  NumericMatrix logpOut(nAll, K);
  for (int r = 0; r < nAll; ++r)
    for (int kc = 0; kc < K; ++kc) logpOut(r, kc) = rpemOp.logp[(size_t)r * K + kc];
  return List::create(_["logn"] = logn, _["lnL"] = lnL, _["eta"] = etaOut,
                      _["lognik"] = lognik, _["logp"] = logpOut);
}

// K=1 M-step: Metropolis-Hastings over the stored E-step samples (design/rpem/05).
// State s=(i,j) indexes a stored sample; propose (i',j') uniformly and accept by
// Eq 32 with w=1: A = min(1, exp((logp' - logp) + (logn_i - logn_i'))).  The
// accepted theta = muIn + eta samples give the conjugate updates (Eq 15-16):
//   mu^(1)    = mean of accepted theta
//   Omega^(1) = covariance of accepted theta about the new mu
// No solving: reuses stored log p / etas.  All randomness is pre-drawn up front
// via the threefry rxRmvn (uniforms via pnorm of standard normals, D19); the
// draw resets rxode2's solve state but the M-step does not solve, so that is fine.
// addSd0 is the additive residual SD used in the E-step; the additive residual is
// updated in closed form (paper Eq 17) from the accepted samples' residual sum of
// squares, backed out of the stored log p (which is exactly the additive-normal
// log-likelihood, C1.1): SS_ij = -2 addSd0^2 (logp_ij + nobs_i (0.5 log 2pi +
// log addSd0)); new addSd = sqrt(sum SS / sum nobs).  (General error structures
// use numeric re-scoring instead -- deferred; see design/rpem/05.)
//[[Rcpp::export]]
List rpemMstepK1(NumericVector muIn, double addSd0, int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, nEta = rpemOp.nEta;
  if ((int)muIn.size() != nEta) stop("muIn must have nEta entries");
  (void)addSd0;

  // Per-subject log n_i (= log of the MC-mean importance weight) from the stored logw.
  std::vector<double> logn(nsub);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) s += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    logn[i] = mx + log(s) - log((double)nG);
  }

  // Pre-draw all MH uniforms: 3 per trial (proposal i', proposal j', accept).
  int total = nTrials + burn;
  size_t nU = (size_t)3 * total;
  NumericVector mu0(1); NumericMatrix s0(1, 1); s0(0, 0) = 1.0;
  NumericVector lo(1, R_NegInf), hi(1, R_PosInf);
  SEXP zr = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create((int)nU)),
                                wrap(mu0), wrap(s0), wrap(lo), wrap(hi),
                                wrap(IntegerVector::create(1)),
                                wrap(LogicalVector::create(false)),
                                wrap(LogicalVector::create(false)),
                                wrap(NumericVector::create(0.4)),
                                wrap(NumericVector::create(2.05)),
                                wrap(NumericVector::create(1e-10)),
                                wrap(IntegerVector::create(100)));
  NumericMatrix zm(zr);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(zm[(int)k], 0.0, 1.0, 1, 0);

  int ci = 0, cj = 0;
  double clogp = rpemOp.logp[0];
  std::vector<double> sumT(nEta, 0.0), sumTT((size_t)nEta * nEta, 0.0);
  double sumSS = 0.0; long sumNobs = 0;
  long m = 0, naccept = 0;
  for (int t = 0; t < total; ++t) {
    double u1 = U[(size_t)3 * t], u2 = U[(size_t)3 * t + 1], u3 = U[(size_t)3 * t + 2];
    int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
    int pjh = (int)(u2 * nG);   if (pjh >= nG)   pjh = nG - 1;
    double plogp = rpemOp.logp[(size_t)pih * nG + pjh];
    double logA = (plogp - clogp) + (logn[ci] - logn[pih]);
    if (log(u3) < logA) { ci = pih; cj = pjh; clogp = plogp; ++naccept; }
    if (t >= burn) {
      size_t off = ((size_t)ci * nG + cj) * nEta;
      for (int a = 0; a < nEta; ++a) {
        double tha = muIn[a] + rpemOp.etaS[off + a];
        sumT[a] += tha;
        for (int b = 0; b < nEta; ++b)
          sumTT[(size_t)a * nEta + b] += tha * (muIn[b] + rpemOp.etaS[off + b]);
      }
      // Additive residual: this sample's raw SS = sum (dv - cp)^2 (stored in ssv;
      // logp now holds the importance weight, so it can no longer back out SS).
      int nobsi = rpemOp.nobs[ci];
      sumSS += rpemOp.ssv[(size_t)ci * nG + cj]; sumNobs += nobsi;
      ++m;
    }
  }
  NumericVector muNew(nEta);
  NumericMatrix omegaNew(nEta, nEta);
  for (int a = 0; a < nEta; ++a) muNew[a] = sumT[a] / m;
  for (int a = 0; a < nEta; ++a)
    for (int b = 0; b < nEta; ++b)
      omegaNew(a, b) = sumTT[(size_t)a * nEta + b] / m - muNew[a] * muNew[b];
  double addNew = sqrt(sumSS / (double)sumNobs);
  return List::create(_["mu"] = muNew, _["omega"] = omegaNew,
                      _["addSd"] = addNew, _["accept"] = (double)naccept / total);
}

// Full C++ EM loop for the additive / proportional, diagonal multi-eta, mu-referenced
// core (design/rpem/12 M5): runs niter iterations of the E-step (threefry eta draw +
// par_solve + log-sum-exp) and the conjugate MH M-step in one C++ call, avoiding the
// per-iteration R round-trip.  The eta draw uses rxode2's per-thread threefry engine
// with a deterministic per-(iteration, subject) seed -- setSeedEng1(seed0 + it*nsub + i)
// after setRxThreadId() -- so it is thread-safe AND reproducible for any core count.
// The ODE solve keeps par_solve (rxode2 owns the parallel solver; a manual OpenMP loop
// over the solve is unsafe).  errType 0 additive, 1 proportional.  Returns per-iteration
// mu / omega (diagonal) / add.sd / lnL traces.
//[[Rcpp::export]]
List rpemEMLoopK1(Environment e, NumericVector base, IntegerVector etaIdx,
                  IntegerVector muIdx, int addSdIdx, int errType,
                  NumericVector mu0, NumericVector omDiag0, double addSd0,
                  IntegerVector resIdx, NumericVector resPar0,
                  int niter, int nGauss, int ncores, int nMH, int mhBurn,
                  unsigned int seed) {
  RObject pred = e["predOnly"];
  List rxControl = as<List>(e["rxControl"]);
  NumericVector param = as<NumericVector>(e["param"]);
  RObject data = e["data"];
  int nEta = etaIdx.size();
  if (rxNormEng == NULL || seedEng == NULL || setSeedEng1 == NULL)
    stop("rxode2 threefry engine not initialized");
  // second residual parameter (index / initial value): [prop.sd, power, lambda];
  // index -1 when the model has no such parameter.  Combined / power / TBS have no
  // closed-form residual, so the M-step accumulates per-(subject, sample) visit counts
  // and re-optimizes (same optimizers as the non-mixture combined/power/TBS M-steps).
  bool doComb = (errType == 2), doPow = (errType == 4), doTbs = (errType == 3);
  int propSdIdx = resIdx[0], powIdx = resIdx[1], lambdaIdx = resIdx[2];
  double propSd = resPar0[0], power = resPar0[1], lambda = resPar0[2];
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  if ((int)muIdx.size() != nEta || (int)mu0.size() != nEta || (int)omDiag0.size() != nEta)
    stop("muIdx / mu0 / omDiag0 must have nEta entries");

  rpemDoSetup(pred, rxControl, param, data);
  int nsub = rpemOp.nsub, nAll = nsub * nGauss;
  rpemOp.nGauss = nGauss; rpemOp.nEta = nEta; rpemOp.nMix = 1;
  rpemOp.logp.assign((size_t)nAll, 0.0);
  rpemOp.etaS.assign((size_t)nAll * nEta, 0.0);
  rpemOp.ssv.assign((size_t)nAll, 0.0);
  rpemOp.wssv.assign((size_t)nAll, 0.0);
  rpemOp.sampObsOff.assign((size_t)nsub, 0);
  long acc = 0;
  for (int id = 0; id < nsub; ++id) { rpemOp.sampObsOff[id] = acc; acc += (long)nGauss * rpemOp.nobs[id]; }
  rpemOp.cpv.assign((size_t)acc, 0.0);
  rpemOp.dvv.assign((size_t)acc, 0.0);
  rpemOp.yjv.assign((size_t)rpemOp.nobsTot, 0);
  rpemOp.lowv.assign((size_t)rpemOp.nobsTot, 0.0);
  rpemOp.hiv.assign((size_t)rpemOp.nobsTot, 1.0);
  rpemOp.censv.assign((size_t)rpemOp.nobsTot, 0);
  rpemOp.limv.assign((size_t)rpemOp.nobsTot, R_NegInf);

  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) baseBuf[i] = base[i];
  std::vector<int> etaIdxBuf(nEta), muIdxBuf(nEta);
  for (int a = 0; a < nEta; ++a) { etaIdxBuf[a] = etaIdx[a]; muIdxBuf[a] = muIdx[a]; }
  std::vector<double> mu(nEta), omDiag(nEta);
  for (int a = 0; a < nEta; ++a) { mu[a] = mu0[a]; omDiag[a] = omDiag0[a]; }
  double addSd = addSd0;

  NumericMatrix muTr(niter, nEta), omTr(niter, nEta);
  NumericVector sdTr(niter), propTr(niter, NA_REAL), powTr(niter, NA_REAL);
  NumericVector lamTr(niter, NA_REAL), llTr(niter);
  std::vector<double> row(rpemOp.ntheta);
  std::vector<long> counts3((doComb || doPow || doTbs) ? (size_t)nsub * nGauss : 0, 0);
  bool doPar = (ncores > 1);
  seedEng(ncores);
  uint32_t seed0 = (uint32_t)seed;

  for (int it = 0; it < niter; ++it) {
    for (int a = 0; a < nEta; ++a) baseBuf[muIdxBuf[a]] = mu[a];
    baseBuf[addSdIdx] = addSd;
    if (doComb) baseBuf[propSdIdx] = propSd;
    if (doPow) baseBuf[powIdx] = power;
    if (doTbs) baseBuf[lambdaIdx] = lambda;
    // re-establish the solve (the previous M-step's rxRmvn resets rxode2 solve state)
    rpemDoSetup(pred, rxControl, param, data);

    // threefry eta draw: eta_ija = z * sqrt(omDiag[a]).  Deterministic per-(iter,subject)
    // seed -> thread-safe and reproducible independent of core count.
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) if(doPar)
#endif
    for (int id = 0; id < nsub; ++id) {
#ifdef _OPENMP
      if (doPar) setRxThreadId(omp_get_thread_num());
#endif
      setSeedEng1(seed0 + (uint32_t)((size_t)it * nsub + id) + 1u);
      for (int j = 0; j < nGauss; ++j)
        for (int a = 0; a < nEta; ++a)
          rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a] = rxNormEng(0.0, 1.0) * sqrt(omDiag[a]);
    }

    // E-step: solve each MC sample's subjects in one par_solve, read the likelihood.
    for (int j = 0; j < nGauss; ++j) {
      for (int id = 0; id < nsub; ++id) {
        for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
        for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
        rpemSetSubject(row.data(), id);
        setIndSolve(getSolvingOptionsInd(rx, id), -1);
      }
      resetRxBadSolve(rx);
      par_solve(rx);
      for (int id = 0; id < nsub; ++id) {
        size_t r = (size_t)id * nGauss + j;
        double ssTmp = 0.0, wssTmp = 0.0;
        long ob = rpemOp.sampObsOff[id] + (long)j * rpemOp.nobs[id];
        long so = rpemOp.idS[id];
        double logLik = -rpemReadSubject(id, &ssTmp, &wssTmp, &rpemOp.cpv[ob], &rpemOp.dvv[ob],
                                         &rpemOp.yjv[so], &rpemOp.lowv[so], &rpemOp.hiv[so],
                                         &rpemOp.censv[so], &rpemOp.limv[so]);
        rpemOp.logp[r] = logLik; rpemOp.ssv[r] = ssTmp; rpemOp.wssv[r] = wssTmp;
      }
    }
    // per-subject log n_i and lnL.
    std::vector<double> logn(nsub); double lnL = 0.0;
    for (int id = 0; id < nsub; ++id) {
      double mx = R_NegInf;
      for (int j = 0; j < nGauss; ++j) { double v = rpemOp.logp[(size_t)id * nGauss + j]; if (v > mx) mx = v; }
      double s = 0.0;
      for (int j = 0; j < nGauss; ++j) s += exp(rpemOp.logp[(size_t)id * nGauss + j] - mx);
      logn[id] = mx + log(s) - log((double)nGauss); lnL += logn[id];
    }

    // M-step: conjugate MH over the joint (subject, sample) posterior (Eq 17).  The MH
    // uniforms come from the SAME threefry engine as the eta draw (serial chain on thread
    // 0) with a per-iteration seed disjoint from the eta-draw seeds, so the whole loop is
    // reproducible run-to-run (no dependence on rxode2's advancing global RNG state).
    int total = nMH + mhBurn;
    size_t nU = (size_t)3 * total;
#ifdef _OPENMP
    setRxThreadId(0);
#endif
    setSeedEng1(seed0 + (uint32_t)((size_t)niter * nsub) + (uint32_t)it + 1u);
    std::vector<double> U(nU);
    for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(rxNormEng(0.0, 1.0), 0.0, 1.0, 1, 0);
    int ci = 0, cj = 0; double clogp = rpemOp.logp[0];
    std::vector<double> sumT(nEta, 0.0), sumTT((size_t)nEta * nEta, 0.0);
    double sumSS = 0.0; long sumNobs = 0, m = 0;
    for (int t = 0; t < total; ++t) {
      double u1 = U[(size_t)3 * t], u2 = U[(size_t)3 * t + 1], u3 = U[(size_t)3 * t + 2];
      int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
      int pjh = (int)(u2 * nGauss); if (pjh >= nGauss) pjh = nGauss - 1;
      double plogp = rpemOp.logp[(size_t)pih * nGauss + pjh];
      double logA = (plogp - clogp) + (logn[ci] - logn[pih]);
      if (log(u3) < logA) { ci = pih; cj = pjh; clogp = plogp; }
      if (t >= mhBurn) {
        size_t off = ((size_t)ci * nGauss + cj) * nEta;
        for (int a = 0; a < nEta; ++a) {
          double tha = mu[a] + rpemOp.etaS[off + a];
          sumT[a] += tha;
          for (int b = 0; b < nEta; ++b) sumTT[(size_t)a * nEta + b] += tha * (mu[b] + rpemOp.etaS[off + b]);
        }
        sumSS += (errType == 1) ? rpemOp.wssv[(size_t)ci * nGauss + cj] : rpemOp.ssv[(size_t)ci * nGauss + cj];
        if (doComb || doPow || doTbs) counts3[(size_t)ci * nGauss + cj]++;
        sumNobs += rpemOp.nobs[ci]; ++m;
      }
    }
    for (int a = 0; a < nEta; ++a) mu[a] = sumT[a] / (double)m;
    for (int a = 0; a < nEta; ++a) omDiag[a] = sumTT[(size_t)a * nEta + a] / (double)m - mu[a] * mu[a];
    // Residual: additive/proportional pooled SS; combined / power / TBS have no closed
    // form, so re-optimize over the accepted (subject, sample) states' stored cp/dv
    // (same optimizers as the non-mixture combined / power / TBS M-steps).
    if (doComb) {
      const double eps = 1e-12;
      auto over = [&](double aa, double bb, double &Q, double *ga, double *gb,
                      double *Haa, double *Hab, double *Hbb) {
        Q = 0.0; if (ga) { *ga = *gb = *Haa = *Hab = *Hbb = 0.0; }
        for (int i = 0; i < nsub; ++i) {
          int nobsi = rpemOp.nobs[i];
          for (int j = 0; j < nGauss; ++j) {
            long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
            long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
            for (int o = 0; o < nobsi; ++o) {
              double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
              double wcp = cp * cp, r2 = rr * rr, V = aa + bb * wcp; if (V < eps) V = eps;
              double invV = 1.0 / V; Q += (double)c * (-0.5 * log(V) - 0.5 * r2 * invV);
              if (ga) {
                double g = -0.5 * invV + 0.5 * r2 * invV * invV;
                double h = 0.5 * invV * invV - r2 * invV * invV * invV;
                *ga += c * g; *gb += c * wcp * g;
                *Haa += c * h; *Hab += c * wcp * h; *Hbb += c * wcp * wcp * h;
              }
            }
          }
        }
      };
      double a = addSd * addSd, bpar = propSd * propSd, Qcur; over(a, bpar, Qcur, 0,0,0,0,0);
      for (int itn = 0; itn < 200; ++itn) {
        double Qtmp, ga, gb, Haa, Hab, Hbb; over(a, bpar, Qtmp, &ga, &gb, &Haa, &Hab, &Hbb);
        if (fabs(ga) + fabs(gb) < 1e-9) break;
        double det = Haa * Hbb - Hab * Hab, da, db;
        if (Haa < 0.0 && det > 0.0) { da = -(Hbb * ga - Hab * gb) / det; db = -(-Hab * ga + Haa * gb) / det; }
        else { double sc = 1.0 / (fabs(Haa) + fabs(Hbb) + 1.0); da = sc * ga; db = sc * gb; }
        double step = 1.0; bool ok = false;
        for (int ls = 0; ls < 40; ++ls) {
          double an = a + step * da, bn = bpar + step * db;
          if (an < eps) an = eps; if (bn < eps) bn = eps;
          double Qn; over(an, bn, Qn, 0,0,0,0,0);
          if (Qn > Qcur) { a = an; bpar = bn; Qcur = Qn; ok = true; break; }
          step *= 0.5;
        }
        if (!ok) break;
      }
      addSd = sqrt(a); propSd = sqrt(bpar);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doPow) {
      auto stat = [&](double cc, double &SSc, double &SL) {
        SSc = 0.0; SL = 0.0;
        for (int i = 0; i < nsub; ++i) {
          int nobsi = rpemOp.nobs[i];
          for (int j = 0; j < nGauss; ++j) {
            long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
            long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
            for (int o = 0; o < nobsi; ++o) {
              double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp, acp = fabs(cp) + 1e-300;
              SSc += (double)c * rr * rr / pow(acp, 2.0 * cc); SL += (double)c * log(acp);
            }
          }
        }
      };
      double N = (double)sumNobs;
      auto f = [&](double cc) { double SSc, SL; stat(cc, SSc, SL);
        if (SSc < 1e-300) SSc = 1e-300; return -0.5 * (N * log(SSc / N) + 2.0 * cc * SL); };
      power = rpemGolden(f, 0.0, 3.0, 100);
      double SSc, SL; stat(power, SSc, SL); addSd = sqrt(SSc / N);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doTbs) {
      auto ssJac = [&](double lam, double &SS, double &Jac) {
        SS = 0.0; Jac = 0.0;
        for (int i = 0; i < nsub; ++i) {
          int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
          for (int j = 0; j < nGauss; ++j) {
            long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
            long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
            for (int o = 0; o < nobsi; ++o) {
              double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
              int yj = rpemOp.yjv[so0 + o]; double low = rpemOp.lowv[so0 + o], hi = rpemOp.hiv[so0 + o];
              double d = _powerD(dv, lam, yj, low, hi) - _powerD(cp, lam, yj, low, hi);
              SS += (double)c * d * d;
              Jac += (double)c * log(fabs(_powerDD(dv, lam, yj, low, hi)) + 1e-300);
            }
          }
        }
      };
      double N = (double)sumNobs;
      auto f = [&](double lam) { double SS, Jac; ssJac(lam, SS, Jac);
        if (SS < 1e-300) SS = 1e-300; return -0.5 * N * log(SS / N) + Jac; };
      lambda = rpemGolden(f, -2.0, 3.0, 100);
      double SS, Jac; ssJac(lambda, SS, Jac); addSd = sqrt(SS / N);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else {
      addSd = sqrt(sumSS / (double)sumNobs);
    }
    for (int a = 0; a < nEta; ++a) { muTr(it, a) = mu[a]; omTr(it, a) = omDiag[a]; }
    sdTr[it] = addSd; llTr[it] = lnL;
    if (doComb) propTr[it] = propSd;
    if (doPow) powTr[it] = power;
    if (doTbs) lamTr[it] = lambda;
  }
  rpemFree();
  return List::create(_["muTrace"] = muTr, _["omegaTrace"] = omTr,
                      _["sdTrace"] = sdTr, _["propTrace"] = propTr,
                      _["powTrace"] = powTr, _["lamTrace"] = lamTr, _["lnL"] = llTr);
}

// Mixture M-step (mix(), split-ETA, single mixed eta).  The MH state is now
// (i, j, k) -- subject, Monte Carlo sample, and component -- with stationary
// distribution proportional to w_k p(Y_i | k, eta_ij) / n_i (the joint posterior
// of subject, sample and component).  Proposing (i',j',k') uniformly accepts by
//   A = min(1, exp[ (logw_k' + logp_i'j'k' - logn_i') - (logw_k + logp_ijk - logn_i) ]).
// Accepted states give the mixture conjugate updates:
//   w_k    = fraction of accepted states in component k
//   mu_k   = mean accepted theta among component-k states (theta = muK_k + eta_ij)
//   Omega  = pooled within-component variance of accepted theta (shared, scalar)
// and the shared additive/proportional residual from the pooled accepted SS.
// muK is the per-component typical value of the mixed eta (paper's mu_k); w the
// current mixture weights.  errType 0=additive, 1=proportional, 6=lognormal.
//[[Rcpp::export]]
List rpemMstepMix(NumericMatrix muK, NumericVector w, IntegerMatrix etaForComp,
                  int errType, double addSd0, double propSd0, int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepMixDraw before rpemMstepMix");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, K = rpemOp.nMix, nEta = rpemOp.nEta;
  int nParam = muK.nrow();                       // one mix() call (mixed parameter) per row
  if (muK.ncol() != K || (int)w.size() != K) stop("muK/w must have K components");
  if (etaForComp.nrow() != nParam || etaForComp.ncol() != K)
    stop("etaForComp must be nParam x K");
  bool doComb = (errType == 2);                  // combined add+prop: no closed form
  bool doPow = (errType == 4);                   // power (scale*cp^exponent)^2: profile
  bool doTbs = (errType == 3);                   // transform-both-sides + dynamic lambda
  std::vector<long> counts3((doComb || doPow || doTbs) ? (size_t)nsub * nG * K : 0, 0);  // per-(i,j,k)

  std::vector<double> logw(K);
  for (int k = 0; k < K; ++k) logw[k] = log(w[k]);
  // Per-subject mixture log n_i from the stored per-component log p.
  std::vector<double> logn(nsub);
  for (int i = 0; i < nsub; ++i) {
    double mxk = R_NegInf; std::vector<double> lnik(K);
    for (int k = 0; k < K; ++k) {
      double mx = R_NegInf;
      for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[((size_t)i * nG + j) * K + k]; if (v > mx) mx = v; }
      double s = 0.0;
      for (int j = 0; j < nG; ++j) s += exp(rpemOp.logp[((size_t)i * nG + j) * K + k] - mx);
      lnik[k] = mx + log(s) - log((double)nG);
      double wk = logw[k] + lnik[k]; if (wk > mxk) mxk = wk;
    }
    double sk = 0.0;
    for (int k = 0; k < K; ++k) sk += exp(logw[k] + lnik[k] - mxk);
    logn[i] = mxk + log(sk);
  }

  int total = nTrials + burn;
  size_t nU = (size_t)4 * total;                 // proposal i', j', k', accept
  NumericVector mu0(1); NumericMatrix s0(1, 1); s0(0, 0) = 1.0;
  NumericVector lo(1, R_NegInf), hi(1, R_PosInf);
  SEXP zr = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create((int)nU)),
                                wrap(mu0), wrap(s0), wrap(lo), wrap(hi),
                                wrap(IntegerVector::create(1)),
                                wrap(LogicalVector::create(false)),
                                wrap(LogicalVector::create(false)),
                                wrap(NumericVector::create(0.4)),
                                wrap(NumericVector::create(2.05)),
                                wrap(NumericVector::create(1e-10)),
                                wrap(IntegerVector::create(100)));
  NumericMatrix zm(zr);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(zm[(int)k], 0.0, 1.0, 1, 0);

  int ci = 0, cj = 0, ck = 0;
  double clogp = rpemOp.logp[0];
  std::vector<double> sumTk((size_t)nParam * K, 0.0), sumTTk((size_t)nParam * K, 0.0);
  std::vector<long> countK(K, 0);
  double sumSS = 0.0; long sumNobs = 0, m = 0, naccept = 0;
  for (int t = 0; t < total; ++t) {
    double u1 = U[(size_t)4 * t], u2 = U[(size_t)4 * t + 1],
           u3 = U[(size_t)4 * t + 2], u4 = U[(size_t)4 * t + 3];
    int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
    int pjh = (int)(u2 * nG);   if (pjh >= nG)   pjh = nG - 1;
    int pkh = (int)(u3 * K);    if (pkh >= K)    pkh = K - 1;
    double plogp = rpemOp.logp[((size_t)pih * nG + pjh) * K + pkh];
    double logA = (logw[pkh] + plogp - logn[pih]) - (logw[ck] + clogp - logn[ci]);
    if (log(u4) < logA) { ci = pih; cj = pjh; ck = pkh; clogp = plogp; ++naccept; }
    if (t >= burn) {
      // one latent class label ck governs every mixed parameter; accumulate each
      // parameter's per-component typical value using that parameter's active eta
      // (shared-eta -> one eta per parameter; split-ETA -> its own eta per component).
      size_t etaBase = ((size_t)ci * nG + cj) * nEta;
      for (int p = 0; p < nParam; ++p) {
        double eta = rpemOp.etaS[etaBase + etaForComp(p, ck)];
        double theta = muK(p, ck) + eta;
        sumTk[(size_t)p * K + ck] += theta;
        sumTTk[(size_t)p * K + ck] += theta * theta;
      }
      ++countK[ck];
      if (doComb || doPow || doTbs) counts3[((size_t)ci * nG + cj) * K + ck]++;
      int nobsi = rpemOp.nobs[ci];
      size_t r = ((size_t)ci * nG + cj) * K + ck;
      if (errType == 6) {
        long ob = (rpemOp.sampObsOff[ci] + (long)cj * nobsi) * K + (long)ck * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
          if (cp > 0.0 && dv > 0.0) { double lr = log(dv) - log(cp); sumSS += lr * lr; }
        }
      } else {
        sumSS += (errType == 1) ? rpemOp.wssv[r] : rpemOp.ssv[r];
      }
      sumNobs += nobsi; ++m;
    }
  }
  NumericMatrix muNew(nParam, K); NumericVector wNew(K), omegaNew(nEta);
  std::vector<double> sumSSk(nEta, 0.0);   // within-component SS grouped by active eta
  std::vector<long> countEta(nEta, 0);
  for (int k = 0; k < K; ++k) {
    for (int p = 0; p < nParam; ++p) {
      if (countK[k] > 0) {
        muNew(p, k) = sumTk[(size_t)p * K + k] / (double)countK[k];
        int a = etaForComp(p, k);        // components sharing an eta pool into its Sigma
        sumSSk[a] += sumTTk[(size_t)p * K + k] - (double)countK[k] * muNew(p, k) * muNew(p, k);
        countEta[a] += countK[k];
      } else {
        muNew(p, k) = muK(p, k);         // component starved this iteration: hold
      }
    }
    wNew[k] = (double)countK[k] / (double)m;
  }
  for (int a = 0; a < nEta; ++a)
    omegaNew[a] = (countEta[a] > 0) ? sumSSk[a] / (double)countEta[a] : 0.0;
  // Empty/collapsing-component guard (design/rpem/07): a weight of exactly 0 makes
  // logw = -Inf next iteration, permanently killing the component (it can never be
  // accepted back).  Floor each weight at a small share and renormalize so every
  // component stays reachable; a genuinely empty component then decays to the floor.
  double wFloor = 1e-3, wsum = 0.0;
  for (int k = 0; k < K; ++k) { if (wNew[k] < wFloor) wNew[k] = wFloor; wsum += wNew[k]; }
  for (int k = 0; k < K; ++k) wNew[k] /= wsum;
  // Residual: additive/proportional/lognormal use the pooled SS; combined (add+prop)
  // has no closed form, so guarded-Newton-maximize the shared (a=add^2, b=prop^2)
  // Gaussian log-likelihood over the accepted (i,j,k) states' stored cp/dv (mixture
  // stride) -- the same optimizer as the non-mixture combined M-step.
  double addNew, propNew = NA_REAL, powNew = NA_REAL, lamNew = NA_REAL;
  if (doComb) {
    const double eps = 1e-12;
    auto over = [&](double aa, double bb, double &Q, double *ga, double *gb,
                    double *Haa, double *Hab, double *Hbb) {
      Q = 0.0; if (ga) { *ga = *gb = *Haa = *Hab = *Hbb = 0.0; }
      for (int i = 0; i < nsub; ++i) {
        int nobsi = rpemOp.nobs[i];
        for (int j = 0; j < nG; ++j) for (int k = 0; k < K; ++k) {
          long c = counts3[((size_t)i * nG + j) * K + k]; if (!c) continue;
          long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
          for (int o = 0; o < nobsi; ++o) {
            double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
            double wcp = cp * cp, r2 = rr * rr, V = aa + bb * wcp; if (V < eps) V = eps;
            double invV = 1.0 / V;
            Q += (double)c * (-0.5 * log(V) - 0.5 * r2 * invV);
            if (ga) {
              double g = -0.5 * invV + 0.5 * r2 * invV * invV;
              double h = 0.5 * invV * invV - r2 * invV * invV * invV;
              *ga += c * g; *gb += c * wcp * g;
              *Haa += c * h; *Hab += c * wcp * h; *Hbb += c * wcp * wcp * h;
            }
          }
        }
      }
    };
    double a = addSd0 * addSd0, bpar = propSd0 * propSd0, Qcur; over(a, bpar, Qcur, 0,0,0,0,0);
    for (int it = 0; it < 200; ++it) {
      double Qtmp, ga, gb, Haa, Hab, Hbb;
      over(a, bpar, Qtmp, &ga, &gb, &Haa, &Hab, &Hbb);
      if (fabs(ga) + fabs(gb) < 1e-9) break;
      double det = Haa * Hbb - Hab * Hab, da, db;
      if (Haa < 0.0 && det > 0.0) { da = -(Hbb * ga - Hab * gb) / det; db = -(-Hab * ga + Haa * gb) / det; }
      else { double sc = 1.0 / (fabs(Haa) + fabs(Hbb) + 1.0); da = sc * ga; db = sc * gb; }
      double step = 1.0; bool ok = false;
      for (int ls = 0; ls < 40; ++ls) {
        double an = a + step * da, bn = bpar + step * db;
        if (an < eps) an = eps; if (bn < eps) bn = eps;
        double Qn; over(an, bn, Qn, 0,0,0,0,0);
        if (Qn > Qcur) { a = an; bpar = bn; Qcur = Qn; ok = true; break; }
        step *= 0.5;
      }
      if (!ok) break;
    }
    addNew = sqrt(a); propNew = sqrt(bpar);
  } else if (doPow) {
    // power residual V = (scale * cp^c)^2: the scale profiles out (scale^2 = SSc/N,
    // SSc = sum count r^2 / cp^(2c)), so golden-section the exponent c over the accepted
    // (i,j,k) states -- same profile as the non-mixture power M-step.  addNew holds the
    // scale (its theta slot), powNew the exponent.
    auto stat = [&](double cc, double &SSc, double &SumLogCp) {
      SSc = 0.0; SumLogCp = 0.0;
      for (int i = 0; i < nsub; ++i) {
        int nobsi = rpemOp.nobs[i];
        for (int j = 0; j < nG; ++j) for (int k = 0; k < K; ++k) {
          long c = counts3[((size_t)i * nG + j) * K + k]; if (!c) continue;
          long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
          for (int o = 0; o < nobsi; ++o) {
            double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp, acp = fabs(cp) + 1e-300;
            SSc += (double)c * rr * rr / pow(acp, 2.0 * cc);
            SumLogCp += (double)c * log(acp);
          }
        }
      }
    };
    double N = (double)sumNobs;
    auto f = [&](double cc) { double SSc, SL; stat(cc, SSc, SL);
      if (SSc < 1e-300) SSc = 1e-300; return -0.5 * (N * log(SSc / N) + 2.0 * cc * SL); };
    powNew = rpemGolden(f, 0.0, 3.0, 100);
    double SSc, SL; stat(powNew, SSc, SL);
    addNew = sqrt(SSc / N);
  } else if (doTbs) {
    // transform-both-sides: additive on the transformed scale + dynamic lambda.  The
    // transformed-scale variance profiles out (add^2 = SS(lambda)/N), so golden-section
    // the profile loglik -0.5 N log(SS/N) + sum log|dt/dDV| over the accepted (i,j,k)
    // states -- same profile as the non-mixture TBS M-step (per-obs yj/low/hi from the
    // E-step).  addNew holds the transformed-scale add.sd, lamNew the lambda.
    auto ssJac = [&](double lam, double &SS, double &Jac) {
      SS = 0.0; Jac = 0.0;
      for (int i = 0; i < nsub; ++i) {
        int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
        for (int j = 0; j < nG; ++j) for (int k = 0; k < K; ++k) {
          long c = counts3[((size_t)i * nG + j) * K + k]; if (!c) continue;
          long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
          for (int o = 0; o < nobsi; ++o) {
            double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
            int yj = rpemOp.yjv[so0 + o]; double low = rpemOp.lowv[so0 + o], hi = rpemOp.hiv[so0 + o];
            double d = _powerD(dv, lam, yj, low, hi) - _powerD(cp, lam, yj, low, hi);
            SS += (double)c * d * d;
            Jac += (double)c * log(fabs(_powerDD(dv, lam, yj, low, hi)) + 1e-300);
          }
        }
      }
    };
    double N = (double)sumNobs;
    auto f = [&](double lam) { double SS, Jac; ssJac(lam, SS, Jac);
      if (SS < 1e-300) SS = 1e-300; return -0.5 * N * log(SS / N) + Jac; };
    lamNew = rpemGolden(f, -2.0, 3.0, 100);
    double SS, Jac; ssJac(lamNew, SS, Jac);
    addNew = sqrt(SS / N);
  } else {
    addNew = sqrt(sumSS / (double)sumNobs);
  }
  return List::create(_["muK"] = muNew, _["omega"] = omegaNew, _["w"] = wNew,
                      _["addSd"] = addNew, _["propSd"] = propNew, _["power"] = powNew,
                      _["lambda"] = lamNew, _["accept"] = (double)naccept / total);
}

// K=1 M-step with a covariate design (mu2, D22).  Generalizes rpemMstepK1's mu
// update to a weighted linear regression of the accepted theta samples on the
// per-subject design matrix (nEta==1).  design: nsub x nCoef with row i =
// [1, cov1_i, cov2_i, ...]; coefs: current [typical, covCoef...] used to
// reconstruct theta_ij = design_i . coefs + eta_ij.  Returns new coefs, omega
// (residual variance of the regression), the additive add.sd, and accept rate.
// The intercept-only case (nCoef==1, design all ones) reduces to rpemMstepK1.
// errType 0 = additive residual, sd is add.sd; errType 1 = proportional, sd is
// prop.sd times cp.  The residual param is updated in closed form from the stored
// per-sample SS for additive or WSS for proportional: sqrt(sum acc / sum nobs).
//[[Rcpp::export]]
List rpemMstepK1Reg(NumericMatrix design, NumericVector coefs, int errType,
                    int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Reg");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Reg currently supports nEta==1");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  int nCoef = design.ncol();
  if (design.nrow() != nsub) stop("design must have one row per subject");
  if ((int)coefs.size() != nCoef) stop("coefs length must match design columns");

  std::vector<double> logn(nsub);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) s += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    logn[i] = mx + log(s) - log((double)nG);
  }
  // per-subject current linear predictor
  std::vector<double> muLin(nsub, 0.0), dmat((size_t)nsub * nCoef);
  for (int i = 0; i < nsub; ++i) {
    double s = 0.0;
    for (int k = 0; k < nCoef; ++k) { double x = design(i, k); dmat[(size_t)i * nCoef + k] = x; s += x * coefs[k]; }
    muLin[i] = s;
  }

  int total = nTrials + burn;
  size_t nU = (size_t)3 * total;
  NumericVector mu0(1); NumericMatrix s0(1, 1); s0(0, 0) = 1.0;
  NumericVector lo(1, R_NegInf), hi(1, R_PosInf);
  SEXP zr = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create((int)nU)),
                                wrap(mu0), wrap(s0), wrap(lo), wrap(hi),
                                wrap(IntegerVector::create(1)),
                                wrap(LogicalVector::create(false)),
                                wrap(LogicalVector::create(false)),
                                wrap(NumericVector::create(0.4)),
                                wrap(NumericVector::create(2.05)),
                                wrap(NumericVector::create(1e-10)),
                                wrap(IntegerVector::create(100)));
  NumericMatrix zm(zr);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(zm[(int)k], 0.0, 1.0, 1, 0);

  int ci = 0, cj = 0;
  double clogp = rpemOp.logp[0];
  std::vector<double> XtX((size_t)nCoef * nCoef, 0.0), Xty(nCoef, 0.0);
  double Stt = 0.0, sumSS = 0.0; long sumNobs = 0, m = 0, naccept = 0;
  for (int t = 0; t < total; ++t) {
    double u1 = U[(size_t)3 * t], u2 = U[(size_t)3 * t + 1], u3 = U[(size_t)3 * t + 2];
    int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
    int pjh = (int)(u2 * nG);   if (pjh >= nG)   pjh = nG - 1;
    double plogp = rpemOp.logp[(size_t)pih * nG + pjh];
    double logA = (plogp - clogp) + (logn[ci] - logn[pih]);
    if (log(u3) < logA) { ci = pih; cj = pjh; clogp = plogp; ++naccept; }
    if (t >= burn) {
      double theta = muLin[ci] + rpemOp.etaS[(size_t)ci * nG + cj];
      const double *xi = &dmat[(size_t)ci * nCoef];
      for (int a = 0; a < nCoef; ++a) {
        Xty[a] += xi[a] * theta;
        for (int b = 0; b < nCoef; ++b) XtX[(size_t)a * nCoef + b] += xi[a] * xi[b];
      }
      Stt += theta * theta;
      int nobsi = rpemOp.nobs[ci];
      size_t r = (size_t)ci * nG + cj;
      if (errType == 6) {                        // lognormal: additive on log scale
        long ob = rpemOp.sampObsOff[ci] + (long)cj * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
          if (cp > 0.0 && dv > 0.0) { double lr = log(dv) - log(cp); sumSS += lr * lr; }
        }
      } else {
        sumSS += (errType == 1) ? rpemOp.wssv[r] : rpemOp.ssv[r];
      }
      sumNobs += nobsi;
      ++m;
    }
  }
  arma::mat A(XtX.data(), nCoef, nCoef);  // symmetric, so column/row-major agree
  arma::vec b(Xty.data(), nCoef);
  arma::vec betaNew = arma::solve(A, b, arma::solve_opts::likely_sympd);
  double omegaNew = (Stt - arma::dot(betaNew, b)) / (double)m;
  double addNew = sqrt(sumSS / (double)sumNobs);
  NumericVector coefOut(nCoef);
  for (int a = 0; a < nCoef; ++a) coefOut[a] = betaNew[a];
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["addSd"] = addNew, _["accept"] = (double)naccept / total);
}

// Number of residual score parameters for a given errType (design/rpem/08).
//   0 add / 1 prop / 6 lnorm : single sd
//   2 combined (add.sd, prop.sd) ; 3 TBS (add.sd, lambda) ; 4 power (prop.sd, power)
static inline int rpemNResPar(int errType) {
  return (errType == 2 || errType == 3 || errType == 4) ? 2 : 1;
}

// Unweighted residual-parameter score for one stored sample (i,j), written into
// sRes (length rpemNResPar(errType)).  For add/prop the closed-form uses the
// precomputed SS/WSS; combined and power have no closed form, so re-score per obs
// from the stored structural prediction cp (rpemOp.cpv) and observation dv
// (rpemOp.dvv).  Score of a Gaussian residual param theta with per-obs variance V:
//   d/dtheta = sum_obs (r^2 - V)/(2 V^2) * dV/dtheta ,  r = dv - cp.
static inline void rpemResidScoreSample(int errType, int i, int j,
                                        const double *resPar, double *sRes) {
  int nG = rpemOp.nGauss, nobsi = rpemOp.nobs[i];
  size_t r = (size_t)i * nG + j;
  if (errType == 2 || errType == 4) {
    long off = rpemOp.sampObsOff[i] + (long)j * nobsi;
    double d0 = 0.0, d1 = 0.0;
    if (errType == 2) {                       // combined: V = add^2 + prop^2 cp^2
      double add = resPar[0], prop = resPar[1], add2 = add * add, prop2 = prop * prop;
      for (int o = 0; o < nobsi; ++o) {
        double cp = rpemOp.cpv[off + o], dv = rpemOp.dvv[off + o];
        double cp2 = cp * cp, V = add2 + prop2 * cp2, rr = dv - cp;
        double g = (rr * rr - V) / (V * V);
        d0 += add * g;                        // d/d add.sd
        d1 += prop * cp2 * g;                 // d/d prop.sd
      }
    } else {                                  // power: V = b^2 cp^(2c)
      double b = resPar[0], c = resPar[1];
      for (int o = 0; o < nobsi; ++o) {
        double cp = rpemOp.cpv[off + o], dv = rpemOp.dvv[off + o];
        double cp2c = pow(cp, 2.0 * c), V = b * b * cp2c, rr = dv - cp;
        double g = (rr * rr - V) / (V * V);
        d0 += b * cp2c * g;                   // d/d prop.sd (b)
        d1 += V * log(cp) * g;                // d/d power (c)
      }
    }
    sRes[0] = d0; sRes[1] = d1;
    return;
  }
  if (errType == 3) {                         // TBS: add.sd (transformed scale) + lambda
    double sd = resPar[0], lam = resPar[1], sd2 = sd * sd, sd3 = sd2 * sd;
    long off = rpemOp.sampObsOff[i] + (long)j * nobsi, so0 = rpemOp.idS[i];
    double h = 1e-4 * (fabs(lam) + 1.0);      // central FD step for the lambda score
    double SSt = 0.0, Lp = 0.0, Lm = 0.0;
    for (int o = 0; o < nobsi; ++o) {
      double cp = rpemOp.cpv[off + o], dv = rpemOp.dvv[off + o];
      int yj = rpemOp.yjv[so0 + o]; double low = rpemOp.lowv[so0 + o], hi = rpemOp.hiv[so0 + o];
      double d = _powerD(dv, lam, yj, low, hi) - _powerD(cp, lam, yj, low, hi);
      SSt += d * d;
      // lambda score by central difference of the lambda-dependent residual + Jacobian
      // loglik (the -0.5 log(2pi sd^2) constant is lambda-free and drops out):
      double dp = _powerD(dv, lam + h, yj, low, hi) - _powerD(cp, lam + h, yj, low, hi);
      double dm = _powerD(dv, lam - h, yj, low, hi) - _powerD(cp, lam - h, yj, low, hi);
      Lp += -0.5 * dp * dp / sd2 + log(fabs(_powerDD(dv, lam + h, yj, low, hi)) + 1e-300);
      Lm += -0.5 * dm * dm / sd2 + log(fabs(_powerDD(dv, lam - h, yj, low, hi)) + 1e-300);
    }
    sRes[0] = -(double)nobsi / sd + SSt / sd3; // d/d add.sd
    sRes[1] = (Lp - Lm) / (2.0 * h);           // d/d lambda
    return;
  }
  double sd = resPar[0], sd3 = sd * sd * sd;  // add (0) / prop (1) / lnorm (6)
  double acc = (errType == 1) ? rpemOp.wssv[r] : rpemOp.ssv[r];
  sRes[0] = -(double)nobsi / sd + acc / sd3;
}

// Fisher-score information for the K=1 regression case (design/rpem/08).  At the
// converged estimates, each subject's marginal-likelihood score is formed by the
// Fisher identity s_i = E_{eta|Y_i}[grad complete-data loglik], with the posterior
// expectation taken over the stored samples via the self-normalized importance
// weights w_ij = softmax_j logp_ij (the samples were drawn from the prior N(0,omega),
// so these weights ARE the posterior).  Score blocks per subject:
//   coef_k  : design_i,k * (theta_ij - design_i.coefs)/omega = design_i,k * eta_ij/omega
//   sd      : -nobs_i/sd + acc_ij/sd^3   (acc = SS additive / WSS proportional)
//   omega   : -1/(2 omega) + eta_ij^2/(2 omega^2)
// Returns the nsub x (nCoef+2) per-subject score matrix S (columns: coefs, sd,
// omega); the caller forms the empirical Fisher information I = S^T S and inverts it.
//[[Rcpp::export]]
NumericMatrix rpemFisherReg(NumericMatrix design, NumericVector coefs, double omega,
                            int errType, NumericVector resPar) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemFisherReg");
  if (rpemOp.nEta != 1) stop("rpemFisherReg currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, nCoef = design.ncol();
  if (design.nrow() != nsub) stop("design must have one row per subject");
  int nRes = rpemNResPar(errType);
  if ((int)resPar.size() != nRes) stop("resPar length must match errType");
  int p = nCoef + nRes + 1;                     // coefs, residual params, omega
  NumericMatrix S(nsub, p);
  double om2 = omega * omega;
  std::vector<double> sRes(nRes), sResAcc(nRes);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double sw = 0.0;
    for (int j = 0; j < nG; ++j) sw += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    double sEta = 0.0, sEta2 = 0.0;              // importance-weighted score sums
    for (int m = 0; m < nRes; ++m) sResAcc[m] = 0.0;
    for (int j = 0; j < nG; ++j) {
      double w = exp(rpemOp.logp[(size_t)i * nG + j] - mx) / sw;
      double eta = rpemOp.etaS[(size_t)i * nG + j];
      sEta  += w * eta;
      sEta2 += w * eta * eta;
      rpemResidScoreSample(errType, i, j, resPar.begin(), sRes.data());
      for (int m = 0; m < nRes; ++m) sResAcc[m] += w * sRes[m];
    }
    for (int k = 0; k < nCoef; ++k) S(i, k) = design(i, k) * sEta / omega;
    for (int m = 0; m < nRes; ++m) S(i, nCoef + m) = sResAcc[m];
    S(i, nCoef + nRes) = -0.5 / omega + 0.5 * sEta2 / om2;
  }
  return S;
}

// Fisher-score information for the general multi-eta case (diagonal Omega, no
// covariates).  Same Fisher-identity construction as rpemFisherReg, but with one
// typical-value (mu_a) and one variance (om_a) score per random effect a:
//   mu_a  : eta_ij,a / om_a
//   sd    : -nobs_i/sd + acc_ij/sd^3   (acc = SS additive / WSS proportional)
//   om_a  : -1/(2 om_a) + eta_ij,a^2/(2 om_a^2)
// Returns the nsub x (2*nEta+1) per-subject score matrix S with columns ordered
// [mu_1..nEta, sd, om_1..nEta]; the caller forms I = S^T S and inverts it.
//[[Rcpp::export]]
NumericMatrix rpemFisherDiag(NumericVector muVec, NumericVector omVec,
                             int errType, NumericVector resPar) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemFisherDiag");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, nEta = rpemOp.nEta;
  if ((int)muVec.size() != nEta || (int)omVec.size() != nEta)
    stop("muVec/omVec must have nEta entries");
  int nRes = rpemNResPar(errType);
  if ((int)resPar.size() != nRes) stop("resPar length must match errType");
  int p = 2 * nEta + nRes;                       // mu_1..nEta, residual params, om_1..nEta
  NumericMatrix S(nsub, p);
  std::vector<double> sMu(nEta), sOm(nEta), sRes(nRes), sResAcc(nRes);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double sw = 0.0;
    for (int j = 0; j < nG; ++j) sw += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    for (int a = 0; a < nEta; ++a) { sMu[a] = 0.0; sOm[a] = 0.0; }
    for (int m = 0; m < nRes; ++m) sResAcc[m] = 0.0;
    for (int j = 0; j < nG; ++j) {
      double w = exp(rpemOp.logp[(size_t)i * nG + j] - mx) / sw;
      for (int a = 0; a < nEta; ++a) {
        double eta = rpemOp.etaS[((size_t)i * nG + j) * nEta + a];
        double om = omVec[a];
        sMu[a] += w * eta / om;
        sOm[a] += w * (-0.5 / om + 0.5 * eta * eta / (om * om));
      }
      rpemResidScoreSample(errType, i, j, resPar.begin(), sRes.data());
      for (int m = 0; m < nRes; ++m) sResAcc[m] += w * sRes[m];
    }
    for (int a = 0; a < nEta; ++a) S(i, a) = sMu[a];
    for (int m = 0; m < nRes; ++m) S(i, nEta + m) = sResAcc[m];
    for (int a = 0; a < nEta; ++a) S(i, nEta + nRes + a) = sOm[a];
  }
  return S;
}

// K=1 combined-error M-step (add + prop, errType 2).  Same regression MH as
// rpemMstepK1Reg for the structural coefs / omega, but the residual has no closed
// form: obs variance V = a + b*cp^2 with a=add.sd^2, b=prop.sd^2.  We accumulate
// per-(i,j) MH visit counts, then Newton-maximize the Gaussian log-likelihood
// sum_visited count * sum_obs [-0.5 log V - r^2/(2V)] over (a,b) using the stored
// per-obs cp^2 (rpemOp.cp2v) and r^2 (rpemOp.r2v).  addSd0/propSd0 seed the Newton.
//[[Rcpp::export]]
List rpemMstepK1Comb(NumericMatrix design, NumericVector coefs, double addSd0,
                     double propSd0, int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Comb");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Comb currently supports nEta==1");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  int nCoef = design.ncol();
  if (design.nrow() != nsub) stop("design must have one row per subject");
  if ((int)coefs.size() != nCoef) stop("coefs length must match design columns");

  std::vector<double> logn(nsub);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) s += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    logn[i] = mx + log(s) - log((double)nG);
  }
  std::vector<double> muLin(nsub, 0.0), dmat((size_t)nsub * nCoef);
  for (int i = 0; i < nsub; ++i) {
    double s = 0.0;
    for (int k = 0; k < nCoef; ++k) { double x = design(i, k); dmat[(size_t)i * nCoef + k] = x; s += x * coefs[k]; }
    muLin[i] = s;
  }

  int total = nTrials + burn;
  size_t nU = (size_t)3 * total;
  NumericVector mu0(1); NumericMatrix s0(1, 1); s0(0, 0) = 1.0;
  NumericVector lo(1, R_NegInf), hi(1, R_PosInf);
  SEXP zr = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create((int)nU)),
                                wrap(mu0), wrap(s0), wrap(lo), wrap(hi),
                                wrap(IntegerVector::create(1)),
                                wrap(LogicalVector::create(false)),
                                wrap(LogicalVector::create(false)),
                                wrap(NumericVector::create(0.4)),
                                wrap(NumericVector::create(2.05)),
                                wrap(NumericVector::create(1e-10)),
                                wrap(IntegerVector::create(100)));
  NumericMatrix zm(zr);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(zm[(int)k], 0.0, 1.0, 1, 0);

  int ci = 0, cj = 0;
  double clogp = rpemOp.logp[0];
  std::vector<double> XtX((size_t)nCoef * nCoef, 0.0), Xty(nCoef, 0.0);
  std::vector<long> counts((size_t)nsub * nG, 0);
  double Stt = 0.0; long m = 0, naccept = 0;
  for (int t = 0; t < total; ++t) {
    double u1 = U[(size_t)3 * t], u2 = U[(size_t)3 * t + 1], u3 = U[(size_t)3 * t + 2];
    int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
    int pjh = (int)(u2 * nG);   if (pjh >= nG)   pjh = nG - 1;
    double plogp = rpemOp.logp[(size_t)pih * nG + pjh];
    double logA = (plogp - clogp) + (logn[ci] - logn[pih]);
    if (log(u3) < logA) { ci = pih; cj = pjh; clogp = plogp; ++naccept; }
    if (t >= burn) {
      double theta = muLin[ci] + rpemOp.etaS[(size_t)ci * nG + cj];
      const double *xi = &dmat[(size_t)ci * nCoef];
      for (int a = 0; a < nCoef; ++a) {
        Xty[a] += xi[a] * theta;
        for (int b = 0; b < nCoef; ++b) XtX[(size_t)a * nCoef + b] += xi[a] * xi[b];
      }
      Stt += theta * theta;
      counts[(size_t)ci * nG + cj]++;
      ++m;
    }
  }
  arma::mat A(XtX.data(), nCoef, nCoef);
  arma::vec bb(Xty.data(), nCoef);
  arma::vec betaNew = arma::solve(A, bb, arma::solve_opts::likely_sympd);
  double omegaNew = (Stt - arma::dot(betaNew, bb)) / (double)m;

  // Guarded optimization of (a=add^2, b=prop^2) over the visited (i,j).  The
  // combined-variance Gaussian log-likelihood is not globally concave, so raw
  // Newton diverges; we take a Newton direction only when the local Hessian is
  // negative-definite, otherwise gradient ascent, and accept a step only if it
  // increases the objective Q (backtracking line search keeps a,b > 0).
  const double eps = 1e-12;
  auto qFun = [&](double aa, double bb) -> double {
    double Q = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j];
        if (c == 0) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double V = aa + bb * cp * cp; if (V < eps) V = eps;
          Q += (double)c * (-0.5 * log(V) - 0.5 * rr * rr / V);
        }
      }
    }
    return Q;
  };
  double a = addSd0 * addSd0, bpar = propSd0 * propSd0;
  double Qcur = qFun(a, bpar);
  for (int it = 0; it < 200; ++it) {
    double ga = 0, gb = 0, Haa = 0, Hab = 0, Hbb = 0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j];
        if (c == 0) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double w = cp * cp;                  // cp^2
          double r2 = rr * rr;
          double V = a + bpar * w; if (V < eps) V = eps;
          double invV = 1.0 / V;
          double g = -0.5 * invV + 0.5 * r2 * invV * invV;         // df/da
          double h = 0.5 * invV * invV - r2 * invV * invV * invV;  // d2f/da2
          ga += c * g;      gb += c * w * g;
          Haa += c * h;     Hab += c * w * h;     Hbb += c * w * w * h;
        }
      }
    }
    double gnorm = fabs(ga) + fabs(gb);
    if (gnorm < 1e-9) break;
    double det = Haa * Hbb - Hab * Hab, da, db;
    if (Haa < 0.0 && det > 0.0) {             // Hessian negative-definite: Newton
      da = -(Hbb * ga - Hab * gb) / det;
      db = -(-Hab * ga + Haa * gb) / det;
    } else {                                  // fall back to gradient ascent
      double sc = 1.0 / (fabs(Haa) + fabs(Hbb) + 1.0);
      da = sc * ga; db = sc * gb;
    }
    double step = 1.0; bool ok = false;
    for (int ls = 0; ls < 40; ++ls) {
      // Project onto the feasible region (a,b >= eps) rather than rejecting a step
      // that violates the boundary, so the optimizer does not stall when one
      // component sits on the boundary.
      double an = a + step * da, bn = bpar + step * db;
      if (an < eps) an = eps;
      if (bn < eps) bn = eps;
      double Qn = qFun(an, bn);
      if (Qn > Qcur) { a = an; bpar = bn; Qcur = Qn; ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;
  }

  NumericVector coefOut(nCoef);
  for (int k = 0; k < nCoef; ++k) coefOut[k] = betaNew[k];
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["addSd"] = sqrt(a), _["propSd"] = sqrt(bpar),
                      _["accept"] = (double)naccept / total);
}

// Shared regression M-step core for the numeric residual updates (TBS lambda,
// power): runs the joint Metropolis-Hastings over the stored E-step samples,
// accumulates per-(i,j) visit counts, and returns the regression coefs / omega.
static void rpemMHReg(NumericMatrix design, NumericVector coefs, int nTrials, int burn,
                      std::vector<long> &counts, NumericVector &coefOut,
                      double &omegaOut, double &acceptOut, long &mOut) {
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  int nCoef = design.ncol();
  std::vector<double> logn(nsub);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) s += exp(rpemOp.logp[(size_t)i * nG + j] - mx);
    logn[i] = mx + log(s) - log((double)nG);
  }
  std::vector<double> muLin(nsub, 0.0), dmat((size_t)nsub * nCoef);
  for (int i = 0; i < nsub; ++i) {
    double s = 0.0;
    for (int k = 0; k < nCoef; ++k) { double x = design(i, k); dmat[(size_t)i * nCoef + k] = x; s += x * coefs[k]; }
    muLin[i] = s;
  }
  int total = nTrials + burn;
  size_t nU = (size_t)3 * total;
  NumericVector mu0(1); NumericMatrix s0(1, 1); s0(0, 0) = 1.0;
  NumericVector lo(1, R_NegInf), hi(1, R_PosInf);
  SEXP zr = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create((int)nU)),
                                wrap(mu0), wrap(s0), wrap(lo), wrap(hi),
                                wrap(IntegerVector::create(1)),
                                wrap(LogicalVector::create(false)),
                                wrap(LogicalVector::create(false)),
                                wrap(NumericVector::create(0.4)),
                                wrap(NumericVector::create(2.05)),
                                wrap(NumericVector::create(1e-10)),
                                wrap(IntegerVector::create(100)));
  NumericMatrix zm(zr);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(zm[(int)k], 0.0, 1.0, 1, 0);
  int ci = 0, cj = 0;
  double clogp = rpemOp.logp[0];
  std::vector<double> XtX((size_t)nCoef * nCoef, 0.0), Xty(nCoef, 0.0);
  counts.assign((size_t)nsub * nG, 0);
  double Stt = 0.0; long m = 0, naccept = 0;
  for (int t = 0; t < total; ++t) {
    double u1 = U[(size_t)3 * t], u2 = U[(size_t)3 * t + 1], u3 = U[(size_t)3 * t + 2];
    int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
    int pjh = (int)(u2 * nG);   if (pjh >= nG)   pjh = nG - 1;
    double plogp = rpemOp.logp[(size_t)pih * nG + pjh];
    double logA = (plogp - clogp) + (logn[ci] - logn[pih]);
    if (log(u3) < logA) { ci = pih; cj = pjh; clogp = plogp; ++naccept; }
    if (t >= burn) {
      double theta = muLin[ci] + rpemOp.etaS[(size_t)ci * nG + cj];
      const double *xi = &dmat[(size_t)ci * nCoef];
      for (int a = 0; a < nCoef; ++a) {
        Xty[a] += xi[a] * theta;
        for (int b = 0; b < nCoef; ++b) XtX[(size_t)a * nCoef + b] += xi[a] * xi[b];
      }
      Stt += theta * theta;
      counts[(size_t)ci * nG + cj]++;
      ++m;
    }
  }
  arma::mat A(XtX.data(), nCoef, nCoef);
  arma::vec bb(Xty.data(), nCoef);
  arma::vec betaNew = arma::solve(A, bb, arma::solve_opts::likely_sympd);
  omegaOut = (Stt - arma::dot(betaNew, bb)) / (double)m;
  coefOut = NumericVector(nCoef);
  for (int k = 0; k < nCoef; ++k) coefOut[k] = betaNew[k];
  acceptOut = (double)naccept / total;
  mOut = m;
}

// Golden-section maximization of f over [a,b]; returns argmax.
static double rpemGolden(const std::function<double(double)> &f, double a, double b, int iter) {
  const double gr = (sqrt(5.0) - 1.0) / 2.0;
  double c = b - gr * (b - a), d = a + gr * (b - a);
  double fc = f(c), fd = f(d);
  for (int it = 0; it < iter; ++it) {
    if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = f(c); }
    else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = f(d); }
    if (b - a < 1e-7) break;
  }
  return 0.5 * (a + b);
}

// K=1 TBS (transform-both-sides) M-step: additive error on the transformed scale
// with a dynamic (estimated) Box-Cox / Yeo-Johnson lambda.  The transformed-scale
// residual variance profiles out in closed form (add.sd^2 = SS(lambda)/N), so we
// golden-section maximize the profile log-likelihood
//   f(lambda) = -0.5 N log(SS(lambda)/N) + sum log|dt/dDV|
// where t = _powerD(., lambda, yj, low, high) and the Jacobian uses _powerDD.
//[[Rcpp::export]]
List rpemMstepK1TBS(NumericMatrix design, NumericVector coefs, double addSd0,
                    double lambda0, int yj, double low, double high,
                    int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1TBS");
  if (rpemOp.nEta != 1) stop("rpemMstepK1TBS currently supports nEta==1");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  if (design.nrow() != nsub) stop("design must have one row per subject");
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, counts, coefOut, omegaNew, accept, m);

  long N = 0;
  for (int i = 0; i < nsub; ++i) {
    int nobsi = rpemOp.nobs[i];
    for (int j = 0; j < nG; ++j) { long c = counts[(size_t)i * nG + j]; if (c) N += c * (long)nobsi; }
  }
  auto ssJac = [&](double lam, double &SS, double &Jac) {
    SS = 0.0; Jac = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
          double d = _powerD(dv, lam, yj, low, high) - _powerD(cp, lam, yj, low, high);
          SS += (double)c * d * d;
          Jac += (double)c * log(fabs(_powerDD(dv, lam, yj, low, high)) + 1e-300);
        }
      }
    }
  };
  auto f = [&](double lam) { double SS, Jac; ssJac(lam, SS, Jac);
    if (SS < 1e-300) SS = 1e-300; return -0.5 * (double)N * log(SS / (double)N) + Jac; };
  double lamHat = rpemGolden(f, -2.0, 3.0, 100);
  double SS, Jac; ssJac(lamHat, SS, Jac);
  double addNew = sqrt(SS / (double)N);
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["addSd"] = addNew, _["lambda"] = lamHat, _["accept"] = accept);
}

// K=1 power-error M-step: variance V = (propSd * cp^power)^2 with an estimated
// exponent.  The scale profiles out (propSd^2 = SSc/N with SSc = sum (DV-cp)^2 /
// cp^(2c)), so we golden-section maximize the profile log-likelihood
//   f(c) = -0.5 [ N log(SSc/N) + 2 c sum log|cp| ]  over the exponent c.
//[[Rcpp::export]]
List rpemMstepK1Pow(NumericMatrix design, NumericVector coefs, double propSd0,
                    double power0, int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Pow");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Pow currently supports nEta==1");
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  if (design.nrow() != nsub) stop("design must have one row per subject");
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, counts, coefOut, omegaNew, accept, m);

  long N = 0;
  for (int i = 0; i < nsub; ++i) {
    int nobsi = rpemOp.nobs[i];
    for (int j = 0; j < nG; ++j) { long c = counts[(size_t)i * nG + j]; if (c) N += c * (long)nobsi; }
  }
  auto stat = [&](double cc, double &SSc, double &SumLogCp) {
    SSc = 0.0; SumLogCp = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double acp = fabs(cp) + 1e-300;
          SSc += (double)c * rr * rr / pow(acp, 2.0 * cc);
          SumLogCp += (double)c * log(acp);
        }
      }
    }
  };
  auto f = [&](double cc) { double SSc, SL; stat(cc, SSc, SL);
    if (SSc < 1e-300) SSc = 1e-300; return -0.5 * ((double)N * log(SSc / (double)N) + 2.0 * cc * SL); };
  double cHat = rpemGolden(f, 0.0, 3.0, 100);
  double SSc, SL; stat(cHat, SSc, SL);
  double propNew = sqrt(SSc / (double)N);
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["propSd"] = propNew, _["power"] = cHat, _["accept"] = accept);
}

// K=1 censored residual M-step (single endpoint, additive or proportional error
// with M2/M3/M4 BLQ observations).  The naive sum-of-squares is biased under
// censoring, so the residual scale is estimated by maximizing the censored
// log-likelihood over the accepted MH samples: for each candidate sd the
// observed observations contribute the Gaussian density and the censored ones the
// CENS probability via doCensNormal1.  1-D golden-section over sd (>0).
// errType 0 = additive (sd_obs = sd), 1 = proportional (sd_obs = sd*|cp|).
//[[Rcpp::export]]
List rpemMstepK1Cens(NumericMatrix design, NumericVector coefs, int errType,
                     double sd0, int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Cens");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Cens currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, counts, coefOut, omegaNew, accept, m);

  auto Q = [&](double sd) -> double {
    double q = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          long so = so0 + o;
          double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
          double sdo = (errType == 1) ? sd * (fabs(cp) + 1e-300) : sd;
          double r = sdo * sdo;
          double gauss = -M_LN_SQRT_2PI - log(sdo) - 0.5 * (dv - cp) * (dv - cp) / r;
          q += (double)c * doCensNormal1((double)rpemOp.censv[so], dv, rpemOp.limv[so],
                                         gauss, cp, r, 0);
        }
      }
    }
    return q;
  };
  std::function<double(double)> f = [&](double sd) { return Q(sd); };
  double sdHat = rpemGolden(f, 1e-4, std::max(5.0 * sd0, 1.0), 120);
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["addSd"] = sdHat, _["accept"] = accept);
}

// Numeric M-step for non-mu-referenced structural fixed effects (the paper's
// beta): parameters with no random effect, which the conjugate mu update cannot
// move (their sampled eta is identically 0).  Maximizes the importance-weighted
// complete-data log-likelihood Q(beta) = sum_i sum_j w_ij log p(Y_i | beta,
// eta_ij), with w_ij the E-step self-normalized weights and eta_ij the stored
// samples -- so this re-solves p(Y_i | .) for candidate beta (the etas are held).
// Takes one damped diagonal-Newton step (finite-difference gradient/curvature)
// with a backtracking line search; the outer EM loop converges beta over
// iterations.  Must be called while the E-step solve struct is still loaded
// (before the MH residual step, whose rxRmvn draw clobbers it).
//[[Rcpp::export]]
NumericVector rpemMstepBeta(NumericVector base, IntegerVector etaIdx,
                            IntegerVector structIdx, NumericVector struct0) {
  if (!rpemOp.loaded) stop("run rpemEstepK1Draw before rpemMstepBeta");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, nEta = rpemOp.nEta;
  int nB = structIdx.size();
  // self-normalized importance weights w_ij from the stored log p (E-step).
  std::vector<double> w((size_t)nsub * nG);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) { double e = exp(rpemOp.logp[(size_t)i * nG + j] - mx); w[(size_t)i * nG + j] = e; s += e; }
    for (int j = 0; j < nG; ++j) w[(size_t)i * nG + j] /= s;
  }
  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int t = 0; t < rpemOp.ntheta; ++t) baseBuf[t] = base[t];
  std::vector<double> row(rpemOp.ntheta);
  std::vector<double> beta(nB); for (int k = 0; k < nB; ++k) beta[k] = struct0[k];
  auto Q = [&](const std::vector<double> &b) -> double {
    double q = 0.0;
    for (int i = 0; i < nsub; ++i) {
      for (int j = 0; j < nG; ++j) {
        for (unsigned int t = 0; t < rpemOp.ntheta; ++t) row[t] = baseBuf[t];
        for (int k = 0; k < nB; ++k) row[structIdx[k]] = b[k];
        for (int a = 0; a < nEta; ++a) row[etaIdx[a]] = rpemOp.etaS[((size_t)i * nG + j) * nEta + a];
        double logp = -rpemSolveSubject(row.data(), i);
        q += w[(size_t)i * nG + j] * logp;
      }
    }
    return q;
  };
  const double h = 1e-3;
  double Q0 = Q(beta);
  std::vector<double> step(nB, 0.0);
  for (int k = 0; k < nB; ++k) {
    std::vector<double> bp = beta, bm = beta;
    bp[k] += h; bm[k] -= h;
    double qp = Q(bp), qm = Q(bm);
    double g = (qp - qm) / (2.0 * h);
    double hess = (qp - 2.0 * Q0 + qm) / (h * h);   // diagonal curvature
    step[k] = (hess < -1e-8) ? -g / hess : g;       // Newton if concave, else gradient
  }
  double t = 1.0; bool ok = false;
  for (int ls = 0; ls < 30; ++ls) {
    std::vector<double> bn = beta;
    for (int k = 0; k < nB; ++k) bn[k] = beta[k] + t * step[k];
    if (Q(bn) > Q0) { beta = bn; ok = true; break; }
    t *= 0.5;
  }
  (void)ok;
  NumericVector out(nB);
  for (int k = 0; k < nB; ++k) out[k] = beta[k];
  return out;
}

// K=1 multiple-endpoint M-step (mirrors SAEM's per-endpoint residual loop).  The
// E-step already computes the joint multi-endpoint log-likelihood, so this shares
// one MH + regression (rpemMHReg) and then updates each endpoint's scale
// separately over the observations that belong to it.  endpt[idS[i]+o] is the
// 0-based endpoint of subject i's o-th observation (in solve order, from R);
// errTypes[b] is 0 additive or 1 proportional for endpoint b.  Returns the new
// per-endpoint residual SDs (add.sd or prop.sd).
//[[Rcpp::export]]
// Guarded (a=add^2,b=prop^2) optimize over endpoint endB's observations (verbatim
// copy of rpemMstepK1Comb's optimizer with an endpoint filter).
static void rpemGuardedComb(const std::vector<long> &counts, const int *endpt,
                            int endB, double a0, double b0, double &aOut, double &bOut) {
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  const double eps = 1e-12;
  auto inB = [&](int i, int o) -> bool { return endB < 0 || endpt[rpemOp.idS[i] + o] == endB; };
  auto qFun = [&](double a, double bp) -> double {
    double Q = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          if (!inB(i, o)) continue;
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double V = a + bp * cp * cp; if (V < eps) V = eps;
          Q += (double)c * (-0.5 * log(V) - 0.5 * rr * rr / V);
        }
      }
    }
    return Q;
  };
  double a = a0, bpar = b0, Qcur = qFun(a, bpar);
  for (int it = 0; it < 200; ++it) {
    double ga = 0, gb = 0, Haa = 0, Hab = 0, Hbb = 0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          if (!inB(i, o)) continue;
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double w = cp * cp, r2 = rr * rr;
          double V = a + bpar * w; if (V < eps) V = eps; double iv = 1.0 / V;
          double g = -0.5 * iv + 0.5 * r2 * iv * iv, hh = 0.5 * iv * iv - r2 * iv * iv * iv;
          ga += c * g; gb += c * w * g; Haa += c * hh; Hab += c * w * hh; Hbb += c * w * w * hh;
        }
      }
    }
    if (fabs(ga) + fabs(gb) < 1e-9) break;
    double det = Haa * Hbb - Hab * Hab, da, db;
    if (Haa < 0.0 && det > 0.0) { da = -(Hbb * ga - Hab * gb) / det; db = -(-Hab * ga + Haa * gb) / det; }
    else { double sc = 1.0 / (fabs(Haa) + fabs(Hbb) + 1.0); da = sc * ga; db = sc * gb; }
    double step = 1.0; bool ok = false;
    for (int ls = 0; ls < 40; ++ls) {
      // Project the trial point onto the feasible region (a,b >= eps) rather than
      // rejecting it -- otherwise, once one component sits on the boundary, a step
      // whose OTHER component would improve Q is thrown away and the optimizer
      // stalls at a non-stationary point.
      double an = a + step * da, bn = bpar + step * db;
      if (an < eps) an = eps;
      if (bn < eps) bn = eps;
      double Qn = qFun(an, bn);
      if (Qn > Qcur) { a = an; bpar = bn; Qcur = Qn; ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;
  }
  aOut = a; bOut = bpar;
}

// Power-error profile over endpoint endB's observations: variance (scale*cp^c)^2,
// scale profiled out (scale^2 = SSc/N), golden-section on the exponent c.
static void rpemGuardedPow(const std::vector<long> &counts, const int *endpt,
                           int endB, double &sclOut, double &powOut) {
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  auto inB = [&](int i, int o) -> bool { return endB < 0 || endpt[rpemOp.idS[i] + o] == endB; };
  long N = 0;
  for (int i = 0; i < nsub; ++i) {
    int nobsi = rpemOp.nobs[i];
    for (int j = 0; j < nG; ++j) { long c = counts[(size_t)i * nG + j]; if (!c) continue;
      for (int o = 0; o < nobsi; ++o) if (inB(i, o)) N += c; }
  }
  auto stat = [&](double cc, double &SSc, double &SumLogCp) {
    SSc = 0.0; SumLogCp = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          if (!inB(i, o)) continue;
          double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp;
          double acp = fabs(cp) + 1e-300;
          SSc += (double)c * rr * rr / pow(acp, 2.0 * cc);
          SumLogCp += (double)c * log(acp);
        }
      }
    }
  };
  std::function<double(double)> f = [&](double cc) { double SSc, SL; stat(cc, SSc, SL);
    if (SSc < 1e-300) SSc = 1e-300; return -0.5 * ((double)N * log(SSc / (double)N) + 2.0 * cc * SL); };
  double cHat = rpemGolden(f, 0.0, 3.0, 100);
  double SSc, SL; stat(cHat, SSc, SL);
  sclOut = sqrt(SSc / (double)N); powOut = cHat;
}

// TBS lambda profile over endpoint endB's observations: additive error on the
// transformed scale, add.sd profiled out (add.sd^2 = SS(lambda)/N), golden-section
// on lambda maximizing -0.5 N log(SS/N) + sum log|dt/dDV|.  Uses the per-obs
// transform code/bounds cached in the E-step (yjv/lowv/hiv) so each endpoint's own
// boxCox/yeoJohnson transform is applied.
static void rpemGuardedTBS(const std::vector<long> &counts, const int *endpt,
                           int endB, double &addOut, double &lamOut) {
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  auto inB = [&](int i, int o) -> bool { return endB < 0 || endpt[rpemOp.idS[i] + o] == endB; };
  long N = 0;
  for (int i = 0; i < nsub; ++i) {
    int nobsi = rpemOp.nobs[i];
    for (int j = 0; j < nG; ++j) { long c = counts[(size_t)i * nG + j]; if (!c) continue;
      for (int o = 0; o < nobsi; ++o) if (inB(i, o)) N += c; }
  }
  auto ssJac = [&](double lam, double &SS, double &Jac) {
    SS = 0.0; Jac = 0.0;
    for (int i = 0; i < nsub; ++i) {
      int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
      for (int j = 0; j < nG; ++j) {
        long c = counts[(size_t)i * nG + j]; if (!c) continue;
        long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
        for (int o = 0; o < nobsi; ++o) {
          if (!inB(i, o)) continue;
          long so = so0 + o; int yj = rpemOp.yjv[so]; double low = rpemOp.lowv[so], hi = rpemOp.hiv[so];
          double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
          double d = _powerD(dv, lam, yj, low, hi) - _powerD(cp, lam, yj, low, hi);
          SS += (double)c * d * d;
          Jac += (double)c * log(fabs(_powerDD(dv, lam, yj, low, hi)) + 1e-300);
        }
      }
    }
  };
  std::function<double(double)> f = [&](double lam) { double SS, Jac; ssJac(lam, SS, Jac);
    if (SS < 1e-300) SS = 1e-300; return -0.5 * (double)N * log(SS / (double)N) + Jac; };
  double lamHat = rpemGolden(f, -2.0, 3.0, 100);
  double SS, Jac; ssJac(lamHat, SS, Jac);
  addOut = sqrt(SS / (double)N); lamOut = lamHat;
}

//[[Rcpp::export]]
List rpemMstepK1Multi(NumericMatrix design, NumericVector coefs, IntegerVector endpt,
                      IntegerVector errTypes, NumericVector add0, NumericVector prop0,
                      int nTrials, int burn) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Multi");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Multi currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  int nEndpt = errTypes.size();
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, counts, coefOut, omegaNew, accept, m);

  // additive/proportional endpoints: closed form over their own observations.
  std::vector<double> sumSS((size_t)nEndpt, 0.0);
  std::vector<long> sumN((size_t)nEndpt, 0);
  for (int i = 0; i < nsub; ++i) {
    int nobsi = rpemOp.nobs[i];
    long base = rpemOp.idS[i];
    for (int j = 0; j < nG; ++j) {
      long c = counts[(size_t)i * nG + j]; if (!c) continue;
      long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
      for (int o = 0; o < nobsi; ++o) {
        int b = endpt[base + o];
        if (errTypes[b] == 2 || errTypes[b] == 3 || errTypes[b] == 4) continue;  // numeric, done below
        double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o], rr = dv - cp;
        double contrib;
        if (errTypes[b] == 6) {                    // lognormal: additive on log scale
          double lr = (cp > 0.0 && dv > 0.0) ? (log(dv) - log(cp)) : 0.0;
          contrib = lr * lr;
        } else if (errTypes[b] == 1) {             // proportional
          contrib = (cp != 0.0) ? (rr / cp) * (rr / cp) : 0.0;
        } else {                                   // additive
          contrib = rr * rr;
        }
        sumSS[b] += (double)c * contrib;
        sumN[b] += c;
      }
    }
  }
  // combined / power endpoints: numeric optimize over just their observations.
  NumericVector sdAdd(nEndpt), sdProp(nEndpt);
  const int *endptp = &endpt[0];
  for (int b = 0; b < nEndpt; ++b) {
    if (errTypes[b] == 2) {
      double aa, bp;
      rpemGuardedComb(counts, endptp, b, add0[b] * add0[b], prop0[b] * prop0[b], aa, bp);
      sdAdd[b] = sqrt(aa); sdProp[b] = sqrt(bp);
    } else if (errTypes[b] == 4) {
      double scl, pw;
      rpemGuardedPow(counts, endptp, b, scl, pw);
      sdAdd[b] = scl; sdProp[b] = pw;   // scale in add slot, exponent in prop slot
    } else if (errTypes[b] == 3) {
      double add, lam;
      rpemGuardedTBS(counts, endptp, b, add, lam);
      sdAdd[b] = add; sdProp[b] = lam;  // add.sd in add slot, lambda in prop slot
    } else {
      sdAdd[b] = sqrt(sumSS[b] / (double)sumN[b]); sdProp[b] = NA_REAL;
    }
  }
  return List::create(_["coefs"] = coefOut, _["omega"] = omegaNew,
                      _["sd"] = sdAdd, _["propSd"] = sdProp, _["accept"] = accept);
}
