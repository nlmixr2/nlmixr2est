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
#include "nmMcmcRng.h" // nmSetSeedEng1 / nmRestoreMcmcSeed: undo a solve's per-subject
                       // threefry re-seed so it never carries into the sampling draws
#include <vector>
#include <functional>

#define _(String) (String)

#include "utilc.h"     // RSprintf, used by scale.h's iteration-print routines
#include "scale.h"     // shared iteration-print + parameter-history machinery (after the
                       // gettext-style _() macro so its warning(_(...)) calls resolve)

// L-BFGS-B with C linkage (src/lbfgsR.c), shared with saem.cpp -- used to refine the
// fixed-effect likelihood parameters of a general log-likelihood RPEM model by a direct
// box-constrained optimization of the importance-weighted observation log-likelihood.
typedef double rpemOptimfn(int n, double *par, void *ex);
typedef void rpemOptimgr(int n, double *par, double *gr, void *ex);
extern "C" void lbfgsbRX(int n, int lmm, double *x, double *lower,
                         double *upper, int *nbd, double *Fmin, rpemOptimfn fn,
                         rpemOptimgr gr, int *fail, void *ex, double factr,
                         double pgtol, int *fncount, int *grcount,
                         int maxit, char *msg, int trace, int nREPORT);

// Box-constrained L-BFGS-B refinement of the ll() structural (likelihood) params over the
// current E-step samples (defined after the RpemLikOpt objective; used by the cLoop).
void rpemLbfgsBeta(const std::vector<double> &baseBuf, const std::vector<int> &structIdxBuf,
                   const std::vector<int> &etaIdxBuf, std::vector<double> &betaIO,
                   const double *lower, const double *upper, const int *nbd,
                   int lmm, double factr, double pgtol, int maxit);

// Per-endpoint residual optimizers (multi-endpoint / errType 5); defined later, used by the
// cLoop's multi-endpoint M-step to re-optimize each endpoint's combined/power/TBS residual
// over just that endpoint's observations.
static void rpemGuardedComb(const std::vector<long> &counts, const int *endpt,
                            int endB, double a0, double b0, double &aOut, double &bOut);
static void rpemGuardedPow(const std::vector<long> &counts, const int *endpt,
                           int endB, double &sclOut, double &powOut);
static void rpemGuardedTBS(const std::vector<long> &counts, const int *endpt,
                           int endB, double &addOut, double &lamOut);

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

// Shared scale.h iteration-print state (setup in R via rpemIterPrintStart_); declared here
// so the C++ E-M loops (rpemEMLoopK1 / rpemEMLoopMix) can stream rows live as they iterate.
static scaling _rpemScale;

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
  // Record the sampling-block engine seed (the post-eta-draw state) so the per-subject
  // re-seeding inside the E-step par_solve -- which under multiple threads leaves the engine
  // in a thread-scheduling-dependent state -- can be undone before the R-side M-step's rxRmvn
  // draws, keeping the M-step deterministic for any core count (nmMcmcRng guard).
  nmSetSeedEng1(getRxSeed1(ncores > 1 ? ncores : 1));
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
  // Undo the solve's per-subject re-seed so the M-step's rxRmvn is deterministic.
  nmRestoreMcmcSeed();

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

  // Record the sampling-block seed so the E-step solves' per-subject re-seed can be undone
  // before the R-side mixture M-step's rxRmvn draws (nmMcmcRng guard; deterministic for any
  // core count).
  nmSetSeedEng1(getRxSeed1(ncores > 1 ? ncores : 1));
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
  // Undo the solve's per-subject re-seed so the mixture M-step's rxRmvn is deterministic.
  nmRestoreMcmcSeed();

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
List rpemMstepK1(NumericVector muIn, double addSd0, int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1");
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

  // Pre-draw all MH uniforms: 3 per trial (proposal i', proposal j', accept).  Drawn
  // directly from rxode2's threefry engine (rxUnifEng), seeded per iteration via
  // nmSetSeedEng1 so the MH stream is reproducible AND independent of the E-step solve's
  // per-subject re-seed (same idiom as saem's _saemFillUnifEng).
  int total = nTrials + burn;
  size_t nU = (size_t)3 * total;
  nmSetSeedEng1(seed);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = rxUnifEng(0.0, 1.0);

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
// with a deterministic per-(iteration, subject) seed on the EVEN threefry stream --
// setSeedEng1(seed0 + (it*nsub + i)*2) after setRxThreadId() -- so it is thread-safe,
// reproducible for any core count, and niter-independent (extending niter shares the
// exact prefix of a shorter run at the same seed; the MH uses the odd stream).
// The ODE solve keeps par_solve (rxode2 owns the parallel solver; a manual OpenMP loop
// over the solve is unsafe).  errType 0 additive, 1 proportional.  Returns per-iteration
// mu / omega (diagonal) / add.sd / lnL traces.
//[[Rcpp::export]]
List rpemEMLoopK1(Environment e, List cfg) {
  // The single-component C++ E-M loop takes one config List (built in R by .rpemFit) rather
  // than a long positional argument list, so new estimation features can add config keys
  // without a signature/registration/caller churn.
  NumericVector base = cfg["base"];
  IntegerVector etaIdx = cfg["etaIdx"], muIdx = cfg["muIdx"];
  int addSdIdx = cfg["addSdIdx"], errType = cfg["errType"];
  NumericVector mu0 = cfg["mu0"], omDiag0 = cfg["omDiag0"];
  double addSd0 = cfg["addSd0"];
  IntegerVector resIdx = cfg["resIdx"];
  NumericVector resPar0 = cfg["resPar0"];
  IntegerVector structIdx = cfg["structIdx"];
  NumericVector struct0 = cfg["struct0"];
  int niter = cfg["niter"], nGauss = cfg["nGauss"], ncores = cfg["ncores"];
  int nMH = cfg["nMH"], mhBurn = cfg["mhBurn"];
  unsigned int seed = (unsigned int)(int)cfg["seed"];
  NumericMatrix design = cfg["design"];
  IntegerVector covCoefIdx = cfg["covCoefIdx"];
  NumericVector structLower = cfg["structLower"], structUpper = cfg["structUpper"];
  IntegerVector structNbd = cfg["structNbd"];
  int likLbfgs = cfg["likLbfgs"], collect = cfg["collect"], lbfgsLmm = cfg["lbfgsLmm"];
  double lbfgsFactr = cfg["lbfgsFactr"], lbfgsPgtol = cfg["lbfgsPgtol"];
  int lbfgsMaxIter = cfg["lbfgsMaxIter"];
  double cInflate = cfg["cInflate"];
  // live iteration printing (design like saem/vae): when printLive, emit each iteration's
  // theta + omega row through the shared scale.h machinery (set up in R via
  // rpemIterPrintStart_) as the loop runs, so the trace appears incrementally rather than
  // all at once after the loop.  nThetaPrint is the number of leading theta columns.
  int printLive = cfg.containsElementNamed("printLive") ? (int)cfg["printLive"] : 0;
  int nThetaPrint = cfg.containsElementNamed("nThetaPrint") ? (int)cfg["nThetaPrint"] : 0;
  // Post-M-step holds/clamps (mirror the R loop) so the cLoop covers every K=1 case:
  //   muRefBuf[a]==0 -> centered eta (typical value pinned at 0)
  //   etaFixBuf[a]==1 -> omega diagonal a held at its initial value omDiag0
  //   omGroupBuf[a]   -> IOV pooling group; omega is averaged within a group
  //   muFixBuf[a]==1  -> typical value a held at its initial value mu0
  //   addSdFix/propSdFix/lambdaFix/powFix -> that residual parameter held at its initial value
  IntegerVector muRefV  = cfg.containsElementNamed("muRef")  ? as<IntegerVector>(cfg["muRef"])  : IntegerVector(0);
  IntegerVector etaFixV = cfg.containsElementNamed("etaFix") ? as<IntegerVector>(cfg["etaFix"]) : IntegerVector(0);
  IntegerVector omGroupV= cfg.containsElementNamed("omGroup")? as<IntegerVector>(cfg["omGroup"]): IntegerVector(0);
  IntegerVector muFixV  = cfg.containsElementNamed("muFix")  ? as<IntegerVector>(cfg["muFix"])  : IntegerVector(0);
  int addSdFix  = cfg.containsElementNamed("addSdFix")  ? (int)cfg["addSdFix"]  : 0;
  int propSdFix = cfg.containsElementNamed("propSdFix") ? (int)cfg["propSdFix"] : 0;
  int lambdaFix = cfg.containsElementNamed("lambdaFix") ? (int)cfg["lambdaFix"] : 0;
  int powFix    = cfg.containsElementNamed("powFix")    ? (int)cfg["powFix"]    : 0;
  // multi-endpoint (errType 5): per-obs endpoint index + per-endpoint residual thetas /
  // types (empty for single-endpoint models).
  IntegerVector endpt = cfg["endpt"], endErrType = cfg["endErrType"];
  IntegerVector endSclIdx = cfg["endSclIdx"], endPropIdx = cfg["endPropIdx"];
  NumericVector endScl0 = cfg["endScl0"], endProp0 = cfg["endProp0"];
  bool doMulti = (errType == 5);
  int nEndpt = endErrType.size();
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
  bool doLik = (errType == 7);   // general log-likelihood endpoint: no residual parameter
  bool doLnorm = (errType == 6); // lognormal: residual SD estimated on the log scale
  int propSdIdx = resIdx[0], powIdx = resIdx[1], lambdaIdx = resIdx[2];
  double propSd = resPar0[0], power = resPar0[1], lambda = resPar0[2];
  int nStruct = structIdx.size();          // non-mu-ref structural fixed effects (beta)
  std::vector<int> structIdxBuf(nStruct);
  std::vector<double> beta(nStruct);
  for (int k = 0; k < nStruct; ++k) { structIdxBuf[k] = structIdx[k]; beta[k] = struct0[k]; }
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

  // mu2 covariate regression (D22): when covariate coefficients are present the mu
  // M-step becomes a weighted linear regression of the accepted theta samples on the
  // per-subject design (nEta==1).  coefs = [typical, covCoef...]; the covariate coefs
  // are written to their theta slots each iteration and theta_ij = design_i.coefs + eta.
  bool doReg = (covCoefIdx.size() > 0);
  int nCoef = doReg ? design.ncol() : 1;
  std::vector<int> covIdxBuf(doReg ? covCoefIdx.size() : 0);
  std::vector<double> coefs(doReg ? nCoef : 0), dmat(doReg ? (size_t)nsub * nCoef : 0);
  if (doReg) {
    if (nEta != 1) stop("rpemEMLoopK1 covariate regression requires nEta==1");
    if (design.nrow() != nsub) stop("design must have one row per subject");
    if (nCoef != (int)covCoefIdx.size() + 1) stop("design columns must be 1 + covCoefIdx");
    for (int k = 0; k < (int)covCoefIdx.size(); ++k) covIdxBuf[k] = covCoefIdx[k];
    coefs[0] = mu[0];
    for (int k = 0; k < (int)covCoefIdx.size(); ++k) coefs[k + 1] = baseBuf[covIdxBuf[k]];
    for (int i = 0; i < nsub; ++i)
      for (int k = 0; k < nCoef; ++k) dmat[(size_t)i * nCoef + k] = design(i, k);
  }

  NumericMatrix muTr(niter, nEta), omTr(niter, nEta), betaTr(niter, nStruct);
  NumericMatrix coefTr(niter, doReg ? nCoef : 0);
  NumericVector sdTr(niter), propTr(niter, NA_REAL), powTr(niter, NA_REAL);
  NumericVector lamTr(niter, NA_REAL), llTr(niter);
  // per-endpoint residual scale/second-param (multi-endpoint) + their per-iteration traces.
  std::vector<double> sdVec(doMulti ? nEndpt : 0), propVec(doMulti ? nEndpt : 0);
  std::vector<int> endSclBuf(doMulti ? nEndpt : 0), endPropBuf(doMulti ? nEndpt : 0);
  std::vector<int> endErrBuf(doMulti ? nEndpt : 0);
  NumericMatrix sdMat(doMulti ? niter : 0, doMulti ? nEndpt : 0);
  NumericMatrix propMat(doMulti ? niter : 0, doMulti ? nEndpt : 0);
  if (doMulti)
    for (int b = 0; b < nEndpt; ++b) {
      sdVec[b] = endScl0[b]; propVec[b] = R_IsNA(endProp0[b]) ? 0.0 : endProp0[b];
      endSclBuf[b] = endSclIdx[b]; endPropBuf[b] = endPropIdx[b]; endErrBuf[b] = endErrType[b];
    }
  std::vector<double> row(rpemOp.ntheta);
  // counts3: accepted-sample visit counts, needed by the residual re-optimizers
  // (combined/power/TBS and, when the data have BLQ records, the censored additive/
  // proportional M-step).  Censoring (doCens) is auto-detected from the E-step's censv.
  std::vector<long> counts3((size_t)nsub * nGauss, 0);
  bool doCens = false, censChecked = false;
  // mode-centered importance sampling (impInflate / cInflate > 0): draw eta ~ N(EBE_i,
  // cInflate*Omega) and importance-weight against the prior N(0,Omega); the per-subject EBE
  // (posterior mean) is updated each iteration and kept entirely in C++.  cInflate == 1 with
  // EBE == 0 (the default) reduces to prior sampling with logRatio == 0 (byte-identical).
  bool doModeIS = (cInflate != 1.0);
  std::vector<double> ebeBuf((size_t)nsub * nEta, 0.0);
  double halfNetaLogC = 0.5 * nEta * log(cInflate);
  bool doPar = (ncores > 1);
  seedEng(ncores);
  uint32_t seed0 = (uint32_t)seed;

  for (int it = 0; it < niter; ++it) {
    if (doReg) {
      mu[0] = coefs[0];
      for (int k = 0; k < (int)covIdxBuf.size(); ++k) baseBuf[covIdxBuf[k]] = coefs[k + 1];
    }
    // centered (non-mu-referenced) etas have muIdx -1: no theta carries their typical value
    for (int a = 0; a < nEta; ++a) if (muIdxBuf[a] >= 0) baseBuf[muIdxBuf[a]] = mu[a];
    if (addSdIdx >= 0) baseBuf[addSdIdx] = addSd;   // LL (doLik): no residual parameter
    if (doComb) baseBuf[propSdIdx] = propSd;
    if (doPow) baseBuf[powIdx] = power;
    if (doTbs) baseBuf[lambdaIdx] = lambda;
    if (doMulti)
      for (int b = 0; b < nEndpt; ++b) {
        baseBuf[endSclBuf[b]] = sdVec[b];
        if ((endErrBuf[b] == 2 || endErrBuf[b] == 3 || endErrBuf[b] == 4) && endPropBuf[b] >= 0)
          baseBuf[endPropBuf[b]] = propVec[b];
      }
    for (int k = 0; k < nStruct; ++k) baseBuf[structIdxBuf[k]] = beta[k];
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
      nmSetSeedEng1(seed0 + (uint32_t)(((size_t)it * nsub + id) * 2));
      for (int j = 0; j < nGauss; ++j)
        for (int a = 0; a < nEta; ++a)
          rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a] =
            ebeBuf[(size_t)id * nEta + a] + rxNormEng(0.0, 1.0) * sqrt(cInflate * omDiag[a]);
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
        // importance weight of the mode-centered proposal vs the prior (0 when !doModeIS).
        double lr = 0.0;
        if (doModeIS) {
          lr = halfNetaLogC;
          for (int a = 0; a < nEta; ++a) {
            double eta = rpemOp.etaS[r * nEta + a], om = omDiag[a];
            double d = eta - ebeBuf[(size_t)id * nEta + a];
            lr += -0.5 * eta * eta / om + 0.5 * d * d / (cInflate * om);
          }
        }
        rpemOp.logp[r] = logLik + lr; rpemOp.ssv[r] = ssTmp; rpemOp.wssv[r] = wssTmp;
      }
    }
    // The par_solve above re-seeds the threefry engine per subject; restore the sampling
    // seed so that leak never carries into the structural / MH draws below (nmMcmcRng guard).
    nmRestoreMcmcSeed();
    // Detect BLQ censoring once from the populated censv (any non-zero cens code): the
    // naive pooled SS is biased under censoring, so the additive/proportional residual is
    // then estimated by maximizing the censored log-likelihood (as rpemMstepK1Cens does).
    if (!censChecked) {
      for (unsigned int k = 0; k < rpemOp.nobsTot; ++k) if (rpemOp.censv[k] != 0) { doCens = true; break; }
      censChecked = true;
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
    // mode-centered IS: refresh each subject's EBE (posterior mean) = sum_j w_ij eta_ij with
    // w_ij = softmax_j(logw), the next iteration's proposal center (kept in C++).
    if (doModeIS) {
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
          ebeBuf[(size_t)id * nEta + a] = acc;
        }
      }
    }

    // Numeric M-step for non-mu-referenced structural fixed effects (the paper's beta):
    // one damped diagonal-Newton step maximizing the importance-weighted complete-data
    // Q(beta) = sum_ij w_ij log p(Y_i | beta, eta_ij) -- re-solving per (i,j) with the
    // held etas.  Runs while the E-step solve is loaded, before the MH clobbers it.
    // ll() likelihood-parameter refinement: during the terminal smoothing window, refine
    // the structural (likelihood) params with the box-constrained L-BFGS-B (respecting
    // their declared bounds), matching the R loop; exploration iterations use the cheaper
    // damped-Newton below.
    if (nStruct > 0 && doLik && likLbfgs && it >= niter - collect) {
      bool hasBnd = ((int)structNbd.size() == nStruct);
      rpemLbfgsBeta(baseBuf, structIdxBuf, etaIdxBuf, beta,
                    hasBnd ? structLower.begin() : nullptr, hasBnd ? structUpper.begin() : nullptr,
                    hasBnd ? structNbd.begin() : nullptr, lbfgsLmm, lbfgsFactr, lbfgsPgtol, lbfgsMaxIter);
    } else if (nStruct > 0) {
      std::vector<double> wij((size_t)nsub * nGauss);
      for (int i = 0; i < nsub; ++i) {
        double mx = R_NegInf;
        for (int j = 0; j < nGauss; ++j) { double v = rpemOp.logp[(size_t)i * nGauss + j]; if (v > mx) mx = v; }
        double sw = 0.0;
        for (int j = 0; j < nGauss; ++j) { double ee = exp(rpemOp.logp[(size_t)i * nGauss + j] - mx); wij[(size_t)i * nGauss + j] = ee; sw += ee; }
        for (int j = 0; j < nGauss; ++j) wij[(size_t)i * nGauss + j] /= sw;
      }
      auto Qbeta = [&](const std::vector<double> &b) -> double {
        double q = 0.0;
        for (int i = 0; i < nsub; ++i)
          for (int j = 0; j < nGauss; ++j) {
            for (unsigned int tt = 0; tt < rpemOp.ntheta; ++tt) row[tt] = baseBuf[tt];
            for (int k = 0; k < nStruct; ++k) row[structIdxBuf[k]] = b[k];
            for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)i * nGauss + j) * nEta + a];
            q += wij[(size_t)i * nGauss + j] * (-rpemSolveSubject(row.data(), i));
          }
        return q;
      };
      const double hb = 1e-3;
      double Q0 = Qbeta(beta);
      std::vector<double> step(nStruct, 0.0);
      for (int k = 0; k < nStruct; ++k) {
        std::vector<double> bp = beta, bm = beta; bp[k] += hb; bm[k] -= hb;
        double qp = Qbeta(bp), qm = Qbeta(bm);
        double g = (qp - qm) / (2.0 * hb), hess = (qp - 2.0 * Q0 + qm) / (hb * hb);
        step[k] = (hess < -1e-8) ? -g / hess : g;
      }
      double tstep = 1.0;
      for (int ls = 0; ls < 30; ++ls) {
        std::vector<double> bn = beta;
        for (int k = 0; k < nStruct; ++k) bn[k] = beta[k] + tstep * step[k];
        if (Qbeta(bn) > Q0) { beta = bn; break; }
        tstep *= 0.5;
      }
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
    // MH uses the ODD threefry stream (eta draw uses the EVEN stream), so the two never
    // collide AND neither seed depends on the total niter -- extending niter reproduces the
    // exact prefix of a shorter run at the same seed (dynamic-iteration stable, imp.cpp style).
    nmSetSeedEng1(seed0 + (uint32_t)it * 2u + 1u);
    std::vector<double> U(nU);
    for (size_t k = 0; k < nU; ++k) U[k] = R::pnorm(rxNormEng(0.0, 1.0), 0.0, 1.0, 1, 0);
    int ci = 0, cj = 0; double clogp = rpemOp.logp[0];
    std::vector<double> sumT(nEta, 0.0), sumTT((size_t)nEta * nEta, 0.0);
    // regression accumulators (doReg): per-subject linear predictor + normal equations
    std::vector<double> muLin(doReg ? nsub : 0, 0.0);
    std::vector<double> XtX(doReg ? (size_t)nCoef * nCoef : 0, 0.0), Xty(doReg ? nCoef : 0, 0.0);
    double Stt = 0.0;
    if (doReg)
      for (int i = 0; i < nsub; ++i) {
        double s = 0.0;
        for (int k = 0; k < nCoef; ++k) s += dmat[(size_t)i * nCoef + k] * coefs[k];
        muLin[i] = s;
      }
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
        if (doReg) {
          double theta = muLin[ci] + rpemOp.etaS[off];
          const double *xi = &dmat[(size_t)ci * nCoef];
          for (int a = 0; a < nCoef; ++a) {
            Xty[a] += xi[a] * theta;
            for (int b = 0; b < nCoef; ++b) XtX[(size_t)a * nCoef + b] += xi[a] * xi[b];
          }
          Stt += theta * theta;
        } else {
          for (int a = 0; a < nEta; ++a) {
            double tha = mu[a] + rpemOp.etaS[off + a];
            sumT[a] += tha;
            for (int b = 0; b < nEta; ++b) sumTT[(size_t)a * nEta + b] += tha * (mu[b] + rpemOp.etaS[off + b]);
          }
        }
        sumSS += (errType == 1) ? rpemOp.wssv[(size_t)ci * nGauss + cj] : rpemOp.ssv[(size_t)ci * nGauss + cj];
        if (doComb || doPow || doTbs || doCens || doMulti || doLnorm) counts3[(size_t)ci * nGauss + cj]++;
        sumNobs += rpemOp.nobs[ci]; ++m;
      }
    }
    if (doReg) {
      arma::mat A(XtX.data(), nCoef, nCoef);   // symmetric -> row/col-major agree
      arma::vec bvec(Xty.data(), nCoef);
      arma::vec betaNew = arma::solve(A, bvec, arma::solve_opts::likely_sympd);
      for (int k = 0; k < nCoef; ++k) coefs[k] = betaNew[k];
      mu[0] = coefs[0];
      omDiag[0] = (Stt - arma::dot(betaNew, bvec)) / (double)m;
    } else {
      for (int a = 0; a < nEta; ++a) mu[a] = sumT[a] / (double)m;
      for (int a = 0; a < nEta; ++a) omDiag[a] = sumTT[(size_t)a * nEta + a] / (double)m - mu[a] * mu[a];
    }
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
    } else if (doCens) {
      // censored additive (errType 0) / proportional (errType 1) residual: maximize the
      // censored log-likelihood over the accepted (subject, sample) states -- observed
      // records contribute the Gaussian density, BLQ records the CENS probability
      // (doCensNormal1) -- by golden-section over sd (as rpemMstepK1Cens).
      auto Qc = [&](double sd) -> double {
        double q = 0.0;
        for (int i = 0; i < nsub; ++i) {
          int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
          for (int j = 0; j < nGauss; ++j) {
            long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
            long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
            for (int o = 0; o < nobsi; ++o) {
              long so = so0 + o;
              double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
              double sdo = (errType == 1) ? sd * (fabs(cp) + 1e-300) : sd;
              double rv = sdo * sdo;
              double gauss = -M_LN_SQRT_2PI - log(sdo) - 0.5 * (dv - cp) * (dv - cp) / rv;
              q += (double)c * doCensNormal1((double)rpemOp.censv[so], dv, rpemOp.limv[so],
                                             gauss, cp, rv, 0);
            }
          }
        }
        return q;
      };
      addSd = rpemGolden(Qc, 1e-4, std::max(5.0 * addSd, 1.0), 120);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doMulti) {
      // multi-endpoint: additive/proportional/lognormal endpoints get a closed-form SS over
      // their own observations; combined/power/TBS endpoints re-optimize over theirs (same
      // per-endpoint logic as rpemMstepK1Multi).
      std::vector<double> ssE(nEndpt, 0.0); std::vector<long> nE(nEndpt, 0);
      const int *ep = &endpt[0];
      for (int i = 0; i < nsub; ++i) {
        int nobsi = rpemOp.nobs[i]; long b0 = rpemOp.idS[i];
        for (int j = 0; j < nGauss; ++j) {
          long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
          long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
          for (int o = 0; o < nobsi; ++o) {
            int b = ep[b0 + o];
            if (endErrBuf[b] == 2 || endErrBuf[b] == 3 || endErrBuf[b] == 4) continue;  // numeric below
            double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o], rr = dv - cp, contrib;
            if (endErrBuf[b] == 6) { double lr = (cp > 0.0 && dv > 0.0) ? (log(dv) - log(cp)) : 0.0; contrib = lr * lr; }
            else if (endErrBuf[b] == 1) { contrib = (cp != 0.0) ? (rr / cp) * (rr / cp) : 0.0; }
            else { contrib = rr * rr; }
            ssE[b] += (double)c * contrib; nE[b] += c;
          }
        }
      }
      for (int b = 0; b < nEndpt; ++b) {
        if (endErrBuf[b] == 2) { double aa, bp; rpemGuardedComb(counts3, ep, b, sdVec[b]*sdVec[b], propVec[b]*propVec[b], aa, bp); sdVec[b] = sqrt(aa); propVec[b] = sqrt(bp); }
        else if (endErrBuf[b] == 4) { double scl, pw; rpemGuardedPow(counts3, ep, b, scl, pw); sdVec[b] = scl; propVec[b] = pw; }
        else if (endErrBuf[b] == 3) { double add, lam; rpemGuardedTBS(counts3, ep, b, add, lam); sdVec[b] = add; propVec[b] = lam; }
        else { sdVec[b] = (nE[b] > 0) ? sqrt(ssE[b] / (double)nE[b]) : sdVec[b]; }
      }
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doLnorm) {
      // lognormal: residual SD on the LOG scale = sqrt(mean (log dv - log cp)^2) over the
      // accepted (subject, sample) states (the pooled natural-scale SS is wrong here).
      double ss = 0.0; long n = 0;
      for (int i = 0; i < nsub; ++i) {
        int nobsi = rpemOp.nobs[i];
        for (int j = 0; j < nGauss; ++j) {
          long c = counts3[(size_t)i * nGauss + j]; if (!c) continue;
          long ob = rpemOp.sampObsOff[i] + (long)j * nobsi;
          for (int o = 0; o < nobsi; ++o) {
            double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
            double lr = (cp > 0.0 && dv > 0.0) ? (log(dv) - log(cp)) : 0.0;
            ss += (double)c * lr * lr; n += c;
          }
        }
      }
      if (n > 0) addSd = sqrt(ss / (double)n);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (!doLik) {                 // LL: no residual sd to update
      addSd = sqrt(sumSS / (double)sumNobs);
    }
    // Post-M-step holds/clamps (mirror the R loop), applied before recording traces so the
    // parameter history reflects them.  Centered etas / held omega / IOV pooling apply to the
    // scalar (non-regression) mu update; fixed typical + residual params apply to both.
    if (!doReg) {
      for (int a = 0; a < nEta; ++a)
        if (a < muRefV.size() && muRefV[a] == 0) mu[a] = 0.0;              // centered eta
      for (int a = 0; a < nEta; ++a)
        if (a < etaFixV.size() && etaFixV[a] != 0) omDiag[a] = omDiag0[a]; // held omega
      if ((int)omGroupV.size() == nEta) {                                  // IOV omega pooling
        std::vector<int> gid; std::vector<double> gsum, gn;
        for (int a = 0; a < nEta; ++a) {
          int g = omGroupV[a], gi = -1;
          for (int q = 0; q < (int)gid.size(); ++q) if (gid[q] == g) { gi = q; break; }
          if (gi < 0) { gid.push_back(g); gsum.push_back(0.0); gn.push_back(0.0); gi = (int)gid.size() - 1; }
          gsum[gi] += omDiag[a]; gn[gi] += 1.0;
        }
        for (int a = 0; a < nEta; ++a)
          for (int q = 0; q < (int)gid.size(); ++q) if (gid[q] == omGroupV[a]) { omDiag[a] = gsum[q] / gn[q]; break; }
      }
    }
    if ((int)muFixV.size() == nEta) {                                      // held typical values
      if (doReg) { if (muFixV[0] != 0) { coefs[0] = mu0[0]; mu[0] = coefs[0]; } }
      else for (int a = 0; a < nEta; ++a) if (muFixV[a] != 0) mu[a] = mu0[a];
    }
    if (addSdFix && addSdIdx >= 0) addSd = addSd0;                         // held residual params
    if (propSdFix && doComb) propSd = resPar0[0];
    if (powFix && doPow) power = resPar0[1];
    if (lambdaFix && doTbs) lambda = resPar0[2];
    for (int a = 0; a < nEta; ++a) { muTr(it, a) = mu[a]; omTr(it, a) = omDiag[a]; }
    for (int k = 0; k < nStruct; ++k) betaTr(it, k) = beta[k];
    if (doReg) for (int k = 0; k < nCoef; ++k) coefTr(it, k) = coefs[k];
    sdTr[it] = addSd; llTr[it] = lnL;
    if (doComb) propTr[it] = propSd;
    if (doPow) powTr[it] = power;
    if (doTbs) lamTr[it] = lambda;
    if (doMulti) for (int b = 0; b < nEndpt; ++b) {
      sdMat(it, b) = sdVec[b];
      bool hasProp = (endErrBuf[b] == 2 || endErrBuf[b] == 3 || endErrBuf[b] == 4);
      propMat(it, b) = hasProp ? propVec[b] : NA_REAL;   // NA for add/prop/lnorm endpoints
    }
    // Live iteration row: assemble the post-M-step theta vector (baseBuf's leading theta
    // slots re-filled with this iteration's estimates) + the omega diagonal, and stream it
    // through the shared scale.h printer.  phase = terminal averaging window (Smooth) vs
    // exploration (EM), matching the R post-loop reconstruction.
    if (printLive && nThetaPrint > 0) {
      std::vector<double> prow((size_t)nThetaPrint + nEta);
      for (int i = 0; i < nThetaPrint; ++i) prow[i] = baseBuf[i];
      if (doReg) {
        prow[muIdxBuf[0]] = coefs[0];
        for (int k = 0; k < (int)covIdxBuf.size(); ++k) prow[covIdxBuf[k]] = coefs[k + 1];
      } else {
        for (int a = 0; a < nEta; ++a) if (muIdxBuf[a] >= 0) prow[muIdxBuf[a]] = mu[a];
      }
      if (addSdIdx >= 0) prow[addSdIdx] = addSd;
      if (doComb) prow[propSdIdx] = propSd;
      if (doPow) prow[powIdx] = power;
      if (doTbs) prow[lambdaIdx] = lambda;
      if (doMulti)
        for (int b = 0; b < nEndpt; ++b) {
          prow[endSclBuf[b]] = sdVec[b];
          if ((endErrBuf[b] == 2 || endErrBuf[b] == 3 || endErrBuf[b] == 4) && endPropBuf[b] >= 0)
            prow[endPropBuf[b]] = propVec[b];
        }
      for (int k = 0; k < nStruct; ++k) prow[structIdxBuf[k]] = beta[k];
      for (int a = 0; a < nEta; ++a) prow[(size_t)nThetaPrint + a] = omDiag[a];
      _rpemScale.phaseLabel = ((it + 1) > (niter - collect)) ? "Smooth" : "EM";
      scalePrintFun(&_rpemScale, prow.data(), llTr[it]);
    }
  }
  rpemFree();
  NumericMatrix ebeOut(nsub, nEta);   // converged proposal center (mode-centered IS)
  for (int id = 0; id < nsub; ++id)
    for (int a = 0; a < nEta; ++a) ebeOut(id, a) = ebeBuf[(size_t)id * nEta + a];
  return List::create(_["muTrace"] = muTr, _["omegaTrace"] = omTr,
                      _["sdTrace"] = sdTr, _["propTrace"] = propTr,
                      _["powTrace"] = powTr, _["lamTrace"] = lamTr,
                      _["betaTrace"] = betaTr, _["coefTrace"] = coefTr,
                      _["ebe"] = ebeOut, _["sdMat"] = sdMat, _["propMat"] = propMat,
                      _["lnL"] = llTr);
}

// Full C++ E-M loop for mixtures (design/rpem/12 M5): the whole mixture E-M in one C++
// call -- threefry eta draw (per-eta omega, split-ETA aware), per-component par_solve
// E-step (setIndMixest), and the joint (subject, sample, component) MH M-step with
// threefry uniforms -- so it is thread-safe and reproducible for any core count.
// Additive / proportional / lognormal residual (combined / power / TBS mixtures keep the
// R loop).  Returns per-iteration traces of the per-(parameter, component) typical
// values, the mixture weights, the per-eta omega, and add.sd.
//[[Rcpp::export]]
List rpemEMLoopMix(Environment e, NumericVector base, IntegerVector etaIdx,
                   NumericMatrix muComp0, IntegerMatrix muCompIdx, IntegerMatrix etaForComp,
                   NumericVector omDiag0, int addSdIdx, int errType, double addSd0,
                   IntegerVector resIdx, NumericVector resPar0,
                   NumericVector w0, int K, int nParam, bool perComp,
                   int niter, int nGauss, int ncores, int nMH, int mhBurn, unsigned int seed) {
  RObject pred = e["predOnly"]; List rxControl = as<List>(e["rxControl"]);
  NumericVector param = as<NumericVector>(e["param"]); RObject data = e["data"];
  int nEta = etaIdx.size();
  if (rxNormEng == NULL || seedEng == NULL || setSeedEng1 == NULL)
    stop("rxode2 threefry engine not initialized");
  // shared residual: additive/prop/lnorm (closed form) vs combined/power/TBS (re-optimize
  // over the accepted (subject, sample, component) states with the mixture stride).
  bool doComb = (errType == 2), doPow = (errType == 4), doTbs = (errType == 3);
  int propSdIdx = resIdx[0], powIdx = resIdx[1], lambdaIdx = resIdx[2];
  double propSd = resPar0[0], power = resPar0[1], lambda = resPar0[2];

  rpemDoSetup(pred, rxControl, param, data);
  int nsub = rpemOp.nsub, nAll = nsub * nGauss;
  rpemOp.nGauss = nGauss; rpemOp.nEta = nEta; rpemOp.nMix = K;
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

  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) baseBuf[i] = base[i];
  std::vector<int> etaIdxBuf(nEta);
  for (int a = 0; a < nEta; ++a) etaIdxBuf[a] = etaIdx[a];
  std::vector<double> muComp((size_t)nParam * K); std::vector<int> muCIdx((size_t)nParam * K), efc((size_t)nParam * K);
  for (int p = 0; p < nParam; ++p) for (int k = 0; k < K; ++k) {
    muComp[(size_t)p * K + k] = muComp0(p, k); muCIdx[(size_t)p * K + k] = muCompIdx(p, k); efc[(size_t)p * K + k] = etaForComp(p, k);
  }
  std::vector<double> omDiag(nEta); for (int a = 0; a < nEta; ++a) omDiag[a] = omDiag0[a];
  std::vector<double> w(K); for (int k = 0; k < K; ++k) w[k] = w0[k];
  double addSd = addSd0;

  NumericVector muTr((size_t)niter * nParam * K), wTr((size_t)niter * K), omTr((size_t)niter * nEta);
  NumericVector sdTr(niter), propTr(niter, NA_REAL), powTr(niter, NA_REAL), lamTr(niter, NA_REAL), llTr(niter);
  std::vector<double> row(rpemOp.ntheta);
  std::vector<long> counts3((doComb || doPow || doTbs) ? (size_t)nsub * nGauss * K : 0, 0);
  bool doPar = (ncores > 1);
  seedEng(ncores);
  uint32_t seed0 = (uint32_t)seed;

  for (int it = 0; it < niter; ++it) {
    for (int p = 0; p < nParam; ++p) for (int k = 0; k < K; ++k) baseBuf[muCIdx[(size_t)p * K + k]] = muComp[(size_t)p * K + k];
    baseBuf[addSdIdx] = addSd;
    if (doComb) baseBuf[propSdIdx] = propSd;
    if (doPow) baseBuf[powIdx] = power;
    if (doTbs) baseBuf[lambdaIdx] = lambda;
    rpemDoSetup(pred, rxControl, param, data);

    // threefry eta draw (per-eta omega); deterministic per-(iter, subject).
#ifdef _OPENMP
#pragma omp parallel for num_threads(ncores) if(doPar)
#endif
    for (int id = 0; id < nsub; ++id) {
#ifdef _OPENMP
      if (doPar) setRxThreadId(omp_get_thread_num());
#endif
      nmSetSeedEng1(seed0 + (uint32_t)(((size_t)it * nsub + id) * 2));
      for (int j = 0; j < nGauss; ++j) for (int a = 0; a < nEta; ++a)
        rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a] = rxNormEng(0.0, 1.0) * sqrt(omDiag[a]);
    }

    // E-step: solve each sample once per component (setIndMixest), read -> logp[(i,j,k)].
    for (int kc = 0; kc < K; ++kc) {
      for (int j = 0; j < nGauss; ++j) {
        for (int id = 0; id < nsub; ++id) {
          for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
          for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[((size_t)id * nGauss + j) * nEta + a];
          rpemSetSubject(row.data(), id);
          setIndMixest(getSolvingOptionsInd(rx, id), kc + 1);
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
      }
    }
    // undo the par_solve's per-subject threefry re-seed before the MH draws (nmMcmcRng).
    nmRestoreMcmcSeed();

    // per-subject mixture log n_i.
    std::vector<double> logw(K), logn(nsub); double lnL = 0.0;
    for (int k = 0; k < K; ++k) logw[k] = log(w[k]);
    for (int id = 0; id < nsub; ++id) {
      double mxk = R_NegInf; std::vector<double> lnik(K);
      for (int k = 0; k < K; ++k) {
        double mx = R_NegInf;
        for (int j = 0; j < nGauss; ++j) { double v = rpemOp.logp[((size_t)id * nGauss + j) * K + k]; if (v > mx) mx = v; }
        double s = 0.0;
        for (int j = 0; j < nGauss; ++j) s += exp(rpemOp.logp[((size_t)id * nGauss + j) * K + k] - mx);
        lnik[k] = mx + log(s) - log((double)nGauss);
        double wk = logw[k] + lnik[k]; if (wk > mxk) mxk = wk;
      }
      double sk = 0.0; for (int k = 0; k < K; ++k) sk += exp(logw[k] + lnik[k] - mxk);
      logn[id] = mxk + log(sk); lnL += logn[id];
    }

    // M-step: joint (subject, sample, component) MH with threefry uniforms.
    int total = nMH + mhBurn; size_t nU = (size_t)4 * total;
#ifdef _OPENMP
    setRxThreadId(0);
#endif
    // MH uses the ODD threefry stream (eta draw uses the EVEN stream), so the two never
    // collide AND neither seed depends on the total niter -- extending niter reproduces the
    // exact prefix of a shorter run at the same seed (dynamic-iteration stable, imp.cpp style).
    nmSetSeedEng1(seed0 + (uint32_t)it * 2u + 1u);
    std::vector<double> U(nU);
    for (size_t kk = 0; kk < nU; ++kk) U[kk] = R::pnorm(rxNormEng(0.0, 1.0), 0.0, 1.0, 1, 0);
    int ci = 0, cj = 0, ck = 0; double clogp = rpemOp.logp[0];
    std::vector<double> sumTk((size_t)nParam * K, 0.0), sumTTk((size_t)nParam * K, 0.0);
    std::vector<long> countK(K, 0);
    double sumSS = 0.0; long sumNobs = 0, m = 0;
    for (int t = 0; t < total; ++t) {
      double u1 = U[(size_t)4 * t], u2 = U[(size_t)4 * t + 1], u3 = U[(size_t)4 * t + 2], u4 = U[(size_t)4 * t + 3];
      int pih = (int)(u1 * nsub); if (pih >= nsub) pih = nsub - 1;
      int pjh = (int)(u2 * nGauss); if (pjh >= nGauss) pjh = nGauss - 1;
      int pkh = (int)(u3 * K); if (pkh >= K) pkh = K - 1;
      double plogp = rpemOp.logp[((size_t)pih * nGauss + pjh) * K + pkh];
      double logA = (logw[pkh] + plogp - logn[pih]) - (logw[ck] + clogp - logn[ci]);
      if (log(u4) < logA) { ci = pih; cj = pjh; ck = pkh; clogp = plogp; }
      if (t >= mhBurn) {
        size_t etaBase = ((size_t)ci * nGauss + cj) * nEta;
        for (int p = 0; p < nParam; ++p) {
          double eta = rpemOp.etaS[etaBase + efc[(size_t)p * K + ck]];
          double theta = muComp[(size_t)p * K + ck] + eta;
          sumTk[(size_t)p * K + ck] += theta; sumTTk[(size_t)p * K + ck] += theta * theta;
        }
        ++countK[ck];
        int nobsi = rpemOp.nobs[ci];
        size_t r = ((size_t)ci * nGauss + cj) * K + ck;
        if (errType == 6) {
          long ob = (rpemOp.sampObsOff[ci] + (long)cj * nobsi) * K + (long)ck * nobsi;
          for (int o = 0; o < nobsi; ++o) { double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o]; if (cp > 0.0 && dv > 0.0) { double lr = log(dv) - log(cp); sumSS += lr * lr; } }
        } else sumSS += (errType == 1) ? rpemOp.wssv[r] : rpemOp.ssv[r];
        if (doComb || doPow || doTbs) counts3[((size_t)ci * nGauss + cj) * K + ck]++;
        sumNobs += nobsi; ++m;
      }
    }
    std::vector<double> sumSSk(nEta, 0.0); std::vector<long> countEta(nEta, 0);
    for (int k = 0; k < K; ++k) {
      for (int p = 0; p < nParam; ++p) {
        if (countK[k] > 0) {
          muComp[(size_t)p * K + k] = sumTk[(size_t)p * K + k] / (double)countK[k];
          int a = efc[(size_t)p * K + k];
          sumSSk[a] += sumTTk[(size_t)p * K + k] - (double)countK[k] * muComp[(size_t)p * K + k] * muComp[(size_t)p * K + k];
          countEta[a] += countK[k];
        }
      }
      w[k] = (double)countK[k] / (double)m;
    }
    for (int a = 0; a < nEta; ++a) omDiag[a] = (countEta[a] > 0) ? sumSSk[a] / (double)countEta[a] : 0.0;
    double wFloor = 1e-3, wsum = 0.0;
    for (int k = 0; k < K; ++k) { if (w[k] < wFloor) w[k] = wFloor; wsum += w[k]; }
    for (int k = 0; k < K; ++k) w[k] /= wsum;
    // residual: add/prop/lnorm closed form; combined/power/TBS re-optimize over the
    // accepted (subject, sample, component) states (mixture stride).
    if (doComb) {
      const double eps = 1e-12;
      auto over = [&](double aa, double bb, double &Q, double *ga, double *gb,
                      double *Haa, double *Hab, double *Hbb) {
        Q = 0.0; if (ga) { *ga = *gb = *Haa = *Hab = *Hbb = 0.0; }
        for (int i = 0; i < nsub; ++i) { int nobsi = rpemOp.nobs[i];
          for (int j = 0; j < nGauss; ++j) for (int k = 0; k < K; ++k) {
            long c = counts3[((size_t)i * nGauss + j) * K + k]; if (!c) continue;
            long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
            for (int o = 0; o < nobsi; ++o) {
              double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp, wcp = cp * cp, r2 = rr * rr, V = aa + bb * wcp; if (V < eps) V = eps;
              double invV = 1.0 / V; Q += (double)c * (-0.5 * log(V) - 0.5 * r2 * invV);
              if (ga) { double g = -0.5 * invV + 0.5 * r2 * invV * invV, h = 0.5 * invV * invV - r2 * invV * invV * invV;
                *ga += c * g; *gb += c * wcp * g; *Haa += c * h; *Hab += c * wcp * h; *Hbb += c * wcp * wcp * h; }
            } } }
      };
      double a = addSd * addSd, bpar = propSd * propSd, Qcur; over(a, bpar, Qcur, 0,0,0,0,0);
      for (int itn = 0; itn < 200; ++itn) {
        double Qtmp, ga, gb, Haa, Hab, Hbb; over(a, bpar, Qtmp, &ga, &gb, &Haa, &Hab, &Hbb);
        if (fabs(ga) + fabs(gb) < 1e-9) break;
        double det = Haa * Hbb - Hab * Hab, da, db;
        if (Haa < 0.0 && det > 0.0) { da = -(Hbb * ga - Hab * gb) / det; db = -(-Hab * ga + Haa * gb) / det; }
        else { double sc = 1.0 / (fabs(Haa) + fabs(Hbb) + 1.0); da = sc * ga; db = sc * gb; }
        double step = 1.0; bool ok = false;
        for (int ls = 0; ls < 40; ++ls) { double an = a + step * da, bn = bpar + step * db;
          if (an < eps) an = eps; if (bn < eps) bn = eps; double Qn; over(an, bn, Qn, 0,0,0,0,0);
          if (Qn > Qcur) { a = an; bpar = bn; Qcur = Qn; ok = true; break; } step *= 0.5; }
        if (!ok) break;
      }
      addSd = sqrt(a); propSd = sqrt(bpar);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doPow) {
      auto stat = [&](double cc, double &SSc, double &SL) { SSc = 0.0; SL = 0.0;
        for (int i = 0; i < nsub; ++i) { int nobsi = rpemOp.nobs[i];
          for (int j = 0; j < nGauss; ++j) for (int k = 0; k < K; ++k) {
            long c = counts3[((size_t)i * nGauss + j) * K + k]; if (!c) continue;
            long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
            for (int o = 0; o < nobsi; ++o) { double cp = rpemOp.cpv[ob + o], rr = rpemOp.dvv[ob + o] - cp, acp = fabs(cp) + 1e-300;
              SSc += (double)c * rr * rr / pow(acp, 2.0 * cc); SL += (double)c * log(acp); } } } };
      double N = (double)sumNobs;
      auto f = [&](double cc) { double SSc, SL; stat(cc, SSc, SL); if (SSc < 1e-300) SSc = 1e-300; return -0.5 * (N * log(SSc / N) + 2.0 * cc * SL); };
      power = rpemGolden(f, 0.0, 3.0, 100);
      double SSc, SL; stat(power, SSc, SL); addSd = sqrt(SSc / N);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else if (doTbs) {
      auto ssJac = [&](double lam, double &SS, double &Jac) { SS = 0.0; Jac = 0.0;
        for (int i = 0; i < nsub; ++i) { int nobsi = rpemOp.nobs[i]; long so0 = rpemOp.idS[i];
          for (int j = 0; j < nGauss; ++j) for (int k = 0; k < K; ++k) {
            long c = counts3[((size_t)i * nGauss + j) * K + k]; if (!c) continue;
            long ob = (rpemOp.sampObsOff[i] + (long)j * nobsi) * K + (long)k * nobsi;
            for (int o = 0; o < nobsi; ++o) { double cp = rpemOp.cpv[ob + o], dv = rpemOp.dvv[ob + o];
              int yj = rpemOp.yjv[so0 + o]; double low = rpemOp.lowv[so0 + o], hi = rpemOp.hiv[so0 + o];
              double d = _powerD(dv, lam, yj, low, hi) - _powerD(cp, lam, yj, low, hi);
              SS += (double)c * d * d; Jac += (double)c * log(fabs(_powerDD(dv, lam, yj, low, hi)) + 1e-300); } } } };
      double N = (double)sumNobs;
      auto f = [&](double lam) { double SS, Jac; ssJac(lam, SS, Jac); if (SS < 1e-300) SS = 1e-300; return -0.5 * N * log(SS / N) + Jac; };
      lambda = rpemGolden(f, -2.0, 3.0, 100);
      double SS, Jac; ssJac(lambda, SS, Jac); addSd = sqrt(SS / N);
      std::fill(counts3.begin(), counts3.end(), 0);
    } else {
      addSd = sqrt(sumSS / (double)sumNobs);
    }
    // label-switching guard: shared-eta components are exchangeable -> order by the
    // first parameter's mu (ascending) each iteration; split-ETA components are tied to
    // distinct symbols, so they are left in place.
    if (!perComp && K > 1) {
      std::vector<int> ord(K); for (int k = 0; k < K; ++k) ord[k] = k;
      std::sort(ord.begin(), ord.end(), [&](int x, int y) { return muComp[x] < muComp[y]; });
      std::vector<double> mc2((size_t)nParam * K), w2(K);
      for (int k = 0; k < K; ++k) { w2[k] = w[ord[k]]; for (int p = 0; p < nParam; ++p) mc2[(size_t)p * K + k] = muComp[(size_t)p * K + ord[k]]; }
      muComp = mc2; w = w2;
    }
    for (int p = 0; p < nParam; ++p) for (int k = 0; k < K; ++k) muTr[(size_t)it * nParam * K + (size_t)p * K + k] = muComp[(size_t)p * K + k];
    for (int k = 0; k < K; ++k) wTr[(size_t)it * K + k] = w[k];
    for (int a = 0; a < nEta; ++a) omTr[(size_t)it * nEta + a] = omDiag[a];
    sdTr[it] = addSd; llTr[it] = lnL;
    if (doComb) propTr[it] = propSd;
    if (doPow) powTr[it] = power;
    if (doTbs) lamTr[it] = lambda;
  }
  rpemFree();
  return List::create(_["muTrace"] = muTr, _["wTrace"] = wTr, _["omegaTrace"] = omTr,
                      _["sdTrace"] = sdTr, _["propTrace"] = propTr, _["powTrace"] = powTr,
                      _["lamTrace"] = lamTr, _["lnL"] = llTr);
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
                  int errType, double addSd0, double propSd0, int nTrials, int burn,
                  unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepMixDraw before rpemMstepMix");
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
  // MH uniforms drawn directly from the threefry engine (rxUnifEng), seeded per iteration
  // via nmSetSeedEng1 -- reproducible and independent of the E-step solve's re-seed.
  nmSetSeedEng1(seed);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = rxUnifEng(0.0, 1.0);

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
                    int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Reg");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Reg currently supports nEta==1");
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
  // MH uniforms drawn directly from the threefry engine (rxUnifEng), seeded per iteration
  // via nmSetSeedEng1 -- reproducible and independent of the E-step solve's re-seed.
  nmSetSeedEng1(seed);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = rxUnifEng(0.0, 1.0);

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
                     double propSd0, int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Comb");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Comb currently supports nEta==1");
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
  // MH uniforms drawn directly from the threefry engine (rxUnifEng), seeded per iteration
  // via nmSetSeedEng1 -- reproducible and independent of the E-step solve's re-seed.
  nmSetSeedEng1(seed);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = rxUnifEng(0.0, 1.0);

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
                      unsigned int seed, std::vector<long> &counts, NumericVector &coefOut,
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
  // MH uniforms drawn directly from the threefry engine (rxUnifEng), seeded per iteration
  // via nmSetSeedEng1 -- reproducible and independent of the E-step solve's re-seed.
  nmSetSeedEng1(seed);
  std::vector<double> U(nU);
  for (size_t k = 0; k < nU; ++k) U[k] = rxUnifEng(0.0, 1.0);
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
                    int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1TBS");
  if (rpemOp.nEta != 1) stop("rpemMstepK1TBS currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  if (design.nrow() != nsub) stop("design must have one row per subject");
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, seed, counts, coefOut, omegaNew, accept, m);

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
                    double power0, int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Pow");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Pow currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  if (design.nrow() != nsub) stop("design must have one row per subject");
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, seed, counts, coefOut, omegaNew, accept, m);

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
                     double sd0, int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Cens");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Cens currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, seed, counts, coefOut, omegaNew, accept, m);

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

// L-BFGS-B refinement of the fixed-effect likelihood parameters of a general log-
// likelihood RPEM model (mirrors saem's refinePhi0Lik / saemix ind.fix10).  Instead of
// rpemMstepBeta's single damped-Newton step, this box-constrained-optimizes the same
// importance-weighted complete-data objective Q(beta) = sum_ij w_ij log p(Y_i|beta,eta_ij)
// (etas held from the E-step) via L-BFGS-B, so bounded parameters are respected natively
// (the box IS the clamping -- lbfgsbRX never leaves [lower, upper]).  A stochastic-
// approximation smoothing gain damps the update: beta <- beta0 + gain*(x_opt - beta0),
// re-clamped to the box.  Must be called while the E-step solve struct is loaded.
struct RpemLikOpt {
  std::vector<double> w;                 // importance weights (nsub*nG)
  std::vector<double> baseBuf;           // full theta row template
  const int *structIdx; int nB;
  const int *etaIdx; int nEta;
  int nsub, nG;
  std::vector<double> row;               // scratch (ntheta)
  long nfe;
};

static double rpemLikNegQ(int n, double *par, void *ex) {
  RpemLikOpt *o = (RpemLikOpt *)ex;
  double q = 0.0;
  for (int i = 0; i < o->nsub; ++i)
    for (int j = 0; j < o->nG; ++j) {
      for (unsigned int t = 0; t < rpemOp.ntheta; ++t) o->row[t] = o->baseBuf[t];
      for (int k = 0; k < n; ++k) o->row[o->structIdx[k]] = par[k];
      for (int a = 0; a < o->nEta; ++a)
        o->row[o->etaIdx[a]] = rpemOp.etaS[((size_t)i * o->nG + j) * o->nEta + a];
      q += o->w[(size_t)i * o->nG + j] * (-rpemSolveSubject(o->row.data(), i));
    }
  o->nfe++;
  if (!std::isfinite(q)) return 1e300;   // out-of-support candidate -> reject
  return -q;                             // minimize -Q
}

static void rpemLikGrad(int n, double *par, double *gr, void *ex) {
  for (int k = 0; k < n; ++k) {
    double p0 = par[k], h = 1e-4 * (std::fabs(p0) + 1e-4);
    par[k] = p0 + h; double fp = rpemLikNegQ(n, par, ex);
    par[k] = p0 - h; double fm = rpemLikNegQ(n, par, ex);
    par[k] = p0;
    gr[k] = (fp - fm) / (2.0 * h);
  }
}

// Shared core: box-constrained L-BFGS-B optimum of Q(beta) over the structural params,
// using the importance weights from the CURRENT rpemOp E-step state.  betaIO seeds the
// search and receives the box-constrained optimum (clamped to the box).  Used by both the
// R-exported rpemMstepBetaLik and the C++ cLoop (rpemEMLoopK1) so ll() bounded-parameter
// refinement is identical whether the loop runs in R or C++.
void rpemLbfgsBeta(const std::vector<double> &baseBuf,
                   const std::vector<int> &structIdxBuf, const std::vector<int> &etaIdxBuf,
                   std::vector<double> &betaIO, const double *lower, const double *upper,
                   const int *nbd, int lmm, double factr, double pgtol, int maxit) {
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss, nB = (int)structIdxBuf.size();
  RpemLikOpt o;
  o.w.assign((size_t)nsub * nG, 0.0);
  for (int i = 0; i < nsub; ++i) {
    double mx = R_NegInf;
    for (int j = 0; j < nG; ++j) { double v = rpemOp.logp[(size_t)i * nG + j]; if (v > mx) mx = v; }
    double s = 0.0;
    for (int j = 0; j < nG; ++j) { double e = exp(rpemOp.logp[(size_t)i * nG + j] - mx); o.w[(size_t)i * nG + j] = e; s += e; }
    for (int j = 0; j < nG; ++j) o.w[(size_t)i * nG + j] /= s;
  }
  o.baseBuf = baseBuf;
  o.row.assign(rpemOp.ntheta, 0.0);
  o.structIdx = structIdxBuf.data(); o.nB = nB; o.etaIdx = etaIdxBuf.data();
  o.nEta = (int)etaIdxBuf.size(); o.nsub = nsub; o.nG = nG; o.nfe = 0;
  std::vector<double> x(betaIO), lo(nB, 0.0), up(nB, 0.0); std::vector<int> bnd(nB, 0);
  for (int k = 0; k < nB; ++k) if (nbd) { bnd[k] = nbd[k]; lo[k] = lower[k]; up[k] = upper[k]; }
  double Fmin; int fail = 0, fncount = 0, grcount = 0; char msg[100];
  lbfgsbRX(nB, lmm, x.data(), lo.data(), up.data(), bnd.data(), &Fmin,
           rpemLikNegQ, rpemLikGrad, &fail, (void *)&o, factr, pgtol,
           &fncount, &grcount, maxit, msg, 0, maxit + 1);
  for (int k = 0; k < nB; ++k) {
    double v = x[k];
    if (nbd) { if ((bnd[k] == 1 || bnd[k] == 2) && v < lo[k]) v = lo[k];
               if ((bnd[k] == 2 || bnd[k] == 3) && v > up[k]) v = up[k]; }
    betaIO[k] = v;
  }
}

//[[Rcpp::export]]
NumericVector rpemMstepBetaLik(NumericVector base, IntegerVector etaIdx,
                               IntegerVector structIdx, NumericVector struct0,
                               NumericVector lower, NumericVector upper, IntegerVector nbd,
                               double gain, int lmm, double factr, double pgtol, int maxit) {
  if (!rpemOp.loaded) stop("run rpemEstepK1Draw before rpemMstepBetaLik");
  int nB = structIdx.size();
  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int t = 0; t < rpemOp.ntheta; ++t) baseBuf[t] = base[t];
  std::vector<int> sIdx(structIdx.begin(), structIdx.end()), eIdx(etaIdx.begin(), etaIdx.end());
  std::vector<double> beta(struct0.begin(), struct0.end());
  bool hasBnd = ((int)nbd.size() == nB);
  rpemLbfgsBeta(baseBuf, sIdx, eIdx, beta,
                hasBnd ? lower.begin() : nullptr, hasBnd ? upper.begin() : nullptr,
                hasBnd ? nbd.begin() : nullptr, lmm, factr, pgtol, maxit);
  // smoothing-gain damped update, re-clamped to the box.
  NumericVector out(nB);
  for (int k = 0; k < nB; ++k) {
    double v = struct0[k] + gain * (beta[k] - struct0[k]);
    if (hasBnd) {
      if ((nbd[k] == 1 || nbd[k] == 2) && v < lower[k]) v = lower[k];
      if ((nbd[k] == 2 || nbd[k] == 3) && v > upper[k]) v = upper[k];
    }
    out[k] = v;
  }
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
                      int nTrials, int burn, unsigned int seed) {
  if (rpemOp.nGauss == 0) stop("run rpemEstepK1Draw before rpemMstepK1Multi");
  if (rpemOp.nEta != 1) stop("rpemMstepK1Multi currently supports nEta==1");
  int nsub = rpemOp.nsub, nG = rpemOp.nGauss;
  int nEndpt = errTypes.size();
  std::vector<long> counts; NumericVector coefOut; double omegaNew, accept; long m;
  rpemMHReg(design, coefs, nTrials, burn, seed, counts, coefOut, omegaNew, accept, m);

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

// RPEM iteration printing + parameter-history capture.  Reuses the SHARED iteration-print
// machinery in scale.h (scaleSetup / scalePrintHeader / scalePrintFun / scaleParHisDf) that
// saem, focei, the nlm family and vae use, so RPEM prints the same iteration table and
// produces parHistData in the standard format.  RPEM never scales its parameters
// (scaleTypeNone drops the redundant "U" rows); the R-side `xform` list drives the "X"
// back-transform row (e.g. exp() typical values), and `phase` labels the algorithm stage
// (EM exploration vs the terminal smoothing/averaging window).
// (_rpemScale is declared near the top of the file so the C++ E-M loops can stream rows.)
static std::vector<double> _rpemIpInit, _rpemIpScaleC;
static std::vector<int> _rpemIpXPar;
static std::string _rpemIpPhase;

//[[Rcpp::export]]
RObject rpemIterPrintStart_(NumericVector initPar, CharacterVector names,
                            List iterPrintControl, RObject xform = R_NilValue) {
  int np = initPar.size();
  _rpemIpInit.assign(initPar.begin(), initPar.end());
  _rpemIpScaleC.assign(np, 1.0);
  scaleSetup(&_rpemScale, _rpemIpInit.data(), _rpemIpScaleC.data(), names,
             /*useColor*/ 0, /*printNcol*/ np, /*print*/ 1,
             normTypeConstant, scaleTypeNone, 0.0, 0.0, 0.0, np);
  if (!Rf_isNull(xform)) scaleAttachXform(&_rpemScale, as<List>(xform));
  if (_rpemScale.xPar == NULL) {
    _rpemIpXPar.assign(np, 0);
    _rpemScale.xPar = _rpemIpXPar.data();
    _rpemScale.probitIdx = NULL;
    _rpemScale.logitThetaLow = _rpemScale.logitThetaHi = NULL;
    _rpemScale.probitThetaLow = _rpemScale.probitThetaHi = NULL;
  }
  _rpemScale.keyExtra = "EM: exploration; Smooth: terminal averaging (collect) window\n";
  scaleApplyIterPrintControl(&_rpemScale, iterPrintControl);
  if (_rpemScale.every > 0) scalePrintHeader(&_rpemScale);   // header only when displaying
  return R_NilValue;
}

//[[Rcpp::export]]
RObject rpemIterPrintRow_(NumericVector x, double f, std::string phase = "") {
  _rpemIpPhase = phase;
  _rpemScale.phaseLabel = _rpemIpPhase.empty() ? NULL : _rpemIpPhase.c_str();
  scalePrintFun(&_rpemScale, &x[0], f);
  return R_NilValue;
}

//[[Rcpp::export]]
RObject rpemIterPrintGet_(bool printLine = true) {
  _rpemScale.save = 0;
  _rpemScale.every = 0;
  if (printLine) scalePrintLine(&_rpemScale, min2(_rpemScale.npars, _rpemScale.ncol));
  return scaleParHisDf(&_rpemScale);
}
