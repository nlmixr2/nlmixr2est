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
#include <vector>

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
  std::vector<double> logp;  // [nsub*nGauss] log p(Y_i | theta_ij)
  std::vector<double> etaS;  // [nsub*nGauss*nEta] drawn etas
};

rpemOptions rpemOp;

//[[Rcpp::export]]
RObject rpemFree() {
  rpemOp.buf.clear();
  rpemOp.nobs = rpemOp.idS = rpemOp.idF = NULL;
  rpemOp.nobsTot = 0;
  rpemOp.nsub = 0;
  rpemOp.loaded = false;
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
static inline double rpemSolveSubject(const double *par, int id) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) {
    setIndParPtr(ind, (int)i, par[i]);
  }
  iniSubjectE(id, 1, ind, op, rx, rxPred.update_inis);
  rpemPredOde(id);
  double s = 0.0;
  double curT;
  int kk;
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
      double v = lhs[0];
      if (ISNA(v)) v = 0.0;
      s += v;
    }
  }
  return s;
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

// K=1 E-step with in-C++ eta sampling via rxode2's thread-safe threefry RNG
// (design/rpem/06, D18). All etas are drawn UP FRONT with the registered
// _rxode2_rxRmvnSEXP_ so no RNG runs inside the (future) parallel solve loop.
//   base:   length-ntheta template param row (ETA slots overwritten per draw)
//   etaIdx: 0-based positions of the ETA params in the row
//   omega:  nEta x nEta between-subject covariance (not Cholesky)
// Draws nGauss etas per subject ~ N(0, omega), solves, accumulates n_i via
// log-sum-exp, stores per-sample log p and etas for the M-step, and returns
// logn, lnL, and the drawn etas.
//[[Rcpp::export]]
List rpemEstepK1Draw(Environment e, NumericVector base, IntegerVector etaIdx,
                     NumericMatrix omega, int nGauss, int ncores) {
  if (_rxode2_rxRmvnSEXP_ == NULL) stop("rxode2 rxRmvn pointer not initialized");
  int nEta = etaIdx.size();
  if (omega.nrow() != nEta || omega.ncol() != nEta) stop("omega must be nEta x nEta");

  RObject pred = e["predOnly"];
  List rxControl = as<List>(e["rxControl"]);
  NumericVector param = as<NumericVector>(e["param"]);
  RObject data = e["data"];

  // Establish the solve once to learn nsub / ntheta.
  rpemDoSetup(pred, rxControl, param, data);
  int nsub = rpemOp.nsub;
  if ((unsigned int)base.size() != rpemOp.ntheta) stop("base must have ntheta entries");
  int nAll = nsub * nGauss;

  // Threefry multivariate-normal draws (thread-safe), untruncated (+/-Inf bounds).
  // NOTE: this resets rxode2's global solve state, so we re-establish `rx` below.
  NumericVector mu(nEta);                       // zeros
  NumericVector lower(nEta, R_NegInf), upper(nEta, R_PosInf);
  SEXP draw = _rxode2_rxRmvnSEXP_(wrap(IntegerVector::create(nAll)),
                                  wrap(mu), wrap(omega), wrap(lower), wrap(upper),
                                  wrap(IntegerVector::create(ncores)),
                                  wrap(LogicalVector::create(false)),  // isChol
                                  wrap(LogicalVector::create(false)),  // keepNames
                                  wrap(NumericVector::create(0.4)),
                                  wrap(NumericVector::create(2.05)),
                                  wrap(NumericVector::create(1e-10)),
                                  wrap(IntegerVector::create(100)));
  NumericMatrix etaAll(draw);                   // nAll x nEta

  // Re-establish the E-step solve struct clobbered by the rxRmvn draw.
  rpemDoSetup(pred, rxControl, param, data);

  rpemOp.nGauss = nGauss;
  rpemOp.nEta = nEta;
  rpemOp.logp.assign((size_t)nAll, 0.0);
  rpemOp.etaS.assign((size_t)nAll * nEta, 0.0);

  // Copy all R objects into plain C++ buffers BEFORE the parallel region; the
  // OpenMP loop must not touch any Rcpp/R object (only these buffers + the
  // rxode2 C solve). etaS doubles as the shared eta buffer for the loop.
  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) rpemOp.etaS[k] = etaAll[k];
  std::vector<double> baseBuf(rpemOp.ntheta);
  for (unsigned int i = 0; i < rpemOp.ntheta; ++i) baseBuf[i] = base[i];
  std::vector<int> etaIdxBuf(nEta);
  for (int a = 0; a < nEta; ++a) etaIdxBuf[a] = etaIdx[a];
  std::vector<double> lognV(nsub, 0.0);

  // Serial per-subject solve.  (OpenMP parallelization is a follow-up: it must go
  // through rxode2's par_solve, which allocates the per-thread ODE workspace --
  // driving ind_solve from a manual OpenMP loop segfaults because that scratch is
  // not set up per thread.  See design/rpem/04, 06.)
  std::vector<double> row(rpemOp.ntheta), lp(nGauss);
  for (int id = 0; id < nsub; ++id) {
    double mx = R_NegInf;
    for (int j = 0; j < nGauss; ++j) {
      size_t r = (size_t)id * nGauss + j;
      for (unsigned int i = 0; i < rpemOp.ntheta; ++i) row[i] = baseBuf[i];
      for (int a = 0; a < nEta; ++a) row[etaIdxBuf[a]] = rpemOp.etaS[r * nEta + a];
      double logp = -rpemSolveSubject(row.data(), id);
      lp[j] = logp;
      rpemOp.logp[r] = logp;
      if (logp > mx) mx = logp;
    }
    double s = 0.0;
    for (int j = 0; j < nGauss; ++j) s += exp(lp[j] - mx);
    lognV[id] = mx + log(s) - log((double)nGauss);
  }

  // Build the R return objects serially, after the parallel region.
  NumericVector logn(nsub);
  NumericMatrix etaOut(nAll, nEta);
  double lnL = 0.0;
  for (int id = 0; id < nsub; ++id) { logn[id] = lognV[id]; lnL += lognV[id]; }
  for (size_t k = 0; k < (size_t)nAll * nEta; ++k) etaOut[k] = rpemOp.etaS[k];
  return List::create(_["logn"] = logn, _["lnL"] = lnL, _["eta"] = etaOut);
}
