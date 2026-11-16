#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "inner.h"
#include "odeSwap.h"
#include "nmParallel.h"
#include "rxomp.h"
#include <R_ext/Rdynload.h>
#include <vector>

// Post-fit table solves through the shared ODE-swap pool (#1356).
//
// The table step solves the fit's inner model (ipred/pred, with the eta
// sensitivities CWRES needs) and its predOnly model (the display columns) over
// the fit's own translated events.  Instead of a fresh rxSolve() per model, R
// builds ONE pool with rxSolve(.setupOnly = 1L) for the larger model and both
// are registered as slots; each solve then writes only the columns asked for.
//
// The read loop follows rxode2_df's: every record calls calc_lhs (a dose first
// advances tlast, so tad/dosenum match), and only observation records (evid 0
// and 2) are output.

static rxSolveF rxTableInner;
static rxSolveF rxTablePred;
static long _tablePoolSolveN = 0;

typedef void (*tablePoolTlast_t)(double *time, rx_solving_options_ind *ind);
static tablePoolTlast_t tablePoolTlast = NULL;

static rxSolveF *tablePoolFns(int slot) {
  if (slot == odeSlotInner) return &rxTableInner;
  if (slot == odeSlotTable) return &rxTablePred;
  return NULL;
}

long tablePoolSolveN() { return _tablePoolSolveN; }

static bool tablePoolCanSolve(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op);

//[[Rcpp::export]]
RObject tablePoolFree_() {
  odeSwapClearAll();
  rxClearFuns(&rxTableInner);
  rxClearFuns(&rxTablePred);
  return R_NilValue;
}

// Declare the table models (metadata only) and return the model the pool must
// be built with: the one with the most states.
//[[Rcpp::export]]
RObject tablePoolDeclare_(RObject inner, RObject predOnly) {
  tablePoolFree_();
  if (!Rf_isNull(inner)) odeSwapDeclare(odeSlotInner, "tableInner", inner);
  if (!Rf_isNull(predOnly)) odeSwapDeclare(odeSlotTable, "tablePred", predOnly);
  return RObject(odeSwapPoolModelSEXP());
}

// Bind the declared models once the pool exists.  Returns, per slot, whether
// it can be solved and read in this pool.
//[[Rcpp::export]]
LogicalVector tablePoolRegister_(RObject inner, RObject predOnly) {
  if (tablePoolTlast == NULL) {
    tablePoolTlast = (tablePoolTlast_t) R_GetCCallable("rxode2", "handleTlast");
  }
  rx_solve *rx = getRxSolve_();
  LogicalVector ret = LogicalVector::create(_["inner"] = false, _["predOnly"] = false);
  if (rx == NULL) return ret;
  rx_solving_options *op = getSolvingOptions(rx);
  if (!Rf_isNull(inner) && odeSwapRegister(odeSlotInner, "tableInner", inner, &rxTableInner)) {
    ret[0] = tablePoolCanSolve(odeSlotInner, &rxTableInner, rx, op);
  }
  if (!Rf_isNull(predOnly) && odeSwapRegister(odeSlotTable, "tablePred", predOnly, &rxTablePred)) {
    ret[1] = tablePoolCanSolve(odeSlotTable, &rxTablePred, rx, op);
  }
  return ret;
}

// odeSwapCanPool() also refuses a model wider than rxode2's lhs slice; here
// every read goes through OdeSwapScope::lhs(), which hands back the private
// buffer odeSwapPlan() sized for it.
static bool tablePoolCanSolve(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op) {
  int deny = odeSwapCanPool(slot);
  if (deny == odeDenyPoolNotSized && getOpNeq(op) >= odeSwapNeq(slot) &&
      odeSwapPlan().scratchNlhs >= odeSwapNlhs(slot)) {
    deny = odeDenyNone;
  }
  return deny == odeDenyNone && odeSwapCheckLhsWidth(slot, fns, rx, op, true);
}

static inline bool tablePoolOutRow(int evid) { return evid == 0 || evid == 2; }

// Solve every subject as table model `which` (0 inner, 1 predOnly) and read
// the requested columns.
//
// parMat:    one row per subject (or one row for all of them); its columns are
//            written to the 0-based par_ptr indices parIdx before solving.
// lhsIdx:    0-based lhs indices of the slot model to output as columns.
// sensIdx:   0-based lhs indices returned as sparse (row, k, value) triplets.
// stateIdx:  0-based state indices to output.
// covIdx:    0-based par_ptr indices to output (covariates, as interpolated).
// tbs:       also output rxLambda/rxYj/rxLow/rxHi.
//[[Rcpp::export]]
List tablePoolSolve_(int which, NumericMatrix parMat, IntegerVector parIdx,
                     IntegerVector lhsIdx, IntegerVector sensIdx,
                     IntegerVector stateIdx, IntegerVector covIdx, bool tbs) {
  int slot = (which == 0) ? (int)odeSlotInner : (int)odeSlotTable;
  rxSolveF *fns = tablePoolFns(slot);
  rx_solve *rx = getRxSolve_();
  if (fns == NULL || rx == NULL || !odeSwapLoaded(slot) || fns->calc_lhs == NULL) {
    stop("table pool slot %d is not set up", slot);
  }
  rx_solving_options *op = getSolvingOptions(rx);
  int nsub = (int)getRxNsub(rx);
  if (parMat.nrow() != 1 && parMat.nrow() != nsub) {
    stop("table pool parameters have %d rows for %d subjects", parMat.nrow(), nsub);
  }
  if (parMat.ncol() != parIdx.size()) stop("table pool parameter names do not match");
  // output rows per subject, and each subject's first output row.  Counted by
  // record, not through ix: the time order is only sorted by the solve.
  std::vector<int> start((size_t)nsub + 1, 0);
  for (int i = 0; i < nsub; ++i) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, i);
    int n = 0, nAll = getIndNallTimes(ind);
    for (int j = 0; j < nAll; ++j) {
      if (tablePoolOutRow(getIndEvid(ind, j))) n++;
    }
    start[(size_t)i + 1] = start[(size_t)i] + n;
  }
  int nr = start[(size_t)nsub];
  int nl = lhsIdx.size(), ns = stateIdx.size(), nc = covIdx.size(),
    nsens = sensIdx.size(), ntbs = tbs ? 4 : 0;
  int neqSlot = odeSwapNeq(slot);
  for (int k = 0; k < ns; ++k) {
    if (stateIdx[k] < 0 || stateIdx[k] >= neqSlot) stop("table pool state index out of range");
  }
  int nlhsSlot = odeSwapNlhs(slot);
  for (int k = 0; k < nl; ++k) {
    if (lhsIdx[k] < 0 || lhsIdx[k] >= nlhsSlot) stop("table pool lhs index out of range");
  }
  for (int k = 0; k < nsens; ++k) {
    if (sensIdx[k] < 0 || sensIdx[k] >= nlhsSlot) stop("table pool lhs index out of range");
  }
  IntegerVector id(nr), rec(nr), evid(nr);
  NumericVector time(nr);
  int ncol = nl + ns + nc + ntbs;
  NumericMatrix vals(nr, ncol);
  double *pv = REAL(vals);
  std::vector<std::vector<int> > sensR((size_t)nsub), sensK((size_t)nsub);
  std::vector<std::vector<double> > sensV((size_t)nsub);
  int *pid = INTEGER(id), *prec = INTEGER(rec), *pevid = INTEGER(evid);
  double *ptime = REAL(time);
  const double *ppar = REAL(parMat);
  int npr = parMat.nrow(), npc = parMat.ncol();
  const int *pparIdx = INTEGER(parIdx);
  std::vector<int> rowBad((size_t)nsub, 0);
  int cores = getOpCores(op);
  bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  // the batch installs (or turns off) this slot's event-sensitivity shape
  OdeSwapEsBatch esBatch(slot);
  nmForEachSubject(rx, nsub, cores, doParallel, [&](int i) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, i);
    int pr = (npr == 1) ? 0 : i;
    for (int k = 0; k < npc; ++k) setIndParPtr(ind, pparIdx[k], ppar[pr + k * npr]);
    OdeSwapScope neqGuard(slot, ind, op);
    OdeSwapCmtScope cmtGuard(slot, op, ind);
    setIndSolve(ind, -1);
    odeSwapSolveInd(slot, i);
    bool bad = odeSwapIndBadSolveSlot(op, ind, slot);
    iniSubjectE(i, 1, ind, op, rx, fns->update_inis);
    double *lhs = neqGuard.lhs();
    int nAll = getIndNallTimes(ind);
    int r = start[(size_t)i], rEnd = start[(size_t)i + 1];
    for (int j = 0; j < nAll; ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, j);
      int ev = getIndEvid(ind, kk);
      if (ev == 9) continue;
      double curT = getTime(kk, ind);
      if (!tablePoolOutRow(ev)) {
        if (ev >= 100 || ev == 3) tablePoolTlast(&curT, ind);
        if (!bad) {
          lhsWalkNewind(ind, j);
          fns->calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
        }
        continue;
      }
      if (r >= rEnd) { rowBad[(size_t)i] = 1; break; }
      pid[r] = i + 1;
      prec[r] = kk;
      pevid[r] = ev;
      ptime[r] = curT;
      if (bad) {
        for (int k = 0; k < ncol; ++k) pv[r + (size_t)k * nr] = NA_REAL;
        // a failed solve's sensitivities are missing, not structural zeros
        for (int k = 0; k < nsens; ++k) {
          sensR[(size_t)i].push_back(r);
          sensK[(size_t)i].push_back(k);
          sensV[(size_t)i].push_back(NA_REAL);
        }
        r++;
        continue;
      }
      lhsWalkNewind(ind, j);
      double *solve = getOpIndSolve(op, ind, j);
      fns->calc_lhs(i, curT, solve, lhs);
      int c = 0;
      for (int k = 0; k < nl; ++k) pv[r + (size_t)(c++) * nr] = lhs[lhsIdx[k]];
      for (int k = 0; k < ns; ++k) pv[r + (size_t)(c++) * nr] = solve[stateIdx[k]];
      for (int k = 0; k < nc; ++k) pv[r + (size_t)(c++) * nr] = getIndParPtr(ind, covIdx[k]);
      if (tbs) {
        pv[r + (size_t)(c++) * nr] = getIndLambda(ind);
        pv[r + (size_t)(c++) * nr] = (double)getIndLambdaYj(ind);
        pv[r + (size_t)(c++) * nr] = getIndLogitLow(ind);
        pv[r + (size_t)(c++) * nr] = getIndLogitHi(ind);
      }
      for (int k = 0; k < nsens; ++k) {
        double v = lhs[sensIdx[k]];
        if (v != 0.0) {
          sensR[(size_t)i].push_back(r);
          sensK[(size_t)i].push_back(k);
          sensV[(size_t)i].push_back(v);
        }
      }
      r++;
    }
    if (r != rEnd) rowBad[(size_t)i] = 1;
  });
  _tablePoolSolveN++;
  bool rowsOk = true;
  for (int i = 0; i < nsub; ++i) if (rowBad[(size_t)i]) rowsOk = false;
  size_t nnz = 0;
  for (int i = 0; i < nsub; ++i) nnz += sensV[(size_t)i].size();
  IntegerVector si(nnz), sk(nnz);
  NumericVector sv(nnz);
  size_t o = 0;
  for (int i = 0; i < nsub; ++i) {
    for (size_t q = 0; q < sensV[(size_t)i].size(); ++q, ++o) {
      si[o] = sensR[(size_t)i][q] + 1;
      sk[o] = sensK[(size_t)i][q] + 1;
      sv[o] = sensV[(size_t)i][q];
    }
  }
  return List::create(_["rowsOk"] = rowsOk, _["id"] = id, _["rec"] = rec, _["evid"] = evid, _["time"] = time,
                      _["vals"] = vals,
                      _["sens"] = List::create(_["i"] = si, _["k"] = sk, _["x"] = sv,
                                               _["nrow"] = nr, _["ncol"] = nsens));
}
