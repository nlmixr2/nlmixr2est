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

// each subject's first output row (start[nsub] = total).  Counted by record,
// not through ix: the time order is only sorted by the solve.
static std::vector<int> tablePoolRowStarts(rx_solve *rx, int nsub) {
  std::vector<int> start((size_t)nsub + 1, 0);
  for (int i = 0; i < nsub; ++i) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, i);
    int n = 0, nAll = getIndNallTimes(ind);
    for (int j = 0; j < nAll; ++j) {
      if (tablePoolOutRow(getIndEvid(ind, j))) n++;
    }
    start[(size_t)i + 1] = start[(size_t)i] + n;
  }
  return start;
}

static void tablePoolCheckIdx(IntegerVector idx, int n, const char *what) {
  for (int k = 0; k < idx.size(); ++k) {
    if (idx[k] < 0 || idx[k] >= n) stop("table pool %s index out of range", what);
  }
}

// output buffers shared by every subject; each subject writes only its own rows
struct TablePoolOut {
  int nr, ncol, nsens;
  int *pid, *prec, *pevid;
  double *ptime, *pv;
  const int *lhsIdx, *stateIdx, *covIdx, *sensIdx;
  int nl, ns, nc;
  bool tbs;
  std::vector<std::vector<int> > sensR, sensK;
  std::vector<std::vector<double> > sensV;
};

static inline void tablePoolSens(TablePoolOut &o, int i, int r, int k, double v) {
  o.sensR[(size_t)i].push_back(r);
  o.sensK[(size_t)i].push_back(k);
  o.sensV[(size_t)i].push_back(v);
}

// a failed solve's values and sensitivities are missing, not structural zeros
static void tablePoolBadRow(TablePoolOut &o, int i, int r) {
  for (int k = 0; k < o.ncol; ++k) o.pv[r + (size_t)k * o.nr] = NA_REAL;
  for (int k = 0; k < o.nsens; ++k) tablePoolSens(o, i, r, k, NA_REAL);
}

static void tablePoolGoodRow(TablePoolOut &o, int i, int r, rx_solving_options_ind *ind,
                             const double *lhs, const double *solve) {
  int c = 0;
  double *pv = o.pv;
  size_t nr = (size_t)o.nr;
  for (int k = 0; k < o.nl; ++k) pv[r + (size_t)(c++) * nr] = lhs[o.lhsIdx[k]];
  for (int k = 0; k < o.ns; ++k) pv[r + (size_t)(c++) * nr] = solve[o.stateIdx[k]];
  for (int k = 0; k < o.nc; ++k) pv[r + (size_t)(c++) * nr] = getIndParPtr(ind, o.covIdx[k]);
  if (o.tbs) {
    pv[r + (size_t)(c++) * nr] = getIndLambda(ind);
    pv[r + (size_t)(c++) * nr] = (double)getIndLambdaYj(ind);
    pv[r + (size_t)(c++) * nr] = getIndLogitLow(ind);
    pv[r + (size_t)(c++) * nr] = getIndLogitHi(ind);
  }
  for (int k = 0; k < o.nsens; ++k) {
    double v = lhs[o.sensIdx[k]];
    if (v != 0.0) tablePoolSens(o, i, r, k, v);
  }
}

// walk subject i's records after its solve, writing rows [r, rEnd); returns
// the row reached, or -1 when the subject has more output rows than counted
static int tablePoolWalkInd(TablePoolOut &o, rxSolveF *fns, rx_solving_options *op,
                            rx_solving_options_ind *ind, int i, int r, int rEnd,
                            double *lhs, bool bad) {
  int nAll = getIndNallTimes(ind);
  for (int j = 0; j < nAll; ++j) {
    setIndIdx(ind, j);
    int kk = getIndIx(ind, j);
    int ev = getIndEvid(ind, kk);
    if (ev == 9) continue;
    double curT = getTime(kk, ind);
    bool out = tablePoolOutRow(ev);
    if (!out && (ev >= 100 || ev == 3)) tablePoolTlast(&curT, ind);
    if (out && r >= rEnd) return -1;
    double *solve = NULL;
    if (!bad) {
      lhsWalkNewind(ind, j);
      solve = getOpIndSolve(op, ind, j);
      fns->calc_lhs(i, curT, solve, lhs);
    }
    if (!out) continue;
    o.pid[r] = i + 1;
    o.prec[r] = kk;
    o.pevid[r] = ev;
    o.ptime[r] = curT;
    if (bad) {
      tablePoolBadRow(o, i, r);
    } else {
      tablePoolGoodRow(o, i, r, ind, lhs, solve);
    }
    r++;
  }
  return r;
}

static List tablePoolSensList(TablePoolOut &o, int nsub) {
  size_t nnz = 0;
  for (int i = 0; i < nsub; ++i) nnz += o.sensV[(size_t)i].size();
  IntegerVector si(nnz), sk(nnz);
  NumericVector sv(nnz);
  size_t q0 = 0;
  for (int i = 0; i < nsub; ++i) {
    for (size_t q = 0; q < o.sensV[(size_t)i].size(); ++q, ++q0) {
      si[q0] = o.sensR[(size_t)i][q] + 1;
      sk[q0] = o.sensK[(size_t)i][q] + 1;
      sv[q0] = o.sensV[(size_t)i][q];
    }
  }
  return List::create(_["i"] = si, _["k"] = sk, _["x"] = sv,
                      _["nrow"] = o.nr, _["ncol"] = o.nsens);
}

// One subject's rows through rxode2's own output-record walk (rxDfWalk*, from
// rxode2's pointer table): the rows, times and dose columns rxSolve() gives for
// the same addDosing/subsetNonmem, and its keep columns.  The rows are known
// only after the solve, so each subject fills its own buffers.
struct TablePoolRxSub {
  std::vector<int> rec, evid, resetno;
  std::vector<double> time, dose, vals, keep;
  std::vector<int> sensR, sensK;
  std::vector<double> sensV;
  int dull[4] = {1, 1, 1, 1};
  bool walkBad = false;
};

// what a walk needs besides the subject
struct TablePoolRxCtx {
  rx_solve *rx;
  rx_solving_options *op;
  rxSolveF *fns;
  int doDose, nkeep, stateSize;
  const TablePoolOut *o;
};

static void tablePoolRxRow(const TablePoolRxCtx &c, TablePoolRxSub &sub, rx_solving_options_ind *ind,
                           int i, int j, const double *wo, const double *lhs, const double *solve,
                           bool bad, int fid, double *keepBuf) {
  const TablePoolOut &o = *c.o;
  int r = (int)sub.time.size();
  sub.rec.push_back(getIndIx(ind, j));
  sub.evid.push_back((int)wo[RXDF_OUT_EVID]);
  sub.resetno.push_back((int)wo[RXDF_OUT_RESETNO]);
  sub.time.push_back(wo[RXDF_OUT_TIME]);
  if (c.doDose > 0) {
    for (int k = RXDF_OUT_CMT; k <= RXDF_OUT_II; ++k) sub.dose.push_back(wo[k]);
  }
  if (bad) {
    for (int k = 0; k < o.ncol; ++k) sub.vals.push_back(NA_REAL);
    for (int k = 0; k < o.nsens; ++k) {
      sub.sensR.push_back(r);
      sub.sensK.push_back(k);
      sub.sensV.push_back(NA_REAL);
    }
  } else {
    for (int k = 0; k < o.nl; ++k) sub.vals.push_back(lhs[o.lhsIdx[k]]);
    for (int k = 0; k < o.ns; ++k) sub.vals.push_back(solve[o.stateIdx[k]]);
    for (int k = 0; k < o.nc; ++k) sub.vals.push_back(getIndParPtr(ind, o.covIdx[k]));
    if (o.tbs) {
      sub.vals.push_back(getIndLambda(ind));
      sub.vals.push_back((double)getIndLambdaYj(ind));
      sub.vals.push_back(getIndLogitLow(ind));
      sub.vals.push_back(getIndLogitHi(ind));
    }
    for (int k = 0; k < o.nsens; ++k) {
      double v = lhs[o.sensIdx[k]];
      if (v != 0.0) {
        sub.sensR.push_back(r);
        sub.sensK.push_back(k);
        sub.sensV.push_back(v);
      }
    }
  }
  if (c.nkeep > 0) {
    nmDfKeepRow(c.rx, ind, i, fid, j, wo[RXDF_OUT_LHS_TIME], keepBuf);
    for (int k = 0; k < c.nkeep; ++k) sub.keep.push_back(keepBuf[k]);
  }
}

// walk subject i after its solve, calling calc_lhs where rxode2_df() does
static void tablePoolRxWalkInd(const TablePoolRxCtx &c, TablePoolRxSub &sub, rx_solving_options_ind *ind,
                               int i, double *lhs, bool bad) {
  std::vector<char> state((size_t)c.stateSize);
  if (nmDfWalkInit(state.data(), c.stateSize, c.rx, ind, c.doDose, 0) != 1) {
    sub.walkBad = true;
    return;
  }
  std::vector<double> keepBuf((size_t)c.nkeep + 1);
  int fid = c.nkeep > 0 ? nmDfKeepFirst(c.rx, i) : 0;
  double wo[RXDF_OUT_N];
  int newind = 1;
  int nAll = getIndNallTimes(ind);
  for (int j = 0; j < nAll; ++j) {
    int wr = nmDfWalkStep(state.data(), j, wo, RXDF_OUT_N);
    if (wr == RXDF_SKIP) continue;
    if (wr == RXDF_BAD) {
      sub.walkBad = true;
      return;
    }
    double *solve = NULL;
    if (!bad) {
      if (nmSetIndNewind != NULL) nmSetIndNewind(ind, newind);
      solve = getOpIndSolve(c.op, ind, j);
      c.fns->calc_lhs(i, wo[RXDF_OUT_LHS_TIME], solve, lhs);
    }
    if (wo[RXDF_OUT_NEWIND] == 1.0) newind = 2;
    if (wr == RXDF_ROW) tablePoolRxRow(c, sub, ind, i, j, wo, lhs, solve, bad, fid, keepBuf.data());
  }
  nmDfWalkDull(state.data(), sub.dull);
}

// the subjects' rows, one after the other, as the solve's result list
static List tablePoolRxCollect(const TablePoolRxCtx &c, std::vector<TablePoolRxSub> &subs) {
  const TablePoolOut &o = *c.o;
  int nsub = (int)subs.size();
  size_t nr = 0, nnz = 0;
  bool rowsOk = true;
  int dull[4] = {1, 1, 1, 1};
  for (int i = 0; i < nsub; ++i) {
    nr += subs[(size_t)i].time.size();
    nnz += subs[(size_t)i].sensV.size();
    if (subs[(size_t)i].walkBad) rowsOk = false;
    for (int k = 0; k < 4; ++k) dull[k] = dull[k] && subs[(size_t)i].dull[k];
  }
  IntegerVector id(nr), rec(nr), evid(nr), resetno(nr), si(nnz), sk(nnz);
  NumericVector time(nr), sv(nnz);
  NumericMatrix vals(nr, o.ncol);
  int ndose = c.doDose > 0 ? 6 : 0;
  NumericMatrix dose(nr, ndose);
  NumericMatrix keep(nr, c.nkeep);
  size_t r0 = 0, q0 = 0;
  for (int i = 0; i < nsub; ++i) {
    TablePoolRxSub &sub = subs[(size_t)i];
    size_t n = sub.time.size();
    for (size_t r = 0; r < n; ++r) {
      size_t rr = r0 + r;
      id[rr] = i + 1;
      rec[rr] = sub.rec[r];
      evid[rr] = sub.evid[r];
      resetno[rr] = sub.resetno[r];
      time[rr] = sub.time[r];
      for (int k = 0; k < o.ncol; ++k) vals[rr + (size_t)k * nr] = sub.vals[r * (size_t)o.ncol + k];
      for (int k = 0; k < ndose; ++k) dose[rr + (size_t)k * nr] = sub.dose[r * 6 + k];
      for (int k = 0; k < c.nkeep; ++k) keep[rr + (size_t)k * nr] = sub.keep[r * (size_t)c.nkeep + k];
    }
    for (size_t q = 0; q < sub.sensV.size(); ++q, ++q0) {
      si[q0] = (int)(sub.sensR[q] + r0) + 1;
      sk[q0] = sub.sensK[q] + 1;
      sv[q0] = sub.sensV[q];
    }
    r0 += n;
  }
  List ret = List::create(_["rowsOk"] = rowsOk, _["id"] = id, _["rec"] = rec, _["evid"] = evid,
                          _["time"] = time, _["vals"] = vals,
                          _["sens"] = List::create(_["i"] = si, _["k"] = sk, _["x"] = sv,
                                                   _["nrow"] = (int)nr, _["ncol"] = o.nsens));
  ret["resetno"] = resetno;
  ret["dose"] = dose;
  ret["dull"] = IntegerVector::create(_["rate"] = dull[0], _["dur"] = dull[1],
                                      _["ss"] = dull[2], _["ii"] = dull[3]);
  ret["keep"] = c.nkeep > 0 ? RObject(nmDfKeepCols(keep)) : RObject(R_NilValue);
  return ret;
}

// tablePoolSolve_() through rxode2's walk: every subject solved, then walked
static List tablePoolSolveRx(int slot, rxSolveF *fns, rx_solve *rx, rx_solving_options *op,
                             const TablePoolOut &o, int doDose, bool keep,
                             const double *ppar, int npr, int npc, const int *pparIdx) {
  int nsub = (int)getRxNsub(rx);
  TablePoolRxCtx c;
  c.rx = rx;
  c.op = op;
  c.fns = fns;
  c.doDose = doDose;
  c.nkeep = keep ? nmDfKeepSetup(rx) : 0;
  c.stateSize = nmDfWalkInit(NULL, 0, rx, NULL, doDose, 0);
  c.o = &o;
  std::vector<TablePoolRxSub> subs((size_t)nsub);
  int cores = getOpCores(op);
  bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
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
    tablePoolRxWalkInd(c, subs[(size_t)i], ind, i, neqGuard.lhs(), bad);
  });
  _tablePoolSolveN++;
  return tablePoolRxCollect(c, subs);
}

// The factor levels rxode2 gives the output columns `names` (rxode2's
// getDfLevels(), which rxSolve() makes a string lhs or a covariate column
// with), NULL for a numeric one; the solve must still be set up.
//[[Rcpp::export]]
List tablePoolLevels_(CharacterVector names) {
  int n = names.size();
  List ret(n);
  ret.names() = names;
  rx_solve *rx = getRxSolve_();
  if (nmDfLevels == NULL || rx == NULL) return ret;
  for (int i = 0; i < n; ++i) {
    SEXP t = PROTECT(nmDfLevels(CHAR(STRING_ELT(names, i)), rx, 0));
    if (Rf_isFactor(t)) ret[i] = Rf_getAttrib(t, R_LevelsSymbol);
    UNPROTECT(1);
  }
  return ret;
}

// Whether rxode2 provides the output-record walk and keep columns
//[[Rcpp::export]]
bool tablePoolHasRxWalk_() {
  return nmDfWalkInit != NULL && nmDfWalkStep != NULL && nmDfWalkDull != NULL &&
    nmDfKeepSetup != NULL && nmDfKeepFirst != NULL && nmDfKeepRow != NULL &&
    nmDfKeepCols != NULL && nmDfLevels != NULL && nmSetIndNewind != NULL;
}

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
// doDose:    NA walks the records itself (observation rows only).  Otherwise
//            rxSolve()'s addDosing code, and the rows, dose columns and (with
//            `keep`) keep columns come from rxode2's own walk; that needs
//            tablePoolHasRxWalk_().
//[[Rcpp::export]]
List tablePoolSolve_(int which, NumericMatrix parMat, IntegerVector parIdx,
                     IntegerVector lhsIdx, IntegerVector sensIdx,
                     IntegerVector stateIdx, IntegerVector covIdx, bool tbs,
                     int doDose, bool keep) {
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
  tablePoolCheckIdx(stateIdx, odeSwapNeq(slot), "state");
  tablePoolCheckIdx(lhsIdx, odeSwapNlhs(slot), "lhs");
  tablePoolCheckIdx(sensIdx, odeSwapNlhs(slot), "lhs");
  bool rxWalk = doDose != NA_INTEGER;
  if (rxWalk && !tablePoolHasRxWalk_()) stop("this rxode2 has no output-record walk");
  std::vector<int> start = rxWalk ? std::vector<int>((size_t)nsub + 1, 0) : tablePoolRowStarts(rx, nsub);
  int nr = start[(size_t)nsub];
  IntegerVector id(nr), rec(nr), evid(nr);
  NumericVector time(nr);
  TablePoolOut o;
  o.nr = nr;
  o.nl = lhsIdx.size();
  o.ns = stateIdx.size();
  o.nc = covIdx.size();
  o.nsens = sensIdx.size();
  o.tbs = tbs;
  o.ncol = o.nl + o.ns + o.nc + (tbs ? 4 : 0);
  NumericMatrix vals(nr, o.ncol);
  o.pv = REAL(vals);
  o.pid = INTEGER(id);
  o.prec = INTEGER(rec);
  o.pevid = INTEGER(evid);
  o.ptime = REAL(time);
  o.lhsIdx = INTEGER(lhsIdx);
  o.stateIdx = INTEGER(stateIdx);
  o.covIdx = INTEGER(covIdx);
  o.sensIdx = INTEGER(sensIdx);
  o.sensR.resize((size_t)nsub);
  o.sensK.resize((size_t)nsub);
  o.sensV.resize((size_t)nsub);
  const double *ppar = REAL(parMat);
  int npr = parMat.nrow(), npc = parMat.ncol();
  const int *pparIdx = INTEGER(parIdx);
  if (rxWalk) return tablePoolSolveRx(slot, fns, rx, op, o, doDose, keep, ppar, npr, npc, pparIdx);
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
    int rEnd = start[(size_t)i + 1];
    int r = tablePoolWalkInd(o, fns, op, ind, i, start[(size_t)i], rEnd, neqGuard.lhs(), bad);
    if (r != rEnd) rowBad[(size_t)i] = 1;
  });
  _tablePoolSolveN++;
  bool rowsOk = true;
  for (int i = 0; i < nsub; ++i) if (rowBad[(size_t)i]) rowsOk = false;
  return List::create(_["rowsOk"] = rowsOk, _["id"] = id, _["rec"] = rec, _["evid"] = evid,
                      _["time"] = time, _["vals"] = vals, _["sens"] = tablePoolSensList(o, nsub));
}
