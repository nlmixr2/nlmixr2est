#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "utilc.h"
#include "censEst.h"
#include "nearPD.h"
#include "shi21.h"
#include "inner.h"

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("nlmixr2est", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#define popedOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)

extern void doAssignFn(void);
extern rxSolveF rxInner;
extern void rxUpdateFuns(SEXP trans, rxSolveF *inner);
extern void rxClearFuns(rxSolveF *inner);
extern ind_solve_t ind_solve;
extern getRxSolve_t getRx;
extern rx_solve *rx;
extern getTime_t getTimeF;
extern iniSubjectI_t iniSubjectI;
extern isRstudio_t isRstudio;


struct popedOptions {
  int ntheta=0;
  int stickyTol=0;
  int stickyRecalcN=1;
  int stickyRecalcN2=0;
  int stickyRecalcN1=0;
  int maxOdeRecalc;
  int reducedTol;
  int reducedTol2;
  int naZero;
  double odeRecalcFactor;
  bool loaded=false;
};

popedOptions popedOp;

//[[Rcpp::export]]
RObject popedFree() {
  return R_NilValue;
}

Environment _popedE;

//[[Rcpp::export]]
RObject popedSetup(Environment e, bool full) {
  doAssignFn();
  popedFree();
  _popedE=e;
  List control = e["control"];
  List rxControl = as<List>(e["rxControl"]);

  RObject model;
  NumericVector p;
  RObject data;
  if (full) {
    model = e["modelF"];
    p = as<NumericVector>(e["paramF"]);
    data = e["dataF"]; //const RObject &events =
  } else {
    model = e["modelMT"];
    p = as<NumericVector>(e["paramMT"]);
    data = e["dataMT"]; //const RObject &events =
  }
  List mvp = rxode2::rxModelVars_(model);
  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxInner);

  // initial value of parameters
  CharacterVector pars = mvp[RxMv_params];
  popedOp.ntheta = pars.size();
  if (p.size() != popedOp.ntheta) {
    Rprintf("pars\n");
    print(pars);
    Rprintf("p\n");
    print(p);
    Rcpp::stop("size mismatch");
  }
  popedOp.stickyRecalcN=as<int>(control["stickyRecalcN"]);
  popedOp.stickyTol=0;
  popedOp.stickyRecalcN2=0;
  popedOp.stickyRecalcN1=0;
  popedOp.reducedTol = 0;
  popedOp.reducedTol2 = 0;
  popedOp.naZero=0;
  popedOp.maxOdeRecalc = as<int>(control["maxOdeRecalc"]);
  popedOp.odeRecalcFactor = as<double>(control["odeRecalcFactor"]);
  rxode2::rxSolve_(model, rxControl,
                   R_NilValue,//const Nullable<CharacterVector> &specParams =
                   R_NilValue,//const Nullable<List> &extraArgs =
                   p,//const RObject &params =
                   data,//const RObject &events =
                   R_NilValue, // inits
                   1);//const int setupOnly = 0
  rx = getRx();
  return R_NilValue;
}

void popedSolve(int &id) {
  rx_solving_options *op = rx->op;
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  popedOde(id);
  int j=0;
  while (popedOp.stickyRecalcN2 <= popedOp.stickyRecalcN &&
         op->badSolve && j < popedOp.maxOdeRecalc) {
    popedOp.stickyRecalcN2++;
    popedOp.reducedTol2 = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(popedOp.odeRecalcFactor);
    ind->solved = -1;
    popedOde(id);
    j++;
  }
  if (j != 0) {
    if (popedOp.stickyRecalcN2 <= popedOp.stickyRecalcN){
      // Not thread safe
      rxode2::atolRtolFactor_(pow(popedOp.odeRecalcFactor, -j));
    } else {
      popedOp.stickyTol=1;
    }
  }
}

static inline rx_solving_options_ind* updateParamRetInd(NumericVector &theta, int &id) {
  rx = getRx();
  rx_solving_options_ind *ind = &(rx->subjects[id]);
  for (int i = popedOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  return ind;
}

// Solve prediction and saved based on modeling time
void popedSolveFid(double *f, double *w, double *t, NumericVector &theta, int id, int totn) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = rx->op;
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < ind->n_all_times; ++j) {
    ind->idx=j;
    kk = ind->ix[j];
    curT = getTimeF(kk, ind);
    if (isDose(ind->evid[kk])) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      continue;
    } else if (ind->evid[kk] == 0) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      if (ISNA(ind->lhs[0])) {
        popedOp.naZero=1;
        ind->lhs[0] = 0.0;
      }
      // ret(k) = ind->lhs[0];
      // k++;
    } else if (ind->evid[kk] >= 10 && ind->evid[kk] <= 99) {
      // mtimes to calculate information
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      f[k] = ind->lhs[0];
      w[k] = sqrt(ind->lhs[1]);
      t[k] = curT;
      k++;
      if (k >= totn) return; // vector has been created, break
    }
  }
}

void popedSolveFid2(double *f, double *w, double *t, NumericVector &theta, int id, int totn) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = rx->op;
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < ind->n_all_times; ++j) {
    ind->idx=j;
    kk = ind->ix[j];
    curT = getTimeF(kk, ind);
    if (isDose(ind->evid[kk])) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      continue;
    } else if (ind->evid[kk] == 0) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      if (ISNA(ind->lhs[0])) {
        popedOp.naZero=1;
        ind->lhs[0] = 0.0;
      }
      // ret(k) = ind->lhs[0];
      // k++;
      f[k] = ind->lhs[0];
      w[k] = sqrt(ind->lhs[1]);
      t[k] = curT;
      k++;
      if (k >= totn) return; // vector has been created, break
    } else if (ind->evid[kk] >= 10 && ind->evid[kk] <= 99) {
      // mtimes to calculate information
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
    }
  }
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdN2(NumericVector &theta, NumericVector &mt, int id, int totn) {
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  popedSolveFid2(&f[0], &w[0], &t[0], theta, id, totn);
  DataFrame ret = DataFrame::create(_["t"]=mt,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;
  return ret;
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdN(NumericVector &theta, NumericVector &mt, int id, int totn) {
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  popedSolveFid(&f[0], &w[0], &t[0], theta, id, totn);
  // arma::uvec m = as<arma::uvec>(match(mt, t))-1;
  // f = f(m);
  // w = w(m);
  DataFrame ret = DataFrame::create(_["t"]=mt,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;
  return ret;
}
