// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "utilc.h"
#include "censEst.h"
#include "nearPD.h"
#include "shi21.h"
#include "inner.h"


#define _(String) (String)

#include "scale.h"


#define nlmOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
#define predOde(id) ind_solve(rx, id, rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)

struct nlmOptions {
  int ntheta=0;
  int *thetaFD=NULL; // theta needs finite difference?
  int *nobs = NULL;
  int *idS  = NULL;
  int *idF  = NULL;
  int *xPar = NULL;
  int nobsTot = 0;
  double *thetahf=NULL; // Shi step size
  double *thetahh=NULL;
  double *initPar= NULL; // initial parameters
  double *thetaSave = NULL;
  double *valSave   = NULL;
  double *grSave    = NULL;
  double *hSave     = NULL;
  double *scaleC    = NULL;
  double *logitThetaLow = NULL;
  double *logitThetaHi  = NULL;

  int eventType=3; // eventType
  int shi21maxFD=1000; //maxiter for shi
  int stickyTol=0;
  int stickyRecalcN=1;
  int stickyRecalcN2=0;
  int stickyRecalcN1=0;
  int maxOdeRecalc;
  int reducedTol;
  int reducedTol2;
  double odeRecalcFactor;
  bool needFD=false;
  int optimHessType=3;
  int shi21maxHess=1000;
  double shiErr;
  double hessErr;
  int solveType;
#define solveType_none 0
#define solveType_pred 1
#define solveType_grad 2
#define solveType_hess 3
#define solveType_nls 10
#define solveType_nls_pred 11
  int saveType;
#define save_pred 1
#define save_grad 2
#define save_hess 3
  int naZero;
  int naGrad;
  scaling scale;
  bool loaded=false;
};

nlmOptions nlmOp;

//[[Rcpp::export]]
RObject nlmFree() {
  if (nlmOp.thetaFD != NULL) R_Free(nlmOp.thetaFD);
  nlmOp.thetaFD = NULL;
  nlmOp.nobs    = NULL;
  nlmOp.idS     = NULL;
  nlmOp.idF     = NULL;
  nlmOp.xPar    = NULL;
  if (nlmOp.thetahf != NULL) R_Free(nlmOp.thetahf);
  nlmOp.thetahf = NULL;
  nlmOp.thetahh = NULL;
  nlmOp.thetaSave = NULL;
  nlmOp.valSave = NULL;
  nlmOp.grSave = NULL;
  nlmOp.hSave = NULL;
  nlmOp.initPar = NULL;
  nlmOp.scaleC  = NULL;
  nlmOp.logitThetaLow = NULL;
  nlmOp.logitThetaHi  = NULL;

  nlmOp.loaded = false;
  return R_NilValue;
}

//[[Rcpp::export]]
RObject nlmSetup(Environment e) {
  nlmFree();
  List control = e["control"];

  RObject pred = e["predOnly"];
  List mvp = rxode2::rxModelVars_(pred);
  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);

  nlmOp.solveType = as<int>(control["solveType"]);
  RObject model;
  if (e.exists("thetaGrad")) {
    model = e["thetaGrad"];
    List mv = rxode2::rxModelVars_(model);
    rxUpdateFuns(as<SEXP>(mv["trans"]), &rxInner);
  } else {
    if (nlmOp.solveType != solveType_nls_pred) {
      nlmOp.solveType = solveType_pred;
    }
    model = pred;
  }

  List rxControl = as<List>(e["rxControl"]);

  NumericVector p = as<NumericVector>(e["param"]);
  nlmOp.ntheta = p.size();

  nlmOp.stickyRecalcN=as<int>(control["stickyRecalcN"]);
  nlmOp.stickyTol=0;
  nlmOp.stickyRecalcN2=0;
  nlmOp.stickyRecalcN1=0;
  nlmOp.reducedTol = 0;
  nlmOp.reducedTol2 = 0;
  nlmOp.naZero=0;
  nlmOp.saveType = 0;
  nlmOp.naGrad=0;
  nlmOp.maxOdeRecalc = as<int>(control["maxOdeRecalc"]);
  nlmOp.odeRecalcFactor = as<double>(control["odeRecalcFactor"]);

  nlmOp.eventType = control["eventType"];
  nlmOp.shi21maxFD = control["shi21maxFD"];
  nlmOp.shiErr = control["shiErr"];

  nlmOp.optimHessType = control["optimHessType"];
  nlmOp.shi21maxHess = control["shi21maxHess"];
  nlmOp.hessErr = control["hessErr"];


  rxode2::rxSolve_(model, rxControl,
                   R_NilValue,//const Nullable<CharacterVector> &specParams =
                   R_NilValue,//const Nullable<List> &extraArgs =
                   p,//const RObject &params =
                   e["data"],//const RObject &events =
                   R_NilValue, // inits
                   1);//const int setupOnly = 0
  rx = getRxSolve_();

  nlmOp.thetaFD = R_Calloc(nlmOp.ntheta*2 + getRxNsub(rx)*3, int); // [ntheta]
  nlmOp.nobs = nlmOp.thetaFD + nlmOp.ntheta; // [nsub]
  nlmOp.idS = nlmOp.nobs + getRxNsub(rx); // [nsub]
  nlmOp.idF = nlmOp.idS + getRxNsub(rx); // [nsub]
  nlmOp.xPar = nlmOp.idF + getRxNsub(rx); // [ntheta]

  // now calculate nobs per id
  nlmOp.nobsTot = 0;
  for (int id = 0; id < getRxNsub(rx); ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
    int no = 0;
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      if (getIndEvid(ind, j) == 0) {
        nlmOp.nobsTot++;
        no++;
      }
    }
    nlmOp.nobs[id]=no;
    if (id == 0) {
      nlmOp.idS[0] = 0;
      nlmOp.idF[0] = no-1;
      continue;
    }
    nlmOp.idS[id] = nlmOp.idF[id-1]+1;
    nlmOp.idF[id] = nlmOp.idS[id] + no - 1;
  }

  IntegerVector needFD = as<IntegerVector>(e["needFD"]);


  // nlmOp.ntheta nlmOp.ntheta+1
  switch(nlmOp.solveType) {
  case solveType_nls:
    nlmOp.thetahf = R_Calloc(nlmOp.ntheta*(5+getRxNsub(rx)) + nlmOp.nobsTot*(1+nlmOp.ntheta), double);// [ntheta*nsub]
    nlmOp.thetaSave = nlmOp.thetahf + nlmOp.ntheta*getRxNsub(rx); // [ntheta]
    nlmOp.initPar = nlmOp.thetaSave + nlmOp.ntheta; // [ntheta]
    nlmOp.scaleC  = nlmOp.initPar   + nlmOp.ntheta; // [ntheta]
    nlmOp.logitThetaLow = nlmOp.scaleC + nlmOp.ntheta; // [ntheta]
    nlmOp.logitThetaHi  = nlmOp.logitThetaLow + nlmOp.ntheta; // [ntheta]
    nlmOp.valSave = nlmOp.logitThetaHi + nlmOp.ntheta; //[nlmOp.nobsTot]
    nlmOp.grSave  = nlmOp.valSave + nlmOp.nobsTot; // [nlmOp.nobsTot*ntheta]
    break;
  case solveType_nls_pred:
    nlmOp.thetahf = R_Calloc(nlmOp.ntheta*(4+getRxNsub(rx)), double);// [ntheta*nsub]
    nlmOp.initPar = nlmOp.thetahf + nlmOp.ntheta*getRxNsub(rx); // [ntheta]
    nlmOp.scaleC  = nlmOp.initPar   + nlmOp.ntheta; // [ntheta]
    nlmOp.logitThetaLow = nlmOp.scaleC + nlmOp.ntheta; // [ntheta]
    nlmOp.logitThetaHi  = nlmOp.logitThetaLow + nlmOp.ntheta; // [ntheta]
    break;
  default:
    // 7*ntheta + nsub*ntheta + 1 + ntheta*ntheta
    // ntheta*(7+nsub+ntheta) + 1
#define ntheta nlmOp.ntheta
#define nsub getRxNsub(rx)
    //nsub*ntheta
    nlmOp.thetahf = R_Calloc(ntheta*(nsub + 7 + ntheta) + 1, double); //[nsub*ntheta]
    nlmOp.thetahh = nlmOp.thetahf   + ntheta*nsub; // [ntheta]
    nlmOp.thetaSave = nlmOp.thetahh + ntheta; // [ntheta]
    nlmOp.valSave = nlmOp.thetaSave + ntheta; // [1]
    nlmOp.grSave = nlmOp.valSave + 1; // [ntheta]
    nlmOp.hSave = nlmOp.grSave + ntheta;// [ntheta*ntheta]
    nlmOp.initPar = nlmOp.hSave + ntheta*ntheta; // [ntheta]
    nlmOp.scaleC  = nlmOp.initPar   + ntheta; // [ntheta]
    nlmOp.logitThetaLow = nlmOp.scaleC + ntheta; // [ntheta]
    nlmOp.logitThetaHi  = nlmOp.logitThetaLow + ntheta; // [ntheta]
#undef ntheta
#undef nsub

    std::fill_n(nlmOp.thetaSave, nlmOp.ntheta, R_PosInf); // not likely to be equal
  }

  std::copy(&p[0], &p[0] + nlmOp.ntheta, nlmOp.initPar);

  scaleSetup(&(nlmOp.scale),
             nlmOp.initPar,
             nlmOp.scaleC,
             nlmOp.xPar,
             nlmOp.logitThetaLow,
             nlmOp.logitThetaHi,
             as<CharacterVector>(e["thetaNames"]) ,
             as<int>(control["useColor"]),
             as<int>(control["printNcol"]),
             as<int>(control["print"]),
             as<int>(control["normType"]),
             as<int>(control["scaleType"]),
             as<double>(control["scaleCmin"]),
             as<double>(control["scaleCmax"]),
             as<double>(control["scaleTo"]),
             nlmOp.ntheta);
  nlmOp.needFD=false;
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    nlmOp.thetaFD[i] = needFD[i];
    if (nlmOp.thetaFD[i]) {
      nlmOp.needFD=true;
    }
  }

  nlmOp.loaded = true;
  return R_NilValue;
}

//[[Rcpp::export]]
RObject nlmScalePar(RObject p0) {
  if (p0.sexp_type() != REALSXP) {
    return p0;
  }
  NumericVector p = as<NumericVector>(p0);
  if (p.size() != nlmOp.ntheta) stop("parameter dimension mismatch");
  NumericVector ret(nlmOp.ntheta);
  for (int i = 0; i < nlmOp.ntheta; i++) {
    ret[i] = scaleScalePar(&(nlmOp.scale), &p[0], i);
  }
  return ret;
}

//[[Rcpp::export]]
NumericVector nlmUnscalePar(NumericVector p) {
  if (p.size() != nlmOp.ntheta) stop("parameter dimension mismatch");
  NumericVector ret(nlmOp.ntheta);
  for (int i = 0; i < nlmOp.ntheta; i++) {
    ret[i] = scaleUnscalePar(&(nlmOp.scale), &p[0], i);
  }
  ret.attr("names") = p.attr("names");
  return ret;
}

void nlmSolveNlm(int id) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  nlmOde(id);
  int j=0;
  while (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN &&
         hasOpBadSolve(op) && j < nlmOp.maxOdeRecalc) {
    nlmOp.stickyRecalcN2++;
    nlmOp.reducedTol  = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(nlmOp.odeRecalcFactor);
    setIndSolve(ind, -1);
    nlmOde(id);
    j++;
  }
  if (j != 0) {
    if (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN){
      // Not thread safe
      rxode2::atolRtolFactor_(pow(nlmOp.odeRecalcFactor, -j));
    } else {
      nlmOp.stickyTol=1;
    }
  }
}

void nlmSolvePred(int &id) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  predOde(id);
  int j=0;
  while (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN &&
         hasOpBadSolve(op) && j < nlmOp.maxOdeRecalc) {
    nlmOp.stickyRecalcN2++;
    nlmOp.reducedTol2 = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(nlmOp.odeRecalcFactor);
    setIndSolve(ind, -1);
    predOde(id);
    j++;
  }
  if (j != 0) {
    if (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN){
      // Not thread safe
      rxode2::atolRtolFactor_(pow(nlmOp.odeRecalcFactor, -j));
    } else {
      nlmOp.stickyTol=1;
    }
  }
}

extern arma::vec calcGradForward(arma::vec &f0, arma::vec &grPH,  double h);
extern arma::vec calcGradCentral(arma::vec &grMH, arma::vec &f0, arma::vec &grPH,  double h);

static inline rx_solving_options_ind* updateParamRetInd(arma::vec &theta, int &id) {
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  for (int i = nlmOp.ntheta; i--;) {
    setIndParPtr(ind, i, scaleUnscalePar(&(nlmOp.scale), &theta[0], i));
  }
  return ind;
}

static inline bool isThetaSame(arma::vec &theta) {
  arma::vec thetaSave(nlmOp.thetaSave, nlmOp.ntheta, false, true);
  if (arma::approx_equal(theta, thetaSave, "absdiff",  std::numeric_limits<double>::epsilon())) return true;
  return false;
}

static inline void saveTheta(arma::vec &theta) {
  arma::vec thetaSave(nlmOp.thetaSave, nlmOp.ntheta, false, true);
  thetaSave = theta;
}

// Solve prediction
void nlmSolveFid(double *retD, int nobs, arma::vec &theta, int id) {
  arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  iniSubjectE(id, 1, ind, op, rx, rxPred.update_inis);
  nlmSolvePred(id);
  int kk, k=0;
  double curT;
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
      if (ISNA(lhs[0])) {
        nlmOp.naZero=1;
        lhs[0] = 0.0;
      }
      ret(k) = lhs[0];
      k++;
    }
  }
}

arma::vec nlmSolveFid(arma::vec &theta, int id) {
  arma::vec ret(nlmOp.nobs[id]);
  double *retD = ret.memptr();
  nlmSolveFid(retD, nlmOp.nobs[id], theta, id);
  return ret;
}


arma::vec nlmSolveF(arma::vec &theta) {
  arma::vec ret(nlmOp.nobsTot);
  double *retD = ret.memptr();
  rx_solving_options *op = getSolvingOptions(rx);
  int cores = getOpCores(op);
  // #ifdef _OPENMP
  // #pragma omp parallel for num_threads(cores)
  // #endif
  for (int id = 0; id < getRxNsub(rx); ++id) {
    nlmSolveFid(retD + nlmOp.idS[id], nlmOp.nobs[id], theta, id);
  }
  return ret;
}

//[[Rcpp::export]]
double nlmSolveR(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  arma::vec ret = nlmSolveF(theta);
  return arma::sum(ret);
}

// Solve gradient with finite difference thetas This is done for every
// observation so it could possibly be used in nls types of models
arma::mat nlmSolveGradId(arma::vec &theta, int id) {
  // first solve the nlm problem
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  arma::mat ret(nlmOp.nobs[id], nlmOp.ntheta+1);
  int kk, k=0;
  double curT;
  iniSubjectE(id, 1, ind, op, rx, rxInner.update_inis);
  nlmSolveNlm(id);
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      for (int kk = 0; kk < getOpNlhs(op); ++kk) {
        if (ISNA(lhs[kk])) {
          lhs[kk] = 0.0;
          nlmOp.naZero=1;
        }
        if (kk == 0) {
          ret(k, kk) = lhs[kk];
        } else {
          ret(k, kk) = scaleAdjustGradScale(&(nlmOp.scale), lhs[kk], &theta[0], kk-1);
        }
      }
      k++;
    }
    if (k >= getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind)) {
      // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
      break;
    }
  }
  double *thetahf = nlmOp.thetahf + id*nlmOp.ntheta;

  arma::vec f0 = ret.col(0);
  arma::vec grTheta(nlmOp.nobs[id]);
  arma::vec grPH(nlmOp.nobs[id]);
  arma::vec grMH(nlmOp.nobs[id]);

  arma::vec hTheta(nlmOp.ntheta);
  arma::vec curTheta = theta;
  for (int ii = 0; ii < nlmOp.ntheta; ++ii) {
    if (nlmOp.thetaFD[ii] == 0) {
      if (!ret.col(ii+1).has_nan()) {
        continue;
      }
      nlmOp.naGrad=1;
    }
    if (thetahf[ii] == 0.0) {
      double h = 0;
      switch(nlmOp.eventType) {
      case 2: // central
        thetahf[ii] = shi21Central(nlmSolveFid, curTheta, h,
                                   f0, grTheta, id, ii,
                                   nlmOp.shiErr, // ef,
                                   1.5,//double rl = 1.5,
                                   4.5,//double ru = 4.5,
                                   3.0,//double nu = 8.0);
                                   nlmOp.shi21maxFD); // maxiter
        ret.col(ii+1) = grTheta;
        break;
      case 1: // forward
        thetahf[ii] = shi21Forward(nlmSolveFid, curTheta, h,
                                   f0, grTheta, id, ii,
                                   nlmOp.shiErr, // ef,
                                   1.5,  //double rl = 1.5,
                                   6.0,  //double ru = 6.0);;
                                   nlmOp.shi21maxFD); // maxiter
      }
      ret.col(ii+1) = grTheta;
      continue;
    }
    // already calculated optimum thetahf, now do the differences
    hTheta = curTheta;
    hTheta[ii] += nlmOp.thetahf[ii];
    grPH = nlmSolveFid(hTheta, id);
    bool useForward = false;
    if (nlmOp.eventType == 1) {
      // if this isn't true try backward
      if (grPH.is_finite()) {
        useForward = true;
        ret.col(ii+1) = calcGradForward(f0, grPH,  nlmOp.thetahf[ii]);
        continue;
      }
    }
    if (!useForward) {
      // stencil or central
      hTheta = curTheta;
      hTheta[ii] -= nlmOp.thetahf[ii];
      grMH = nlmSolveFid(hTheta, id);
      // central
      ret.col(ii+1) = calcGradCentral(grMH, f0, grPH,  nlmOp.thetahf[ii]);
    }
  }
  // restore save (may not be needed)
  //  }
  return ret;
}

arma::mat nlmSolveGrad(arma::vec &theta) {
  arma::mat ret(nlmOp.nobsTot, nlmOp.ntheta+1);
  rx_solving_options *op = getSolvingOptions(rx);
  int cores = getOpCores(op);
  // #ifdef _OPENMP
  // #pragma omp parallel for num_threads(cores)
  // #endif
  for (int id = 0; id < getRxNsub(rx); ++id) {
    ret.rows(nlmOp.idS[id], nlmOp.idF[id]) = nlmSolveGradId(theta, id);
  }
  return ret;
}

//[[Rcpp::export]]
RObject nlmSetScaleC(NumericVector scaleC) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (scaleC.size() != nlmOp.ntheta) {
    REprintf("ntheta %d\n", nlmOp.ntheta);
    stop("scaleC size mismatch");
  }
  std::copy(scaleC.begin(), scaleC.end(), nlmOp.scaleC);
  return R_NilValue;
}

//[[Rcpp::export]]
NumericVector nlmGetScaleC(arma::vec &theta, double to) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.solveType == solveType_pred)  return NumericVector::create();
  if (to <= 0) return NumericVector::create();
  std::fill_n(nlmOp.scaleC, nlmOp.ntheta, 1.0);
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  NumericVector scaleC(nlmOp.ntheta);
  for(int i = 0; i < nlmOp.ntheta; ++i) {
    // grad*C = to
    //
    scaleC[i] = fabs(to/cs(i+1));
  }
  std::copy(scaleC.begin(), scaleC.end(), nlmOp.scaleC);
  return scaleC;
}



//[[Rcpp::export]]
RObject nlmSolveGradR(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.solveType == solveType_pred) stop("incorrect solve type");
  int ntheta = theta.size();
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  NumericVector ret(1);
  NumericVector grad(ntheta);
  ret[0] = cs[0];
  grad = wrap(cs(span(1, ntheta)));
  ret.attr("gradient") = grad;
  scalePrintFun(&(nlmOp.scale), &theta[0], cs[0]);
  return ret;
}

arma::vec nlmSolveGrad1(arma::vec &theta, int id) {
  arma::mat ret0 = nlmSolveGrad(theta);
  ret0 = ret0.cols(1, nlmOp.ntheta);
  return (arma::sum(ret0, 0)).t();
}

//[[Rcpp::export]]
NumericVector solveGradNls(arma::vec &theta, int returnType) {
  if (!nlmOp.loaded) stop("'nls' problem not loaded");
  if (nlmOp.solveType == solveType_nls_pred) {
    arma::vec resid = nlmSolveF(theta);
    arma::vec r2 = resid % resid;
    double rss = arma::sum(r2);
    scalePrintFun(&(nlmOp.scale), &theta[0], rss);
    return wrap(resid);
  }
  if (nlmOp.solveType != solveType_nls) {
    stop(_("incorrect solve type"));
  }
  if (!isThetaSame(theta)) {
    arma::mat ret0(nlmOp.valSave, nlmOp.nobsTot, nlmOp.ntheta+1, false, true);
    ret0 = nlmSolveGrad(theta);
    if (ret0.has_nan()) {
      nlmOp.naZero=1;
      ret0.replace(datum::nan, 0);
    }
    double llik;
    arma::vec resid =ret0.col(0);
    resid = resid % resid;
    double rss = arma::sum(resid);
    scalePrintFun(&(nlmOp.scale), &theta[0], rss);
    saveTheta(theta);
  }
  if (returnType == 1) {
    // vector with gradient attached
    NumericVector ret(nlmOp.nobsTot);
    NumericVector grad(nlmOp.nobsTot*nlmOp.ntheta);
    std::copy(nlmOp.valSave, nlmOp.valSave + nlmOp.nobsTot, ret.begin());
    std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta*nlmOp.nobsTot, grad.begin());
    grad.attr("dim") = IntegerVector::create(nlmOp.nobsTot, nlmOp.ntheta);
    ret.attr("gradient") = grad;
    return ret;
  } else if (returnType == 2) {
    // vector only
    NumericVector ret(nlmOp.nobsTot);
    std::copy(nlmOp.valSave, nlmOp.valSave + nlmOp.nobsTot, ret.begin());
    return ret;
  } else if (returnType == 3) {
    // gradient only
    NumericVector grad(nlmOp.nobsTot*nlmOp.ntheta);
    std::copy(nlmOp.grSave,
              nlmOp.grSave + nlmOp.ntheta*nlmOp.nobsTot,
              grad.begin());
    grad.attr("dim") = IntegerVector::create(nlmOp.nobsTot, nlmOp.ntheta);
    return grad;
  }
  return NumericVector::create();
}

arma::mat nlmCalcHessian(arma::vec &gr0, arma::vec &theta) {
  int id = 0; // dummy id
  mat H(nlmOp.ntheta, nlmOp.ntheta);
  H.zeros();
  arma::vec grPH(nlmOp.ntheta);
  arma::vec grMH(nlmOp.ntheta);
  double h;
  double *thetahh = nlmOp.thetahh;
  for (int k = nlmOp.ntheta; k--;) {
    h = thetahh[k];
    if (nlmOp.optimHessType == 1 && h <= 0) {
      arma::vec t = theta;
      thetahh[k] = shi21Forward(nlmSolveGrad1, theta, h,
                                gr0, grPH, id, k,
                                nlmOp.hessErr, //double ef = 7e-7,
                                1.5,  //double rl = 1.5,
                                6.0,  //double ru = 6.0);;
                                nlmOp.shi21maxHess);  //maxiter=15
      H.col(k) = grPH;
      continue;
    }
    if (nlmOp.optimHessType == 2 && h <= 0) {
      // Central
      arma::vec t = theta;
      thetahh[k] = shi21Central(nlmSolveGrad1, t, h,
                                gr0, grPH, id, k,
                                nlmOp.hessErr, // ef,
                                1.5,//double rl = 1.5,
                                4.5,//double ru = 4.5,
                                3.0,//double nu = 8.0);
                                nlmOp.shi21maxHess); // maxiter
      H.col(k) = grPH;
      continue;
    }
    // x + h
    theta[k] += h;
    grPH = nlmSolveGrad1(theta, id);
    bool forwardFinite =  grPH.is_finite();
    if (nlmOp.optimHessType == 1 && forwardFinite) { // forward
      H.col(k) = (grPH-gr0)/h;
      theta[k] -= h;
      continue;
    }
    // x - h
    theta[k] -= 2*h;
    grMH = nlmSolveGrad1(theta, 0);
    bool backwardFinite = grMH.is_finite();
    if (nlmOp.optimHessType == 2 &&
        forwardFinite && backwardFinite) {
      // central
      theta[k] += h;
      H.col(k) = (grPH-grMH)/(2.0*h);
      continue;
    }
    if (forwardFinite && !backwardFinite) {
      // forward difference
      H.col(k) = (grPH-gr0)/h;
      theta[k] += h;
      continue;
    }
    if (!forwardFinite && backwardFinite) {
      // backward difference
      H.col(k) = (gr0-grMH)/h;
      theta[k] += h;
      continue;
    }
  }
  // symmetrize
  H = 0.5*(H + H.t());
  return H;
}

//[[Rcpp::export]]
RObject nlmSolveGradHess(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.solveType == solveType_pred) stop("incorrect solve type");
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  double ll = cs[0];
  arma::vec gr0 = cs(span(1, nlmOp.ntheta));
  mat H = nlmCalcHessian(gr0, theta);
  NumericVector ret(1);
  ret[0] = ll;
  NumericVector grad = wrap(gr0(span(0, nlmOp.ntheta-1)));
  ret.attr("gradient") = grad;
  ret.attr("hessian") = wrap(H);
  scalePrintFun(&(nlmOp.scale), &theta[0], ll);
  scalePrintGrad(&(nlmOp.scale), &grad[0], iterTypeSens);
  return ret;
}

//[[Rcpp::export]]
RObject nlmSolveSwitch(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  NumericVector ret;
  switch(nlmOp.solveType) {
  case solveType_pred:
    ret = wrap(nlmSolveR(theta));
    scalePrintFun(&(nlmOp.scale), &theta[0], ret[0]);
    return ret;
  case solveType_grad:
    return nlmSolveGradR(theta);
  case solveType_hess:
    return nlmSolveGradHess(theta);
  }
  return R_NilValue;
}


//[[Rcpp::export]]
NumericVector optimFunC(arma::vec &theta, bool grad=false) {
  if (!nlmOp.loaded) stop("'optim' problem not loaded");
  if (nlmOp.solveType == solveType_pred) {
    if (grad) stop(_("incorrect solve type"));
    NumericVector ret(1);
    ret[0] = nlmSolveR(theta);
    scalePrintFun(&(nlmOp.scale), &theta[0], ret[0]);
    return ret;
  }
  if (isThetaSame(theta)) {
    if (grad) {
      NumericVector ret(nlmOp.ntheta);
      std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, ret.begin());
      return ret;
    }
    NumericVector ret(1);
    ret[0] = nlmOp.valSave[0];
    return ret;
  }
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec saveVec(nlmOp.valSave, nlmOp.ntheta + 1, false, true);
  saveVec = (arma::sum(ret0, 0)).t();
  saveTheta(theta);
  if (grad) {
    NumericVector ret(nlmOp.ntheta);
    std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, ret.begin());
    scalePrintFun(&(nlmOp.scale), &theta[0], nlmOp.valSave[0]);
    scalePrintGrad(&(nlmOp.scale), nlmOp.grSave, iterTypeSens);
    return ret;
  }
  NumericVector ret(1);
  ret[0] = nlmOp.valSave[0];
  scalePrintFun(&(nlmOp.scale), &theta[0], ret[0]);
  scalePrintGrad(&(nlmOp.scale), nlmOp.grSave, iterTypeSens);
  return ret;
}

//[[Rcpp::export]]
NumericVector nlminbFunC(arma::vec &theta, int type) {
  if (!nlmOp.loaded) stop("'nlminb' problem not loaded");
  // restore saved values
  bool isSame = isThetaSame(theta);
  if (isSame) {
    switch (type) {
    case solveType_pred:
      if (nlmOp.saveType >= solveType_pred) {
        NumericVector ret(1);
        ret[0] = nlmOp.valSave[0];
        return ret;
      }
      break;
    case solveType_grad:
      if (nlmOp.saveType >= solveType_grad) {
        NumericVector ret(nlmOp.ntheta);
        std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, ret.begin());
        return ret;
      }
      break;
    case solveType_hess:
      if (nlmOp.saveType == solveType_hess) {
        NumericVector ret(nlmOp.ntheta*nlmOp.ntheta);
        std::copy(nlmOp.hSave, nlmOp.hSave + nlmOp.ntheta * nlmOp.ntheta, ret.begin());
        ret.attr("dim") = IntegerVector::create(nlmOp.ntheta, nlmOp.ntheta);
        return ret;
      }
      break;
    }
  }
  // calculate saved values
  switch (type) {
  case solveType_pred: {
    NumericVector ret(1);
    nlmOp.valSave[0] = ret[0] = nlmSolveR(theta);
    nlmOp.saveType = solveType_pred;
    saveTheta(theta);
    scalePrintFun(&(nlmOp.scale), &theta[0], ret[0]);
    return ret;
  }
    break;
  case solveType_grad: {
    // You have to solve the full system for the grad anyway
    arma::mat ret0 = nlmSolveGrad(theta);
    arma::vec saveVec(nlmOp.valSave, nlmOp.ntheta + 1, false, true);
    saveVec = (arma::sum(ret0, 0)).t();
    nlmOp.saveType = solveType_grad;
    saveTheta(theta);
    NumericVector ret(nlmOp.ntheta);
    std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, ret.begin());
    scalePrintGrad(&(nlmOp.scale), &ret[0], iterTypeSens);
    return ret;
  }
    break;
  case solveType_hess: {
    if (isSame && nlmOp.saveType == solveType_grad) {
      // Just add hessian
      arma::vec gr0(nlmOp.ntheta);
      std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, gr0.begin());
      mat H = nlmCalcHessian(gr0, theta);
      std::copy(H.begin(), H.end(), nlmOp.hSave);
      nlmOp.saveType = solveType_hess;
    } else {
      // Calculate everything
      NumericVector attrRet = as<NumericVector>(nlmSolveGradHess(theta));
      saveTheta(theta);
      nlmOp.saveType = solveType_hess;
      nlmOp.valSave[0] = attrRet[0];
      NumericVector grad = as<NumericVector>(attrRet.attr("gradient"));
      std::copy(grad.begin(), grad.end(), nlmOp.grSave);
      NumericVector hess = as<NumericVector>(attrRet.attr("hessian"));
      std::copy(hess.begin(), hess.end(), nlmOp.hSave);
    }
    NumericVector ret(nlmOp.ntheta*nlmOp.ntheta);
    std::copy(nlmOp.hSave, nlmOp.hSave + nlmOp.ntheta * nlmOp.ntheta, ret.begin());
    ret.attr("dim") = IntegerVector::create(nlmOp.ntheta, nlmOp.ntheta);
    return ret;
  }
    break;
  }
  stop("Couldn't find solution type: %d", type);
  return NumericVector::create(NA_REAL);
}

//[[Rcpp::export]]
RObject nlmPrintHeader() {
  scalePrintHeader(&(nlmOp.scale));
  return R_NilValue;
}

//[[Rcpp::export]]
RObject nlmWarnings() {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.naGrad) {
    warning(_("NaN symbolic gradients were resolved with finite differences"));
  }
  if (nlmOp.naZero) {
    warning(_("solved items that were NaN/NA were replaced with 0.0"));
  }
  if (nlmOp.reducedTol){
    if (nlmOp.stickyTol){
      warning(_("tolerances (atol/rtol) were increased (after %d bad solves) for some difficult ODE solving during the optimization.\ncan control with foceiControl(stickyRecalcN=)\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"), nlmOp.stickyRecalcN);
    } else {
      warning(_("tolerances (atol/rtol) were temporarily increased for some difficult ODE solving during the optimization.\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"));
    }
  }
  return R_NilValue;
}

//[[Rcpp::export]]
RObject nlmGetParHist(bool p=true) {
  nlmOp.scale.save = 0;
  nlmOp.scale.print = 0;
  if (p) {
    scalePrintLine(min2(nlmOp.scale.npars, nlmOp.scale.printNcol));
  }
  return scaleParHisDf(&(nlmOp.scale));
}


//[[Rcpp::export]]
RObject nlmAdjustHessian(RObject Hin, arma::vec theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  arma::mat J(nlmOp.ntheta, nlmOp.ntheta);
  arma::mat H = as<arma::mat>(Hin);
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    J(i, i) =1.0/scaleAdjustGradScale(&(nlmOp.scale), 1.0, &theta[0], i);
  }
  H = J * H * J;
  RObject ret = wrap(H);
  ret.attr("dimnames") = Hin.attr("dimnames");
  return ret;
}

//[[Rcpp::export]]
RObject nlmAdjustCov(RObject CovIn, arma::vec theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  arma::mat J(nlmOp.ntheta, nlmOp.ntheta);
  arma::mat Cov = as<arma::mat>(CovIn);
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    J(i, i) = scaleAdjustGradScale(&(nlmOp.scale), 1.0, &theta[0], i);
  }
  Cov = J * Cov * J;
  RObject ret = wrap(Cov);
  ret.attr("dimnames") = CovIn.attr("dimnames");
  return ret;
}
