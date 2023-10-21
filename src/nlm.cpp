// [[Rcpp::plugins(openmp)]]
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

#define nlmOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
#define predOde(id) ind_solve(rx, id, rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)

extern void doAssignFn(void);
extern rxSolveF rxInner;
extern rxSolveF rxPred;
extern void rxUpdateFuns(SEXP trans, rxSolveF *inner);
extern void rxClearFuns(rxSolveF *inner);
extern ind_solve_t ind_solve;
extern getRxSolve_t getRx;
extern rx_solve *rx;
extern getTime_t getTimeF;
extern iniSubjectI_t iniSubjectI;


struct nlmOptions {
  int ntheta=0;
  int *thetaFD=NULL; // theta needs finite difference?
  int *nobs = NULL;
  int *idS  = NULL;
  int *idF  = NULL;
  int nobsTot = 0;
  double *thetahf=NULL; // Shi step size
  double *thetahh=NULL;
  double *thetaSave = NULL;
  double *valSave   = NULL;
  double *grSave    = NULL;

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
  int naZero;
  int naGrad;
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
  if (nlmOp.thetahf != NULL) R_Free(nlmOp.thetahf);
  nlmOp.thetahf = NULL;
  nlmOp.thetahh = NULL;
  nlmOp.loaded = false;
  return R_NilValue;
}


//[[Rcpp::export]]
RObject nlmSetup(Environment e) {
  doAssignFn();
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
    nlmOp.solveType = 1;
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
                   e["param"],//const RObject &params =
                   e["data"],//const RObject &events =
                   R_NilValue, // inits
                   1);//const int setupOnly = 0
  rx = getRx();

  nlmOp.thetaFD = R_Calloc(nlmOp.ntheta + rx->nsub*3, int);
  nlmOp.nobs = nlmOp.thetaFD + nlmOp.ntheta;
  nlmOp.idS = nlmOp.nobs + rx->nsub;
  nlmOp.idF = nlmOp.idS + rx->nsub;

  IntegerVector needFD = as<IntegerVector>(e["needFD"]);

  // nlmOp.ntheta nlmOp.ntheta+1
  nlmOp.thetahf = R_Calloc(nlmOp.ntheta*(rx->nsub+3) + 1, double);
  nlmOp.thetahh = nlmOp.thetahf   + nlmOp.ntheta*rx->nsub; // [nlmOp.ntheta]
  nlmOp.thetaSave = nlmOp.thetahh + nlmOp.ntheta; // [nlmOp.ntheta]
  nlmOp.valSave = nlmOp.thetaSave + nlmOp.ntheta; // [1]
  nlmOp.grSave = nlmOp.valSave + 1; // nlmOp.ntheta
  std::fill_n(nlmOp.thetaSave, nlmOp.ntheta, R_PosInf); // not likely to be equal

  nlmOp.needFD=false;
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    nlmOp.thetaFD[i] = needFD[i];
    if (nlmOp.thetaFD[i]) {
      nlmOp.needFD=true;
    }
  }
  // now calculate nobs per id
  nlmOp.nobsTot = 0;
  for (int id = 0; id < rx->nsub; ++id) {
    rx_solving_options_ind *ind = &(rx->subjects[id]);
    int no = 0;
    for (int j = 0; j < ind->n_all_times; ++j) {
      if (ind->evid[j] == 0) {
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
  // NumericVector nobsa(rx->nsub);
  // // save and restore memory pointer
  // // std::copy(ind->solve, ind->solve + nsolve, solveSave.memptr());
  // std::copy(nlmOp.nobs, nlmOp.nobs + rx->nsub, nobsa.begin());
  // IntegerVector idSa(rx->nsub);
  // std::copy(nlmOp.idS, nlmOp.idS + rx->nsub, idSa.begin());
  // IntegerVector idFa(rx->nsub);
  // std::copy(nlmOp.idF, nlmOp.idF + rx->nsub, idFa.begin());
  // List df = List::create(_["nobs"]=nobsa,
  //                        _["idS"]=idSa,
  //                        _["idF"]=idFa);
  // df.attr("class") = "data.frame";
  // df.attr("row.names") = IntegerVector::create(NA_INTEGER, -rx->nsub);
  // print(wrap(df));
  nlmOp.loaded = true;
  return R_NilValue;
}

void nlmSolveNlm(int id) {
  rx_solving_options *op = rx->op;
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  nlmOde(id);
  int j=0;
  while (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN &&
         op->badSolve && j < nlmOp.maxOdeRecalc) {
    nlmOp.stickyRecalcN2++;
    nlmOp.reducedTol  = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(nlmOp.odeRecalcFactor);
    ind->solved = -1;
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

void nlmSolvePred(int id) {
  rx_solving_options *op = rx->op;
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  predOde(id);
  int j=0;
  while (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN &&
         op->badSolve && j < nlmOp.maxOdeRecalc) {
    nlmOp.stickyRecalcN2++;
    nlmOp.reducedTol2 = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(nlmOp.odeRecalcFactor);
    ind->solved = -1;
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

// Solve prediction
void nlmSolveFid(double *retD, int nobs, arma::vec &theta, int id) {
  arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  rx_solving_options *op = rx->op;
  iniSubjectI(id, 1, ind, op, rx, rxPred.update_inis);
  nlmSolvePred(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < ind->n_all_times; ++j) {
    ind->idx=j;
    kk = ind->ix[j];
    curT = getTimeF(kk, ind);
    if (isDose(ind->evid[kk])) {
      rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
      continue;
    } else if (ind->evid[kk] == 0) {
      double *d = getSolve(j);
      rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
      if (ISNA(ind->lhs[0])) {
        nlmOp.naZero=1;
        ind->lhs[0] = 0.0;
      }
      ret(k) = ind->lhs[0];
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
  rx_solving_options *op = rx->op;
  int cores = op->cores;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(cores)
// #endif
  for (int id = 0; id < rx->nsub; ++id) {
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
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  rx_solving_options *op = rx->op;
  arma::mat ret(nlmOp.nobs[id], nlmOp.ntheta+1);
  int kk, k=0;
  double curT;
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
  nlmSolveNlm(id);
  for (int j = 0; j < ind->n_all_times; ++j) {
    ind->idx=j;
    kk = ind->ix[j];
    curT = getTimeF(kk, ind);
    if (isDose(ind->evid[kk])) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      continue;
    } else if (ind->evid[kk] == 0) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      for (int kk = 0; kk < op->nlhs; ++kk) {
        if (ISNA(ind->lhs[kk])) {
          ret(k, kk) = 0.0;
          nlmOp.naZero=1;
          continue;
        }
        ret(k, kk) = ind->lhs[kk];
      }
      k++;
    }
    if (k >= ind->n_all_times - ind->ndoses - ind->nevid2) {
      // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
      break;
    }
  }
  //if (nlmOp.needFD) {
    // Save solve; not needed can corrupt memory space with preds
    //int nsolve = (op->neq + op->nlin)*ind->n_all_times;
    //arma::vec solveSave(nsolve);
    // save and restore memory pointer
    // std::copy(ind->solve, ind->solve + nsolve, solveSave.memptr());
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
    //std::copy(solveSave.begin(), solveSave.end(), ind->solve);
    //  }
  return ret;
}

arma::mat nlmSolveGrad(arma::vec &theta) {
  arma::mat ret(nlmOp.nobsTot, nlmOp.ntheta+1);
  rx_solving_options *op = rx->op;
  int cores = op->cores;
// #ifdef _OPENMP
// #pragma omp parallel for num_threads(cores)
// #endif
  for (int id = 0; id < rx->nsub; ++id) {
    ret.rows(nlmOp.idS[id], nlmOp.idF[id]) = nlmSolveGradId(theta, id);
  }
  return ret;
}

//[[Rcpp::export]]
RObject nlmSolveGradR(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.solveType == 1) stop("incorrect solve type");
  int ntheta = theta.size();
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  NumericVector ret(1);
  NumericVector grad(ntheta);
  ret[0] = cs[0];
  ret.attr("gradient") = wrap(cs(span(1, ntheta)));
  return ret;
}

arma::vec nlmSolveGrad1(arma::vec &theta, int id) {
  arma::mat ret0 = nlmSolveGrad(theta);
  ret0 = ret0.cols(1, nlmOp.ntheta);
  return (arma::sum(ret0, 0)).t();
}

//[[Rcpp::export]]
RObject nlmSolveGradHess(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  if (nlmOp.solveType == 1) stop("incorrect solve type");
  int id = 0; // dummy id
  arma::mat ret0 = nlmSolveGrad(theta);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  double ll = cs[0];
  arma::vec gr0 = cs(span(1, nlmOp.ntheta));
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
  NumericVector ret(1);
  ret[0] = ll;
  ret.attr("gradient") = wrap(gr0(span(0, nlmOp.ntheta-1)));
  ret.attr("hessian") = wrap(H);
  return ret;
}

//[[Rcpp::export]]
RObject nlmSolveSwitch(arma::vec &theta) {
  if (!nlmOp.loaded) stop("'nlm' problem not loaded");
  switch(nlmOp.solveType) {
  case 1:
    return wrap(nlmSolveR(theta));
  case 2:
    return nlmSolveGradR(theta);
  case 3:
    return nlmSolveGradHess(theta);
  }
  return R_NilValue;
}

//[[Rcpp::export]]
NumericVector optimFunC(arma::vec &theta, bool grad=false) {
  if (!nlmOp.loaded) stop("'optim' problem not loaded");
  if (nlmOp.solveType == 1) {
    if (grad) stop(_("incorrect solve type"));
    NumericVector ret(1);
    ret[0] = nlmSolveR(theta);
    return ret;
  }
  bool isSame = true;
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    if (nlmOp.thetaSave[i] != theta[i]) {
      isSame = false;
      break;
    }
  }
  if (isSame) {
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
  arma::vec thetaSave(nlmOp.thetaSave, nlmOp.ntheta, false, true);
  thetaSave = theta;
  if (grad) {
    NumericVector ret(nlmOp.ntheta);
    std::copy(nlmOp.grSave, nlmOp.grSave + nlmOp.ntheta, ret.begin());
    return ret;
  }
  NumericVector ret(1);
  ret[0] = nlmOp.valSave[0];
  return ret;
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
