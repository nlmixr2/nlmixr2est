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
  double *thetahf=NULL; // Shi step size
  double *thetahh=NULL;

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
};

nlmOptions nlmOp;

//[[Rcpp::export]]
RObject nlmFree() {
  if (nlmOp.thetaFD != NULL) R_Free(nlmOp.thetaFD);
  nlmOp.thetaFD = NULL;
  if (nlmOp.thetahf != NULL) R_Free(nlmOp.thetahf);
  nlmOp.thetahf = NULL;
  nlmOp.thetahh = NULL;
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

  nlmOp.thetaFD = R_Calloc(nlmOp.ntheta, int);
  IntegerVector needFD = as<IntegerVector>(e["needFD"]);

  nlmOp.needFD=false;
  for (int i = 0; i < nlmOp.ntheta; ++i) {
    nlmOp.thetaFD[i] = needFD[i];
    if (nlmOp.thetaFD[i]) {
      nlmOp.needFD=true;
    }
  }

  nlmOp.thetahf = R_Calloc(nlmOp.ntheta*2, double);
  nlmOp.thetahh = nlmOp.thetahf+nlmOp.ntheta;

  nlmOp.stickyRecalcN=as<int>(control["stickyRecalcN"]);
  nlmOp.stickyTol=0;
  nlmOp.stickyRecalcN2=0;
  nlmOp.stickyRecalcN1=0;
  nlmOp.reducedTol = 0;
  nlmOp.reducedTol2 = 0;
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
arma::vec calcGradCentral(arma::vec &grMH, arma::vec &f0, arma::vec &grPH,  double h);

// Solve prediction
arma::vec nlmSolveF(arma::vec &theta, int id) {
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  nlmSolvePred(id);
  rx_solving_options *op = rx->op;
  arma::vec ret(ind->n_all_times);
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
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
      rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
      ret(k) = ind->lhs[0];
      k++;
    }
    if (k >= ind->n_all_times - ind->ndoses - ind->nevid2) {
      // With moving doses this may be at the very end, so drop out
      // now if all the observations were accounted for
      break;
    }
  }
  return ret(span(0, k-1));
}

//[[Rcpp::export]]
double nlmSolveR(arma::vec &theta) {
  int id = 0;
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  nlmSolvePred(id);
  rx_solving_options *op = rx->op;
  double ret = 0.0;
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
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
      rxPred.calc_lhs(id, curT, getSolve(j), ind->lhs);
      ret += ind->lhs[0];
      k++;
    }
    if (k >= ind->n_all_times - ind->ndoses - ind->nevid2) {
      // With moving doses this may be at the very end, so drop out
      // now if all the observations were accounted for
      break;
    }
  }
  return ret;
}

// Solve gradient with finite difference thetas This is done for every
// observation so it could possibly be used in nls types of models
arma::mat nlmSolveGrad(arma::vec &theta, int id) {
  // first solve the nlm problem
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  nlmSolveNlm(id);
  rx_solving_options *op = rx->op;
  arma::mat ret(ind->n_all_times, nlmOp.ntheta+1);
  int kk, k=0;
  double curT;
  iniSubjectI(id, 1, ind, op, rx, rxInner.update_inis);
  for (int j = 0; j < ind->n_all_times; ++j) {
    ind->idx=j;
    kk = ind->ix[j];
    curT = getTimeF(kk, ind);
    if (isDose(ind->evid[kk])) {
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      continue;
    } else if (ind->evid[kk] == 0) {
      double *solveCur = getSolve(j);
      rxInner.calc_lhs(id, curT, getSolve(j), ind->lhs);
      for (int kk = 0; kk < op->nlhs; ++kk) {
        ret(k, kk) = ind->lhs[kk];
      }
      k++;
    }
    if (k >= ind->n_all_times - ind->ndoses - ind->nevid2) {
      // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
      break;
    }
  }
  ret = ret.rows(0, k-1);
  int nObs = k;
  if (nlmOp.needFD) {
    // Save solve; not needed can corrupt memory space with preds
    //int nsolve = (op->neq + op->nlin)*ind->n_all_times;
    //arma::vec solveSave(nsolve);
    // save and restore memory pointer
    // std::copy(ind->solve, ind->solve + nsolve, solveSave.memptr());
    arma::vec f0 = ret.col(0);
    arma::vec grTheta(nObs);
    arma::vec grPH(nObs);
    arma::vec grMH(nObs);

    arma::vec hTheta(nlmOp.ntheta);
    arma::vec curTheta = theta;
    for (int ii = 0; ii < nlmOp.ntheta; ++ii) {
      if (nlmOp.thetaFD[ii] == 0) continue;
      if (nlmOp.thetahf[ii] == 0.0) {
        double h = 0;
        switch(nlmOp.eventType) {
        case 2: // central
          nlmOp.thetahf[ii] = shi21Central(nlmSolveF, curTheta, h,
                                           f0, grTheta, id, ii,
                                           nlmOp.shiErr, // ef,
                                           1.5,//double rl = 1.5,
                                           4.5,//double ru = 4.5,
                                           3.0,//double nu = 8.0);
                                           nlmOp.shi21maxFD); // maxiter
          break;
        case 1: // forward
          nlmOp.thetahf[ii] = shi21Forward(nlmSolveF, curTheta, h,
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
      grPH = nlmSolveF(hTheta, id);
      bool useForward = false;
      if (nlmOp.eventType == 1) {
        // if this isn't true try backward
        if (grPH.is_finite()) {
          useForward = true;
          ret.col(ii) = calcGradForward(f0, grPH,  nlmOp.thetahf[ii]);
          continue;
        }
      }
      if (!useForward) {
        // stencil or central
        hTheta = curTheta;
        hTheta[ii] -= nlmOp.thetahf[ii];
        grMH = nlmSolveF(hTheta, id);
        // central
        ret.col(ii) = calcGradCentral(grMH, f0, grPH,  nlmOp.thetahf[ii]);
      }
    }
    // restore save (may not be needed)
    //std::copy(solveSave.begin(), solveSave.end(), ind->solve);
  }
  return ret;
}

//[[Rcpp::export]]
RObject nlmSolveGradR(arma::vec &theta) {
  int ntheta = theta.size();
  arma::mat ret0 = nlmSolveGrad(theta, 0);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  NumericVector ret(1);
  NumericVector grad(ntheta);
  ret[0] = cs[0];
  ret.attr("gradient") = wrap(cs(span(1, ntheta)));
  return ret;
}

arma::vec nlmSolveGrad1(arma::vec &theta, int id) {
  int ntheta = theta.size();
  arma::mat ret0 = nlmSolveGrad(theta, 0);
  ret0 = ret0.cols(1, ntheta);
  return (arma::sum(ret0, 0)).t();
}
//[[Rcpp::export]]
NumericVector nlmSolveGradOnly(arma::vec &theta) {
  arma::vec gr=nlmSolveGrad1(theta, 0);
  return wrap(gr(span(0, nlmOp.ntheta-1)));
}

//[[Rcpp::export]]
RObject nlmSolveGradHess(arma::vec &theta) {
  int id = 0;
  int ntheta = theta.size();
  arma::mat ret0 = nlmSolveGrad(theta, 0);
  arma::vec cs = (arma::sum(ret0, 0)).t();
  double ll = cs[0];
  arma::vec gr0 = cs(span(1, ntheta));
  mat H(nlmOp.ntheta, nlmOp.ntheta);
  H.zeros();
  arma::vec grPH(nlmOp.ntheta);
  arma::vec grMH(nlmOp.ntheta);
  double h;
  for (int k = nlmOp.ntheta; k--;) {
    h = nlmOp.thetahh[k];
    if (nlmOp.optimHessType == 1 && h <= 0) {
      arma::vec t = theta;
      nlmOp.thetahh[k] = shi21Forward(nlmSolveGrad1, theta, h,
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
      nlmOp.thetahh[k] = shi21Central(nlmSolveGrad1, t, h,
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
    grMH = nlmSolveGrad1(theta, id);
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
  ret.attr("gradient") = wrap(gr0(span(0, ntheta-1)));
  ret.attr("hessian") = wrap(H);
  return ret;
}

//[[Rcpp::export]]
RObject nlmSolveSwitch(arma::vec &theta) {
  switch(nlmOp.solveType) {
  case 1:
    return wrap(nlmSolveR(theta));
  case 2:
    return nlmSolveGradR(theta);
  default:
    return nlmSolveGradHess(theta);
  }
}
