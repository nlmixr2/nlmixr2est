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

#define innerOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
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
  // calc/get
  bool needFD=false;
  double shiErr;

};

nlmOptions nlmOp;

void nlmFree(void) {
  if (nlmOp.thetaFD != NULL) R_Free(nlmOp.thetaFD);
  nlmOp.thetaFD = NULL;
  if (nlmOp.thetahf != NULL) R_Free(nlmOp.thetahf);
  nlmOp.thetahf = NULL;
}

RObject nlmSetup(Environment e) {
  doAssignFn();
  List control = e["control"];
  List model = e["thetaGrad"];

  List mv = rxode2::rxModelVars_(model);
  rxUpdateFuns(as<SEXP>(mv["trans"]), &rxPred);

  List pred = e["predOnly"];
  List mvp = rxode2::rxModelVars_(pred);
  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);

  List rxControl = as<List>(e["rxControl"]);

  NumericVector p = as<NumericVector>(e["param"]);
  nlmOp.ntheta = p.size();

  nlmOp.thetaFD = R_Calloc(nlmOp.ntheta, int);
  nlmOp.thetahf = R_Calloc(nlmOp.ntheta, double);

  nlmOp.stickyRecalcN=as<int>(control["stickyRecalcN"]);
  nlmOp.stickyTol=0;
  nlmOp.stickyRecalcN2=0;
  nlmOp.stickyRecalcN1=0;
  nlmOp.maxOdeRecalc = as<int>(control["maxOdeRecalc"]);
  nlmOp.odeRecalcFactor = as<double>(control["odeRecalcFactor"]);

  nlmOp.eventType = control["eventType"];
  nlmOp.shi21maxFD = control["shi21maxFD"];

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

void nlmSolveInner(int id) {
  rx_solving_options *op = rx->op;
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  innerOde(id);
  int j=0;
  while (nlmOp.stickyRecalcN2 <= nlmOp.stickyRecalcN &&
         op->badSolve && j < nlmOp.maxOdeRecalc) {
    nlmOp.stickyRecalcN2++;
    nlmOp.reducedTol  = 1;
    nlmOp.reducedTol2 = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(nlmOp.odeRecalcFactor);
    ind->solved = -1;
    innerOde(id);
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
    nlmOp.reducedTol  = 1;
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
    } else {
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

// Solve gradient with finite difference thetas
arma::mat nlmSolveGrad(arma::vec &theta, int id) {
  // first solve the inner problem
  rx_solving_options_ind *ind =  &(rx->subjects[id]);
  for (int i = nlmOp.ntheta; i--;) {
    ind->par_ptr[i]=theta[i];
  }
  nlmSolveInner(id);
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
    } else {
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
        case 3: // forward
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
      if (nlmOp.eventType == 3) {
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
