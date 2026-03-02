#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "npde.h"
#include "cwres.h"
#include "res.h"
#include "ires.h"
#include "utilc.h"
#include "shrink.h"

/* Internal C calls, should not be called outside of C code. */
typedef void (*S_fp) (double *, double *);
extern void nelder_fn(S_fp func, int n, double *start, double *step,
                      int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
                      int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
                      int *iprint);

/* .Call calls */
extern SEXP neldermead_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
/* extern SEXP n1qn1_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); */
extern SEXP _nlmixr2est_llik_binomial_c(SEXP, SEXP, SEXP);

extern SEXP _nlmixr2est_llik_poisson(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_normal(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_betabinomial(SEXP, SEXP, SEXP);
extern SEXP _nlmixr2est_llik_student_t(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_beta(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_neg_binomial(SEXP, SEXP);

// FOCEi
extern SEXP _nlmixr2est_nlmixr2Parameters(SEXP, SEXP);

SEXP _nlmixr2est_foceiInnerLp(SEXP, SEXP);
SEXP _nlmixr2est_likInner(SEXP, SEXP);
SEXP _nlmixr2est_cholSE_(SEXP, SEXP);
SEXP _nlmixr2est_foceiLik(SEXP);
SEXP _nlmixr2est_foceiOfv(SEXP);
SEXP _nlmixr2est_foceiLik(SEXP);
SEXP _nlmixr2est_foceiOfv(SEXP);
SEXP _nlmixr2est_foceiNumericGrad(SEXP);

SEXP _nlmixr2est_foceiSetup_(SEXP, SEXP, SEXP, SEXP, SEXP,
                             SEXP, SEXP, SEXP, SEXP, SEXP,
                             SEXP);

SEXP _nlmixr2est_foceiOuterF(SEXP);
SEXP _nlmixr2est_foceiOuterG(SEXP);
SEXP _nlmixr2est_foceiOuter(SEXP);
SEXP _nlmixr2est_sqrtm(SEXP);
SEXP _nlmixr2est_foceiCalcCov(SEXP);
SEXP _nlmixr2est_foceiFitCpp_(SEXP);
SEXP _nlmixr2est_boxCox_(SEXP, SEXP, SEXP);
SEXP _nlmixr2est_iBoxCox_(SEXP, SEXP, SEXP);
SEXP _nlmixr2est_freeFocei(void);
SEXP _nlmixr2est_nlmixr2Gill83_(SEXP, SEXP, SEXP, SEXP, SEXP,
                                SEXP, SEXP, SEXP, SEXP);

SEXP _nlmixr2est_rxode2version4(SEXP);
SEXP _nlmixr2est_rxode2stateOde(SEXP);

SEXP _nlmixr2est_nlmixr2Grad_(SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2Eval_(SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2ParHist_(SEXP);
SEXP _nlmixr2est_nlmixr2Hess_(SEXP, SEXP, SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2Unscaled_(SEXP, SEXP);

SEXP _nlmixr2est_saem_fit(SEXP);
SEXP _nlmixr2est_saem_do_pred(SEXP, SEXP, SEXP);

SEXP _nlmixr2est_augPredTrans(SEXP, SEXP, SEXP, SEXP, SEXP,
                              SEXP);


SEXP _nlmixr2est_uninformativeEta(SEXP);

static const R_CMethodDef CEntries[] = {
  {NULL, NULL, 0}
};

SEXP _nlmixr2est_powerD(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _saemResidF(SEXP v);

SEXP _nlmixr2est_nlmixrExpandFdParNlme_(SEXP, SEXP);

//SEXP _nlmixr2est_nmNearPD_()
SEXP _nlmixr2est_nmNearPD_(SEXP, SEXP, SEXP, SEXP, SEXP,
                           SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _nlmixr2est_filterNormalLikeAndDoses(SEXP, SEXP, SEXP);

SEXP _nlmixr2est_rxode2hasLlik(void);

SEXP _nlmixr2est_RcppExport_registerCCallable(void);

SEXP _nlmixr2est_nlmSetup(SEXP);
SEXP _nlmixr2est_nlmFree(void);
SEXP _nlmixr2est_nlmSolveGradHess(SEXP);
SEXP _nlmixr2est_nlmSolveGradR(SEXP);
SEXP _nlmixr2est_nlmSolveR(SEXP);
SEXP _nlmixr2est_nlmSolveSwitch(SEXP);
SEXP _nlmixr2est_optimFunC(SEXP, SEXP);
SEXP _nlmixr2est_nlminbFunC(SEXP, SEXP);
SEXP _nlmixr2est_nlmWarnings(void);
SEXP _nlmixr2est_nlmCensInfo(void);

SEXP _nlmixr2est_nlmScalePar(SEXP);
SEXP _nlmixr2est_nlmUnscalePar(SEXP);

SEXP _nlmixr2est_solveGradNls(SEXP, SEXP);
SEXP _nlmixr2est_nlmGetScaleC(SEXP, SEXP);

SEXP _nlmixr2est_nlmAdjustHessian(SEXP, SEXP);
SEXP _nlmixr2est_nlmAdjustCov(SEXP, SEXP);
SEXP _nlmixr2est_nlmSetScaleC(SEXP);
SEXP _nlmixr2est_nlmPrintHeader(void);
SEXP _nlmixr2est_nlmGetParHist(SEXP);

SEXP _nlmixr2est_iniLotriPtr(SEXP ptr);

SEXP _nlmixr2est_iniRxodePtrs(SEXP ptr);

SEXP _nlmixr2est_iniN1qn1cPtrs(SEXP ptr);

SEXP _nlmixr2est_iniLbfgsb3ptr(SEXP ptr);

SEXP _rxode2version4(SEXP);
SEXP _nlmixr2est_rxode2stateOde(SEXP);
SEXP _rxode2rxFixRes(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_rxode2rxFixRes", (DL_FUNC) &_rxode2rxFixRes, 2},
  {"_rxode2version4", (DL_FUNC) &_rxode2version4, 1},
  {"_nlmixr2est_rxode2stateOde", (DL_FUNC) &_nlmixr2est_rxode2stateOde, 1},
  {"_nlmixr2est_iniLbfgsb3ptr", (DL_FUNC) &_nlmixr2est_iniLbfgsb3ptr, 1},
  {"_nlmixr2est_iniN1qn1cPtrs", (DL_FUNC) &_nlmixr2est_iniN1qn1cPtrs, 1},
  {"_nlmixr2est_iniRxodePtrs", (DL_FUNC) &_nlmixr2est_iniRxodePtrs, 1},
  {"_nlmixr2est_iniLotriPtr", (DL_FUNC) &_nlmixr2est_iniLotriPtr, 1},
  {"_nlmixr2est_uninformativeEta", (DL_FUNC) &_nlmixr2est_uninformativeEta, 1},
  {"_nlmixr2est_nlmGetParHist", (DL_FUNC) &_nlmixr2est_nlmGetParHist, 1},
  {"_nlmixr2est_nlmPrintHeader", (DL_FUNC) &_nlmixr2est_nlmPrintHeader, 0},
  {"_nlmixr2est_nlmSetScaleC", (DL_FUNC) &_nlmixr2est_nlmSetScaleC, 1},
  {"_nlmixr2est_nlmAdjustCov", (DL_FUNC) &_nlmixr2est_nlmAdjustCov, 2},
  {"_nlmixr2est_nlmAdjustHessian", (DL_FUNC) &_nlmixr2est_nlmAdjustHessian, 2},
  {"_nlmixr2est_nlmGetScaleC", (DL_FUNC) &_nlmixr2est_nlmGetScaleC, 2},
  {"_nlmixr2est_nlmScalePar", (DL_FUNC) &_nlmixr2est_nlmScalePar, 1},
  {"_nlmixr2est_nlmUnscalePar", (DL_FUNC) &_nlmixr2est_nlmUnscalePar, 1},
  {"_nlmixr2est_solveGradNls", (DL_FUNC) &_nlmixr2est_solveGradNls, 2},
  {"_nlmixr2est_nlminbFunC", (DL_FUNC) &_nlmixr2est_nlminbFunC, 2},
  {"_nlmixr2est_nlmWarnings", (DL_FUNC) &_nlmixr2est_nlmWarnings, 0},
  {"_nlmixr2est_nlmCensInfo", (DL_FUNC) &_nlmixr2est_nlmCensInfo, 0},
  {"_nlmixr2est_optimFunC", (DL_FUNC) &_nlmixr2est_optimFunC, 2},
  {"_nlmixr2est_nlmSolveSwitch", (DL_FUNC) &_nlmixr2est_nlmSolveSwitch, 1},
  {"_nlmixr2est_nlmSolveR", (DL_FUNC) &_nlmixr2est_nlmSolveR, 1},
  {"_nlmixr2est_nlmSetup", (DL_FUNC) &_nlmixr2est_nlmSetup, 1},
  {"_nlmixr2est_nlmSolveGradR", (DL_FUNC) &_nlmixr2est_nlmSolveGradR, 1},
  {"_nlmixr2est_nlmSolveGradHess", (DL_FUNC) &_nlmixr2est_nlmSolveGradHess, 1},
  {"_nlmixr2est_nlmFree", (DL_FUNC) &_nlmixr2est_nlmFree, 0},
  {"_nlmixr2est_RcppExport_registerCCallable", (DL_FUNC) &_nlmixr2est_RcppExport_registerCCallable, 0},
  {"_nlmixr2est_rxode2hasLlik", (DL_FUNC) &_nlmixr2est_rxode2hasLlik, 0},
  {"_nlmixr2est_freeFocei", (DL_FUNC) &_nlmixr2est_freeFocei, 0},
  {"_nlmixr2est_filterNormalLikeAndDoses", (DL_FUNC) &_nlmixr2est_filterNormalLikeAndDoses, 3},
  {"neldermead_wrap",      (DL_FUNC) &neldermead_wrap,      11},
  /* {"n1qn1_wrap",           (DL_FUNC) &n1qn1_wrap,           13}, */
  {"_nlmixr2est_nlmixr2Parameters", (DL_FUNC) &_nlmixr2est_nlmixr2Parameters, 2},
  // FOCEi
  {"_nlmixr2est_foceiInnerLp", (DL_FUNC) &_nlmixr2est_foceiInnerLp, 2},
  {"_nlmixr2est_cholSE_", (DL_FUNC) &_nlmixr2est_cholSE_, 2},
  {"_nlmixr2est_likInner", (DL_FUNC) &_nlmixr2est_likInner, 2},
  {"_nlmixr2est_foceiLik", (DL_FUNC) &_nlmixr2est_foceiLik, 1},
  {"_nlmixr2est_foceiOfv", (DL_FUNC) &_nlmixr2est_foceiOfv, 1},
  {"_nlmixr2est_foceiNumericGrad", (DL_FUNC) &_nlmixr2est_foceiNumericGrad, 1},
  {"_nlmixr2est_foceiSetup_", (DL_FUNC) &_nlmixr2est_foceiSetup_, 11},
  {"_nlmixr2est_foceiOuterF", (DL_FUNC) &_nlmixr2est_foceiOuterF, 1},
  {"_nlmixr2est_foceiOuterG", (DL_FUNC) &_nlmixr2est_foceiOuterG, 1},
  {"_nlmixr2est_foceiOuter", (DL_FUNC) &_nlmixr2est_foceiOuter, 1},
  {"_nlmixr2est_sqrtm", (DL_FUNC) &_nlmixr2est_sqrtm, 1},
  {"_nlmixr2est_foceiCalcCov", (DL_FUNC) &_nlmixr2est_foceiCalcCov, 1},
  {"_nlmixr2est_foceiFitCpp_", (DL_FUNC) &_nlmixr2est_foceiFitCpp_, 1},
  {"_nlmixr2est_boxCox_", (DL_FUNC) &_nlmixr2est_boxCox_, 3},
  {"_nlmixr2est_iBoxCox_", (DL_FUNC) &_nlmixr2est_iBoxCox_, 3},
  {"_nlmixr2est_nlmixr2Gill83_", (DL_FUNC) &_nlmixr2est_nlmixr2Gill83_, 9},
  {"_nlmixr2est_nlmixr2Grad_", (DL_FUNC) &_nlmixr2est_nlmixr2Grad_, 2},
  {"_nlmixr2est_nlmixr2Eval_", (DL_FUNC) &_nlmixr2est_nlmixr2Eval_, 2},
  {"_nlmixr2est_nlmixr2ParHist_", (DL_FUNC) &_nlmixr2est_nlmixr2ParHist_, 1},
  {"_nlmixr2est_nlmixr2Hess_", (DL_FUNC) &_nlmixr2est_nlmixr2Hess_, 4},
  {"_nlmixr2est_augPredTrans", (DL_FUNC) &_nlmixr2est_augPredTrans, 6},
  {"_nlmixr2est_nlmixr2Unscaled_", (DL_FUNC) &_nlmixr2est_nlmixr2Unscaled_, 2},
  {"_nlmixr2est_setSilentErr", (DL_FUNC) &_nlmixr2est_setSilentErr, 1},
  {"_nlmixr2est_saem_fit", (DL_FUNC) &_nlmixr2est_saem_fit, 1},
  {"_nlmixr2est_saem_do_pred", (DL_FUNC) &_nlmixr2est_saem_do_pred, 3},
  {"_nlmixr2est_powerD", (DL_FUNC) &_nlmixr2est_powerD, 5},
  {"_nlmixr2est_powerL", (DL_FUNC) &_nlmixr2est_powerL, 5},
  {"_saemResidF", (DL_FUNC) &_saemResidF, 1},
  {"_nlmixr2est_npdeCalc", (DL_FUNC) &_nlmixr2est_npdeCalc, 6},
  {"_nlmixr2est_cwresCalc",  (DL_FUNC) &_nlmixr2est_cwresCalc, 12},
  {"_nlmixr2est_resCalc",  (DL_FUNC) &_nlmixr2est_resCalc, 12},
  {"_nlmixr2est_iresCalc", (DL_FUNC) &_nlmixr2est_iresCalc, 10},
  {"_nlmixr2est_calcShrinkOnly", (DL_FUNC) &_nlmixr2est_calcShrinkOnly, 3},
  {"_nlmixr2est_popResFinal", (DL_FUNC) &_nlmixr2est_popResFinal, 1},
  {"_nlmixr2est_nlmixrExpandFdParNlme_", (DL_FUNC) &_nlmixr2est_nlmixrExpandFdParNlme_, 2},
  {"_nlmixr2est_nmNearPD_", (DL_FUNC) &_nlmixr2est_nmNearPD_, 10},
  {NULL, NULL, 0}
};

void R_init_nlmixr2est(DllInfo *dll)
{
  R_RegisterCCallable("nlmixr2est","nelder_fn", (DL_FUNC) &nelder_fn);
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll,FALSE);
}

void rxOptionsFreeFocei(void);
void R_unload_nlmixr2est(DllInfo *info){
  rxOptionsFreeFocei();
}
