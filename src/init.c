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
extern SEXP slice_wrap(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

extern SEXP _nlmixr2est_llik_poisson(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_normal(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_betabinomial(SEXP, SEXP, SEXP);
extern SEXP _nlmixr2est_llik_student_t(SEXP, SEXP);
extern SEXP _nlmixr2est_llik_beta(SEXP, SEXP);
extern SEXP _nlmixr2est_lin_cmt_stan(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP);
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
			 SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _nlmixr2est_foceiOuterF(SEXP);
SEXP _nlmixr2est_foceiOuterG(SEXP);
SEXP _nlmixr2est_foceiOuter(SEXP);
SEXP _nlmixr2est_sqrtm(SEXP);
SEXP _nlmixr2est_foceiCalcCov(SEXP);
SEXP _nlmixr2est_foceiFitCpp_(SEXP);
SEXP _nlmixr2est_boxCox_(SEXP, SEXP, SEXP);
SEXP _nlmixr2est_iBoxCox_(SEXP, SEXP, SEXP);
SEXP _nlmixr2est_freeFocei();
SEXP _nlmixr2est_nlmixr2Gill83_(SEXP, SEXP, SEXP, SEXP, SEXP,
			   SEXP, SEXP, SEXP, SEXP);

SEXP _nlmixr2est_nlmixr2Grad_(SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2Eval_(SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2ParHist_(SEXP);
SEXP _nlmixr2est_nlmixr2Hess_(SEXP, SEXP, SEXP, SEXP);
SEXP _nlmixr2est_nlmixr2Unscaled_(SEXP, SEXP);

SEXP _nlmixr2est_saem_fit(SEXP);
SEXP _nlmixr2est_saem_do_pred(SEXP, SEXP, SEXP);

SEXP _nlmixr2est_augPredTrans(SEXP, SEXP, SEXP, SEXP, SEXP,
			  SEXP);

static const R_CMethodDef CEntries[] = {
    {NULL, NULL, 0}
};

SEXP _nlmixr2est_powerD(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP _saemResidF(SEXP v);

SEXP _nlmixr2est_nlmixrExpandFdParNlme_(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_nlmixr2est_freeFocei", (DL_FUNC) &_nlmixr2est_freeFocei, 0},
  {"neldermead_wrap",      (DL_FUNC) &neldermead_wrap,      11},
  /* {"n1qn1_wrap",           (DL_FUNC) &n1qn1_wrap,           13}, */
  {"_nlmixr2est_lin_cmt_stan",  (DL_FUNC) &_nlmixr2est_lin_cmt_stan,   9},
  {"_nlmixr2est_llik_binomial_c", (DL_FUNC) &_nlmixr2est_llik_binomial_c,  3},
  {"_nlmixr2est_llik_poisson",  (DL_FUNC) &_nlmixr2est_llik_poisson,   2},
  {"_nlmixr2est_llik_normal",   (DL_FUNC) &_nlmixr2est_llik_normal,    2},
  {"_nlmixr2est_llik_betabinomial", (DL_FUNC) &_nlmixr2est_llik_betabinomial, 3},
  {"_nlmixr2est_llik_student_t",  (DL_FUNC) &_nlmixr2est_llik_student_t, 2},
  {"_nlmixr2est_llik_beta",     (DL_FUNC) &_nlmixr2est_llik_beta, 2},
  {"_nlmixr2est_llik_neg_binomial", (DL_FUNC) &_nlmixr2est_llik_neg_binomial, 2},
  {"slice_wrap",           (DL_FUNC) &slice_wrap,            7},
  {"_nlmixr2est_nlmixr2Parameters", (DL_FUNC) &_nlmixr2est_nlmixr2Parameters, 2},
  // FOCEi
  {"_nlmixr2est_foceiInnerLp", (DL_FUNC) &_nlmixr2est_foceiInnerLp, 2},
  {"_nlmixr2est_cholSE_", (DL_FUNC) &_nlmixr2est_cholSE_, 2},
  {"_nlmixr2est_likInner", (DL_FUNC) &_nlmixr2est_likInner, 2},
  {"_nlmixr2est_foceiLik", (DL_FUNC) &_nlmixr2est_foceiLik, 1},
  {"_nlmixr2est_foceiOfv", (DL_FUNC) &_nlmixr2est_foceiOfv, 1},
  {"_nlmixr2est_foceiNumericGrad", (DL_FUNC) &_nlmixr2est_foceiNumericGrad, 1},
  {"_nlmixr2est_foceiSetup_", (DL_FUNC) &_nlmixr2est_foceiSetup_, 10},
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
  {NULL, NULL, 0}
};

void R_init_nlmixr2est(DllInfo *dll)
{
  R_RegisterCCallable("nlmixr2est","nelder_fn", (DL_FUNC) &nelder_fn);
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll,FALSE);
}

void rxOptionsFreeFocei();
void R_unload_nlmixr2est(DllInfo *info){
  rxOptionsFreeFocei();
}
