#ifndef __INNER_H__
#define __INNER_H__
struct rxSolveF {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  t_dydt dydt = NULL;
  t_calc_jac calc_jac = NULL;
  t_calc_lhs calc_lhs = NULL;
  t_update_inis update_inis = NULL;
  t_dydt_lsoda_dum dydt_lsoda_dum = NULL;
  t_dydt_liblsoda dydt_liblsoda = NULL;
  t_jdum_lsoda jdum_lsoda = NULL;
  t_set_solve set_solve = NULL;
  t_get_solve get_solve = NULL;
  int global_jt = 2;
  int global_mf = 22;
  int global_debug = 0;
  // No neq here: the state count lives in the odeSwap registry (src/odeSwap.h),
  // which is what the pool decision reads.  The field that used to sit here was
  // written in four places and never read.
};

extern rxSolveF rxInner;
extern rxSolveF rxPred;
extern rxSolveF rxThetaSens;
extern rxSolveF rxHess2;
extern rxSolveF rxVaeOuter;
extern rxSolveF rxOuterNode;
extern rxSolveF rxOuterCov;
extern void rxUpdateFuns(SEXP trans, rxSolveF *inner);
extern void rxClearFuns(rxSolveF *inner);
extern rx_solve *rx;
#endif
