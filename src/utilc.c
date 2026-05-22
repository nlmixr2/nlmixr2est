#define STRICT_R_HEADER
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>   /* dj: import intptr_t */
#include <errno.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>
#include <rxode2.h>
#define _(String) (String)

#include "utilc.h"
#include "rxProtect.h"

int _setSilentErr=0;
extern void setSilentErr(int silent){
  _setSilentErr = silent;
}

SEXP _nlmixr2est_setSilentErr(SEXP in) {
  rxProtectGuard;
  SEXP ret = rxP(Rf_allocVector(LGLSXP, 1));
  int t = TYPEOF(in);
  if (Rf_length(in) > 0) {
    if (t == INTSXP) {
      if (INTEGER(in)[0] > 0) {
        _setSilentErr = 1;
        INTEGER(ret)[0] = 1;
        rxUPAll();
        return ret;
      } else {
        _setSilentErr = 0;
        INTEGER(ret)[0] = 0;
        rxUPAll();
        return ret;
      }
    } else if (t == LGLSXP) {
      if (INTEGER(in)[0] > 0) {
        _setSilentErr = 1;
        INTEGER(ret)[0] = 1;
        rxUPAll();
        return ret;
      } else {
        _setSilentErr = 0;
        INTEGER(ret)[0] = 0;
        rxUPAll();
        return ret;
      }
    } else if (t == REALSXP) {
      if (REAL(in)[0] > 0) {
        _setSilentErr = 1;
        INTEGER(ret)[0] = 1;
        rxUPAll();
        return ret;
      } else {
        _setSilentErr = 0;
        INTEGER(ret)[0] = 0;
        rxUPAll();
        return ret;
        return R_NilValue;
      }
    }
  }
  _setSilentErr = 0;
  INTEGER(ret)[0] = 0;
  rxUPAll();
  return ret;
}

void RSprintf(const char *format, ...) {
  if (_setSilentErr == 0) {
    va_list args;
    va_start(args, format);
    Rvprintf(format, args);
    va_end(args);
  }
}
// double x, double lambda, int yj, double low, double high
SEXP _nlmixr2est_powerD(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS) {
  int t = TYPEOF(xS);
  int len = Rf_length(xS);
  if (t != REALSXP) {
    Rf_errorcall(R_NilValue, _("'x' must be a real number"));
  }
  double *x =REAL(xS);
  if (len != Rf_length(lambdaS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(yjS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(lowS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(hiS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  t = TYPEOF(lambdaS);
  if (t != REALSXP) {
    Rf_errorcall(R_NilValue, _("'lambda' must be a real number"));
  }
  double *lambda = REAL(lambdaS);
  int *yj;
  t = TYPEOF(yjS);
  if (t == INTSXP) {
    yj = INTEGER(yjS);
  } else {
    Rf_errorcall(R_NilValue, _("'yj' must be an integer number"));
  }
  double *hi;
  t = TYPEOF(hiS);
  if (t == REALSXP) {
    hi = REAL(hiS);
  } else {
    Rf_errorcall(R_NilValue, _("'hi' must be a real number"));
  }
  t = TYPEOF(lowS);
  double *low;
  if (t == REALSXP) {
    low = REAL(lowS);
  } else {
    Rf_errorcall(R_NilValue, _("'low' must be a real number"));
  }
  rxProtectGuard;
  SEXP retS = rxP(Rf_allocVector(REALSXP, len));
  double *ret = REAL(retS);
  for (int i = len; i--;) {
    ret[i] = _powerD(x[i], lambda[i], yj[i], low[i], hi[i]);
  }
  rxUPAll();
  return retS;
}

// double x, double lambda, int yj, double low, double high
SEXP _nlmixr2est_powerL(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS) {
  int t = TYPEOF(xS);
  int len = Rf_length(xS);
  if (t != REALSXP) {
    Rf_errorcall(R_NilValue, _("'x' must be a real number"));
  }
  double *x =REAL(xS);
  if (len != Rf_length(lambdaS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(yjS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(lowS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  if (len != Rf_length(hiS)) {
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  }
  t = TYPEOF(lambdaS);
  if (t != REALSXP) {
    Rf_errorcall(R_NilValue, _("'lambda' must be a real number"));
  }
  double *lambda = REAL(lambdaS);
  int *yj;
  t = TYPEOF(yjS);
  if (t == INTSXP) {
    yj = INTEGER(yjS);
  } else {
    Rf_errorcall(R_NilValue, _("'yj' must be an integer number"));
  }
  double *hi;
  t = TYPEOF(hiS);
  if (t == REALSXP) {
    hi = REAL(hiS);
  } else {
    Rf_errorcall(R_NilValue, _("'hi' must be a real number"));
  }
  t = TYPEOF(lowS);
  double *low;
  if (t == REALSXP) {
    low = REAL(lowS);
  } else {
    Rf_errorcall(R_NilValue, _("'low' must be a real number"));
  }
  rxProtectGuard;
  SEXP retS = rxP(Rf_allocVector(REALSXP, 1));
  double *ret = REAL(retS);
  ret[0] = 0;
  for (int i = len; i--;) {
    ret[0] += _powerL(x[i], lambda[i], yj[i], low[i], hi[i]);
  }
  rxUPAll();
  return retS;
}

SEXP getDfSubsetVars(SEXP ipred, SEXP lhs) {
  int type = TYPEOF(lhs);
  if (type != STRSXP) return R_NilValue;
  if (Rf_length(lhs) == 0) return R_NilValue;
  rxProtectGuard;
  SEXP ipredNames = rxP(Rf_getAttrib(ipred, R_NamesSymbol));
  int *keepVals = R_Calloc((size_t)Rf_length(ipredNames), int);
  R_xlen_t k = 0;
  for (R_xlen_t i = 0; i < Rf_length(ipredNames); ++i) {
    for (R_xlen_t j = 0; j < Rf_length(lhs); ++j) {
      if (!strcmp(CHAR(STRING_ELT(ipredNames, i)), CHAR(STRING_ELT(lhs, j)))) {
        keepVals[k++] = i;
        break;
      }
    }
  }
  if (k == 0) {
    R_Free(keepVals);
    rxUPAll();
    return R_NilValue;
  }
  SEXP ret = rxP(Rf_allocVector(VECSXP, k));
  SEXP nm = rxP(Rf_allocVector(STRSXP, k));
  for (R_xlen_t i = 0; i < k; ++i) {
    SET_VECTOR_ELT(ret,i,VECTOR_ELT(ipred, keepVals[i]));
    SET_STRING_ELT(nm,i,STRING_ELT(ipredNames, keepVals[i]));
  }
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  SEXP cls = rxP(allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, mkChar("data.frame"));
  Rf_setAttrib(ret, R_ClassSymbol, cls);
  SEXP rn = rxP(Rf_allocVector(INTSXP, 2));
  int *rni =INTEGER(rn);
  rni[0] = NA_INTEGER;
  rni[1] = -Rf_length(VECTOR_ELT(ret,0));
  Rf_setAttrib(ret, R_RowNamesSymbol, rn);
  R_Free(keepVals);
  rxUPAll();
  return ret;
}


SEXP dfCbindList(SEXP lst) {
  int type = TYPEOF(lst);
  if (type != VECSXP) return R_NilValue;
  int totN=0;
  rxProtectGuard;
  SEXP curS;
  SEXP curN;
  SEXP curElt;
  for (int i = 0; i < Rf_length(lst); ++i) {
    curS = rxP(VECTOR_ELT(lst, i));
    if (TYPEOF(curS) == VECSXP) {
      totN += Rf_length(curS);
    }
  }
  if (totN == 0) {
    rxUPAll();
    return R_NilValue;
  }
  SEXP ret = rxP(Rf_allocVector(VECSXP, totN));
  SEXP nm  = rxP(Rf_allocVector(STRSXP, totN));
  int k=0;
  for (int i = 0; i < Rf_length(lst); ++i) {
    curS = rxP(VECTOR_ELT(lst, i));
    if (TYPEOF(curS) == VECSXP) {
      curN = rxP(Rf_getAttrib(curS, R_NamesSymbol));
      for (int j = 0; j < Rf_length(curN); ++j) {
        curElt = VECTOR_ELT(curS, j);
        Rf_setAttrib(curElt, R_DimSymbol, R_NilValue);
        SET_VECTOR_ELT(ret, k, curElt);
        SET_STRING_ELT(nm, k++, STRING_ELT(curN, j));
      }
    }
  }
  Rf_setAttrib(ret, R_NamesSymbol, nm);
  SEXP rn = rxP(Rf_allocVector(INTSXP, 2));
  int *rni = INTEGER(rn);
  rni[0] = NA_INTEGER;
  rni[1] = -Rf_length(VECTOR_ELT(ret, 0));
  Rf_setAttrib(ret, R_RowNamesSymbol, rn);
  SEXP cls = rxP(allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, mkChar("data.frame"));
  Rf_setAttrib(ret, R_ClassSymbol, cls);
  rxUPAll();
  return ret;
}
