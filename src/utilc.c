#define STRICT_R_HEADER
#define USE_FC_LEN_T
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
#include <R_ext/BLAS.h>
#include <Rmath.h>
#include <rxode2.h>
#define _(String) (String)

#include "utilc.h"
#include "rxProtect.h"
#ifndef FCONE
#define FCONE
#endif

int _setSilentErr=0;

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

#ifdef _OPENMP
#include <omp.h>
#endif
void RSprintf(const char *format, ...) {
#ifdef _OPENMP
  // R's print API is single-threaded; only emit from the master thread.
  if (omp_get_thread_num() != 0) return;
#endif
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

// Shared arg-extraction + per-obs apply for the transform derivative wrappers
// (dy'/dlambda, d2y'/dlambda2).  Same (x,lambda,yj,low,hi)
// contract as _nlmixr2est_powerD; returns a length-`len` vector of fn() per obs.
static SEXP _nlmixr2estPowerApply(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS,
                                  double (*fn)(double, double, int, double, double)) {
  int t = TYPEOF(xS);
  int len = Rf_length(xS);
  if (t != REALSXP) Rf_errorcall(R_NilValue, _("'x' must be a real number"));
  double *x = REAL(xS);
  if (len != Rf_length(lambdaS) || len != Rf_length(yjS) ||
      len != Rf_length(lowS) || len != Rf_length(hiS))
    Rf_errorcall(R_NilValue, _("all arguments must be the same length"));
  if (TYPEOF(lambdaS) != REALSXP) Rf_errorcall(R_NilValue, _("'lambda' must be a real number"));
  double *lambda = REAL(lambdaS);
  if (TYPEOF(yjS) != INTSXP) Rf_errorcall(R_NilValue, _("'yj' must be an integer number"));
  int *yj = INTEGER(yjS);
  if (TYPEOF(hiS) != REALSXP) Rf_errorcall(R_NilValue, _("'hi' must be a real number"));
  double *hi = REAL(hiS);
  if (TYPEOF(lowS) != REALSXP) Rf_errorcall(R_NilValue, _("'low' must be a real number"));
  double *low = REAL(lowS);
  rxProtectGuard;
  SEXP retS = rxP(Rf_allocVector(REALSXP, len));
  double *ret = REAL(retS);
  for (int i = len; i--;) ret[i] = fn(x[i], lambda[i], yj[i], low[i], hi[i]);
  rxUPAll();
  return retS;
}

// dy'/dlambda (transform value lambda-derivative), per observation
SEXP _nlmixr2est_powerDLambda(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS) {
  return _nlmixr2estPowerApply(xS, lambdaS, yjS, lowS, hiS, _powerDLambda);
}
// d2y'/dlambda2, per observation
SEXP _nlmixr2est_powerDLambda2(SEXP xS, SEXP lambdaS, SEXP yjS, SEXP lowS, SEXP hiS) {
  return _nlmixr2estPowerApply(xS, lambdaS, yjS, lowS, hiS, _powerDLambda2);
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
  SEXP cls = rxP(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
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
  SEXP cls = rxP(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(cls, 0, Rf_mkChar("data.frame"));
  Rf_setAttrib(ret, R_ClassSymbol, cls);
  rxUPAll();
  return ret;
}

// sum(x * y) as R computes it: in order, in long double
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("fp-contract=off")))
#endif
static double nmRsum(const double *x, const double *y, int n) {
  long double t = 0.0;
  for (int i = 0; i < n; ++i) {
    double p = x[i] * y[i];
    t += p;
  }
  return (double)t;
}

// Damped-BFGS update of the n x n column-major H from the secant pair (s, y),
// computed as R's .trustOuterBfgs() did, so the two agree bitwise: H %*% s
// through R's BLAS, sum() in long double and no fused multiply-adds.  Hs and r
// are n-long work vectors.  Returns 1 when H was updated.
#if defined(__GNUC__) && !defined(__clang__)
__attribute__((optimize("fp-contract=off")))
#endif
int nmTrustBfgsUpdate(int n, double *H, const double *s, const double *y,
                      double *Hs, double *r) {
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif
  if (n < 1) return 0;
  // R's %*% takes another route through a non-finite operand, but then s'Hs
  // is not finite either and the update is skipped
  for (int i = 0; i < n; ++i) {
    if (!R_FINITE(s[i]) || !R_FINITE(y[i])) return 0;
  }
  for (int i = 0; i < n * n; ++i) {
    if (!R_FINITE(H[i])) return 0;
  }
  const char *tr = "N";
  double one = 1.0, zero = 0.0;
  int inc = 1;
  F77_CALL(dgemv)(tr, &n, &n, &one, H, &n, s, &inc, &zero, Hs, &inc FCONE);
  double sBs = nmRsum(s, Hs, n);
  if (!R_FINITE(sBs) || !(sBs > 0)) return 0;
  double sy = nmRsum(s, y, n);
  // Nocedal & Wright, Numerical Optimization 2nd ed, Procedure 18.2
  if (sy >= 0.2 * sBs) {
    for (int i = 0; i < n; ++i) r[i] = y[i];
  } else {
    double th = 0.8 * sBs / (sBs - sy);
    double th1 = 1 - th;
    for (int i = 0; i < n; ++i) {
      double a = th * y[i];
      double b = th1 * Hs[i];
      r[i] = a + b;
    }
  }
  double sr = nmRsum(s, r, n);
  double sNorm = sqrt(nmRsum(s, s, n));
  double rNorm = sqrt(nmRsum(r, r, n));
  // a reject-then-shrink step gives a near-zero denominator
  double lim = 1e-10 * sNorm;
  lim = lim * rNorm;
  if (!R_FINITE(sr) || !(sr > lim)) return 0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double a = Hs[i] * Hs[j];
      a = a / sBs;
      double b = r[i] * r[j];
      b = b / sr;
      double h = H[i + n * j] - a;
      H[i + n * j] = h + b;
    }
  }
  return 1;
}

SEXP _nlmixr2est_trustBfgsUpdate(SEXP bS, SEXP sS, SEXP yS) {
  rxProtectGuard;
  int n = Rf_length(sS);
  if (TYPEOF(bS) != REALSXP || TYPEOF(sS) != REALSXP || TYPEOF(yS) != REALSXP ||
      Rf_length(yS) != n || Rf_length(bS) != n * n) {
    rxUPAll();
    Rf_errorcall(R_NilValue, _("trustBfgsUpdate: bad arguments"));
  }
  SEXP ret = rxP(Rf_duplicate(bS));
  double *Hs = (double*)R_alloc(2 * (size_t)n + 1, sizeof(double));
  nmTrustBfgsUpdate(n, REAL(ret), REAL(sS), REAL(yS), Hs, Hs + n);
  rxUPAll();
  return ret;
}
