#define STRICT_R_HEADER
#include <stdio.h>
#include <stdarg.h>
#include <thread>
#include <chrono>
#include <vector>
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <rxode2ptr.h>
#include "nmMcmcRng.h"
#include "utilc.h"
#include "censEst.h"
#include "nearPD.h"
#include "inner.h"
#include "odeSwap.h"
#include "rxomp.h"
#include <memory>
#include "solveWarnHelper.h"
#include "truncNorm.h"

#define _(String) (String)

#define PHI(x) 0.5*(1.0+erf((x)/M_SQRT2))


#ifndef __SAEM_CLASS_RCPP_HPP__
#define __SAEM_CLASS_RCPP_HPP__
#define MAXENDPNT 40
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
using namespace std;
using namespace arma;
using namespace Rcpp;

// scale.h needs Rcpp:: types in scope (CharacterVector, RObject, warning, stop)
// -- must be included AFTER the `using namespace Rcpp;` above.
#include "scale.h"

typedef void (*fn_ptr) (double *, double *);

extern "C" void nelder_fn(fn_ptr func, int n, double *start, double *step,
			  int itmax, double ftol_rel, double rcoef, double ecoef, double ccoef,
			  int *iconv, int *it, int *nfcall, double *ynewlo, double *xmin,
			  int *iprint);

double *_saemYptr;
double *_saemFptr;
int _saemLen;
int _saemYj;
int _saemAddProp;
double _saemLambda;
double _saemLow;
double _saemHi;
fn_ptr _saemFn;
double *_saemStart;
double *_saemStep;
double _saemLambdaR;
double _saemPowR;
int _saemPropT=0;
bool _warnAtolRtol=false;
// ODE-freeze for the general-likelihood phi0 optimization: when true, user_function
// skips par_solve and reuses the states already in the solve buffers, recomputing
// only the per-observation log-likelihood (phi0 params are independent of the ODE).
static bool _saemFreezeOde=false;
static std::vector<double> _saemFtCache;
static std::vector<double> _saemYtrCache;
static std::vector<double> _saemFaCache;
static std::vector<double> _saemFaAdjustCache;
static double* _saemCacheYptr = nullptr;
static double* _saemCacheFptr = nullptr;
static int _saemCacheLen = -1;
static int _saemCacheYj = -1;
static int _saemCachePropT = -1;
static double _saemCacheLambda = std::numeric_limits<double>::quiet_NaN();
static double _saemCacheLow = std::numeric_limits<double>::quiet_NaN();
static double _saemCacheHi = std::numeric_limits<double>::quiet_NaN();
struct saem_state_t {
  int _saemIncreaseTol=0;
  int _saemIncreasedTol2=0;
  double _saemOdeRecalcFactor = 1.0;
  int _saemMaxOdeRecalc = 0;
  bool _saemIndTolRelax = true;
  mat _saemUE;
  int _saemMixest = 0;
};
static saem_state_t* current_saem_state = nullptr;

int _saemFixedIdx[4] = {0, 0, 0, 0};
double _saemFixedValue[4] = {0.0, 0.0, 0.0, 0.0};

// res_mod defines
#define rmAdd 1
#define rmProp 2
#define rmPow 3
#define rmAddProp 4
#define rmAddPow 5
#define rmAddLam 6
#define rmPropLam 7
#define rmPowLam 8
#define rmAddPropLam 9
#define rmAddPowLam 10

static inline double handleF(int powt, double &ft, double &f, bool trunc, bool adjustF) {
  double xmin = 1.0e-200, xmax=1e300;
  double fa = powt ? ft : f;
  if (adjustF && fa == 0.0) {
    fa = 1.0;
  }
  if (trunc){
    if (fa < xmin) fa = xmin;
    else if (fa > xmax) fa = xmax;
  }
  return fa;
}

// Per-observation combined-error SD for the E-step/simulation, matching the
// per-endpoint combined1/combined2 branch the M-step objective functions use
// (objC()/objD(): combined1 g = a + b*|f|^c, combined2 g = sqrt(a^2+b^2*f^(2c))).
// c (cres) defaults to 1 for every non-pow() residual model (rmAdd/rmProp/
// rmAddProp and their lambda siblings), so this reduces to the pre-#972
// formula there; only rmPow/rmAddPow/rmPowLam/rmAddPowLam actually estimate
// c != 1 (#972).
// Fills g in place (a scalar loop, no temporaries) so it can run inside the
// pre-allocated per-chain E-step scratch buffers (_scratch_g) without
// defeating their point.
static inline void saemFormG(vec &g, const vec &a, const vec &b, const vec &ft, const vec &c, const uvec &addPropVec) {
  const arma::uword n = ft.n_elem;
  for (arma::uword i = 0; i < n; ++i) {
    double fa = std::fabs(ft[i]);
    // c[i]==1 is every non-pow() endpoint (the overwhelming common case) --
    // skip pow() there so g is bit-identical to the pre-#972 formula instead
    // of merely numerically equal (SAEM's MCMC acceptance is sensitive
    // enough to a ULP-level g difference that it is not a no-op in practice).
    double fac = (c[i] == 1.0) ? fa : std::pow(fa, c[i]);
    if (addPropVec[i] == 1) {
      g[i] = a[i] + b[i]*fac;
    } else {
      g[i] = std::sqrt(a[i]*a[i] + b[i]*b[i]*fac*fac);
    }
  }
}

static inline void ensureSaemFixedTransformCache() {
  if (_saemCacheYptr == _saemYptr &&
      _saemCacheFptr == _saemFptr &&
      _saemCacheLen == _saemLen &&
      _saemCacheYj == _saemYj &&
      _saemCachePropT == _saemPropT &&
      _saemCacheLambda == _saemLambda &&
      _saemCacheLow == _saemLow &&
      _saemCacheHi == _saemHi) {
    return;
  }

  _saemFtCache.resize(_saemLen);
  _saemYtrCache.resize(_saemLen);
  _saemFaCache.resize(_saemLen);
  _saemFaAdjustCache.resize(_saemLen);
  for (int i = 0; i < _saemLen; ++i) {
    double ft = _powerD(_saemFptr[i], _saemLambda, _saemYj, _saemLow, _saemHi);
    double f = _saemFptr[i];
    _saemFtCache[i] = ft;
    _saemYtrCache[i] = _powerD(_saemYptr[i], _saemLambda, _saemYj, _saemLow, _saemHi);
    _saemFaCache[i] = handleF(_saemPropT, ft, f, false, false);
    _saemFaAdjustCache[i] = handleF(_saemPropT, ft, f, false, true);
  }
  _saemCacheYptr = _saemYptr;
  _saemCacheFptr = _saemFptr;
  _saemCacheLen = _saemLen;
  _saemCacheYj = _saemYj;
  _saemCachePropT = _saemPropT;
  _saemCacheLambda = _saemLambda;
  _saemCacheLow = _saemLow;
  _saemCacheHi = _saemHi;
}

#define toLambda(x) _powerDi(x, 1.0, 4, -_saemLambdaR, _saemLambdaR)
#define toLambdaEst(x) _powerD((x < -0.99*_saemLambdaR ? -0.99*_saemLambdaR : (x > 0.99*_saemLambdaR ? 0.99*_saemLambdaR : x)), 1.0, 4, -_saemLambdaR, _saemLambdaR)

#define toPow(x) _powerDi(x, 1.0, 4, -_saemPowR, _saemPowR)
#define toPowEst(x) _powerD((x < -0.99*_saemPowR ? -0.99*_saemPowR : (x > 0.99*_saemPowR ? 0.99*_saemPowR : x)), 1.0, 4, -_saemPowR, _saemPowR)

// add+prop
void obj(double *ab, double *fx) {
  ensureSaemFixedTransformCache();
  int i;
  double g, sum, cur, fa;
  double xmin = 1.0e-200, xmax=1e300, ft, ytr;
  double ab02;
  double ab12;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  ab02 = ab02*ab02;
  ab12 = ab12*ab12;
  for (i=0, sum=0; i<_saemLen; ++i) {
    ft = _saemFtCache[i];
    ytr = _saemYtrCache[i];
    fa = _saemFaCache[i];
    if (_saemAddProp == 1) {
      g = ab02 + ab12*fa;
    } else {
      g = sqrt(ab02*ab02 + ab12*ab12*fa*fa);
    }
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur*cur + 2*log(g);
  }
  *fx = sum;
}

// add + pow
void objC(double *ab, double *fx) {
  ensureSaemFixedTransformCache();
  int i;
  double g, sum, cur, ft, ytr, fa=1.0;
  double xmin = 1.0e-200, xmax = 1e300;
  double ab02, ab12, ab22;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  if (_saemFixedIdx[2] == 1) {
    ab22 = _saemFixedValue[2];
  } else {
    ab22 = ab[curi++];
  }
  double pw = toPow(ab22);
  for (i=0, sum=0; i<_saemLen; ++i) {
    ft = _saemFtCache[i];
    ytr = _saemYtrCache[i];
    fa = _saemFaCache[i];
    if (_saemAddProp == 1){
      g = ab02*ab02 + ab12*ab12*pow(fa, pw);
    } else {
      double ab0 = ab02*ab02;
      double ab1 = ab12*ab12;
      g = sqrt(ab0*ab0 + ab1*ab1*pow(fa, 2*pw));
    }
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur * cur + 2*log(g);
  }
  *fx = sum;
}

// Power only
void objD(double *ab, double *fx) {
  ensureSaemFixedTransformCache();
  int i;
  double g, sum, cur, ft, ytr;
  double xmin = 1.0e-200, xmax = 1e300;
  double ab02, ab12;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  double pw = toPow(ab12);
  double fa;
  for (i=0, sum=0; i<_saemLen; ++i) {
    ft = _saemFtCache[i];
    ytr = _saemYtrCache[i];
    fa = _saemFaAdjustCache[i];
    g = ab02*ab02*pow(fa, pw);
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur * cur + 2*log(g);
  }
  *fx = sum;
}

// add+_saemLambda only
void objE(double *ab, double *fx) {
  int i;
  double g, sum, cur, ft, ytr;
  double xmin = 1.0e-200, xmax = 1e300;
  double ab02, ab12;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  double lambda = toLambda(ab12);
  for (i=0, sum=0; i<_saemLen; ++i) {
    // nelder_() does not al_saemLow _saemLower bounds; we force ab[] be positive here
    ft = _powerD(_saemFptr[i],  lambda, _saemYj, _saemLow, _saemHi);
    ytr = _powerD(_saemYptr[i], lambda, _saemYj, _saemLow, _saemHi);
    g = ab02*ab02;
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur * cur + 2*log(g);
  }
  *fx = sum;
}

// prop+_saemLambda only
void objF(double *ab, double *fx) {
  int i;
  double g, sum, cur, ft, ytr, fa;
  double xmin = 1.0e-200, xmax = 1e300;
  double ab02, ab12;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  double lambda = toLambda(ab12);
  for (i=0, sum=0; i<_saemLen; ++i) {
    // nelder_() does not al_saemLow _saemLower bounds; we force ab[] be positive here
    ft = _powerD(_saemFptr[i],  lambda, _saemYj, _saemLow, _saemHi);
    ytr = _powerD(_saemYptr[i], lambda, _saemYj, _saemLow, _saemHi);
    fa = handleF(_saemPropT, ft, _saemFptr[i], false, true);
    g = ab02*ab02*fa;
    if (g == 0) g = 1;
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur * cur + 2*log(g);
  }
  *fx = sum;
}

// pow+_saemLambda only
void objG(double *ab, double *fx) {
  int i;
  double g, sum, cur, ft, ytr, fa;
  double xmin = 1.0e-200, xmax = 1e300;
  double ab02, ab12, ab22;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  if (_saemFixedIdx[2] == 1) {
    ab22 = _saemFixedValue[2];
  } else {
    ab22 = ab[curi++];
  }
  double lambda = toLambda(ab22);
  double pw = toPow(ab12);
  for (i=0, sum=0; i<_saemLen; ++i) {
    // nelder_() does not al_saemLow _saemLower bounds; we force ab[] be positive here
    ft = _powerD(_saemFptr[i],  lambda, _saemYj, _saemLow, _saemHi);
    ytr = _powerD(_saemYptr[i], lambda, _saemYj, _saemLow, _saemHi);
    fa = handleF(_saemPropT, ft, _saemFptr[i], false, true);
    g = ab02*ab02*pow(fa, pw);
    if (g == 0) g = 1.0;
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur * cur + 2*log(g);
  }
  *fx = sum;
}

// add + prop + _saemLambda
void objH(double *ab, double *fx) {
  int i;
  double g, sum, cur, fa;
  double xmin = 1.0e-200, xmax = 1e300, ft, ytr;
  double ab02, ab12, ab22;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  if (_saemFixedIdx[2] == 1) {
    ab22 = _saemFixedValue[2];
  } else {
    ab22 = ab[curi++];
  }
  double lambda = toLambda(ab22);
  for (i=0, sum=0; i<_saemLen; ++i) {
    // nelder_() does not al_saemLow _saemLower bounds; we force ab[] be positive here
    ft = _powerD(_saemFptr[i],  lambda, _saemYj, _saemLow, _saemHi);
    ytr = _powerD(_saemYptr[i], lambda, _saemYj, _saemLow, _saemHi);
    // focei: rx_r_ = eff^2 * prop.sd^2 + add_sd^2
    // focei g = sqrt(eff^2*prop.sd^2 + add.sd^2)
    fa = handleF(_saemPropT, ft, _saemFptr[i], false, false);
    if (_saemAddProp == 1) {
      g = ab02*ab02 + ab12*ab12*fa;
    } else {
      double ab0 = ab02*ab02;
      double ab1 = ab12*ab12;
      g = sqrt(ab0*ab0 + ab1*ab1*fa*fa);
    }
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur*cur + 2*log(g);
  }
  *fx = sum;
}

// add + pow + _saemLambda
void objI(double *ab, double *fx) {
  int i;
  double g, sum, cur, fa=1.0;
  double xmin = 1.0e-200, xmax = 1e300, ft, ytr;
  double ab02, ab12, ab22, ab32;
  int curi = 0;
  if (_saemFixedIdx[0] == 1) {
    ab02 = _saemFixedValue[0];
  } else {
    ab02 = ab[curi++];
  }
  if (_saemFixedIdx[1] == 1) {
    ab12 = _saemFixedValue[1];
  } else {
    ab12 = ab[curi++];
  }
  if (_saemFixedIdx[2] == 1) {
    ab22 = _saemFixedValue[2];
  } else {
    ab22 = ab[curi++];
  }
  if (_saemFixedIdx[3] == 1) {
    ab32 = _saemFixedValue[3];
  } else {
    ab32 = ab[curi++];
  }

  double lambda = toLambda(ab32);
  double pw = toPow(ab22);
  for (i=0, sum=0; i<_saemLen; ++i) {
    // nelder_() does not al_saemLow _saemLower bounds; we force ab[] be positive here
    ft = _powerD(_saemFptr[i],  lambda, _saemYj, _saemLow, _saemHi);
    ytr = _powerD(_saemYptr[i], lambda, _saemYj, _saemLow, _saemHi);
    fa = handleF(_saemPropT, ft, _saemFptr[i], false, false);
    if (_saemAddProp == 1) {
      g = ab02*ab02 + ab12*ab12*pow(fa, pw);
    } else {
      double ab0 = ab02*ab02;
      double ab1 = ab12*ab12;
      fa = pow(fa, pw);
      g = sqrt(ab0*ab0 + ab1*ab1*fa*fa);
    }
    if (g < xmin) g = xmin;
    if (g > xmax) g = xmax;
    cur = (ytr-ft)/g;
    sum += cur*cur + 2*log(g);
  }
  *fx = sum;
}

int _saemItmax = 100;
double _saemTol = 1e-4;
int _saemType = 1;

static inline void _saemOpt(int n, double *pxmin) {
  if (n == 0) return;
  if (n == 1) {
    // Use R's optimize for unidimensional optimization
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    Environment nlmixr2 = loadNamespace("nlmixr2est");
    Function optimize1 = nlmixr2[".saemOpt1"];
    NumericVector par0(1);
    par0[0] = _saemStart[0];
    double x0 = as<double>(optimize1(par0));
    pxmin[0] = x0;
  } else {
    if (_saemType == 1) {
      int iconv, it, nfcall, iprint=0, itmax=_saemItmax*n;
      double ynewlo;
      nelder_fn(_saemFn, n, _saemStart, _saemStep, itmax, _saemTol, 1.0, 2.0, .5,
                &iconv, &it, &nfcall, &ynewlo, pxmin, &iprint);
    } else if (_saemType == 2) {
      // Try Newoua
      Function loadNamespace("loadNamespace", R_BaseNamespace);
      Environment nlmixr2 = loadNamespace("nlmixr2est");
      Function newuoa = nlmixr2[".newuoa"];
      NumericVector par0(n);
      for (int i = n; i--;) {
        par0[i] = _saemStart[i];
      }
      List ret = newuoa(_["par"] = par0, _["fn"] = nlmixr2[".saemResidF"],
                        _["control"]=List::create(_["rhoend"]=_saemTol,
                                                  _["maxfun"]=_saemItmax*n*n));
      double f = as<double>(ret["value"]);
      if (ISNA(f)) {
        RSprintf("newoua failed, switch to nelder-mead\n");
        int iconv, it, nfcall, iprint=0, itmax=_saemItmax*n;
        double ynewlo;
        nelder_fn(_saemFn, n, _saemStart, _saemStep, itmax, _saemTol, 1.0, 2.0, .5,
                  &iconv, &it, &nfcall, &ynewlo, pxmin, &iprint);
      } else {
        NumericVector x = ret["x"];
        for (int i = n; i--;) {
          pxmin[i] = x[i];
        }
      }
    }
  }
}

extern "C" SEXP _saemResidF(SEXP v) {
  SEXP ret = PROTECT(Rf_allocVector(REALSXP, 1));
  _saemFn(REAL(v),REAL(ret));
  UNPROTECT(1);
  return ret;
}


struct mcmcphi {
  int nphi;
  uvec i;
  mat Gamma_phi;
  mat Gdiag_phi;
  mat IGamma2_phi;
  mat mprior_phiM;
};

struct mcmcaux {
  int nM;
  uvec indio;
  //double sigma2;  //not needed?
  vec y;
  mat evtM;
  List optM;
};


uvec getObsIdx(umat m) {
  uvec x;
  x.set_size(0);

  for (unsigned int b=0; b<m.n_rows; ++b) {
    uvec i=linspace<uvec>(m(b,0), m(b,1), m(b,1) - m(b,0) + 1);
    x = join_cols(x, i);
  }
  return x;
}


// phi0 objective for the general-likelihood direct optimization (bounded bobyqa
// via .boundedResidOpt); gPhi0Self is the active SAEM object.  R-callable through
// Rcpp::InternalFunction.  Defined after the class body.
class SAEM;
static SAEM* gPhi0Self = nullptr;
static double gPhi0ObjR(Rcpp::NumericVector p);
// Coordinate-descent state for the normal-model phi0 direct optimization:
// gPhi0Work is the full phi0 vector, gPhi0Coord the coordinate optimize() varies.
static arma::vec gPhi0Work;
static int gPhi0Coord = 0;
// Free (non-FIXED) phi0 coordinates the bounded optimizer varies; gPhi0Full holds
// the full phi0 vector so FIXED coordinates keep their ini value in the objective.
static arma::vec gPhi0Full;
static std::vector<int> gPhi0FreeIx;
static double gPhi0Obj1DR(double x);
// Shared state for the multivariate phi0 refinements (nelder-mead and newuoa).
// Both are unbounded, so the objective clamps each candidate into the trust
// region before evaluating it, and both need their evaluation budget enforced
// here: nelder_fn's itmax counts IMPROVING iterations rather than evaluations,
// and newuoa's own maxfun stop still returns whatever point it was holding.  The
// objective therefore tracks the best point it actually saw, and the caller
// takes that rather than whatever the optimizer reports.
static arma::vec gPhi0Lo, gPhi0Hi, gPhi0RefBest;
static int gPhi0RefEvalMax = 0, gPhi0RefEvalN = 0;
static double gPhi0RefBestF = 0.0;
static double gPhi0RefObj(const double *p);
static void gPhi0NmFn(double *p, double *fx);
static double gPhi0RefObjR(Rcpp::NumericVector p);

// phi1 objective for the general-likelihood Laplace-corrected direct
// optimization (Phase 4, SAEM general-likelihood theta plan) -- the phi1
// sibling of gPhi0Self/gPhi0ObjR above, same bounded-bobyqa/.boundedResidOpt
// wiring.  Defined after the class body.
static SAEM* gPhi1Self = nullptr;
static double gPhi1ObjR(Rcpp::NumericVector p);
static arma::vec gPhi1Full;
static std::vector<int> gPhi1FreeIx;
// Diagnostic: how many times refinePhi1Lik actually ran, so a test can prove
// the mechanism executed rather than infer it from matching estimates alone
// (evaluation criterion #2 in the plan).
static long _saemPhi1RefineN = 0;

// Fill an armadillo mat/vec from rxode2's threefry engine (the current seeded
// stream).  Used for the MCMC proposals; the saem ODE solve does not draw from
// the engine, so these do not interfere with user_fn.
static inline void _saemFillNormEng(arma::mat &m) {
  for (arma::uword ii = 0; ii < m.n_elem; ++ii) m(ii) = rxNormEng(0.0, 1.0);
}
static inline void _saemFillNormEng(arma::vec &v) {
  for (arma::uword ii = 0; ii < v.n_elem; ++ii) v(ii) = rxNormEng(0.0, 1.0);
}
static inline void _saemFillUnifEng(arma::vec &v) {
  for (arma::uword ii = 0; ii < v.n_elem; ++ii) v(ii) = rxUnifEng(0.0, 1.0);
}
// Iteration-indexed threefry stream seed for a do_mcmc proposal block: a pure
// mixing function of (baseSeed, kiter, method, u, k1) so every iteration's RNG
// is independent of the total iteration count (dynamic phase extension) and
// independent of thread scheduling.  Pins thread 0 (serial draw).  `mixIdx`
// (0 for non-mixture / a 1-based component index for parallel per-component
// chains) offsets the seed so each mixture component draws an independent
// threefry stream instead of the identical proposals that collapse the mixture;
// a bare add suffices since threefry streams for distinct seeds are independent.
static inline void _saemSeedDoMcmc(uint32_t baseSeed, int kiter, int method, int u, int k1,
                                   int mixIdx = 0) {
  setRxThreadId(0);
  uint32_t s = baseSeed;
  s = s * 2654435761u + 0x00006D0Cu;   // "do_mcmc" namespace tag
  s = s * 2654435761u + (uint32_t)kiter;
  s = s * 2654435761u + (uint32_t)method;
  s = s * 2654435761u + (uint32_t)u;
  s = s * 2654435761u + (uint32_t)k1;
  s += (uint32_t)mixIdx;               // per-component stream offset
  nmSetSeedEng1(s);
}

// Iteration-indexed threefry stream seed for the censored-DV data-augmentation
// draw (see simCensDv()/augmentCensY() below), independent of the do_mcmc
// proposal streams above.  Keyed by (kiter, chain k, endpoint b, mixIdx) so
// every (iteration, chain, endpoint, component) combination gets its own
// stream; a fit is reproducible given saemControl(seed=).
static inline void _saemSeedCensAug(uint32_t baseSeed, int kiter, int k, int b,
                                    int mixIdx) {
  setRxThreadId(0);
  uint32_t s = baseSeed;
  s = s * 2654435761u + 0x63656E73u;   // "cens" namespace tag
  s = s * 2654435761u + (uint32_t)kiter;
  s = s * 2654435761u + (uint32_t)k;
  s = s * 2654435761u + (uint32_t)b;
  // per-component stream offset (mixIdx can be -1); multiply-folded like every
  // other field above -- a bare `+=` here collides whenever b and mixIdx+1
  // sum to the same value (e.g. b=0,mixIdx=1 vs b=1,mixIdx=0).
  s = s * 2654435761u + (uint32_t)(mixIdx + 1);
  nmSetSeedEng1(s);
}

// Simulate the "true" value of a censored (M3/M4) observation from the
// truncated normal implied by the current transformed prediction/residual SD
// -- data augmentation (Samson, Lavielle & Mentre 2006) so the M-step
// residual SSR sees a draw from the censored region instead of the recorded
// LOQ/limit.  M2 rows carry a real measurement and are returned unchanged.
// cens/limDv/lim follow the doCensNormal1() convention (all on the
// transformed scale here); sd is the endpoint's current residual SD.  The
// draw itself is rxTruncNorm() (truncNorm.h) -- the same Botev (2015)
// algorithm CWRES's censored-observation simulation uses (censResid.h's
// truncnorm(), via rxode2's rxRmvn) -- rather than a plain inverse-CDF draw,
// which loses precision once the truncation bounds are a few SDs from the
// mean (the regime a BQL row's bound often sits in).
static inline double simCensDv(double cens, double limDv, double lim, double f,
                               double sd) {
  if (!(cens == 1.0 || cens == -1.0)) return limDv;
  double lo = R_NegInf, hi = R_PosInf;
  if (R_FINITE(lim) && !ISNA(lim)) {
    // M4: truncation interval is between limDv (the LOQ) and lim (the other,
    // informative bound); cens picks which side is which.
    if (cens > 0) { lo = lim;   hi = limDv; } else { lo = limDv; hi = lim; }
  } else if (cens > 0) {
    // M3, left-censored: y <= limDv
    hi = limDv;
  } else {
    // M3, right-censored: y >= limDv
    lo = limDv;
  }
  double zl = R_FINITE(lo) ? (lo - f) / sd : R_NegInf;
  double zu = R_FINITE(hi) ? (hi - f) / sd : R_PosInf;
  if (!(zu > zl)) return limDv;        // degenerate/inverted bound: keep the historical value
  return f + sd * rxTruncNorm(zl, zu);
}

// Phase 4 (SAEM general-likelihood theta plan): phi1Objective (a class
// method, below) solves the odeSlotHess2 peer directly, so it needs the
// process rx_solve* -- defined further down in this same translation unit
// (after the class body), forward-declared here so the class can reference
// it. Every field/function this pulls in (rxHess2, OdeSwapScope/CmtScope,
// odeSwapSolveInd, ...) already comes from inner.h/odeSwap.h, included above.
extern rx_solve* _rx;

// Phase 4: THETA[k]/ETA[k] -> phi column maps and pool-readiness state,
// shared between phi1Objective (a method, uses odeSlotHess2/odeSlotPred) and
// user_function (a free function, uses odeSlotPred for its own per-row
// solve) -- see the fuller doc comment and real definitions after the class
// body, alongside _saemPhi1PoolActive.
extern bool _saemPhi1PoolReady;
extern bool _saemPhi1UseAnalyticHess;
extern arma::ivec _saemPhi1H2ThetaKind;
extern arma::ivec _saemPhi1H2ThetaCol;
extern arma::vec _saemPhi1H2ThetaFixedVal;
extern arma::ivec _saemPhi1H2EtaCol;
extern bool _saemPhi1WantHessian;
extern int _saemPhi1PredOffset;
// 0-based DV parameter-slot position, resolved by .saemPhi1TargetMap
// (R/saemPhi1Inner.R) purely as a readiness check -- confirms the compiled
// general-likelihood model actually declares a DV parameter. DV itself is
// supplied by the ordinary solve setup (see rxUiGet.saemInParsAndMuRefCovariates,
// R/saemRxUiGet.R), not written here. -1 when not resolved.
extern int _saemPhi1DvCol;
extern int _saemPhi1DvColHess2;
extern int _saemPhi1H2PredOffset;
extern int _saemPhi1H2HessOffset;
extern arma::uvec _saemPhi1I0;
extern arma::uvec _saemPhi1I1;

// class def starts
class SAEM {
  typedef mat (*user_funct) (const mat&, const mat&, const List&);

public:

  SAEM() {
    user_fn = NULL;
  }

  ~SAEM() {}

  // Total observation -log-likelihood at candidate fixed-effect (phi0) values p,
  // holding the current phi1 samples fixed (general-likelihood / distribution==4:
  // the model prediction column is the per-observation log-likelihood).  Summed
  // over all chains, which is the SAEM stochastic-approximation objective.
  double phi0Objective(double *p) {
    mat phiCand = phiM;
    for (int c = 0; c < nphi0; c++) {
      phiCand.col(i0(c)).fill(p[c]);
    }
    mat fMat = user_fn(phiCand, evt, optM);
    double v;
    if (distribution == 4) {
      // general log-likelihood: the prediction column IS the per-obs loglik
      v = -accu(fMat.col(0));
    } else {
      // normal models (nonMuTheta="regress"): residual objective from f
      v = phi0NormalSSR(fMat.col(0));
    }
    if (!std::isfinite(v)) return 1e300;
    return v;
  }

  // Side-effect-free observation -log-likelihood objective for the normal-model
  // phi0 direct optimization (nonMuTheta="regress").  This is the SAME
  // per-observation Gaussian -loglik the MCMC acceptance uses (arDYF /
  // independent branch): 0.5*((yt-ft)/g)^2 + log(g), with residual SD
  // g = ares(b) + bres(b)*|ft| (saem.cpp:1457) and TBS transform ft via
  // _powerD.  Properly normalized across error models and endpoints (unlike a
  // raw SSR), summed over MCMC chains.  Does NOT touch the AR accumulators /
  // M-step state (unlike arResk), so it is safe to call inside the optimizer.
  // The TBS data jacobian is omitted: it is constant in phi0, so it does not
  // change the argmin.  AR correlation is left to the M-step.
  double phi0NormalSSR(const vec &fall) {
    const double double_xmin = 1.0e-200, xmax = 1e300;
    double v = 0.0;
    for (int k = 0; k < nmc; k++) {
      vec fk = fall.subvec(k * ntotal, (k + 1) * ntotal - 1);
      fk = fk(ix_sorting);
      for (int b = 0; b < nendpnt; b++) {
        int nb = (int)(y_offset(b + 1) - y_offset(b));
        for (int i = 0; i < nb; i++) {
          double fi = fk(y_offset(b) + i);
          double ft = _powerD(fi, lambda(b), yj(b), low(b), hi(b));
          double yt = hasFixedObsTransform ? ysTrans(y_offset(b) + i)
            : _powerD(ys(y_offset(b) + i), lambda(b), yj(b), low(b), hi(b));
          double g = ares(b) + bres(b) * std::fabs(ft);
          if (g == 0.0) g = 1.0;
          else if (g < double_xmin) g = double_xmin;
          else if (g > xmax) g = xmax;
          double e = yt - ft;
          v += 0.5 * (e / g) * (e / g) + std::log(g);
        }
      }
    }
    return v;
  }

  // saemix "ind.fix10" step: refine the fixed-effect-only (phi0) parameters of a
  // general-likelihood model by a direct optimization of the observation
  // log-likelihood (seeded at the current phi(theta) values), then apply a
  // stochastic-approximation update and keep MCOV0 consistent so the next
  // iteration's mprior_phi0 = COV0*MCOV0 reproduces it.  Only the intercept-only
  // (no phi0 covariate) case is handled.
  //
  // phi0 does not enter the ODE, so the states are solved once and then held
  // fixed (ODE-freeze) while phi0 is optimized -- each objective evaluation only
  // recomputes the log-likelihood.  The model emits no analytic d(ll)/d(phi0), so
  // the optimization is derivative-free (nelder-mead / newuoa), reusing the shared
  // _saemOpt driver selected by the `type` control (_saemType).
  // Detect (once) whether any phi0 param changes the structural prediction f.
  // Perturb each phi0 column and re-solve; if f never moves, phi0 touches only
  // the residual/likelihood and the ODE can be frozen during its optimization.
  bool phi0AffectsOde() {
    bool savedFreeze = _saemFreezeOde;
    _saemFreezeOde = false;
    vec f0 = user_fn(phiM, evt, optM).col(0);
    double maxd = 0.0;
    for (int c = 0; c < nphi0; c++) {
      mat pp = phiM;
      pp.col(i0(c)) += 0.1;
      vec f1 = user_fn(pp, evt, optM).col(0);
      double d = arma::abs(f1 - f0).max();
      if (std::isfinite(d) && d > maxd) maxd = d;
    }
    // restore states at the true phiM
    { mat _t = user_fn(phiM, evt, optM); (void)_t; }
    _saemFreezeOde = savedFreeze;
    return maxd > 1e-8;
  }

  // Does the phi0 refinement need a LIVE re-solve?
  //
  // Freezing the ODE (saemix ind.fix10) is only valid for a phi0 parameter the
  // SOLVE does not see -- a likelihood SD, a residual parameter.  An IOV
  // magnitude theta is a phi0 parameter that drives the STRUCTURAL model, and
  // with the solve frozen the objective is exactly constant in it, so the
  // bounded optimizer walks to its upper bound and the SA update turns that
  // into an unbounded runaway (#1000).  Answer by measuring the frozen/live
  // discrepancy on the free phi0 columns, never by assuming from
  // `distribution`.  Must run after gPhi0FreeIx is filled.
  bool phi0NeedsLiveSolve() {
    bool savedFreeze = _saemFreezeOde;
    bool needs = false;
    _saemFreezeOde = false;
    { mat _t = user_fn(phiM, evt, optM); (void)_t; }   // states at phiM
    for (size_t fi = 0; fi < gPhi0FreeIx.size() && !needs; ++fi) {
      mat pp = phiM;
      pp.col(i0(gPhi0FreeIx[fi])) += 0.1;
      _saemFreezeOde = true;
      vec ff = user_fn(pp, evt, optM).col(0);   // frozen: states still at phiM
      _saemFreezeOde = false;
      vec fl = user_fn(pp, evt, optM).col(0);   // live
      double d = arma::abs(ff - fl).max();
      if (!std::isfinite(d) || d > 1e-8) needs = true;
      { mat _t = user_fn(phiM, evt, optM); (void)_t; } // restore states at phiM
    }
    _saemFreezeOde = savedFreeze;
    return needs;
  }

  // Impose the two-level (IOV) equality constraint on an Omega-shaped matrix.
  //
  // Columns sharing an omegaPool group id are the same occasion parameter at
  // different occasion levels, so their variances are one parameter.  Replace
  // each group's diagonal entries by their mean.  A user-FIXED occasion
  // variance needs no special case here: Gamma2_phi1fixedIx is restored after
  // this runs, so the pin wins.
  void poolOmegaGroups(mat &G) {
    if (omegaPool.n_elem != (unsigned int)nphi1) return;
    if (omegaPool.n_elem == 0) return;
    unsigned int maxg = omegaPool.max();
    for (unsigned int g = 1; g <= maxg; ++g) {
      uvec ix = find(omegaPool == g);
      if (ix.n_elem < 2) continue;
      double m = 0.0;
      for (unsigned int j = 0; j < ix.n_elem; ++j) m += G(ix(j), ix(j));
      m /= (double)ix.n_elem;
      for (unsigned int j = 0; j < ix.n_elem; ++j) G(ix(j), ix(j)) = m;
    }
  }

  // Re-run the uninformative-eta test at the current estimates.
  //
  // The test (R/uninformativeEtas.R) perturbs each eta and asks whether the prediction
  // moves, and it runs ONCE, at the INITIAL estimates.  Poor initial estimates make it
  // answer a question about the estimates rather than about the eta, and whatever it
  // decides is carried for the whole fit.  Probing again at the end of burn-in, with
  // theta and Omega where SAEM has moved them, gets the verdict the data supports.
  //
  // Predictions come from user_fn -- the same evaluation the loop already does every
  // iteration -- so this adds no rxSolve of its own and no R callback: the cached _rx
  // and rxode2's solve state are untouched.  user_fn draws no random numbers, so the
  // RNG stream (and hence the rest of the fit) is unchanged.
  // The revisit writes into the mask and reads the flat prediction vector through
  // ix_idM, so both have to have the layout it assumes; anything else means the two
  // sides disagree and writing would be writing somewhere else entirely.
  bool ueRevisitLayoutOk() {
    if (ueRevisitCols.n_elem == 0 || N <= 0 || nmc <= 0 || ntotal <= 0) return false;
    if (current_saem_state->_saemUE.n_rows != (unsigned int)(N * nmc) ||
        current_saem_state->_saemUE.n_cols != (unsigned int)nphi) return false;
    // ix_idM is indexed by PSEUDO-subject (chain c, subject i) at row c*N+i and already
    // carries that chain's offset into the flat prediction vector -- do not re-derive it
    if (ix_idM.n_rows != (unsigned int)(N * nmc)) return false;
    return ueDelta.n_elem == ueRevisitCols.n_elem;
  }

  // Read one solved perturbation level back into the accumulators: chain `c` of the
  // solve carries level `l` of column `jj`.
  void ueProbeAccum(const vec &g, unsigned int jj, int c, int l,
                    mat &ret, mat &scl, umat &bad) {
    const double sgn[3] = {1.0, -2.0, 1.0};  // pred(-) + pred(+) - 2*pred(0)
    for (int i = 0; i < N; ++i) {
      unsigned int row = (unsigned int)(c * N + i);
      unsigned int st = ix_idM(row, 0), en = ix_idM(row, 1);
      // en < st is how a subject with no observations arrives (end = start - 1,
      // which wraps); either way there is nothing to read and nothing to judge
      if (en < st || en >= g.n_elem) { bad(i, jj) = 1; continue; }
      for (unsigned int idx = st; idx <= en; ++idx) {
        double p = g(idx);
        // user_fn substitutes 1e99 for a NaN prediction, so a failed solve
        // arrives finite and huge rather than as NaN -- catch the sentinel.
        if (!R_finite(p) || std::abs(p) >= 1.0e99) { bad(i, jj) = 1; continue; }
        ret(i, jj) += sgn[l] * p;
        double a = std::abs(p);
        if (a > scl(i, jj)) scl(i, jj) = a;
      }
    }
  }

  // Accumulate the second-difference statistic for one eta column, spreading the three
  // perturbation levels across the MCMC chains so each solve carries up to three of them.
  void ueProbeCol(unsigned int jj, unsigned int col, const mat &phiBase,
                  mat &ret, mat &scl, umat &bad) {
    // The probe half-width comes from the INITIAL Omega, exactly as the first test's
    // does, so the two evaluations differ only in theta -- which is the whole point.
    // Reading the CURRENT Omega instead would make the revisit fight itself: a frozen
    // eta contributes nothing to its own variance, so Omega shrinks, the probe narrows,
    // and the eta looks even less informative the longer it has been frozen.
    double delta = ueDelta(jj);
    if (!R_finite(delta) || delta <= 0) return;
    const int nLev = std::min(nmc, 3);      // perturbation levels carried per solve
    const double lev[3] = {-1.0, 0.0, 1.0};
    for (int l0 = 0; l0 < 3; l0 += nLev) {
      mat phiProbe(N * nmc, nphi);
      for (int c = 0; c < nmc; ++c) {
        int l = std::min(l0 + (c % nLev), 2);   // pad any unused tail chains
        phiProbe.rows(c * N, c * N + N - 1) = phiBase;
        phiProbe.submat(c * N, col, c * N + N - 1, col) += lev[l] * delta;
      }
      vec g = user_fn(phiProbe, evt, optM).col(0);
      for (int c = 0; c < nLev && l0 + c < 3; ++c) {
        ueProbeAccum(g, jj, c, l0 + c, ret, scl, bad);
      }
    }
  }

  // Same verdict as _nlmixr2est_uninformativeEta: only take "uninformative" when the
  // predictions it is built from are real.
  void ueApplyMask(unsigned int jj, unsigned int col,
                   const mat &ret, const mat &scl, const umat &bad) {
    for (int i = 0; i < N; ++i) {
      bool havePred = !bad(i, jj) && R_finite(ret(i, jj)) &&
        R_finite(scl(i, jj)) && scl(i, jj) > ueTol;
      double m = ((std::abs(ret(i, jj)) > ueTol) || !havePred) ? 1.0 : 0.0;
      double was = current_saem_state->_saemUE(i, col);
      if (was == 0.0 && m == 1.0) ueRevisitUnfroze++;
      else if (was == 1.0 && m == 0.0) ueRevisitFroze++;
      for (int c = 0; c < nmc; ++c) {
        current_saem_state->_saemUE(c * N + i, col) = m;
      }
    }
  }

  void revisitUninformativeEtas() {
    if (!ueRevisitLayoutOk()) return;
    ueRevisitRan = 1;
    const unsigned int nc = ueRevisitCols.n_elem;

    // base phi at eta = 0: each subject's mu-referenced population value
    mat phiBase(N, nphi, fill::zeros);
    phiBase.cols(i1) = mprior_phi1;
    if (nphi0 > 0) phiBase.cols(i0) = mprior_phi0;

    mat ret(N, nc, fill::zeros);            // the second-difference statistic
    mat scl(N, nc, fill::zeros);            // largest |pred| it is built from
    umat bad(N, nc, fill::zeros);           // a solve that failed anywhere in the cell

    bool savedFreeze = _saemFreezeOde;
    _saemFreezeOde = false;                 // the probe needs a live re-solve
    for (unsigned int jj = 0; jj < nc; ++jj) {
      unsigned int col = ueRevisitCols(jj);
      if (col < (unsigned int)nphi) ueProbeCol(jj, col, phiBase, ret, scl, bad);
    }
    { mat _t = user_fn(phiM, evt, optM); (void)_t; }  // restore states at the true phiM
    _saemFreezeOde = savedFreeze;

    for (unsigned int jj = 0; jj < nc; ++jj) {
      unsigned int col = ueRevisitCols(jj);
      if (col < (unsigned int)nphi) ueApplyMask(jj, col, ret, scl, bad);
    }
  }

  void refinePhi0Lik(unsigned int kiter, const vec &pas) {
    if (nphi0 <= 0) return;
    // A user-FIXED phi0 theta must not be touched here.  Once this refinement
    // owns phi0 the stochastic update that restores MCOV0(fixedIx0) is skipped
    // (skipStochPhi0), so anything this function writes to a fixed entry is
    // never put back and the theta drifts off its ini value.  fixedIx0 indexes
    // phi0 columns (refinePhi0Lik only handles the intercept-only phi0 case).
    std::vector<bool> phi0Fix((size_t)nphi0, false);
    for (unsigned int j = 0; j < fixedIx0.n_elem; ++j) {
      if (fixedIx0(j) < (unsigned int)nphi0) phi0Fix[(size_t)fixedIx0(j)] = true;
    }
    gPhi0FreeIx.clear();
    for (int c = 0; c < nphi0; ++c) {
      if (!phi0Fix[(size_t)c]) gPhi0FreeIx.push_back(c);
    }
    if (gPhi0FreeIx.empty()) return;
    // Snapshot the fixed MCOV0 entries: the closing least-squares update rewrites
    // all of MCOV0 from mprior_phi0, which redistributes a fixed theta's value
    // across the design even when its coordinate never moved.
    vec mcov0Fixed;
    if (fixedIx0.n_elem > 0) mcov0Fixed = vec(MCOV0(jcov0(fixedIx0)));
    // Decide whether to freeze the ODE during the phi0 optimization.  General-
    // likelihood phi0 params (a likelihood SD) never enter the ODE.  For a
    // normal model under nonMuTheta="regress", phi0 thetas that drive the ODE
    // (ka, V, ...) need a LIVE re-solve each evaluation, but phi0 params that
    // touch only the residual/likelihood (not f) should be frozen -- solve once,
    // recompute only the objective, exactly like npag's ELS residual step.  The
    // f-sensitivity is detected once (perturb each phi0, see if f moves).
    _saemFreezeOde = false;
    { mat _tmp = user_fn(phiM, evt, optM); (void)_tmp; }  // establish states
    bool doFreeze;
    if (distribution == 4) {
      // NOT unconditionally frozen: a general-likelihood model can still carry a
      // phi0 that drives the solve (an IOV magnitude theta), and freezing makes
      // the objective constant in it -- see phi0NeedsLiveSolve() and #1000.
      if (_phi0NeedsLive < 0) _phi0NeedsLive = phi0NeedsLiveSolve() ? 1 : 0;
      doFreeze = (_phi0NeedsLive == 0);
    } else if (nonMuThetaRegress) {
      if (_phi0OdeSensitive < 0) _phi0OdeSensitive = phi0AffectsOde() ? 1 : 0;
      doFreeze = (_phi0OdeSensitive == 0);
    } else {
      doFreeze = false;
    }
    _saemFreezeOde = doFreeze;
    // optimize phi0 with the BOUNDED bobyqa (.boundedResidOpt), honoring the
    // ini-block bounds of the phi0 thetas.  An unbounded method (newuoa/
    // nelder-mead) could push a phi0 like a likelihood SD into an invalid region.
    gPhi0Self = this;
    Rcpp::NumericVector par0(nphi0), lo(nphi0), hi(nphi0);
    // Normal-model phi0 thetas (nonMuTheta="regress") drive the ODE, so the
    // observation objective has huge NaN plateaus far from the current value;
    // an unbounded bobyqa span breaks there.  Constrain each step to a LOCAL
    // trust region around the current value (intersected with any ini bounds)
    // so the SA iteration refines it gradually, like a clamped regression step.
    // localTrust also selects the optimizer; trustBounds only clamps the step.
    bool localTrust = (distribution != 4) && nonMuThetaRegress;
    // A general-likelihood phi0 that drives the solve (an IOV magnitude) gets the
    // same absolute trust region: its ini bounds are typically (0, Inf), and an
    // unbounded span over an ODE-driven objective is what let it run away (#1000).
    // A general-likelihood phi0 the solve never sees keeps the original wide bounds.
    bool trustBounds = localTrust || (distribution == 4 && _phi0NeedsLive == 1);
    for (int c = 0; c < nphi0; c++) {
      par0[c] = mprior_phi0(0, c);
      double userLo = ((int)phi0Lower.n_elem == nphi0) ? phi0Lower(c) : R_NegInf;
      double userHi = ((int)phi0Upper.n_elem == nphi0) ? phi0Upper(c) : R_PosInf;
      if (trustBounds) {
        // ABSOLUTE trust radius (not relative to par0): a relative radius lets a
        // param that starts to drift grow its own step and run away.
        double trust = 0.75;
        double tl = par0[c] - trust, th = par0[c] + trust;
        lo[c] = std::isfinite(userLo) ? std::max(userLo, tl) : tl;
        hi[c] = std::isfinite(userHi) ? std::min(userHi, th) : th;
      } else {
        lo[c] = userLo;
        hi[c] = userHi;
      }
    }
    Rcpp::NumericVector xmin(nphi0);
    if (localTrust) {
      // Normal-model phi0 objective is extremely ill-conditioned (a tiny
      // proportional-error SD makes it change by orders of magnitude over a
      // small phi0 step), which breaks bobyqa's quadratic model.  Two
      // derivative-free options within the local trust bounds:
      //   nonMuThetaOpt="optimize"   -- coordinate descent with R's golden-section
      //     optimize(); exact per coordinate but costs sweeps*nphi0*~20 objective
      //     evaluations, each a full ODE re-solve when phi0 drives the model.
      //   nonMuThetaOpt="nelderMead" -- one clamped nelder-mead over all free
      //     coordinates; a fixed, much smaller budget, and it sees the coupling
      //     between coordinates that coordinate descent cannot.
      //   nonMuThetaOpt="newuoa"     -- the same, with newuoa's quadratic model
      //     built from those evaluations instead of a simplex.
      // Both spend at most nonMuThetaMaxEval objective evaluations (enforced by
      // gPhi0RefObj) and both are clamped into the trust region.
      gPhi0Self = this;
      gPhi0Work.set_size(nphi0);
      for (int c = 0; c < nphi0; c++) gPhi0Work[c] = par0[c];
      int nFree = (int)gPhi0FreeIx.size();
      if (nonMuThetaOptType > 0 && nFree > 1) {
        gPhi0Lo.set_size(nphi0);
        gPhi0Hi.set_size(nphi0);
        for (int c = 0; c < nphi0; c++) { gPhi0Lo(c) = lo[c]; gPhi0Hi(c) = hi[c]; }
        gPhi0RefBest.set_size(nFree);
        gPhi0RefEvalN = 0;
        gPhi0RefBestF = 0.0;
        gPhi0RefEvalMax = (nonMuThetaMaxEval > 0) ? nonMuThetaMaxEval : 10*nFree;
        if (nonMuThetaOptType == 1) {
          std::vector<double> st((size_t)nFree), stp((size_t)nFree), xm((size_t)nFree);
          for (int fi = 0; fi < nFree; fi++) {
            int c = gPhi0FreeIx[(size_t)fi];
            st[(size_t)fi] = par0[c];
            xm[(size_t)fi] = par0[c];
            double span = hi[c] - lo[c];
            stp[(size_t)fi] = (R_finite(span) && span > 0.0) ? 0.1*span : 0.1;
          }
          int iconv, it, nfcall, iprint = 0;
          double ynewlo;
          // itmax is deliberately generous; gPhi0RefObj owns the real budget
          nelder_fn(gPhi0NmFn, nFree, st.data(), stp.data(), 100*nFree,
                    nonMuThetaTol, 1.0, 2.0, 0.5,
                    &iconv, &it, &nfcall, &ynewlo, xm.data(), &iprint);
        } else {
          // newuoa needs npt = 2n+1 interpolation points before it can move, so
          // a budget below that leaves it no working evaluations at all.
          int npt = 2*nFree + 1;
          if (gPhi0RefEvalMax < npt + 2) gPhi0RefEvalMax = npt + 2;
          Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
          Rcpp::Function phi0Newuoa = nlmixr2[".saemPhi0Newuoa"];
          Rcpp::InternalFunction fnRef(&gPhi0RefObjR);
          Rcpp::NumericVector parFree(nFree);
          double rhobeg = R_PosInf;
          for (int fi = 0; fi < nFree; fi++) {
            int c = gPhi0FreeIx[(size_t)fi];
            parFree[fi] = par0[c];
            double span = hi[c] - lo[c];
            if (R_finite(span) && span > 0.0 && 0.2*span < rhobeg) rhobeg = 0.2*span;
          }
          if (!R_finite(rhobeg) || rhobeg <= 0.0) rhobeg = 0.2;
          phi0Newuoa(Rcpp::_["par"] = parFree, Rcpp::_["fn"] = fnRef,
                     Rcpp::_["maxfun"] = gPhi0RefEvalMax,
                     Rcpp::_["rhobeg"] = rhobeg,
                     Rcpp::_["rhoend"] = nonMuThetaTol,
                     Rcpp::_["npt"] = npt);
        }
        for (int fi = 0; fi < nFree; fi++) {
          gPhi0Work[gPhi0FreeIx[(size_t)fi]] = gPhi0RefBest(fi);
        }
      } else {
        Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
        Rcpp::Function optimize = stats["optimize"];
        Rcpp::InternalFunction fn1d(&gPhi0Obj1DR);
        for (int sweep = 0; sweep < nonMuThetaSweeps; sweep++) {
          for (int c = 0; c < nphi0; c++) {
            if (phi0Fix[(size_t)c]) continue;
            if (!(hi[c] > lo[c])) continue;
            gPhi0Coord = c;
            Rcpp::List o = optimize(Rcpp::_["f"] = fn1d,
                                    Rcpp::_["lower"] = lo[c],
                                    Rcpp::_["upper"] = hi[c],
                                    Rcpp::_["tol"] = nonMuThetaTol);
            double xm = Rcpp::as<double>(o["minimum"]);
            gPhi0Work[c] = xm;
          }
        }
      }
      for (int c = 0; c < nphi0; c++) xmin[c] = gPhi0Work[c];
    } else {
      Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
      Rcpp::Function boundedOpt = nlmixr2[".boundedResidOpt"];
      Rcpp::InternalFunction fn(&gPhi0ObjR);
      // Optimize only the free coordinates: gPhi0ObjR expands them back into
      // gPhi0Full, which holds the FIXED coordinates at their ini values.
      gPhi0Full.set_size(nphi0);
      for (int c = 0; c < nphi0; c++) gPhi0Full[c] = par0[c];
      int nFree = (int)gPhi0FreeIx.size();
      Rcpp::NumericVector parFree(nFree), loFree(nFree), hiFree(nFree);
      for (int fi = 0; fi < nFree; fi++) {
        int c = gPhi0FreeIx[(size_t)fi];
        parFree[fi] = par0[c];
        loFree[fi] = lo[c];
        hiFree[fi] = hi[c];
      }
      Rcpp::List ret = boundedOpt(Rcpp::_["par"] = parFree, Rcpp::_["fn"] = fn,
                                  Rcpp::_["lower"] = loFree, Rcpp::_["upper"] = hiFree);
      Rcpp::NumericVector rx = ret["x"];
      for (int c = 0; c < nphi0; c++) xmin[c] = par0[c];
      for (int fi = 0; fi < nFree; fi++) xmin[gPhi0FreeIx[(size_t)fi]] = rx[fi];
    }
    _saemFreezeOde = false;
    for (int c = 0; c < nphi0; c++) {
      double cur = mprior_phi0(0, c);
      mprior_phi0.col(c).fill(cur + pas(kiter) * (xmin[c] - cur));
    }
    // MCOV0 is BLOCK structured by LCOV0 -- each lambda row belongs to exactly one
    // phi0 column -- so a single least-squares against all of COV0 is rank
    // deficient whenever nphi0 > 1: with no phi0 covariate every column of COV0 is
    // the same intercept column, and arma warns "solve(): system is singular"
    // every iteration from niter_phi0 on.  It also fills MCOV0 off-structure, so a
    // FIXED phi0 no longer reproduces its value through COV0*MCOV0.  Back-solve
    // each phi0 column against only its own design columns instead.
    for (int c = 0; c < nphi0; c++) {
      uvec li = arma::find(LCOV0.col(c) == 1);
      if (li.n_elem == 0) continue;
      mat Xc = COV0.cols(li);
      vec bc;
      if (arma::solve(bc, Xc.t() * Xc, Xc.t() * mprior_phi0.col(c))) {
        for (unsigned int j = 0; j < li.n_elem; ++j) MCOV0(li(j), c) = bc(j);
      }
    }
    if (fixedIx0.n_elem > 0) MCOV0(jcov0(fixedIx0)) = mcov0Fixed;
  }

  // Phase 4 (SAEM general-likelihood theta plan): Laplace-corrected objective
  // for the phi1 (mu-referenced theta) direct optimization -- the phi1
  // sibling of phi0Objective.  No EBE mode-search: for candidate theta p
  // (nphi1 free coordinates, intercept-only -- guaranteed by _saemPhi1PoolReady/
  // .saemPhi1TargetMap, which declines any phi1 covariate), every one of the
  // nM=N*nmc chain-replicated rows is scored by solving innerHess2
  // (odeSlotHess2) at phi_i_NEW = mu_i(p) + eta_i, where
  // eta_i = phiM(i, i1) - mprior_phi1(subject, .) is that row's CURRENT,
  // UNCHANGED deviation from whatever mu was active when phiM was last
  // sampled by do_mcmc.  Reads back rx_pred_ (the log-density) and
  // rx__d2pred_i_j__ (the exact 2nd-order eta-Hessian AT THAT SUPPLIED
  // POINT, not a re-optimized mode) and scores -2*rx_pred_ + log|H|, with
  // H = -d2pred + Omega^-1 (matching calcEtaHessian's own sign convention,
  // src/inner.cpp).  Summed over every row -- no division by nmc, since this
  // only ever feeds an argmin (matching phi0Objective's own convention).  A
  // row whose solve fails or whose H is not positive-definite makes the
  // whole candidate return phi0Objective's own 1e300 sentinel.  Parallelized
  // the same way inner.cpp's own per-subject loops are.
  // One row's rx_pred_ (summed over its observations) at a supplied eta
  // vector, for the FD-fallback path -- solves odeSlotPred (no sensitivities)
  // fresh at etaVec, under a scope the caller already has open for this row.
  // Sets `bad` on solve failure.
  double phi1PredAt(int i, rx_solving_options_ind *ind, rx_solving_options *op,
                     OdeSwapScope &neqGuard, int nH2Theta,
                     const arma::vec &etaVec, bool &bad) {
    int nEta = (int)etaVec.n_elem;
    for (int k = 0; k < nEta; ++k) setIndParPtr(ind, nH2Theta + k, etaVec(k));
    setIndSolve(ind, -1);
    resetOpBadSolve(op);  // courtesy only; racy under cores>1, not relied on below
    odeSwapSolveInd(odeSlotPred, i);
    if (odeSwapIndBadSolveSlot(op, ind, odeSlotPred)) { bad = true; return 0.0; }
    iniSubjectE(i, 1, ind, op, _rx, rxPred.update_inis);
    double *lhs = neqGuard.lhs();
    double pred = 0.0;
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, j);
      if (getIndEvid(ind, kk) != 0) continue;
      double curT = getTime(kk, ind);
      rxPred.calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
      pred += lhs[_saemPhi1PredOffset];
    }
    return pred;
  }

  // Solve odeSlotHess2 (the exact analytic eta-Hessian model) at eta0 for row
  // i, accumulating rx_pred_ into rowPred and the packed lower-triangular
  // d2pred into H.  Returns false (row is bad) on solve failure, matching the
  // inline `rowBad[i] = 1; continue;` this replaced.
  bool phi1AnalyticHessAt(int i, rx_solving_options_ind *ind, rx_solving_options *op,
                          OdeSwapScope &neqGuard, int nH2Theta,
                          const arma::vec &eta0, double &rowPred, arma::mat &H) {
    int nEta = (int)eta0.n_elem;
    for (int k = 0; k < nEta; ++k) setIndParPtr(ind, nH2Theta + k, eta0(k));
    setIndSolve(ind, -1);
    resetOpBadSolve(op);
    odeSwapSolveInd(odeSlotHess2, i);
    if (odeSwapIndBadSolveSlot(op, ind, odeSlotHess2)) return false;
    iniSubjectE(i, 1, ind, op, _rx, rxHess2.update_inis);
    double *lhs = neqGuard.lhs();
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, j);
      if (getIndEvid(ind, kk) != 0) continue;
      double curT = getTime(kk, ind);
      rxHess2.calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
      rowPred += lhs[_saemPhi1H2PredOffset];
      int r = 0;
      for (int jc = 0; jc < nphi1; ++jc)
        for (int ic = 0; ic <= jc; ++ic) {
          H(ic, jc) += -lhs[_saemPhi1H2HessOffset + r];
          ++r;
        }
    }
    return true;
  }

  // Finite-difference eta-Hessian fallback at eta0 for row i, re-solving
  // odeSlotPred (via phi1PredAt) at small perturbations of eta0 -- the same
  // Shi(2021)-style fallback discipline calcEtaHessian uses, simplified to a
  // fixed step (see phi1Objective's own docs on fdH).  Sets `bad` on any
  // perturbed solve's failure; H accumulates in place.
  void phi1FDHessAt(int i, rx_solving_options_ind *ind, rx_solving_options *op,
                     OdeSwapScope &neqGuard, int nH2Theta, double fdH,
                     const arma::vec &eta0, double rowPred, arma::mat &H, bool &bad) {
    int nEta = (int)eta0.n_elem;
    for (int k = 0; k < nEta && !bad; ++k) {
      arma::vec ep = eta0, em = eta0;
      ep(k) += fdH; em(k) -= fdH;
      double fp = phi1PredAt(i, ind, op, neqGuard, nH2Theta, ep, bad);
      double fm = bad ? 0.0 : phi1PredAt(i, ind, op, neqGuard, nH2Theta, em, bad);
      if (!bad) H(k, k) += -(fp - 2.0 * rowPred + fm) / (fdH * fdH);
    }
    for (int jc = 1; jc < nEta && !bad; ++jc) {
      for (int ic = 0; ic < jc && !bad; ++ic) {
        arma::vec epp = eta0, epm = eta0, emp = eta0, emm = eta0;
        epp(ic) += fdH; epp(jc) += fdH;
        epm(ic) += fdH; epm(jc) -= fdH;
        emp(ic) -= fdH; emp(jc) += fdH;
        emm(ic) -= fdH; emm(jc) -= fdH;
        double fpp = phi1PredAt(i, ind, op, neqGuard, nH2Theta, epp, bad);
        double fpm = bad ? 0.0 : phi1PredAt(i, ind, op, neqGuard, nH2Theta, epm, bad);
        double fmp = bad ? 0.0 : phi1PredAt(i, ind, op, neqGuard, nH2Theta, emp, bad);
        double fmm = bad ? 0.0 : phi1PredAt(i, ind, op, neqGuard, nH2Theta, emm, bad);
        if (!bad) {
          double hij = -(fpp - fpm - fmp + fmm) / (4.0 * fdH * fdH);
          H(ic, jc) += hij;
          H(jc, ic) += hij;
        }
      }
    }
  }

  // Score row i's Laplace objective from its rowPred (log-density, already
  // accumulated by phi1AnalyticHessAt/phi1FDHessAt) and eta-Hessian H
  // (accumulated in place, not yet including the Omega prior): adds the
  // Omega^-1 (IGamma2_phi1) prior term, mirrors it into the lower triangle,
  // computes log|H|, and validates it (finite, positive-definite -- H must
  // be a valid precision matrix for -2*loglik + log|H| to be a real Laplace
  // score). Returns false (row is bad) on any failure.
  bool phi1LaplaceScore(double rowPred, arma::mat &H, double &score) {
    for (int jc = 0; jc < nphi1; ++jc)
      for (int ic = 0; ic <= jc; ++ic) {
        H(ic, jc) += IGamma2_phi1(ic, jc);
        H(jc, ic) = H(ic, jc);
      }
    double logdetH, sgnH;
    if (!arma::log_det(logdetH, sgnH, H) || !(sgnH > 0) || !R_finite(logdetH)) return false;
    score = -2.0 * rowPred + logdetH;
    return true;
  }

  double phi1Objective(double *p) {
    rx_solving_options *op = getSolvingOptions(_rx);
    int cores = getOpCores(op);
    bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
    int nH2Theta = (int)_saemPhi1H2ThetaKind.n_elem;
    int nEta = (int)_saemPhi1H2EtaCol.n_elem;
    int slot = _saemPhi1UseAnalyticHess ? odeSlotHess2 : odeSlotPred;
    // A fixed relative-or-absolute step for the FD Hessian fallback -- a
    // simpler, deliberately non-adaptive alternative to calcEtaHessian's own
    // Shi(2021) stepping (src/shi21.h); see R/saemPhi1Inner.R's own docs on
    // this being a documented v1 simplification, not a load-bearing choice.
    const double fdH = 1e-4;
    std::vector<double> rowScore((size_t)nM, 0.0);
    std::vector<int> rowBad((size_t)nM, 0);
    // The event-sensitivity shape is a process global, installed/cleared only by
    // OdeSwapEsBatch, which MUST run OUTSIDE the OpenMP region below.  odeSlotHess2
    // carries its OWN shape (odeEsHess2); odeSlotPred carries none, but per
    // OdeSwapEsBatch's own contract (src/odeSwap.cpp) a "no ES" slot STILL needs a
    // batch constructed for it, so a shape left installed by some earlier solve
    // (this session's own odeSlotHess2 call, or a prior FOCEi fit's fit-wide inner
    // load) gets explicitly deactivated rather than silently reused with the wrong
    // dimensions on this solve's dosing events.  Always construct one, for `slot`.
    std::unique_ptr<OdeSwapEsBatch> phi1EsBatch(new OdeSwapEsBatch(slot));
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
    for (int i = 0; i < nM; ++i) {
#ifdef _OPENMP
      if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
      int subj = i % N;
      rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
      OdeSwapScope neqGuard(slot, ind, op);
      OdeSwapCmtScope cmtGuard(slot, op, ind);
      for (int k = 0; k < nH2Theta; ++k) {
        int kind = _saemPhi1H2ThetaKind(k);
        int col = _saemPhi1H2ThetaCol(k);
        double val = (kind == 1) ? p[col] :
          ((kind == 0) ? mprior_phi0(0, col) : _saemPhi1H2ThetaFixedVal(k));
        setIndParPtr(ind, k, val);
      }
      arma::vec eta0(nEta);
      for (int k = 0; k < nEta; ++k) {
        int col = _saemPhi1H2EtaCol(k);
        eta0(k) = phiM(i, i1(col)) - mprior_phi1(subj, col);
      }
      bool bad = false;
      double rowPred = 0.0;
      arma::mat H(nphi1, nphi1, arma::fill::zeros);
      if (!_saemPhi1WantHessian) {
        // saemControl(phi1Hessian=FALSE), the default: plain -2*loglik, no
        // Hessian/Omega-prior term at all -- see the field's own docs.
        rowPred = phi1PredAt(i, ind, op, neqGuard, nH2Theta, eta0, bad);
      } else if (_saemPhi1UseAnalyticHess) {
        if (!phi1AnalyticHessAt(i, ind, op, neqGuard, nH2Theta, eta0, rowPred, H)) {
          rowBad[i] = 1; continue;
        }
      } else {
        rowPred = phi1PredAt(i, ind, op, neqGuard, nH2Theta, eta0, bad);
        // restore of the row's solve/lhs state to eta0 is unnecessary --
        // calcEtaHessian's own FD fallback leaves the last perturbation's
        // solve behind too (its caller re-solves before any further read),
        // but SAEM's OWN model (odeSlotSaem) is what user_function reads
        // next, an independent peer/solve buffer.
        if (!bad) phi1FDHessAt(i, ind, op, neqGuard, nH2Theta, fdH, eta0, rowPred, H, bad);
      }
      if (bad) { rowBad[i] = 1; continue; }
      if (!_saemPhi1WantHessian) {
        if (!R_finite(rowPred)) { rowBad[i] = 1; continue; }
        rowScore[i] = -2.0 * rowPred;
        continue;
      }
      if (!phi1LaplaceScore(rowPred, H, rowScore[i])) { rowBad[i] = 1; continue; }
    }
    double total = 0.0;
    for (int i = 0; i < nM; ++i) {
      if (rowBad[i]) return 1e300;
      total += rowScore[i];
    }
    if (!std::isfinite(total)) return 1e300;
    return total;
  }

  // Phase 4: bobyqa-driven refinement of phi1's mu (mirrors refinePhi0Lik's
  // own .boundedResidOpt/minqa::bobyqa wiring exactly -- see the user's own
  // direction that the optimizer is bobyqa and the objective evaluation is
  // pure C++, no R calls). Intercept-only (no phi1 covariate, guaranteed by
  // _saemPhi1PoolReady): each free phi1 column's current mu (mprior_phi1(0,c),
  // identical across all N subjects) seeds bobyqa, unbounded (the Laplace
  // objective is well-behaved -- no ini() bounds apply to a mu-ref theta the
  // way phi0Lower/phi0Upper do to fixed-effect-only ones). The result is
  // damped by the SA step pas(kiter) exactly like refinePhi0Lik's own
  // mprior_phi0 update, and MCOV1 is back-solved so next iteration's
  // mprior_phi1=COV1*MCOV1 reproduces it -- the row-level phiM/eta values
  // do_mcmc already holds are left completely untouched.
  // Back-solve MCOV1's per-column covariate loadings against mprior_phi1,
  // restricted to each phi1 column's own LCOV1 design columns -- mirrors
  // phi0Objective's own MCOV0 back-solve (src/saem.cpp) for the identical
  // rank-deficiency reason: a single least-squares against all of COV1 is
  // rank deficient whenever nphi1 > 1 with no phi1 covariate (every column
  // of COV1 is then the same intercept column).
  void phi1BackSolveMCOV() {
    for (int c = 0; c < nphi1; ++c) {
      uvec li = arma::find(LCOV1.col(c) == 1);
      if (li.n_elem == 0) continue;
      mat Xc = COV1.cols(li);
      vec bc;
      if (arma::solve(bc, Xc.t() * Xc, Xc.t() * mprior_phi1.col(c))) {
        for (unsigned int j = 0; j < li.n_elem; ++j) MCOV1(li(j), c) = bc(j);
      }
    }
  }

  void refinePhi1Lik(unsigned int kiter, const vec &pas) {
    if (!_saemPhi1PoolReady || nphi1 <= 0) return;
    std::vector<bool> phi1Fix((size_t)nphi1, false);
    for (unsigned int j = 0; j < fixedIx1.n_elem; ++j) {
      if (fixedIx1(j) < (unsigned int)nphi1) phi1Fix[(size_t)fixedIx1(j)] = true;
    }
    gPhi1FreeIx.clear();
    for (int c = 0; c < nphi1; ++c) {
      if (!phi1Fix[(size_t)c]) gPhi1FreeIx.push_back(c);
    }
    if (gPhi1FreeIx.empty()) return;

    gPhi1Self = this;
    Rcpp::NumericVector par0(nphi1);
    for (int c = 0; c < nphi1; ++c) par0[c] = mprior_phi1(0, c);
    gPhi1Full.set_size(nphi1);
    for (int c = 0; c < nphi1; ++c) gPhi1Full[c] = par0[c];

    // Bounded local trust region around the current mu (mirrors
    // refinePhi0Lik's own localTrust clamp, src/saem.cpp) -- an UNBOUNDED
    // search lets bobyqa wander toward a point where the eta-Hessian
    // degenerates (log|H| -> -Inf reads as a spuriously low objective,
    // "rewarding" the wander), which then destabilizes the EM's own
    // Gamma2_phi1 update on the resulting garbage mprior_phi1 (measured:
    // an unbounded search on a single-eta TTE model eventually hit
    // "problem with matrix inverse" after enough burn-in iterations).
    const double trust = 0.75;
    int nFree = (int)gPhi1FreeIx.size();
    Rcpp::NumericVector parFree(nFree), loFree(nFree), hiFree(nFree);
    for (int fi = 0; fi < nFree; fi++) {
      int c = gPhi1FreeIx[(size_t)fi];
      parFree[fi] = par0[c];
      loFree[fi] = par0[c] - trust;
      hiFree[fi] = par0[c] + trust;
    }
    Rcpp::Environment nlmixr2 = Rcpp::Environment::namespace_env("nlmixr2est");
    Rcpp::Function boundedOpt = nlmixr2[".boundedResidOpt"];
    Rcpp::InternalFunction fn(&gPhi1ObjR);
    Rcpp::List ctl = Rcpp::List::create(Rcpp::_["maxfun"] = phi1ThetaMaxEval);
    Rcpp::List ret = boundedOpt(Rcpp::_["par"] = parFree, Rcpp::_["fn"] = fn,
                                Rcpp::_["lower"] = loFree, Rcpp::_["upper"] = hiFree,
                                Rcpp::_["control"] = ctl);
    Rcpp::NumericVector rxOpt = ret["x"];
    Rcpp::NumericVector xmin(nphi1);
    for (int c = 0; c < nphi1; ++c) xmin[c] = par0[c];
    for (int fi = 0; fi < nFree; ++fi) xmin[gPhi1FreeIx[(size_t)fi]] = rxOpt[fi];

    for (int c = 0; c < nphi1; ++c) {
      double cur = mprior_phi1(0, c);
      mprior_phi1.col(c).fill(cur + pas(kiter) * (xmin[c] - cur));
    }
    phi1BackSolveMCOV();
    _saemPhi1RefineN++;
  }

  void set_fn(user_funct f) {
    user_fn = f;
  }

  // Continuous-time AR(1) whitening of a residual/SD pair, in place.  e and g are
  // indexed in the original 1-chain observation order (length ntotal).  For each
  // AR endpoint obs with a previous same-subject-same-endpoint record:
  //   phi_i = cor^dt_i, eps_i = e_i - phi_i*e_{prev}, gstar_i = g_i*sqrt(1-phi_i^2).
  // The first record of each subject/endpoint (arPrev<0) is left marginal.  The
  // previous residual is the RAW (pre-whitening) residual, so snapshot e first.
  void arWhiten(vec &e, vec &g) {
    if (!hasAr) return;
    vec e0 = e;
    for (arma::uword i = 0; i < e.n_elem; ++i) {
      int b = (int)ix_endpnt(i);
      if (!arActive(b)) continue;
      arma::sword p = arPrev(i);
      if (p < 0) continue;
      double phi = std::pow(arCor(b), arDt(i));
      double om = 1.0 - phi*phi;
      if (om < 1e-8) om = 1e-8;
      e(i) = e0(i) - phi*e0((arma::uword)p);
      g(i) *= std::sqrt(om);
    }
  }

  // Per-observation Gaussian -LL contribution with the AR(1) whitening applied
  // (reduces to the independent 0.5*((yt-ft)/g)^2 + log(g) when no AR).  Also
  // writes the whitened (conditional) prediction/SD into _scratch_ftAr/_scratch_gAr
  // so a censored row on the SAME chain can be scored against the AR(1)
  // conditional distribution, not the marginal (ft, g) -- see #918.
  vec arDYFhyp(const vec &yt, const vec &ft, const vec &g) {
    vec e = yt - ft;
    vec gg = g;
    arWhiten(e, gg);
    _scratch_ftAr = yt - e;
    _scratch_gAr = gg;
    return 0.5*(e/gg)%(e/gg) + log(gg);
  }

  // Replace the normal per-observation loss in `DYFm` with the censored one
  // (#876).  doCensNormal1 speaks the FOCEi inner's language -- it takes and
  // returns a LOG-LIKELIHOOD, wants the VARIANCE, and reads the DV on the
  // TRANSFORMED scale -- while the SAEM chain carries the NEGATED
  // log-likelihood (a loss), the residual SD `g`, and the untransformed y.
  // Translating in both directions is what makes a censored row score the
  // same here as in likInner0: without it an M3/M4 row's censored term
  // arrives with its sign flipped (a log-likelihood stored back as a loss)
  // and its scale wrong (SD passed where a variance is wanted). An
  // uncensored row comes back untouched.
  inline void applyCensLoss(mat &DYFm, const uvec &indioK, const vec &censk,
                            const vec &ytk, const vec &limT, const vec &ft,
                            const vec &g) const {
    for (int j = ntotal; j--;) {
      DYFm(indioK(j)) = -doCensNormal1(censk[j], ytk[j], limT[j],
                                       -DYFm(indioK(j)), ft[j], g[j]*g[j], 0);
    }
  }

  // Final per-endpoint estimated AR(1) correlation (0 for non-AR endpoints).
  vec get_arCor() { return arCor; }

  // Reset the AR(1) M-step accumulators (called with the statr reset each iter).
  void arResetMstep() {
    for (int b = 0; b < nendpnt; ++b) {
      arPairR[b].clear(); arPairP[b].clear();
      arPairDt[b].clear(); arPairW[b].clear();
      arFirstSSR[b] = 0.0; arNobs[b] = 0;
    }
  }

  // AR(1) correlation M-step for endpoint b: grid-search the profiled negative
  // log-likelihood g(cor) = n*log(WSSR(cor)) + sum w*log(1-phi^2) over cor in
  // [0, ~1), then take a stochastic-approximation step toward the maximizer.
  // WSSR(cor) = firstSSR + sum w*eps(cor)^2/(1-phi^2), eps = r_i - cor^dt*r_prev.
  void arUpdateCor(int b, int kiter, const vec &pas) {
    const double double_xmin = 1.0e-200;
    size_t np = arPairR[b].size();
    if (np == 0 || arNobs[b] == 0) return;
    double best = std::numeric_limits<double>::infinity();
    double corHat = arCor(b);
    for (int gi = 0; gi <= 99; ++gi) {
      double c = gi*0.0099;
      double wssr = arFirstSSR[b];
      double logdet = 0.0;
      for (size_t j = 0; j < np; ++j) {
        double phi = std::pow(c, arPairDt[b][j]);
        double om = 1.0 - phi*phi; if (om < 1e-8) om = 1e-8;
        double eps = arPairR[b][j] - phi*arPairP[b][j];
        wssr += arPairW[b][j]*eps*eps/om;
        logdet += arPairW[b][j]*std::log(om);
      }
      if (wssr < double_xmin) wssr = double_xmin;
      double g = arNobs[b]*std::log(wssr) + logdet;
      if (g < best) { best = g; corHat = c; }
    }
    arCor(b) = arCor(b) + pas(kiter)*(corHat - arCor(b));
    if (arCor(b) < 0.0) arCor(b) = 0.0;
    if (arCor(b) > 0.999) arCor(b) = 0.999;
  }

  // One endpoint's residual SSR contribution for one MCMC chain / mixture
  // component (weightCol into mixWeights).  Independent path = original sum of
  // standardized r^2.  AR path accumulates the whitened SSR (eps^2/(1-phi^2))
  // and stores the (r, r_prev, dt, w) pairs + first-obs SSR for arUpdateCor.
  double arResk(int b, const vec &f_cur, const vec &y_cur, int weightCol) {
    const double double_xmin = 1.0e-200, xmax = 1e300;
    double resk = 0.0;
    if (arActive(b)) {
      for (int i = 0; i < (int)y_cur.size(); i++) {
        int idx_orig = ix_sorting(y_offset(b) + i);
        double ft = _powerD(f_cur[i], lambda(b), yj(b), low(b), hi(b));
        double r_ji = y_cur[i] - ft;
        if (res_mod(b) == rmProp) {
          double fci = f_cur[i];
          double fa = handleF(propT(b), ft, fci, true, true);
          if (fa <= double_xmin) fa = 1.0;
          r_ji /= fa;
        }
        _arRorig(idx_orig) = r_ji;
      }
      for (int i = 0; i < (int)y_cur.size(); i++) {
        int idx_orig = ix_sorting(y_offset(b) + i);
        int i_subj = obs_subject(idx_orig);
        double w = (weightCol < 0) ? 1.0 : mixWeights(i_subj, weightCol);
        double r_ji = _arRorig(idx_orig);
        arma::sword p = arPrev(idx_orig);
        double contrib;
        if (p < 0) {
          contrib = r_ji*r_ji;
          arFirstSSR[b] += w*contrib;
        } else {
          double phi = std::pow(arCor(b), arDt(idx_orig));
          double om = 1.0 - phi*phi; if (om < 1e-8) om = 1e-8;
          double rp = _arRorig((arma::uword)p);
          double eps = r_ji - phi*rp;
          contrib = eps*eps/om;
          arPairR[b].push_back(r_ji); arPairP[b].push_back(rp);
          arPairDt[b].push_back(arDt(idx_orig)); arPairW[b].push_back(w);
        }
        arNobs[b] += 1;
        if (contrib > xmax) contrib = xmax;
        else if (contrib < double_xmin) contrib = double_xmin;
        resk += w * contrib;
      }
    } else {
      for (int i = 0; i < (int)y_cur.size(); i++) {
        int idx_orig = ix_sorting(y_offset(b) + i);
        int i_subj = obs_subject(idx_orig);
        double ft = _powerD(f_cur[i], lambda(b), yj(b), low(b), hi(b));
        double r_ji = y_cur[i] - ft;
        if (res_mod(b) == rmProp) {
          double fci = f_cur[i];
          double fa = handleF(propT(b), ft, fci, true, true);
          if (fa <= double_xmin) fa = 1.0;
          r_ji /= fa;
        }
        double r_ji_sq = r_ji * r_ji;
        if (r_ji_sq > xmax) r_ji_sq = xmax;
        else if (r_ji_sq < double_xmin) r_ji_sq = double_xmin;
        double w = (weightCol < 0) ? 1.0 : mixWeights(i_subj, weightCol);
        resk += w * r_ji_sq;
      }
    }
    return resk;
  }

  // Data augmentation for the M-step residual SSR (see #916): return a copy of
  // this chain/endpoint's transformed observation vector with every censored
  // (M3/M4) row replaced by a draw from the truncated normal implied by this
  // chain's prediction (simCensDv()), so arResk()/the direct SSR see a value
  // consistent with the censoring instead of the recorded LOQ/limit.
  // cens_cur/limit_cur are RAW (untransformed), same length/order as y_cur/
  // f_cur (chain-sliced, ix_sorting-applied, endpoint span).
  //
  // y_cur is on whatever scale the caller's hasFixedObsTransform branch put
  // it on (ysTrans, i.e. already transformed, when the TBS transform is
  // fixed; the raw ys when it is estimated -- see arResk()'s own read of
  // y_cur a few lines below every call site).  simCensDv() itself works
  // entirely on the TRANSFORMED scale (that is what f/sd/the truncation
  // bounds are in), so both the limDv fed to it and the value it returns are
  // converted between that scale and y_cur's ambient one via _powerD()/
  // _powerDi() -- a no-op when hasFixedObsTransform is true.  This keeps a
  // censored row's simulated replacement on the SAME scale as its
  // uncensored neighbors in the returned vector; it does not touch the
  // separate, pre-existing mismatch between y_cur and f_cur that arResk()
  // itself has for EVERY row (censored or not) when the transform is
  // estimated, which traces to #914 (the saem lambda member never actually
  // updates from its initial value) and is out of scope here.
  vec augmentCensY(int b, const vec &f_cur, const vec &y_cur,
                    const vec &cens_cur, const vec &limit_cur,
                    int kiter, int k, int mixIdx) {
    if (!arma::any(cens_cur != 0.0)) return y_cur;
    vec y_aug = y_cur;
    _saemSeedCensAug((uint32_t)saemSeed, kiter, k, b, mixIdx);
    const double double_xmin = 1.0e-200, xmax = 1e300;
    for (unsigned int i = 0; i < y_aug.n_elem; i++) {
      if (cens_cur[i] == 0.0) continue;
      double fci = f_cur[i];
      double ft = _powerD(fci, lambda(b), yj(b), low(b), hi(b));
      double ftT = handleF(propT(b), ft, fci, false, true);
      double sd = ares(b) + bres(b) * std::fabs(ftT);
      if (sd == 0.0) sd = 1.0;
      else if (sd < double_xmin) sd = double_xmin;
      else if (sd > xmax) sd = xmax;
      double limT = _powerD(limit_cur[i], lambda(b), yj(b), low(b), hi(b));
      double limDvT = hasFixedObsTransform ? y_cur[i] :
        _powerD(y_cur[i], lambda(b), yj(b), low(b), hi(b));
      double simT = simCensDv(cens_cur[i], limDvT, limT, ft, sd);
      y_aug[i] = hasFixedObsTransform ? simT :
        _powerDi(simT, lambda(b), yj(b), low(b), hi(b));
    }
    setRxThreadId(-1);
    return y_aug;
  }

  // Fill this chain's per-endpoint residual log-sigma2 score/Hessian entries
  // from resy(b,k) (endpoint b's residual SSR for MCMC chain k).  Only a pure
  // additive endpoint (res_mod==rmAdd) has a valid single-parameter
  // log-sigma2 score; every other endpoint's slot is held at exactly 0 (see
  // the nb_param comment in inits()).  d1_logsigma2 must already be sized
  // nResidEp; d2logk is nb_param x nb_param.
  void fillResidLogSigma2(int k, const mat &resy, vec &d1_logsigma2, mat &d2logk) {
    int resBase = nlambda + nphi1;
    for (int b = 0; b < nendpnt; b++) {
      int idx = residEpIdx[b];
      if (idx < 0) continue;
      int col = resBase + idx;
      if (res_mod(b) == rmAdd) {
        double nb = (double)(y_offset(b + 1) - y_offset(b));
        d1_logsigma2[idx] = 0.5 * resy(b, k) / sigma2[b] - 0.5 * nb;
        d2logk(col, col) = -0.5 * resy(b, k) / sigma2[b];
      } else {
        d1_logsigma2[idx] = 0.0;
        d2logk(col, col) = 0.0;
      }
    }
  }

  mat get_resMat() {
    mat m(nendpnt,4);
    m.col(0) = ares;
    m.col(1) = bres;
    m.col(2) = cres;
    m.col(3) = lres;
    return m;
  }

  mat get_trans() {
    mat m(nendpnt, 4);
    m.col(0) = lambda;
    m.col(1) = conv_to<vec>::from(yj); // convert uvec to double mat
    m.col(2) = low;
    m.col(3) = hi;
    return m;
  }

  mat get_mprior_phi() {
    mat m = mpost_phi;
    m.cols(i1) = mprior_phi1;
    return m;
  }

  mat get_mpost_phi() {
    return mpost_phi;
  }

  mat get_Plambda() {
    return Plambda;
  }

  mat get_Gamma2_phi1() {
    return Gamma2_phi1;
  }

  // Reporting-only pooled BSV for split ETAs; falls back to the live matrix if no pooling was ever applied.
  mat get_Gamma2_phi1Report() {
    if (Gamma2_phi1Report.n_elem == Gamma2_phi1.n_elem) return Gamma2_phi1Report;
    return Gamma2_phi1;
  }

  mat get_Ha() {
    return Ha;
  }

  vec get_sig2() {
    return vcsig2;                                       //FIXME: regression due to multiple endpnts?
  }

  List get_resInfo() {
    vec sig2(bres.size());
    std::copy(sigma2, sigma2+bres.size(), &sig2[0]);
    return List::create(_["sigma2"]  = wrap(sig2),
			_["ares"]    = wrap(ares),
			_["bres"]    = wrap(bres),
			_["cres"]    = wrap(cres),
			_["lres"]    = wrap(lres),
			_["res_mod"] = wrap(res_mod));
  }

  mat get_par_hist() {
    return par_hist;
  }

  IntegerVector get_ueRevisitInfo() {
    return IntegerVector::create(_["ran"] = ueRevisitRan,
                                 _["unfroze"] = ueRevisitUnfroze,
                                 _["froze"] = ueRevisitFroze);
  }
  mat get_HaSa() {
    return HaSa;
  }

  vec get_mixProb() {
    return mixProb;
  }

  mat get_mixWeights() {
    return mixWeights;
  }

  mat get_eta() {
    mat eta = mpost_phi.cols(i1);
    eta -= mprior_phi1;
    mat ue = current_saem_state->_saemUE.rows(0, eta.n_rows - 1);
    ue = ue.cols(i1);
    eta = eta % ue;
    return eta;
  }

  void inits(List x) {
    _saemItmax = as<int>(x["itmax"]);
    _saemTol = as<double>(x["tol"]);
    _saemType = as<int>(x["type"]);
    _saemLambdaR = fabs(as<double>(x["lambdaRange"]));
    _saemPowR = fabs(as<double>(x["powRange"]));
    current_saem_state->_saemIncreaseTol=0;
    current_saem_state->_saemIncreasedTol2=0;
    current_saem_state->_saemMaxOdeRecalc = abs(as<int>(x["maxOdeRecalc"]));
    current_saem_state->_saemOdeRecalcFactor = fabs(as<double>(x["odeRecalcFactor"]));
    current_saem_state->_saemIndTolRelax = as<bool>(x["indTolRelax"]);
    current_saem_state->_saemUE = as<mat>(x["ue"]);

    nmc = as<int>(x["nmc"]);
    nu = as<uvec>(x["nu"]);
    niter = as<int>(x["niter"]);
    saemSeed = x.containsElementNamed("seed") ? as<int>(x["seed"]) : 99;
    nb_correl = as<int>(x["nb_correl"]);
    nb_fixOmega = as<int>(x["nb_fixOmega"]);
    nb_fixResid = as<int>(x["nb_fixResid"]);
    resValue = as<vec>(x["resValue"]);
    resFixed = as<uvec>(x["resFixed"]);
    resKeep = find(resFixed==0);
    niter_phi0 = as<int>(x["niter_phi0"]);
    coef_phi0 = as<double>(x["coef_phi0"]);
    // ini-block bounds of the phi0 thetas (i0 column order) for the bounded
    // general-likelihood phi0 optimization.
    if (x.containsElementNamed("phi0Lower")) {
      phi0Lower = as<vec>(x["phi0Lower"]);
      phi0Upper = as<vec>(x["phi0Upper"]);
    }
    nb_sa = as<int>(x["nb_sa"]);
    // Phase-1 (SA/burn) iteration count for the print's SA/EM row tag; -1
    // (missing, e.g. a saved cfg from an older version) disables the tag.
    nPhase1 = x.containsElementNamed("nPhase1") ? as<int>(x["nPhase1"]) : -1;
    // uninformative-eta revisit (absent in a cfg saved by an older version -> off)
    ueRevisitIter = x.containsElementNamed("ueRevisitIter") ? as<int>(x["ueRevisitIter"]) : -1;
    if (x.containsElementNamed("ueRevisitCols")) ueRevisitCols = as<uvec>(x["ueRevisitCols"]);
    if (x.containsElementNamed("ueDelta")) ueDelta = as<vec>(x["ueDelta"]);
    if (x.containsElementNamed("ueTol")) ueTol = as<double>(x["ueTol"]);
    coef_sa = as<double>(x["coef_sa"]);
    rmcmc = as<double>(x["rmcmc"]);
    pas = as<vec>(x["pas"]);
    pash = as<vec>(x["pash"]);
    // SA (stochastic-approximation) covariance phase: after the niter estimation
    // iterations, run nSaCov extra iterations with the gain frozen at zero so the
    // parameters stay at the converged estimate (theta_hat).  Only the MCMC E-step
    // resimulates phi ~ p(phi|y,theta_hat); the per-iteration Louis observed-information
    // integrand (DDa) is Monte-Carlo averaged into HaSa, giving a converged Fisher
    // information decoupled from the cooling schedule (cf. Monolix "stochastic
    // approximation" standard errors; Kuhn & Lavielle 2005).  nSaCov==0 is a no-op.
    nSaCov = x.containsElementNamed("nSaCov") ? as<int>(x["nSaCov"]) : 0;
    if (nSaCov > 0) {
      pas  = join_cols(pas,  zeros<vec>(nSaCov));   // freeze theta during the cov phase
      pash = join_cols(pash, zeros<vec>(nSaCov));
    }
    minv = as<vec>(x["minv"]);

    N = as<int>(x["N"]);
    ntotal = as<int>(x["ntotal"]);
    mlen = as<int>(x["mlen"]);
    y  = as<vec>(x["y"]);
    evt  = as<mat>(x["evt"]);
    phiM = as<mat>(x["phiM"]);
    indio = as<uvec>(x["indio"]);
    nM = N*nmc;

    nMix = x.containsElementNamed("nMix") ? as<int>(x["nMix"]) : 1;
    if (nMix > 1) {
      mixProb = as<vec>(x["mixProb"]);
      mixProbInit = mixProb;
      mixWeights.zeros(N, nMix);
      priorWeights.zeros(N, nMix);
      mixProbMethod = x.containsElementNamed("mixProbMethod") ? as<int>(x["mixProbMethod"]) : 0;
      pasMix = x.containsElementNamed("pasMix") ? as<vec>(x["pasMix"]) : pas;
      mixProbPriorN = x.containsElementNamed("mixProbPriorN") ? as<double>(x["mixProbPriorN"]) : 0.0;
      mixSampleMethod = x.containsElementNamed("mixSampleMethod") ? as<int>(x["mixSampleMethod"]) : 0;
    } else {
      mixProb.clear();
    }

    opt = as<List>(x["opt"]);                                      //CHECKME
    optM = as<List>(x["optM"]);                                    //CHECKME

    pc1 = as<uvec>(x["pc1"]);
    covstruct1 = as<mat>(x["covstruct1"]);
    Mcovariables = as<mat>(x["Mcovariables"]);

    nphi1 = as<int>(x["nphi1"]);
    i1 = as<uvec>(x["i1"]);
    Gamma2_phi1 = as<mat>(x["Gamma2_phi1"]);
    Gamma2_phi1Init = Gamma2_phi1;
    Gamma2_phi1fixedIxIn = as<umat>(x["Gamma2_phi1fixedIx"]);
    Gamma2_phi1fixedIx = find(Gamma2_phi1fixedIxIn);
    Gamma2_phi1fixed = as<int>(x["Gamma2_phi1fixed"]);
    if (Gamma2_phi1fixed==1) {
      Gamma2_phi1fixedValues = as<mat>(x["Gamma2_phi1fixedValues"]);
    }
    mprior_phi1 = as<mat>(x["mprior_phi1"]);
    COV1 = as<mat>(x["COV1"]);
    LCOV1 = as<mat>(x["LCOV1"]);
    COV21 = as<mat>(x["COV21"]);
    MCOV1 = as<mat>(x["MCOV1"]);
    jcov1 = as<uvec>(x["jcov1"]);
    ind_cov1 = as<uvec>(x["ind_cov1"]);
    statphi11 = as<mat>(x["statphi11"]);
    statphi12 = as<mat>(x["statphi12"]);
    omegaShare = x.containsElementNamed("omegaShare") ? as<uvec>(x["omegaShare"]) : uvec();
    omegaShareSubpop = x.containsElementNamed("omegaShareSubpop") ? as<uvec>(x["omegaShareSubpop"]) : uvec();
    omegaPool = x.containsElementNamed("omegaPool") ? as<uvec>(x["omegaPool"]) : uvec();
    statphi11_mix.set_size(std::max(nMix, 1));
    for (int _j = 0; _j < std::max(nMix, 1); _j++) {
      statphi11_mix(_j) = statphi11;
    }

    nphi0 = as<int>(x["nphi0"]);
    if (nphi0>0) {
      i0 = as<uvec>(x["i0"]);
      Gamma2_phi0 = as<mat>(x["Gamma2_phi0"]);
      // C_i mu_{k+1} = mprior_phi0
      mprior_phi0 = as<mat>(x["mprior_phi0"]);
      COV0 = as<mat>(x["COV0"]);
      LCOV0 = as<mat>(x["LCOV0"]);
      COV20 = as<mat>(x["COV20"]);
      MCOV0 = as<mat>(x["MCOV0"]);
      jcov0 = as<uvec>(x["jcov0"]);
      ind_cov0 = as<uvec>(x["ind_cov0"]);
      statphi01 = as<mat>(x["statphi01"]);
      statphi02 = as<mat>(x["statphi02"]);
    }
    fixedIx0 = as<uvec>(x["fixed.i0"]);
    fixedIx1 = as<uvec>(x["fixed.i1"]);

    nlambda1 = as<int>(x["nlambda1"]);
    nlambda0 = as<int>(x["nlambda0"]);
    nlambda = nlambda1 + nlambda0;
    nphi = nphi1+nphi0;
    Plambda.zeros(nlambda);
    ilambda1 = as<uvec>(x["ilambda1"]);
    ilambda0 = as<uvec>(x["ilambda0"]);

    DYF = zeros<mat>(mlen, nM);
    phi.set_size(N, nphi, nmc);

    //FIXME
    nendpnt=as<int>(x["nendpnt"]);
    distribution=as<int>(x["distribution"]);
    // One FIM residual slot per endpoint that carries a residual parameter --
    // none when the whole model is a general log-likelihood (distribution==4;
    // "any LL endpoint" forces the WHOLE model to distribution==4, so no
    // endpoint has a real residual in that case).  Sizing to nendpnt (not
    // nres, the total residual PARAMETER count) is deliberate: the analytic
    // Louis FIM only ever tracks a single log-sigma2 score/Hessian per
    // endpoint, valid only for a pure additive residual (res_mod==rmAdd);
    // any other endpoint's slot exists (so nb_param has a fixed layout) but
    // its d1_logsigma2/d2logk entries are held at exactly 0 in every
    // iteration, so its row/col of Ha/HaSa stays exactly 0 and
    // .saemFimToCov (R/saem.R) can drop it and fall back to the linFim
    // splice for that endpoint's residual SE.
    nResidEp = (distribution == 4) ? 0 : nendpnt;
    for (int b = 0; b < MAXENDPNT; ++b) residEpIdx[b] = -1;
    for (int b = 0; b < nResidEp; ++b) residEpIdx[b] = b;
    nb_param = nphi1 + nlambda + nResidEp;
    ix_sorting=as<uvec>(x["ix_sorting"]);
    ys = y(ix_sorting);    //ys: obs sorted by endpnt
    y_offset=as<uvec>(x["y_offset"]);
    res_mod = as<uvec>(x["res.mod"]);
    ares = as<vec>(x["ares"]);
    bres = as<vec>(x["bres"]);
    cres = as<vec>(x["cres"]);
    lres = as<vec>(x["lres"]);
    yj = as<uvec>(x["yj"]);
    propT=as<uvec>(x["propT"]);
    // lambda mirrors lres (the M-step's working boxCox/yeoJohnson estimate);
    // seed it from lres, not x["lambda"] (which the R side always ships as 1),
    // and keep it synced wherever lres is updated below (#914).
    lambda = lres;
    low = as<vec>(x["low"]);
    hi = as<vec>(x["hi"]);

    ix_endpnt=as<uvec>(x["ix_endpnt"]);
    ix_idM=as<umat>(x["ix_idM"]);
    res_offset=as<uvec>(x["res_offset"]);
    addProp=as<uvec>(x["addProp"]);
    arCor=as<vec>(x["arCor"]);
    arActive=as<uvec>(x["arActive"]);
    arPrev=as<arma::ivec>(x["arPrev"]);
    arDt=as<vec>(x["arDt"]);
    hasAr = (int)accu(arActive);
    hasFixedObsTransform = true;
    for (unsigned int b = 0; b < res_mod.n_elem; ++b) {
      if (res_mod[b] >= rmAddLam && res_mod[b] <= rmAddPowLam) {
        hasFixedObsTransform = false;
        break;
      }
    }
    if (hasFixedObsTransform) {
      // Compute yTrans for N subjects only (not repeated nmc times)
      yTrans = y;
      for (unsigned int i = 0; i < yTrans.n_elem; ++i) {
        int cur = ix_endpnt(i);
        yTrans[i] = _powerD(y[i], lambda(cur), yj(cur), low(cur), hi(cur));
      }
      ysTrans = ys;
      for (int b = 0; b < nendpnt; ++b) {
        for (unsigned int i = y_offset(b); i < y_offset(b + 1); ++i) {
          ysTrans[i] = _powerD(ys[i], lambda(b), yj(b), low(b), hi(b));
        }
      }
    }
    nres = res_offset.max();
    vcsig2.set_size(nres);
    vecares = ares(ix_endpnt);
    vecbres = bres(ix_endpnt);
    veccres = cres(ix_endpnt);
    veclres = lres(ix_endpnt);
    vecaddProp = addProp(ix_endpnt);
    // Pre-allocate per-chain scratch buffers for the distribution==1 hot loops
    _scratch_ft.set_size(ntotal);
    _scratch_limitT.set_size(ntotal);
    _scratch_ftT.set_size(ntotal);
    _scratch_g.set_size(ntotal);
    _scratch_ftAr.set_size(ntotal);
    _scratch_gAr.set_size(ntotal);
    _scratch_indio = indio;  // same length as indio, initialise from it
    _arRorig.set_size(ntotal);
    for (int b=0; b<nendpnt; ++b) {
      sigma2[b] = 10;
      if (res_mod(b) == rmAdd) {
        sigma2[b] = max(ares(b)*ares(b), 10.0);
      }
      if (res_mod(b) == rmProp) {
        sigma2[b] = max(bres(b)*bres(b), 1.0);
      }
      statrese[b] = 0.0;
    }

    par_hist = as<mat>(x["par.hist"]);
    parHistThetaKeep=as<uvec>(x["parHistThetaKeep"]);
    parHistThetaKeep = find(parHistThetaKeep);
    parHistOmegaKeep=as<uvec>(x["parHistOmegaKeep"]);
    parHistOmegaKeep = find(parHistOmegaKeep);
    // off-diagonal Omega covariances recorded in the iteration history: 0-indexed
    // (row, col) pairs into Gamma2_phi1, appended after the diagonal variances
    parHistOmegaOffPairs = x.containsElementNamed("parHistOmegaOffPairs") ?
      as<umat>(x["parHistOmegaOffPairs"]) : umat(0, 2);

    obs_subject.set_size(ntotal);
    for (int i = 0; i < N; i++) {
      int start = ix_idM(i, 0);
      int end = ix_idM(i, 1);
      for (int idx = start; idx <= end; idx++) {
        obs_subject(idx) = i;
      }
    }

    // Set up the shared scale.h iteration-print struct (scaleTypeNone; xform
    // sub-list wired via scaleAttachXform, same path as other estimators).
    scaleNames = as<CharacterVector>(x["parHistNames"]);
    // Off-diagonal Omega covariances add parHistOmegaOffPairs.n_rows rows;
    // mixture models add (nMix - 1) mixture-weight rows.
    int nprint = parHistThetaKeep.n_elem + parHistOmegaKeep.n_elem + parHistOmegaOffPairs.n_rows + resKeep.n_elem + (nMix > 1 ? nMix - 1 : 0);
    scaleInitPar.assign(std::max(nprint, 1), 0.0);
    scaleC.assign(std::max(nprint, 1), NA_REAL);
    scaleSetup(&scale,
               scaleInitPar.data(),
               scaleC.data(),
               scaleNames,
               /*useColor*/0, /*printNcol*/1, /*print*/0,
               normTypeConstant,
               scaleTypeNone,
               1e-7, 1e7, 0.0,
               nprint);
    scaleAttachXform(&scale, as<List>(x["xform"]));
    scaleApplyIterPrintControl(&scale, as<List>(x["iterPrintControl"]));
    // saem has no per-iteration objective function; suppress the Function
    // Val column entirely so users don't see "nan" in every iteration row.
    scale.showOfv = 0;
    scale.save = 0; // par_hist already records the iteration history
    if (nPhase1 >= 0) {
      scale.keyExtra = "SA: Stochastic-approximation (burn-in) phase; EM: EM phase\n";
    }

    L  = zeros<vec>(nb_param);
    Ha = zeros<mat>(nb_param,nb_param);
    Hb = zeros<mat>(nb_param,nb_param);
    mpost_phi = zeros<mat>(N, nphi);
    cpost_phi = zeros<mat>(N, nphi);

    //handle situation when nphi0=0
    mprior_phi0.set_size(N, nphi0);
    statphi01.set_size(N, nphi0);

    mx.nM     = nM;
    mx.y      = y;
    mx.indio  = indio;
    mx.evtM   = evt;
    mx.optM   = optM;

    nonMuThetaRegress = x.containsElementNamed("nonMuThetaRegress") ?
      as<int>(x["nonMuThetaRegress"]) : 0;
    nonMuThetaOptType = x.containsElementNamed("nonMuThetaOptType") ?
      as<int>(x["nonMuThetaOptType"]) : 0;
    nonMuThetaMaxEval = x.containsElementNamed("nonMuThetaMaxEval") ?
      as<int>(x["nonMuThetaMaxEval"]) : 25;
    nonMuThetaSweeps = x.containsElementNamed("nonMuThetaSweeps") ?
      as<int>(x["nonMuThetaSweeps"]) : 2;
    nonMuThetaEvery = x.containsElementNamed("nonMuThetaEvery") ?
      as<int>(x["nonMuThetaEvery"]) : 1;
    if (nonMuThetaEvery < 1) nonMuThetaEvery = 1;
    nonMuThetaTol = x.containsElementNamed("nonMuThetaTol") ?
      as<double>(x["nonMuThetaTol"]) : 1.0e-4;
    phi1ThetaEvery = x.containsElementNamed("phi1ThetaEvery") ?
      as<int>(x["phi1ThetaEvery"]) : 1;
    if (phi1ThetaEvery < 1) phi1ThetaEvery = 1;
    phi1ThetaMaxEval = x.containsElementNamed("phi1ThetaMaxEval") ?
      as<int>(x["phi1ThetaMaxEval"]) : 50;
    // Phase 4 (SAEM general-likelihood theta plan): THETA[k]/ETA[k] -> phi
    // column maps, present only for a general-lik fit whose
    // saemPhi1TargetMap resolved (R/saemPhi1Inner.R).  Two-tier fallback
    // discipline, matching calcEtaHessian's own (src/inner.cpp): the exact
    // analytic Hessian (innerHess2/odeSlotHess2) when it built, else a
    // finite-difference Hessian over the bare predNoLhs/odeSlotPred model
    // (e.g. linCmt(), where innerHess2 is always NULL) -- _saemPhi1UseAnalyticHess
    // picks which.  _saemPhi1PoolReady requires only the pred offset (always
    // present); a fit whose map did not resolve at all keeps the historic
    // SA-recursion for phi1, unchanged.
    //
    // _saemPhi1PoolReady/_saemPhi1UseAnalyticHess are process-wide globals
    // that outlive a fit (mirrors why _saemPhi1PoolActive/odeSwapClearAll
    // needed the same discipline in setupRx): a LATER fit whose own map did
    // NOT resolve (or is not general-lik at all) never enters the branch
    // below, so without this reset it silently inherits a PRIOR pooled
    // fit's _saemPhi1PoolReady=true and stale THETA/ETA maps sized for a
    // DIFFERENT model -- refinePhi1Lik/phi1Objective then run for a plain
    // normal fit using garbage offsets (measured: segfaulted inside
    // gPhi1ObjR, confirmed via valgrind's backtrace).
    _saemPhi1PoolReady = false;
    _saemPhi1UseAnalyticHess = false;
    _saemPhi1WantHessian = x.containsElementNamed("phi1Hessian") &&
      as<int>(x["phi1Hessian"]) != 0;
    if (opt.containsElementNamed("saemPhi1ThetaKind")) {
      _saemPhi1H2ThetaKind = as<ivec>(opt["saemPhi1ThetaKind"]);
      _saemPhi1H2ThetaCol = as<ivec>(opt["saemPhi1ThetaCol"]);
      _saemPhi1H2ThetaFixedVal = as<vec>(opt["saemPhi1ThetaFixedVal"]);
      _saemPhi1H2EtaCol = as<ivec>(opt["saemPhi1EtaCol"]);
      _saemPhi1DvCol = opt.containsElementNamed("saemPhi1DvCol") ?
        as<int>(opt["saemPhi1DvCol"]) : -1;
      _saemPhi1DvColHess2 = opt.containsElementNamed("saemPhi1DvColHess2") ?
        as<int>(opt["saemPhi1DvColHess2"]) : -1;
      _saemPhi1PredOffset = odeSwapLhsIndex(odeSlotPred, "rx_pred_");
      bool haveHess2 = odeSwapLoaded(odeSlotHess2);
      _saemPhi1H2PredOffset = haveHess2 ? odeSwapLhsIndex(odeSlotHess2, "rx_pred_") : -1;
      _saemPhi1H2HessOffset = haveHess2 ? odeSwapLhsIndex(odeSlotHess2, "rx__d2pred_1_1__") : -1;
      _saemPhi1UseAnalyticHess = haveHess2 && _saemPhi1H2PredOffset >= 0 &&
        _saemPhi1H2HessOffset >= 0 && _saemPhi1DvColHess2 >= 0;
      _saemPhi1PoolReady = _saemPhi1PredOffset >= 0 && _saemPhi1DvCol >= 0 &&
        _saemPhi1H2EtaCol.n_elem == (unsigned int)nphi1 &&
        (_saemPhi1UseAnalyticHess || odeSwapLoaded(odeSlotPred));
      // user_function (a free function) needs i0/i1 too, to map
      // _saemPhi1H2ThetaCol/_saemPhi1H2EtaCol (0-based within phi0's/phi1's
      // own subset) into _phi's nphi-wide column space.
      _saemPhi1I0 = i0;
      _saemPhi1I1 = i1;
    }
    residWarmStart = x.containsElementNamed("residWarmStart") ?
      as<int>(x["residWarmStart"]) : 1;
    mixProbRegress = x.containsElementNamed("mixProbRegress") ?
      as<int>(x["mixProbRegress"]) : 0;
    // Mixtures: the direct phi0 optimizer does not partition a per-component
    // structural theta (tcl1/tcl2) by subject membership, so it would leave an
    // under-populated component's theta unconstrained (runaway).  Fall back to
    // the stochastic phi0 block, which respects membership via the per-chain
    // mixest regressor.  (nMix is read earlier in inits.)  Likewise the residual
    // warm-start (formed at the population eta=0 prediction) is especially
    // unreliable for a mixture -- the poor initial fit inflates the residual and
    // flattens the likelihood, preventing the components from separating.  Both
    // gates must follow the reads above so they are not overwritten.
    if (nMix > 1) { nonMuThetaRegress = 0; residWarmStart = 0; }
    // (mixProbMethod="regress" degeneracy fallback is decided in saem_fit AFTER
    // the initial hard classification -- see below.)
    DEBUG=as<int>(x["DEBUG"]);
    phiMFile=as<std::vector< std::string > >(x["phiMFile"]);
    //Rcout << phiMFile[0];

  }

  // Warm-start the residual-error parameters (ares/bres) from the observed
  // per-endpoint moments at the initial predictions, like npag's npResidMoments:
  // additive SD = sqrt(mean(err^2)), proportional SD = sqrt(mean((err/f)^2)),
  // on the transform-both-sides scale.  Only seeds the scale(s) an endpoint's
  // error model actually uses; a combined (add+prop) endpoint splits the
  // variance half/half as a warm start.  fsaveFull is chain-major f (>= ntotal).
  void warmStartResid(const vec &fsaveFull) {
    if (!residWarmStart || distribution == 4) return;
    if ((int)fsaveFull.n_elem < ntotal) return;
    vec fk = fsaveFull.subvec(0, ntotal - 1);
    fk = fk(ix_sorting);
    for (int b = 0; b < nendpnt; b++) {
      int rm = (int)res_mod(b);
      bool hasAdd = (rm==rmAdd || rm==rmAddProp || rm==rmAddPow ||
                     rm==rmAddLam || rm==rmAddPropLam || rm==rmAddPowLam);
      bool hasProp = (rm==rmProp || rm==rmAddProp || rm==rmPropLam || rm==rmAddPropLam);
      // Max |ft| in this endpoint.  SAEM computes this at the population (eta=0)
      // prediction, so at small |ft| the residual is dominated by between-subject
      // variability, not residual noise -- exclude those (< 5% of max) from the
      // proportional moment so BSV does not inflate the warm-started prop SD.
      double fmax = 0.0;
      for (unsigned int i = y_offset(b); i < y_offset(b+1); i++) {
        double aft = std::fabs(_powerD(fk(i), lambda(b), yj(b), low(b), hi(b)));
        if (std::isfinite(aft) && aft > fmax) fmax = aft;
      }
      double fthr = 0.05 * fmax;
      double sAdd = 0.0, sProp = 0.0; int n = 0, nProp = 0;
      for (unsigned int i = y_offset(b); i < y_offset(b+1); i++) {
        double ft = _powerD(fk(i), lambda(b), yj(b), low(b), hi(b));
        double yt = hasFixedObsTransform ? ysTrans(i)
          : _powerD(ys(i), lambda(b), yj(b), low(b), hi(b));
        double err = ft - yt;
        if (!std::isfinite(err)) continue;
        sAdd += err * err; n++;
        if (std::fabs(ft) <= fthr) continue;
        // proportional guard for f==0: denominator 1 when |f| is tiny so a
        // near-zero prediction does not blow up the proportional moment.
        double denom = (std::fabs(ft) <= 1e-6) ? 1.0 : ft;
        double ratio = err / denom;
        if (std::isfinite(ratio)) { sProp += ratio * ratio; nProp++; }
      }
      if (n == 0) continue;
      if (nProp == 0) nProp = 1;
      double split = (hasAdd && hasProp) ? std::sqrt(0.5) : 1.0;
      double mAdd = std::sqrt(sAdd / n) * split;
      double mProp = std::sqrt(sProp / nProp) * split;
      // SAEM computes this moment at the UNCONVERGED initial phi (npag/npb use
      // converged posterior etas), so it is contaminated by structural misfit
      // and can be far too large -- clamp the warm-started value to [0.2x, 5x]
      // the user's ini value so a poor initial fit cannot push the residual into
      // a runaway basin (an over-large residual flattens the S-step likelihood,
      // which then keeps the residual large).  The f==0 guard above is the
      // shared fix; this clamp is SAEM-specific.
      if (hasAdd && std::isfinite(mAdd) && mAdd > 0.0 && ares(b) > 0.0)
        ares(b) = std::min(std::max(mAdd, 0.2*ares(b)), 5.0*ares(b));
      if (hasProp && std::isfinite(mProp) && mProp > 0.0 && bres(b) > 0.0)
        bres(b) = std::min(std::max(mProp, 0.2*bres(b)), 5.0*bres(b));
    }
    // keep the per-observation SD vectors (used in the S-step g = vecares +
    // vecbres*|ft|) consistent with the warm-started scalars.
    vecares = ares(ix_endpnt);
    vecbres = bres(ix_endpnt);
  }

  void saem_fit() {
    // (arma_rng::set_seed removed -- arma::randn uses R's RNG under
    // RcppArmadillo, seeded by set.seed(seed); the explicit arma seed was a
    // no-op that is under test as a determinism-cycle suspect.)
    double double_xmin = 1.0e-200; //FIXME hard-coded xmin, also in neldermean.hpp
    double xmax = 1e300;
    ofstream phiFile;
    _warnAtolRtol = false;
    phiFile.open(phiMFile[0].c_str());

    if (DEBUG>0) {
      RSprintf("initialization successful\n");
    }
    // Emit the column header once; periodic re-emits are handled in scalePrintFun.
    scalePrintHeader(&scale);
    if (nMix > 1) {
      phiM_mix.set_size(nMix);
      fsave_mix.set_size(nMix);
      cens_mix.set_size(nMix);
      limit_mix.set_size(nMix);
      for (int jMix = 0; jMix < nMix; jMix++) {
        phiM_mix(jMix) = phiM;
        current_saem_state->_saemMixest = jMix + 1;
        mat initMat = user_fn(phiM, evt, optM);
        fsave_mix(jMix) = initMat.col(0);
        cens_mix(jMix) = initMat.col(1);
        limit_mix(jMix) = initMat.col(2);
      }
      current_saem_state->_saemMixest = 0;
    } else {
      fsaveMat = user_fn(phiM, evt, optM);
      limit = fsaveMat.col(2);
      cens = fsaveMat.col(1);
      fsave = fsaveMat.col(0);
    }
    if (DEBUG>0){
      RSprintf("initial user_fn successful\n");
    }
    warmStartResid(nMix > 1 ? fsave_mix(0) : fsave);
    if (nMix > 1 && mixProbRegress) {
      // mixProbMethod="regress": classify each subject to its best component
      // once (hard) and hold membership fixed -- mixWeights becomes a 0/1
      // indicator, so the existing responsibility-weighted machinery (arResk,
      // sufficient stats, phiM_weighted) collapses to a hard assignment.  The
      // soft-EM E-step and the mixProb SA update are skipped below.
      mixFixedAssign = mixNaiveClassify(0.0);
      // Fall back to soft-EM (regularized) when fixed membership cannot work:
      //  (a) SPLIT-ETA mixtures -- each component owns a distinct eta, so
      //      omegaShareSubpop has >= 2 distinct non-zero subpop values (a
      //      shared-eta mixture has just one).  Components start identical and
      //      must differentiate during the fit; a hard split at the symmetric
      //      init is arbitrary and never separates.
      //  (b) a degenerate classification that leaves a component empty (its
      //      theta would be unconstrained).
      uvec _nzSub = omegaShareSubpop.elem(find(omegaShareSubpop > 0));
      uvec _uniqSub = unique(_nzSub);
      bool isSplitEta = (_uniqSub.n_elem >= 2);
      bool allPop = true;
      for (int j = 0; j < nMix; j++) {
        int cnt = 0;
        for (int i = 0; i < N; i++) if ((int)mixFixedAssign(i) == j + 1) cnt++;
        if (cnt == 0) { allPop = false; break; }
      }
      if (isSplitEta || !allPop) {
        mixProbRegress = 0;
        mixProbMethod = 1; // regularized
      } else {
      mixWeights.zeros(N, nMix);
      for (int i = 0; i < N; i++) {
        unsigned int a = mixFixedAssign(i);
        if (a >= 1 && a <= (unsigned int)nMix) mixWeights(i, a - 1) = 1.0;
      }
      mixProb = mean(mixWeights, 0).t();
      // Supply the fixed assignment as a per-subject mixest regressor so the
      // S-step can solve each subject once under its OWN component (mixest=0 in
      // user_function -> per-subject indMixest), instead of nMix per-component
      // chains.  Tile over the nmc MCMC chains (phiM row k*N+i is subject i).
      Rcpp::IntegerVector mixestFull(N * nmc);
      for (int k = 0; k < nmc; k++)
        for (int i = 0; i < N; i++)
          mixestFull[k * N + i] = (int)mixFixedAssign(i);
      mx.optM["mixest"] = mixestFull;
      } // end allPop else
    }
    if (nMix > 1 && mixSampleMethod == 1 && omegaShareSubpop.n_elem == (unsigned int)nphi1) {
      // MSAEM stratified init: nudge each subject's MCMC draw/prior mean toward its
      // best-fitting hypothesis (mixNaiveClassify) so iteration 0 isn't symmetric.
      uvec cls = mixNaiveClassify(1.5);
      for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
        unsigned int subpop = omegaShareSubpop(c);
        if (subpop < 1) continue;
        double sdCol = std::sqrt(Gamma2_phi1(c, c));
        if (!std::isfinite(sdCol) || sdCol <= 0) continue;
        for (int i = 0; i < N; i++) {
          double nudge = (cls(i) == subpop ? 1.0 : -1.0) * 2.0 * sdCol;
          mprior_phi1(i, c) += nudge;
          for (int k = 0; k < nmc; k++) {
            phiM(i + k * N, i1(c)) += nudge;
          }
        }
      }
    }
    if (nSaCov > 0) { HaSa = zeros<mat>(nb_param, nb_param); covCount = 0; }
    for (unsigned int kiter=0; kiter<(unsigned int)(niter + nSaCov); kiter++) {
      // entering the SA covariance phase: snapshot the converged estimate so it can be
      // restored afterward (the cov-phase iterations fluctuate the parameters).
      if (nSaCov > 0 && kiter == (unsigned int)niter) {
        _savPlambda = Plambda; _savGamma2_phi1 = Gamma2_phi1; _savGamma2_phi0 = Gamma2_phi0;
        _savGamma2_phi1Report = Gamma2_phi1Report; _savMprior_phi1 = mprior_phi1;
        _savMprior_phi0 = mprior_phi0; _savAres = ares; _savBres = bres; _savCres = cres;
        _savLres = lres; _savVcsig2 = vcsig2; _savPhiM = phiM; _savHa = Ha;
        if (nMix > 1) { _savMixProb = mixProb; _savMixWeights = mixWeights; }
      }
      // End of burn-in: re-decide which etas the data informs, now that theta and Omega
      // have moved off the initial estimates the first test was run at.
      if (ueRevisitIter >= 0 && kiter == (unsigned int)ueRevisitIter) {
        revisitUninformativeEtas();
      }
      IGamma2_phi1=invSympdNearPd(Gamma2_phi1, "Gamma2_phi1 (Omega)");
      gamma2_phi1=Gamma2_phi1.diag();
      D1Gamma21=LCOV1*IGamma2_phi1;
      D2Gamma21=D1Gamma21*LCOV1.t();
      CGamma21=COV21%D2Gamma21;

      IGamma2_phi0=invSympdNearPd(Gamma2_phi0, "Gamma2_phi0 (Omega)");
      gamma2_phi0=Gamma2_phi0.diag();
      D1Gamma20=LCOV0*IGamma2_phi0;
      D2Gamma20=D1Gamma20*LCOV0.t();
      CGamma20=COV20%D2Gamma20;

      //    MCMC
      mcmcphi mphi1, mphi0;
      set_mcmcphi(mphi1, i1, nphi1, Gamma2_phi1, IGamma2_phi1, mprior_phi1);
      set_mcmcphi(mphi0, i0, nphi0, Gamma2_phi0, IGamma2_phi0, mprior_phi0);

      // CHG hard coded 20
      int nu1, nu2, nu3;
      if (kiter==0) {
        nu1=20*nu(0);
        nu2=20*nu(1);
        nu3=20*nu(2);
      } else {
        nu1=nu(0);
        nu2=nu(1);
        nu3=nu(2);
      }

      mat Statphi11 = zeros<mat>(N, nphi1);
      mat Statphi01 = zeros<mat>(N, nphi0);
      mat Statphi12 = zeros<mat>(nphi1, nphi1);
      mat Statphi02 = zeros<mat>(nphi0, nphi0);
      field<mat> Statphi11_mix(std::max(nMix, 1));
      for (int _j = 0; _j < std::max(nMix, 1); _j++) {
        Statphi11_mix(_j) = zeros<mat>(N, nphi1);
      }
      double statr[MAXENDPNT];
      for (int b = 0; b < nendpnt; b++) {
        statr[b] = 0.0;
      }
      arResetMstep();
      double resk = 0.0;
      vec D1 = zeros<vec>(nb_param);
      mat D11 = zeros<mat>(nb_param, nb_param);
      mat D2 = zeros<mat>(nb_param, nb_param);
      mat d2logk = zeros<mat>(nb_param, nb_param);
      mat resy(nendpnt, nmc);  // resy(b, k): endpoint b's residual SSR for chain k
      vec fsM;
      fsM.set_size(0);

      if (nMix > 1 && mixSampleMethod == 1) {
        // MSAEM (Lavielle & Mbogning 2014): S-step uses a single MCMC trajectory per
        // subject (phiM) targeting the marginal mixture density; no label is simulated/masked.
        vec U_y = mixObsLoss(phiM, mx);
        if (nphi1 > 0) {
          vec U_phi;
          do_mcmc_msaem(1, nu1, mx, mphi1, phiM, U_y, U_phi, (int)kiter);
          mat dphi = phiM.cols(i1) - mphi1.mprior_phiM;
          U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi1), 1);
          do_mcmc_msaem(2, nu2, mx, mphi1, phiM, U_y, U_phi, (int)kiter);
          do_mcmc_msaem(3, nu3, mx, mphi1, phiM, U_y, U_phi, (int)kiter);
        }
        if (nphi0 > 0) {
          vec U_phi;
          do_mcmc_msaem(1, nu1, mx, mphi0, phiM, U_y, U_phi, (int)kiter);
          mat dphi = phiM.cols(i0) - mphi0.mprior_phiM;
          U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi0), 1);
          do_mcmc_msaem(2, nu2, mx, mphi0, phiM, U_y, U_phi, (int)kiter);
          do_mcmc_msaem(3, nu3, mx, mphi0, phiM, U_y, U_phi, (int)kiter);
        }
        if (DEBUG > 0) Rcout << "mcmc successful (msaem)\n";
        if (kiter < (unsigned int)niter) phiFile << phiM;

        // E-step: posterior responsibility gamma_{i,m} (softmax, see mixWeights below) from
        // the one simulated phi, plus per-hypothesis predictions for the residual term below.
        field<vec> fsave_hyp(nMix);
        field<vec> cens_hyp(nMix);
        field<vec> limit_hyp(nMix);
        mat Ly(N, nMix, fill::zeros);
        for (int mHyp = 0; mHyp < nMix; mHyp++) {
          current_saem_state->_saemMixest = mHyp + 1;
          mat hypMat = user_fn(phiM, evt, optM);
          fsave_hyp(mHyp) = hypMat.col(0);
          vec fHyp = hypMat.col(0);
          vec censHyp = hypMat.col(1);
          vec limitHyp = hypMat.col(2);
          cens_hyp(mHyp) = censHyp;
          limit_hyp(mHyp) = limitHyp;
          mat DYFhyp = zeros<mat>(mlen, nM);
          if (distribution == 1) {
            vec yt = hasFixedObsTransform ? yTrans : y;
            if (!hasFixedObsTransform) {
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                yt(i) = _powerD(y(i), lambda(cur), yj(cur), low(cur), hi(cur));
              }
            }
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            for (int k = 0; k < nmc; k++) {
              int obs_start = k * ntotal;
              vec fk = fHyp.subvec(obs_start, obs_start + ntotal - 1);
              const vec censk = censHyp.subvec(obs_start, obs_start + ntotal - 1);
              const vec limitk = limitHyp.subvec(obs_start, obs_start + ntotal - 1);
              _scratch_ft = fk;
              _scratch_limitT = limitk;
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ft(i) = _powerD(fk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fk(i), false, true);
              }
              saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = indio + (arma::uword)k * stride;
              DYFhyp(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              applyCensLoss(DYFhyp, _scratch_indio, censk, yt, _scratch_limitT,
                            _scratch_ftAr, _scratch_gAr);
            }
          } else if (distribution == 2) {
            for (int k = 0; k < nmc; k++) {
              vec fk = fHyp.subvec(k * ntotal, (k + 1) * ntotal - 1);
              uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
              DYFhyp(indio_k) = -y % log(fk) + fk;
            }
          } else if (distribution == 3) {
            for (int k = 0; k < nmc; k++) {
              vec fk = fHyp.subvec(k * ntotal, (k + 1) * ntotal - 1);
              uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
              DYFhyp(indio_k) = -y % log(fk) - (1 - y) % log(1 - fk);
            }
          }
          vec U_y_hyp = sum(DYFhyp, 0).t();
          for (int i = 0; i < N; i++) {
            double sumL = 0.0;
            for (int k = 0; k < nmc; k++) sumL += U_y_hyp(i + k * N);
            Ly(i, mHyp) = sumL / nmc;
          }
        }
        current_saem_state->_saemMixest = 0;

        for (int i = 0; i < N; i++) {
          double minL = Ly(i, 0);
          for (int j = 1; j < nMix; j++) {
            if (Ly(i, j) < minL) minL = Ly(i, j);
          }
          rowvec w_i(nMix);
          double sumW = 0.0;
          for (int j = 0; j < nMix; j++) {
            w_i(j) = mixProb(j) * exp(minL - Ly(i, j));
            sumW += w_i(j);
          }
          if (!mixProbRegress) {
            if (sumW > 0.0) {
              mixWeights.row(i) = w_i / sumW;
            } else {
              mixWeights.row(i) = mixProb.t();
            }
          }
        }

        // Leaspy-style prior-only responsibility (cf. mixWeights' joint obs+prior NLL): weights
        // the variance update by prior-deviation alone, since an unmatched column stays near its prior.
        if (omegaShareSubpop.n_elem == (unsigned int)nphi1) {
          mat priorPenalty(N, nMix, fill::zeros);
          bool anyOwned = false;
          for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
            unsigned int subpop = omegaShareSubpop(c);
            if (subpop < 1 || subpop > (unsigned int)nMix) continue;
            anyOwned = true;
            double v = Gamma2_phi1(c, c);
            if (!std::isfinite(v) || v <= 0) continue;
            for (int i = 0; i < N; i++) {
              double dsum = 0.0;
              for (int k = 0; k < nmc; k++) {
                dsum += phiM(i + k * N, i1(c)) - mprior_phi1(i, c);
              }
              double d = dsum / nmc;
              priorPenalty(i, subpop - 1) += 0.5 * d * d / v;
            }
          }
          if (anyOwned) {
            // Split-ETA columns' genuinely-owned column shows a *larger* prior deviation
            // (unlike leaspy's smaller-is-better case), so softmax uses max, not min.
            for (int i = 0; i < N; i++) {
              double maxP = priorPenalty(i, 0);
              for (int j = 1; j < nMix; j++) {
                if (priorPenalty(i, j) > maxP) maxP = priorPenalty(i, j);
              }
              rowvec pw(nMix);
              double sumP = 0.0;
              for (int j = 0; j < nMix; j++) {
                pw(j) = exp(priorPenalty(i, j) - maxP);
                sumP += pw(j);
              }
              if (sumP > 0.0) {
                priorWeights.row(i) = pw / sumP;
              } else {
                priorWeights.row(i).fill(1.0 / nMix);
              }
            }
          }
        }

        // AE-step: unweighted accumulation (phi doesn't vary by hypothesis, so
        // sum_m(gamma_{i,m}*phi_i) collapses to phi_i); only the residual term needs gamma weighting.
        for (int k = 0; k < nmc; k++) {
          phi.slice(k) = phiM.rows(span(k * N, (k + 1) * N - 1));

          Statphi11 += phi.slice(k).cols(i1);
          Statphi01 += phi.slice(k).cols(i0);
          mat phik = phi.slice(k);
          mat phi1k = phik.cols(i1);
          mat phi0k = phik.cols(i0);
          Statphi12 = Statphi12 + phi1k.t() * phi1k;
          Statphi02 = Statphi02 + phi0k.t() * phi0k;

          for (int b = 0; b < nendpnt; b++) {
            double resk = 0.0;
            for (int mHyp = 0; mHyp < nMix; mHyp++) {
              vec fk = fsave_hyp(mHyp).subvec(k * ntotal, (k + 1) * ntotal - 1);
              fk = fk(ix_sorting);
              vec f_cur = fk(span(y_offset(b), y_offset(b+1)-1));
              vec y_cur;
              if (hasFixedObsTransform) {
                y_cur = ysTrans(span(y_offset(b), y_offset(b+1)-1));
              } else {
                y_cur = ys(span(y_offset(b), y_offset(b+1)-1));
              }
              vec censK = cens_hyp(mHyp).subvec(k * ntotal, (k + 1) * ntotal - 1);
              censK = censK(ix_sorting);
              vec limitK = limit_hyp(mHyp).subvec(k * ntotal, (k + 1) * ntotal - 1);
              limitK = limitK(ix_sorting);
              // #916: data augmentation -- replace censored (M3/M4) rows with a
              // simulated draw before building the residual SSR below.
              y_cur = augmentCensY(b, f_cur, y_cur,
                                    censK(span(y_offset(b), y_offset(b+1)-1)),
                                    limitK(span(y_offset(b), y_offset(b+1)-1)),
                                    (int)kiter, k, mHyp);
              resk += arResk(b, f_cur, y_cur, mHyp);
            }
            statr[b] += resk;
            resy(b, k) = resk;
          }

          mat dphi1k = phi1k - mprior_phi1;
          mat dphi0k = phi0k - mprior_phi0;
          vec sdg1 = sum(dphi1k % dphi1k, 0).t() / gamma2_phi1;
          mat Md1 = (IGamma2_phi1 * (dphi1k.t() * Mcovariables)).t();
          mat Md0 = (IGamma2_phi0 * (dphi0k.t() * Mcovariables)).t();
          vec d1_mu_phi1 = Md1(ind_cov1);
          vec d1_mu_phi0 = Md0(ind_cov0);
          vec d1_loggamma2_phi1 = 0.5 * sdg1 - 0.5 * N;
          vec d1_logsigma2(nResidEp);
          fillResidLogSigma2(k, resy, d1_logsigma2, d2logk);
          vec d1logk = join_cols(d1_mu_phi1, join_cols(d1_mu_phi0, join_cols(d1_loggamma2_phi1, d1_logsigma2)));
          D1 = D1 + d1logk;
          D11 = D11 + d1logk * d1logk.t();

          vec w2phi = -0.5 * sdg1;
          for (int j = 0, l = 0; j < nphi1; j++) {
            for (unsigned int jj = 0; jj < pc1(j); jj++) {
              double temp = -dot(COV1.col(l), dphi1k.col(j)) / gamma2_phi1(j);
              d2logk(l, nlambda + j) = temp;
              d2logk(nlambda + j, l) = temp;
              l = l + 1;
            }
            d2logk(nlambda + j, nlambda + j) = w2phi(j);
          }
          D2 = D2 + d2logk;
        }
      } else if (nMix > 1) {
        field<mat> DYF_mix(nMix);
        field<vec> U_y_mix(nMix);

        // mixProbMethod="regress": membership is fixed, so run the S-step ONCE
        // (jMix=0) with each subject solved under its own component via the
        // per-subject mixest regressor (mx.optM["mixest"], _saemMixest kept 0),
        // then fan the result out to every component below.  The soft-EM path
        // runs one full per-component chain per component.
        int nMixLoop = mixProbRegress ? 1 : nMix;
        for (int jMix = 0; jMix < nMixLoop; jMix++) {
          if (!mixProbRegress) current_saem_state->_saemMixest = jMix + 1;

          mat &cur_phiM = phiM_mix(jMix);
          vec &cur_fsave = fsave_mix(jMix);
          vec &cur_cens = cens_mix(jMix);
          vec &cur_limit = limit_mix(jMix);

          mat &cur_DYF = DYF_mix(jMix);
          cur_DYF.zeros(mlen, nM);

          vec f = cur_fsave;
          cur_fsave = f;

          if (distribution == 1) {
            vec yt = hasFixedObsTransform ? yTrans : y;
            if (!hasFixedObsTransform) {
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                yt(i) = _powerD(y(i), lambda(cur), yj(cur), low(cur), hi(cur));
              }
            }
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            for (int k = 0; k < nmc; k++) {
              int obs_start = k * ntotal;
              vec fk = f.subvec(obs_start, obs_start + ntotal - 1);
              const vec censk = cur_cens.subvec(obs_start, obs_start + ntotal - 1);
              const vec limitk = cur_limit.subvec(obs_start, obs_start + ntotal - 1);
              _scratch_ft = fk;
              _scratch_limitT = limitk;
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ft(i) = _powerD(fk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fk(i), false, true);
              }
              saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = indio + (arma::uword)k * stride;
              cur_DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              applyCensLoss(cur_DYF, _scratch_indio, censk, yt, _scratch_limitT,
                            _scratch_ftAr, _scratch_gAr);
            }
          } else if (distribution == 2) {
            for (int k = 0; k < nmc; k++) {
              vec fk = f.subvec(k * ntotal, (k + 1) * ntotal - 1);
              uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
              cur_DYF(indio_k) = -y % log(fk) + fk;
            }
          } else if (distribution == 3) {
            for (int k = 0; k < nmc; k++) {
              vec fk = f.subvec(k * ntotal, (k + 1) * ntotal - 1);
              uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
              cur_DYF(indio_k) = -y % log(fk) - (1 - y) % log(1 - fk);
            }
          }

          vec U_y = sum(cur_DYF, 0).t();

          if (nphi1 > 0) {
            vec U_phi;
            do_mcmc(1, nu1, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
            mat dphi = cur_phiM.cols(i1) - mphi1.mprior_phiM;
            U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi1), 1);
            do_mcmc(2, nu2, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
            do_mcmc(3, nu3, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
          }
          if (nphi0 > 0) {
            vec U_phi;
            do_mcmc(1, nu1, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
            mat dphi = cur_phiM.cols(i0) - mphi0.mprior_phiM;
            U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi0), 1);
            do_mcmc(2, nu2, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
            do_mcmc(3, nu3, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1);
          }

          // Joint NLL (U_y + U_phi) for mixture weights: U_y alone is insufficient since MCMC
          // adapts eta to fit every component equally; the prior penalty U_phi (Mahalanobis
          // distance from the component mean) is what discriminates the true component.
          {
            vec joint_nll(N * nmc);
            for (int k = 0; k < nmc; k++) {
              // Observation loss for this MCMC sample (rows i + k*N for each subject)
              for (int i = 0; i < N; i++) {
                joint_nll(i + k * N) = 0.0;
                for (int m = 0; m < mlen; m++) {
                  joint_nll(i + k * N) += cur_DYF(m, i + k * N);
                }
              }
              // mprior_phiM is already restricted to i1/i0 columns (see set_mcmcphi()), so it
              // must not be column-indexed again; slice the k-th row-block via .rows() to avoid
              // a full deep copy of cur_phiM.
              if (nphi1 > 0) {
                mat block1     = cur_phiM.rows(k * N, (k + 1) * N - 1);
                mat phi1_k     = block1.cols(i1);
                mat prior1_k   = mphi1.mprior_phiM.rows(k * N, (k + 1) * N - 1);
                mat dphi1_k    = phi1_k - prior1_k;
                vec uphi1_k    = 0.5 * sum(dphi1_k % (dphi1_k * IGamma2_phi1), 1);
                for (int i = 0; i < N; i++) joint_nll(i + k * N) += uphi1_k(i);
              }
              if (nphi0 > 0) {
                mat block0     = cur_phiM.rows(k * N, (k + 1) * N - 1);
                mat phi0_k     = block0.cols(i0);
                mat prior0_k   = mphi0.mprior_phiM.rows(k * N, (k + 1) * N - 1);
                mat dphi0_k    = phi0_k - prior0_k;
                vec uphi0_k    = 0.5 * sum(dphi0_k % (dphi0_k * IGamma2_phi0), 1);
                for (int i = 0; i < N; i++) joint_nll(i + k * N) += uphi0_k(i);
              }
            }
            U_y_mix(jMix) = joint_nll;
          }
          if (!mixProbRegress) current_saem_state->_saemMixest = jMix + 1;
          cur_fsave = user_fn(cur_phiM, mixProbRegress ? mx.evtM : evt,
                              mixProbRegress ? mx.optM : optM).col(0);
        }
        current_saem_state->_saemMixest = 0;

        // Fixed membership: fan the single solved chain out to every component
        // so the downstream AE-step / sufficient-stat accumulation (which sums
        // mixWeights(i,j)*phiM_mix(j)) sees each subject's one solve under a 0/1
        // weight.  Copy phi/fsave/cens/limit from component 0.
        if (mixProbRegress) {
          for (int j = 1; j < nMix; j++) {
            phiM_mix(j) = phiM_mix(0);
            fsave_mix(j) = fsave_mix(0);
            cens_mix(j) = cens_mix(0);
            limit_mix(j) = limit_mix(0);
          }
        }

        // Compute posterior weights a_ji (soft-EM only; regress keeps the fixed
        // 0/1 mixWeights and only ran one component's chain, so U_y_mix(1..) are
        // empty and must not be read).
        if (!mixProbRegress) {
        mat L_ji(N, nMix);
        for (int j = 0; j < nMix; j++) {
          const vec &U_y_j = U_y_mix(j);
          for (int i = 0; i < N; i++) {
            double sum_L = 0.0;
            for (int k = 0; k < nmc; k++) {
              sum_L += U_y_j(i + k * N);
            }
            L_ji(i, j) = sum_L / nmc;
          }
        }

        for (int i = 0; i < N; i++) {
          double min_L = L_ji(i, 0);
          for (int j = 1; j < nMix; j++) {
            if (L_ji(i, j) < min_L) min_L = L_ji(i, j);
          }
          rowvec w_i(nMix);
          double sum_w = 0.0;
          for (int j = 0; j < nMix; j++) {
            w_i(j) = mixProb(j) * exp(min_L - L_ji(i, j));
            sum_w += w_i(j);
          }
          if (sum_w > 0.0) {
            mixWeights.row(i) = w_i / sum_w;
          } else {
            mixWeights.row(i) = mixProb.t();
          }
        }
        }

        if (DEBUG > 0) Rcout << "mcmc successful (mixture)\n";

        mat phiM_weighted(N * nmc, nphi, fill::zeros);
        for (int jMix = 0; jMix < nMix; jMix++) {
          for (int row = 0; row < N * nmc; row++) {
            int i_subj = row % N;
            phiM_weighted.row(row) += mixWeights(i_subj, jMix) * phiM_mix(jMix).row(row);
          }
        }
        if (kiter < (unsigned int)niter) phiFile << phiM_weighted;
        for (int k = 0; k < nmc; k++) {
          phi.slice(k) = phiM_weighted.rows(span(k * N, (k + 1) * N - 1));
        }

        // Integration / sufficient stats accumulation for nMix > 1
        field<mat> phi1_mix(nMix);
        field<mat> phi0_mix(nMix);
        for (int j = 0; j < nMix; j++) {
          phi1_mix(j) = phiM_mix(j).cols(i1);
          phi0_mix(j) = phiM_mix(j).cols(i0);
        }

        for (int k = 0; k < nmc; k++) {
          mat dphi1k_w(N, nphi1, fill::zeros);
          mat dphi0k_w(N, nphi0, fill::zeros);
          vec sdg1_w(nphi1, fill::zeros);
          vec sdg0_w(nphi0, fill::zeros);

          for (int i = 0; i < N; i++) {
            int row_idx = i + k * N;
            rowvec phi1_w = zeros<rowvec>(nphi1);
            rowvec phi0_w = zeros<rowvec>(nphi0);

            for (int jMix = 0; jMix < nMix; jMix++) {
              rowvec phi1_ji = phi1_mix(jMix).row(row_idx);
              rowvec phi0_ji = phi0_mix(jMix).row(row_idx);

              // Unblended per-component accumulation (no mixWeights factor) so the theta
              // M-step for mixture-owned columns isn't diluted by the other component.
              Statphi11_mix(jMix).row(i) += phi1_ji;

              phi1_w += mixWeights(i, jMix) * phi1_ji;
              phi0_w += mixWeights(i, jMix) * phi0_ji;

              Statphi12 += mixWeights(i, jMix) * phi1_ji.t() * phi1_ji;
              Statphi02 += mixWeights(i, jMix) * phi0_ji.t() * phi0_ji;

              rowvec dphi1_ji = phi1_ji - mprior_phi1.row(i);
              rowvec dphi0_ji = phi0_ji - mprior_phi0.row(i);

              dphi1k_w.row(i) += mixWeights(i, jMix) * dphi1_ji;
              dphi0k_w.row(i) += mixWeights(i, jMix) * dphi0_ji;

              sdg1_w += mixWeights(i, jMix) * (dphi1_ji % dphi1_ji).t();
              sdg0_w += mixWeights(i, jMix) * (dphi0_ji % dphi0_ji).t();
            }

            Statphi11.row(i) += phi1_w;
            Statphi01.row(i) += phi0_w;
          }

          for (int b = 0; b < nendpnt; b++) {
            double resk = 0.0;
            for (int jMix = 0; jMix < nMix; jMix++) {
              vec fk = fsave_mix(jMix).subvec(k * ntotal, (k + 1) * ntotal - 1);
              fk = fk(ix_sorting);
              vec f_cur = fk(span(y_offset(b), y_offset(b+1)-1));
              vec y_cur;
              if (hasFixedObsTransform) {
                y_cur = ysTrans(span(y_offset(b), y_offset(b+1)-1));
              } else {
                y_cur = ys(span(y_offset(b), y_offset(b+1)-1));
              }
              vec censK = cens_mix(jMix).subvec(k * ntotal, (k + 1) * ntotal - 1);
              censK = censK(ix_sorting);
              vec limitK = limit_mix(jMix).subvec(k * ntotal, (k + 1) * ntotal - 1);
              limitK = limitK(ix_sorting);
              // #916: data augmentation -- replace censored (M3/M4) rows with a
              // simulated draw before building the residual SSR below.
              y_cur = augmentCensY(b, f_cur, y_cur,
                                    censK(span(y_offset(b), y_offset(b+1)-1)),
                                    limitK(span(y_offset(b), y_offset(b+1)-1)),
                                    (int)kiter, k, jMix);
              resk += arResk(b, f_cur, y_cur, jMix);
            }
            statr[b] += resk;
            resy(b, k) = resk;
          }

          vec sdg1 = sdg1_w / gamma2_phi1;
          mat Md1 = (IGamma2_phi1 * (dphi1k_w.t() * Mcovariables)).t();
          mat Md0 = (IGamma2_phi0 * (dphi0k_w.t() * Mcovariables)).t();
          vec d1_mu_phi1 = Md1(ind_cov1);
          vec d1_mu_phi0 = Md0(ind_cov0);
          vec d1_loggamma2_phi1 = 0.5 * sdg1 - 0.5 * N;
          vec d1_logsigma2(nResidEp);
          fillResidLogSigma2(k, resy, d1_logsigma2, d2logk);
          vec d1logk = join_cols(d1_mu_phi1, join_cols(d1_mu_phi0, join_cols(d1_loggamma2_phi1, d1_logsigma2)));
          D1 = D1 + d1logk;
          D11 = D11 + d1logk * d1logk.t();

          vec w2phi = -0.5 * sdg1;
          for (int j = 0, l = 0; j < nphi1; j++) {
            for (unsigned int jj = 0; jj < pc1(j); jj++) {
              double temp = -dot(COV1.col(l), dphi1k_w.col(j)) / gamma2_phi1(j);
              d2logk(l, nlambda + j) = temp;
              d2logk(nlambda + j, l) = temp;
              l = l + 1;
            }
            d2logk(nlambda + j, nlambda + j) = w2phi(j);
          }
          D2 = D2 + d2logk;
        }
        for (int k = 0; k < nmc; k++) {
          vec fk_w(ntotal, fill::zeros);
          for (int jMix = 0; jMix < nMix; jMix++) {
            vec fk_j = fsave_mix(jMix).subvec(k * ntotal, (k + 1) * ntotal - 1);
            fk_j = fk_j(ix_sorting);
            for (int i = 0; i < ntotal; i++) {
              int idx_orig = ix_sorting(i);
              int i_subj = obs_subject(idx_orig);
              fk_w[i] += mixWeights(i_subj, jMix) * fk_j[i];
            }
          }
          fsM = join_cols(fsM, fk_w);
        }
      } else {
        vec f = fsave;
        fsave = f;
        if (distribution == 1){
          // Build yt once: does not depend on chain index k
          vec yt = hasFixedObsTransform ? yTrans : y;
          if (!hasFixedObsTransform) {
            for (int i = ntotal; i--;) {
              int cur = ix_endpnt(i);
              yt(i) = _powerD(y(i), lambda(cur), yj(cur), low(cur), hi(cur));
            }
          }
          const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
          for (int k = 0; k < nmc; k++) {
            int obs_start = k * ntotal;
            vec fk = f.subvec(obs_start, obs_start + ntotal - 1);
            const vec censk = cens.subvec(obs_start, obs_start + ntotal - 1);
            const vec limitk = limit.subvec(obs_start, obs_start + ntotal - 1);
            _scratch_ft = fk;
            _scratch_limitT = limitk;
            for (int i = ntotal; i--;) {
              int cur = ix_endpnt(i);
              _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
              _scratch_ft(i) = _powerD(fk(i), lambda(cur), yj(cur), low(cur), hi(cur));
              _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fk(i), false, true);
            }
            saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
            _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
            _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
            _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
            _scratch_indio = indio + (arma::uword)k * stride;
            DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
            applyCensLoss(DYF, _scratch_indio, censk, yt, _scratch_limitT,
                          _scratch_ftAr, _scratch_gAr);
          }
        } else if (distribution == 2){
          for (int k = 0; k < nmc; k++) {
            vec fk = f.subvec(k * ntotal, (k + 1) * ntotal - 1);
            uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
            DYF(indio_k) = -y % log(fk) + fk;
          }
        } else if (distribution == 3) {
          for (int k = 0; k < nmc; k++) {
            vec fk = f.subvec(k * ntotal, (k + 1) * ntotal - 1);
            uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
            DYF(indio_k) = -y % log(fk) - (1 - y) % log(1 - fk);
          }
        } else if (distribution == 4) {
          // General log-likelihood endpoint (ll() ~ expr): the model returns the
          // per-observation log-likelihood as its prediction (rx_pred_ ~ <ll>), so
          // the observation loss is simply -ll and the standard RWM kernels run
          // unchanged.  Reachable for est="saem" when .saemGeneralLik(ui) is true
          // (the transform-normal assertion is skipped for such a model).
          // This computes U_y, the CURRENT state's reference likelihood every
          // do_mcmc call below is scored against -- it needs the identical
          // sentinel-inversion/ceiling clamp do_mcmc's own case 4 applies (see
          // its comment), or a bad/NaN solve at the CURRENT phi (e.g. right at
          // kiter=0, before any exploration has adapted the state) makes U_y
          // itself an artificially huge NEGATIVE ("great") reference, which
          // either traps the chain there or makes every subsequent acceptance
          // test's deltu meaningless (both sides dominated by the same
          // ~-1e99 sentinel).
          for (int k = 0; k < nmc; k++) {
            vec fk = f.subvec(k * ntotal, (k + 1) * ntotal - 1);
            fk.elem(find(fk >= 1.0e99)).fill(_saemGenLikBadSolvePenalty);
            fk.elem(find(fk > _saemGenLikCeiling)).fill(_saemGenLikCeiling);
            uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
            DYF(indio_k) = -fk;
          }
        }
        else {
          RSprintf("unknown distribution (id=%d)\n", distribution);
          return;
        }
        //U_y is a vec of subject llik; summed over obs for each subject
        vec U_y=sum(DYF, 0).t();

        if(nphi1>0) {
          vec U_phi;
          do_mcmc(1, nu1, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
          mat dphi = phiM.cols(i1)-mphi1.mprior_phiM;
          U_phi    = 0.5*sum(dphi%(dphi*IGamma2_phi1),1);
          do_mcmc(2, nu2, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
          do_mcmc(3, nu3, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
        }
        if(nphi0>0) {
          vec U_phi;
          do_mcmc(1, nu1, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
          mat dphi = phiM.cols(i0)-mphi0.mprior_phiM;
          U_phi    = 0.5*sum(dphi%(dphi*IGamma2_phi0),1);
          do_mcmc(2, nu2, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
          do_mcmc(3, nu3, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter);
        }
        if (DEBUG>0) Rcout << "mcmc successful\n";
        if (kiter < (unsigned int)niter) phiFile << phiM;

        //integration
        for(int k=0; k<nmc; k++) {
          phi.slice(k)=phiM.rows(span(k*N, (k+1)*N-1));

          Statphi11 += phi.slice(k).cols(i1);
          Statphi01 += phi.slice(k).cols(i0);
          mat phik=phi.slice(k);
          mat phi1k=phik.cols(i1);
          mat phi0k=phik.cols(i0);
          Statphi12=Statphi12+phi1k.t()*phi1k;
          Statphi02=Statphi02+phi0k.t()*phi0k;

          vec fk = fsave(span(k*ntotal, (k+1)*ntotal-1));
          fk = fk(ix_sorting);    //sorted by endpnt
          fsM = join_cols(fsM, fk);
          // vec resid_all(ys.size());// = ys - fk;
          vec censK = cens.subvec(k*ntotal, (k+1)*ntotal-1);
          censK = censK(ix_sorting);
          vec limitK = limit.subvec(k*ntotal, (k+1)*ntotal-1);
          limitK = limitK(ix_sorting);
          vec gk, y_cur, f_cur;
          double ft, fa;
          //loop thru endpoints here
          // general log-likelihood (distribution==4) has no residual error, so
          // skip the residual SSR accumulation entirely (fsM above is kept for
          // downstream predictions)
          if (distribution != 4)
          for(int b=0; b<nendpnt; ++b) {
            if (hasFixedObsTransform) {
              y_cur = ysTrans(span(y_offset(b), y_offset(b+1)-1));
            } else {
              y_cur = ys(span(y_offset(b), y_offset(b+1)-1));
            }
            f_cur = fk(span(y_offset(b), y_offset(b+1)-1));
            // #916: data augmentation -- replace censored (M3/M4) rows with a
            // simulated draw before building the residual SSR below.
            y_cur = augmentCensY(b, f_cur, y_cur,
                                  censK(span(y_offset(b), y_offset(b+1)-1)),
                                  limitK(span(y_offset(b), y_offset(b+1)-1)),
                                  (int)kiter, k, -1);
            vec resid(y_cur.size());
            for (int i = y_cur.size(); i--;){
              resid(i) = y_cur[i];
              if (std::isnan(resid(i))) {
                Rcpp::stop(_("NaN in data or transformed data; please check transformation/data"));
              }
              ft = _powerD(f_cur[i], lambda(b), yj(b), low(b), hi(b));
              resid(i) -=  ft;
              if (res_mod(b) == rmProp) {
                fa = handleF(propT(b), ft, f_cur[i], true, true);
                if (fa <= double_xmin) {
                  fa = 1;
                }
                resid(i) = resid(i)/fa;
              }
            }

            if (arActive(b)) {
              // AR(1): whitened SSR + accumulate residual pairs for arUpdateCor
              resk = arResk(b, f_cur, y_cur, -1);
            } else if (res_mod(b) <= rmProp) {
              resk = dot(resid, resid);
              if (resk > xmax) {
                resk = xmax;
              } else if (resk < double_xmin) {
                resk = double_xmin;
              }
            }
            else {
              resk = 1;                                              //FIXME
            }

            statr[b]=statr[b]+resk;
            resy(b, k) = resk;
          }
          if (DEBUG>1) Rcout << "star[] successful\n";

          mat dphi1k=phi1k-mprior_phi1;
          mat dphi0k=phi0k-mprior_phi0;
          vec sdg1=sum(dphi1k%dphi1k,0).t()/gamma2_phi1;
          mat Md1=(IGamma2_phi1*(dphi1k.t()*Mcovariables)).t();
          mat Md0=(IGamma2_phi0*(dphi0k.t()*Mcovariables)).t();
          vec d1_mu_phi1=Md1(ind_cov1);                              //CHK!! vec or mat
          vec d1_mu_phi0=Md0(ind_cov0);                              //CHK!! vec or mat
          vec d1_loggamma2_phi1=0.5*sdg1-0.5*N;
          // general log-likelihood (distribution==4): no residual param, so
          // nResidEp==0 and this block is empty
          vec d1_logsigma2(nResidEp);
          if (distribution != 4) fillResidLogSigma2(k, resy, d1_logsigma2, d2logk);
          vec d1logk=join_cols(d1_mu_phi1, join_cols(d1_mu_phi0, join_cols(d1_loggamma2_phi1, d1_logsigma2)));
          D1 = D1+d1logk;
          D11= D11+d1logk*d1logk.t();

          vec w2phi=-0.5*sdg1;                                       //CHK!!!
          for(int j=0, l=0; j<nphi1; j++) {
            for(unsigned int jj=0; jj<pc1(j); jj++) {
              double temp=-dot(COV1.col(l),dphi1k.col(j))/gamma2_phi1(j);
              d2logk(l,nlambda+j)=temp;
              d2logk(nlambda+j,l)=temp;
              l=l+1;
            }
            d2logk(nlambda+j,nlambda+j)=w2phi(j);
          }
          D2=D2+d2logk;
        }
      }//k
      if (DEBUG>0) Rcout << "integration successful\n";

      if (nMix > 1 && mixSampleMethod == 1 && omegaShareSubpop.n_elem == (unsigned int)nphi1) {
        // Split-ETA columns: separation can still be developing after pas(kiter) decays
        // (same freeze-before-signal issue as mixProb), so give these columns their own
        // extended-flat schedule; non-split columns and "parallel" fits are unaffected.
        double pasMsaemSplit = (kiter < (unsigned int)(0.8 * niter)) ? 1.0 : pas(kiter);
        for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
          double pasUse = (omegaShareSubpop(c) >= 1) ? pasMsaemSplit : pas(kiter);
          statphi11.col(c) = statphi11.col(c) + pasUse * (Statphi11.col(c) / nmc - statphi11.col(c));
        }
      } else {
        statphi11=statphi11+pas(kiter)*(Statphi11/nmc-statphi11);
      }
      if (nMix > 1 && mixSampleMethod == 0) {
        // Only "parallel" populates the unblended Statphi11_mix accumulator; "msaem" has a
        // single trajectory, so statphi11 above is already clean and statphi11_mix is unused.
        for (int _j = 0; _j < nMix; _j++) {
          statphi11_mix(_j) = statphi11_mix(_j) + pas(kiter)*(Statphi11_mix(_j)/nmc - statphi11_mix(_j));
        }
      }
      statphi12=statphi12+pas(kiter)*(Statphi12/nmc-statphi12);
      statphi01=statphi01+pas(kiter)*(Statphi01/nmc-statphi01);
      // s_{2, k} = statphi02
      statphi02=statphi02+pas(kiter)*(Statphi02/nmc-statphi02);
      for(int b=0; b<nendpnt; ++b) {
        statrese[b]=statrese[b]+pas(kiter)*(statr[b]/nmc-statrese[b]);
      }

      // update parameters
      vec Plambda1, Plambda0;
      Plambda1=inv_sympd(CGamma21)*sum((D1Gamma21%(COV1.t()*statphi11)),1);
      // Split-ETA mixture-owned columns (e.g. eta.cl1/eta.cl2): blended statphi11 dilutes each
      // column with the other component's uninformative draws, coupling tcl1/tcl2 and shrinking
      // apparent BSV. Override these columns' theta via responsibility-weighted regression
      // against the clean statphi11_mix -- safe only when the column has no assumed correlation
      // with any other phi1 column (else Gamma2_phi1's off-diagonal couples it back in).
      if (nMix > 1 && omegaShareSubpop.n_elem == (unsigned int)nphi1) {
        for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
          unsigned int subpop = omegaShareSubpop(c);
          if (subpop < 1 || subpop > (unsigned int)nMix) continue;
          bool diagOnly = true;
          for (unsigned int cc = 0; cc < (unsigned int)nphi1; cc++) {
            if (cc != c && covstruct1(c, cc) != 0) { diagOnly = false; break; }
          }
          if (!diagOnly) continue;
          uvec lambdaIdx = find(LCOV1.col(c) == 1);
          if (lambdaIdx.n_elem == 0) continue;
          vec w = mixWeights.col(subpop - 1);
          mat Xl = COV1.cols(lambdaIdx);
          mat XtW = Xl.t() * diagmat(w);
          mat XtWX = XtW * Xl;
          // "msaem": statphi11 is already the clean single-trajectory statistic. "parallel":
          // fall back to the unblended per-component accumulator (statphi11 there is diluted).
          vec y = (mixSampleMethod == 1) ? statphi11.col(c) : statphi11_mix(subpop - 1).col(c);
          vec XtWy = XtW * y;
          bool fixedCol = false;
          if (fixedIx1.n_elem > 0) {
            for (unsigned int li = 0; li < lambdaIdx.n_elem; li++) {
              if (any(fixedIx1 == lambdaIdx(li))) { fixedCol = true; break; }
            }
          }
          if (fixedCol) continue;
          mat XtWXi;
          if (inv_sympd(XtWXi, XtWX)) {
            Plambda1(lambdaIdx) = XtWXi * XtWy;
          }
        }
      }
      if (fixedIx1.n_elem>0) {
        Plambda1(fixedIx1) = MCOV1(jcov1(fixedIx1));
      }
      MCOV1(jcov1)=Plambda1;
      if (nphi0>0) {
        Plambda0=inv_sympd(CGamma20)*sum((D1Gamma20%(COV0.t()*statphi01)),1);
        if (fixedIx0.n_elem>0) {
          Plambda0(fixedIx0) = MCOV0(jcov0(fixedIx0));
        }
        MCOV0(jcov0)=Plambda0;
      }
      // Phase 4 (SAEM general-likelihood theta plan): once the direct phi1
      // Laplace optimizer owns phi1 (kiter>=niter_phi0 -- the same
      // past-burn-in threshold phi0's own regress mode uses), do NOT
      // overwrite mprior_phi1 with the stochastic sampled-mean update --
      // mirrors skipStochPhi0 exactly.  Only for a general-lik fit whose
      // phi1 map resolved AND innerHess2 built (_saemPhi1PoolReady); every other
      // fit's mprior_phi1 update is completely unchanged.  Row-level
      // phiM/eta values do_mcmc already holds are untouched either way.
      bool skipStochPhi1 = _saemPhi1PoolReady && (kiter >= (unsigned int)niter_phi0);
      if (!skipStochPhi1) {
        mprior_phi1=COV1*MCOV1;
      }
      // nonMuTheta="regress": once the direct phi0 optimizer owns phi0
      // (kiter>=niter_phi0), do NOT overwrite mprior_phi0 with the stochastic
      // sampled-mean update -- that fights the optimizer and lets phi0 drift.
      // The optimizer result from the previous iteration persists as the seed.
      bool skipStochPhi0 = nonMuThetaRegress && (distribution != 4) &&
        (kiter >= (unsigned int)niter_phi0);
      if (!skipStochPhi0) {
        mprior_phi0=COV0*MCOV0;
      }
      if (_saemPhi1PoolReady && kiter >= (unsigned int)niter_phi0 &&
          (kiter - (unsigned int)niter_phi0) % (unsigned int)phi1ThetaEvery == 0) {
        refinePhi1Lik(kiter, pas);
      }
      // The sampled-mean update above only weakly informs fixed-effect-only
      // (phi0) parameters, so once the SA/variance-shrinkage phase has begun,
      // refine them by a direct bounded optimization with the ODE states frozen
      // (saemix ind.fix10).  Enabled for general log-likelihood models
      // (distribution==4) and, via nonMuTheta="regress", for normal models --
      // keeping non-mu thetas as directly-optimized, bound-respecting regressors
      // instead of stochastic phi0 draws.
      // nonMuThetaEvery>1 refines only every k-th iteration; mprior_phi0 holds its
      // last refined value in between (the SA step pas(kiter) moves it a fraction
      // of one optimizer step anyway, so refining every iteration buys little).
      if ((distribution == 4 || nonMuThetaRegress) &&
          nphi0 > 0 && kiter >= (unsigned int)niter_phi0 &&
          (kiter - (unsigned int)niter_phi0) % (unsigned int)nonMuThetaEvery == 0) {
        refinePhi0Lik(kiter, pas);
      }
      mprior_phi0.set_size(N, nphi0);                              // deal w/ nphi0=0
      if (nphi0 > 0) {
        phiM.cols(i0) = repmat(mprior_phi0, nmc, 1);
        if (nMix > 1) {
          for (int jMix = 0; jMix < nMix; jMix++) {
            phiM_mix(jMix).cols(i0) = repmat(mprior_phi0, nmc, 1);
          }
        }
      }

      mat G1=(statphi12+mprior_phi1.t()*mprior_phi1- statphi11.t()*mprior_phi1 - mprior_phi1.t()*statphi11)/N;

      // "msaem" only: G1 above divides by all N uniformly, diluting a split-ETA column's
      // BSV with subjects who rarely belong to that column's component. Override with a
      // responsibility-weighted sample variance around the prior mean (same diagOnly
      // separability argument as the theta override above).
      if (nMix > 1 && mixSampleMethod == 1 && omegaShareSubpop.n_elem == (unsigned int)nphi1) {
        for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
          unsigned int subpop = omegaShareSubpop(c);
          if (subpop < 1 || subpop > (unsigned int)nMix) continue;
          bool diagOnly = true;
          for (unsigned int cc = 0; cc < (unsigned int)nphi1; cc++) {
            if (cc != c && covstruct1(c, cc) != 0) { diagOnly = false; break; }
          }
          if (!diagOnly) continue;
          if (Gamma2_phi1fixed == 1 && any(Gamma2_phi1fixedIx == c * nphi1 + c)) continue;
          vec w = priorWeights.col(subpop - 1);
          double sumW = arma::sum(w);
          if (sumW <= 0.0) continue;
          vec dev = statphi11.col(c) - mprior_phi1.col(c);
          G1(c, c) = arma::sum(w % (dev % dev)) / sumW;
        }
      }

      // Two-level (IOV): the per-occasion columns of one occasion parameter
      // estimate a single Psi, so pool their moments.  Under the equality
      // constraint the maximizer of the complete-data likelihood is the plain
      // mean of the per-occasion moments -- each column carries the same N
      // deviations -- so this stays a closed-form M-step, not a projection.
      poolOmegaGroups(G1);

      if (kiter<=(unsigned int)(nb_sa)) {
        Gamma2_phi1=max(Gamma2_phi1*coef_sa, diagmat(G1));
      } else {
        Gamma2_phi1=G1;
      }
      Gamma2_phi1=Gamma2_phi1%covstruct1;
      // the SA floor above is per-element, so it can pull a pooled group apart
      // again; restore the constraint after it
      poolOmegaGroups(Gamma2_phi1);
      // Split-ETA components sharing an omegaShare group are pooled into a single BSV term
      // (law of total variance) for *reporting only*, into Gamma2_phi1Report; the live
      // Gamma2_phi1 feeding IGamma2_phi1/D1Gamma21 stays untouched so tcl1/tcl2 stay uncoupled.
      Gamma2_phi1Report = Gamma2_phi1;
      if (nMix > 1 && omegaShare.n_elem == (unsigned int)nphi1) {
        unsigned int max_group = 0;
        for (unsigned int i = 0; i < omegaShare.n_elem; ++i) {
          if (omegaShare(i) > max_group) max_group = omegaShare(i);
        }
        for (unsigned int g = 1; g <= max_group; ++g) {
          double sum_weighted_var = 0.0;
          double sum_weights = 0.0;
          int count = 0;
          std::vector<double> weights;
          std::vector<double> means;
          std::vector<unsigned int> indices;
          for (unsigned int i = 0; i < omegaShare.n_elem; ++i) {
            if (omegaShare(i) == g) {
              double w = 1.0;
              if (nMix > 1 && omegaShareSubpop.n_elem == omegaShare.n_elem) {
                unsigned int subpop = omegaShareSubpop(i);
                if (subpop >= 1 && subpop <= (unsigned int)nMix) {
                  w = arma::sum(mixWeights.col(subpop - 1));
                }
              }
              weights.push_back(w);
              double mu = 0.0;
              if (mprior_phi1.n_rows > 0) {
                mu = arma::mean(mprior_phi1.col(i));
              }
              means.push_back(mu);
              indices.push_back(i);
              sum_weighted_var += w * Gamma2_phi1(i, i);
              sum_weights += w;
              count++;
            }
          }
          if (count > 1 && sum_weights > 0.0) {
            double mean_of_vars = sum_weighted_var / sum_weights;
            double mu_total = 0.0;
            for (size_t k = 0; k < weights.size(); ++k) {
              mu_total += weights[k] * means[k];
            }
            mu_total /= sum_weights;
            double var_of_means = 0.0;
            for (size_t k = 0; k < weights.size(); ++k) {
              double diff = means[k] - mu_total;
              var_of_means += weights[k] * diff * diff;
            }
            var_of_means /= sum_weights;
            double total_var = mean_of_vars + var_of_means;
            for (unsigned int i : indices) {
              Gamma2_phi1Report(i, i) = total_var;
            }
          }
        }
      }
      // "msaem" split-ETA columns: generic Gmin/minv floor (1e-20) isn't tight enough to stop
      // IGamma2_phi1 exploding and locking MCMC proposals to zero; floor at a fraction of ini() variance instead.
      if (nMix > 1 && mixSampleMethod == 1 && omegaShareSubpop.n_elem == (unsigned int)nphi1) {
        for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
          if (omegaShareSubpop(c) < 1) continue;
          double floorVar = 0.1 * Gamma2_phi1Init(c, c);
          if (std::isfinite(floorVar) && floorVar > 0 && Gamma2_phi1(c, c) < floorVar) {
            Gamma2_phi1(c, c) = floorVar;
          }
        }
      }
      vec Gmin=minv(i1);
      uvec jDmin=find(Gamma2_phi1.diag()<Gmin);
      for(unsigned int jm=0; jm<jDmin.n_elem; jm++) {
        Gamma2_phi1(jDmin(jm),jDmin(jm))=Gmin(jDmin(jm));
      }
      // fix before diagonals are enforced
      if (Gamma2_phi1fixed==1 && kiter > (unsigned int)(nb_fixOmega)) {
        Gamma2_phi1.elem(Gamma2_phi1fixedIx) = Gamma2_phi1fixedValues(Gamma2_phi1fixedIx);
      }

      if (kiter<=(unsigned int)(nb_correl)) {
        Gamma2_phi1 = diagmat(Gamma2_phi1);
      }

      if (nphi0>0) {
        if (kiter<=(unsigned int)(niter_phi0)) {
          // omega estimation
          Gamma2_phi0=(statphi02 + mprior_phi0.t()*mprior_phi0 - statphi01.t()*mprior_phi0 - mprior_phi0.t()*statphi01)/N;
          Gmin=minv(i0);
          jDmin=find(Gamma2_phi0.diag()<Gmin);
          for(unsigned int jm=0; jm<jDmin.n_elem; jm++) {
            Gamma2_phi0(jDmin(jm),jDmin(jm))=Gmin(jDmin(jm));
          }
          dGamma2_phi0=Gamma2_phi0.diag();
        } else {
          dGamma2_phi0=dGamma2_phi0*coef_phi0;
        }
        Gamma2_phi0=diagmat(dGamma2_phi0);                         //CHK
      }
      //CHECK the following seg on b & yptr & fptr
      // general log-likelihood (distribution==4): no residual error params to update
      if (distribution != 4)
      for(int b=0; b<nendpnt; ++b) {
        // AR(1): update the correlation from this iteration's residual pairs
        // (grid-search profile likelihood + stochastic approximation).  statrese
        // already holds the whitened SSR, so sig2 below is the marginal variance.
        if (arActive(b)) arUpdateCor(b, kiter, pas);
        double sig2=statrese[b]/(y_offset(b+1)-y_offset(b));       //CHK: range
        int offsetR = res_offset[b];
        _saemFixedIdx[0] = _saemFixedIdx[1] = _saemFixedIdx[2] = _saemFixedIdx[3] = 0;
        switch (res_mod(b)) {
        case rmAdd:
          {
            if (resFixed[offsetR] == 1 && kiter > (unsigned int)(nb_fixResid)) {
              ares(b) = resValue[offsetR];
            } else {
              ares(b) = sqrt(sig2);
            }
          }
          break;
        case rmProp:
          {
            if (resFixed[offsetR] == 1 && kiter > (unsigned int)(nb_fixResid)) {
              bres(b) = resValue[offsetR];
            } else {
              if (sig2 == 0) sig2 = 1;
              bres(b) = sqrt(sig2);
            }
          }
          break;
        case rmAddProp:
          {
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            // yptr = ysb.memptr();
            // fptr = fsb.memptr();
            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(2);
            double *pxmin = xmin.memptr();
            int n=2;
            double start[2]={sqrt(fabs(ares(b))), sqrt(fabs(bres(b)))};                  //force are & bres to be positive
            double step[2]={-.2, -.2};

            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                ares(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(ares(b)));
              } else {
                start[curi++] = sqrt(fabs(ares(b)));
                n++;
              }
              if (resFixed[offsetR + 1] == 1) {
                bres(b) = resValue[offsetR + 1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
            }

            // f = sum((ytr-ft)/g);
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp=addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemFn = obj;
            _saemStep = step;
            _saemStart=start;
            _saemOpt(n, pxmin);
            // Adjust back
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 1] == 0) {
                double ab12 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab12*ab12 - bres(b));    //force are & bres to be positive
              }
            } else {
              double ab02 = pxmin[0];
              double ab12 = pxmin[1];
              ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              bres(b) = bres(b) + pas(kiter)*(ab12*ab12 - bres(b));    //force are & bres to be positive
            }
          }
          break;
        case rmAddPow:
          { // add + pow
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            // yptr = ysb.memptr();
            // fptr = fsb.memptr();
            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(3);
            double *pxmin = xmin.memptr();
            int n=3;

            // REprintf("ares: %f bres: %f cres: %f\n", ares(b), bres(b), cres(b));
            double start[3]={sqrt(fabs(ares(b))), sqrt(fabs(bres(b))), toPowEst(cres(b))}; //force are & bres to be positive
            double step[3]={-.2, -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                ares(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(ares(b)));
              } else {
                start[curi++] = sqrt(fabs(ares(b)));
                n++;
              }
              if (resFixed[offsetR + 1] == 1) {
                bres(b) = resValue[offsetR + 1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR + 2] == 1) {
                cres(b) = resValue[offsetR + 2];
                _saemFixedIdx[2] = 1;
                _saemFixedValue[2] = toPowEst(cres(b));
              } else {
                start[curi++] = toPowEst(cres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objC;
            _saemOpt(n, pxmin);
            // REprintf("\tares: %f bres: %f cres: %f\n", pxmin[0], pxmin[1], pxmin[2]);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 1] == 0) {
                double ab12 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab12*ab12 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 2] == 0) {
                cres(b) = cres(b) + pas(kiter)*(toPow(pxmin[curi++]) - cres(b));    //force are & bres to be positive
              }
            } else {
              ares(b) = ares(b) + pas(kiter)*(pxmin[0]*pxmin[0] - ares(b)); //force ares & bres to be positive
              bres(b) = bres(b) + pas(kiter)*(pxmin[1]*pxmin[1] - bres(b)); //force ares & bres to be positive
              cres(b) = cres(b) + pas(kiter)*(toPow(pxmin[2]) - cres(b));
            }
          }
          break;
        case rmPow:
          { // power
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(2);
            double *pxmin = xmin.memptr();
            int n=2;
            double start[2]={sqrt(fabs(bres(b))), toPowEst(cres(b))};                  //force are & bres to be positive
            double step[2]={ -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                bres(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR + 1] == 1) {
                cres(b) = resValue[offsetR + 1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = toPowEst(cres(b));
              } else {
                start[curi++] = toPowEst(cres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp =addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objD;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab02*ab02 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 1] == 0) {
                double ab12 = pxmin[curi++];
                cres(b) = cres(b) + pas(kiter)*(toPow(ab12) - cres(b));    //force are & bres to be positive
              }
            } else {
              bres(b) = bres(b) + pas(kiter)*(pxmin[0]*pxmin[0] - bres(b));    //force are & bres to be positive
              cres(b) = cres(b) + pas(kiter)*(toPow(pxmin[1]) - cres(b));      //force are & bres to be positive
            }
          }
          break;
        case rmAddLam:
          { // additive + lambda
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(2);
            double *pxmin = xmin.memptr();
            int n=2;
            double start[2]={sqrt(fabs(ares(b))), toLambdaEst(lres(b))};                  //force are & bres to be positive
            double step[2]={ -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                ares(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(ares(b)));
              } else {
                start[curi++] = sqrt(fabs(ares(b)));
                n++;
              }
              if (resFixed[offsetR + 1] == 1) {
                lres(b) = resValue[offsetR + 1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = toLambdaEst(lres(b));
              } else {
                start[curi++] = toLambdaEst(lres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objE;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 1] == 0) {
                double ab12 = pxmin[curi++];
                lres(b) = lres(b) + pas(kiter)*(toLambda(ab12) - lres(b));    //force are & bres to be positive
              }
            } else {
              ares(b) = ares(b) + pas(kiter)*(pxmin[0]*pxmin[0] - ares(b));    //force are & bres to be positive
              lres(b) = lres(b) + pas(kiter)*(toLambda(pxmin[1]) - lres(b));   //force are & bres to be positive
            }
            lambda(b) = lres(b);
          }
          break;
        case rmPropLam:
          { // prop + lambda
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(2);
            double *pxmin = xmin.memptr();
            int n=2;
            double start[2]={sqrt(fabs(bres(b))), toLambdaEst(lres(b))};                  //force are & bres to be positive
            double step[2]={ -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                bres(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR + 1] == 1) {
                lres(b) = resValue[offsetR + 1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = toLambdaEst(lres(b));
              } else {
                start[curi++] = toLambdaEst(lres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objF;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab02*ab02 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 1] == 0) {
                double ab12 = pxmin[curi++];
                lres(b) = lres(b) + pas(kiter)*(toLambda(ab12) - lres(b));    //force are & bres to be positive
              }
            } else {
              bres(b) = bres(b) + pas(kiter)*(pxmin[0]*pxmin[0] - bres(b));    //force are & bres to be positive
              lres(b) = lres(b) + pas(kiter)*(toLambda(pxmin[1]) - lres(b));            //force are & bres to be positive
            }
            lambda(b) = lres(b);
          }
          break;
        case rmPowLam:
          { // pow + lambda
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(3);
            double *pxmin = xmin.memptr();
            int n=3;
            double start[3]={sqrt(fabs(bres(b))), toPowEst(cres(b)), toLambdaEst(lres(b))};                  //force are & bres to be positive
            double step[3]={ -.2, -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                bres(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR+1] == 1) {
                cres(b) = resValue[offsetR+1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = toPowEst(cres(b));
              } else {
                start[curi++] = toPowEst(cres(b));
                n++;
              }
              if (resFixed[offsetR + 2] == 1) {
                lres(b) = resValue[offsetR + 2];
                _saemFixedIdx[2] = 1;
                _saemFixedValue[2] = toLambdaEst(lres(b));
              } else {
                start[curi++] = toLambdaEst(lres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objG;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab02*ab02 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR+1] == 0) {
                double ab02 = pxmin[curi++];
                cres(b) = cres(b) + pas(kiter)*(toPow(ab02) - cres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 2] == 0) {
                double ab12 = pxmin[curi++];
                lres(b) = lres(b) + pas(kiter)*(toLambda(ab12) - lres(b));    //force are & bres to be positive
              }
            } else {
              bres(b) = bres(b) + pas(kiter)*(pxmin[0]*pxmin[0] - bres(b));    //force are & bres to be positive
              cres(b) = cres(b) + pas(kiter)*(toPow(pxmin[1]) - cres(b));    //force are & bres to be positive
              lres(b) = lres(b) + pas(kiter)*(toLambda(pxmin[2]) - lres(b));            //force are & bres to be positive
            }
            lambda(b) = lres(b);
          }
          break;
        case rmAddPropLam:
          { // add + prop + lambda
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(3);
            double *pxmin = xmin.memptr();
            int n=3;
            double start[3]={sqrt(fabs(ares(b))), sqrt(fabs(bres(b))), toLambdaEst(lres(b))};                  //force are & bres to be positive
            double step[3]={ -.2, -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                ares(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(ares(b)));
              } else {
                start[curi++] = sqrt(fabs(ares(b)));
                n++;
              }
              if (resFixed[offsetR+1] == 1) {
                bres(b) = resValue[offsetR+1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR + 2] == 1) {
                lres(b) = resValue[offsetR + 2];
                _saemFixedIdx[2] = 1;
                _saemFixedValue[2] = toLambdaEst(lres(b));
              } else {
                start[curi++] = toLambdaEst(lres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objH;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR+1] == 0) {
                double ab02 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab02*ab02 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 2] == 0) {
                double ab12 = pxmin[curi++];
                lres(b) = lres(b) + pas(kiter)*(toLambda(ab12) - lres(b));    //force are & bres to be positive
              }
            } else {
              ares(b) = ares(b) + pas(kiter)*(pxmin[0]*pxmin[0] - ares(b));    //force are & bres to be positive
              bres(b) = bres(b) + pas(kiter)*(pxmin[1]*pxmin[1] - bres(b));    //force are & bres to be positive
              lres(b) = lres(b) + pas(kiter)*(toLambda(pxmin[2]) - lres(b));            //force are & bres to be positive
            }
            lambda(b) = lres(b);
          }
          break;
        case rmAddPowLam:
          { // add + pow + lambda
            uvec idx;
            idx = find(ix_endpnt==b);
            vec ysb, fsb;

            buildFsbYsb(idx, fsM, fsb, ysb);

            //len = ysb.n_elem;                                        //CHK: needed by nelder
            vec xmin(4);
            double *pxmin = xmin.memptr();
            int n=4;
            double start[4]={sqrt(fabs(ares(b))), sqrt(fabs(bres(b))), toPowEst(cres(b)), toLambdaEst(lres(b))};                  //force are & bres to be positive
            double step[4]={ -.2, -.2, -.2, -.2};
            if (kiter > (unsigned int)(nb_fixResid)) {
              n = 0;
              int curi=0;
              if (resFixed[offsetR] == 1) {
                ares(b) = resValue[offsetR];
                _saemFixedIdx[0] = 1;
                _saemFixedValue[0] = sqrt(fabs(ares(b)));
              } else {
                start[curi++] = sqrt(fabs(ares(b)));
                n++;
              }
              if (resFixed[offsetR+1] == 1) {
                bres(b) = resValue[offsetR+1];
                _saemFixedIdx[1] = 1;
                _saemFixedValue[1] = sqrt(fabs(bres(b)));
              } else {
                start[curi++] = sqrt(fabs(bres(b)));
                n++;
              }
              if (resFixed[offsetR+2] == 1) {
                cres(b) = resValue[offsetR+2];
                _saemFixedIdx[2] = 1;
                _saemFixedValue[2] = toPowEst(cres(b));
              } else {
                start[curi++] = toPowEst(cres(b));
                n++;
              }
              if (resFixed[offsetR + 3] == 1) {
                lres(b) = resValue[offsetR + 3];
                _saemFixedIdx[3] = 1;
                _saemFixedValue[3] = toLambdaEst(lres(b));
              } else {
                start[curi++] = toLambdaEst(lres(b));
                n++;
              }
            }
            _saemYptr = ysb.memptr();
            _saemFptr = fsb.memptr();
            _saemLen  = ysb.n_elem;
            _saemYj   = yj(b);
            _saemPropT = propT(b);
            _saemAddProp = addProp(b);
            _saemLambda = lambda(b);
            _saemLow = low(b);
            _saemHi = hi(b);
            _saemStep = step;
            _saemStart = start;
            _saemFn = objI;
            _saemOpt(n, pxmin);
            if (kiter > (unsigned int)(nb_fixResid)) {
              int curi = 0;
              if (resFixed[offsetR] == 0) {
                double ab02 = pxmin[curi++];
                ares(b) = ares(b) + pas(kiter)*(ab02*ab02 - ares(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR+1] == 0) {
                double ab02 = pxmin[curi++];
                bres(b) = bres(b) + pas(kiter)*(ab02*ab02 - bres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 2] == 0) {
                double ab12 = pxmin[curi++];
                cres(b) = cres(b) + pas(kiter)*(toPow(ab12) - cres(b));    //force are & bres to be positive
              }
              if (resFixed[offsetR + 3] == 0) {
                double ab12 = pxmin[curi++];
                lres(b) = lres(b) + pas(kiter)*(toLambda(ab12) - lres(b));    //force are & bres to be positive
              }
            } else {
              ares(b) = ares(b) + pas(kiter)*(pxmin[0]*pxmin[0] - ares(b));    //force are & bres to be positive
              bres(b) = bres(b) + pas(kiter)*(pxmin[1]*pxmin[1] - bres(b));    //force are & bres to be positive
              cres(b) = cres(b) + pas(kiter)*(toPow(pxmin[2]) - cres(b));    //force are & bres to be positive
              lres(b) = lres(b) + pas(kiter)*(toLambda(pxmin[3]) - lres(b));            //force are & bres to be positive
            }
            lambda(b) = lres(b);
          }
          break;
        }
        sigma2[b] = sig2;                                          //CHK: sigma2[] use
        if (sigma2[b]>1.0e99) sigma2[b] = 1.0e99;
        if (std::isnan(sigma2[b])) sigma2[b] = 1.0e99;
      }
      vecares = ares(ix_endpnt);
      vecbres = bres(ix_endpnt);
      veccres = cres(ix_endpnt);
      if (DEBUG>0) Rcout << "par update successful\n";

      //    Fisher information
      DDa=(D1/nmc)*(D1/nmc).t()-D11/nmc-D2/nmc;
      DDb=-D11/nmc-D2/nmc;
      // d2logk omits the deterministic mu-block complete Hessian (-M'Omega^{-1}M), leaving
      // DDa's fixed-effect block equal to -Var[score] (indefinite -> solve(Ha) yields NaN
      // SEs for covMethod="fim").  Add that block (= CGamma2, the mu Fisher information the
      // M-step forms) so DDa = -E[Hessian] - Var[score] is the true observed information.
      {
        unsigned int nl1 = CGamma21.n_rows, nl0 = CGamma20.n_rows;
        if (nl1 > 0) DDa.submat(0, 0, nl1 - 1, nl1 - 1) += CGamma21;
        if (nl0 > 0) DDa.submat(nl1, nl1, nl1 + nl0 - 1, nl1 + nl0 - 1) += CGamma20;
      }
      L=L+pash(kiter)*(D1/nmc-L);
      Ha=Ha+pash(kiter)*(DDa- Ha);
      Hb=Hb+pash(kiter)*(DDb- Hb);
      // SA covariance phase (covMethod="sa"): theta is frozen (pas==pash==0 above), so the
      // now-corrected DDa is a Monte-Carlo draw of the observed information at theta_hat.
      // Monte-Carlo average it (after a short burn-in) into HaSa; the covariance is
      // solve(HaSa).
      if (nSaCov > 0 && kiter >= (unsigned int)niter) {
        unsigned int cc = kiter - (unsigned int)niter;
        if (cc >= (unsigned int)(0.1 * nSaCov)) {
          covCount++;
          HaSa += (DDa - HaSa) / (double)covCount;
        }
      }
      cube phi2 = phi%phi;
      mat sphi1 = sum(phi ,2);
      mat sphi2 = sum(phi2,2);
      mpost_phi=mpost_phi+pash(kiter)*(sphi1/nmc-mpost_phi);
      cpost_phi=cpost_phi+pash(kiter)*(sphi2/nmc-cpost_phi);
      mpost_phi.cols(i0)=mprior_phi0;

      //FIXME: chg according to multiple endpnts; need to chg dim(par_hist)
      for (int b=0; b<nendpnt; ++b) {
        int offset = res_offset[b];
        switch ((int)(res_mod(b))) {
        case rmAdd:
          vcsig2[offset] = ares(b);//sigma2[b];
          // because of old translation use variance
          break;
        case rmProp:
          vcsig2[offset] = bres(b);
          break;
        case rmPow:
          vcsig2[offset]     = bres(b);
          vcsig2[offset + 1] = cres(b);
          break;
        case rmAddProp:
          vcsig2[offset]   = ares(b);
          vcsig2[offset+1] = bres(b);
          break;
        case rmAddPow:
          vcsig2[offset]   = ares(b);
          vcsig2[offset+1] = bres(b);
          vcsig2[offset+2] = cres(b);
          break;
        case rmAddLam:
          vcsig2[offset]   = ares(b);
          vcsig2[offset+1] = lres(b);
          break;
        case rmPropLam:
          vcsig2[offset]   = bres(b);
          vcsig2[offset+1] = lres(b);
          break;
        case rmPowLam:
          vcsig2[offset]   = bres(b);
          vcsig2[offset+1] = cres(b);
          vcsig2[offset+2] = lres(b);
          break;
        case rmAddPropLam:
          vcsig2[offset]   = ares(b);
          vcsig2[offset+1] = bres(b);
          vcsig2[offset+2] = lres(b);
          break;
        case rmAddPowLam:
          vcsig2[offset]   = ares(b);
          vcsig2[offset+1] = bres(b);
          vcsig2[offset+2] = cres(b);
          vcsig2[offset+3] = lres(b);
          break;
        }
      }
      Plambda(ilambda1) = Plambda1;
      Plambda(ilambda0) = Plambda0;
      if (nMix > 1 && mixProbRegress) {
        // Fixed membership: the mixing proportions are just the (constant) hard
        // assignment fractions; no soft-EM SA update.
        mixProb = mean(mixWeights, 0).t();
      } else if (nMix > 1) {
        vec mean_aji = mean(mixWeights, 0).t();
        if (mixProbMethod == 1) {
          // Dirichlet-style regularization: blend in mixProbPriorN pseudo-subjects from the
          // initial mixing distribution before the SA step, damping the responsibility average
          // even during burn-in's flat pas(kiter)==1.
          vec mean_aji_reg = (mean_aji * N + mixProbPriorN * mixProbInit) / (N + mixProbPriorN);
          // "msaem" separation develops more slowly than "parallel", so stay at full-replacement
          // step for 80% of the fit to give mixProb time to track the slower signal.
          double pasMsaem = (mixSampleMethod == 1 && kiter < (unsigned int)(0.8 * niter)) ? 1.0 : pas(kiter);
          mixProb = mixProb + pasMsaem * (mean_aji_reg - mixProb);
        } else {
          // Annealed step-size: mixProb gets its own decaying schedule instead of the shared
          // burn-in pas(kiter)==1, which would let one noisy iteration fully replace it (a
          // runaway feedback loop, since a smaller mixProb(k) shrinks its own responsibility further).
          mixProb = mixProb + pasMix(kiter) * (mean_aji - mixProb);
        }
      }
      vec pl = Plambda.elem(parHistThetaKeep);
      vec g2 = Gamma2_phi1.diag();
      g2 = g2.elem(parHistOmegaKeep);
      pl = join_cols(pl, g2);
      if (parHistOmegaOffPairs.n_rows > 0) {                    // off-diagonal Omega covariances
        vec offv(parHistOmegaOffPairs.n_rows);
        for (unsigned int p = 0; p < parHistOmegaOffPairs.n_rows; ++p) {
          offv(p) = Gamma2_phi1(parHistOmegaOffPairs(p, 0), parHistOmegaOffPairs(p, 1));
        }
        pl = join_cols(pl, offv);
      }
      g2 = vcsig2.elem(resKeep);
      pl = join_cols(pl, g2);
      if (nMix > 1) {
        vec mixP = mixProb.head(nMix - 1);
        pl = join_cols(pl, mixP);
      }
      if (kiter < (unsigned int)niter) {
        par_hist.row(kiter) = pl.t();
        // saem has no per-iteration objective function; scale.showOfv=0 so the
        // `f` argument is ignored here.  scalePrintFun gates printing on
        // (cn % every == 0) and runs the user-interrupt check internally.
        // The X row's tag carries the phase: "SA: X" (burn) / "EM: X".
        if (nPhase1 >= 0) {
          scale.phaseLabel = (kiter < (unsigned int)nPhase1) ? "SA" : "EM";
        }
        scalePrintFun(&scale, pl.memptr(), NA_REAL);
        // One summary line per printed iteration for any ODE-solve warnings
        // accumulated since the last flush (see inner.cpp foceiOfvOptim).
        nmFlushRxSolveWarn(5);
      }
      // SA covariance phase (kiter >= niter): theta is frozen and HaSa is accumulated
      // above; nothing is recorded to par_hist and no printing happens.
    }//kiter
    // restore the converged estimate after the SA covariance phase (the reported fit
    // must be the converged value, not a cov-phase iterate)
    if (nSaCov > 0) {
      Plambda = _savPlambda; Gamma2_phi1 = _savGamma2_phi1; Gamma2_phi0 = _savGamma2_phi0;
      Gamma2_phi1Report = _savGamma2_phi1Report; mprior_phi1 = _savMprior_phi1;
      mprior_phi0 = _savMprior_phi0; ares = _savAres; bres = _savBres; cres = _savCres;
      lres = _savLres; vcsig2 = _savVcsig2; phiM = _savPhiM; Ha = _savHa;
      lambda = lres;
      if (nMix > 1) { mixProb = _savMixProb; mixWeights = _savMixWeights; }
    }
    phiFile.close();
  }


private:

  user_funct user_fn;

  uvec nu;
  int niter;
  int saemSeed = 99;
  int nPhase1;
  // uninformative-eta revisit: re-run the informativeness test at the end of burn-in
  // (see revisitUninformativeEtas).  ueRevisitIter < 0 disables it.
  int ueRevisitIter = -1;
  uvec ueRevisitCols;
  vec ueDelta;                            // probe half-width per revisited column
  double ueTol = 1e-7;
  // diagnostics: whether the revisit ran, and how many (subject, eta) verdicts it
  // changed in each direction.  Tests assert on these -- a mask that is unchanged is
  // indistinguishable from a revisit that never happened.
  int ueRevisitRan = 0, ueRevisitUnfroze = 0, ueRevisitFroze = 0;
  int nb_sa;
  int nb_correl;
  int nb_fixOmega;
  int nb_fixResid;
  int niter_phi0;
  double coef_phi0;
  vec phi0Lower;
  vec phi0Upper;
  double rmcmc;
  double coef_sa;
  vec pas, pash;
  vec minv;
  int nmc;
  int nM;

  int ntotal, N, mlen;
  vec y, ys;    //ys is y sorted by endpnt
  mat evt;
  mat phiM;
  uvec indio;
  mat Mcovariables;
  List opt, optM;

  int nphi0, nphi1, nphi;
  mat covstruct1;
  uvec i1, i0, fixedIx1, fixedIx0;
  umat Gamma2_phi1fixedIxIn;
  uvec Gamma2_phi1fixedIx;
  int Gamma2_phi1fixed;
  mat Gamma2_phi1fixedValues;
  uvec pc1;
  mat COV1, COV0, LCOV1, LCOV0, COV21, COV20, MCOV1, MCOV0;
  mat Gamma2_phi1, Gamma2_phi0, mprior_phi1, mprior_phi0;
  mat Gamma2_phi1Report; // reporting-only pooled BSV for split ETAs sharing an omegaShare group; never fed back into estimation
  mat Gamma2_phi1Init; // initial (ini()) Gamma2_phi1, captured once; used as an msaem-only exploration floor for split-ETA columns (see saem_fit())
  mat IGamma2_phi1, D1Gamma21, D2Gamma21, CGamma21;
  mat IGamma2_phi0, D1Gamma20, D2Gamma20, CGamma20;
  mat Gamma_phi1, Gdiag_phi1, Gamma_phi0, Gdiag_phi0;
  vec gamma2_phi1, gamma2_phi0;
  uvec ind_cov1, ind_cov0, jcov1, jcov0;
  vec dGamma2_phi0;
  vec Plambda;

  int nlambda1, nlambda0, nlambda, nb_param;
  uvec ilambda1, ilambda0;
  // one FIM residual slot per endpoint (see the nb_param comment in inits());
  // residEpIdx[b] is that endpoint's compacted slot, or -1 when the whole
  // model is a general log-likelihood (nResidEp==0)
  int nResidEp;
  int residEpIdx[MAXENDPNT];

  mat statphi01, statphi02, statphi11, statphi12;
  // Per-component, unblended sufficient statistic (never mixed across components); used to fix
  // tcl1/tcl2-style split-ETA fixed effects via weighted regression (see omegaShareSubpop block).
  field<mat> statphi11_mix;
  double statrese[MAXENDPNT];
  double sigma2[MAXENDPNT];
  vec ares, bres, cres, lres, lambda, low, hi;
  vec vecares, vecbres, veccres, veclres;
  uvec res_mod, yj, propT, addProp, vecaddProp;

  mat DYF;
  cube phi;
  bool hasFixedObsTransform = false;
  vec yTrans;
  vec ysTrans;

  vec L;
  mat Ha, Hb, DDa, DDb;
  mat mpost_phi, cpost_phi;

  vec resValue;
  uvec resFixed;
  uvec resKeep;

  mcmcaux mx;

  // Iteration-print formatting shared with focei/nlm via src/scale.h (scaleTypeNone,
  // since Plambda is already on the model scale); transform fields set by scaleAttachXform.
  scaling scale;
  std::vector<double> scaleInitPar;
  std::vector<double> scaleC;
  CharacterVector scaleNames;
  mat par_hist;
  uvec parHistThetaKeep;
  uvec parHistOmegaKeep;
  umat parHistOmegaOffPairs;

  // SA (stochastic-approximation) covariance phase
  int nSaCov;
  mat HaSa;                      // converged Louis FIM (nb_param x nb_param) at theta_hat
  int covCount;                  // averaged cov-phase iterations
  // saved converged state, restored after the cov phase so the fit is unperturbed
  vec _savPlambda, _savAres, _savBres, _savCres, _savLres, _savVcsig2, _savMixProb;
  mat _savGamma2_phi1, _savGamma2_phi0, _savGamma2_phi1Report, _savMprior_phi1, _savMprior_phi0, _savPhiM, _savHa, _savMixWeights;

  int distribution;
  // nonMuTheta="regress": estimate the fixed-effect-only (phi0) parameters by
  // the bounded direct optimizer (refinePhi0Lik) for NORMAL models too, instead
  // of only the stochastic phi0 block.  Keeps such thetas as plain regressors
  // with directly-optimized, bound-respecting values (no shrinking phi0
  // variance).  0 = classic phi0, 1 = direct-optimize.
  int nonMuThetaRegress;
  // Cost controls for that refinement (it is the dominant per-iteration cost when
  // the phi0 thetas drive the ODE): nonMuThetaOptType picks the optimizer
  // (0 = coordinate descent, 1 = nelder-mead, 2 = newuoa -- the latter two over
  // all free coordinates at once), nonMuThetaMaxEval is their objective-
  // evaluation budget (0 means 10 per free coordinate), nonMuThetaSweeps the
  // coordinate-descent sweep count,
  // nonMuThetaTol the inner convergence tolerance and nonMuThetaEvery how often
  // (in iterations) the refinement runs at all.
  int nonMuThetaOptType;
  int nonMuThetaMaxEval;
  int nonMuThetaSweeps;
  int nonMuThetaEvery;
  double nonMuThetaTol;
  // Cached: does any phi0 param change the structural prediction f?  -1 unknown,
  // 0 no (residual/likelihood only -> freeze the ODE during the phi0 opt like
  // npag), 1 yes (structural, e.g. ka/V -> keep the ODE live).
  int _phi0OdeSensitive = -1;
  // Same question for a general-likelihood fit, answered by comparing a frozen
  // and a live evaluation rather than by watching the prediction move: there the
  // prediction IS the log-density, so every phi0 changes it.  -1 not yet probed.
  int _phi0NeedsLive = -1;
  // Phase 4 (SAEM general-likelihood theta plan): the THETA[k]/ETA[k] -> phi
  // column maps, pool-readiness flags, and lhs offsets (_saemPhi1H2ThetaKind,
  // _saemPhi1H2ThetaCol, _saemPhi1H2ThetaFixedVal, _saemPhi1H2EtaCol, _saemPhi1PredOffset,
  // _saemPhi1H2PredOffset, _saemPhi1H2HessOffset, _saemPhi1PoolReady, _saemPhi1UseAnalyticHess)
  // are FILE-SCOPE GLOBALS, not members -- user_function (a free function,
  // set via set_fn/user_fn, not a class method) needs them too, and setupRx
  // (also free) is what populates them from opt.  See their definitions
  // and doc comment near _saemPhi1PoolActive, after this class body.
  int phi1ThetaEvery = 1;
  int phi1ThetaMaxEval = 50;
  // Warm-start residual params from observed per-endpoint moments (npag-style).
  int residWarmStart;
  // mixProbMethod="regress": fix per-subject mixture membership (hard classify
  // once) instead of the soft-EM responsibility step.  mixFixedAssign holds the
  // 1-based assigned component per subject.
  int mixProbRegress;
  uvec mixFixedAssign;

  int nendpnt;
  uvec ix_endpnt;
  umat ix_idM;
  uvec y_offset;
  uvec res_offset;
  vec vcsig2;
  int nres;
  uvec ix_sorting;
  mat fsaveMat;
  vec cens;
  vec limit;
  vec fsave;

  int nMix;
  vec mixProb;
  vec mixProbInit;
  int mixProbMethod; // 0 = annealed step-size, 1 = Dirichlet-style regularization
  vec pasMix;
  double mixProbPriorN;
  int mixSampleMethod; // 0 = parallel per-component chains (NONMEM-style), 1 = MSAEM (Lavielle & Mbogning 2014)
  mat mixWeights;
  mat priorWeights; // leaspy-style prior-only responsibility, msaem only (see saem_fit()'s E-step)
  field<mat> phiM_mix;
  field<vec> fsave_mix;
  field<vec> limit_mix;
  field<vec> cens_mix;
  uvec omegaShare;
  uvec omegaShareSubpop;
  // Two-level (IOV) models: phi1 columns sharing a non-zero group id are one
  // occasion parameter observed at different levels, so they estimate ONE
  // variance (Psi).  0 means the column has its own.  See R/saemIov.R.
  uvec omegaPool;

  // Per-chain scratch buffers pre-allocated in inits() to avoid repeated heap
  // allocation in the hot distribution==1 loops in saem_fit() and do_mcmc().
  vec _scratch_ft;      // transformed-f per chain (replaces ftk/fck)
  vec _scratch_limitT;  // transformed-limit per chain (replaces limitTk)
  vec _scratch_ftT;     // handleF output per chain (replaces ftTk/fcTk)
  vec _scratch_g;       // residual SD per chain (replaces gk/gck)
  uvec _scratch_indio;  // DYF row indices per chain (replaces indio_k)
  vec _scratch_ftAr;    // AR(1)-conditional prediction, filled by arDYFhyp()
  vec _scratch_gAr;     // AR(1)-conditional SD, filled by arDYFhyp()

  uvec obs_subject;

  // AR(1) autocorrelated residuals (continuous-time; phi = cor^dt).  arActive[b]
  // flags an endpoint with ar(); arCor[b] is its correlation (stochastic-approx
  // updated in the M-step).  Per-observation (ORIGINAL order, one chain) arPrev
  // is the previous same-subject-same-endpoint observation index in time order
  // (-1 = first of its subject/endpoint) and arDt the time gap to it.
  vec arCor;
  uvec arActive;
  arma::ivec arPrev;
  vec arDt;
  int hasAr;
  // M-step scratch: residual indexed by original obs (_arRorig), and per-endpoint
  // (r_i, r_prev, dt, weight) pairs + first-obs SSR accumulated over MCMC chains,
  // used by the AR(1) correlation grid-search M-step (arUpdateCor).
  vec _arRorig;
  std::vector<double> arPairR[MAXENDPNT], arPairP[MAXENDPNT], arPairDt[MAXENDPNT], arPairW[MAXENDPNT];
  double arFirstSSR[MAXENDPNT];
  int arNobs[MAXENDPNT];

  int DEBUG;
  std::vector< std::string > phiMFile;

  // Build fsb (predictions across all nmc chains) and ysb (observations
  // replicated nmc times) for a single endpoint b, used by the residual
  // parameter estimation switch cases.
  void buildFsbYsb(const uvec &idx, const vec &fsM, vec &fsb, vec &ysb) const {
    int nb_b = (int)idx.n_elem;
    uvec fsb_idx((arma::uword)(nmc * nb_b));
    for (int k = 0; k < nmc; k++) {
      fsb_idx.subvec((arma::uword)(k * nb_b), (arma::uword)((k + 1) * nb_b - 1)) =
        idx + (arma::uword)(k * ntotal);
    }
    fsb = fsM(fsb_idx);
    ysb = arma::repmat(ys(idx), (arma::uword)nmc, 1);
  }

  // Invert a symmetric covariance (omega) matrix.  If it is not positive
  // definite, project it to the nearest positive-definite matrix (in place, so
  // downstream chol()/set_mcmcphi() see the corrected matrix), warn the user
  // once, and return the inverse of the corrected matrix.
  bool _nearPdWarned = false;
  mat invSympdNearPd(mat &G, const char *what) {
    mat out;
    if (inv_sympd(out, G)) return out;
    mat pd;
    if (nmNearPD(pd, G)) {
      G = pd;
      if (!_nearPdWarned) {
        Rcpp::warning(std::string("SAEM: ") + what +
                      " was not positive definite; projected to the nearest positive-definite matrix (results may be affected)");
        _nearPdWarned = true;
      }
      if (inv_sympd(out, G)) return out;
    }
    // last resort: general inverse of the symmetrized matrix
    return inv(0.5 * (G + G.t()));
  }

  void set_mcmcphi(mcmcphi &mphi1,
		   const uvec i1,
		   const int nphi1,
		   const mat Gamma2_phi1,
		   const mat IGamma2_phi1,
		   const mat mprior_phi1) {
    mphi1.i = i1;
    mphi1.nphi = nphi1;
    mphi1.Gamma_phi=chol(Gamma2_phi1);
    mphi1.IGamma2_phi = IGamma2_phi1;
    mphi1.Gdiag_phi.zeros(nphi1, nphi1);
    mphi1.Gdiag_phi.diag() = sqrt(Gamma2_phi1.diag())*rmcmc;
    mphi1.mprior_phiM = repmat(mprior_phi1,nmc,1);
  }

  // do_mcmc's distribution==4 (general-likelihood) branch clamps: a legitimate
  // per-observation log-likelihood is capped at _saemGenLikCeiling (mirroring
  // case 1's own clamp of its scale g to [double_xmin, xmax] -- a genuine
  // log-density is unbounded above as the effective residual scale collapses
  // toward zero), and the shared +1e99 bad/NaN-solve sentinel is inverted to
  // _saemGenLikBadSolvePenalty (a strongly NEGATIVE value) rather than left
  // as +1e99, which DYF=-fck would otherwise reward as an almost-certain MCMC
  // accept.  Values, not derived constants: large enough to swamp any
  // realistic per-observation contribution in either direction.
  static constexpr double _saemGenLikCeiling = 700.0;
  static constexpr double _saemGenLikBadSolvePenalty = -1.0e10;

  void do_mcmc(const int method,
               const int nu,
               const mcmcaux &mx,
               const mcmcphi &mphi,
               mat &DYF,
               mat &phiM,
               vec &U_y,
               vec &U_phi,
               vec &cur_fsave,
               vec &cur_cens,
               vec &cur_limit,
               int kiter,
               int mixIdx = 0) {
    mat fcMat;
    vec fc, fs, Uc_y, Uc_phi, deltu;
    uvec ind;
    arma::vec accU(mx.nM);   // per-block threefry acceptance uniforms

    uvec i=mphi.i;
    double double_xmin = 1.0e-200;                               //FIXME hard-coded xmin, also in neldermean.hpp
    double xmax = 1e300;
    for (int u=0; u<nu; u++)
      for (int k1=0; k1<mphi.nphi; k1++) {
        mat phiMc=phiM;
        // iteration-indexed threefry stream: proposal noise + the acceptance
        // uniform are drawn here (before the solve) from one seeded stream
        _saemSeedDoMcmc((uint32_t)saemSeed, kiter, method, u, k1, mixIdx);
        switch (method) {
        case 1: {
          mat noise(mx.nM, mphi.nphi); _saemFillNormEng(noise);
          phiMc.cols(i)=noise*mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          break;
        }
        case 2: {
          mat noise(mx.nM, mphi.nphi); _saemFillNormEng(noise);
          phiMc.cols(i)=phiM.cols(i) +
            noise*mphi.Gdiag_phi % current_saem_state->_saemUE.cols(i);
          break;
        }
        case 3: {
          vec noise(mx.nM); _saemFillNormEng(noise);
          phiMc.col(i(k1))=phiM.col(i(k1))+
            noise*mphi.Gdiag_phi(k1,k1) % current_saem_state->_saemUE.col(i(k1));
          break;
        }
        }
        _saemFillUnifEng(accU);   // acceptance uniforms from the same stream

        fcMat = nmRngGuard([&]{ return user_fn(phiMc, mx.evtM, mx.optM); });
        cur_limit = fcMat.col(2);
        cur_cens = fcMat.col(1);

        fc = fcMat.col(0);
        fs = fc;
        switch (distribution) {
        case 1:
          {
            // Build yt once: does not depend on chain index k
            vec yt = hasFixedObsTransform ? yTrans : mx.y;
            if (!hasFixedObsTransform) {
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                yt(i) = _powerD(mx.y(i), lambda(cur), yj(cur), low(cur), hi(cur));
              }
            }
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            for (int k = 0; k < nmc; k++) {
              int obs_start = k * ntotal;
              vec fsk = fs.subvec(obs_start, obs_start + ntotal - 1);
              const vec limitk = cur_limit.subvec(obs_start, obs_start + ntotal - 1);
              const vec censk = cur_cens.subvec(obs_start, obs_start + ntotal - 1);
              _scratch_ft = fsk;
              _scratch_limitT = limitk;
              for (int i = ntotal; i--;) {
                int cur = ix_endpnt(i);
                _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ft(i) = _powerD(fsk(i), lambda(cur), yj(cur), low(cur), hi(cur));
                _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fsk(i), false, true);
              }
              saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = mx.indio + (arma::uword)k * stride;
              DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              applyCensLoss(DYF, _scratch_indio, censk, yt, _scratch_limitT,
                            _scratch_ftAr, _scratch_gAr);
            }
          }
          break;
        case 2:
          {
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            for (int k = 0; k < nmc; k++) {
              vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
              _scratch_indio = mx.indio + (arma::uword)k * stride;
              DYF(_scratch_indio) = -mx.y % log(fck) + fck;
            }
          }
          break;
        case 3:
          {
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            for (int k = 0; k < nmc; k++) {
              vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
              _scratch_indio = mx.indio + (arma::uword)k * stride;
              DYF(_scratch_indio) = -mx.y % log(fck) - (1 - mx.y) % log(1 - fck);
            }
          }
          break;
        case 4:
          {
            // general log-likelihood: model prediction is the per-obs loglik
            const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
            // Case 1 (closed-form normal/etc.) clamps its own scale (g) to
            // [double_xmin, xmax] before scoring, specifically so a candidate
            // that collapses the residual scale toward zero cannot blow up
            // that observation's contribution -- log N(x;mu,sigma) -> +Inf as
            // sigma -> 0.  A general-likelihood fck IS the log-likelihood
            // already (no separate scale to clamp), so it needs the same kind
            // of finite ceiling directly, or a scale-collapsing candidate
            // (e.g. driving Vc, and hence cp=centr/Vc, toward 0 in a
            // proportional-error MM model) gets an unboundedly large reward
            // and the acceptance test (deltu < -log(accU)) takes it almost
            // unconditionally.  (Checked directly for the U009 MM corpus
            // model's own divergence: this clamp never actually fires there,
            // so it is not what drives that specific case -- it is still a
            // real, principled gap for models where a candidate genuinely
            // does collapse the residual scale, and mirrors case 1's own
            // clamp on structural grounds, not as an attempted fix for U009.)
            // saemReadRowsPooled also flags a failed/NaN solve with the
            // shared +1e99 sentinel (see its own comment): for every OTHER
            // distribution fc is a predicted VALUE compared against data, so
            // a huge fc naturally penalizes a bad solve via the residual --
            // here fc IS the log-likelihood, so left unrecognized the
            // sentinel flows into DYF=-fck as an unboundedly NEGATIVE (i.e.
            // rewarded) contribution, exactly backwards. Detect and invert it
            // to a strongly penalized value before applying the ceiling.
            for (int k = 0; k < nmc; k++) {
              vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
              fck.elem(find(fck >= 1.0e99)).fill(_saemGenLikBadSolvePenalty);
              fck.elem(find(fck > _saemGenLikCeiling)).fill(_saemGenLikCeiling);
              _scratch_indio = mx.indio + (arma::uword)k * stride;
              DYF(_scratch_indio) = -fck;
            }
          }
          break;
        }

        Uc_y=sum(DYF,0).t();
        if (method==1) {
          deltu=Uc_y-U_y;
        }
        else {
          mat dphic=phiMc.cols(i)-mphi.mprior_phiM;
          Uc_phi=0.5*sum(dphic%(dphic*mphi.IGamma2_phi),1);
          deltu=Uc_y-U_y+Uc_phi-U_phi;
        }

        ind=find( deltu < -log(accU) );
        phiM(ind,i)=phiMc(ind,i);
        U_y(ind)=Uc_y(ind);
        if (method>1) {
          U_phi(ind)=Uc_phi(ind);
        }
        ind = getObsIdx(ix_idM.rows(ind));
        cur_fsave(ind)=fs(ind);
        if (method<3) {
          break;
        }
      }
    setRxThreadId(-1);
  }

  // MSAEM (Lavielle & Mbogning 2014): mixture-weighted (log-sum-exp) observation loss for a
  // candidate phi, evaluated under every component. mix() only lets the component affect which
  // observation columns are read, so only this loss needs mixture treatment; the prior penalty
  // (in do_mcmc_msaem()) stays a standard single-Gaussian quadratic form.
  vec mixObsLoss(const mat &phiC, const mcmcaux &mx) {
    double double_xmin = 1.0e-200;
    double xmax = 1e300;
    mat lossByM(mx.nM, nMix);
    for (int mHyp = 0; mHyp < nMix; mHyp++) {
      current_saem_state->_saemMixest = mHyp + 1;
      mat fcMat = nmRngGuard([&]{ return user_fn(phiC, mx.evtM, mx.optM); });
      vec curLimit = fcMat.col(2);
      vec curCens = fcMat.col(1);
      vec fc = fcMat.col(0);
      vec fs = fc;
      mat DYFm = zeros<mat>(mlen, mx.nM);
      switch (distribution) {
      case 1:
        {
          vec yt = hasFixedObsTransform ? yTrans : mx.y;
          if (!hasFixedObsTransform) {
            for (int i = ntotal; i--;) {
              int cur = ix_endpnt(i);
              yt(i) = _powerD(mx.y(i), lambda(cur), yj(cur), low(cur), hi(cur));
            }
          }
          const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
          for (int k = 0; k < nmc; k++) {
            int obs_start = k * ntotal;
            vec fsk = fs.subvec(obs_start, obs_start + ntotal - 1);
            const vec limitk = curLimit.subvec(obs_start, obs_start + ntotal - 1);
            const vec censk = curCens.subvec(obs_start, obs_start + ntotal - 1);
            _scratch_ft = fsk;
            _scratch_limitT = limitk;
            for (int i = ntotal; i--;) {
              int cur = ix_endpnt(i);
              _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
              _scratch_ft(i) = _powerD(fsk(i), lambda(cur), yj(cur), low(cur), hi(cur));
              _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fsk(i), false, true);
            }
            saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
            _scratch_g.elem(find(_scratch_g == 0.0)).fill(1);
            _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
            _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
            _scratch_indio = mx.indio + (arma::uword)k * stride;
            DYFm(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
            applyCensLoss(DYFm, _scratch_indio, censk, yt, _scratch_limitT,
                          _scratch_ftAr, _scratch_gAr);
          }
        }
        break;
      case 2:
        {
          const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
          for (int k = 0; k < nmc; k++) {
            vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
            _scratch_indio = mx.indio + (arma::uword)k * stride;
            DYFm(_scratch_indio) = -mx.y % log(fck) + fck;
          }
        }
        break;
      case 3:
        {
          const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
          for (int k = 0; k < nmc; k++) {
            vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
            _scratch_indio = mx.indio + (arma::uword)k * stride;
            DYFm(_scratch_indio) = -mx.y % log(fck) - (1 - mx.y) % log(1 - fck);
          }
        }
        break;
      case 4:
        {
          // general log-likelihood: model prediction is the per-obs loglik
          const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
          for (int k = 0; k < nmc; k++) {
            vec fck = fc.subvec(k * ntotal, (k + 1) * ntotal - 1);
            _scratch_indio = mx.indio + (arma::uword)k * stride;
            DYFm(_scratch_indio) = -fck;
          }
        }
        break;
      }
      lossByM.col(mHyp) = sum(DYFm, 0).t();
    }
    current_saem_state->_saemMixest = 0;

    vec result(mx.nM);
    for (int row = 0; row < mx.nM; row++) {
      double minL = lossByM(row, 0);
      for (int m = 1; m < nMix; m++) {
        if (lossByM(row, m) < minL) minL = lossByM(row, m);
      }
      double sumExp = 0.0;
      for (int m = 0; m < nMix; m++) {
        sumExp += mixProb(m) * exp(minL - lossByM(row, m));
      }
      if (sumExp <= 0.0) sumExp = 1e-300;
      result(row) = minL - log(sumExp);
    }
    return result;
  }

  // Model-aware naive classification for MSAEM's stratified init (run once before iteration 1):
  // for each subject/hypothesis, shift phi1's owned columns toward/away from that hypothesis
  // (pertSd BSV-SD) and evaluate the compiled model's actual fit. Returns the per-subject
  // argmin-hypothesis classification (1-indexed).
  uvec mixNaiveClassify(double pertSd) {
    double double_xmin = 1.0e-200;
    double xmax = 1e300;
    mat lossByHyp(N, nMix, fill::zeros);
    for (int mHyp = 0; mHyp < nMix; mHyp++) {
      mat cand = phiM;
      if (omegaShareSubpop.n_elem == (unsigned int)nphi1) {
        for (unsigned int c = 0; c < (unsigned int)nphi1; c++) {
          unsigned int subpop = omegaShareSubpop(c);
          if (subpop < 1) continue;
          double sdCol = std::sqrt(Gamma2_phi1(c, c));
          if (!std::isfinite(sdCol) || sdCol <= 0) continue;
          double shift = (subpop == (unsigned int)(mHyp + 1)) ? pertSd * sdCol : -pertSd * sdCol;
          cand.col(i1(c)) += shift;
        }
      }
      current_saem_state->_saemMixest = mHyp + 1;
      mat hypMat = user_fn(cand, evt, optM);
      vec fHyp = hypMat.col(0);
      vec censHyp = hypMat.col(1);
      vec limitHyp = hypMat.col(2);
      mat DYFhyp = zeros<mat>(mlen, nM);
      if (distribution == 1) {
        vec yt = hasFixedObsTransform ? yTrans : y;
        if (!hasFixedObsTransform) {
          for (int i = ntotal; i--;) {
            int cur = ix_endpnt(i);
            yt(i) = _powerD(y(i), lambda(cur), yj(cur), low(cur), hi(cur));
          }
        }
        const arma::uword stride = (arma::uword)N * (arma::uword)mlen;
        for (int k = 0; k < nmc; k++) {
          int obs_start = k * ntotal;
          vec fk = fHyp.subvec(obs_start, obs_start + ntotal - 1);
          const vec censk = censHyp.subvec(obs_start, obs_start + ntotal - 1);
          const vec limitk = limitHyp.subvec(obs_start, obs_start + ntotal - 1);
          _scratch_ft = fk;
          _scratch_limitT = limitk;
          for (int i = ntotal; i--;) {
            int cur = ix_endpnt(i);
            _scratch_limitT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
            _scratch_ft(i) = _powerD(fk(i), lambda(cur), yj(cur), low(cur), hi(cur));
            _scratch_ftT(i) = handleF(propT(cur), _scratch_ft(i), fk(i), false, true);
          }
          saemFormG(_scratch_g, vecares, vecbres, _scratch_ftT, veccres, vecaddProp);
          _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
          _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
          _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
          _scratch_indio = indio + (arma::uword)k * stride;
          DYFhyp(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
          applyCensLoss(DYFhyp, _scratch_indio, censk, yt, _scratch_limitT,
                        _scratch_ftAr, _scratch_gAr);
        }
      } else if (distribution == 2) {
        for (int k = 0; k < nmc; k++) {
          vec fk = fHyp.subvec(k * ntotal, (k + 1) * ntotal - 1);
          uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
          DYFhyp(indio_k) = -y % log(fk) + fk;
        }
      } else if (distribution == 3) {
        for (int k = 0; k < nmc; k++) {
          vec fk = fHyp.subvec(k * ntotal, (k + 1) * ntotal - 1);
          uvec indio_k = indio + (arma::uword)k * (arma::uword)(N * mlen);
          DYFhyp(indio_k) = -y % log(fk) - (1 - y) % log(1 - fk);
        }
      }
      vec U_y_hyp = sum(DYFhyp, 0).t();
      for (int i = 0; i < N; i++) {
        double sumL = 0.0;
        for (int k = 0; k < nmc; k++) sumL += U_y_hyp(i + k * N);
        lossByHyp(i, mHyp) = sumL / nmc;
      }
    }
    current_saem_state->_saemMixest = 0;

    uvec cls(N);
    for (int i = 0; i < N; i++) {
      unsigned int best = 0;
      double bestL = lossByHyp(i, 0);
      for (int m = 1; m < nMix; m++) {
        if (lossByHyp(i, m) < bestL) { bestL = lossByHyp(i, m); best = (unsigned int)m; }
      }
      cls(i) = best + 1;
    }
    return cls;
  }

  // MSAEM's S-step MCMC kernel: same proposal/prior machinery as do_mcmc(), but acceptance
  // uses mixture-weighted mixObsLoss() instead of a single hypothesis's DYF. No column is ever
  // masked or frozen -- every proposal is judged against the full mixture density, avoiding the
  // label-simulation instability the paper documents.
  void do_mcmc_msaem(const int method,
                      const int nu,
                      const mcmcaux &mx,
                      const mcmcphi &mphi,
                      mat &phiM,
                      vec &U_y,
                      vec &U_phi,
                      int kiter) {
    mat phiMc;
    vec Uc_y, Uc_phi, deltu;
    uvec ind;
    uvec i = mphi.i;
    arma::vec accU(mx.nM);

    for (int u = 0; u < nu; u++)
      for (int k1 = 0; k1 < mphi.nphi; k1++) {
        phiMc = phiM;
        // iteration-indexed threefry stream (msaem tag = method + 16 so it does
        // not collide with the plain do_mcmc streams)
        _saemSeedDoMcmc((uint32_t)saemSeed, kiter, method + 16, u, k1);
        switch (method) {
        case 1: {
          mat noise(mx.nM, mphi.nphi); _saemFillNormEng(noise);
          phiMc.cols(i) = noise * mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          break;
        }
        case 2: {
          mat noise(mx.nM, mphi.nphi); _saemFillNormEng(noise);
          phiMc.cols(i) = phiM.cols(i) +
            noise * mphi.Gdiag_phi % current_saem_state->_saemUE.cols(i);
          break;
        }
        case 3: {
          vec noise(mx.nM); _saemFillNormEng(noise);
          phiMc.col(i(k1)) = phiM.col(i(k1)) +
            noise * mphi.Gdiag_phi(k1, k1) % current_saem_state->_saemUE.col(i(k1));
          break;
        }
        }
        _saemFillUnifEng(accU);

        Uc_y = mixObsLoss(phiMc, mx);

        if (method == 1) {
          deltu = Uc_y - U_y;
        } else {
          mat dphic = phiMc.cols(i) - mphi.mprior_phiM;
          Uc_phi = 0.5 * sum(dphic % (dphic * mphi.IGamma2_phi), 1);
          deltu = Uc_y - U_y + Uc_phi - U_phi;
        }

        ind = find(deltu < -log(accU));
        phiM(ind, i) = phiMc(ind, i);
        U_y(ind) = Uc_y(ind);
        if (method > 1) {
          U_phi(ind) = Uc_phi(ind);
        }
        if (method < 3) {
          break;
        }
      }
    setRxThreadId(-1);
  }
};

// phi0 objective trampoline for the general-likelihood direct optimization.
static double gPhi0ObjR(Rcpp::NumericVector p) {
  for (size_t i = 0; i < gPhi0FreeIx.size(); ++i) {
    gPhi0Full[gPhi0FreeIx[i]] = p[i];
  }
  return gPhi0Self->phi0Objective(gPhi0Full.memptr());
}

static double gPhi1ObjR(Rcpp::NumericVector p) {
  for (size_t i = 0; i < gPhi1FreeIx.size(); ++i) {
    gPhi1Full[gPhi1FreeIx[i]] = p[i];
  }
  return gPhi1Self->phi1Objective(gPhi1Full.memptr());
}

//[[Rcpp::export]]
long saemPhi1RefineN_() { return _saemPhi1RefineN; }

static double gPhi0Obj1DR(double x) {
  gPhi0Work[gPhi0Coord] = x;
  return gPhi0Self->phi0Objective(gPhi0Work.memptr());
}

static double gPhi0RefObj(const double *p) {
  if (gPhi0RefEvalN >= gPhi0RefEvalMax) {
    // Out of budget: report a value worse than anything seen so the optimizer
    // collapses onto the best point instead of wandering.
    return gPhi0RefBestF + 1.0e10;
  }
  for (size_t i = 0; i < gPhi0FreeIx.size(); ++i) {
    int c = gPhi0FreeIx[i];
    double v = p[i];
    if (v < gPhi0Lo(c)) v = gPhi0Lo(c);
    if (v > gPhi0Hi(c)) v = gPhi0Hi(c);
    gPhi0Work[c] = v;
  }
  double f = gPhi0Self->phi0Objective(gPhi0Work.memptr());
  gPhi0RefEvalN++;
  if (gPhi0RefEvalN == 1 || f < gPhi0RefBestF) {
    gPhi0RefBestF = f;
    for (size_t i = 0; i < gPhi0FreeIx.size(); ++i) {
      gPhi0RefBest(i) = gPhi0Work[gPhi0FreeIx[i]];
    }
  }
  return f;
}

static void gPhi0NmFn(double *p, double *fx) { *fx = gPhi0RefObj(p); }

static double gPhi0RefObjR(Rcpp::NumericVector p) {
  return gPhi0RefObj(&(p[0]));
}


// closing for #ifndef __SAEM_CLASS_RCPP_HPP__
#endif

using namespace Rcpp;



t_calc_lhs saem_lhs = NULL;
t_update_inis saem_inis = NULL;

// Phase 4 (SAEM general-likelihood theta plan): true when setupRx registered
// predNoLhs/innerHess2 (odeSlotPred/odeSlotHess2) instead of the historic
// direct rxUpdateFuns(..., &rxInner) bind -- i.e. this fit is a general-
// likelihood fit whose saemPhi1TargetMap resolved (see R/saemPhi1Inner.R).
// SAEM's own per-iteration likelihood read is then reparameterized through
// predNoLhs/innerHess2 too (user_function's own phi1RowSolve), rather than
// pooling a separately-parameterized "SAEM's own model" peer -- see
// setupRx's own comment for why. false for EVERY normal/poisson/binomial
// fit and for a general-lik fit whose map did not resolve (covariate on a
// phi0/phi1 theta, unpaired eta, ...) -- both keep the original
// single-model rxInner path, byte-identically.
static bool _saemPhi1PoolActive = false;
bool _saemPhi1PoolReady = false;
bool _saemPhi1UseAnalyticHess = false;
arma::ivec _saemPhi1H2ThetaKind;
arma::ivec _saemPhi1H2ThetaCol;
arma::vec _saemPhi1H2ThetaFixedVal;
arma::ivec _saemPhi1H2EtaCol;
// saemControl(phi1Hessian=FALSE) is the default -- see its own docs (an
// ablation check found the Laplace log|H| correction was not what fixed a
// diverging Gaussian twin, and it can dominate/diverge for a heavy-tailed
// t()/cauchy() endpoint, nlmixr2/nlmixr2est#999).  When false, phi1Objective
// scores plain -2*loglik with no Hessian term at all (not even Omega^-1) --
// solving only odeSlotPred, never odeSlotHess2/the FD-Hessian fallback.
bool _saemPhi1WantHessian = false;
int _saemPhi1PredOffset = -1;
int _saemPhi1H2PredOffset = -1;
int _saemPhi1H2HessOffset = -1;
int _saemPhi1DvCol = -1;
int _saemPhi1DvColHess2 = -1;
arma::uvec _saemPhi1I0;
arma::uvec _saemPhi1I1;

rx_solve* _rx = NULL;

RObject mat2NumMat(const mat &m) {
  RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
  x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
  return x;
}

CharacterVector parNames;

// Phase 4 (SAEM general-likelihood theta plan): set every row's THETA[k]/
// ETA[k] for the odeSlotPred solve, from _phi (already mu+eta summed) --
// THETA[k]=this row's phi value, ETA[k]=0, so THETA[k]+ETA[k] reproduces the
// combined phi exactly (see setupRx's own comment).  Called ONCE, mirroring
// the original "fill in subject parameter information" loop's own timing --
// parameters do not change across the bad-solve retry loop, only tolerance
// does, so this must not be re-run per retry.
static void saemSetRowsPooled(const mat &_phi) {
  int nInd = (int)_phi.n_rows;
  int nH2Theta = (int)_saemPhi1H2ThetaKind.n_elem;
  int nEta = (int)_saemPhi1H2EtaCol.n_elem;
  for (int i = 0; i < nInd; ++i) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
    for (int k = 0; k < nH2Theta; ++k) {
      int kind = _saemPhi1H2ThetaKind(k);
      int col = _saemPhi1H2ThetaCol(k);
      double val = (kind == 1) ? _phi(i, _saemPhi1I1(col)) :
        ((kind == 0) ? _phi(i, _saemPhi1I0(col)) : _saemPhi1H2ThetaFixedVal(k));
      setIndParPtr(ind, k, val);
    }
    for (int k = 0; k < nEta; ++k) setIndParPtr(ind, nH2Theta + k, 0.0);
  }
}

// Phase 4: solve odeSlotPred (cheap, no sensitivities -- a Hessian is never
// needed for a plain likelihood read) for every row, under an
// OdeSwapScope/OdeSwapCmtScope (the pool is sized for the larger
// odeSlotHess2 peer when it exists).  par_solve is itself just this same
// per-individual ind_solve loop with rxode2's own bookkeeping around it, so
// this changes nothing structurally except the guards and the model solved.
// Parallelized the same way inner.cpp's own per-subject loops are; called
// once per (re)solve, same as the original par_solve(_rx) call it replaces
// (including from inside the bad-solve retry loop).
static void saemSolveIndividualsPooled(int nInd) {
  rx_solving_options *op = getSolvingOptions(_rx);
  int cores = getOpCores(op);
  bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  // odeSlotPred carries no event sensitivities of its own, but the ES shape is
  // still a process global: without a batch here, a shape some earlier solve
  // left installed (this fit's own odeSlotHess2 call, or a prior FOCEi fit's
  // fit-wide inner load) stays live and handle_evid injects jumps sized for
  // that OTHER model into this (smaller, compacted) solve's dosing events.
  // OdeSwapEsBatch's own constructor is what deactivates it for a "no ES" slot
  // -- see its contract in src/odeSwap.cpp.  Must be constructed outside the
  // OpenMP region below, matching every other peer solve's own batch.
  OdeSwapEsBatch predEsBatch(odeSlotPred);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int i = 0; i < nInd; ++i) {
#ifdef _OPENMP
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
    OdeSwapScope neqGuard(odeSlotPred, ind, op);
    OdeSwapCmtScope cmtGuard(odeSlotPred, op, ind);
    setIndSolve(ind, -1);
    odeSwapSolveInd(odeSlotPred, i);
  }
}

// Phase 4: read rx_pred_ back from the odeSlotPred solve above, ONE g row
// per evid==0 observation (rx_pred_ is a PER-OBSERVATION value, summed by
// the CALLER -- e.g. phi0Objective's accu() -- not pre-summed here),
// matching the original per-timepoint saem_lhs loop exactly.  Called once,
// after the bad-solve retry loop finishes, mirroring where the original
// read loop sits.  Parallelized the same way inner.cpp's own per-subject
// loops are; the sequential second pass only copies already-computed
// per-observation values into g, so elt's running count stays deterministic
// regardless of thread scheduling.
static void saemReadRowsPooled(mat &g, int &elt, bool &hasNan, int nInd) {
  rx_solving_options *op = getSolvingOptions(_rx);
  int cores = getOpCores(op);
  bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  bool hasCens = hasRxCens(_rx), hasLimit = hasRxLimit(_rx);
  std::vector<std::vector<double> > rowObs((size_t)nInd);
  std::vector<int> rowNan((size_t)nInd, 0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int i = 0; i < nInd; ++i) {
#ifdef _OPENMP
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
    OdeSwapScope neqGuard(odeSlotPred, ind, op);
    OdeSwapCmtScope cmtGuard(odeSlotPred, op, ind);
    int nAll = getIndNallTimes(ind);
    std::vector<double> &obs = rowObs[(size_t)i];
    obs.reserve((size_t)nAll);
    if (odeSwapIndBadSolveSlot(op, ind, odeSlotPred)) {
      for (int j = 0; j < nAll; ++j) {
        if (getIndEvid(ind, getIndIx(ind, j)) == 0) obs.push_back(1.0e99);
      }
      rowNan[i] = 1;
      continue;
    }
    iniSubjectE(i, 1, ind, op, _rx, rxPred.update_inis);
    double *lhs = neqGuard.lhs();
    for (int j = 0; j < nAll; ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, j);
      if (getIndEvid(ind, kk) != 0) continue;
      double curT = getTime(kk, ind);
      rxPred.calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
      double cur = lhs[_saemPhi1PredOffset];
      if (std::isnan(cur)) { cur = 1.0e99; rowNan[i] = 1; }
      obs.push_back(cur);
    }
  }
  for (int i = 0; i < nInd; ++i) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
    const std::vector<double> &obs = rowObs[(size_t)i];
    size_t oi = 0;
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      int kk = getIndIx(ind, j);
      if (getIndEvid(ind, kk) != 0) continue;
      g(elt, 0) = obs[oi++];
      g(elt, 1) = hasCens ? getIndCens(ind, kk) : 0.0;
      g(elt, 2) = hasLimit ? getIndLimit(ind, kk) : R_NegInf;
      elt++;
    }
    if (rowNan[i]) hasNan = true;
  }
}

mat user_function(const mat &_phi, const mat &_evt, const List &_opt) {
  // yp has all the observations in the dataset
  rx_solving_options_ind *ind;
  rx_solving_options *op = getSolvingOptions(_rx);
  // _phi has N*nmc rows (all chains); _evt has only N subjects (chain 0 template)
  int _Nnlmixr2 = (int)_phi.n_rows;
  SEXP paramUpdate = _opt["paramUpdate"];
  int *doParam = INTEGER(paramUpdate);
  int nPar = Rf_length(paramUpdate);

  int *indMixest = nullptr;
  if (_opt.hasAttribute("names")) {
    if (_opt.containsElementNamed("mixest")) {
      SEXP indMixestR = _opt["mixest"];
      if (indMixestR != R_NilValue) {
        indMixest = INTEGER(indMixestR);
      }
    }
  }


  // Fill in subject parameter information.  Phase 4 (SAEM general-lik theta
  // plan): a fit whose phi1 map resolved sets THETA[k]/ETA[k] instead (see
  // saemSetRowsPooled), reparameterizing the SAME combined phi value rather
  // than using SAEM's own native paramUpdate/doParam positions -- every
  // other fit takes the unchanged path below.
  if (_saemPhi1PoolActive) saemSetRowsPooled(_phi);
  for (int _i = 0; _i < _Nnlmixr2; ++_i) {
    ind = getSolvingOptionsInd(_rx, _i);
    setIndSolve(ind, -1);
    if (current_saem_state != nullptr && current_saem_state->_saemMixest != 0) {
      setIndMixest(ind, current_saem_state->_saemMixest);
    } else if (indMixest != nullptr) {
      setIndMixest(ind, indMixest[_i]);
    }
    if (!_saemPhi1PoolActive) {
      int k=0;
      for (int _j = 0; _j < nPar; _j++){
        if (doParam[_j] == 1) {
          setIndParPtr(ind, _j, _phi(_i, k++));
        }
      }
    }
  }
  // ODE-freeze: phi0 (general-likelihood fixed effects) do not change the ODE, so
  // during the phi0 optimization the states are reused and only the log-likelihood
  // is recomputed below -- skip the (costly) re-integration entirely.
  if (!_saemFreezeOde) {
  resetRxBadSolve(_rx);
  // Phase 4: a general-lik fit whose phi1 map resolved routes SAEM's own
  // solve through odeSlotPred (see setupRx) instead of the bulk par_solve --
  // every other fit is byte-identical to before.
  if (_saemPhi1PoolActive) {
    saemSolveIndividualsPooled(_Nnlmixr2);
  } else {
    par_solve(_rx); // Solve the complete system (possibly in parallel)
  }
  int j=0;
  while (hasRxBadSolve(_rx) && j < current_saem_state->_saemMaxOdeRecalc){
    current_saem_state->_saemIncreaseTol=1;
    if (current_saem_state->_saemIndTolRelax) {
      // Only loosen tolerance for subjects whose ODE solve produced NaN/Inf.
      // Tolerance is sticky via ind->tolFactor so iniSubject reapplies it on
      // subsequent SAEM iterations; genuinely stiff subjects stay loosened.
      if (getOpNeq(op) > 0) {
        for (int _i = 0; _i < _Nnlmixr2; _i++) {
          rx_solving_options_ind *_indI = getSolvingOptionsInd(_rx, _i);
          double *_solveI = getIndSolve(_indI);
          int _nsolveI = getOpNeq(op) * getIndNallTimes(_indI);
          for (int _ns = 0; _ns < _nsolveI; _ns++) {
            if (ISNA(_solveI[_ns]) || std::isnan(_solveI[_ns]) || std::isinf(_solveI[_ns])) {
              setIndTolFactor(_indI, getIndTolFactor(_indI) * current_saem_state->_saemOdeRecalcFactor);
              break;
            }
          }
        }
      }
    } else {
      // Loosen all subjects uniformly; reset after retry.
      for (int _i = 0; _i < _Nnlmixr2; _i++) {
        rx_solving_options_ind *_indI = getSolvingOptionsInd(_rx, _i);
        setIndTolFactor(_indI, getIndTolFactor(_indI) * current_saem_state->_saemOdeRecalcFactor);
      }
    }
    resetRxBadSolve(_rx);
    if (_saemPhi1PoolActive) {
      saemSolveIndividualsPooled(_Nnlmixr2);
    } else {
      par_solve(_rx);
    }
    j++;
  }
  if (!current_saem_state->_saemIndTolRelax && j != 0) {
    // Reset all subjects' tolFactor after the non-selective retry.
    for (int _i = 0; _i < _Nnlmixr2; _i++) {
      setIndTolFactor(getSolvingOptionsInd(_rx, _i), 1.0);
    }
  }
  } // end (!_saemFreezeOde) integration guard
  // indTolRelax=TRUE: stiff subjects retain their loosened tolFactor across iterations.
  mat g(getRxNsim(_rx) * getRxNobs2(_rx), 3); // nobs across all chains
  int elt=0;
  bool hasNan = false;
  if (_saemPhi1PoolActive) {
    saemReadRowsPooled(g, elt, hasNan, _Nnlmixr2);
  } else {
  for (int id = 0; id < _Nnlmixr2; ++id) {
    ind = getSolvingOptionsInd(_rx, id);
    double *lhs = getIndLhs(ind);
    iniSubjectE(getOpNeq(op), 1, ind, op, _rx, saem_inis);
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, getIndIdx(ind));
      double curT = getTime(kk, ind);
      if (isDose(getIndEvid(ind, kk))) {
        // Need to calculate for advan sensitivities
        saem_lhs((int)id, curT,
                 getOpIndSolve(op, ind, j), lhs);
      } else if (getIndEvid(ind,kk) == 0) {
        saem_lhs((int)id, curT,
                 getOpIndSolve(op, ind, j), lhs);
        double cur = lhs[0];
        if (std::isnan(cur)) {
          cur = 1.0e99;
          hasNan = true;
        }
        g(elt, 0) = cur;
        if (hasRxCens(_rx)) {
          g(elt, 1) = getIndCens(ind, kk);
        } else {
          g(elt, 1) = 0;
        }
        if (hasRxLimit(_rx)) {
          g(elt, 2) = getIndLimit(ind, kk);
        } else {
          g(elt, 2) = R_NegInf;
        }
        elt++;
      } // evid=2 does not need to be calculated
    }
  }
  }
  if (solveMethodThreadSafe(op)) { // liblsoda
    // Order by the overall solve time
    // Should it be done every time? Every x times?
    sortIds(_rx, 0);
  }
  if (hasNan && !_warnAtolRtol) {
    RSprintf("NaN in prediction; Consider: relax atol & rtol; change initials; change seed; change structural model\n  warning only issued once per problem\n");
    _warnAtolRtol = true;
  }
  return g;
}

// Set up the rxode2 solve structure for N subjects across nmc chains.
// Passing an N*nmc-row params matrix triggers rxode2's nsim mechanism:
//   nsim = nPopPar / nsub = (N*nmc) / N = nmc
// Chains 1..nmc-1 automatically share chain 0's event data pointers
// (all_times, evid, dose, ii, idose, cov_ptr) while each subject retains
// its own solve/ix/tolFactor buffers, reducing event-table memory by ~nmc times.
void setupRx(List &opt, SEXP evt, int nmc, int N) {
  RObject obj = opt[".rx"];
  List mv = _rxode2_rxModelVars_(obj);
  parNames = mv[RxMv_params];

  // Phase 4 (SAEM general-likelihood theta plan): a general-lik fit whose
  // saemPhi1TargetMap resolved (R/saemPhi1Inner.R) attaches saemPhi1Pred
  // (always) and saemPhi1Hess2 (NULL for some model shapes, e.g. linCmt()).
  // Rather than pooling SAEM's own (differently-parameterized, native-name)
  // model alongside these, solve EVERYTHING -- SAEM's own per-iteration
  // likelihood read AND the phi1 theta step's Laplace objective -- through
  // predNoLhs/innerHess2 directly: rx_pred_ there IS the same log-density
  // SAEM's own flattened model computes, just parameterized as THETA[k]
  // (mu) + ETA[k] (deviation) instead of one pre-summed phi value. Setting
  // THETA[k]=phi_value, ETA[k]=0 reproduces the combined phi exactly, so no
  // separate "SAEM's own model" peer/parameter-order reconciliation is
  // needed at all -- there is only ONE parameterization in play, matching
  // predNoLhs's/innerHess2's own (see phi1RowSolve, user_function). EVERY
  // other fit (normal/poisson/binomial, or a general-lik fit whose map did
  // not resolve) takes the unchanged rxUpdateFuns(..., &rxInner) path below.
  // The odeSwap registry is a global that outlives a fit (mirrors
  // foceiFitCpp_'s own "start every fit from empty" comment, src/inner.cpp)
  // -- otherwise a LATER non-pooled fit's own rxSolve_(obj,...) call still
  // sees a PRIOR pooled fit's declared Hess2/Pred peers in the registry,
  // corrupting its own pool sizing even though _saemPhi1PoolActive is false
  // for it (measured: a normal-model fit run after a general-lik pooled fit
  // in the same R session segfaulted from exactly this).
  odeSwapClearAll();
  _saemPhi1PoolActive = opt.containsElementNamed("saemPhi1Pred") &&
    !Rf_isNull(opt["saemPhi1Pred"]);
  if (_saemPhi1PoolActive) {
    bool haveHess2 = opt.containsElementNamed("saemPhi1Hess2") &&
      !Rf_isNull(opt["saemPhi1Hess2"]);
    if (haveHess2) odeSwapDeclare(odeSlotHess2, "hess2", opt["saemPhi1Hess2"]);
    odeSwapDeclare(odeSlotPred, "pred", opt["saemPhi1Pred"]);
    // rxSolve_ on whichever peer has the most states (innerHess2's extra
    // eta-sensitivity states when it built, else predNoLhs) -- matches the
    // number-of-ODEs sizing rule odeSwap already uses for neq; both share
    // the identical THETA[]/ETA[]/DV parameter declaration (verified by
    // .saemPhi1TargetMap), so either works as the params-matrix source.
    RObject widePar = haveHess2 ? opt["saemPhi1Hess2"] : opt["saemPhi1Pred"];
    List odeO = opt["rxControl"];
    List wideMv = _rxode2_rxModelVars_(widePar);
    CharacterVector wideParNames = wideMv[RxMv_params];
    int npars = wideParNames.size();
    int nrows = N * nmc;
    NumericMatrix parsM(nrows, npars);
    // Placeholder only -- every real THETA[k]/ETA[k] value is written fresh
    // by phi1RowSolve before every solve; DV is a per-observation covariate
    // the shared event table supplies (rxSolve_ still requires a params-
    // matrix entry for every declared parameter before it can match evt's
    // own covariate columns).
    for (int k = 0; k < nrows; k++)
      for (int j = 0; j < npars; j++) parsM(k, j) = 1.1;
    parsM.attr("dimnames") = List::create(R_NilValue, wideParNames);
    rxode2::rxSolve_(widePar, odeO,
                     R_NilValue, R_NilValue,
                     parsM, evt, R_NilValue, 1);
    // Register AFTER rxSolve_ has sized/built the pool -- registering (which
    // rxDynLoad's) before it exists rebinds rxode2's event-sensitivity
    // globals and corrupts the solve (see odeSwap.h).
    if (haveHess2) odeSwapRegister(odeSlotHess2, "hess2", opt["saemPhi1Hess2"], &rxHess2);
    odeSwapRegister(odeSlotPred, "pred", opt["saemPhi1Pred"], &rxPred);
    return;
  }

  rxUpdateFuns(mv["trans"], &rxInner);
  if (!Rf_isNull(obj)){
    RObject pars0 = opt[".pars"];
    List odeO = opt["rxControl"];
    if (Rf_isNull(pars0)) {
      stop("params must be non-nil");
    }
    NumericVector parsV = as<NumericVector>(pars0);
    int npars = parsV.size();
    int nrows = N * nmc;
    NumericMatrix parsM(nrows, npars);
    CharacterVector parsNames = parsV.names();
    for (int k = 0; k < nrows; k++) {
      for (int j = 0; j < npars; j++) parsM(k, j) = parsV[j];
    }
    parsM.attr("dimnames") = List::create(R_NilValue, parsNames);
    rxode2::rxSolve_(obj, odeO,
                     R_NilValue, R_NilValue,
                     parsM, evt, R_NilValue, 1);
  } else {
    stop("cannot find rxode2 model");
  }
}

// Unused when _saemPhi1PoolActive (user_function's own phi1RowSolve reads
// predNoLhs/innerHess2 directly instead) -- harmless to still point these at
// rxInner, which every fit's setupRx (this function) populates regardless.
static inline void saemSetLhsInis() {
  saem_lhs = rxInner.calc_lhs;
  saem_inis = rxInner.update_inis;
}

//[[Rcpp::export]]
SEXP saem_do_pred(SEXP in_phi, SEXP in_evt, SEXP in_opt) {
  List opt = List(in_opt);
  mat phi = as<mat>(in_phi);
  setupRx(opt, in_evt, 1, (int)phi.n_rows);
  saemSetLhsInis();
  _rx=getRxSolve_();
  mat evt = as<mat>(in_evt);
  saem_state_t dummy_st;
  if (opt.containsElementNamed("maxOdeRecalc")) dummy_st._saemMaxOdeRecalc = abs(as<int>(opt["maxOdeRecalc"]));
  if (opt.containsElementNamed("odeRecalcFactor")) dummy_st._saemOdeRecalcFactor = fabs(as<double>(opt["odeRecalcFactor"]));
  if (opt.containsElementNamed("indTolRelax")) dummy_st._saemIndTolRelax = as<bool>(opt["indTolRelax"]);
  if (opt.containsElementNamed("ue")) dummy_st._saemUE = as<mat>(opt["ue"]);
  current_saem_state = &dummy_st;
  mat gMat = user_function(phi, evt, opt);
  current_saem_state = nullptr;
  vec g = gMat.col(0);
  return wrap(g);
}


//[[Rcpp::export]]
SEXP saem_fit(SEXP xSEXP) {
  List x(xSEXP);
  List opt = x["opt"];
  setupRx(opt, x["evt"], as<int>(x["nmc"]), as<int>(x["N"]));

  // if (rxSingleSolve == NULL) rxSingleSolve = (rxSingleSolve_t) R_GetCCallable("rxode2","rxSingleSolve");
  saemSetLhsInis();
  _rx=getRxSolve_();

  saem_state_t dummy_st;
  if (opt.containsElementNamed("maxOdeRecalc")) dummy_st._saemMaxOdeRecalc = abs(as<int>(opt["maxOdeRecalc"]));
  if (opt.containsElementNamed("odeRecalcFactor")) dummy_st._saemOdeRecalcFactor = fabs(as<double>(opt["odeRecalcFactor"]));
  if (opt.containsElementNamed("indTolRelax")) dummy_st._saemIndTolRelax = as<bool>(opt["indTolRelax"]);
  if (opt.containsElementNamed("ue")) dummy_st._saemUE = as<mat>(opt["ue"]);
  current_saem_state = &dummy_st;

  SAEM saem;
  saem.inits(x);
  saem.set_fn(user_function);

  saem.saem_fit();

  int _saemNsub = (int)getRxNsub(_rx);
  NumericVector _saemTf(_saemNsub);
  for (int _i = 0; _i < _saemNsub; _i++) {
    _saemTf[_i] = getIndTolFactor(getSolvingOptionsInd(_rx, _i));
  }
  List out = List::create(
    Named("resMat") = saem.get_resMat(),
    Named("arCor") = saem.get_arCor(),
    Named("transMat") = saem.get_trans(),
    Named("mprior_phi") = saem.get_mprior_phi(),
    Named("mpost_phi") = saem.get_mpost_phi(),
    Named("Gamma2_phi1") = saem.get_Gamma2_phi1(),
    Named("Gamma2_phi1Report") = saem.get_Gamma2_phi1Report(),
    Named("Plambda") = saem.get_Plambda(),
    Named("Ha") = saem.get_Ha(),
    Named("sig2") = saem.get_sig2(),
    Named("eta") = saem.get_eta(),
    Named("par_hist") = saem.get_par_hist(),
    Named("ueRevisitInfo") = saem.get_ueRevisitInfo(),
    Named("HaSa") = saem.get_HaSa(),
    Named("res_info") = saem.get_resInfo(),
    Named("tolFactor") = _saemTf,
    Named("mixProb") = wrap(saem.get_mixProb()),
    Named("mixWeights") = wrap(saem.get_mixWeights())
  );
  current_saem_state = nullptr;
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
}

// Test-only wrapper: exposes the E-step's per-observation combined-error SD
// (saemFormG(), used at every _scratch_g site) so it can be pinned against
// the M-step's combined1/combined2 formulas without running a full fit.
//[[Rcpp::export]]
SEXP saemFormGTest(SEXP inA, SEXP inB, SEXP inFt, SEXP inC, SEXP inAddProp) {
  vec a = as<vec>(inA);
  vec b = as<vec>(inB);
  vec ft = as<vec>(inFt);
  vec c = as<vec>(inC);
  uvec addPropVec = as<uvec>(inAddProp);
  vec g(a.n_elem);
  saemFormG(g, a, b, ft, c, addPropVec);
  return wrap(g);
}
