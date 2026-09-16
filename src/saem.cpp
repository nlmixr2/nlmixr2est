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
#include "nmSeqSeed.h"
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
#include "nonMuThetaGrad.h"
#include "shi21.h"
#include <n1qn1c.h>

// The declared-distribution family dispatch and the ODE-free maximum-
// likelihood M-step live in their own translation unit so saem and imp share
// one implementation -- see src/etaDistFam.h.
#include "etaDistFam.h"
// eta-scale primitives: bounds, bijectors, the copula joint density.  Only the
// direct parameterization reaches these; the cdf route never does.
#include "etaDistEtaScale.h"

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
// The residual-step transform cache used to be keyed on the ADDRESS of ysb/fsb.
// Those are per-iteration locals, freed and reallocated every M-step, so the
// allocator routinely hands back the same address with completely different
// CONTENTS -- the guard then reported "unchanged" and the residual optimizer
// scored every later iteration against the FIRST iteration's predictions.  Key
// it on a counter bumped wherever the data behind the cache is set instead.
static unsigned long _saemResidGen = 0;
static unsigned long _saemCacheGen = (unsigned long)-1;
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
  if (_saemCacheGen == _saemResidGen &&
      _saemCacheYptr == _saemYptr &&
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
  _saemCacheGen = _saemResidGen;
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
  int block;  // 0 phi1, 1 phi0 (seed layout)
  uvec i;
  mat Gamma_phi;
  mat Gdiag_phi;
  // saemControl(rwOmega=): NONMEM's mode-2 random walk proposes from
  // Z = lambda*Omega (technical guide eq. 1.139) -- the FULL covariance, so the
  // walk moves along the posterior's correlated directions.  Gdiag_phi is the
  // diagonal saemix and nlmixr2 have always used; this is chol(Omega)*rmcmc,
  // the Omega-shaped alternative.
  mat Gfull_phi;
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
// The iteration whose COMPLETE-SYSTEM solve is currently established, or -1.
//
// With the sensitivity peer live, one solve carries rx_pred_ AND every
// d(f)/d(theta) together (_saemSolveCompleteOnce), so it serves the gradient
// step and the single-solve pred both -- there is no reason to solve the
// system more than once per step unless something MOVES phi0 underneath it,
// which is why the writers invalidate rather than the readers re-solving.
static int _saemCompleteSolveIter = -1;
// ...and whether that established solve carries the SENSITIVITIES.  The
// sensitivity solve is the expensive one and only the gradient step needs it;
// the derivative-free search wants the plain prediction solve, which is
// faster.  A sensitivity solve satisfies a later plain request (it is a
// superset); the reverse re-solves.
static bool _saemCompleteSolveHasSens = false;
static double gPhi0Obj1DR(double x);
// ---- Shi-difference fallback for the non-mu theta gradient ----------------
//
// When the sensitivity peer's bad-solve ladder is exhausted the analytic
// d(f)/d(theta) is unavailable.  focei's answer in that situation is to fall
// back to the ORIGINAL model -- no sensitivities -- and difference it with the
// Shi (2021) step search; this is the same fallback for saem's non-mu step.
//
// gShiPredFn is the vector-valued function shi21Forward differentiates: it
// writes the candidate free phi0 values into the phi0 columns and returns the
// population prediction column.  It solves through user_fn with
// _saemSolveCompleteOnce left at 0, which routes to _saemOwnSolveSlot
// (odeSlotPred) rather than the sensitivity peer -- i.e. exactly "the original
// model without the theta gradients", and the same solve phi0Objective uses,
// which is known to succeed where the sensitivity read did not.
//
// shi21fn_type is a bare function pointer, so the instance is reached through a
// file-static the way gPhi0ObjR already does.
static SAEM *gShiSelf = nullptr;
static std::vector<int> gShiFreeIx;
// ---- general declared-distribution M-step objective ------------------------
// nelder_fn takes a bare function pointer, so the problem is reached through
// file statics the way gPhi0ObjR already is.  Only ever touched from the
// serial part of the iteration.
static int gEdFam = -1;
static const std::vector< std::vector<etaDistTok> > *gEdRpn = NULL;
static int gEdNth = 0, gEdNSym = 0, gEdNRec = 0;
static const double *gEdRec = NULL, *gEdEta = NULL, *gEdWt = NULL;
// newuoa's own maxfun stop returns whatever point it was holding rather than
// the best one it saw, so the objective tracks the best itself -- the same
// convention gPhi0RefObjR uses in this file.
static double gEdBest = R_PosInf;
static std::vector<double> gEdBestPar;
// PAIR mode: when gEdRpn2 is set the objective is the JOINT copula density of
// two declarations rather than one marginal.  Under a copula the marginals'
// joint MLE is not the two separate marginal MLEs, so a correlated pair has to
// be one optimization problem, not two.  gEdRhoIdx >= 0 additionally puts
// atanh(rho) in the parameter vector, which is what makes the correlation
// estimated by the copula's own density instead of by a moment statistic.
static int gEdFam2 = -1;
static const std::vector< std::vector<etaDistTok> > *gEdRpn2 = NULL;
static const double *gEdEta2 = NULL;
static double gEdRho = 0.0;
static int gEdRhoIdx = -1;
static inline bool gEdIsPair() { return gEdRpn2 != NULL; }
static double gEdObj(double *p) {
  double v = 0.0;
  if (gEdRpn == NULL) return 1e300;
  bool ok = gEdIsPair() ?
    rxEtaDistPairLoglikObj(gEdFam, *gEdRpn, gEdFam2, *gEdRpn2, gEdNth, gEdNSym,
                           p, gEdRec, gEdEta, gEdEta2, gEdRho, gEdRhoIdx,
                           gEdWt, gEdNRec, &v) :
    rxEtaDistLoglikObj(gEdFam, *gEdRpn, gEdNth, gEdNSym, p,
                       gEdRec, gEdEta, gEdWt, gEdNRec, &v);
  if (!ok) return 1e300;
  return -v;   // the optimizers here minimize
}

// n1qn1 cost: value and gradient together, and NO R API anywhere in it -- the
// whole reason this objective is in C++ is to be callable from the OpenMP
// regions the focei family and imp use, and n1qn1_ is a plain C function
// pointer.  (newuoa is not thread-safe and would have to be reached through
// Rcpp::Function, which is doubly disqualifying here.)
static int gEdN1Bad = 0;
static int gEdN1Evals = 0;
static void gEdN1Cost(int *ind, int *n, double *x, double *f, double *g,
                      int *ti, float *tr, double *td, int *id) {
  (void)ti; (void)tr; (void)td; (void)id; (void)n;
  double v = 0.0;
  std::vector<double> gr((size_t)gEdNth, 0.0);
  bool okg = gEdRpn != NULL &&
    (gEdIsPair() ?
     rxEtaDistPairLoglikGrad(gEdFam, *gEdRpn, gEdFam2, *gEdRpn2, gEdNth, gEdNSym,
                             x, gEdRec, gEdEta, gEdEta2, gEdRho, gEdRhoIdx,
                             gEdWt, gEdNRec, &v, gr.data()) :
     rxEtaDistLoglikGrad(gEdFam, *gEdRpn, gEdNth, gEdNSym, x,
                         gEdRec, gEdEta, gEdWt, gEdNRec, &v, gr.data()));
  if (!okg) {
    gEdN1Bad = 1;
    if (*ind == 2 || *ind == 4) *f = 1e300;
    if (*ind == 3 || *ind == 4) for (int t = 0; t < gEdNth; ++t) g[t] = 0.0;
    return;
  }
  gEdN1Evals++;
  // maximizing the log-likelihood, so minimize its negative -- gradient too
  if (*ind == 2 || *ind == 4) *f = -v;
  if (*ind == 3 || *ind == 4) for (int t = 0; t < gEdNth; ++t) g[t] = -gr[(size_t)t];
  if (-v < gEdBest) {
    gEdBest = -v;
    gEdBestPar.assign(x, x + gEdNth);
  }
}

static arma::vec gShiPredFn(arma::vec &t, int id);

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
// n1qn1 driver for the non-mu (phi0) refinement.  Available only with the exact
// gradient: a quasi-Newton is pointless without one, and with one each iteration
// costs a single COMPLETE solve (objective and gradient together) where the
// derivative-free search spends one solve per objective evaluation and never
// sees a derivative.
static void gPhi0N1Cost(int *ind, int *n, double *x, double *f, double *g,
                        int *ti, float *tr, double *td, int *id);
static int gPhi0N1Bad = 0;
static int gPhi0N1Evals = 0;
static void gPhi0NmFn(double *p, double *fx);
static double gPhi0RefObjR(Rcpp::NumericVector p);

// saemControl(zeroOmegaDirect=): direct maximization of the observation
// likelihood in the phi1 columns whose declared omega was zero.  Same
// file-scope-callback shape as the gPhi0 block above, because nelder_fn()
// takes a plain function pointer.
static SAEM* gZeroOmSelf = nullptr;
static std::vector<int> gZeroOmIx;   // phi1 columns being optimized
static arma::vec gZeroOmLo, gZeroOmHi, gZeroOmBest;
static double gZeroOmBestF = 0.0;
static int gZeroOmEvalN = 0, gZeroOmEvalMax = 0;
static double gZeroOmObj(const double *p);
static void gZeroOmNmFn(double *p, double *fx);

// phi1 objective for the general-likelihood Laplace-corrected direct
// optimization (Phase 4, SAEM general-likelihood theta plan) -- the phi1
// sibling of gPhi0Self/gPhi0ObjR above, same bounded-bobyqa/.boundedResidOpt
// wiring.  Defined after the class body.
static SAEM* gPhi1Self = nullptr;
// Reached from the solve passes, which are static free functions, to
// harvest the declared-distribution argument anchors out of each
// individual's lhs.  Same pattern as gPhi0Self/gPhi1Self above.
static SAEM* gAnchorSelf = nullptr;
// Set per solve request: is this the state solve, or a candidate evaluation
// from the phi0 search?  The pooled read path is a separate function and cannot
// see user_function's local, so the decision is shared here.
static bool gAnchorAtState = false;
static double gPhi1ObjR(Rcpp::NumericVector p);
static arma::vec gPhi1Full;
static std::vector<int> gPhi1FreeIx;
// Diagnostic: how many times refinePhi1Lik actually ran, so a test can prove
// the mechanism executed rather than infer it from matching estimates alone
// (evaluation criterion #2 in the plan).
static long _saemPhi1RefineN = 0;
// Same idea for the declared-distribution M-step: every way it can decline to
// run is silent, and a fit whose estimates look reasonable is no evidence that
// it engaged.  Counted so a test can prove the mechanism executed rather than
// infer it from matching estimates -- which it does not: with the M-step inert,
// etaDistMstep=TRUE and FALSE give bit-identical fits.
static long _saemEtaDistN = 0;
// Whether THIS fit asked for the declared-distribution M-step at all, so
// "never ran" can be told apart from "was never requested".
static int _saemEtaDistOn = 0;
// etaDistLoglik: the declared thetas are estimated from the OBSERVATION
// likelihood instead, so the family M-step standing down is the intended
// behaviour and must not be reported as a no-op.
static int _saemEtaDistObsLik = 0;
// declared correlations the sufficient statistic says the data do not identify
static int _saemEtaDistCorNotEst = 0;

// Closed-form sequential seed layout of one SAEM fit (nmSeqSeed.h).  Per
// iteration: every mixture component's MCMC steps -- phi1 then phi0, methods 1,
// 2 and 3 (3 once per phi column), nu times each, then mode 1B nu1B times -- one
// seed per chain row, then every component's censored-value draws, one seed per
// observation.  With nu1B = 0 this is the layout of the branch without mode 1B.
struct saemSeedLayout {
  uint64_t nu[3] = {0, 0, 0}, nu1B = 0, nphi[2] = {0, 0};
  uint64_t nComp = 1, nM = 0, nmc = 0, ntotal = 0;
  // MCMC steps of one phi block; nu (not nu1B) is 20x at kiter 0
  uint64_t blockSteps(int kiter, int block) const {
    uint64_t f = kiter == 0 ? 20u : 1u;
    return f * (nu[0] + nu[1] + nu[2] * nphi[block]) + nu1B;
  }
  uint64_t steps(int kiter) const {
    return blockSteps(kiter, 0) + blockSteps(kiter, 1);
  }
  uint64_t stride(int kiter) const {
    return nComp * (steps(kiter) * nM + nmc * ntotal);
  }
  uint64_t iterBase(int kiter) const {
    return kiter == 0 ? 0u : stride(0) + (uint64_t)(kiter - 1) * stride(1);
  }
  // first seed of a step; block 0 is phi1, 1 is phi0; method 4 is mode 1B
  uint64_t step(int kiter, int comp, int block, int method, int u, int k1) const {
    uint64_t f = kiter == 0 ? 20u : 1u;
    uint64_t s = (uint64_t)comp * steps(kiter);
    if (block == 1) s += blockSteps(kiter, 0);
    if (method == 1) s += (uint64_t)u;
    else if (method == 2) s += f * nu[0] + (uint64_t)u;
    else if (method == 3) s += f * (nu[0] + nu[1]) + (uint64_t)u * nphi[block] + (uint64_t)k1;
    else s += f * (nu[0] + nu[1] + nu[2] * nphi[block]) + (uint64_t)u;
    return iterBase(kiter) + s * nM;
  }
  // first seed of chain k's observations
  uint64_t cens(int kiter, int comp, int k) const {
    return iterBase(kiter) + nComp * steps(kiter) * nM +
      ((uint64_t)comp * nmc + (uint64_t)k) * ntotal;
  }
};

// Draw one MCMC step, chain row r (k*N + subject) from seed + offset + r: its
// noise columns, any `extra` normals, then its acceptance uniform.
static inline void _saemDrawRows(int seed, uint64_t offset, arma::mat &noise,
                                 arma::vec &accU, arma::mat *extra = nullptr) {
  for (arma::uword r = 0; r < noise.n_rows; ++r) {
    nmSeqSeedSet(seed, offset, r);
    for (arma::uword c = 0; c < noise.n_cols; ++c) noise(r, c) = rxNormEng(0.0, 1.0);
    if (extra != nullptr) {
      for (arma::uword c = 0; c < extra->n_cols; ++c) (*extra)(r, c) = rxNormEng(0.0, 1.0);
    }
    accU(r) = rxUnifEng(0.0, 1.0);
  }
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
// Theta-sensitivity peer state for the non-mu gradient refinement
// (src/nonMuThetaGrad.h).  Separate from the phi1 pool state: it is declared
// for any model shape that produced a sensitivity model, general-likelihood or
// not.
extern bool _saemThetaSensAnalytic;
extern int _saemEtaDistCppMap;
extern int _saemEtaDistRMap;
extern bool _saemThetaSensActive;
extern arma::ivec _saemThetaSensPhi0Col;  // sens output -> phi0 column, -1 = none
extern arma::ivec _saemThetaSensTheta;    // sens output -> 1-based ntheta
extern int _saemNonMuGradEvery;
extern int _saemThetaSensPredOffset;      // lhs index of rx_pred_ in the peer
extern int _saemThetaSensROffset;         // lhs index of rx_r_, -1 when absent
// THETA[k]/ETA[k] -> phi column translation for driving the peer from SAEM's
// own phi matrix, the same shape the phi1 peers use.
extern arma::ivec _saemThetaSensThetaKind;     // 1 = phi1, 0 = phi0, -1 = fixed
extern arma::ivec _saemThetaSensThetaCol;      // column within that group
extern arma::vec _saemThetaSensThetaFixedVal;  // value for kind -1
extern arma::ivec _saemThetaSensEtaCol;        // ETA[k] -> phi1 column
extern int _saemThetaSensDvCol;                // -1 when DV is not a parameter
extern int _saemThetaSensSensOffset;           // lhs index of the FIRST d(f)/d(theta)
// lhs index of EACH d(f)/d(theta) output, resolved by name rather than assumed
// contiguous from the first -- the outputs are named by ntheta, and the map's
// order is .impmapEstTheta()$all, which need not be the lhs order.
extern arma::ivec _saemThetaSensSensIx;
extern int _saemThetaSensNlhs;                 // peer's lhs width (sizes the buffer)

extern bool _saemPhi1PoolReady;
extern bool _saemPhi1UseAnalyticHess;
extern arma::ivec _saemPhi1H2ThetaKind;
extern arma::ivec _saemPhi1H2ThetaCol;
extern arma::vec _saemPhi1H2ThetaFixedVal;
extern arma::ivec _saemPhi1H2EtaCol;
arma::ivec _saemPhi1EtaNonMu;
// 1 when ETA[k]'s parameter has no THETA[] of its own (a nonMuEta, e.g. a
// dist()-declared eta): the pooled setter must put the phi VALUE there, not 0.
extern arma::ivec _saemPhi1EtaNonMu;
extern bool _saemPhi1WantHessian;
extern int _saemPhi1PredOffset;
// Slot SAEM's OWN per-iteration solve goes through, and rx_pred_'s lhs index in
// it.  See the definitions further down for why this is the sensitivity slot
// when that peer is live.
extern int _saemOwnSolveSlot;
extern int _saemOwnPredOffset;
extern int _saemSolveCompleteOnce;
extern int _saemReadSlot;
static void saemPickOwnSolveSlot();
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

// A C++ exception must NEVER escape an OpenMP region: the runtime cannot unwind
// out of one, so it calls std::terminate and the whole R process aborts.
//
// rxode2ll's likelihood functions throw std::domain_error for an out-of-domain
// argument -- `student_t_lpdf: Degrees of freedom parameter is -598.433, but
// must be positive finite` -- and a general-likelihood SAEM fit reaches exactly
// that as a matter of course: the MCMC proposes a negative df, the model's
// calc_lhs evaluates ll() there, and it throws from inside the threaded pooled
// solve.  Serially that would be a caught error; in the parallel region it was
// an abort with no R-level message.
//
// A parameter the likelihood is not defined at is precisely what a bad solve
// already means, and every caller here handles that by rejecting the proposal.
// So catch it and report it as one, rather than letting it kill the session.
template <typename F>
static inline bool saemNoThrow(F &&f) {
  try { f(); return true; } catch (...) { return false; }
}

// class def starts
class SAEM {
  typedef mat (*user_funct) (const mat&, const mat&, const List&);

public:

  SAEM() {
    // Point the harvest's file-static at this fit.  It used to be set inside the
    // declared-distribution branch of inits(), so a later fit on a model with no
    // declaration left it pointing at a destroyed object -- measured as
    // "double free or corruption (out)" part-way through a test file.
    gAnchorSelf = this;
    user_fn = NULL;
  }

  ~SAEM() {
    if (gAnchorSelf == this) gAnchorSelf = nullptr;
  }

  // Total observation -log-likelihood at candidate fixed-effect (phi0) values p,
  // holding the current phi1 samples fixed (general-likelihood / distribution==4:
  // the model prediction column is the per-observation log-likelihood).  Summed
  // over all chains, which is the SAEM stochastic-approximation objective.
  // The copula correlation's SUFFICIENT STATISTIC, every iteration.
  //
  // S_n = L S_z L' with L from the current rho, so the normalized off-diagonal
  //
  //     rho_hat = S_n[j,k] / sqrt(S_n[j,j] * S_n[k,k])
  //
  // is the M-step for a unit-diagonal correlation.  Uncentred: the latent's
  // PRIOR mean is zero, so a shifted sample is information, not nuisance.
  //
  // ESTIMABILITY.  All of the information about rho lives in the departure of
  // S_z from the identity -- if S_z = I then S_n = L L' = R_old and the update
  // returns whatever it was handed, for any estimator.  So the off-diagonal is
  // compared against the sampling noise of a zero correlation, 2/sqrt(N):
  // below that, the data do not identify this correlation and the fit should
  // say so rather than report whatever the search drifted to.
  void etaDistCorSuffStat(unsigned int kiter, const vec &pas) {
    if (etaDistNdist <= 0 || (int)etaDistCorSuff.n_elem != etaDistNdist) return;
    bool tr = (getenv("NLMIXR2_ETADIST_OPT") != NULL);
    for (int k = 0; k < etaDistNdist; ++k) {
      int j = etaDistCorWith(k);
      if (j < 0) { etaDistCorEstim(k) = -1; continue; }
      int cj = etaDistLatent(j), ck = etaDistLatent(k);
      if (cj < 0 || ck < 0 || cj >= (int)phiM.n_cols || ck >= (int)phiM.n_cols) {
        etaDistCorEstim(k) = -1; continue;
      }
      double zjj = 0, zkk = 0, zjk = 0;
      unsigned int nr = phiM.n_rows;
      // On the DIRECT route phiM holds the ETAS, not the latents, and this
      // statistic is written for latents throughout -- it standardizes them and
      // compares S_z against the identity.  Two gammas have neither unit
      // variance nor zero mean, so the raw second moments would report
      // structure that is the MARGINALS and not the dependence.
      //
      // z = qnorm(F_k(eta_k)) is the latent the copula is defined on, which is
      // exactly what rxEtaDistCopulaZ() computes for the prior term, so the
      // same transform serves both and they cannot drift apart.
      const bool edEtaScale = etaDistDirectOn();
      const int naC = (int)etaDistArgs.n_cols;
      std::vector<double> aJ((size_t)naC), aK((size_t)naC);
      if (edEtaScale) {
        for (int t = 0; t < naC; ++t) {
          aJ[(size_t)t] = etaDistArgs(j, t); aK[(size_t)t] = etaDistArgs(k, t);
        }
      }
      for (unsigned int r = 0; r < nr; ++r) {
        double a = phiM(r, cj), b = phiM(r, ck);
        if (edEtaScale) {
          a = rxEtaDistCopulaZ(etaDistFam(j), a, &aJ[0]);
          b = rxEtaDistCopulaZ(etaDistFam(k), b, &aK[0]);
        }
        if (!std::isfinite(a) || !std::isfinite(b)) continue;
        zjj += a*a; zkk += b*b; zjk += a*b;
      }
      if (nr < 2 || !(zjj > 0) || !(zkk > 0)) { etaDistCorEstim(k) = -1; continue; }
      double nn = (double)nr;
      double rho0 = etaDistRho(k);
      if (!std::isfinite(rho0)) rho0 = 0.0;
      double l21 = rho0;
      // STANDARDIZE S_z first.  The raw second moments carry the latent's
      // over-dispersion (diag ~2-3 rather than 1), and feeding that through
      // L S_z L' biases the ratio: njj uses zjj alone while nkk is a MIXTURE of
      // zjj, zjk and zkk, so unequal inflated diagonals inflate rho_hat.  On
      // Bauer's g1 that drove rho to +0.999, which collapses the partner's
      // latent onto its partner's, fails the family M-step's own spread guard,
      // and froze every family theta at its starting value (MARE 41.2% against
      // the baseline's 18.3%).
      //
      // Using the CORRELATION of z instead leaves a well-behaved update with
      // the right fixed point:
      //
      //     rho_hat = (l21 + l22*rz) / sqrt(1 + 2*l21*l22*rz)
      //
      // rz = 0 returns rho unchanged -- which is exactly correct, since S_z = I
      // means the current rho already explains the draws.  rz > 0 moves it up,
      // rz < 0 down, and neither can be driven by the diagonal any more.
      double rz = zjk/std::sqrt(zjj*zkk);
      // ONE copy of this update, shared with imp and the FOCEi family
      // (rxEtaDistCorFromRz, src/etaDistFam.h).  It used to live only here,
      // while imp and focei correlated the COMBINED latents instead -- a
      // different estimator with a fixed point at the current rho, which is
      // how their copula came back frozen at its ini() value.
      double rh = rxEtaDistCorFromRz(l21, rz);
      if (std::isfinite(rh)) {
        etaDistCorSuff(k) = rh;
      }
      double off = std::fabs(zjk/nn);
      etaDistCorOffMag(k) = off;
      // 2 SE of a zero correlation on the standardized scale
      double thresh = 2.0/std::sqrt(nn);
      if ((int)etaDistCorOffMax.n_elem == etaDistNdist && off > etaDistCorOffMax(k)) {
        etaDistCorOffMax(k) = off;
      }
      // estimable if it was EVER informative, not if it is informative right now
      double peak = ((int)etaDistCorOffMax.n_elem == etaDistNdist) ?
        etaDistCorOffMax(k) : off;
      etaDistCorEstim(k) = (peak > thresh) ? 1 : 0;
      // THE DAMPED ADJUSTMENT, every iteration, on the SAME stochastic-
      // approximation series pas(kiter) every other parameter uses.  The
      // statistic is second moments over phiM and costs nothing, so there is no
      // reason to hold the adjustment back to the M-step's cadence -- and a
      // parameter that moves once every 20 iterations while the SA weight is
      // shrinking barely moves at all.
      //
      // Only for method 3, which IS this statistic; the other estimators are
      // applied in the copula loop on their own terms.  And only when the data
      // identify it -- adjusting toward a statistic that is within noise of zero
      // is how a correlation ends up reporting the search's drift.
      if (etaDistCorMethod == 3 && etaDistCorEstim(k) == 1 &&
          std::isfinite(etaDistCorSuff(k))) {
        // note this keeps adjusting after the off-diagonal has shrunk: that is
        // correct, the statistic is still the M-step, and pas(kiter) is what
        // makes the adjustment fade rather than an estimability test
        double cur = etaDistRho(k);
        if (!std::isfinite(cur)) cur = 0.0;
        double v = cur + pas(kiter) * (etaDistCorSuff(k) - cur);
        if (std::isfinite(v)) {
          if (v > 0.99) v = 0.99; else if (v < -0.99) v = -0.99;
          etaDistRho(k) = v;
          int cc = corCol(k);
          if (cc >= 0) {
            double a = std::atanh(v);
            if (std::isfinite(a)) {
              mprior_phi0.col(cc).fill(a);
              std::vector<int> one(1, cc);
              writeBackPhi0(one);
            }
          }
          etaDistCorFired = true;
        }
      }
      if (tr && (kiter % 20 == 0)) {
        RSprintf("[suff] it=%d k=%d S_z/N jj=%.4f kk=%.4f jk=%+.4f rz=%+.4f | "
                 "rho cur=%+.4f suff=%+.4f | off=%.4f thresh=%.4f estimable=%d\n",
                 (int)kiter, k, zjj/nn, zkk/nn, zjk/nn, rz, rho0,
                 etaDistCorSuff(k), off, thresh, (int)etaDistCorEstim(k));
      }
    }
  }

  // Has family k's pooled latent spread STOPPED CHANGING between M-step
  // attempts?  Returns false on the first attempt (nothing to compare with) --
  // one skipped attempt, and it is what makes a chain that never settles never
  // update.
  //
  // MUST be called for every family on every attempt, or the comparison stops
  // spanning one etaDistEvery gap.  The caller therefore may not reach it
  // through a short-circuiting && (see the correlation loop).
  bool etaDistSpreadSettled(int k, double lsd) {
    if (k >= (int)etaDistSdPrev.n_elem) return false;
    // stage this attempt's value; the baseline advances once per iteration, at
    // the end of the M-step.  Advancing it inside a loop instead would let the
    // correlation loop compare an attempt against ITSELF -- relative change
    // zero, so always "settled", which is the opposite of the intent.
    etaDistSdCur(k) = lsd;
    if (!std::isfinite(lsd)) {
      // a measurement that failed is not a gap to be spanned: drop the
      // baseline so the next attempt declines for want of one, rather than
      // silently comparing across two gaps or more
      etaDistSdPrev(k) = NA_REAL;
      return false;
    }
    if (!(lsd >= etaDistSdLo && lsd <= etaDistSdHi)) return false;
    if (!(etaDistSdTol > 0.0)) return true;          // cap only
    double p = etaDistSdPrev(k);
    if (!std::isfinite(p) || !(p > 0.0)) return false;
    return std::fabs(lsd - p) <= etaDistSdTol * p;
  }

  // how many declared correlations the data do not identify
  int etaDistCorNotEstimable() const {
    int n = 0;
    for (int k = 0; k < (int)etaDistCorEstim.n_elem; ++k)
      if (etaDistCorEstim(k) == 0) n++;
    return n;
  }

  // ONE complete-system solve per step, shared.
  //
  // _saemSolveCompleteOnce makes this solve carry rx_pred_ and every
  // d(f)/d(theta) together, so the same one serves the gradient step and the
  // single-solve pred.  Re-solving is only needed when something has MOVED
  // phi0 since -- the writers call invalidateCompleteSolve() -- or when the
  // regressor search runs, which evaluates candidates and therefore solves per
  // candidate by nature.
  void ensureCompleteSolve(unsigned int kiter, bool wantSens) {
    if (_saemCompleteSolveIter == (int)kiter &&
        (_saemCompleteSolveHasSens || !wantSens)) return;
    if (nphi0 > 0) phiM.cols(i0) = repmat(mprior_phi0, nmc, 1);
    _saemSolveCompleteOnce = wantSens ? 1 : 0;
    { mat _tmp = user_fn(phiM, evt, optM); (void)_tmp; }
    _saemSolveCompleteOnce = 0;
    _saemCompleteSolveIter = (int)kiter;
    _saemCompleteSolveHasSens = wantSens;
  }
  // Brent's method: parabolic interpolation with golden-section fallback,
  // the algorithm behind R's optimize().  Self-contained because Brent_fmin is
  // not in R's public headers.
  //
  // Minimizes the OBSERVATION objective directly rather than root-finding on
  // its score.  The score route was tried and is not usable here: its values
  // came back at 1e+04 to 1e+09 on this parameter and the sign change it
  // bracketed sat at rho = -0.92 against a truth of +0.5.  Minimizing the
  // objective cannot pick a minimum of the log-likelihood by accident, and does
  // not depend on the score's scale or sign convention at all.
  double brentMinPhi0Col(int col, double ax, double bx, double tol,
                         int maxIt, std::vector<double> &pv, bool &ok) {
    const double gold = 0.5*(3.0 - std::sqrt(5.0));
    double a = ax, b = bx;
    double x = a + gold*(b - a), w = x, v = x;
    double d = 0.0, e = 0.0;
    pv[(size_t)col] = x;
    double fx = phi0Objective(pv.data());
    ok = std::isfinite(fx) && fx < 1e299;
    if (!ok) return x;
    double fw = fx, fv = fx;
    for (int it = 0; it < maxIt; ++it) {
      double xm = 0.5*(a + b);
      double tol1 = tol*std::fabs(x) + 1e-10, tol2 = 2.0*tol1;
      if (std::fabs(x - xm) <= tol2 - 0.5*(b - a)) break;
      bool useGold = true;
      if (std::fabs(e) > tol1) {
        double r = (x - w)*(fx - fv), q = (x - v)*(fx - fw);
        double pq = (x - v)*q - (x - w)*r;
        q = 2.0*(q - r);
        if (q > 0) pq = -pq; else q = -q;
        if (std::fabs(pq) < std::fabs(0.5*q*e) && pq > q*(a - x) && pq < q*(b - x)) {
          double etmp = e; e = d; d = pq/q; useGold = false;
          if (d - (b - x) > 0 || x + d - a < tol2) d = (xm >= x) ? tol1 : -tol1;
          (void)etmp;
        }
      }
      if (useGold) { e = (x >= xm) ? (a - x) : (b - x); d = gold*e; }
      double u = x + ((std::fabs(d) >= tol1) ? d : ((d > 0) ? tol1 : -tol1));
      pv[(size_t)col] = u;
      double fu = phi0Objective(pv.data());
      if (!std::isfinite(fu) || fu >= 1e299) break;
      if (fu <= fx) {
        if (u >= x) a = x; else b = x;
        v = w; fv = fw; w = x; fw = fx; x = u; fx = fu;
      } else {
        if (u < x) a = u; else b = u;
        if (fu <= fw || w == x) { v = w; fv = fw; w = u; fw = fu; }
        else if (fu <= fv || v == x || v == w) { v = u; fv = fu; }
      }
    }
    return x;
  }

  // persist a phi0 move through MCOV0, per column against its own design
  // block: one least squares over all of COV0 is rank deficient when nphi0 > 1
  void writeBackPhi0(const std::vector<int> &ix) {
    for (size_t fi = 0; fi < ix.size(); ++fi) {
      int c = ix[fi];
      uvec li = arma::find(LCOV0.col(c) == 1);
      if (li.n_elem == 0) continue;
      mat Xc = COV0.cols(li);
      vec bc;
      if (arma::solve(bc, Xc.t() * Xc, Xc.t() * mprior_phi0.col(c))) {
        for (unsigned int j = 0; j < li.n_elem; ++j) MCOV0(li(j), c) = bc(j);
      }
    }
  }
  void invalidateCompleteSolve() {
    _saemCompleteSolveIter = -1; _saemCompleteSolveHasSens = false;
  }

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

  // Population prediction vector at candidate free phi0 values, for the
  // Shi-difference fallback (gShiPredFn).  Writes the candidates into the phi0
  // columns of mprior_phi0/phiM and solves; the caller restores phi0.
  arma::vec shiPredAt(const arma::vec &t, const std::vector<int> &freeIx) {
    for (size_t j = 0; j < freeIx.size(); ++j) {
      int c = freeIx[j];
      if (c < 0 || c >= nphi0) continue;
      if (!std::isfinite(t((arma::uword)j))) return arma::vec();
      mprior_phi0.col(c).fill(t((arma::uword)j));
    }
    if (nphi0 > 0) phiM.cols(i0) = repmat(mprior_phi0, nmc, 1);
    mat fMat = user_fn(phiM, evt, optM);
    return arma::vec(fMat.col(0));
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
      // and the WITHIN-group off-diagonals, which carry Omega in the collapsed
      // form; equal diagonals plus equal off-diagonals is what makes the block
      // compound-symmetric
      double o = 0.0;
      unsigned int no = 0;
      for (unsigned int a = 0; a < ix.n_elem; ++a) {
        for (unsigned int b = a + 1; b < ix.n_elem; ++b) {
          o += G(ix(a), ix(b));
          no++;
        }
      }
      if (no == 0) continue;
      o /= (double)no;
      for (unsigned int a = 0; a < ix.n_elem; ++a) {
        for (unsigned int b = a + 1; b < ix.n_elem; ++b) {
          G(ix(a), ix(b)) = o;
          G(ix(b), ix(a)) = o;
        }
      }
    }
  }

  // Impose the collapsed form's shared-mean constraint on the GLS solution.
  // Exact rather than a projection: see omegaPoolMean.
  void poolLambdaGroups(vec &P) {
    if (!omegaPoolMean) return;
    if (omegaPool.n_elem != (unsigned int)nphi1) return;
    if (lambdaCol1.n_elem != P.n_elem) return;
    unsigned int maxg = omegaPool.max();
    for (unsigned int g = 1; g <= maxg; ++g) {
      uvec cols = find(omegaPool == g);
      if (cols.n_elem < 2) continue;
      std::vector<unsigned int> li;
      for (unsigned int l = 0; l < P.n_elem; ++l) {
        if (lambdaCol1(l) < (unsigned int)nphi1 &&
            omegaPool(lambdaCol1(l)) == g) li.push_back(l);
      }
      if (li.size() < 2) continue;
      double m = 0.0;
      for (size_t k = 0; k < li.size(); ++k) m += P(li[k]);
      m /= (double)li.size();
      for (size_t k = 0; k < li.size(); ++k) P(li[k]) = m;
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

  // Exact-gradient warm start for the non-mu (phi0) thetas.
  //
  // refinePhi0Lik's search has no derivative information and a budget of
  // nonMuThetaMaxEval full-population solves.  The peer model emits exact
  // symbolic d(f)/d(theta) for every estimated theta in ONE solve per row, so a
  // Gauss-Newton step off it is cheaper than a handful of search evaluations
  // and far better directed.  This does NOT replace the search: it moves
  // mprior_phi0 along the local quadratic model, and the search then runs
  // warm-started from there and corrects wherever that model was poor
  // (src/nonMuThetaGrad.h -- the hedge between linearity and non-linearity).
  //
  // The objective differentiated here is exactly phi0NormalSSR's, so the two
  // halves agree on what they are minimizing:
  //     0.5*((y - f)/g)^2 + log(g),   g = ares(b) + bres(b)*|f|
  // g and dg/df come from SAEM's OWN live ares/bres, not from the peer's
  // rx_r_.  SAEM keeps the residual error outside phi, so the peer's residual
  // THETA is pinned at its ini() value and its rx_r_ would weight the score
  // with a stale sd while SAEM's actual one moves.
  //
  // Returns true when it actually moved something.
  // xEval != nullptr puts this in EVALUATE mode: write those free-coordinate
  // values into mprior_phi0, take ONE solve of the complete system, and hand
  // back both the objective and its exact gradient without stepping.  That is
  // precisely what a quasi-Newton wants, and it is only affordable because the
  // complete system emits rx_pred_ alongside d(f)/d(theta) -- a derivative-free
  // search pays a solve per objective evaluation and gets no gradient at all.
  // Shi (2021) finite-difference gradient of the ORIGINAL (no-sensitivity)
  // model, used when the analytic path's bad-solve ladder is exhausted.
  //
  // Differentiates the population prediction vector with shi21Forward, one free
  // phi0 coordinate at a time, then walks the observations and feeds the SAME
  // accumulator the analytic path uses (nonMuGradAccumObs).  That is what makes
  // the two paths interchangeable: score and information come out with
  // identical semantics, including the per-subject BHHH form, so nothing
  // downstream needs to know which one produced them.
  //
  // Returns false when the fallback itself cannot produce a usable gradient, in
  // which case the caller leaves its thetas to the derivative-free search.
  bool shiGradPhi0(nonMuObjKind objKind, int nFree,
                   const std::vector<int> &obsOff, const arma::uvec &invSort,
                   arma::vec &score, arma::mat &info,
                   std::vector<double> &rowScore, std::vector<double> &rowInfo,
                   int nRow) {
    if (nFree <= 0 || nphi0 <= 0) return false;
    gShiSelf = this;
    gShiFreeIx = gPhi0FreeIx;
    arma::vec t((arma::uword)nFree);
    for (int fi = 0; fi < nFree; ++fi) t(fi) = mprior_phi0(0, gPhi0FreeIx[(size_t)fi]);
    // Snapshot phi0 so a probe cannot leave the model at a perturbed value.
    arma::rowvec phi0Save = mprior_phi0.row(0);
    bool frz = _saemFreezeOde;
    _saemFreezeOde = false;
    arma::vec f0 = gShiPredFn(t, 0);
    bool ok = f0.is_finite() && f0.n_elem == (arma::uword)(nmc * ntotal);
    std::vector< arma::vec > gr((size_t)nFree);
    for (int fi = 0; fi < nFree && ok; ++fi) {
      double h = 0.0;
      arma::vec g1;
      // shi21Forward picks the step from the objective's own noise floor; the
      // defaults are focei's.  A non-finite result for any coordinate makes the
      // whole fallback unusable, for the same reason a partial population is:
      // a gradient missing one direction is not a descent direction.
      double rc = shi21Forward(gShiPredFn, t, h, f0, g1, 0, fi);
      (void)rc;
      if (!g1.is_finite() || g1.n_elem != f0.n_elem) ok = false;
      else gr[(size_t)fi] = g1;
    }
    // Restore phi0 and the caller's solve state before returning either way.
    for (int c = 0; c < nphi0; ++c) mprior_phi0.col(c).fill(phi0Save(c));
    if (nphi0 > 0) phiM.cols(i0) = repmat(mprior_phi0, nmc, 1);
    { mat _t = user_fn(phiM, evt, optM); (void)_t; }
    _saemFreezeOde = frz;
    if (!ok) return false;

    score.zeros(nFree);
    info.zeros(nFree, nFree);
    std::fill(rowScore.begin(), rowScore.end(), 0.0);
    std::fill(rowInfo.begin(), rowInfo.end(), 0.0);
    std::vector<double> dfdth((size_t)nFree);
    // phiM stacks the chains, so row r = k*N + i (the same layout phi.slice(k)
    // reads).  Within one chain the prediction vector is in SOLVE order, which
    // obsOff indexes by subject and invSort maps to ys/ix_endpnt order --
    // exactly the convention the analytic loop established.
    for (int k = 0; k < nmc; ++k) {
      for (int i = 0; i < N; ++i) {
        int r = k * N + i;
        if (r >= nRow) continue;
        arma::vec sc((arma::uword)nFree, fill::zeros);
        arma::mat inf((arma::uword)nFree, (arma::uword)nFree, fill::zeros);
        for (int q = obsOff[(size_t)i]; q < obsOff[(size_t)i + 1]; ++q) {
          arma::uword idx = (arma::uword)(k * ntotal + q);
          double f = f0(idx);
          if (!std::isfinite(f)) return false;
          double y = 0.0, gsd = 0.0, dgsdf = 0.0;
          if (objKind != nonMuObjLl) {
            if (q >= ntotal) return false;
            arma::uword tt = invSort((arma::uword)q);
            y = ys(tt);
            int b = (nendpnt == 1) ? 0 : (int)ix_endpnt(tt);
            if (!std::isfinite(y)) return false;
            gsd = ares(b) + bres(b) * std::fabs(f);
            if (!(gsd > 0.0) || !std::isfinite(gsd)) continue;
            dgsdf = bres(b) * ((f < 0.0) ? -1.0 : 1.0);
          }
          for (int fi = 0; fi < nFree; ++fi) dfdth[(size_t)fi] = gr[(size_t)fi](idx);
          nonMuGradAccumObs(objKind, y, f, gsd, dgsdf,
                            dfdth.data(), nFree, 1.0, sc, inf);
        }
        for (int a = 0; a < nFree; ++a) {
          rowScore[(size_t)r * (size_t)nFree + (size_t)a] = sc(a);
          score(a) += sc(a);
        }
        if (nonMuThetaBhhh) {
          // same per-subject outer product the analytic path forms
          for (int a = 0; a < nFree; ++a)
            for (int bb = 0; bb < nFree; ++bb)
              info(a, bb) += sc(a) * sc(bb);
        } else {
          for (int a = 0; a < nFree; ++a)
            for (int bb = 0; bb < nFree; ++bb) {
              rowInfo[((size_t)r * (size_t)nFree + (size_t)a) * (size_t)nFree + (size_t)bb] =
                inf(a, bb);
              info(a, bb) += inf(a, bb);
            }
        }
      }
    }
    return score.is_finite() && info.is_finite();
  }

  bool nonMuGradPhi0(unsigned int kiter, const vec &pas,
                     const double *xEval = nullptr,
                     double *fOut = nullptr, double *gOut = nullptr) {
    const bool gchk = (getenv("NLMIXR2_SAEM_GRADCHECK") != NULL);
    if (gchk) Rprintf("gradPhi0 kiter=%u active=%d nFreeIx=%d nphi0=%d dist=%d nendpnt=%d\n",
                      kiter, (int)_saemThetaSensActive, (int)gPhi0FreeIx.size(),
                      nphi0, distribution, nendpnt);
    if (!_saemThetaSensActive || gPhi0FreeIx.empty()) return false;
    if (_saemNonMuGradEvery > 1 &&
        ((int)kiter % _saemNonMuGradEvery) != 0) return false;
    // Which objective this is differentiating MUST match the one the search
    // minimizes (phi0Objective), or the warm start pulls phi0 toward the argmin
    // of a different function.  For a general-likelihood model the prediction IS
    // the per-observation log-likelihood and phi0Objective is -sum(f); for a
    // normal model it is the Gaussian phi0NormalSSR.
    const nonMuObjKind objKind = (distribution == 4) ? nonMuObjLl : nonMuObjGauss;
    if (objKind == nonMuObjGauss) {
      // Transform-both-sides changes the objective's shape (phi0NormalSSR
      // applies _powerD to both f and y); the peer's sensitivities are on the
      // untransformed scale, so the chain rule below would be wrong.  Fall back
      // to the search.
      for (int b = 0; b < nendpnt; ++b) {
        if (yj(b) != 2 || lambda(b) != 1.0) return false;
      }
    }
    const int nFree = (int)gPhi0FreeIx.size();
    const int nTheta = (int)_saemThetaSensThetaKind.n_elem;
    const int nEta = (int)_saemThetaSensEtaCol.n_elem;
    const int nSens = (int)_saemThetaSensPhi0Col.n_elem;
    // sens output -> free-index position, -1 when that output is not a free
    // phi0 column (fixed, M-step-owned, or not a phi0 column at all)
    std::vector<int> sensFree((size_t)nSens, -1);
    bool any = false;
    for (int s = 0; s < nSens; ++s) {
      int c = _saemThetaSensPhi0Col(s);
      if (c < 0) continue;
      for (int fi = 0; fi < nFree; ++fi) {
        if (gPhi0FreeIx[(size_t)fi] == c) { sensFree[(size_t)s] = fi; any = true; break; }
      }
    }
    if (gchk) Rprintf("  nSens=%d anyFree=%d ownSlot=%d (thetaSens=%d pred=%d) ownPredOff=%d\n",
                      nSens, (int)any, _saemOwnSolveSlot, (int)odeSlotThetaSens,
                      (int)odeSlotPred, _saemOwnPredOffset);
    if (!any) return false;

    // Observation identity, in SOLVE order.
    //
    // This walk visits observations subject-by-subject in solve order, but SAEM
    // holds ys/ix_endpnt in endpoint-blocked ix_sorting order (phi0NormalSSR
    // reorders f the same way before it can index ys).  getIndDv() is NOT the
    // observation here -- it reads 0 in this path, which silently made every
    // residual -f and is what made the analytic score disagree with a finite
    // difference by seven orders of magnitude.
    //
    // ix_sorting maps sorted position t -> solve-order position q, so invert it
    // once to go the other way.  With the endpoint recovered per observation
    // there is no need to restrict this to a single endpoint either.
    arma::uvec invSort;
    if (objKind == nonMuObjGauss) {
      if ((int)ix_sorting.n_elem != ntotal || (int)ys.n_elem < ntotal) return false;
      invSort.set_size(ntotal);
      for (int t = 0; t < ntotal; ++t) {
        unsigned int q = ix_sorting(t);
        if ((int)q >= ntotal) return false;
        invSort(q) = (unsigned int)t;
      }
    }

    rx_solving_options *op = getSolvingOptions(_rx);
    int cores = getOpCores(op);
    bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
    const int nRow = N * nmc;
    // Where each subject's observations start within one chain's solve-order
    // block.  Chains replicate the same subjects, so this is computed once.
    std::vector<int> obsOff((size_t)N + 1, 0);
    if (objKind == nonMuObjGauss) {
      int acc = 0;
      for (int i = 0; i < N; ++i) {
        obsOff[(size_t)i] = acc;
        rx_solving_options_ind *indI = getSolvingOptionsInd(_rx, i);
        for (int j = 0; j < getIndNallTimes(indI); ++j) {
          if (getIndEvid(indI, getIndIx(indI, j)) == 0) acc++;
        }
      }
      obsOff[(size_t)N] = acc;
      if (acc != ntotal) return false;   // layout is not what this assumes
    }
    // Per-row score/information, reduced serially afterwards -- accumulating
    // into shared arma objects inside the parallel region would race.
    if (xEval != nullptr) {
      for (int fi = 0; fi < nFree; ++fi) {
        int c = gPhi0FreeIx[(size_t)fi];
        if (!std::isfinite(xEval[fi])) return false;
        mprior_phi0.col(c).fill(xEval[fi]);
      }
      if (nphi0 > 0) phiM.cols(i0) = repmat(mprior_phi0, nmc, 1);
      bool frz = _saemFreezeOde;
      _saemFreezeOde = false;
      _saemSolveCompleteOnce = 1;
      mat fMat = user_fn(phiM, evt, optM);
      _saemSolveCompleteOnce = 0;
      _saemFreezeOde = frz;
      double v = (distribution == 4) ? -accu(fMat.col(0))
                                     : phi0NormalSSR(fMat.col(0));
      if (!std::isfinite(v)) return false;
      if (fOut != nullptr) *fOut = v;
    }
    std::vector<double> rowScore((size_t)nRow * (size_t)nFree, 0.0);
    std::vector<double> rowInfo((size_t)nRow * (size_t)nFree * (size_t)nFree, 0.0);
    std::vector<int> rowBad((size_t)nRow, 0);
    // The event-sensitivity shape is a process global installed only by
    // OdeSwapEsBatch, which MUST be constructed outside the OpenMP region --
    // and a "no ES" slot still needs one built, so a shape left installed by an
    // earlier solve is deactivated rather than reused with the wrong dimensions
    // (OdeSwapEsBatch's own contract, src/odeSwap.cpp).
    // Only when we are going to SOLVE.  This object installs/deactivates the
    // process-wide event-sensitivity shape for a solve; constructing it when we
    // are merely reading a solve someone else took can only disturb it.
    std::unique_ptr<OdeSwapEsBatch> tsEsBatch;
    if (_saemOwnSolveSlot != odeSlotThetaSens) {
      tsEsBatch.reset(new OdeSwapEsBatch(odeSlotThetaSens));
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
    for (int r = 0; r < nRow; ++r) {
#ifdef _OPENMP
      if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
      int subj = r % N;
      (void)subj;
      rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, r);
      OdeSwapScope neqGuard(odeSlotThetaSens, ind, op);
      OdeSwapCmtScope cmtGuard(odeSlotThetaSens, op, ind);
      // Set the peer exactly the way saemSetRowsPooled sets its own: the WHOLE
      // combined phi value goes in THETA[k] and ETA[k] is 0.  SAEM's phi is
      // already theta+eta, and the model depends on the pair only through their
      // sum, so this gives the same f AND the same d(f)/d(THETA[k]) as splitting
      // it back into a mu and a deviation -- with none of the mprior_phi1
      // bookkeeping that split needed.
      //
      // The exception is a nonMuEta: no THETA[] refers to it, the parameter IS
      // the eta, so its phi value belongs in ETA[k] (same rule, and same flag,
      // as the pooled setter).
      for (int k = 0; k < nTheta; ++k) {
        double v;
        int kind = _saemThetaSensThetaKind(k), col = _saemThetaSensThetaCol(k);
        if (kind == 1) v = phiM(r, i1(col));
        else if (kind == 0) v = mprior_phi0(0, col);
        else v = _saemThetaSensThetaFixedVal(k);
        setIndParPtr(ind, k, v);
      }
      for (int k = 0; k < nEta; ++k) {
        double v = 0.0;
        if (k < (int)_saemPhi1EtaNonMu.n_elem && _saemPhi1EtaNonMu(k) != 0) {
          v = phiM(r, i1(_saemThetaSensEtaCol(k)));
        }
        setIndParPtr(ind, nTheta + k, v);
      }
      // REUSE the solve the caller just took.
      //
      // refinePhi0Lik establishes states with a full population solve
      // immediately before calling this, and when the sensitivity peer is live
      // that solve goes through odeSlotThetaSens -- the complete system, whose
      // lhs carries rx_pred_ AND every d(f)/d(theta) column.  Its parameters
      // are the ones set just above, because saemSetRowsPooled uses the very
      // same convention this loop does (THETA = combined phi, ETA = 0 except a
      // nonMuEta).  So re-solving here would integrate the identical system a
      // second time and double the per-iteration cost.  imp already takes this
      // route -- impThetaSensCollect's reuseSolve path (src/inner.cpp).
      //
      // Only solve when the caller's solve did NOT come through this slot.
      // DISABLED: reuse is measurably wrong as written -- the finite-difference
      // check went from agreeing to 1e-8 to reporting two columns as exactly 0
      // and the rest sign-flipped.  Something between the caller's solve and
      // this read does not survive, so the saving is not free the way the
      // parameter-convention argument suggested.  Left in place, and off, until
      // that is understood; re-solving is correct but pays for the system twice.
      // Reuse the caller's solve.  refinePhi0Lik establishes states with a full
      // population solve immediately before this, through _saemOwnSolveSlot --
      // which is odeSlotThetaSens whenever the peer is registered, i.e. the
      // COMPLETE system, whose lhs carries rx_pred_ and every d(f)/d(theta)
      // column.  Its parameters are the ones set just above, because
      // saemSetRowsPooled uses this loop's own convention (THETA = combined
      // phi, ETA = 0 except a nonMuEta).  Re-integrating would solve the
      // identical system twice per refinement iteration.
      bool reuse = (_saemReadSlot == odeSlotThetaSens);
      int _tolTries = 0;
      double _tol0 = getIndTolFactor(ind);
      // setIndSolve() selects which solve buffer this individual reads through,
      // so it is needed whether or not we re-integrate -- skipping it with the
      // solve is what made the reused read return zeros for some columns.
      setIndSolve(ind, -1);
      if (!reuse) {
        // Bad-solve ladder, following focei: retry with loosened tolerances up
        // to maxOdeRecalc times, each by odeRecalcFactor.  Skipped entirely
        // when the peer is analytic -- there is no integration for a tolerance
        // to affect, so a retry would just repeat the same answer.  The
        // subject's tolFactor is restored afterwards so a row that needed
        // loosening here does not silently loosen SAEM's own later solves.
        int maxTry = _saemThetaSensAnalytic ? 0 :
          ((current_saem_state != nullptr) ? current_saem_state->_saemMaxOdeRecalc : 0);
        double fac = (current_saem_state != nullptr) ?
          current_saem_state->_saemOdeRecalcFactor : 1.0;
        bool solved = false;
        while (true) {
          if (saemNoThrow([&]{ odeSwapSolveInd(odeSlotThetaSens, r); }) &&
              !odeSwapIndBadSolveSlot(op, ind, odeSlotThetaSens)) { solved = true; break; }
          if (_tolTries >= maxTry || !(fac > 1.0)) break;
          setIndTolFactor(ind, getIndTolFactor(ind) * fac);
          _tolTries++;
          setIndSolve(ind, -1);
        }
        if (_tolTries > 0) setIndTolFactor(ind, _tol0);
        if (!solved) { rowBad[(size_t)r] = 1; continue; }
      } else if (odeSwapIndBadSolveSlot(op, ind, odeSlotThetaSens)) {
        rowBad[(size_t)r] = 1; continue;
      }
      iniSubjectE(r, 1, ind, op, _rx, rxThetaSens.update_inis);
      double *lhs = neqGuard.lhs();
      arma::vec sc((int)nFree, fill::zeros);
      arma::mat inf((int)nFree, (int)nFree, fill::zeros);
      std::vector<double> dfdth((size_t)nFree);
      int nObs = 0;
      for (int j = 0; j < getIndNallTimes(ind); ++j) {
        setIndIdx(ind, j);
        int kk = getIndIx(ind, j);
        if (getIndEvid(ind, kk) != 0) continue;
        double curT = getTime(kk, ind);
        if (gchk && r == 0 && nObs == 0) {
          // Print the states this row actually HAS.  A fixed 10 walks off the
          // end of the solve buffer on any model with fewer (a 2-compartment
          // linCmt has 2), which aborted the check with an out-of-bounds throw
          // -- debug-only, but a trap for whoever next turns the check on.
          double *st = getOpIndSolve(op, ind, j);
          int nSt = getOpNeq(op);
          if (nSt > 10) nSt = 10;
          Rprintf("    [row0 obs0 reuse=%d states:", (int)reuse);
          for (int q = 0; q < nSt; ++q) Rprintf(" %.4g", st[q]);
          Rprintf("]\n");
        }
        if (!saemNoThrow([&]{
              rxThetaSens.calc_lhs(r, curT, getOpIndSolve(op, ind, j), lhs); })) {
          rowBad[(size_t)r] = 1; break;
        }
        double f = lhs[_saemThetaSensPredOffset];
        if (!std::isfinite(f)) { rowBad[(size_t)r] = 1; break; }
        double y = 0.0, gsd = 0.0, dgsdf = 0.0;
        if (objKind == nonMuObjGauss) {
          int q = obsOff[(size_t)subj] + nObs;   // solve-order position in the chain
          if (q >= ntotal) { rowBad[(size_t)r] = 1; break; }
          unsigned int t = invSort(q);           // its position in ys/ix_endpnt
          y = ys(t);
          int b = (nendpnt == 1) ? 0 : (int)ix_endpnt(t);
          if (!std::isfinite(y)) { rowBad[(size_t)r] = 1; break; }
          gsd = ares(b) + bres(b) * std::fabs(f);
          if (!(gsd > 0.0) || !std::isfinite(gsd)) { nObs++; continue; }
          dgsdf = bres(b) * ((f < 0.0) ? -1.0 : 1.0);
        }
        std::fill(dfdth.begin(), dfdth.end(), 0.0);
        bool okObs = true;
        for (int sIx = 0; sIx < nSens; ++sIx) {
          int fi = sensFree[(size_t)sIx];
          if (fi < 0) continue;
          int lix = (sIx < (int)_saemThetaSensSensIx.n_elem) ?
            _saemThetaSensSensIx(sIx) : -1;
          if (lix < 0) { okObs = false; break; }
          double d = lhs[lix];
          if (!std::isfinite(d)) { okObs = false; break; }
          dfdth[(size_t)fi] = d;
        }
        if (!okObs) { rowBad[(size_t)r] = 1; break; }
        nonMuGradAccumObs(objKind, y, f, gsd, dgsdf,
                          dfdth.data(), nFree, 1.0, sc, inf);
        nObs++;
      }
      if (rowBad[(size_t)r]) continue;
      for (int a = 0; a < nFree; ++a) {
        rowScore[(size_t)r * (size_t)nFree + (size_t)a] = sc(a);
        for (int bb = 0; bb < nFree; ++bb)
          rowInfo[((size_t)r * (size_t)nFree + (size_t)a) * (size_t)nFree + (size_t)bb] =
            inf(a, bb);
      }
    }
    tsEsBatch.reset();
    arma::vec score(nFree, fill::zeros);
    arma::mat info(nFree, nFree, fill::zeros);
    int nGood = 0;
    for (int r = 0; r < nRow; ++r) {
      if (rowBad[(size_t)r]) continue;
      nGood++;
      for (int a = 0; a < nFree; ++a) {
        score(a) += rowScore[(size_t)r * (size_t)nFree + (size_t)a];
      }
      if (nonMuThetaBhhh) {
        // NONMEM eq. 1.51: H = sum_i g_i g_i', the outer product of each
        // SUBJECT's score.  rowScore is already per (subject, chain), so the
        // subject-level gradient is in hand.  The default path instead sums
        // per-OBSERVATION outer products (rowInfo), which is a different matrix
        // -- it drops the within-subject correlation of the score and so
        // overstates the information.
        for (int a = 0; a < nFree; ++a) {
          double ga = rowScore[(size_t)r * (size_t)nFree + (size_t)a];
          for (int bb = 0; bb < nFree; ++bb) {
            info(a, bb) += ga * rowScore[(size_t)r * (size_t)nFree + (size_t)bb];
          }
        }
      } else {
        for (int a = 0; a < nFree; ++a)
          for (int bb = 0; bb < nFree; ++bb)
            info(a, bb) +=
              rowInfo[((size_t)r * (size_t)nFree + (size_t)a) * (size_t)nFree + (size_t)bb];
      }
    }
    // A partial population would bias the step toward whoever happened to
    // solve; the search alone is better than a skewed Newton step.
    if (gchk) Rprintf("  nGood=%d / nRow=%d\n", nGood, nRow);
    // NLMIXR2_SAEM_SHICHECK: compute the Shi fallback ALONGSIDE a healthy
    // analytic gradient and print both.  The fallback is otherwise reached only
    // when the sensitivity ladder is exhausted, which is rare and not
    // reproducible on demand -- so without this it ships untested.  The two are
    // meant to be the same gradient by different means; large disagreement is a
    // wrong index, a wrong solve slot, or a step search that never converged.
    if (getenv("NLMIXR2_SAEM_SHICHECK") != NULL && nGood == nRow && xEval == nullptr) {
      arma::vec sA = score;
      arma::mat iA = info;
      std::vector<double> rsA = rowScore, riA = rowInfo;
      arma::vec sS(nFree, fill::zeros); arma::mat iS(nFree, nFree, fill::zeros);
      std::vector<double> rsS = rowScore, riS = rowInfo;
      if (shiGradPhi0(objKind, nFree, obsOff, invSort, sS, iS, rsS, riS, nRow)) {
        Rprintf("saem shi-vs-analytic (kiter=%u)\n", kiter);
        for (int fi = 0; fi < nFree; ++fi) {
          double rel = (std::fabs(sA(fi)) > 1e-8) ?
            std::fabs(sS(fi) - sA(fi)) / std::fabs(sA(fi)) :
            std::fabs(sS(fi) - sA(fi));
          Rprintf("  phi0[%d] analytic=% .8e  shi=% .8e  rel=%.3e\n",
                  gPhi0FreeIx[(size_t)fi], sA(fi), sS(fi), rel);
        }
      } else {
        Rprintf("saem shi-vs-analytic (kiter=%u): fallback declined\n", kiter);
      }
      score = sA; info = iA; rowScore = rsA; rowInfo = riA;
    }
    if (nGood < nRow) {
      // The analytic sensitivities did not survive for the whole population,
      // and a partial one biases the step toward whoever happened to solve.
      // Fall back to the original model differenced with Shi (2021) steps --
      // focei's own answer in this situation.  Costs nFree+1 population solves,
      // paid only after the ladder has already failed.
      if (xEval != nullptr) return false;   // evaluate-at-x mode has no fallback
      if (!shiGradPhi0(objKind, nFree, obsOff, invSort, score, info,
                       rowScore, rowInfo, nRow)) return false;
      if (gchk) Rprintf("  shi fallback supplied the gradient\n");
      _saemShiFallbackN++;
    }
    if (xEval != nullptr) {
      if (gOut != nullptr) for (int fi = 0; fi < nFree; ++fi) gOut[fi] = score(fi);
      return true;
    }
    // Finite-difference verification of the exact gradient, off by default.
    // The analytic score above and a central difference of phi0Objective --
    // the very objective the search minimizes -- must agree; anything else is
    // a wrong index, a wrong sign, or a stale residual sd.  Kept behind an
    // env var so a normal fit never pays for the 2*nFree extra population
    // solves it costs.
    if (gchk) {
      // phi0Objective() re-solves only when the ODE is not frozen, and
      // refinePhi0Lik() establishes the solve states before it starts
      // optimizing.  This check runs BEFORE both of those, so without doing the
      // same here it would difference a stale or frozen solve and report a
      // mismatch that says nothing about the gradient.
      bool _frz = _saemFreezeOde;
      _saemFreezeOde = false;
      { mat _tmp = user_fn(phiM, evt, optM); (void)_tmp; }
      std::vector<double> pv((size_t)nphi0);
      for (int c = 0; c < nphi0; ++c) pv[(size_t)c] = mprior_phi0(0, c);
      Rprintf("saem non-mu gradient check (kiter=%u)\n", kiter);
      for (int fi = 0; fi < nFree; ++fi) {
        int c = gPhi0FreeIx[(size_t)fi];
        double x0 = pv[(size_t)c];
        double h = 1e-5 * std::max(1.0, std::fabs(x0));
        pv[(size_t)c] = x0 + h;
        double fp = phi0Objective(pv.data());
        pv[(size_t)c] = x0 - h;
        double fm = phi0Objective(pv.data());
        pv[(size_t)c] = x0;
        double fd = (fp - fm) / (2.0 * h);
        double rel = (std::fabs(fd) > 1e-8) ?
          std::fabs(score(fi) - fd) / std::fabs(fd) : std::fabs(score(fi) - fd);
        Rprintf("  phi0[%d] analytic=% .8e  fd=% .8e  rel=%.3e\n",
                c, score(fi), fd, rel);
      }
      _saemFreezeOde = _frz;
    }
    // Damped by the SA step exactly like the M-step's own update, and clamped
    // to refinePhi0Lik's local trust radius so a poorly conditioned information
    // matrix cannot throw a phi0 outside the region the search then works in.
    arma::vec cur(nFree), step;
    for (int fi = 0; fi < nFree; ++fi) cur(fi) = mprior_phi0(0, gPhi0FreeIx[(size_t)fi]);
    if (!nonMuGradStep(score, info, cur, 0.75, step)) return false;
    double damp = (kiter < pas.n_elem) ? pas(kiter) : 1.0;
    if (!std::isfinite(damp) || damp <= 0.0) return false;
    // NONMEM's alpha acceptance test (the text after eq. 1.46): try alpha = 1,
    // evaluate the objective at the proposed point, and if it did not improve
    // shrink alpha by sqrt(2) and try again.  This is the guard the default
    // path does not have: nonMuGradStep() gates on the CONDITIONING of the
    // information matrix, which says the step is numerically trustworthy, not
    // that it goes downhill.  On a nearly flat or degenerate direction those
    // are different questions, and only the second one stops a march to the
    // boundary.
    if (nonMuThetaBhhh) {
      std::vector<double> pv((size_t)nphi0);
      for (int c = 0; c < nphi0; ++c) pv[(size_t)c] = mprior_phi0(0, c);
      bool frz = _saemFreezeOde;
      _saemFreezeOde = false;
      double f0 = phi0Objective(pv.data());
      double alpha = 1.0;
      bool accepted = false;
      if (std::isfinite(f0)) {
        for (int t = 0; t < 8 && !accepted; ++t) {
          for (int fi = 0; fi < nFree; ++fi) {
            int c = gPhi0FreeIx[(size_t)fi];
            double v = cur(fi) + alpha * step(fi);
            double lo = ((int)phi0Lower.n_elem == nphi0) ? phi0Lower(c) : R_NegInf;
            double hi = ((int)phi0Upper.n_elem == nphi0) ? phi0Upper(c) : R_PosInf;
            if (std::isfinite(lo) && v < lo) v = lo;
            if (std::isfinite(hi) && v > hi) v = hi;
            pv[(size_t)c] = v;
          }
          double ft = phi0Objective(pv.data());
          if (std::isfinite(ft) && ft < f0) accepted = true;
          else alpha /= std::sqrt(2.0);
        }
      }
      // restore, so the damped write below is the only thing that moves phi0
      for (int c = 0; c < nphi0; ++c) pv[(size_t)c] = mprior_phi0(0, c);
      _saemFreezeOde = frz;
      if (!accepted) {
        if (gchk) Rprintf("  bhhh: no alpha improved the objective; theta held\n");
        return false;
      }
      step *= alpha;
    }
    bool moved = false;
    for (int fi = 0; fi < nFree; ++fi) {
      int c = gPhi0FreeIx[(size_t)fi];
      double v = cur(fi) + damp * step(fi);
      double lo = ((int)phi0Lower.n_elem == nphi0) ? phi0Lower(c) : R_NegInf;
      double hi = ((int)phi0Upper.n_elem == nphi0) ? phi0Upper(c) : R_PosInf;
      if (std::isfinite(lo) && v < lo) v = lo;
      if (std::isfinite(hi) && v > hi) v = hi;
      if (std::isfinite(v) && v != mprior_phi0(0, c)) {
        mprior_phi0(0, c) = v;
        moved = true;
      }
    }
    return moved;
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
    // A theta the distribution M-step owns must NOT also be optimized here --
    // two optimizers on the same parameter, against different objectives, is
    // exactly the fight this change exists to end.
    // Columns another mechanism owns, held out of this refinement.  Split the
    // same way the M-step is: the family thetas only when etaDistOn, the copula
    // theta whenever the closed form is running.  The closed form has ALREADY
    // written its damped value into mprior_phi0 this iteration, so excluding it
    // here is what makes that value the warm start rather than something the
    // search immediately overwrites.
    std::vector<bool> phi0Dist((size_t)nphi0, false);
    // Held out in BOTH modes, only the owner differs: the family M-step by
    // default, the one-solve gradient step (etaDistGradStep) under
    // etaDistLoglik.  Either way this SEARCH must not also move them -- it
    // pays a population solve per candidate to rediscover a derivative the
    // sensitivity peer already emits.
    if (etaDistOn && etaDistNdist > 0 &&
        (int)etaDistThetaPhi0.n_rows == etaDistNdist) {
      for (int k = 0; k < etaDistNdist; ++k) {
        // Only cede family k's thetas to the M-step once it has actually moved
        // them.  Otherwise this search keeps them: an M-step that never fires
        // must not leave its parameters unowned, which returns the ini()
        // values as if they were estimates.
        if ((int)etaDistFiredK.size() == etaDistNdist &&
            etaDistFiredK[(size_t)k] == 0) continue;
        for (int t = 0; t < etaDistNth(k) && t < (int)etaDistThetaPhi0.n_cols; ++t) {
          int c = etaDistPhi0Col(k, t);
          if (c >= 0 && c < nphi0) phi0Dist[(size_t)c] = true;
        }
      }
    }
    // the copula theta keyed off its OWN flag, not the family M-step's -- and,
    // like them, only once its closed form has actually written a value.  Mode
    // independent: the closed form runs in the observation-likelihood mode too
    // (a latent correlation is not a property of the mean function), so it owns
    // this column either way.
    // Held out whenever something else owns it: the closed form once it has
    // fired, or -- in the observation-likelihood mode -- the gradient step,
    // which now carries rxCor alongside the family thetas.
    if ((etaDistObsLik() && etaDistCorMethod == 2) ||
        (etaDistCorOn && etaDistCorFired)) {
      for (int k = 0; k < etaDistNdist; ++k) {
        int c = corCol(k);
        if (c >= 0) phi0Dist[(size_t)c] = true;
      }
    }
    // SCOPE, and a correction to why this is here.
    //
    // The reasoning it was written for was wrong: "unrestricted, this hands
    // refinePhi0Lik prop.sd as well".  It does not.  On Bauer's model phi0 is
    // {lclm, lv1m, lclrv, lv1rv, rxCor} -- nphi0 == 5, all five free --
    // because a residual-error parameter lives in ares/bres, not in phi0.  So
    // the free set was already the declared thetas plus the copula.
    //
    // It is also nearly unreachable: nonMuTheta defaults to "regress", so
    // nonMuThetaRegress is 1 on essentially every fit and this block is
    // skipped.  It applies only to nonMuTheta="eta", where refinePhi0Lik would
    // otherwise not run at all and etaDistLoglik is the only thing asking for
    // it -- there, restricting the free set to the declared thetas is what
    // keeps this control from quietly becoming "regress".
    //
    // Kept, narrowly, for that case.  Do not read it as the mechanism behind
    // any measured difference: with the default control it never executes.
    if (etaDistObsLik() && !nonMuThetaRegress && distribution != 4 &&
        (int)etaDistThetaPhi0.n_rows == etaDistNdist) {
      std::vector<bool> phi0Decl((size_t)nphi0, false);
      for (int k = 0; k < etaDistNdist; ++k) {
        for (int t = 0; t < etaDistNth(k) && t < (int)etaDistThetaPhi0.n_cols; ++t) {
          int c = etaDistPhi0Col(k, t);
          if (c >= 0 && c < nphi0) phi0Decl[(size_t)c] = true;
        }
      }
      for (int c = 0; c < nphi0; ++c) {
        if (!phi0Decl[(size_t)c]) phi0Dist[(size_t)c] = true;
      }
    }
    gPhi0FreeIx.clear();
    for (int c = 0; c < nphi0; ++c) {
      if (!phi0Fix[(size_t)c] && !phi0Dist[(size_t)c]) gPhi0FreeIx.push_back(c);
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
    // Refresh phiM's phi0 columns from mprior_phi0 before establishing states.
    // phiM is only rebuilt from mprior_phi0 at the END of the iteration, so
    // anything that moved phi0 earlier in THIS one (the distribution M-step,
    // say) leaves phiM carrying the previous value.  Without this the solve is
    // taken at a slightly stale phi0 while the gradient reads parameters at the
    // current one, and the two disagree: the finite-difference check sat at
    // ~5e-4 instead of ~1e-8, which is exactly that inconsistency and not
    // solver noise.
    // Sensitivities only when nonMuGradPhi0() below can actually consume them.
    // Without the peer it declines outright, and the SEARCH that follows is
    // derivative-free -- so the plain prediction solve is the right one and is
    // faster.
    ensureCompleteSolve(kiter, _saemThetaSensActive != 0);
    // Gauss-Newton warm start off the exact sensitivities, then the search
    // below refines from there (src/nonMuThetaGrad.h).  Placed AFTER
    // gPhi0FreeIx so it moves exactly the columns the search owns -- never one
    // the distribution M-step owns or the user fixed -- and AFTER the
    // establish-states solve above so it can READ that solve instead of taking
    // its own.  When the sensitivity peer is live that solve IS the complete
    // system (rx_pred_ and d(f)/d(theta) together, _saemOwnSolveSlot), and it
    // was going to happen anyway, so the gradient costs no solve at all.
    bool gradMoved = nonMuGradPhi0(kiter, pas);
    // The BHHH arm replaces the search, so it may only do that when the
    // gradient actually exists.  nonMuGradPhi0() declines outright whenever the
    // sensitivity peer is not registered (_saemThetaSensActive == 0) -- which
    // is the case for a declared-distribution linCmt model, among others -- and
    // returning here on that path would leave phi0 refined by NOTHING at all
    // rather than by one Newton step.  Measured: on Bauer's gamma model
    // nonMuThetaBhhh=TRUE was silently disabling phi0 refinement entirely, and
    // the numbers it produced were that, not a BHHH step.
    if (nonMuThetaBhhh && gradMoved) {
      // The BHHH step IS the update (NONMEM eqs. 1.47-1.52): one accepted
      // Newton step per iteration, averaged across iterations by the SA gain
      // (eq. 1.152).  No search follows it, which is the whole point -- what
      // the default path damps is the argmax of a full maximization, and a
      // full maximization is what walks a degenerate direction to its bound.
      // Still close the same way, so mprior_phi0 = COV0*MCOV0 next iteration.
      _saemFreezeOde = false;
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
      return;
    }
    bool doFreeze;
    if (distribution == 4) {
      // NOT unconditionally frozen: a general-likelihood model can still carry a
      // phi0 that drives the solve (an IOV magnitude theta), and freezing makes
      // the objective constant in it -- see phi0NeedsLiveSolve() and #1000.
      if (_phi0NeedsLive < 0) _phi0NeedsLive = phi0NeedsLiveSolve() ? 1 : 0;
      doFreeze = (_phi0NeedsLive == 0);
    } else if (phi0ObsLikRoute()) {
      if (_phi0OdeSensitive < 0) _phi0OdeSensitive = phi0AffectsOde() ? 1 : 0;
      doFreeze = (_phi0OdeSensitive == 0);
    } else {
      doFreeze = false;
    }
    _saemFreezeOde = doFreeze;
    if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
      static int _edOnce = 0;
      if (_edOnce++ < 2) {
        std::string fx;
        for (size_t q = 0; q < gPhi0FreeIx.size(); ++q)
          fx += std::to_string(gPhi0FreeIx[q]) + " ";
        RSprintf("[phi0] nphi0=%d nFree=%d free={%s} obsLikRoute=%d regress=%d "
                 "dist4=%d doFreeze=%d optType=%d thetaSensActive=%d\n",
                 nphi0, (int)gPhi0FreeIx.size(), fx.c_str(),
                 (int)phi0ObsLikRoute(), nonMuThetaRegress, (int)(distribution == 4),
                 (int)doFreeze, nonMuThetaOptType, (int)_saemThetaSensActive);
        // Which phi0 COLUMN each declared theta writes into.  The family M-step
        // writes back through etaDistThetaPhi0(k,t), so a wrong entry here puts
        // a family argument into a theta the declaration has no claim on -- and
        // that is invisible in the parameter table, which just shows an
        // unrelated theta having moved.
      }
    }
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
    bool localTrust = (distribution != 4) && phi0ObsLikRoute();
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
        bool n1qn1Done = false;
        if (nonMuThetaOptType == 3) {
          // n1qn1 (BFGS) on the EXACT gradient.  Each iteration costs one
          // complete solve returning objective and gradient together; the
          // derivative-free alternatives pay a solve per objective evaluation
          // and never see a derivative.  Requires the gradient -- without it
          // there is nothing to hand n1qn1, so fall through to newuoa.
          if (_saemThetaSensActive) {
            std::vector<double> x((size_t)nFree), gg((size_t)nFree, 0.0);
            for (int fi = 0; fi < nFree; ++fi) x[(size_t)fi] = par0[gPhi0FreeIx[(size_t)fi]];
            // n1qn1 is unbounded; clamp each accepted iterate back into the
            // trust region the same way the searches are clamped.
            std::vector<double> zm((size_t)(nFree*(nFree+13)/2 + 1), 0.0);
            std::vector<double> var((size_t)nFree, 0.1);
            double f = 0.0, eps = nonMuThetaTol;
            int nn = nFree, mode = 1, niter = gPhi0RefEvalMax,
              nsim = gPhi0RefEvalMax, impr = 0, izs = 0; float rzs = 0;
            double dzs = 0; int idz = 0;
            gPhi0Self = this;
            gPhi0N1Bad = 0;
            gPhi0N1Evals = 0;
            if (n1qn1_ != NULL) {
              n1qn1_(gPhi0N1Cost, &nn, x.data(), &f, gg.data(), var.data(), &eps,
                     &mode, &niter, &nsim, &impr, zm.data(), &izs, &rzs, &dzs, &idz);
            }
            if (!gPhi0N1Bad && gPhi0N1Evals > 0) {
              for (int fi = 0; fi < nFree; ++fi) {
                int c = gPhi0FreeIx[(size_t)fi];
                double v = x[(size_t)fi];
                if (v < lo[c]) v = lo[c];
                if (v > hi[c]) v = hi[c];
                if (std::isfinite(v)) xmin[c] = v; else xmin[c] = par0[c];
              }
              for (int c = 0; c < nphi0; c++) {
                bool free = false;
                for (size_t q = 0; q < gPhi0FreeIx.size(); ++q)
                  if (gPhi0FreeIx[q] == c) { free = true; break; }
                if (!free) xmin[c] = par0[c];
              }
              n1qn1Done = true;
            }
            // n1qn1 got nowhere; fall through to the derivative-free search
            for (int c = 0; c < nphi0; c++) gPhi0Work[c] = par0[c];
          }
        }
        if (n1qn1Done) {
          // n1qn1 already wrote xmin; skip the derivative-free search entirely
        } else if (nonMuThetaOptType == 1) {
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
        for (int fi = 0; n1qn1Done ? false : (fi < nFree); fi++) {
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
      Rcpp::Function boundedOpt = nlmixr2[".saemBoundedResidOpt"];
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
    // the search moved phi0 underneath the established solve
    invalidateCompleteSolve();
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
    iniSubjectE(i, 1, ind, op, _rx,
                (_saemReadSlot == odeSlotThetaSens) ? rxThetaSens.update_inis
                                                    : rxPred.update_inis);
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



  // ONE-SOLVE gradient step for the DECLARED thetas.
  //
  // This is what etaDistLoglik should have been.  The observation likelihood is
  // the right objective, but handing it to refinePhi0Lik's SEARCH pays a full
  // population solve per candidate -- 25 of them per firing -- to rediscover a
  // derivative nlmixr2 already emits.
  //
  // Instead: the solve is FROZEN at the current parameters and taken once, the
  // theta-sensitivity peer gives d(f)/d(theta) for every declared theta out of
  // that same solve (linCmt promoted to linCmtB through odeSwap, so an analytic
  // model costs no integration at all), and nonMuGradPhi0() turns that into a
  // damped step.  The declared thetas reach the likelihood ONLY through
  // eta = Q(phiU(z); args(theta)), and the sensitivity model differentiates
  // through gammapInv/phiU exactly, so the chain rule is already in that one
  // derivative.  No re-solve, and the steps are only loosely coupled -- which
  // is fine, because SA damps each one anyway.
  //
  // Scheduled on its OWN cadence from iteration 0.  nonMuGradPhi0 is otherwise
  // called from inside refinePhi0Lik, so it inherits the SEARCH's schedule and
  // cannot run before nonMuThetaStart (half of nBurn+nEm by default) -- the
  // opposite of what src/nonMuThetaGrad.h prescribes: "run the cheap directed
  // step often, the expensive undirected one rarely".
  bool etaDistGradStep(unsigned int kiter, const vec &pas) {
    // The declared thetas' OBSERVATION-likelihood owner, which does not depend
    // on the family MLE being switched on.  `etaDistObsLik()` requires
    // `etaDistOn`, and that flag means "the FAMILY M-step is wanted" (see
    // R/saem.R) -- a different route entirely, on the prior side.  Requiring it
    // here made the gradient step unreachable under `etaDistMstep = FALSE`,
    // which is both the default and the setting the documentation recommends:
    // measured, `loglik=1 on=0 ndist=1 nphi0=4` and the step never ran, leaving
    // a declared theta on the cdf route with no owner that carries derivative
    // information.
    if (!etaDistLoglik || etaDistNdist <= 0) return false;
    // ...but ONLY for declarations whose thetas are in the observation path.
    // A declaration whose thetas are prior-only is owned by the Q2 step, and
    // letting this one at them too is two owners on one parameter -- the
    // failure this area keeps relearning.  Measured on the direct route, where
    // every declared theta is Q2: with this step also running, bWT went 0.7099
    // -> 0.5854 and prop.sd 0.1519 -> 0.2936 against truths of 0.75 and 0.15.
    if (!etaDistAnyQ1()) return false;
    if (!_saemThetaSensActive || nphi0 <= 0) return false;
    if (etaDistThetaPhi0.n_rows != (unsigned int)etaDistNdist) return false;
    // CADENCE.  This is the declared-distribution M-step, so it runs on that
    // step's own schedule (etaDistStart / etaDistEvery, default 20) rather
    // than every iteration.
    //
    // Every iteration is what the first implementation did, and it is not
    // affordable: each firing needs a COMPLETE-SYSTEM solve to produce the
    // sensitivity columns, which the MCMC's own solves do not carry.  At
    // etaDistEvery = 1 that is ~300 extra population solves and the fit ran
    // 7x the baseline without finishing.  The step is SA-damped anyway, so
    // firing it every iteration bought little even in principle.
    if (kiter < (unsigned int)etaDistStart) return false;
    if (((int)(kiter - (unsigned int)etaDistStart) % etaDistEvery) != 0) return false;
    if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
      RSprintf("[edGrad] it=%d gate passed (start=%d every=%d)\n",
               (int)kiter, etaDistStart, etaDistEvery);
    }
    if (_saemNonMuGradEvery > 1 &&
        ((int)kiter % _saemNonMuGradEvery) != 0) return false;
    // the free set is exactly the declared thetas, minus anything the user
    // fixed.  NOT every phi0 column: taking the rest is nonMuTheta="regress",
    // which has its own control.  The copula keeps its own closed form.
    std::vector<bool> isFix((size_t)nphi0, false);
    for (unsigned int j = 0; j < fixedIx0.n_elem; ++j) {
      if (fixedIx0(j) < (unsigned int)nphi0) isFix[(size_t)fixedIx0(j)] = true;
    }
    std::vector<int> saveFree = gPhi0FreeIx;
    gPhi0FreeIx.clear();
    for (int k = 0; k < etaDistNdist; ++k) {
      for (int t = 0; t < etaDistNth(k) && t < (int)etaDistThetaPhi0.n_cols; ++t) {
        int c = etaDistPhi0Col(k, t);
        if (c < 0 || c >= nphi0 || isFix[(size_t)c]) continue;
        bool dup = false;
        for (size_t q = 0; q < gPhi0FreeIx.size(); ++q)
          if (gPhi0FreeIx[q] == c) { dup = true; break; }
        if (!dup) gPhi0FreeIx.push_back(c);
      }
    }
    // The copula is stepped SEPARATELY below, not added here -- see there.
    if (gPhi0FreeIx.empty()) { gPhi0FreeIx = saveFree; return false; }
    std::vector<int> famIx = gPhi0FreeIx;
    // The copula gets its OWN step, over its own information matrix.
    //
    // Not because the observation likelihood does not identify it -- it does,
    // and rx__sens_rx_pred__BY_THETA_8___ is exactly that derivative -- but
    // because nonMuGradStep() gates on the CONDITIONING of the information
    // matrix and escalates Levenberg-Marquardt damping when rcond says it is
    // untrustworthy.  rxCor is a nearly flat direction, so putting it in the
    // same matrix as the family thetas made the whole 5x5 ill-conditioned and
    // the damping collapsed EVERY column: the family thetas went from
    // moved=0.694 in their own 4x4 step to moved=0.001 sharing one with it.
    // One Newton step couples its parameters through that matrix; two steps do
    // not.
    std::vector<int> corIx;
    for (int k = 0; k < etaDistNdist; ++k) {
      int cc = corCol(k);
      if (cc < 0 || isFix[(size_t)cc]) continue;
      bool dup = false;
      for (size_t q = 0; q < famIx.size(); ++q)
        if (famIx[q] == cc) { dup = true; break; }
      for (size_t q = 0; q < corIx.size() && !dup; ++q)
        if (corIx[q] == cc) dup = true;
      if (!dup) corIx.push_back(cc);
    }
    bool frz = _saemFreezeOde;
    _saemFreezeOde = false;
    bool moved = false;
    bool edTr = (getenv("NLMIXR2_ETADIST_OPT") != NULL);

    // PASS 1 -- the non-correlation thetas, by the damped Newton step.
    gPhi0FreeIx = famIx;
    if (!gPhi0FreeIx.empty()) {
      ensureCompleteSolve(kiter, true);
      bool m = nonMuGradPhi0(kiter, pas);
      if (edTr) RSprintf("[edGrad] it=%d family nFree=%d moved=%d\n",
                         (int)kiter, (int)gPhi0FreeIx.size(), (int)m);
      if (m) { moved = true; writeBackPhi0(gPhi0FreeIx); invalidateCompleteSolve(); }
    }

    // PASS 2 -- the correlation, with pass 1's values now FIXED.
    //
    // A ROOT FIND, not a Newton step.  With one free parameter the stationarity
    // condition is scalar -- score(rho) = 0 -- and that is a uniroot problem,
    // which is robust exactly where the Newton step is not: nonMuGradStep()
    // gates on the conditioning of the information matrix, and rxCor is a
    // nearly flat direction, so Levenberg-Marquardt damping eats the step and
    // returns "moved" without moving.  Measured: the family thetas stepped
    // moved=0.694 in their own 4x4, and 0.001 sharing a 5x5 with rxCor;
    // splitting into two Newton steps left it at 0.003, because a 1x1 Newton
    // step on a flat direction has the same problem the 5x5 did.
    //
    // nonMuGradPhi0()'s xEval mode returns the score at an arbitrary point,
    // which is the only thing a bracketing root find needs.  Bisection rather
    // than boost's toms748: each evaluation costs a population solve, so the
    // budget is tens of evaluations and the extra order of convergence buys
    // less than the guarantee of never leaving the bracket.
    gPhi0FreeIx = corIx;
    if (etaDistCorMethod != 2) gPhi0FreeIx.clear();   // a closed form owns it
    if (!gPhi0FreeIx.empty()) {
      ensureCompleteSolve(kiter, true);
      int c = gPhi0FreeIx[0];
      double cur = mprior_phi0(0, c);
      // Only ONE correlation may take this route.  With one free parameter the
      // problem is a scalar minimization; with two or more it is not, and
      // minimizing one coordinate at a time is not solving it.  Anything else
      // is left to the damped Newton step.
      bool oneCor = (corIx.size() == 1);
      bool rooted = false;
      double best = cur;
      if (oneCor) {
        std::vector<double> pv((size_t)nphi0);
        for (int q = 0; q < nphi0; ++q) pv[(size_t)q] = mprior_phi0(0, q);
        bool ok = false;
        // A LOCAL TRUST REGION, not the whole range.
        //
        // Maximizing the conditional observation likelihood over rho has no
        // interior optimum: collapsing the two latents onto one always fits the
        // CURRENT draws better, so the objective is monotone toward rho -> 1.
        // Searched over the full (-3, 3) this lands on the bound every time --
        // measured, argmin 2.9965, rho 0.995 -- which is the degeneracy
        // etaDistCorMstep's own documentation warns about.
        //
        // A trust region makes it a LOCAL problem, which is why the ordinary
        // bounded search keeps rho interior where an unbounded one does not.
        // Wider than the 0.75 phi0 uses, because the correlation's derivative
        // information is weak -- that is what ruled out both the Newton step
        // and the score root find -- so it needs room to move per firing rather
        // than sharper local curvature.  Derivative-free for the same reason.
        double r = (etaDistCorTrust > 0.0) ? etaDistCorTrust : 1.5;
        double lo2 = cur - r, hi2 = cur + r;
        if (lo2 < -3.0) lo2 = -3.0;
        if (hi2 > 3.0) hi2 = 3.0;
        double xm = brentMinPhi0Col(c, lo2, hi2, 1e-3, 40, pv, ok);
        if (ok && std::isfinite(xm)) { best = xm; rooted = true; }
        mprior_phi0.col(c).fill(cur);   // restore; the damped move is below
      }
      if (rooted) {
        double v = cur + pas(kiter) * (best - cur);
        if (std::isfinite(v)) {
          mprior_phi0.col(c).fill(v);
          moved = true; writeBackPhi0(gPhi0FreeIx); invalidateCompleteSolve();
        }
      }
      if (edTr) RSprintf("[edGrad] it=%d copula nCor=%d solved=%d argmin=%.4f cur=%.4f -> %.4f\n",
                         (int)kiter, (int)corIx.size(), (int)rooted, best, cur,
                         mprior_phi0(0, c));
    }
    _saemFreezeOde = frz;
    if (moved && (int)etaDistFiredK.size() == etaDistNdist) {
      for (int k = 0; k < etaDistNdist; ++k) etaDistFiredK[(size_t)k] = 1;
    }
    gPhi0FreeIx = saveFree;
    return moved;
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
        if (!saemNoThrow([&]{
              rowPred = phi1PredAt(i, ind, op, neqGuard, nH2Theta, eta0, bad); })) {
          rowBad[i] = 1; continue;
        }
      } else if (_saemPhi1UseAnalyticHess) {
        bool okH = false;
        if (!saemNoThrow([&]{
              okH = phi1AnalyticHessAt(i, ind, op, neqGuard, nH2Theta, eta0,
                                       rowPred, H); }) || !okH) {
          rowBad[i] = 1; continue;
        }
      } else {
        if (!saemNoThrow([&]{
              rowPred = phi1PredAt(i, ind, op, neqGuard, nH2Theta, eta0, bad); })) {
          rowBad[i] = 1; continue;
        }
        // restore of the row's solve/lhs state to eta0 is unnecessary --
        // calcEtaHessian's own FD fallback leaves the last perturbation's
        // solve behind too (its caller re-solves before any further read),
        // but SAEM's OWN model (odeSlotSaem) is what user_function reads
        // next, an independent peer/solve buffer.
        if (!bad && !saemNoThrow([&]{
              phi1FDHessAt(i, ind, op, neqGuard, nH2Theta, fdH, eta0,
                           rowPred, H, bad); })) {
          rowBad[i] = 1; continue;
        }
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

  // saemControl(zeroOmegaAnneal=).  saemix decays the variance of a parameter
  // WITHOUT IIV geometrically across the annealing phase
  // (diag.omega[i0.omega2] *= alpha0.sa every iteration, with
  // alpha0.sa = 10^(-3/nbiter.sa), R/main_mstep.R:91) rather than holding it:
  // wide exploration early, convergence late.  nlmixr2 pinned the placeholder
  // at zeroOmegaTune forever, which makes it a permanent noise floor -- no
  // single constant can be both, which is why widening zeroOmegaTune degrades
  // every parameter monotonically.
  //
  // Called ONLY from the fixed-restore branch, which has just reset these
  // diagonals to the placeholder, so the result is the absolute value
  // tune*coef^k rather than a compounding one.  It also deliberately runs
  // after the Gmin floor, which would otherwise fight the decay.
  void zeroOmegaAnneal(unsigned int kiter) {
    if (zeroOmegaAnnealCoef >= 1.0 || saemZeroOmegaPhi1.n_elem == 0) return;
    unsigned int kd = (kiter <= (unsigned int)nb_sa) ? kiter : (unsigned int)nb_sa;
    double fac = std::pow(zeroOmegaAnnealCoef, (double)kd);
    for (unsigned int f = 0; f < saemZeroOmegaPhi1.n_elem; ++f) {
      unsigned int c = saemZeroOmegaPhi1(f);
      if (c >= Gamma2_phi1.n_rows) continue;
      double v = Gamma2_phi1(c, c) * fac;
      if (v < 1e-12) v = 1e-12;   // keep IGamma2_phi1 invertible (saemix floors at double eps)
      Gamma2_phi1(c, c) = v;
    }
  }

  // Objective for zeroOmegaDirectStep(): the observation -log-likelihood with
  // the named phi1 columns overwritten by the candidate mu and every other
  // column left at its sampled value -- saemix's compute.Uy, which does
  // exactly this for i0.omega2 (R/func_aux.R:357-359).
  double zeroOmegaObjective(const double *p) const {
    mat phiCand = phiM;
    for (size_t j = 0; j < gZeroOmIx.size(); ++j) {
      phiCand.col(i1(gZeroOmIx[j])).fill(p[j]);
    }
    return computeUy(phiCand);
  }

  // saemix's ind.fix10 branch (R/main_mstep.R:63) for the phi1 columns the GLS
  // above cannot move: a mu-referenced random effect whose DECLARED variance
  // was zero.  Plambda1 is an Omega^-1-weighted normal equation, so such a
  // column's update collapses to "reproduce its sampled mean" -- and the
  // sampler cannot move it off its prior mean, which IS the current theta.
  // The M-step's fixed point is therefore the ini() value.
  //
  // NONMEM has the same two routes -- technical guide eqs. 1.45/1.46 through
  // mu, 1.47-1.52 by differentiating the entire joint density for a theta that
  // is not reachable through mu -- and demonstrably takes the second one here:
  // in ~/src/gamma_indpar/gamma_clv1_saem.ctl all seven thetas are MU_
  // referenced and five sit on $OMEGA (0.0 FIXED), so every one of those phi
  // is deterministic, yet the .ext still moves them (THETA5 -3.0 -> -2.42126).
  // A sampled-mean update returns the starting value by construction.
  //
  // So: maximize the observation likelihood in those coordinates directly,
  // then take the SAME damped stochastic-approximation step refinePhi1Lik()
  // takes, and keep MCOV1 consistent so the next iteration's COV1*MCOV1
  // reproduces it.  Intercept-only columns only -- the shape mu-referencing
  // produces, and the only one whose mu IS the theta.
  void zeroOmegaDirectStep(unsigned int kiter, const vec &pas) {
    if (nphi1 <= 0) return;
    gZeroOmIx.clear();
    for (unsigned int f = 0; f < saemZeroOmegaPhi1.n_elem; ++f) {
      unsigned int c = saemZeroOmegaPhi1(f);
      if (c >= (unsigned int)nphi1) continue;
      // Intercept-only columns only -- the shape mu-referencing produces, and
      // the only one whose mu IS the theta.  A covariate column's theta is a
      // regression coefficient and belongs to the GLS.
      uvec li = arma::find(LCOV1.col(c) == 1);
      if (li.n_elem != 1) continue;
      // fixedIx1 indexes LAMBDA rows, not phi1 columns.  For an intercept-only
      // column the two coincide, but map through LCOV1 rather than lean on
      // that: a model with a covariate on any OTHER column shifts the lambda
      // numbering and the coincidence stops holding.
      bool isFixed = false;
      for (unsigned int j = 0; j < fixedIx1.n_elem; ++j) {
        if (fixedIx1(j) == li(0)) { isFixed = true; break; }
      }
      if (isFixed) continue;
      gZeroOmIx.push_back((int)c);
    }
    if (gZeroOmIx.empty()) return;
    int nFree = (int)gZeroOmIx.size();

    gZeroOmSelf = this;
    gZeroOmLo.set_size(nFree);
    gZeroOmHi.set_size(nFree);
    std::vector<double> st((size_t)nFree), stp((size_t)nFree), xm((size_t)nFree);
    for (int j = 0; j < nFree; ++j) {
      double cur = mprior_phi1(0, gZeroOmIx[(size_t)j]);
      // ABSOLUTE local trust region around the current mu, for the same reason
      // refinePhi0Lik() uses one: the observation objective has NaN plateaus
      // far from the current value, and a relative radius lets a coordinate
      // that starts to drift grow its own step.
      gZeroOmLo(j) = cur - 0.75;
      gZeroOmHi(j) = cur + 0.75;
      st[(size_t)j] = cur;
      xm[(size_t)j] = cur;
      stp[(size_t)j] = 0.15;
    }
    gZeroOmBest.set_size(nFree);
    for (int j = 0; j < nFree; ++j) gZeroOmBest(j) = st[(size_t)j];
    gZeroOmEvalN = 0;
    gZeroOmBestF = 0.0;
    gZeroOmEvalMax = (nonMuThetaMaxEval > 0) ? nonMuThetaMaxEval : 20*nFree;
    int iconv, it, nfcall, iprint = 0;
    double ynewlo;
    // itmax generous; gZeroOmObj owns the real evaluation budget
    nelder_fn(gZeroOmNmFn, nFree, st.data(), stp.data(), 100*nFree,
              nonMuThetaTol, 1.0, 2.0, 0.5,
              &iconv, &it, &nfcall, &ynewlo, xm.data(), &iprint);

    for (int j = 0; j < nFree; ++j) {
      int c = gZeroOmIx[(size_t)j];
      double cur = mprior_phi1(0, c);
      double tgt = gZeroOmBest(j);
      if (!std::isfinite(tgt)) continue;
      mprior_phi1.col(c).fill(cur + pas(kiter) * (tgt - cur));
    }
    phi1BackSolveMCOV();
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
    Rcpp::Function boundedOpt = nlmixr2[".saemBoundedResidOpt"];
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
  void arWhiten(vec &e, vec &g) const {
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
    return arDYFinto(yt, ft, g, _scratch_ftAr, _scratch_gAr);
  }

  // The same computation with the AR(1) conditional prediction/SD written to
  // CALLER-OWNED buffers instead of the shared _scratch_ members.  const, so
  // the compiler enforces that a caller inside an OpenMP region (computeUy())
  // cannot corrupt the MCMC's own working state.
  vec arDYFinto(const vec &yt, const vec &ft, const vec &g,
                vec &ftAr, vec &gAr) const {
    vec e = yt - ft;
    vec gg = g;
    arWhiten(e, gg);
    ftAr = yt - e;
    gAr = gg;
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
    // one seed per observation: this chain's rows, endpoint b's span
    const uint64_t seedOff = _seedLayout.cens(kiter, mixIdx < 0 ? 0 : mixIdx, k) +
      (uint64_t)y_offset(b);
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
      nmSeqSeedSet(saemSeed, seedOff, i);
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

  int get_etaDistMapFail() { return _saemEtaDistMapFail; }
  // declared correlations the data do not identify, for the $runInfo warning
  int get_etaDistCorNotEst() { return etaDistCorNotEstimable(); }
  // The declared copula correlations as the sampler last had them.
  //
  // On the cdf route these also live in `rxCor.*` thetas and the fit rebuilds
  // them from there; on the DIRECT route there is no such theta, so without
  // this the value the sampler used is simply lost -- the fit reported a
  // correlation it never used, and .etaDistWarnCorFrozen() had to say so.
  // Returned alongside `etaDistCorWith` so R knows which pair each one joins.
  vec get_etaDistRho()     { return etaDistRho; }
  ivec get_etaDistCorWith(){ return etaDistCorWith; }
  mat get_mcmcAccTrace()   { return mcmcAccTrace; }
  mat get_mcmcStuckTrace() { return mcmcStuckTrace; }
  mat get_phiSdTrace()     { return phiSdTrace; }
  mat get_phiAcfTrace()    { return phiAcfTrace; }

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
    if (x.containsElementNamed("iaccept")) iaccept = as<double>(x["iaccept"]);
    if (!std::isfinite(iaccept) || iaccept < 0.0 || iaccept >= 1.0) iaccept = 0.0;
    if (x.containsElementNamed("iacceptSingle")) iacceptSingle = as<double>(x["iacceptSingle"]);
    if (x.containsElementNamed("iacceptPerId")) iacceptPerId = as<int>(x["iacceptPerId"]);
    if (x.containsElementNamed("rwOmega")) rwOmega = as<int>(x["rwOmega"]);
    if (x.containsElementNamed("nonMuThetaBhhh")) nonMuThetaBhhh = as<int>(x["nonMuThetaBhhh"]);
    if (x.containsElementNamed("etaDistLoglik")) etaDistLoglik = as<int>(x["etaDistLoglik"]);
    if (!std::isfinite(iacceptSingle) || iacceptSingle < 0.0 || iacceptSingle >= 1.0) iacceptSingle = 0.0;
    if (x.containsElementNamed("etaDistOn")) etaDistOn = as<int>(x["etaDistOn"]);
    if (x.containsElementNamed("etaDistCorOn")) etaDistCorOn = as<int>(x["etaDistCorOn"]);
    if (x.containsElementNamed("etaDistDebug")) etaDistDebug = as<int>(x["etaDistDebug"]);
    if (x.containsElementNamed("etaDistStart")) etaDistStart = as<int>(x["etaDistStart"]);
    if (x.containsElementNamed("etaDistEvery")) etaDistEvery = as<int>(x["etaDistEvery"]);
    if (x.containsElementNamed("etaDistCorTrust")) etaDistCorTrust = as<double>(x["etaDistCorTrust"]);
    if (x.containsElementNamed("etaDistCorMethod")) {
      etaDistCorMethod = as<int>(x["etaDistCorMethod"]);
    } else if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
      // Worth saying out loud.  A control that never arrives is invisible: two
      // estimators came back byte identical on six quantities because both ran
      // as method 0, and nothing in the output said so.
      RSprintf("[etaDist] etaDistCorMethod ABSENT from the control -- using %d\n",
               etaDistCorMethod);
    }
    // Argument expressions + their theta names, for the C++ native->theta map.
    etaDistExprs.clear(); etaDistExprThetas.clear();
    if (x.containsElementNamed("etaDistExprs") && !Rf_isNull(x["etaDistExprs"]) &&
        x.containsElementNamed("etaDistExprThetas") && !Rf_isNull(x["etaDistExprThetas"])) {
      List le = x["etaDistExprs"], lt = x["etaDistExprThetas"];
      for (int k = 0; k < le.size(); ++k) {
        CharacterVector ce = le[k], ct = lt[k];
        std::vector<std::string> e, t;
        for (int i = 0; i < ce.size(); ++i) e.push_back(as<std::string>(ce[i]));
        for (int i = 0; i < ct.size(); ++i) t.push_back(as<std::string>(ct[i]));
        etaDistExprs.push_back(e); etaDistExprThetas.push_back(t);
      }
    }
    if (etaDistEvery < 1) etaDistEvery = 1;
    if (x.containsElementNamed("etaDistSpreadGuard")) etaDistSpreadGuard = as<int>(x["etaDistSpreadGuard"]);
    if (x.containsElementNamed("etaDistSdLo")) etaDistSdLo = as<double>(x["etaDistSdLo"]);
    if (x.containsElementNamed("etaDistSdHi")) etaDistSdHi = as<double>(x["etaDistSdHi"]);
    if (x.containsElementNamed("etaDistSdTol")) etaDistSdTol = as<double>(x["etaDistSdTol"]);
    // AFTER every read above.  Printed before them it reported the C++
    // defaults no matter what the control carried, which is worse than no
    // trace at all: it says a control did not arrive when it did, and a valid
    // measurement gets thrown away on its word.
    if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
      RSprintf("[etaDist] etaDistCorMethod=%d etaDistCorTrust=%.3g "
               "spreadGuard=%d sdLo=%.3g sdHi=%.3g sdTol=%.3g etaDistEvery=%d "
               "etaDistStart=%d\n",
               etaDistCorMethod, etaDistCorTrust, etaDistSpreadGuard,
               etaDistSdLo, etaDistSdHi, etaDistSdTol, etaDistEvery,
               etaDistStart);
    }
    // per fit, not per session: the question this answers is "did THIS fit's
    // M-step run", so it cannot accumulate across fits the way
    // _saemPhi1RefineN does
    _saemEtaDistN = 0;
    _saemEtaDistOn = etaDistOn;
    _saemEtaDistObsLik = (etaDistLoglik && etaDistOn && etaDistNdist > 0) ? 1 : 0;
    if (etaDistNdist > 0) etaDistFiredK.assign((size_t)etaDistNdist, 0);
    etaDistCorFired = false;
    // Ingest the declaration metadata whenever R SENT it.  Not gated on either
    // M-step flag: those say which OWNER is wanted, and the declarations exist
    // regardless of who updates their parameters.
    //
    // It used to read `(etaDistOn || etaDistCorOn) && ...`, which made the whole
    // feature ride on `etaDistCorMstep` defaulting TRUE: on the direct route
    // `etaDistOn` is 0 by construction (the family MLE has nothing to do
    // there), so `saemControl(etaDistParam="direct", etaDistCorMstep=FALSE)`
    // ingested no latent, no family and no route, and the sampler fell back to
    // a standard normal without saying so.  Measured on the covariate arm:
    // bWT -0.2081 against a truth of 0.75 where the baseline is 0.7099, and
    // sampled etas at -0.4370 -- outside a gamma's support.
    if (x.containsElementNamed("etaDistLatent")) {
      etaDistLatent  = as<ivec>(x["etaDistLatent"]);
      etaDistFam     = as<ivec>(x["etaDistFam"]);
      etaDistDirect  = x.containsElementNamed("etaDistDirect") ?
        as<int>(x["etaDistDirect"]) : 0;
      etaDistQ2 = x.containsElementNamed("etaDistQ2") ?
        as<imat>(x["etaDistQ2"]) : imat();
      etaDistCorWith = as<ivec>(x["etaDistCorWith"]);
      if (x.containsElementNamed("etaDistUsable")) {
        etaDistUsable = as<ivec>(x["etaDistUsable"]);
      } else {
        etaDistUsable = ivec();
      }
      etaDistCov.clear(); etaDistCovNames.clear();
      if (x.containsElementNamed("etaDistCov")) {
        Rcpp::List cvl(x["etaDistCov"]);
        for (int k = 0; k < cvl.size(); ++k) {
          if (Rf_isNull(cvl[k])) {
            etaDistCov.push_back(arma::mat());
            etaDistCovNames.push_back(std::vector<std::string>());
            continue;
          }
          Rcpp::NumericMatrix m(cvl[k]);
          etaDistCov.push_back(as<arma::mat>(cvl[k]));
          std::vector<std::string> nms;
          Rcpp::List dn(m.attr("dimnames"));
          if (dn.size() == 2 && !Rf_isNull(dn[1])) {
            Rcpp::CharacterVector cn(dn[1]);
            for (int c = 0; c < cn.size(); ++c) nms.push_back(Rcpp::as<std::string>(cn[c]));
          }
          etaDistCovNames.push_back(nms);
        }
      }
      etaDistArgs    = as<mat>(x["etaDistArgs"]);
      // WHERE the model computes each argument: `rxEdA.<eta>.<role>`, one model
      // line per family argument, in the same column order as etaDistArgs.  The
      // compiled model evaluates these per observation -- covariates included,
      // through rxode2's ordinary covariate machinery, inside the ODE model
      // pool -- so the arguments are READ from the solve rather than evaluated a
      // second time here.
      //
      // R resolves the names to lhs indices, because that is where saem's model
      // is in hand.  Resolving them here through odeSwapLhsIndex() does not
      // work: saem drives its own solve (saem_lhs = rxInner.calc_lhs) and, at
      // the point the sampler first asks, no odeSwap slot is loaded at all --
      // every name came back -1.
      //
      // -1 means no anchor: the family emitted no line for that argument (a
      // normal-based family collapses onto its latent), and it keeps its
      // population value.
      etaDistAnchorNlhs = x.containsElementNamed("etaDistAnchorNlhs") ?
        as<int>(x["etaDistAnchorNlhs"]) : 0;
      etaDistAnchorIdx.clear();
      if (x.containsElementNamed("etaDistAnchorIdx")) {
        Rcpp::IntegerMatrix ai(x["etaDistAnchorIdx"]);
        int nr = ai.nrow(), nc = ai.ncol();
        etaDistAnchorIdx.resize((size_t)nr);
        for (int r = 0; r < nr; ++r) {
          etaDistAnchorIdx[(size_t)r].assign((size_t)nc, -1);
          for (int c = 0; c < nc; ++c) {
            int v = ai(r, c);
            etaDistAnchorIdx[(size_t)r][(size_t)c] =
              (v == NA_INTEGER || v < 0) ? -1 : v;
          }
        }
      }
      etaDistRho     = as<vec>(x["etaDistRho"]);
      etaDistThetaPhi0 = as<imat>(x["etaDistThetaPhi0"]);
      etaDistNth     = as<ivec>(x["etaDistNth"]);
      etaDistCorPhi0 = as<arma::ivec>(x["etaDistCorPhi0"]);
      if (x.containsElementNamed("etaDistMapFn")) etaDistMapR = x["etaDistMapFn"];
      etaDistNdist   = (int)etaDistLatent.n_elem;
      // sized HERE, not with the other resets above: etaDistNdist is assigned on
      // this line, so anything sized before it gets length zero -- and
      // etaDistCorSuffStat() then returns immediately on every iteration, which
      // reads exactly like the statistic being uninformative.
      if (etaDistNdist > 0) {
        etaDistCorSuff = arma::vec((unsigned int)etaDistNdist, arma::fill::zeros);
        etaDistCorEstim = arma::ivec((unsigned int)etaDistNdist);
        etaDistCorEstim.fill(-1);
        etaDistCorOffMag = arma::vec((unsigned int)etaDistNdist, arma::fill::zeros);
        etaDistCorOffMax = arma::vec((unsigned int)etaDistNdist, arma::fill::zeros);
        etaDistSdPrev = arma::vec((unsigned int)etaDistNdist);
        etaDistSdPrev.fill(NA_REAL);
        etaDistSdCur = arma::vec((unsigned int)etaDistNdist);
        etaDistSdCur.fill(NA_REAL);
      }
    }
    if (x.containsElementNamed("nu1B")) nu1B = as<int>(x["nu1B"]);
    if (nu1B < 0) nu1B = 0;
    if (x.containsElementNamed("nb1B")) nb1B = as<int>(x["nb1B"]);
    if (nb1B < 1) nb1B = 1;
    if (x.containsElementNamed("stepsizeRw")) stepsizeRw = as<double>(x["stepsizeRw"]);
    if (!std::isfinite(stepsizeRw) || stepsizeRw <= 0.0) stepsizeRw = 0.4;
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
    // A MIXTURE plus the direct parameterization is refused rather than run.
    //
    // The mixture E-step has its own kernel (do_mcmc_msaem), and it scores every
    // proposal with the Gaussian quadratic.  Teaching it the family prior is not
    // hard, but nothing exercises the combination and an untested prior in an
    // MCMC acceptance is exactly the kind of defect that shows up as slightly
    // wrong estimates rather than as an error.  So say so instead: on the cdf
    // route a mixture works as it always has, because there the sampled column
    // genuinely IS standard normal.
    if (nMix > 1 && etaDistDirect != 0 && etaDistNdist > 0) {
      Rf_error("a mixture model cannot use etaDistParam=\"direct\" yet\n"
               "  the mixture MCMC kernel scores its proposals with the Gaussian "
               "prior, which is not this model's prior\n"
               "  use etaDistParam=\"cdf\", where the sampled random effect really "
               "is standard normal");
    }
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
    // Gamma2_phi1 columns carrying no between-subject variability -- a
    // mu-referenced parameter whose omega was declared zero.  Excluded from
    // the omega sufficient statistics (it does not contribute to them) and
    // given a flat prior in the sampler, since it is not part of Omega.
    if (x.containsElementNamed("saemFlatPhi1")) {
      saemFlatPhi1 = as<uvec>(x["saemFlatPhi1"]);
    } else {
      saemFlatPhi1 = uvec();
    }
    if (x.containsElementNamed("saemZeroOmegaPhi1")) {
      saemZeroOmegaPhi1 = as<uvec>(x["saemZeroOmegaPhi1"]);
    } else {
      saemZeroOmegaPhi1 = uvec();
    }
    if (x.containsElementNamed("zeroOmegaAnnealCoef")) {
      zeroOmegaAnnealCoef = as<double>(x["zeroOmegaAnnealCoef"]);
    }
    if (!std::isfinite(zeroOmegaAnnealCoef) || zeroOmegaAnnealCoef <= 0.0) {
      zeroOmegaAnnealCoef = 1.0;
    }
    if (x.containsElementNamed("zeroOmegaDirect")) {
      zeroOmegaDirect = as<int>(x["zeroOmegaDirect"]);
    }
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
    omegaPoolMean = x.containsElementNamed("omegaPoolMean") ? as<int>(x["omegaPoolMean"]) : 0;
    _buildLambdaCol1 = true;
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
    nonMuThetaStart = x.containsElementNamed("nonMuThetaStart") ?
      as<int>(x["nonMuThetaStart"]) : -1;
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
      if (opt.containsElementNamed("saemPhi1EtaNonMu") &&
          !Rf_isNull(opt["saemPhi1EtaNonMu"])) {
        _saemPhi1EtaNonMu = as<ivec>(opt["saemPhi1EtaNonMu"]);
      } else {
        _saemPhi1EtaNonMu = arma::ivec(_saemPhi1H2EtaCol.n_elem, arma::fill::zeros);
      }
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
    _seedLayout.nu[0] = (uint64_t)nu(0);
    _seedLayout.nu[1] = (uint64_t)nu(1);
    _seedLayout.nu[2] = (uint64_t)nu(2);
    _seedLayout.nu1B = (uint64_t)nu1B;
    _seedLayout.nphi[0] = (uint64_t)nphi1;
    _seedLayout.nphi[1] = (uint64_t)nphi0;
    _seedLayout.nComp = (uint64_t)std::max(nMix, 1);
    _seedLayout.nM = (uint64_t)nM;
    _seedLayout.nmc = (uint64_t)nmc;
    _seedLayout.ntotal = (uint64_t)ntotal;
    // the draws own every iteration's seeds; the solves take the sequence after
    const uint64_t seedReserved = _seedLayout.iterBase(niter + nSaCov);
    if (seedReserved > 0xFFFFFFFFull) {
      Rcpp::warning("random draws exceed the seed range; some seeds repeat");
    }
    nmSeqSeedStart(saemSeed, seedReserved);
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
      //      shared-eta mixture has none).  Components start identical and
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
      // A flat column is not part of Omega, so it gets no prior: zeroing its
      // row and column of the inverse is the same as leaving it out of the
      // Cholesky (it cannot be correlated with anything -- a zero variance
      // carries no covariance -- so Omega is block diagonal about it and this
      // IS the smaller inverse embedded back).  The MCMC prior term and the
      // objective's quadratic form both lose it, which is what a parameter
      // with no between-subject variability should contribute.
      for (unsigned int _f = 0; _f < saemFlatPhi1.n_elem; ++_f) {
        unsigned int _c = saemFlatPhi1(_f);
        if (_c < IGamma2_phi1.n_rows) {
          IGamma2_phi1.row(_c).zeros();
          IGamma2_phi1.col(_c).zeros();
        }
      }
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
      set_mcmcphi(mphi1, i1, nphi1, Gamma2_phi1, IGamma2_phi1, mprior_phi1, rwScale1, rwScale1b);
      set_mcmcphi(mphi0, i0, nphi0, Gamma2_phi0, IGamma2_phi0, mprior_phi0, rwScale0, rwScale0b);
      mphi1.block = 0;
      mphi0.block = 1;

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
          do_mcmc_msaem(1, nu1, mx, mphi1, phiM, U_y, U_phi, (int)kiter, &rwScale1, &rwLam1);
          U_phi = etaDistSeedUphi(mphi1, phiM);
          do_mcmc_msaem(2, nu2, mx, mphi1, phiM, U_y, U_phi, (int)kiter, &rwScale1, &rwLam1);
          do_mcmc_msaem(3, nu3, mx, mphi1, phiM, U_y, U_phi, (int)kiter, &rwScale1, &rwLam1, &rwScale1b, &rwLam1b);
        }
        if (nphi0 > 0) {
          vec U_phi;
          do_mcmc_msaem(1, nu1, mx, mphi0, phiM, U_y, U_phi, (int)kiter, &rwScale0, &rwLam0);
          U_phi = etaDistSeedUphi(mphi0, phiM);
          do_mcmc_msaem(2, nu2, mx, mphi0, phiM, U_y, U_phi, (int)kiter, &rwScale0, &rwLam0);
          do_mcmc_msaem(3, nu3, mx, mphi0, phiM, U_y, U_phi, (int)kiter, &rwScale0, &rwLam0, &rwScale0b, &rwLam0b);
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
            do_mcmc(1, nu1, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale1, &rwLam1);
            U_phi = etaDistSeedUphi(mphi1, cur_phiM);
            do_mcmc(2, nu2, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale1, &rwLam1);
            do_mcmc(3, nu3, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale1, &rwLam1, &rwScale1b, &rwLam1b);
          }
          if (nphi0 > 0) {
            vec U_phi;
            do_mcmc(1, nu1, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale0, &rwLam0);
            U_phi = etaDistSeedUphi(mphi0, cur_phiM);
            do_mcmc(2, nu2, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale0, &rwLam0);
            do_mcmc(3, nu3, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit, (int)kiter, jMix + 1, &rwScale0, &rwLam0, &rwScale0b, &rwLam0b);
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
          do_mcmc(1, nu1, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale1, &rwLam1);
          U_phi    = etaDistSeedUphi(mphi1, phiM);
          // NONMEM runs mode 1B directly after mode 1, once each subject's
          // conditional moments have had time to accumulate
          if (buildMode1B(i1, kiter)) {
            do_mcmc(4, nu1B, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale1, &rwLam1);
          }
          do_mcmc(2, nu2, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale1, &rwLam1);
          do_mcmc(3, nu3, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale1, &rwLam1, &rwScale1b, &rwLam1b);
        }
        if(nphi0>0) {
          vec U_phi;
          do_mcmc(1, nu1, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale0, &rwLam0);
          U_phi    = etaDistSeedUphi(mphi0, phiM);
          do_mcmc(2, nu2, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale0, &rwLam0);
          do_mcmc(3, nu3, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit, (int)kiter, 0, &rwScale0, &rwLam0, &rwScale0b, &rwLam0b);
        }
        if (DEBUG>0) Rcout << "mcmc successful\n";
        // Close this iteration's mixing diagnostics.  Placed here, after every
        // kernel and before the sufficient statistics are accumulated, so the
        // traces describe exactly the draws the M-step is about to use.
        mcmcCloseIter(kiter, (unsigned int)niter, phiM);
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
      // A flat column has had its row and column of IGamma2_phi1 zeroed (it is
      // not part of Omega, so it gets no prior), and that zero propagates into
      // CGamma21 -- whose inverse then fails outright.  It also should: the
      // Omega^-1-weighted normal equation is the wrong update for a column with
      // no Omega.  With a flat prior the GLS reduces to ordinary least squares
      // on that column alone -- it cannot be correlated with any other, a zero
      // variance carrying no covariance -- and for the intercept-only design
      // that mu-referencing produces, OLS is just the mean of the sampled phi.
      //
      // So the flat rows are decoupled from the system before the solve and
      // given that mean afterwards.  This is the same "theta += mean(eta)"
      // update imp gets from impMuInterceptStep(), written in saem's
      // parameterization, where the theta IS the column's location.
      mat CG21 = CGamma21;
      std::vector<unsigned int> flatRow;
      std::vector<unsigned int> flatCol;
      for (unsigned int _f = 0; _f < saemFlatPhi1.n_elem; ++_f) {
        unsigned int _c = saemFlatPhi1(_f);
        if (_c >= (unsigned int)LCOV1.n_cols) continue;
        uvec _li = find(LCOV1.col(_c) == 1);
        if (_li.n_elem != 1) continue;   // covariate design: leave it alone
        unsigned int _l = _li(0);
        if (_l >= CG21.n_rows) continue;
        CG21.row(_l).zeros();
        CG21.col(_l).zeros();
        CG21(_l, _l) = 1.0;
        flatRow.push_back(_l);
        flatCol.push_back(_c);
      }
      Plambda1=inv_sympd(CG21)*sum((D1Gamma21%(COV1.t()*statphi11)),1);
      for (size_t _k = 0; _k < flatRow.size(); ++_k) {
        Plambda1(flatRow[_k]) = arma::mean(statphi11.col(flatCol[_k]));
      }
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
      // collapsed IOV: the group's columns share one mean, so average their
      // solutions before anything downstream reads them
      if (omegaPoolMean) {
        if (_buildLambdaCol1) {
          lambdaCol1.set_size(LCOV1.n_rows);
          lambdaCol1.fill((unsigned int)nphi1);
          for (unsigned int l = 0; l < LCOV1.n_rows; ++l) {
            for (unsigned int c = 0; c < LCOV1.n_cols; ++c) {
              if (LCOV1(l, c) == 1) { lambdaCol1(l) = c; break; }
            }
          }
          _buildLambdaCol1 = false;
        }
        poolLambdaGroups(Plambda1);
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
        // A theta the declared-distribution M-step owns must not ALSO be moved
        // by this GLS.  The two answer to different objectives and pull in
        // opposite directions: measured on Bauer's gamma model, the GLS alone
        // drives lclm from 6.686 up to 7.51 while the M-step's own target is
        // 4.5, and the spread guard means the M-step only fires on a minority
        // of iterations -- so the GLS wins the rest and the family's variance
        // absorbs the conflict (rvCL 51.6 against a truth of 0.086).
        //
        // Holding them here is what makes the M-step the SOLE owner, which is
        // the whole point of taking them out of refinePhi0Lik as well: on an
        // iteration the guard blocks, the right behavior is for these thetas
        // to stay put, not to drift.
        //
        // Conditional on the M-step having ACTUALLY FIRED for that family, for
        // the same reason refinePhi0Lik's hold-out is.  An M-step that never
        // fires must not leave its parameters unowned: the GLS is then the only
        // thing that would move them, and holding them out of it returns the
        // ini() values as if they were estimates.
        //
        // That is not hypothetical.  On Bauer's g1 arm the spread guard
        // rejected every one of 400 iterations, and lclm/lv1m/lclrv/lv1rv and
        // rxCor all came back at exactly their ini values while tq/tv2/prop.sd
        // -- the thetas this block does not hold -- moved normally.  The
        // etaDistFiredK guard added to refinePhi0Lik did not catch it: that
        // refinement is gated on (distribution == 4 || nonMuThetaRegress), and
        // Bauer's model has a prop() endpoint, so it never runs there at all.
        if (etaDistOn && etaDistNdist > 0) {
          for (int k = 0; k < etaDistNdist; ++k) {
            // In the observation-likelihood mode refinePhi0Lik owns these, so
            // the GLS holds them out unconditionally -- there is no "has the
            // M-step fired yet" question, because the owner is not the family
            // M-step.
            if (!etaDistObsLik() &&
                (int)etaDistFiredK.size() == etaDistNdist &&
                etaDistFiredK[(size_t)k] == 0) continue;   // GLS: see below
            for (int t = 0; t < etaDistNth(k); ++t) {
              int c = etaDistPhi0Col(k, t);
              if (c < 0 || c >= nphi0) continue;
              uvec li = arma::find(LCOV0.col(c) == 1);
              for (unsigned int q = 0; q < li.n_elem; ++q) {
                uvec hit = arma::find(jcov0 == (li(q) + c*(unsigned int)LCOV0.n_rows));
                if (hit.n_elem == 1) Plambda0(hit(0)) = MCOV0(jcov0(hit(0)));
              }
            }
          }
          if ((etaDistObsLik() && etaDistCorMethod == 2) ||
              (etaDistCorOn && etaDistCorFired)) {
            for (int kc = 0; kc < etaDistNdist; ++kc) {
              int cc = corCol(kc);
              if (cc < 0) continue;
              uvec li = arma::find(LCOV0.col(cc) == 1);
              for (unsigned int q = 0; q < li.n_elem; ++q) {
                uvec hit = arma::find(jcov0 == (li(q) + cc*(unsigned int)LCOV0.n_rows));
                if (hit.n_elem == 1) Plambda0(hit(0)) = MCOV0(jcov0(hit(0)));
              }
            }
          }
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
      // Only nonMuThetaRegress suppresses the stochastic phi0 update.  The
      // etaDist-only route takes just the DECLARED columns, and the GLS still
      // owns the rest -- suppressing it there would freeze prop.sd and every
      // other non-mu theta at whatever the last GLS pass left.
      bool skipStochPhi0 = nonMuThetaRegress && (distribution != 4) &&
        (kiter >= (unsigned int)niter_phi0);
      if (!skipStochPhi0) {
        mprior_phi0=COV0*MCOV0;
      }
      if (_saemPhi1PoolReady && kiter >= (unsigned int)niter_phi0 &&
          (kiter - (unsigned int)niter_phi0) % (unsigned int)phi1ThetaEvery == 0) {
        refinePhi1Lik(kiter, pas);
      }
      // saemix gates its ind.fix10 branch to kiter >= nbiter.sa
      // (R/main_mstep.R:57): through the annealing phase everything goes via
      // the GLS, and only once the variances have settled does the direct
      // maximization take over.  nonMuThetaEvery thins it further -- each step
      // is damped by pas(kiter) anyway, so refining every iteration buys
      // little against a whole extra optimization (nonMuThetaMaxEval solves of
      // the full population) per iteration.
      // ODE-free distribution M-step.  The step itself reads only the sampled
      // etas -- no solve -- which is why it originally ran every iteration.
      // That reasoning is incomplete: the step is cheap, but the parameters it
      // moves to are not.  Measured on Bauer's gamma model, 80 iterations, the
      // copula M-step alone (etaDistCorMstep, which defaults on) took the fit
      // from 184.9s to 700.8s -- 3.8x -- while improving both CL (3.722 ->
      // 5.963, truth 5.03) and the correlation (0.344 -> 0.384, truth 0.438).
      // The cost is not the M-step arithmetic; it is that the parameters it
      // reaches make gammapInv's iterative inversion work much harder.
      //
      // So it gets a cadence, exactly like the non-mu theta refinement's
      // nonMuThetaEvery: each step is damped by pas(kiter) anyway, so running
      // it every k-th iteration keeps most of the benefit for a fraction of the
      // cost.
      // The cadence exists because the FAMILY MLE is expensive -- the comment
      // above measures 184.9s -> 700.8s -- and because each step is damped by
      // pas(kiter) anyway.  Neither applies to the eta-density objective: it
      // does no ODE work at all, which is the same reason etaDistCorSuffStat()
      // runs every iteration.  Starving it is not free: at etaDistEvery=20 it
      // fires ~15 times across 300 iterations and lclrv crawled from -0.2231 to
      // -0.3349 against a truth of -0.6931 (MARE 19.26%), where the MLE route
      // reached -0.6169 (3.71%).
      //
      // So a model whose declared thetas are prior-only runs every iteration;
      // everything else keeps the cadence it had.
      const bool edQ2Every = etaDistAnyQ2();
      if ((etaDistOn || etaDistCorOn || edQ2Every) && etaDistNdist > 0 && nphi0 > 0 &&
          kiter >= (unsigned int)etaDistStart &&
          (edQ2Every ||
           ((int)(kiter - (unsigned int)etaDistStart) % etaDistEvery) == 0)) {
        if (etaDistMstep(kiter, pas)) {
          _saemEtaDistN++;
          // Map the updated NATIVE parameters back onto the user's thetas.
          // This is the ONE place the loop touches R -- once per iteration, not
          // per objective evaluation -- because a declared family's arguments
          // are arbitrary expressions over thetas that C++ cannot evaluate.
          if (getenv("NLMIXR2_ETADIST_OPT") != NULL && kiter % 20 == 0)
        Rprintf("etaDist map: C++=%d R=%d\n", _saemEtaDistCppMap, _saemEtaDistRMap);
      if (etaDistOn && !etaDistMapR.isNULL()) {
            Rcpp::Function mapFn(etaDistMapR);
            for (int k = 0; k < etaDistNdist; ++k) {
              int na = rxEtaDistNarg(etaDistFam(k));
              if (na <= 0) continue;
              NumericVector av(na);
              for (int i = 0; i < na; ++i) av[i] = etaDistArgs(k, i);
              // C++ first.  This is the map that used to be the ONE place this
              // loop touched R; it now only falls back when an argument
              // expression is outside the C++ grammar, which the parse reports
              // rather than guesses at.
              int nth = etaDistNth(k);
              if (k < (int)etaDistExprs.size() &&
                  (int)etaDistExprThetas[(size_t)k].size() == nth && nth > 0) {
                std::vector<double> st((size_t)nth), got2((size_t)nth);
                for (int t = 0; t < nth; ++t) {
                  int c = etaDistPhi0Col(k, t);
                  st[(size_t)t] = (c >= 0 && c < nphi0) ? mprior_phi0(0, c) : 0.0;
                }
                if (rxEtaDistArgsToThetas(etaDistExprs[(size_t)k],
                                          etaDistExprThetas[(size_t)k],
                                          st.data(), av.begin(), got2.data())) {
                if (getenv("NLMIXR2_ETADIST_OPT") != NULL) _saemEtaDistCppMap++;
                  for (int t = 0; t < nth; ++t) {
                    int c = etaDistPhi0Col(k, t);
                    if (c < 0 || c >= nphi0 || !std::isfinite(got2[(size_t)t])) continue;
                    mprior_phi0.col(c).fill(got2[(size_t)t]);
                  }
                  continue;
                }
              }
              if (getenv("NLMIXR2_ETADIST_OPT") != NULL) _saemEtaDistRMap++;
              RObject got = mapFn(k + 1, av);
              if (got.isNULL()) { _saemEtaDistMapFail++; continue; }
              NumericVector th(got);
              if ((int)th.size() != nth) { _saemEtaDistMapFail++; continue; }
              for (int t = 0; t < nth; ++t) {
                int c = etaDistPhi0Col(k, t);
                if (c < 0 || c >= nphi0 || !std::isfinite(th[t])) continue;
                mprior_phi0.col(c).fill(th[t]);
              }
            }
          }
          // Each correlation to its OWN column.  This used to AVERAGE every
          // rho into a single column, which blends a +0.6 pair and a -0.4 pair
          // into 0.1; it was unreachable only because the R side disabled the
          // whole path whenever there was more than one.
          for (int k = 0; k < etaDistNdist; ++k) {
            int cc = corCol(k);
            if (cc < 0 || etaDistCorWith(k) < 0) continue;
            // the copula theta is atanh(rho): the expansion writes the
            // correlation as tanh() of it
            double r = etaDistRho(k);
            if (!std::isfinite(r)) continue;
            double a = std::atanh(std::max(std::min(r, 0.999), -0.999));
            if (std::isfinite(a)) mprior_phi0.col(cc).fill(a);
          }
          // keep MCOV0 consistent so the next COV0*MCOV0 reproduces this
          for (int c = 0; c < nphi0; c++) {
            uvec li = arma::find(LCOV0.col(c) == 1);
            if (li.n_elem == 0) continue;
            mat Xc = COV0.cols(li);
            vec bc;
            if (arma::solve(bc, Xc.t() * Xc, Xc.t() * mprior_phi0.col(c))) {
              for (unsigned int j = 0; j < li.n_elem; ++j) MCOV0(li(j), c) = bc(j);
            }
          }
        }
      }
      if (zeroOmegaDirect && saemZeroOmegaPhi1.n_elem > 0 &&
          kiter >= (unsigned int)nb_sa &&
          (kiter - (unsigned int)nb_sa) % (unsigned int)nonMuThetaEvery == 0) {
        zeroOmegaDirectStep(kiter, pas);
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
      // The start point is its OWN quantity, not niter_phi0.  niter_phi0 is half
      // of nBurn+nEm, so on a default fit this refinement -- nlmixr2's analogue
      // of the NONMEM technical guide's non-mu theta route (eqs 1.47-1.52) --
      // is barred until the whole burn-in is over.  Measured on Bauer's gamma
      // model: the phi0 GLS parks CL at 7.51 within 24 iterations, this cannot
      // run until 200, and by then pas(kiter) ~ 1/k is far too small to travel
      // back to the 4.79 NONMEM's SAEM reaches.  The thetas are frozen for 376
      // of 400 iterations, which is why a 12x evaluation budget moved CL by
      // 0.004: the search never gets a chance to matter.
      unsigned int phi0Start = (nonMuThetaStart >= 0) ?
        (unsigned int)nonMuThetaStart : (unsigned int)niter_phi0;
      // The copula's sufficient statistic, EVERY iteration.  Second moments
      // over phiM -- no solve, so no reason to put it on a cadence -- and it is
      // what the parameter adjustment below reads.
      // the prior needs the parsed covariate expressions as much as the M-step
      // does; parsing is arithmetic, so once per iteration costs nothing
      etaDistBuildCovRpn();
      etaDistCorSuffStat(kiter, pas);
      // The declared thetas' own step: ONE solve, exact gradient, one damped
      // move -- from iteration 0, on nonMuThetaGradEvery, independent of the
      // search's nonMuThetaStart.  Placed BEFORE the search so that when both
      // run the search sees the gradient's result, never the other way round.
      etaDistGradStep(kiter, pas);
      if ((distribution == 4 || phi0ObsLikRoute()) &&
          nphi0 > 0 && kiter >= phi0Start &&
          (kiter - phi0Start) % (unsigned int)nonMuThetaEvery == 0) {
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
      // The flat columns contribute nothing to the sufficient statistics that
      // estimate Omega, so they are held at their placeholder rather than
      // taking a value out of G1.  (The placeholder only keeps the matrix
      // invertible; IGamma2_phi1 has already had the column zeroed, so it
      // reaches neither the sampler nor the objective.)
      for (unsigned int _f = 0; _f < saemFlatPhi1.n_elem; ++_f) {
        unsigned int _c = saemFlatPhi1(_f);
        if (_c < Gamma2_phi1.n_rows) {
          Gamma2_phi1.row(_c).zeros();
          Gamma2_phi1.col(_c).zeros();
          Gamma2_phi1(_c, _c) = 1.0;
        }
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
        // Gamma2_phi1Report is what the fit REPORTS, and it was snapshotted
        // above -- before this restore -- so a fix()ed variance came back as
        // the M-step's unconstrained estimate: fix(0.3) reported as 0.318
        // while the sampler correctly used 0.3 (#1073).  Only the fixed cells
        // are touched; every other reported value stays exactly as it was.
        Gamma2_phi1Report.elem(Gamma2_phi1fixedIx) = Gamma2_phi1fixedValues(Gamma2_phi1fixedIx);
        zeroOmegaAnneal(kiter);
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
            _saemResidGen++;
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
      // Full per-subject second moment for mode 1B (see xpost_phi).  Same
      // Robbins-Monro weighting as mpost_phi/cpost_phi, so the three stay
      // consistent with one another.
      if (nu1B > 0) {
        if (xpost_phi.n_slices != (unsigned int)N) {
          xpost_phi.zeros(nphi, nphi, N);
        }
        for (int _i = 0; _i < N; ++_i) {
          mat S(nphi, nphi, fill::zeros);
          for (int _k = 0; _k < nmc; ++_k) {
            vec v = phi.slice(_k).row(_i).t();
            S += v * v.t();
          }
          S /= (double)nmc;
          xpost_phi.slice(_i) += pash(kiter) * (S - xpost_phi.slice(_i));
        }
      }
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
    // Report the proportion at the score-zero point: one exact M-step at the final
    // responsibilities, so sum_i (r_i - p) == 0 holds for every mixProbMethod.  The
    // annealed step size and the Dirichlet-style regularization above stabilize the
    // TRAJECTORY; leaving their shrinkage in the reported value makes the proportion
    // disagree with the fit's own per-subject probabilities (#1058).
    if (nMix > 1 && arma::accu(mixWeights) > 0.0) mixProb = mean(mixWeights, 0).t();
    phiFile.close();
  }


private:

  user_funct user_fn;

  uvec nu;
  int niter;
  int saemSeed = 99;
  saemSeedLayout _seedLayout;
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
  // saemControl(iaccept=): TARGET Metropolis acceptance rate for the
  // random-walk kernels (0 = do not adapt, the historical behaviour).
  // saemix calls this proba.mcmc (default 0.4); NONMEM calls it IACCEPT, and
  // Bauer's gamma control streams set it to 0.3.
  double iaccept = 0.0;
  // saemControl(iacceptSingle=): the SAME target for the coordinate-wise
  // kernel.  The optimal acceptance rate depends on the proposal's dimension:
  // ~0.234 for a multidimensional symmetric random walk (kernel 2, all
  // coordinates at once) and ~0.44 for a one-at-a-time
  // Metropolis-within-Gibbs update (kernel 3).  saemix uses ONE proba.mcmc
  // (0.4) for both, which is about right for kernel 3 and too HIGH for
  // kernel 2 -- and targeting too high an acceptance rate forces the step too
  // small, which is the direction that under-disperses the latent normals.
  double iacceptSingle = 0.0;
  // saemControl(stepsizeRw=): Robbins-Monro rate for that adaptation
  // (saemix stepsize.rw, default 0.4).
  double stepsizeRw = 0.4;
  // Per-column multiplier on the random-walk step, adapted toward `iaccept`.
  // rmcmc is only the STARTING value of this -- it is saemix's rw.init, which
  // nlmixr2 froze because the adaptation was never ported.
  vec rwScale1, rwScale0;
  // Kernel 3 (Metropolis-within-Gibbs) adapts toward iacceptSingle, kernel 2
  // toward iaccept.  Separate vectors: sharing one is what froze the chain.
  vec rwScale1b, rwScale0b;
  // saemControl(iacceptPerId=): one random-walk scale per SUBJECT, the way
  // NONMEM tunes its lambda, rather than one per coordinate pooled over the
  // population.  Sized to N on first use; left empty (and so inert) when off.
  vec rwLam1, rwLam0;      // per-subject scale for kernel 2 (target iaccept)
  vec rwLam1b, rwLam0b;    // ... and for kernel 3 (target iacceptSingle)

  // ---- MCMC mixing diagnostics -------------------------------------------
  //
  // saem computed its acceptance rate every iteration, used it to adapt the
  // random-walk scale, and threw it away.  Nothing about the chain reached the
  // fit, so a chain that had stopped moving was indistinguishable from one
  // exploring properly -- and "the chain is not mixing" is the diagnosis for
  // the declared-distribution M-step collapsing, so it had to become
  // measurable before it could be fixed.
  //
  // Four traces, one row per iteration:
  //   mcmcAccTrace   pooled acceptance rate, one column per kernel (1,2,3,1B)
  //   mcmcStuckTrace fraction of SUBJECTS that accepted NOTHING that iteration
  //                  -- the number a pooled rate cannot show, and the one that
  //                  says whether a healthy-looking 0.3 is everyone at 0.3 or
  //                  half the population never moving
  //   phiSdTrace     pooled SD of each phi column across subjects x chains.
  //                  For a declared distribution's latent this is 1 BY
  //                  CONSTRUCTION, so any departure is mixing, not signal.
  //   phiAcfTrace    lag-1 autocorrelation of each phi column against the
  //                  PREVIOUS iteration's draws.  This is the direct measure:
  //                  1.0 means the chain did not move at all.
  mat mcmcAccTrace, mcmcStuckTrace, phiSdTrace, phiAcfTrace;
  vec mcmcAccNum, mcmcAccDen;      // per-kernel, accumulated within an iteration
  vec mcmcAccById;                 // per-subject accepts within an iteration
  double mcmcAccByIdTrials = 0.0;
  // Per-subject acceptances, counted SEPARATELY per random-walk kernel.
  // Kernel 2 and kernel 3 are different proposals with different optimal
  // acceptance rates (iaccept ~0.234 multivariate, iacceptSingle ~0.44
  // one-at-a-time), so pooling them into one rate and comparing it against one
  // target drives the scale toward a compromise that is right for neither --
  // the same defect that froze the pooled path.
  vec mcmcAccByIdK2, mcmcAccByIdK3;
  double mcmcAccByIdK2Trials = 0.0, mcmcAccByIdK3Trials = 0.0;
  mat phiMprevIter;                // last iteration's phiM, for the lag-1 acf
  int iacceptPerId = 0;
  // saemControl(rwOmega=): propose the mode-2 random walk from lambda*Omega
  // (NONMEM eq. 1.139) rather than from a diagonal.
  int rwOmega = 0;
  // saemControl(nonMuThetaBhhh=): update the non-mu thetas with ONE per-subject
  // BHHH Newton step subject to an alpha acceptance test (NONMEM eqs.
  // 1.47-1.52 and the text after 1.46), instead of damping the argmax of a
  // full derivative-free maximization.  A full maximization of a nearly flat
  // or degenerate direction lands on the boundary, and the SA damping then
  // only sets how fast the theta marches there; a step that has to prove it
  // improved the objective cannot.
  int nonMuThetaBhhh = 0;
  // saemControl(etaDistLoglik=): use the general log-likelihood objective for
  // the declared-distribution M-step instead of fitting native parameters and
  // inverting them.  Opt-in while it is measured against the inversion.
  int etaDistLoglik = 0;

  // saemControl(etaDistLoglik=TRUE): estimate the DECLARED thetas from the
  // OBSERVATION likelihood rather than by fitting the family to the sampled
  // etas.
  //
  // This is the M-step the construction actually implies.  The complete data
  // is (y, z) with z the latent standard normal, so
  //
  //   log p(y, z | theta) = log p(y | z, theta) + log p(z)
  //
  // and log p(z) is theta-free.  eta = Q(phiU(z); args(theta)) is a
  // DETERMINISTIC transform, not observed data, so the family density never
  // enters the Q-function at all -- the declared parameters are structural
  // parameters of the mean function and belong to the observation likelihood
  // like any other non-mu theta.
  //
  // Fitting the family to the eta sample -- what both the MLE route and the
  // peer-density route do -- is a heuristic.  It measurably helps (the gamma
  // arms), but it is not this, and it is what needs the spread guard: an
  // over-dispersed latent makes the eta sample look over-dispersed and the
  // family widens to cover it.  The observation likelihood has no such
  // failure mode, because a theta that makes the PREDICTIONS worse is
  // rejected whatever the latents look like.
  //
  // A covariate on a distribution parameter then needs no machinery at all:
  // the model itself recomputes eta per record from the candidate thetas.
  bool etaDistObsLik() const {
    return etaDistLoglik && etaDistOn && etaDistNdist > 0;
  }

  // "phi0 is estimated by direct optimization of the observation likelihood".
  //
  // ONE predicate, because nonMuThetaRegress gates FOUR things -- whether
  // refinePhi0Lik runs, whether the stochastic phi0 update is skipped, whether
  // the ODE may be frozen during the search, and whether the search gets a
  // local trust region -- and a second entry added to only some of them would
  // leave the mode running that optimizer without the trust region a normal
  // model needs.
  //
  // In practice the second entry changes nothing: nonMuTheta DEFAULTS to
  // "regress", so nonMuThetaRegress is already 1 and this refinement already
  // runs on every saem fit.  What is actually restricted is WHEN
  // (nonMuThetaStart defaults to half of nBurn+nEm) and HOW OFTEN
  // (nonMuThetaEvery), and that the search is derivative-free unless
  // nonMuThetaOpt="n1qn1" and nonMuThetaGrad=TRUE.  So etaDistLoglik is a
  // SCHEDULING/ownership change to machinery that already runs -- not a new
  // estimation method -- and any comparison of it has to be against
  // nonMuThetaStart/nonMuThetaEvery, or the two are conflated.
  bool phi0ObsLikRoute() const {
    // etaDistObsLik() is deliberately NOT here any more.  It used to widen the
    // SEARCH's gate, which is the expensive way to use the observation
    // likelihood -- a population solve per candidate.  The declared thetas now
    // get etaDistGradStep() instead: one solve, exact gradient, one damped
    // step.  nonMuThetaRegress still governs the search for everything else.
    return nonMuThetaRegress;
  }
  // How often the Shi-difference fallback supplied the non-mu gradient because
  // the analytic sensitivity path's bad-solve ladder was exhausted.  Reported,
  // not hidden: a fit that spends most of its refinements on a finite
  // difference is a fit whose sensitivity peer is failing, and that is worth
  // knowing rather than inferring from the runtime.
  int _saemShiFallbackN = 0;
  // Times the native-parameters -> thetas map could not produce thetas, so the
  // family's fitted parameters were computed and then thrown away.  Counted
  // because the failure is otherwise invisible: the fit converges, looks
  // entirely normal, and the M-step simply had no effect.  The commonest cause
  // is an argument expression the map cannot invert -- a COVARIATE on a
  // distribution parameter, which the model-block dist() form now allows and
  // this inversion cannot represent (it solves for one POPULATION-level set of
  // native parameters, and a covariate gives every subject their own).
  int _saemEtaDistMapFail = 0;
  // Whether the M-step has ever actually moved family k's parameters.
  // refinePhi0Lik() hands those thetas over to the M-step and stops optimizing
  // them; if the M-step then never fires -- the spread guard can legitimately
  // hold it back for a whole fit -- NOBODY moves them and they are silently
  // returned as the ini() values.  Ownership is therefore conditional on the
  // step having actually fired for that family.
  std::vector<int> etaDistFiredK;
  // Same idea for the copula: its hold-out from the GLS is only legitimate
  // once its closed form has actually written a value.
  bool etaDistCorFired = false;
  // NONMEM's proposal kernel "mode 1B" (technical guide, "The MCMC method of
  // Expectation in SAEM"): after the first few iterations, propose from a
  // Gaussian built out of each subject's OWN accumulated conditional mean and
  // variance, rather than from the population prior (mode 1) or a random walk
  // around the current point (modes 2/3).  Bauer describes it as "a type of
  // importance sampling kernel for SAEM".
  //
  // It is an INDEPENDENCE sampler, so unlike modes 2/3 its acceptance ratio
  // carries a proposal-density correction; and unlike mode 1 the proposal is
  // not the prior, so that correction does not cancel.  See do_mcmc case 4.
  //
  // nu1B = 0 disables it (the default), which is exactly the historical
  // behaviour.  nb1B is the iteration it starts at -- the moments have to
  // accumulate first (NONMEM starts after the 10th).
  int nu1B = 0;
  int nb1B = 10;
  // Per-(subject x chain) proposal mean and SD, rebuilt each iteration from
  // mpost_phi / cpost_phi.  Empty when mode 1B is off or not yet started.
  mat m1bMean, m1bSd;
  // ---- declared-distribution M-step (see rxEtaDistMle) --------------------
  // saemControl(etaDistMstep=): estimate a declared family's parameters by a
  // distribution fit to the sampled etas rather than through the data
  // likelihood.  0 = off (the historical behaviour).
  int etaDistOn = 0;
  // saemControl(etaDistCorMstep=): update a declared Gaussian copula's
  // correlation from its CLOSED FORM -- the sample correlation of the latent
  // pair -- rather than leaving it to the general theta refinement.
  //
  // Independent of etaDistOn: the correlation needs this even when the family
  // parameters are estimated the ordinary way.  rxEtaDistExpand() gives the
  // rxCor.* theta an atanh scale with lower=-Inf/upper=Inf (rxode2
  // R/etaDist.R), so tanh() of it can reach 1 and refinePhi0Lik -- which
  // maximizes the observation likelihood CONDITIONAL on the current draws, with
  // no prior term to penalize degeneracy -- walks there: at rho=1 the copula
  // partner's latent collapses onto its partner's and two random effects become
  // one.  Measured on Bauer's gamma data, rho pinned at 1.000 in 3 of 7 fits
  // across seeds and start points, contributing 128% of one of the eight MARE
  // terms on its own.
  //
  // The closed form cannot do that: it is a correlation coefficient, bounded by
  // construction.  Damped like every other M-step here, and the damped value is
  // what the refinement then warm-starts from.
  int etaDistCorOn = 0;
  // getOption("nlmixr2.etaDistDebug"): 0 off; 1 traces the M-step (the first two
  // iterations and every tenth thereafter) while it acts; 2 traces the same but
  // does NOT apply the update, so the trajectory shown is an ordinary fit's,
  // unperturbed.  Level 2 is what established that the early wide latent spread
  // is a mixing transient rather than a broken construction.
  int etaDistDebug = 0;
  // First iteration the M-step is allowed to run.  It must NOT run from
  // iteration 0: the whole signal it reads is the pooled latents' departure
  // from N(0,1), so an unmixed chain is indistinguishable from a badly wrong
  // family, and it acts on the difference.  Measured on Bauer's gamma model
  // (etaDistDebug=2, which traces without acting), the pooled latent spread
  // starts at 2.15 and relaxes through 1.58 (it=10) and 1.25 (it=20) to 1.04
  // (it=30); acting on the it=0 draws drove the fitted shape from 7.39 to 1.74
  // and then to 0.51, each shrink widening the mapped etas and feeding the
  // next -- a runaway, with the copula correlation pinned at its clamp on the
  // way.  saemix gates the analogous ind.fix10 step to kiter >= nbiter.sa for
  // the same reason (R/main_mstep.R:57).
  int etaDistStart = 0;
  // How often the declared-distribution M-step runs.  1 = every iteration.
  int etaDistEvery = 1;
  // saemControl(etaDistCorTrust=): half-width of the local trust region the
  // copula's bounded 1-D search works in, on the atanh scale.  Wider than the
  // 0.75 the phi0 search uses -- see the search itself for why.
  double etaDistCorTrust = 1.5;
  // saemControl(etaDistCor=): HOW the copula correlation is updated, always
  // AFTER the distributional thetas have moved.
  //   0 "observed"  the product-moment correlation of the latent pair -- what
  //                 this has always done.  Sensitive to the latent SPREAD: an
  //                 over-dispersed latent drives it to its clamp.
  //   1 "analytic"  the Gaussian-copula identity rho = 2*sin(pi*rho_S/6) from
  //                 the RANKS.  Exact for the copula, and rank-based, so the
  //                 spread cannot reach it.
  //   3 "posterior" the EM M-step for a unit-diagonal correlation: second
  //                 moments about ZERO rather than about the sample mean, since
  //                 the latent's PRIOR mean is zero and a shifted sample is
  //                 information rather than nuisance.
  //   2 "optimize"  a bounded 1-D search of the observation objective inside a
  //                 local trust region (etaDistCorTrust).  Only for a model
  //                 with ONE correlation; a system is not a scalar problem.
  int etaDistCorMethod = 0;
  // Sufficient statistic for the copula correlation, refreshed EVERY iteration.
  // R enters the complete-data likelihood only through S_n = sum(n n'), and
  // n = L z, so S_n = L S_z L' -- three scalars per pair.  It is second moments
  // over phiM, so it costs nothing next to a solve and there is no reason to
  // put it on the M-step's cadence.
  arma::vec etaDistCorSuff;     // the normalized suff-stat rho, per family
  arma::ivec etaDistCorEstim;   // 1 estimable, 0 not, -1 no partner
  arma::vec etaDistCorOffMag;   // |S_z off-diagonal|/N, for the report
  // the LARGEST |S_z off-diagonal| seen over the whole fit.  Estimability is a
  // property of the trajectory, not of the endpoint: S_z -> I is what CONVERGENCE
  // looks like (S_n = L S_z L' -> L L' = R), so the instantaneous off-diagonal
  // being small says rho is settled, not that it was never identified.  A
  // correlation the data cannot inform is one whose off-diagonal NEVER left the
  // noise floor.
  arma::vec etaDistCorOffMax;
  // Per declared distribution: its argument expressions and the theta names
  // they are written over, so the map back onto thetas stays in C++.
  std::vector<std::vector<std::string> > etaDistExprs;
  std::vector<std::vector<std::string> > etaDistExprThetas;
  // Acceptable pooled spread for the latent normals; outside this the draws are
  // not yet a sample from anything worth fitting.
  //
  // The upper bound is 1 on principle, not by tuning: the latent's PRIOR is
  // exactly N(0,1), and a posterior is not wider than its prior, so a pooled
  // spread above 1 means the chain has not settled -- it is never evidence that
  // the family is too narrow.  Measured on Bauer's gamma model the spread
  // starts at 2.15, is still 1.43 at iteration 10 and 1.25 at 20, and only
  // settles (1.04 at 30, 0.93 at 50, ~0.85 thereafter) once the chain has
  // mixed.  An earlier version of this bound was 1.5, which admitted exactly
  // those transients: acting on the it=10 draws took the gamma's shape from
  // 7.389 to 3.696, and ten iterations later it was 0.0885 with the mapped etas
  // averaging 176 -- each widening feeding the next.
  // A pure DIVERGENCE cap, not a calibration.  The real test is
  // etaDistSdTol below; these only exclude a latent that has run away or
  // collapsed outright, so they are deliberately loose.
  double etaDistSdLo = 0.2, etaDistSdHi = 5.0;
  // saemControl(etaDistSdTol=): the M-step waits until the pooled latent spread
  // STOPS CHANGING between attempts, rather than until it reaches some level.
  // Level cannot do this job.  A settled chain under a wrong family sits at
  // sd 1.40 on Bauer's g1 and that spread IS the information the step needs; a
  // still-burning chain passes through 1.40 on its way down from 2.6 and acting
  // there is a runaway.  The two are indistinguishable by value and obvious by
  // trajectory.  <= 0 disables the test and falls back to the cap alone.
  double etaDistSdTol = 0.10;
  // spread at each family's previous M-step attempt, and this attempt's,
  // for that comparison.  Prev advances only at the end of the M-step.
  arma::vec etaDistSdPrev, etaDistSdCur;
  // saemControl(etaDistSpreadGuard=FALSE): skip the spread test entirely.
  // Setting the bounds wide is NOT the same thing -- the lower bound still
  // applies, and a sample whose spread has collapsed is rejected just as a
  // widened upper bound stops rejecting an inflated one.  Off means off.
  int etaDistSpreadGuard = 1;
  int etaDistNdist = 0;          // number of declared random effects
  ivec etaDistLatent;            // phi column of each one's OWN latent normal
  ivec etaDistFam;               // family code (rxEtaDistQ/rxEtaDistLogD)
  // saemControl(etaDistParam="direct"): the declared eta IS the random effect
  // and carries its family as its prior, rather than being a standard normal
  // latent that a decoder turns into one.  The two routes fit the SAME model by
  // different parameterizations, so they are comparable -- which is the point
  // of having both -- but the MCMC has to score them differently: on the cdf
  // route the sampled phi column really is N(0,1) and the Gaussian quadratic is
  // exact, while here it is (say) Gamma(0.5, 0.5) and the quadratic is simply a
  // different prior.  Reading one as the other is what "fits the wrong model
  // silently" means, so the flag comes from what the EXPANSION built, not from
  // what the control asked for (see .etaDistIsDirect()).
  int etaDistDirect = 0;
  // Per declaration x theta slot: 1 where the PRIOR identifies that theta and
  // the observation likelihood does not.  NoLimits.jl's Q1/Q2 partition
  // (_partition_q1_q2_names, src/estimation/common.jl:4557), computed on the R
  // side by .etaDistThetaSplit() from the model text.
  //
  // All zero on the cdf route: there the declared thetas sit inside the decoder
  // and DO reach the prediction, so the observation likelihood identifies them
  // and the existing steps own them.  Non-zero only where an `rxEdA.*` anchor
  // is read by nothing, which is the direct route's signature.
  imat etaDistQ2;
  // Is theta t of declaration k in the Q2 set?  Bounds-checked for the same
  // reason etaDistPhi0Col() is: this package builds with -DNDEBUG, so armadillo
  // does not check, and an empty matrix reads out of bounds and returns garbage
  // that happens to pass a 0/1 test.
  bool etaDistIsQ2(int k, int t) const {
    if (etaDistNdist <= 0) return false;
    if ((int)etaDistQ2.n_rows != etaDistNdist) return false;
    if (k < 0 || k >= (int)etaDistQ2.n_rows) return false;
    if (t < 0 || t >= (int)etaDistQ2.n_cols) return false;
    return etaDistQ2(k, t) != 0;
  }
  // Are ALL of declaration k's thetas prior-only?
  //
  // All, not any.  A declaration mixing prior-only and observation-path thetas
  // is the NON-SEPARABLE case: the complete-data likelihood does not split into
  // a piece each owner can maximize alone, so it stays with the
  // observation-likelihood steps entirely.  NoLimits does the same thing, from
  // the same reasoning -- it empties its whole Q2 set when `extra_objective`
  // couples the two (saem.jl:3204-3211).
  //
  // It also keeps the optimization vector whole: Q2 can vary every theta of the
  // declaration rather than a subset with the rest pinned.
  bool etaDistAllQ2(int k) const {
    int n = etaDistNth(k);
    if (n <= 0) return false;
    for (int t = 0; t < n; ++t) if (!etaDistIsQ2(k, t)) return false;
    return true;
  }
  // Does any declaration qualify?  Cheap gate so a model with none pays nothing.
  // Is any declaration's theta set in the OBSERVATION path?  The complement of
  // etaDistAllQ2() over declarations: cdf puts every declared theta there (the
  // decoder reads the anchors), direct puts none.
  bool etaDistAnyQ1() const {
    for (int k = 0; k < etaDistNdist; ++k) if (!etaDistAllQ2(k)) return true;
    return false;
  }

  bool etaDistAnyQ2() const {
    for (int k = 0; k < etaDistNdist; ++k) if (etaDistAllQ2(k)) return true;
    return false;
  }
  // Per-declaration covariate values, one row per SUBJECT in saem's own subject
  // order (phiM row r is subject r % N -- see the mixture weighting at
  // phiM_weighted).  Empty for a declaration with no covariate, which is the
  // nSym == 0 case and behaves exactly as before.
  std::vector<arma::mat> etaDistCov;
  // and their names, in the same order as the matrix columns -- the parser
  // needs the symbol, the objective needs the value.
  std::vector< std::vector<std::string> > etaDistCovNames;
  // The lhs index of each declaration argument's `rxEdA.<eta>.<role>` line,
  // resolved in R at setup against saem's own model.  -1 = no anchor.
  std::vector< std::vector<int> > etaDistAnchorIdx;
  // Harvested anchor values, one row per SUBJECT and one column per (k, arg).
  // This is the subject's FIRST record, which is the whole story when the
  // declaration's covariates do not vary within a subject.
  mutable arma::mat etaDistAnchorVal;
  mutable std::vector<char> etaDistAnchorHave;
  // EVERY record's anchors, per subject: a flat nrec x (ndist*na) block.  Kept
  // because a declaration whose covariates vary WITHIN a subject has no single
  // argument set, and its prior is a weighted per-observation likelihood rather
  // than one density at one record.
  mutable std::vector< std::vector<double> > etaDistAnchorRows;
  mutable std::vector<int> etaDistAnchorNrec;
  // Per declaration: do this declaration's anchors vary within any subject?
  // 0 = not yet computed for this pass, 1 = varies, 2 = constant.
  mutable std::vector<char> etaDistVariesK;
  // 1 where this M-step owns the declaration's family; 0 where it stands down
  // for that declaration alone -- an argument that varies by subject has no
  // single population value to fit.  Empty means every declaration is usable,
  // so metadata from an older R side behaves as it always did.
  ivec etaDistUsable;
  ivec etaDistCorWith;           // declared eta it is copula-correlated with, or -1
  mat etaDistArgs;               // current NATIVE parameters, ndist x maxNarg
  vec etaDistRho;                // current copula correlation per declared eta
  imat etaDistThetaPhi0;         // phi0 COLUMN of each declared theta
  // The ONLY way to read etaDistThetaPhi0.  It is a 0x0 matrix whenever the R
  // side could not match a declared theta to a phi0 column (saem_fit.R leaves
  // the `matrix(-1L, 0, 0)` it starts with), which is a real state and not an
  // error -- but this package builds with -DNDEBUG, so armadillo's bounds
  // checking is OFF and `etaDistThetaPhi0(k, t)` on an empty matrix reads out
  // of bounds and returns garbage.  That garbage then passed the
  // `0 <= c < nphi0` test at the call sites and a family argument was written
  // into whichever theta the garbage named -- on a 3-theta model, a 1-in-3
  // chance of landing on a theta the declaration has no claim on, showing up
  // only as an unrelated parameter having moved.
  int etaDistPhi0Col(int k, int t) const {
    if (etaDistNdist <= 0) return -1;
    if ((int)etaDistThetaPhi0.n_rows != etaDistNdist) return -1;
    if (k < 0 || k >= (int)etaDistThetaPhi0.n_rows) return -1;
    if (t < 0 || t >= (int)etaDistThetaPhi0.n_cols) return -1;
    return etaDistThetaPhi0(k, t);
  }
  ivec etaDistNth;               // how many thetas each declared eta has
  // phi0 column of each declared family's copula correlation theta, -1 where
  // that family has no partner.  A VECTOR, not a scalar: a model with two
  // correlated pairs has two correlations, and collapsing them to one column
  // either averaged them (blending +0.6 and -0.4 into 0.1) or -- as it actually
  // behaved -- disabled the machinery outright and left the second one with no
  // owner at all.
  arma::ivec etaDistCorPhi0;
  // the phi0 column family k's correlation lives in, or -1
  int corCol(int k) const {
    if (k < 0 || k >= (int)etaDistCorPhi0.n_elem) return -1;
    int c = etaDistCorPhi0(k);
    return (c >= 0 && c < nphi0) ? c : -1;
  }
  bool anyCorCol() const {
    for (int k = 0; k < (int)etaDistCorPhi0.n_elem; ++k)
      if (corCol(k) >= 0) return true;
    return false;
  }
  RObject etaDistMapR;           // R closure: (k, args) -> thetas (once/iter)
  // Running per-subject SECOND MOMENT E[phi phi'] (nphi x nphi x N).  cpost_phi
  // only keeps the elementwise E[phi^2], i.e. the diagonal; NONMEM's mode 1B
  // proposal density uses the full individual conditional variance --- its
  // B_i (technical guide eq. 1.147) is an OUTER product, so it carries the
  // off-diagonals.  For a model whose individual posterior is correlated
  // across coordinates (a copula/inverse-CDF model certainly is) a diagonal
  // proposal is the materially weaker version.  Only accumulated when mode 1B
  // is actually on, so an ordinary fit pays nothing.
  cube xpost_phi;
  // Lower Cholesky of each subject's conditional covariance, restricted to the
  // sampled columns; empty when the covariance was not usable this iteration.
  field<mat> m1bChol;
  bool m1bFull = false;
  vec pas, pash;
  vec minv;
  int nmc;
  int nM;

  int ntotal, N, mlen;
  vec y, ys;    //ys is y sorted by endpnt
  mat evt;
  mat phiM;
  uvec indio;
  uvec saemFlatPhi1;
  // Mu-referenced phi1 columns whose DECLARED omega was zero and which
  // .preProcessZeroOmegaMuRef() rewrote to the zeroOmegaTune placeholder.
  // Distinct from saemFlatPhi1 (a genuine zero): these carry a real, fixed
  // variance whose only job is to let the sampler move the column.
  uvec saemZeroOmegaPhi1;
  // saemControl(zeroOmegaAnneal=): per-iteration multiplier applied to that
  // placeholder over the simulated-annealing phase.  saemix decays the
  // variance of a parameter WITHOUT IIV by alpha0.sa = 10^(-3/nbiter.sa) every
  // SA iteration (R/main_mstep.R:91), a 1000x shrink across the phase, rather
  // than holding it -- exploration early, convergence late.  1.0 = hold (the
  // historical behaviour).
  double zeroOmegaAnnealCoef = 1.0;
  // saemControl(zeroOmegaDirect=): update those columns' thetas by directly
  // maximizing the observation likelihood (computeUy) instead of by the
  // Omega^-1-weighted GLS, which cannot move them.  saemix's ind.fix10 branch;
  // NONMEM technical guide eqs. 1.47-1.52.
  int zeroOmegaDirect = 0;
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
int nonMuThetaStart = -1;  // first iteration refinePhi0Lik may run; -1 = niter_phi0
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
  // Collapsed (Panhard & Samson) IOV: the pooled columns also share ONE mean --
  // the single theta the user declared.  Because the group's Gamma block is
  // compound-symmetric, 1 is an eigenvector of it, so Gamma^-1 1 is proportional
  // to 1 and the exact constrained GLS for that shared mean is the equal-weight
  // average of the group's unconstrained solutions.  0 leaves the means alone
  // (the two-level form pins them at 0 instead).
  int omegaPoolMean = 0;
  // lambda -> phi1 column, from LCOV1's single 1 per row; empty when unused
  uvec lambdaCol1;
  bool _buildLambdaCol1 = false;

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
		   const mat mprior_phi1,
		   vec &rwScale,
		   vec &rwScale3) {
    mphi1.i = i1;
    mphi1.nphi = nphi1;
    mphi1.Gamma_phi=chol(Gamma2_phi1);
    mphi1.IGamma2_phi = IGamma2_phi1;
    mphi1.Gdiag_phi.zeros(nphi1, nphi1);
    // rmcmc is the INITIAL scale (saemix rw.init); rwScale1 carries the
    // acceptance-rate adaptation from iteration to iteration.  With
    // iaccept == 0 rwScale1 stays all-ones and this is bit-identical to
    // before.
    if (rwScale.n_elem != (unsigned int)nphi1) rwScale.ones(nphi1);
    if (rwScale3.n_elem != (unsigned int)nphi1) rwScale3.ones(nphi1);
    if (iacceptPerId) {
      if ((int)rwLam1.n_elem != N) rwLam1.ones(N);
      if ((int)rwLam0.n_elem != N) rwLam0.ones(N);
      if ((int)rwLam1b.n_elem != N) rwLam1b.ones(N);
      if ((int)rwLam0b.n_elem != N) rwLam0b.ones(N);
    }
    // BASE step only -- rmcmc * the prior SD.  The acceptance-rate scale is
    // applied per KERNEL in do_mcmc, not folded in here.
    //
    // Folding it in here is what made kernels 2 and 3 share one rwScale vector
    // while adapting it toward DIFFERENT targets (0.234 over all coordinates
    // for the multivariate block walk, 0.44 one coordinate at a time for the
    // Metropolis-within-Gibbs sweep).  They fought, the scale ran to its clamp,
    // a wild proposal produced a non-finite objective, and the chain latched
    // (see the U_y rescue in do_mcmc).  Measured on theo_sd: with both
    // adaptations on, acceptance was exactly 0.000 on all three kernels and
    // every subject was frozen, with the objective 60 units worse than
    // origin/main, which has no acceptance adaptation at all.  Either
    // adaptation ALONE was fine.
    mphi1.Gdiag_phi.diag() = sqrt(Gamma2_phi1.diag())*rmcmc;
    // Omega-shaped step (saemControl(rwOmega=)).  chol() already succeeded
    // above for Gamma_phi, so this cannot throw where that did not.  No
    // per-coordinate rwScale: an Omega-shaped proposal has ONE scale, which is
    // what NONMEM's lambda is, and per-subject when iacceptPerId is on.
    mphi1.Gfull_phi = mphi1.Gamma_phi * rmcmc;
    mphi1.mprior_phiM = repmat(mprior_phi1,nmc,1);
  }

  // Build mode 1B's per-subject proposal moments from the accumulated
  // conditional mean (mpost_phi) and second moment (cpost_phi):
  //   var_ic = E[phi^2] - E[phi]^2
  // replicated across the nmc chains so it lines up with phiM's row layout.
  // Returns false when the moments are not usable yet (too early, or a
  // degenerate/negative variance), in which case the kernel is skipped for
  // this iteration rather than proposing from a broken density.
  bool buildMode1B(const uvec &cols, unsigned int kiter) {
    if (nu1B <= 0 || (int)kiter < nb1B) return false;
    unsigned int nc = cols.n_elem;
    if (nc == 0) return false;
    mat m = mpost_phi.cols(cols);              // N x nc
    m1bMean = repmat(m, nmc, 1);
    if (!m1bMean.is_finite()) return false;

    // FULL individual conditional covariance (NONMEM's B_i, eq. 1.147):
    //   V_i = E[phi phi'] - E[phi] E[phi]'
    // Symmetrized, with a jitter escalation so a not-yet-moved or numerically
    // indefinite V_i still factors instead of aborting the kernel.
    m1bFull = false;
    if (xpost_phi.n_slices == (unsigned int)N) {
      field<mat> L(N);
      bool ok = true;
      for (int i = 0; i < N && ok; ++i) {
        mat Vi = xpost_phi.slice(i).submat(cols, cols);
        vec mi = m.row(i).t();
        Vi -= mi * mi.t();
        Vi = 0.5 * (Vi + Vi.t());
        double sc = Vi.diag().max();
        if (!std::isfinite(sc) || sc <= 0.0) sc = 1.0;
        mat Li;
        bool got = false;
        double jit = 0.0;
        for (int t = 0; t < 8 && !got; ++t) {
          mat Vj = Vi;
          if (jit > 0.0) Vj.diag() += jit * sc;
          if (arma::chol(Li, Vj, "lower") && Li.is_finite()) got = true;
          else jit = (jit == 0.0) ? 1e-8 : jit * 100.0;
        }
        if (!got) ok = false; else L(i) = Li;
      }
      if (ok) { m1bChol = L; m1bFull = true; }
    }
    if (m1bFull) return true;

    // Fallback: diagonal-only proposal from cpost_phi.  Correct, just the
    // weaker kernel -- it cannot exploit correlation between coordinates.
    mat v = cpost_phi.cols(cols) - (m % m);
    for (unsigned int c = 0; c < nc; ++c) {
      for (unsigned int r = 0; r < v.n_rows; ++r) {
        if (!std::isfinite(v(r, c)) || v(r, c) < 1e-8) v(r, c) = 1e-8;
      }
    }
    m1bSd = sqrt(repmat(v, nmc, 1));
    return m1bSd.is_finite();
  }

  // Correlated mode 1B noise: L_i z per (subject x chain) row.
  mat mode1BNoise(const mat &noise) const {
    mat out(noise.n_rows, noise.n_cols);
    for (unsigned int r = 0; r < noise.n_rows; ++r) {
      unsigned int subj = r % (unsigned int)N;
      out.row(r) = (m1bChol(subj) * noise.row(r).t()).t();
    }
    return out;
  }

  // Q(x) = 0.5 (x-m)' V^-1 (x-m), the mode 1B proposal's exponent.  The
  // -0.5*log|V| normalizer is identical at the current and proposed points and
  // cancels from the Metropolis-Hastings ratio, so it is deliberately omitted.
  vec mode1BQ(const mat &x) const {
    mat d = x - m1bMean;
    vec q(d.n_rows);
    if (m1bFull) {
      for (unsigned int r = 0; r < d.n_rows; ++r) {
        unsigned int subj = r % (unsigned int)N;
        vec z = arma::solve(arma::trimatl(m1bChol(subj)), d.row(r).t());
        q(r) = 0.5 * arma::dot(z, z);
      }
    } else {
      mat vv = m1bSd % m1bSd;
      q = 0.5 * sum((d % d) / vv, 1);
    }
    return q;
  }

  // The ODE-free distribution M-step.  Runs EVERY iteration, like the residual
  // step and for the same reason: it reads only the sampled etas, so it costs
  // nothing beyond one simplex over a handful of parameters.
  //
  // Returns true when anything moved, in which case the caller maps the new
  // NATIVE parameters back onto the user's thetas.
  bool etaDistMstep(unsigned int kiter, const vec &pas) {
    // Q2 counts as a reason to be here.  This function hosts THREE owners --
    // the family MLE (etaDistOn), the copula closed form (etaDistCorOn) and the
    // prior-only Q2 step -- and bailing on the first two left Q2, which needs
    // neither, unreachable.  Measured on the direct route, where every declared
    // theta is Q2-owned and etaDistOn is 0 by construction:
    // `etaDistCorMstep=FALSE` returned bWT 0.1382 against a truth of 0.75
    // where the default path gives 0.7099, because nothing owned the thetas.
    if ((!etaDistOn && !etaDistCorOn && !etaDistAnyQ2()) ||
        etaDistNdist <= 0) return false;
    if (etaDistArgs.n_rows != (unsigned int)etaDistNdist) return false;
    // Per ATTEMPT, not per fit.  Left standing from the previous attempt, an
    // entry for a family neither loop visits this time would be copied into the
    // baseline below as though it had just been measured.
    if (etaDistSdCur.n_elem == (unsigned int)etaDistNdist) etaDistSdCur.fill(NA_REAL);
    bool moved = false;
    // In the observation-likelihood mode the declared thetas are estimated by
    // refinePhi0Lik against the observation likelihood, so this step has
    // nothing to do for the FAMILIES.  The copula closed form still runs: the
    // correlation is a property of the latent block, not of the mean function,
    // and no observation-likelihood term identifies it.
    bool famOff = etaDistObsLik();
    // latent normals actually seen by each declared eta's quantile: its own
    // sampled column, or -- for a copula member -- the correlated combination
    // the model forms (rho*z_j + sqrt(1-rho^2)*z_k).
    std::vector< std::vector<double> > w((size_t)etaDistNdist);
    for (int k = 0; k < etaDistNdist; ++k) {
      int ck = etaDistLatent(k);
      if (ck < 0 || ck >= (int)phiM.n_cols) return false;
      int j = etaDistCorWith(k);
      w[(size_t)k].resize(phiM.n_rows);
      // On the DIRECT route phiM ALREADY holds the eta, and there is no latent
      // to combine: the correlation lives in the prior's copula term
      // (rxEtaDistPairLogD), not in a linear combination of normals.  Forming
      // rho*z_j + sqrt(1-rho^2)*z_k out of two gamma draws would hand the MLE a
      // mixture of two subjects' clearances and call it one subject's.  The
      // MARGINAL is what the family MLE wants, and a copula's marginal is the
      // marginal.
      const bool edEtaScale = etaDistDirectOn();
      for (unsigned int r = 0; r < phiM.n_rows; ++r) {
        double zk = phiM(r, ck);
        if (j < 0 || edEtaScale) { w[(size_t)k][r] = zk; continue; }
        int cj = etaDistLatent(j);
        if (cj < 0 || cj >= (int)phiM.n_cols) return false;
        double rho = etaDistRho(k);
        if (!std::isfinite(rho)) rho = 0.0;
        double s2 = 1.0 - rho*rho;
        w[(size_t)k][r] = rho*phiM(r, cj) + (s2 > 0 ? std::sqrt(s2) : 0.0)*zk;
      }
    }
    // each declared family: latent -> eta via the CURRENT parameters, then MLE
    // NOT gated on famOff.  The general objective below lives in this loop and
    // is enabled by exactly the flag that sets famOff, so gating the loop on it
    // made that block unreachable -- etaDistLoglik=TRUE turned the family MLE
    // off and then skipped the thing meant to replace it, leaving every
    // declared theta at its ini() value while the step still counted as having
    // run.  famOff now guards only the MLE-and-inversion half, which is what
    // the comment above it always described.
    // `etaDistOn` is the etaDistMstep CONTROL (R/saem.R:589), so this loop used
    // to vanish entirely when the family M-step was switched off -- taking the
    // eta-density objective below with it, even though that objective is not
    // the family M-step and is the ONLY thing that identifies a prior-only
    // theta.  Measured: on the direct route with etaDistMstep=FALSE, lclm came
    // back 1.3901 from a start of 1.3863, never having moved, while the rest of
    // the fit converged normally.
    const bool edQ2 = etaDistAnyQ2();
    etaDistBuildCovRpn();
    for (int k = 0; (etaDistOn || edQ2) && k < etaDistNdist; ++k) {
      // Per DECLARATION, and it guards the MLE-and-inversion half ONLY.
      //
      // It used to `continue` here, which also skipped the general objective
      // below -- the one route that CAN fit a covariate-carrying declaration,
      // because it maximizes over the thetas directly instead of fitting
      // population native parameters and inverting them (there is no single
      // population `a` to invert when an argument varies by subject).  Skipping
      // the declaration outright therefore stood down the only thing able to
      // help it.
      //
      // The copula loop further down is deliberately not gated on this either:
      // the correlation is a property of the raw latent block, which a
      // covariate on a family argument does not touch.
      // ...and the family MLE half still obeys its own control: entering the
      // loop for Q2's sake must not switch the MLE back on.
      //
      // It also stands down for a declaration Q2 OWNS, whatever the control
      // says.  The PARTITION decides ownership, not the flag -- that is the
      // whole point -- and two owners on one parameter is named at the bottom
      // of this file as "the fight this whole area keeps losing".  Measured
      // when both were allowed to run: direct went from MARE 3.71% (MLE alone,
      // cadence 20) to 17.98%, while Q2 alone reaches 3.27%.
      //
      // This is not a loss of capability.  On the direct route the family MLE
      // fits NATIVE parameters and inverts them onto the thetas; Q2 maximizes
      // the same likelihood over the thetas directly, with no inversion and no
      // spread guard to satisfy.  It is the better instrument for the same job.
      bool famUsable = etaDistOn && !etaDistAllQ2(k) &&
        !(etaDistUsable.n_elem == (unsigned int)etaDistNdist &&
          etaDistUsable(k) == 0);
      int fam = etaDistFam(k);
      int na = rxEtaDistNarg(fam);
      if (na <= 0) continue;
      double a0[4];
      for (int i = 0; i < na; ++i) a0[i] = etaDistArgs(k, i);
      // This declaration's covariate columns, if any.  nSym == 0 is the
      // no-covariate case and everything below collapses to what it was.
      const arma::mat *cvK = (k < (int)etaDistCov.size()) ? &etaDistCov[(size_t)k] : NULL;
      int nSym = (cvK != NULL && cvK->n_rows > 0) ? (int)cvK->n_cols : 0;
      std::vector<double> ev; ev.reserve(w[(size_t)k].size());
      // `rec` is nRec x nSym ROW-MAJOR, which is the layout rxEtaDistLoglikObj
      // indexes as rec[r*nSym + c].
      std::vector<double> rec;
      if (nSym > 0) rec.reserve(w[(size_t)k].size() * (size_t)nSym);
      // A covariate declaration has NO population argument set -- `rate` here
      // is 1/(exp(lclrv)*exp(lclm + bWT*log(WT/70))) and the R side leaves it
      // NA precisely because there is no single value to give it.  So the eta
      // SAMPLE cannot be decoded with etaDistArgs either: rxEtaDistQ(fam, u,
      // NA) is non-finite for every draw, every one is dropped, and `ev` comes
      // back EMPTY -- which silently skipped the whole objective below for the
      // one declaration it exists to serve.  Decode per record instead, from
      // the current thetas and that subject's covariate row.
      std::vector< std::vector<etaDistTok> > rpnQ;
      std::vector<double> stQ, valsQ;
      bool perRec = false;
      if (nSym > 0) {
        int nthQ = etaDistNth(k);
        if (nthQ > 0 && k < (int)etaDistExprs.size() &&
            (int)etaDistExprThetas[(size_t)k].size() == nthQ &&
            k < (int)etaDistCovNames.size() &&
            (int)etaDistCovNames[(size_t)k].size() == nSym) {
          std::vector<std::string> pv = etaDistExprThetas[(size_t)k];
          for (int c = 0; c < nSym; ++c) pv.push_back(etaDistCovNames[(size_t)k][(size_t)c]);
          if (rxEtaDistLoglikParse(etaDistExprs[(size_t)k], pv, rpnQ)) {
            stQ.assign((size_t)nthQ, 0.0);
            for (int t = 0; t < nthQ; ++t) {
              int c = etaDistPhi0Col(k, t);
              stQ[(size_t)t] = (c >= 0 && c < nphi0) ? mprior_phi0(0, c) : 0.0;
            }
            valsQ.assign((size_t)(nthQ + nSym), 0.0);
            for (int t = 0; t < nthQ; ++t) valsQ[(size_t)t] = stQ[(size_t)t];
            perRec = (rpnQ.size() == (size_t)na);
          }
        }
        // Without a usable parse there is nothing to decode the sample with, so
        // leave it empty rather than fabricating one from NA arguments.
        if (!perRec) nSym = 0;
      }
      // On the DIRECT route the sample is ALREADY on the eta scale -- there is
      // no latent and no decoder -- so `w` must be taken as it stands.
      //
      // Decoding it anyway is not a small error.  Q(Phi(eta)) on a gamma sample
      // whose mean is 5.1 evaluates Phi(5.1) = 1 - 1.7e-7, which the guard below
      // clamps to 1 - 1e-15, and the inverse CDF there is the extreme upper
      // tail: measured, the M-step came back with lclm = 16.93 against a truth
      // of 1.63 while the eta SAMPLE it was fitted to had mean 5.098 against a
      // truth of 5.104.  The sample was right and the fit to it was not.
      const bool edSampleIsEta = etaDistDirectOn();
      for (size_t r = 0; r < w[(size_t)k].size(); ++r) {
        // same boundary guard phiU() applies: pnorm saturates to 0/1 in double
        // precision and an inverse CDF there is +/-Inf
        double u = R::pnorm(w[(size_t)k][r], 0.0, 1.0, 1, 0);
        if (u < 1e-15) u = 1e-15; else if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
        double aR[4];
        const double *aUse = a0;
        if (perRec) {
          unsigned int sj = (unsigned int)(r % (size_t)N);
          if (sj >= cvK->n_rows) break;
          for (int c = 0; c < nSym; ++c) {
            valsQ[(size_t)((int)stQ.size() + c)] = (*cvK)(sj, (unsigned int)c);
          }
          bool okA = true;
          for (int q = 0; q < na; ++q) {
            aR[q] = etaDistExprEval(rpnQ[(size_t)q], valsQ.data(), (int)valsQ.size());
            if (!std::isfinite(aR[q])) { okA = false; break; }
          }
          if (!okA) continue;
          aUse = aR;
        }
        double e = edSampleIsEta ? w[(size_t)k][r] : rxEtaDistQ(fam, u, aUse);
        if (!std::isfinite(e)) continue;
        // ...and, on the DIRECT route, in the family's SUPPORT.
        //
        // A decoded eta is in support by construction -- rxEtaDistQ is the
        // family's own quantile function -- so the finiteness test above was
        // sufficient while the cdf route was the only one.  A sampled eta is
        // not: phiM is initialized on the Gaussian scale, so a positive-support
        // family starts some rows at or below zero, and a NEGATIVE FINITE value
        // passes the test above unchanged.  It then reaches rxode2ll's density,
        // which THROWS rather than returning -Inf:
        //
        //   gamma_lpdf: Random variable is -0.951673, but must be positive finite!
        //
        // killing the fit.  Dropping the record is the same treatment a
        // non-finite one already gets, and it is the right one: a draw outside
        // the support carries no information about the family's parameters.
        if (edSampleIsEta) {
          double lo, hi;
          rxEtaDistBounds(fam, aUse, &lo, &hi);
          if ((R_finite(lo) && e <= lo) || (R_finite(hi) && e >= hi)) continue;
        }
        // ONE loop, dropped TOGETHER.  Pushing the covariate in a second pass
        // over the same range would keep every record that this one skips and
        // shift the whole column by however many etas came back non-finite --
        // silently, since both vectors would still look well formed.
        ev.push_back(e);
        if (nSym > 0) {
          // phiM stacks chains: row r is subject r % N (see phiM_weighted).
          unsigned int subj = (unsigned int)(r % (size_t)N);
          if (subj >= cvK->n_rows) { rec.clear(); nSym = 0; break; }
          for (int c = 0; c < nSym; ++c) rec.push_back((*cvK)(subj, (unsigned int)c));
        }
      }
      // Guard on the assumption the whole step rests on: the latent is
      // standard normal BY CONSTRUCTION, so the only reason the pooled draws
      // depart from that is information -- or a chain that has not mixed.  The
      // step cannot tell those apart, and acting on the second is a runaway
      // (see etaDistStart).  A spread far from 1 is therefore not signal to be
      // fitted but a sample not yet worth fitting, so skip it.  Self-tuning,
      // unlike an iteration count, and it degrades the right way: a chain that
      // never settles simply never updates, which the M-step counter reports
      // rather than hides.
      double lsd = NA_REAL;
      bool famTr = (getenv("NLMIXR2_ETADIST_OPT") != NULL);
      // Measure unconditionally (lo=0, hi=inf just fills lsd), so the trace
      // reports the spread even on an iteration that rejects, and so the
      // trajectory is recorded on every attempt rather than only on the ones
      // that pass.  Running unguarded should not also mean running blind.
      //
      // ON THE SCALE THE BAND IS WRITTEN FOR.  etaDistSdLo/etaDistSdHi bound a
      // LATENT -- the comment above says so: "the latent is standard normal BY
      // CONSTRUCTION".  On the direct route `ev` holds the ETAS, which are
      // standard normal by no construction at all, so a positive-support
      // family with a mean of 5 measures 35-45 against a band that ends at 5.
      // Every attempt then rejects, and what rejects with it is not just the
      // family MLE but etaDistQ2PairStep -- the joint maximization over the
      // pair's marginals AND atanh(rho) together, which is the only step that
      // estimates a declared copula properly.  Measured on Bauer's gamma4:
      // 250 iterations, okK=0 okJ=0 every one, rxCor returned its ini().
      //
      // z = qnorm(F(eta)) is the latent, and it is the same transform
      // rxEtaDistPairLogD uses for the copula term.  The copula loop below
      // measures the same way, which matters: both read etaDistSdPrev, so a
      // baseline set on one scale and compared on the other spans nothing.
      rxEtaDistSpreadOk(w[(size_t)k], 0.0, R_PosInf, &lsd);
      double sdPrevWas = (k < (int)etaDistSdPrev.n_elem) ? etaDistSdPrev(k) : NA_REAL;
      bool spreadOk = (etaDistSpreadGuard == 0) ||
        etaDistSpreadSettled(k, lsd);
      double aNew[4];
      for (int i = 0; i < na; ++i) aNew[i] = a0[i];
      // famOff: in the observation-likelihood mode the families are estimated
      // by the objective below (or by refinePhi0Lik), not by this MLE.
      bool mleOk = !famOff && famUsable && spreadOk && rxEtaDistMle(fam, ev, aNew);
      if (famTr) {
        double relCh = (std::isfinite(sdPrevWas) && sdPrevWas > 0) ?
          std::fabs(lsd - sdPrevWas)/sdPrevWas : NA_REAL;
        RSprintf("[fam] it=%d k=%d latentSd=%.4f prev=%.4f relCh=%.4f tol=%.3g "
                 "cap=[%.2f,%.2f] guard=%d spreadOk=%d mleOk=%d\n",
                 (int)kiter, k, lsd, sdPrevWas, relCh, etaDistSdTol,
                 etaDistSdLo, etaDistSdHi,
                 etaDistSpreadGuard, (int)spreadOk, (int)mleOk);
        // The phi0 column each of THIS declaration's thetas writes into, printed
        // where the loop actually runs.  An earlier version of this sat beside
        // the [phi0] trace, which is emitted BEFORE the control is read, so it
        // reported zeros for everything and said the machinery was off while
        // this loop was plainly running.
        {
          std::string tp;
          for (int t = 0; t < etaDistNth(k); ++t) {
            char b[32];
            snprintf(b, sizeof(b), "%d ", etaDistPhi0Col(k, t));
            tp += b;
          }
          RSprintf("[phi0map] k=%d nth=%d thetaPhi0Rows=%d cols=%d cols={%s} "
                   "famUsable=%d\n",
                   k, etaDistNth(k), (int)etaDistThetaPhi0.n_rows,
                   (int)etaDistThetaPhi0.n_cols, tp.c_str(), (int)famUsable);
        }
      }
      if (etaDistDebug && (kiter % 10 == 0 || kiter < 2)) {
        double ws = 0, ws2 = 0, es = 0, es2 = 0;
        size_t nw = w[(size_t)k].size();
        for (size_t r = 0; r < nw; ++r) { ws += w[(size_t)k][r]; ws2 += w[(size_t)k][r]*w[(size_t)k][r]; }
        for (size_t r = 0; r < ev.size(); ++r) { es += ev[r]; es2 += ev[r]*ev[r]; }
        double wm = nw ? ws/nw : 0, wsd = nw ? std::sqrt(std::max(0.0, ws2/nw - wm*wm)) : 0;
        double em = ev.size() ? es/ev.size() : 0;
        double esd = ev.size() ? std::sqrt(std::max(0.0, es2/ev.size() - em*em)) : 0;
        RSprintf("[etaDist k=%d it=%d] nw=%d latent m=%.4f sd=%.4f | eta n=%d m=%.4g sd=%.4g | a0=(%.4g,%.4g) -> aNew=(%.4g,%.4g) ok=%d\n",
                 k, (int)kiter, (int)nw, wm, wsd, (int)ev.size(), em, esd,
                 a0[0], na>1?a0[1]:0.0, aNew[0], na>1?aNew[1]:0.0, (int)mleOk);
        if (!spreadOk) RSprintf("[etaDist k=%d it=%d] SKIPPED: latent sd %.4f outside [%.2f, %.2f]\n",
                                k, (int)kiter, lsd, etaDistSdLo, etaDistSdHi);
      }
      // saemControl(etaDistLoglik=): the general objective (see the design in
      // R/etaDistMstep.R).  Maximizes the family log-likelihood over the THETAS
      // directly, so there is no population-level native parameter set to fit
      // and none to invert -- which is what lets a covariate on a distribution
      // parameter be represented at all.
      //
      // This first increment covers the nSym == 0 case: no covariate, one
      // record per sampled draw, equal weights.  That is design test T1 -- it
      // must reproduce the MLE-plus-inversion answer, which is the check that
      // the objective is right before the per-record plumbing is added on top.
      // A copula-linked PAIR is one joint problem, taken once.
      //
      // Scoring the two marginals separately and carrying a correlation
      // alongside them is a different estimator -- the marginals' MLE from the
      // joint likelihood equals the separate marginal MLEs only at rho == 0 --
      // and it leaves rho to a moment statistic rather than to the copula's own
      // density.  etaDistQ2PairStep() does both together, with atanh(rho) in
      // the same parameter vector.
      //
      // Entered from the LOWER member so the pair is handled once whichever
      // order the loop reaches them in, and `continue` because the joint step
      // has already written back everything the per-declaration step below
      // would have.
      {
        double rhoPair = 0.0;
        int part = etaDistPartnerOf(k, &rhoPair);
        // NOT gated on spreadOk.  That guard exists for ONE consumer: the
        // family MLE, which fits native parameters to the sample and then
        // INVERTS them onto the user's thetas, consulting no objective on the
        // way -- a bad sample there maps straight to an extreme theta.  This
        // step is not that.  It is an n1qn1 maximization of the declared
        // prior that accepts only on `gEdBest < f0` and then damps by
        // pas(kiter), so a sample not worth fitting simply fails to improve
        // and nothing moves.  NoLimits guards its equivalent step with
        // nothing but isfinite on the optimizer's result, for the same reason.
        //
        // And on the direct route the family MLE cannot run at all --
        // famUsable is `etaDistOn && !etaDistAllQ2(k)`, and every declared
        // theta there is Q2 -- so spreadOk was gating only steps that do not
        // need it.  That is what blocked Bauer's gamma2 and gamma4 for entire
        // fits (lclm 1.244 and 5.269 against a truth of 1.630).
        if (part > k && etaDistAllQ2(k) && etaDistAllQ2(part)) {
          if (etaDistQ2PairStep(k, part, kiter, pas)) { moved = true; continue; }
        } else if (part >= 0 && part < k && etaDistAllQ2(k) && etaDistAllQ2(part)) {
          continue;   // already done from the lower member
        }
      }
      // Q2: the eta-density objective.  Reached either because the user asked
      // for it (etaDistLoglik) or because the PARTITION says this declaration's
      // thetas are prior-only and nothing else can identify them.
      // spreadOk dropped here too, and for the same reason: this block
      // accepts only on `ynew < f0` and damps the accepted step.
      if ((etaDistLoglik || etaDistAllQ2(k)) && !ev.empty()) {
        int nth = etaDistNth(k);
        if (nth > 0 && k < (int)etaDistExprs.size() &&
            (int)etaDistExprThetas[(size_t)k].size() == nth) {
          // Parse against thetas THEN the covariate symbols.  That order is
          // the contract: rxEtaDistLoglikObj lays out vals[0..nth) as the
          // candidate thetas and vals[nth..nth+nSym) as this record's symbols,
          // and the parser stores each symbol's index into that same flat
          // vector.  Without the covariate names the parse simply DECLINES on
          // the unknown symbol, which is how a covariate silently cost the
          // whole objective before this.
          std::vector<std::string> pvars = etaDistExprThetas[(size_t)k];
          if (nSym > 0 && k < (int)etaDistCovNames.size()) {
            for (size_t c = 0; c < etaDistCovNames[(size_t)k].size(); ++c) {
              pvars.push_back(etaDistCovNames[(size_t)k][c]);
            }
          }
          std::vector< std::vector<etaDistTok> > rpn;
          if ((int)pvars.size() == nth + nSym &&
              rxEtaDistLoglikParse(etaDistExprs[(size_t)k], pvars, rpn)) {
            std::vector<double> wt(ev.size(), 1.0);
            gEdFam = fam; gEdRpn = &rpn; gEdNth = nth; gEdNSym = nSym;
            gEdNRec = (int)ev.size();
            gEdRec = (nSym > 0) ? rec.data() : NULL;
            gEdEta = ev.data(); gEdWt = wt.data();
            std::vector<double> st((size_t)nth), step((size_t)nth),
              xmin((size_t)nth);
            for (int t = 0; t < nth; ++t) {
              int c = etaDistPhi0Col(k, t);
              st[(size_t)t] = (c >= 0 && c < nphi0) ? mprior_phi0(0, c) : 0.0;
              double a = std::fabs(st[(size_t)t]);
              step[(size_t)t] = (a > 1e-8) ? 0.2 * a : 0.1;
            }
            double f0 = gEdObj(st.data());
            if (f0 < 1e299) {
              // n1qn1: quasi-Newton on the exact gradient.  NOT nelder_fn
              // (measured badly enough elsewhere here that it is not a
              // defensible default) and NOT newuoa, which is not thread-safe
              // and would need Rcpp::Function to reach -- both disqualifying
              // for a step whose point is to run inside an OpenMP region.
              gEdBest = R_PosInf;
              gEdBestPar.assign(st.begin(), st.end());
              gEdN1Bad = 0; gEdN1Evals = 0;
              std::vector<double> gg((size_t)nth, 0.0);
              std::vector<double> zm((size_t)(nth*(nth+13)/2 + 1), 0.0);
              std::vector<double> var((size_t)nth, 0.1);
              double fN = 0.0, eps = 1e-8;
              int nn = nth, mode = 1, niter = 200, nsim = 200, impr = 0,
                izs = 0, idz = 0; float rzs = 0; double dzs = 0;
              if (n1qn1_ != NULL) {
                n1qn1_(gEdN1Cost, &nn, st.data(), &fN, gg.data(), var.data(),
                       &eps, &mode, &niter, &nsim, &impr, zm.data(),
                       &izs, &rzs, &dzs, &idz);
              }
              double ynew = gEdBest;
              for (int t = 0; t < nth; ++t) xmin[(size_t)t] = gEdBestPar[(size_t)t];
              if (!gEdN1Bad && gEdN1Evals > 0 && ynew < f0) {

                for (int t = 0; t < nth; ++t) {
                  int c = etaDistPhi0Col(k, t);
                  if (c < 0 || c >= nphi0 || !std::isfinite(xmin[(size_t)t])) continue;
                  double cur = mprior_phi0(0, c);
                  double v = cur + pas(kiter) * (xmin[(size_t)t] - cur);
                  if (std::isfinite(v)) { mprior_phi0.col(c).fill(v); moved = true; }
                }
                gEdRpn = NULL; gEdEta = NULL; gEdWt = NULL;
                if ((int)etaDistFiredK.size() == etaDistNdist) etaDistFiredK[(size_t)k] = 1;
                continue;   // thetas set directly; no MLE, no inversion
              }
            }
            gEdRpn = NULL; gEdEta = NULL; gEdWt = NULL;
          }
        }
      }
      if (!mleOk) continue;
      // etaDistDebug >= 2 observes without acting: the trace above then shows
      // what the M-step WOULD have seen over an otherwise ordinary fit, which
      // is the only way to watch the latent spread settle without the M-step's
      // own updates perturbing it.
      if (etaDistDebug > 1) continue;
      // stochastic-approximation damping, as every other M-step here does
      for (int i = 0; i < na; ++i) {
        double cur = etaDistArgs(k, i);
        double v = cur + pas(kiter) * (aNew[i] - cur);
        if (std::isfinite(v)) { etaDistArgs(k, i) = v; moved = true; }
      }
      if ((int)etaDistFiredK.size() == etaDistNdist) etaDistFiredK[(size_t)k] = 1;
    }
    // copula correlations: closed form, no search.  Runs independently of the
    // family M-step -- see etaDistCorOn.
    // The copula closed form stands down when the gradient step owns rxCor:
    // two owners on one parameter, against different objectives, is the fight
    // this whole area keeps losing.
    // Stands down only when the SEARCH owns rxCor.  "observed" and "analytic"
    // are closed forms and run here in every mode -- including the
    // observation-likelihood mode, where they are the only thing that keeps the
    // correlation off the boundary the search walks to.
    // method 2 is owned by the search; method 3 by the every-iteration damped
    // adjustment above.  Only 0 and 1 are applied here.
    bool corBySearch = (etaDistObsLik() && etaDistCorMethod == 2) ||
      etaDistCorMethod == 3;
    for (int k = 0; etaDistCorOn && !corBySearch && k < etaDistNdist; ++k) {
      int j = etaDistCorWith(k);
      if (j < 0) continue;
      // Same spread guard the family fits get.  This loop used to bypass it
      // entirely, so the correlation was updated from exactly the unmixed draws
      // the guard exists to reject -- and being pinned at its clamp makes the
      // partner's latent numerically equal to its partner's, which breaks BOTH
      // family fits, not just the correlation.
      // EVERY firing, no modulus.  The M-step already has its own cadence
      // (etaDistStart + n*etaDistEvery); an independent `kiter % 20` cannot be
      // relied on to coincide with it, and when it does not the trace is simply
      // absent -- which reads as "nothing happened" rather than "not printed".
      // That has now cost two debugging rounds in this file.
      // The spread guard exists because the RAW product-moment estimator runs
      // away when the latents are over-dispersed.  The normalized sufficient
      // statistic cannot: it is a ratio bounded in [-1, 1] by construction, and
      // an inflated diagonal ATTENUATES it toward zero rather than inflating it.
      // So method 3 is exempt -- and that matters, because the guard rejecting
      // every iteration is what made three estimators produce byte-identical
      // fits while none of them ran.
      // Same test the family fits get, against the same baseline: both loops
      // read etaDistSdPrev, which does not move until the end of the M-step.
      double lsdK = NA_REAL, lsdJ = NA_REAL;
      rxEtaDistSpreadOk(w[(size_t)k], 0.0, R_PosInf, &lsdK);
      rxEtaDistSpreadOk(w[(size_t)j], 0.0, R_PosInf, &lsdJ);
      // Both, ALWAYS, before combining.  Through a short-circuiting && a false
      // from k would skip j entirely, j's spread would never be staged for this
      // attempt, and the advance below would carry j's value from some EARLIER
      // attempt into the baseline -- so j's next comparison spans two gaps or
      // more.  That mis-rejects a settled j, and worse, an oscillating j whose
      // frozen baseline it keeps matching passes as "settled".  Reachable
      // whenever the family loop above does not run: etaDistMstep=FALSE with
      // etaDistCorMstep=TRUE, and the observation-likelihood mode (famOff).
      bool okK = etaDistSpreadSettled(k, lsdK);
      bool okJ = etaDistSpreadSettled(j, lsdJ);
      // The DIRECT route is exempt, because on it this guard protects nothing.
      //
      // The guard exists for the family MLE -- fit native parameters to the
      // sample, then invert them onto the user's thetas, consulting no
      // objective on the way.  That step cannot run here at all: famUsable is
      // `etaDistOn && !etaDistAllQ2(k)` and on the direct route every declared
      // theta is Q2.  There is no MLE to protect, and NoLimits, which has no
      // fit-then-invert step anywhere, has no spread guard anywhere either --
      // its M-step skips on isfinite alone.
      //
      // What the guard DID do here was block the estimators, all of which are
      // already safe by construction: each closed form is clamped into
      // [-0.999, 0.999] and the accepted value is damped by pas(kiter).  The
      // runaway that motivated it -- Bauer's g1 reaching 0.999 -- was the RAW
      // product moment on over-dispersed LATENTS, on the cdf route, which is
      // also why method 3 was already exempt on the same reasoning.
      //
      // Measured: `w` holds the etas here, so a positive-support family with a
      // mean of 5 reports a spread of 35-45 against a band ending at 5, okK
      // and okJ were 0 on every one of 250 iterations, and rxCor returned its
      // ini() value of 0.600 with truth 0.500.
      bool corSpreadOk = (etaDistSpreadGuard == 0) || etaDistDirectOn() ||
        (etaDistCorMethod == 3) || (okK && okJ);
      if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
        int cj2 = etaDistLatent(j), ck2 = etaDistLatent(k);
        if (cj2 >= 0 && ck2 >= 0 && cj2 < (int)phiM.n_cols &&
            ck2 < (int)phiM.n_cols) {
          double zjj = 0, zkk = 0, zjk = 0, njj = 0, nkk = 0, njk = 0;
          unsigned int nr2 = phiM.n_rows;
          for (unsigned int rr = 0; rr < nr2; ++rr) {
            double zj = phiM(rr, cj2), zk2 = phiM(rr, ck2);
            zjj += zj*zj; zkk += zk2*zk2; zjk += zj*zk2;
          }
          for (size_t rr = 0; rr < w[(size_t)k].size() && rr < w[(size_t)j].size(); ++rr) {
            double a = w[(size_t)j][rr], b = w[(size_t)k][rr];
            njj += a*a; nkk += b*b; njk += a*b;
          }
          double nn = (double)nr2;
          double rho0 = etaDistRho(k);
          double l21 = std::isfinite(rho0) ? rho0 : 0.0;
          double l22 = std::sqrt(std::max(0.0, 1.0 - l21*l21));
          // S_n = L S_z L' for the 2x2 block, rebuilt from S_z
          double bjj = zjj/nn;
          double bkk = l21*l21*(zjj/nn) + 2*l21*l22*(zjk/nn) + l22*l22*(zkk/nn);
          double bjk = l21*(zjj/nn) + l22*(zjk/nn);
          RSprintf("[Sz] it=%d k=%d  S_z/N: jj=%.4f kk=%.4f jk=%+.4f (I would be 1,1,0)\n",
                   (int)kiter, k, zjj/nn, zkk/nn, zjk/nn);
          RSprintf("[Sz] it=%d k=%d  S_n/N direct jj=%.4f kk=%.4f jk=%+.4f | via L*Sz*L' jj=%.4f kk=%.4f jk=%+.4f\n",
                   (int)kiter, k, njj/nn, nkk/nn, njk/nn, bjj, bkk, bjk);
          RSprintf("[Sz] it=%d k=%d  rho: current=%+.4f  suff-stat=%+.4f\n",
                   (int)kiter, k, rho0,
                   (njj > 0 && nkk > 0) ? njk/std::sqrt(njj*nkk) : NA_REAL);
        }
      }
      if (getenv("NLMIXR2_ETADIST_OPT") != NULL) {
        RSprintf("[Sz] it=%d k=%d  spreadOk=%d  (0 means this loop returns "
                 "WITHOUT running any estimator)\n",
                 (int)kiter, k, (int)corSpreadOk);
      }
      // The guard stays, but note what its failing MEANS: no estimator runs, and
      // rxCor is then left to whatever else owns it -- which is why three
      // different estimators produced byte-identical fits.
      if (!corSpreadOk) continue;

      // "observed" is the product-moment correlation of the latent pair;
      // "analytic" the rank-based Gaussian-copula identity.  Both are closed
      // forms evaluated after the distributional thetas moved; "optimize" is
      // handled in etaDistGradStep() instead and skips this loop.
      // SUFFICIENT-STATISTIC DIAGNOSTIC.
      //
      // R enters the complete-data likelihood only through S_n = sum(n n'),
      // and n = L z is not sampled -- z is (phiM's rxz.* columns).  So
      //
      //     S_n = L * S_z * L'
      //
      // and if S_z is the identity then S_n = L L' = R_old and the M-step
      // returns exactly what it was handed, whatever estimator is used.  ALL
      // the information about R lives in the departure of S_z from I, so print
      // it: three numbers per pair settle whether the update has anything to
      // work with, or whether the chain simply is not moving z.
      //
      // Also prints S_n both ways -- accumulated directly from w, and rebuilt
      // as L S_z L' -- because those must agree, and agreeing is what shows
      // "posterior" IS the normalized sufficient statistic rather than merely
      // resembling it.
      double r;
      if (etaDistCorMethod == 1) {
        r = rxEtaDistCorSpearman(w[(size_t)j], w[(size_t)k]);
      } else if (etaDistCorMethod == 3) {
        // straight from the statistic refreshed this iteration
        r = ((int)etaDistCorSuff.n_elem == etaDistNdist) ?
          etaDistCorSuff(k) : rxEtaDistCorPost(w[(size_t)j], w[(size_t)k]);
      } else {
        r = rxEtaDistCorMle(w[(size_t)j], w[(size_t)k]);
      }
      if (etaDistDebug && (kiter % 10 == 0 || kiter < 2)) {
        RSprintf("[etaDist cor k=%d<-j=%d it=%d] rhoCur=%.4f rhoNew=%.4f\n",
                 k, j, (int)kiter, etaDistRho(k), r);
      }
      if (!std::isfinite(r)) continue;
      if (etaDistDebug > 1) continue;
      double cur = etaDistRho(k);
      double v = cur + pas(kiter) * (r - cur);
      if (std::isfinite(v)) { etaDistRho(k) = v; moved = true; etaDistCorFired = true; }
    }
    // ONE advance per M-step attempt, after every loop that consulted the
    // baseline, so the next attempt compares across a full etaDistEvery gap.
    for (int k = 0; k < (int)etaDistSdCur.n_elem; ++k) {
      if (std::isfinite(etaDistSdCur(k))) etaDistSdPrev(k) = etaDistSdCur(k);
    }
    return moved;
  }

  // Robbins-Monro adaptation of the random-walk scale toward the target
  // acceptance rate, exactly saemix's rule (R/main_estep.R:59 and :95):
  //
  //   domega2 <- domega2 * (1 + stepsize.rw*(nbc2/nt2 - proba.mcmc))
  //
  // NONMEM does the same thing through IACCEPT (Bauer's gamma streams set
  // IACCEPT=0.3).  nlmixr2 had neither: it carried saemix's INITIAL value
  // (rw.init = 0.5, here `rmcmc`) frozen for the whole fit, so a chain whose
  // acceptance was far from target never corrected, and the step could not
  // grow to reach the tails of the latent normal.
  //
  // The per-iteration factor is clamped to [0.5, 2] and the accumulated scale
  // to [1e-3, 1e3] so one unlucky iteration cannot collapse or explode the
  // proposal.
  // Accumulate one proposal block's acceptances into this iteration's
  // diagnostic counters.  `acc` holds the accepted ROW indices; phiM stacks the
  // nmc chains, so row r belongs to subject r % N.  Kernels are indexed 1,2,3
  // and 4 (mode 1B) -> columns 0..3.
  void mcmcRecordAccept(int method, const uvec &acc, int nM) {
    if (method < 1 || method > 4 || nM <= 0) return;
    if (mcmcAccNum.n_elem != 4) { mcmcAccNum.zeros(4); mcmcAccDen.zeros(4); }
    mcmcAccNum(method - 1) += (double)acc.n_elem;
    mcmcAccDen(method - 1) += (double)nM;
    if (N > 0) {
      if ((int)mcmcAccById.n_elem != N) mcmcAccById.zeros(N);
      if ((int)mcmcAccByIdK2.n_elem != N) mcmcAccByIdK2.zeros(N);
      if ((int)mcmcAccByIdK3.n_elem != N) mcmcAccByIdK3.zeros(N);
      for (unsigned int j = 0; j < acc.n_elem; ++j) {
        arma::uword sIdx = acc(j) % (arma::uword)N;
        mcmcAccById(sIdx) += 1.0;
        if (method == 2) mcmcAccByIdK2(sIdx) += 1.0;
        else if (method == 3) mcmcAccByIdK3(sIdx) += 1.0;
      }
      mcmcAccByIdTrials += 1.0;
      // Trials per subject in this block is the chain count.  Only the
      // random-walk kernels are counted: kernel 1 draws from the prior, so its
      // acceptance is not a function of any step size and NONMEM does not
      // adapt on it either.
      double perSubj = (double)(nM / (N > 0 ? N : 1));
      if (method == 2) mcmcAccByIdK2Trials += perSubj;
      else if (method == 3) mcmcAccByIdK3Trials += perSubj;
    }
  }

  // Close out one iteration's diagnostics and start the next.  Called once per
  // iteration, after the whole MCMC block, so it sees the chain as the M-step
  // will see it.
  void mcmcCloseIter(unsigned int kiter, unsigned int niterTot, const mat &phiCur) {
    if ((int)mcmcAccTrace.n_rows != (int)niterTot) {
      mcmcAccTrace.set_size(niterTot, 4);    mcmcAccTrace.fill(NA_REAL);
      mcmcStuckTrace.set_size(niterTot, 1);  mcmcStuckTrace.fill(NA_REAL);
      phiSdTrace.set_size(niterTot, phiCur.n_cols);   phiSdTrace.fill(NA_REAL);
      phiAcfTrace.set_size(niterTot, phiCur.n_cols);  phiAcfTrace.fill(NA_REAL);
    }
    if (kiter >= niterTot) return;
    if (mcmcAccNum.n_elem == 4) {
      for (int k = 0; k < 4; ++k) {
        if (mcmcAccDen(k) > 0) mcmcAccTrace(kiter, k) = mcmcAccNum(k) / mcmcAccDen(k);
      }
    }
    if ((int)mcmcAccById.n_elem == N && N > 0) {
      double nStuck = 0.0;
      for (int i = 0; i < N; ++i) if (mcmcAccById(i) <= 0.0) nStuck += 1.0;
      mcmcStuckTrace(kiter, 0) = nStuck / (double)N;
    }
    for (unsigned int c = 0; c < phiCur.n_cols; ++c) {
      vec col = phiCur.col(c);
      if (col.n_elem > 1) phiSdTrace(kiter, c) = arma::stddev(col);
      if (phiMprevIter.n_rows == phiCur.n_rows &&
          phiMprevIter.n_cols == phiCur.n_cols) {
        vec pv = phiMprevIter.col(c);
        double sd1 = arma::stddev(col), sd0 = arma::stddev(pv);
        if (sd1 > 1e-12 && sd0 > 1e-12) {
          phiAcfTrace(kiter, c) = arma::as_scalar(arma::cor(pv, col));
        }
      }
    }
    phiMprevIter = phiCur;
    // saemControl(iacceptPerId=) adapts here, on the pooled-over-blocks rate,
    // and must run before the counters are cleared.
    adaptRwPerIdIter();
    mcmcAccNum.zeros(4); mcmcAccDen.zeros(4);
    if (N > 0) { mcmcAccById.zeros(N); mcmcAccByIdK2.zeros(N); mcmcAccByIdK3.zeros(N); }
    mcmcAccByIdTrials = 0.0;
    mcmcAccByIdK2Trials = 0.0;
    mcmcAccByIdK3Trials = 0.0;
  }

  // Per-subject form of adaptRw (saemControl(iacceptPerId=)).  `acc` holds the
  // ACCEPTED ROW indices of this proposal block; phiM stacks the nmc chains
  // vertically, so row r belongs to subject r % N and each subject gets nmc
  // trials per block.  Same Robbins-Monro rule, same clamps, applied to one
  // scalar per subject rather than one per coordinate.
  // Per-subject form of adaptRw (saemControl(iacceptPerId=)), applied ONCE PER
  // ITERATION against acceptances pooled over every random-walk block.
  //
  // Pooling is the whole point.  One block gives a subject only `nmc` trials,
  // so its acceptance rate is one of {0, 1/nmc, ..., 1} -- far too coarse to
  // drive a multiplicative update.  Over a whole iteration a subject sees
  // (1 + nphi) * nu blocks, which is a usable estimate.  Same Robbins-Monro
  // rule and same clamps as adaptRw, applied to one scalar per subject.
  void adaptRwPerIdIter() {
    if (!iacceptPerId || N <= 0) return;
    // One scale per SUBJECT and per KERNEL, each against that kernel's own
    // target.  Pooling the two kernels into one rate would compare a mixture of
    // a ~0.234-optimal and a ~0.44-optimal proposal against a single number and
    // satisfy neither -- the per-subject form of the defect that froze the
    // pooled path.
    struct { vec *lam; const vec *acc; double trials; double target; } arms[4] = {
      { &rwLam1,  &mcmcAccByIdK2, mcmcAccByIdK2Trials, iaccept       },
      { &rwLam0,  &mcmcAccByIdK2, mcmcAccByIdK2Trials, iaccept       },
      { &rwLam1b, &mcmcAccByIdK3, mcmcAccByIdK3Trials, iacceptSingle },
      { &rwLam0b, &mcmcAccByIdK3, mcmcAccByIdK3Trials, iacceptSingle }
    };
    for (int a = 0; a < 4; ++a) {
      vec *lam = arms[a].lam;
      const vec *acc = arms[a].acc;
      if (lam == nullptr || (int)lam->n_elem != N) continue;
      if ((int)acc->n_elem != N || !(arms[a].trials > 0.0)) continue;
      if (!(arms[a].target > 0.0)) continue;
      for (int i = 0; i < N; ++i) {
        double rate = (*acc)(i) / arms[a].trials;
        double f = 1.0 + stepsizeRw * (rate - arms[a].target);
        if (f < 0.5) f = 0.5;
        else if (f > 2.0) f = 2.0;
        double v = (*lam)(i) * f;
        if (v < 1e-3) v = 1e-3;
        else if (v > 1e3) v = 1e3;
        (*lam)(i) = v;
      }
    }
  }

  void adaptRw(vec *rwScale, const uvec &cols, double accRate, double target) {
    if (rwScale == nullptr || target <= 0.0) return;
    double f = 1.0 + stepsizeRw * (accRate - target);
    if (f < 0.5) f = 0.5;
    else if (f > 2.0) f = 2.0;
    for (unsigned int j = 0; j < cols.n_elem; ++j) {
      unsigned int c = cols(j);
      if (c >= rwScale->n_elem) continue;
      double v = (*rwScale)(c) * f;
      if (v < 1e-3) v = 1e-3;
      else if (v > 1e3) v = 1e3;
      (*rwScale)(c) = v;
    }
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

  // Phi(z), clamped off 0 and 1 exactly as the decoder's phiU() is: a draw in
  // the far tail otherwise saturates in double precision and the quantile
  // function returns the support endpoint for a z the model is happy with.
  //
  // Named for the normal CDF explicitly.  etaDistKernel.h already has an
  // rxEtaDistUPhi(), and it is the FAMILY cdf F(x) -- a different function with
  // a compatible-looking call.
  static double etaDistNormCdf(double z) {
    double u = R::pnorm5(z, 0.0, 1.0, 1, 0);
    if (u < 1e-15) u = 1e-15; else if (u > 1.0 - 1e-15) u = 1.0 - 1e-15;
    return u;
  }

  // -------------------------------------------------------------------------
  // Q2 for a copula-linked PAIR: one joint problem, not two marginal ones.
  // -------------------------------------------------------------------------

  // The union of two declarations' theta names, and where each one writes back.
  // Union rather than concatenation so two declarations may SHARE a theta -- a
  // common dispersion, say -- and still be one optimization vector.
  bool etaDistPairUnion(int k, int j, std::vector<std::string> &uni,
                        std::vector<int> &col) const {
    uni.clear(); col.clear();
    if (k >= (int)etaDistExprThetas.size() || j >= (int)etaDistExprThetas.size()) {
      return false;
    }
    for (int pass = 0; pass < 2; ++pass) {
      int d = (pass == 0) ? k : j;
      const std::vector<std::string> &nm = etaDistExprThetas[(size_t)d];
      for (size_t t = 0; t < nm.size(); ++t) {
        bool dup = false;
        for (size_t q = 0; q < uni.size(); ++q) if (uni[q] == nm[t]) { dup = true; break; }
        if (dup) continue;
        int c = etaDistPhi0Col(d, (int)t);
        if (c < 0 || c >= nphi0) return false;   // cannot write it back
        uni.push_back(nm[t]);
        col.push_back(c);
      }
    }
    return !uni.empty();
  }

  // One joint step for the pair (k, j).  Returns true when it moved something.
  //
  // COVARIATES INCLUDED.  A coefficient on a declared distribution is a
  // parameter of that joint distribution exactly as rho is, so it belongs in
  // the same vector: the two declarations' covariate symbols are unioned,
  // built once into a shared per-subject record matrix, and indexed
  // consistently from both expression sets.  rxEtaDistPairLoglikObj already
  // took `nSym`/`rec` per record -- only this caller declined to fill them,
  // which left a covariate coefficient and the copula estimated by two
  // different owners against two different objectives.
  bool etaDistQ2PairStep(int k, int j, unsigned int kiter, const vec &pas) {
    if (!etaDistDirectOn()) return false;          // eta sample, not a latent
    if (!etaDistAllQ2(k) || !etaDistAllQ2(j)) return false;
    if (k >= (int)etaDistExprs.size() || j >= (int)etaDistExprs.size()) return false;
    int ck = etaDistLatent(k), cj = etaDistLatent(j);
    if (ck < 0 || cj < 0 || ck >= (int)phiM.n_cols || cj >= (int)phiM.n_cols) return false;
    std::vector<std::string> uni; std::vector<int> col;
    if (!etaDistPairUnion(k, j, uni, col)) return false;
    // The covariate symbols of BOTH members, deduplicated by name, each
    // remembering which declaration's matrix supplies its column.
    std::vector<std::string> covUni;
    std::vector<int> covDecl, covCol;
    for (int pass = 0; pass < 2; ++pass) {
      int d = (pass == 0) ? k : j;
      if (d >= (int)etaDistCovNames.size() || d >= (int)etaDistCov.size()) continue;
      const std::vector<std::string> &nm = etaDistCovNames[(size_t)d];
      if (nm.size() != (size_t)etaDistCov[(size_t)d].n_cols) return false;
      for (size_t c = 0; c < nm.size(); ++c) {
        bool dup = false;
        for (size_t q = 0; q < covUni.size(); ++q) {
          if (covUni[q] == nm[c]) { dup = true; break; }
        }
        if (dup) continue;
        covUni.push_back(nm[c]);
        covDecl.push_back(d);
        covCol.push_back((int)c);
      }
    }
    const int nSymP = (int)covUni.size();
    // rho joins the vector on the atanh scale, LAST
    double rho0 = ((int)etaDistRho.n_elem == etaDistNdist) ? etaDistRho(k) : 0.0;
    if (!std::isfinite(rho0)) rho0 = 0.0;
    if (rho0 > 0.99) rho0 = 0.99; else if (rho0 < -0.99) rho0 = -0.99;
    int nUni = (int)uni.size();
    int nth = nUni + 1;
    int rhoIdx = nUni;
    std::vector<std::string> pvars = uni;
    pvars.push_back("rxEdRho");     // never appears in an expression; a slot only
    // Symbols come AFTER every theta, because that is the layout the objective
    // reads: vals[0..nth) are the candidate parameters -- rho among them -- and
    // vals[nth..nth+nSym) this record's symbols.
    for (size_t c = 0; c < covUni.size(); ++c) pvars.push_back(covUni[c]);
    std::vector< std::vector<etaDistTok> > rpn1, rpn2;
    if (!rxEtaDistLoglikParse(etaDistExprs[(size_t)k], pvars, rpn1)) return false;
    if (!rxEtaDistLoglikParse(etaDistExprs[(size_t)j], pvars, rpn2)) return false;
    // ONE loop, dropped TOGETHER: a record is usable only when BOTH etas are,
    // or the two columns silently shift relative to each other.
    int fk = etaDistFam(k), fj = etaDistFam(j);
    const int nak = (int)etaDistArgs.n_cols;
    std::vector<double> ak((size_t)nak), aj((size_t)nak);
    for (int t = 0; t < nak; ++t) { ak[(size_t)t] = etaDistArgs(k, t); aj[(size_t)t] = etaDistArgs(j, t); }
    double lok, hik, loj, hij;
    rxEtaDistBounds(fk, &ak[0], &lok, &hik);
    rxEtaDistBounds(fj, &aj[0], &loj, &hij);
    std::vector<double> e1, e2, recP;
    e1.reserve(phiM.n_rows); e2.reserve(phiM.n_rows);
    if (nSymP > 0) recP.reserve(phiM.n_rows * (size_t)nSymP);
    for (unsigned int r = 0; r < phiM.n_rows; ++r) {
      double x1 = phiM(r, (unsigned int)ck), x2 = phiM(r, (unsigned int)cj);
      if (!std::isfinite(x1) || !std::isfinite(x2)) continue;
      if ((R_finite(lok) && x1 <= lok) || (R_finite(hik) && x1 >= hik)) continue;
      if ((R_finite(loj) && x2 <= loj) || (R_finite(hij) && x2 >= hij)) continue;
      // The covariate row FIRST, so a record with an unusable one is dropped
      // before either eta is pushed -- the objective refuses the whole sample
      // on a single non-finite symbol, and a partial push would misalign every
      // column after it.
      if (nSymP > 0) {
        if (N <= 0) return false;
        unsigned int subj = (unsigned int)(r % (size_t)N);
        bool okc = true;
        double buf[32];
        if (nSymP > 32) return false;
        for (int c = 0; c < nSymP; ++c) {
          const arma::mat &cm = etaDistCov[(size_t)covDecl[(size_t)c]];
          if (subj >= cm.n_rows) { okc = false; break; }
          double v = cm(subj, (unsigned int)covCol[(size_t)c]);
          if (!std::isfinite(v)) { okc = false; break; }
          buf[c] = v;
        }
        if (!okc) continue;
        for (int c = 0; c < nSymP; ++c) recP.push_back(buf[c]);
      }
      e1.push_back(x1); e2.push_back(x2);
    }
    if (e1.size() < 2) return false;
    std::vector<double> wt(e1.size(), 1.0);
    gEdFam = fk; gEdRpn = &rpn1; gEdFam2 = fj; gEdRpn2 = &rpn2;
    gEdNth = nth; gEdNSym = nSymP; gEdNRec = (int)e1.size();
    gEdRec = (nSymP > 0) ? recP.data() : NULL;
    gEdEta = e1.data(); gEdEta2 = e2.data(); gEdWt = wt.data();
    gEdRho = rho0; gEdRhoIdx = rhoIdx;
    std::vector<double> st((size_t)nth);
    for (int t = 0; t < nUni; ++t) st[(size_t)t] = mprior_phi0(0, col[(size_t)t]);
    st[(size_t)rhoIdx] = std::atanh(rho0);
    bool moved = false;
    double f0 = gEdObj(st.data());
    if (f0 < 1e299) {
      gEdBest = R_PosInf; gEdBestPar.assign(st.begin(), st.end());
      gEdN1Bad = 0; gEdN1Evals = 0;
      std::vector<double> gg((size_t)nth, 0.0);
      std::vector<double> zm((size_t)(nth*(nth+13)/2 + 1), 0.0);
      std::vector<double> var((size_t)nth, 0.1);
      double fN = 0.0, eps = 1e-8;
      int nn = nth, mode = 1, niter = 200, nsim = 200, impr = 0,
        izs = 0, idz = 0; float rzs = 0; double dzs = 0;
      if (n1qn1_ != NULL) {
        n1qn1_(gEdN1Cost, &nn, st.data(), &fN, gg.data(), var.data(),
               &eps, &mode, &niter, &nsim, &impr, zm.data(),
               &izs, &rzs, &dzs, &idz);
      }
      if (!gEdN1Bad && gEdN1Evals > 0 && gEdBest < f0) {
        for (int t = 0; t < nUni; ++t) {
          double xv = gEdBestPar[(size_t)t];
          if (!std::isfinite(xv)) continue;
          int c = col[(size_t)t];
          double cur = mprior_phi0(0, c);
          double v = cur + pas(kiter) * (xv - cur);
          if (std::isfinite(v)) { mprior_phi0.col(c).fill(v); moved = true; }
        }
        // rho, damped on the atanh scale like every other coordinate
        double an = gEdBestPar[(size_t)rhoIdx];
        if (std::isfinite(an)) {
          double a0 = std::atanh(rho0);
          double av = a0 + pas(kiter) * (an - a0);
          double rv = std::tanh(av);
          if (std::isfinite(rv)) {
            if (rv > 0.99) rv = 0.99; else if (rv < -0.99) rv = -0.99;
            // BOTH slots, and that is not belt-and-braces.
            //
            // `etaDistCorWith` is recorded only on the HIGHER-indexed member of
            // a pair, so that is the slot every reader looks in -- while this
            // step is entered from the LOWER member.  Writing only etaDistRho(k)
            // put the estimate somewhere nothing reads: the fit reported the
            // ini() value 0.3 while the sampler had moved on.
            if ((int)etaDistRho.n_elem == etaDistNdist) {
              etaDistRho(k) = rv;
              if (j >= 0 && j < etaDistNdist) etaDistRho(j) = rv;
            }
            // ...and to its own theta where one exists.  On the direct route it
            // does not (the correlation stays in the omega), which is what
            // .etaDistWarnCorFrozen() reports.
            int cc = corCol(k);
            if (cc >= 0 && cc < nphi0) {
              double aa = std::atanh(rv);
              if (std::isfinite(aa)) {
                mprior_phi0.col(cc).fill(aa);
                std::vector<int> one(1, cc);
                writeBackPhi0(one);
              }
            }
            etaDistCorFired = true;
            moved = true;
          }
        }
        if ((int)etaDistFiredK.size() == etaDistNdist) {
          etaDistFiredK[(size_t)k] = 1; etaDistFiredK[(size_t)j] = 1;
        }
      }
    }
    gEdRpn = NULL; gEdRpn2 = NULL; gEdEta = NULL; gEdEta2 = NULL; gEdWt = NULL;
    gEdRhoIdx = -1; gEdFam2 = -1;
    return moved;
  }

  // -------------------------------------------------------------------------
  // The DIRECT route's prior, on the eta scale.
  //
  // Everything here is a no-op unless rxEtaDistExpand(param="direct") built the
  // model, so the cdf route's arithmetic is untouched to the bit.
  // -------------------------------------------------------------------------

  // Per-declaration parsed argument expressions, for a declaration whose
  // arguments depend on a COVARIATE.
  //
  // etaDistArgs holds one population argument set per declaration, and for a
  // covariate declaration there is no such thing -- the R side leaves it NA on
  // purpose.  The M-step already recomputes per record; the PRIOR did not, so
  // every draw scored -Inf and a covariate on a declared eta could not be
  // sampled at all on this route.  Parsed once per iteration (pure arithmetic,
  // no solve) and reused across every row.
  mutable std::vector< std::vector< std::vector<etaDistTok> > > etaDistRpnCov;
  mutable std::vector<char> etaDistRpnCovOk;
  void etaDistBuildCovRpn() const {
    etaDistRpnCov.assign((size_t)(etaDistNdist > 0 ? etaDistNdist : 0),
                         std::vector< std::vector<etaDistTok> >());
    etaDistRpnCovOk.assign((size_t)(etaDistNdist > 0 ? etaDistNdist : 0), 0);
    for (int k = 0; k < etaDistNdist; ++k) {
      if (k >= (int)etaDistCov.size() || etaDistCov[(size_t)k].n_cols == 0) continue;
      int nth = etaDistNth(k);
      if (nth <= 0 || k >= (int)etaDistExprs.size()) continue;
      if ((int)etaDistExprThetas[(size_t)k].size() != nth) continue;
      std::vector<std::string> pv = etaDistExprThetas[(size_t)k];
      for (size_t c = 0; c < etaDistCovNames[(size_t)k].size(); ++c) {
        pv.push_back(etaDistCovNames[(size_t)k][c]);
      }
      std::vector< std::vector<etaDistTok> > rpn;
      if (rxEtaDistLoglikParse(etaDistExprs[(size_t)k], pv, rpn)) {
        etaDistRpnCov[(size_t)k] = rpn;
        etaDistRpnCovOk[(size_t)k] = 1;
      }
    }
  }

  // This declaration's arguments for the subject owning phiM row `r`.
  //
  // Falls back to the population set when there is no covariate, so every
  // caller can use one path.  Returns false when the arguments cannot be built,
  // which the callers treat the way they treat a non-finite density: reject,
  // never guess.
  // Does this declaration have covariate-dependent arguments at all?
  bool etaDistArgsHaveCov(int k) const {
    return k < (int)etaDistRpnCovOk.size() && etaDistRpnCovOk[(size_t)k] == 1 &&
      k < (int)etaDistCov.size() && etaDistCov[(size_t)k].n_cols > 0;
  }

  // The population argument set, which is what a declaration without a
  // covariate has and always had.
  bool etaDistArgsPop(int k, double *a, int na) const {
    for (int t = 0; t < na; ++t) {
      a[t] = etaDistArgs(k, t);
      if (!std::isfinite(a[t])) return false;
    }
    return true;
  }

  // The argument set for the SUBJECT owning phiM row `r`: the declaration's
  // thetas read from phi0, this subject's covariates, and the parsed argument
  // expressions evaluated on the two together.
  bool etaDistArgsCovRow(int k, unsigned int r, double *a, int na) const {
    const arma::mat &cv = etaDistCov[(size_t)k];
    unsigned int subj = (N > 0) ? (r % (unsigned int)N) : 0;
    if (subj >= cv.n_rows) return false;
    int nth = etaDistNth(k);
    int nSym = (int)cv.n_cols;
    std::vector<double> vals((size_t)(nth + nSym), 0.0);
    for (int t = 0; t < nth; ++t) {
      int c = etaDistPhi0Col(k, t);
      if (c < 0 || c >= nphi0) return false;
      vals[(size_t)t] = mprior_phi0(0, c);
      if (!std::isfinite(vals[(size_t)t])) return false;
    }
    for (int c = 0; c < nSym; ++c) {
      vals[(size_t)(nth + c)] = cv(subj, (unsigned int)c);
      if (!std::isfinite(vals[(size_t)(nth + c)])) return false;
    }
    const std::vector< std::vector<etaDistTok> > &rpn = etaDistRpnCov[(size_t)k];
    if ((int)rpn.size() != na) return false;
    for (int t = 0; t < na; ++t) {
      a[t] = etaDistExprEval(rpn[(size_t)t], vals.data(), nth + nSym);
      if (!std::isfinite(a[t])) return false;
    }
    return true;
  }

  // Does the model compute ANY of this declaration's arguments?
  bool etaDistAnchorOn(int k) const {
    if (k < 0 || k >= (int)etaDistAnchorIdx.size()) return false;
    for (size_t t = 0; t < etaDistAnchorIdx[(size_t)k].size(); ++t) {
      if (etaDistAnchorIdx[(size_t)k][t] >= 0) return true;
    }
    return false;
  }

  // This declaration's arguments for the subject owning phiM row `r`.
  //
  // Preference order, and each step is a fallback for exactly one reason:
  //
  //   1. the values HARVESTED from the solve -- the model computed them per
  //      observation, so a covariate is already in them;
  //   2. the population set, for an argument the model emits no line for.
  //
  // Reading the solve rather than re-evaluating the argument expressions is
  // what keeps one source of truth: the same arithmetic, in the pooled model,
  // done once.
  bool etaDistArgsFor(int k, unsigned int r, double *a, int na) const {
    // No anchors -- the model does not compute this declaration's arguments, or
    // their lhs indices did not resolve -- so use the expression evaluator.
    if (!etaDistAnchorOn(k)) {
      if (!etaDistArgsHaveCov(k)) return etaDistArgsPop(k, a, na);
      return etaDistArgsCovRow(k, r, a, na);
    }
    unsigned int subj = (N > 0) ? (r % (unsigned int)N) : 0;
    bool harvested = (subj < etaDistAnchorHave.size() && etaDistAnchorHave[subj]);
    if (!harvested) {
      // the harvest has not run for this subject yet this pass
      if (etaDistArgsHaveCov(k)) return etaDistArgsCovRow(k, r, a, na);
      return etaDistArgsPop(k, a, na);
    }
    // Seed from the population set WITHOUT requiring it to be finite, then let
    // the harvested values overwrite.
    //
    // Requiring it finite first is wrong and was measured to be: a
    // covariate-carrying argument has NO population value -- the R side leaves
    // it NA on purpose, and the trace shows `pop(shape,rate) = 1.25 nan`.
    // `etaDistArgsPop()` therefore failed, this returned false, the prior became
    // +Inf and the sampler degenerated: the gamma's etas went NEGATIVE at every
    // seed and bWT scattered (0.1234 / 0.4962 / 0.4425 against a truth of 0.75,
    // where the evaluator path gives 0.7441 / 0.7432 / 0.7618).  It looked for a
    // long time like "any second lhs in saem's model breaks the fit"; it was
    // this, tripped by the anchors resolving.
    for (int t = 0; t < na; ++t) a[t] = etaDistArgs(k, t);
    const std::vector<int> &ix = etaDistAnchorIdx[(size_t)k];
    for (int t = 0; t < na && t < (int)ix.size(); ++t) {
      if (ix[(size_t)t] < 0) continue;            // no line; keep population
      a[t] = etaDistAnchorVal(subj, (unsigned int)etaDistAnchorCol(k, t));
    }
    // now every argument must be usable, whatever supplied it
    for (int t = 0; t < na; ++t) if (!std::isfinite(a[t])) return false;
    return true;
  }

  // Flat column for (declaration, argument) in etaDistAnchorVal
  int etaDistAnchorCol(int k, int t) const {
    int na = (int)etaDistArgs.n_cols;
    return k*na + t;
  }

public:
  // Is this solve at the CURRENT STATE, or is it a candidate evaluation?
  //
  // saem solves many times per iteration and not all of them are the state: the
  // phi0 search (newuoa/bobyqa, reached through R) re-solves once per candidate,
  // with phi0 set to the candidate rather than to the fit's current value.
  // Harvesting indiscriminately captured those -- measured, a solve arrived with
  // phi0 = 2.29335, 3.65121, 0.358568, -1.36372 where the state was 1.38629,
  // 3.91202, -0.223144, 0.3, and `rxd.eta.cl` came through NEGATIVE, which a
  // gamma random effect cannot be.
  //
  // The test is the thing we actually care about rather than a proxy for it: a
  // state solve carries the state's phi0.  A candidate differs in at least the
  // coordinate being searched, so it is rejected by construction.
  bool etaDistAnchorAtState(const arma::mat &phi) const {
    if (nphi0 <= 0 || phi.n_rows == 0) return false;
    if ((int)i0.n_elem != nphi0 || (int)mprior_phi0.n_cols < nphi0) return false;
    for (int c = 0; c < nphi0; ++c) {
      unsigned int col = (unsigned int)i0(c);
      if (col >= phi.n_cols) return false;
      double a = phi(0, col), b = mprior_phi0(0, c);
      if (!std::isfinite(a) || !std::isfinite(b)) return false;
      if (fabs(a - b) > 1e-8*(1.0 + fabs(b))) return false;
    }
    return true;
  }

  // Start a fresh harvest.  Called once per solve pass: the anchors are
  // functions of the thetas and the subject's covariates, so they change when
  // the thetas move and are constant across the chains and across the MCMC
  // sweep that follows.
  void etaDistAnchorReset() const {
    int na = (int)etaDistArgs.n_cols;
    if (etaDistNdist <= 0 || etaDistAnchorIdx.empty() || N <= 0 || na <= 0) {
      etaDistAnchorHave.clear();
      return;
    }
    if ((int)etaDistAnchorVal.n_rows != N ||
        (int)etaDistAnchorVal.n_cols != etaDistNdist*na) {
      etaDistAnchorVal.set_size((unsigned int)N, (unsigned int)(etaDistNdist*na));
    }
    etaDistAnchorVal.fill(NA_REAL);
    etaDistAnchorHave.assign((size_t)N, 0);
    etaDistAnchorRows.assign((size_t)N, std::vector<double>());
    etaDistAnchorNrec.assign((size_t)N, 0);
    etaDistVariesK.clear();
  }

  // Take one individual's anchor values out of its lhs buffer.
  //
  // Only `id < N` -- the FIRST chain -- writes.  phiM stacks the chains, so
  // every subject appears nmc times and the anchors are identical across them
  // (they depend on the thetas and the covariates, not on eta); restricting to
  // one chain makes the writer unique, which matters because the pooled solve
  // runs this under OpenMP.
  // `nlhs` of the model actually being solved.  The harvest is called from more
  // than one solve path and they do NOT all solve the same model, so an index
  // resolved against saem's own model can point past a shorter model's lhs
  // buffer.  Measured: reading lhs[1] on such a path segfaulted.
  mutable int etaDistAnchorNlhs = 0;
  // One record's anchor values, appended to this subject's block.  The first
  // record also fills the per-subject row, which is what a declaration with no
  // within-subject variation reads.
  void etaDistAnchorAppend(int id, const double *lhs, std::vector<double> &rows,
                           size_t base, bool first) const {
    int na = (int)etaDistArgs.n_cols;
    for (int k = 0; k < etaDistNdist && k < (int)etaDistAnchorIdx.size(); ++k) {
      const std::vector<int> &ix = etaDistAnchorIdx[(size_t)k];
      for (int t = 0; t < na && t < (int)ix.size(); ++t) {
        if (ix[(size_t)t] < 0 || ix[(size_t)t] >= etaDistAnchorNlhs) continue;
        double v = lhs[ix[(size_t)t]];
        rows[base + (size_t)etaDistAnchorCol(k, t)] = v;
        if (first) {
          etaDistAnchorVal((unsigned int)id,
                           (unsigned int)etaDistAnchorCol(k, t)) = v;
        }
      }
    }
  }

  void etaDistAnchorHarvest(int id, const double *lhs) const {
    if (etaDistAnchorHave.empty() || lhs == NULL) return;
    if (id < 0 || id >= N || etaDistAnchorNlhs <= 0) return;
    int ncol = etaDistNdist*(int)etaDistArgs.n_cols;
    if (ncol <= 0) return;
    std::vector<double> &rows = etaDistAnchorRows[(size_t)id];
    size_t base = rows.size();
    rows.resize(base + (size_t)ncol, NA_REAL);
    etaDistAnchorAppend(id, lhs, rows, base, !etaDistAnchorHave[(size_t)id]);
    etaDistAnchorNrec[(size_t)id]++;
    etaDistAnchorHave[(size_t)id] = 1;
  }

  // Does this declaration's argument set vary WITHIN a subject?
  //
  // Decided from the harvested records rather than declared up front: a
  // covariate that is constant in the data behaves as constant whatever its
  // column could in principle do, and a declaration with no covariate at all
  // must not pay for the weighted path.  Computed once per harvest pass.
  // Does declaration kk's argument set differ across any subject's records?
  bool etaDistVariesOne(int kk, int na, int ncol) const {
    for (size_t sj = 0; sj < etaDistAnchorRows.size(); ++sj) {
      int nrec = etaDistAnchorNrec[sj];
      if (nrec <= 1) continue;
      const std::vector<double> &rows = etaDistAnchorRows[sj];
      for (int t = 0; t < na; ++t) {
        double v0 = rows[(size_t)etaDistAnchorCol(kk, t)];
        for (int r = 1; r < nrec; ++r) {
          double v = rows[(size_t)r*(size_t)ncol + (size_t)etaDistAnchorCol(kk, t)];
          if (v != v0 && (std::isfinite(v) || std::isfinite(v0))) return true;
        }
      }
    }
    return false;
  }

  // Does this declaration's argument set vary WITHIN a subject?
  //
  // Decided from the harvested records rather than declared up front: a
  // covariate that is constant in the data behaves as constant whatever its
  // column could in principle do, and a declaration with no covariate at all
  // must not pay for the weighted path.  Computed once per harvest pass.
  bool etaDistVaries(int k) const {
    if (k < 0 || k >= etaDistNdist) return false;
    int na = (int)etaDistArgs.n_cols;
    int ncol = etaDistNdist*na;
    if (ncol <= 0) return false;
    if ((int)etaDistVariesK.size() != etaDistNdist) {
      etaDistVariesK.assign((size_t)etaDistNdist, 2);
      for (int kk = 0; kk < etaDistNdist; ++kk) {
        if (etaDistVariesOne(kk, na, ncol)) etaDistVariesK[(size_t)kk] = 1;
      }
    }
    return etaDistVariesK[(size_t)k] == 1;
  }

  // This declaration's arguments at ONE harvested record of a subject.
  bool etaDistArgsAtRec(int k, unsigned int subj, int rec,
                        double *a, int na) const {
    if (subj >= etaDistAnchorRows.size()) return false;
    int ncol = etaDistNdist*na;
    const std::vector<double> &rows = etaDistAnchorRows[subj];
    if (rec < 0 || (size_t)(rec + 1)*(size_t)ncol > rows.size()) return false;
    for (int t = 0; t < na; ++t) a[t] = etaDistArgs(k, t);
    const std::vector<int> &ix = etaDistAnchorIdx[(size_t)k];
    for (int t = 0; t < na && t < (int)ix.size(); ++t) {
      if (ix[(size_t)t] < 0) continue;
      a[t] = rows[(size_t)rec*(size_t)ncol + (size_t)etaDistAnchorCol(k, t)];
    }
    for (int t = 0; t < na; ++t) if (!std::isfinite(a[t])) return false;
    return true;
  }

  int etaDistNrec(unsigned int subj) const {
    if (subj >= etaDistAnchorNrec.size()) return 0;
    return etaDistAnchorNrec[subj];
  }

private:

  bool etaDistDirectOn() const {
    return etaDistDirect != 0 && etaDistNdist > 0 &&
      (int)etaDistLatent.n_elem == etaDistNdist &&
      (int)etaDistFam.n_elem == etaDistNdist &&
      (int)etaDistArgs.n_rows == etaDistNdist;
  }

  // The declaration k is copula-paired with, and their rho, looked up in BOTH
  // directions.
  //
  // `etaDistCorWith` is ASYMMETRIC by construction: the R side records the pair
  // on the higher-indexed member only (`.cw[.hi] <- .lo - 1L`), leaving the
  // lower one at -1 with rho 0.  A loop that asks only `etaDistCorWith(k)`
  // therefore sees the lower member as unpaired -- and since columns are walked
  // in order, it scores that member's MARGINAL first and then scores the pair,
  // counting the lower marginal twice.
  int etaDistPartnerOf(int k, double *rho) const {
    *rho = 0.0;
    if (k < 0 || k >= etaDistNdist ||
        (int)etaDistCorWith.n_elem != etaDistNdist) return -1;
    int j = etaDistCorWith(k);
    if (j >= 0 && j < etaDistNdist) {
      if ((int)etaDistRho.n_elem == etaDistNdist) *rho = etaDistRho(k);
      return j;
    }
    // k is the LOWER member: find whoever names it
    for (int q = 0; q < etaDistNdist; ++q) {
      if (etaDistCorWith(q) == k) {
        if ((int)etaDistRho.n_elem == etaDistNdist) *rho = etaDistRho(q);
        return q;
      }
    }
    return -1;
  }

  // Which declaration, if any, owns each column of a sampled block.  `i` holds
  // the block's PHI columns; the answer is indexed by the block's own column
  // position, which is what phiM.cols(i) and dphi are indexed by.
  void etaDistDirectLocal(const uvec &i, std::vector<int> &declOf) const {
    declOf.assign(i.n_elem, -1);
    if (!etaDistDirectOn()) return;
    for (unsigned int c = 0; c < i.n_elem; ++c) {
      for (int k = 0; k < etaDistNdist; ++k) {
        if ((int)i(c) == etaDistLatent(k)) { declOf[c] = k; break; }
      }
    }
  }

  // Is any column of this block declared?  Cheap guard so a model that mixes a
  // declared eta with ordinary ones still pays nothing on the ordinary blocks.
  static bool etaDistAnyLocal(const std::vector<int> &declOf) {
    for (size_t c = 0; c < declOf.size(); ++c) if (declOf[c] >= 0) return true;
    return false;
  }

  // -log p(eta) for the declared columns of one block, per row.
  //
  // A correlated PAIR is scored ONCE, jointly, through the Gaussian copula --
  // scoring each marginal separately would drop the dependence entirely, and
  // scoring it twice would double-count both marginals.  `seen` is what keeps
  // that straight when both members are in the same block, which is the normal
  // case (they are correlated, so they share an omega block).
  //
  // A pair whose partner is NOT in this block falls back to its marginal: that
  // is the correct conditional up to a factor not depending on this column only
  // when rho == 0, so it is also flagged -- see etaDistDirectSplitPair.
  arma::vec etaDistDirectU(const mat &phiCols,
                           const std::vector<int> &declOf) const {
    arma::vec out(phiCols.n_rows, arma::fill::zeros);
    if (!etaDistDirectOn()) return out;
    const int na = (int)etaDistArgs.n_cols;
    std::vector<char> seen(declOf.size(), 0);
    // local column of each declaration, so a partner can be found
    std::vector<int> colOf((size_t)etaDistNdist, -1);
    for (size_t c = 0; c < declOf.size(); ++c) {
      if (declOf[c] >= 0) colOf[(size_t)declOf[c]] = (int)c;
    }
    std::vector<double> a1((size_t)na), a2((size_t)na);
    for (size_t c = 0; c < declOf.size(); ++c) {
      int k = declOf[c];
      if (k < 0 || seen[c]) continue;
      double rho = 0.0;
      int j = etaDistPartnerOf(k, &rho);
      int cj = (j >= 0 && j < etaDistNdist) ? colOf[(size_t)j] : -1;
      // The arguments are resolved PER ROW, because a covariate on a
      // declaration makes them subject-specific.  With no covariate
      // etaDistArgsFor() returns the population set and this is the same
      // arithmetic it always was.
      if (cj >= 0 && rho != 0.0) {
        seen[c] = 1; seen[(size_t)cj] = 1;
        etaDistDirectUPair(phiCols, k, (unsigned int)c, j, (unsigned int)cj,
                           rho, na, &a1[0], &a2[0], out);
      } else {
        seen[c] = 1;
        etaDistDirectUOne(phiCols, k, (unsigned int)c, na, &a1[0], out);
      }
    }
    return out;
  }

  // One declaration's contribution to the negative log prior, accumulated over
  // every row.  A row whose arguments cannot be built contributes +Inf: the
  // same treatment a non-finite density gets, so an unbuildable argument set is
  // rejected rather than guessed at.
  void etaDistDirectUOne(const mat &phiCols, int k, unsigned int c, int na,
                         double *a1, arma::vec &out) const {
    bool varies = etaDistVaries(k);
    for (unsigned int r = 0; r < phiCols.n_rows; ++r) {
      double l;
      if (varies) {
        l = etaDistWeightedLogD(k, r, phiCols(r, c), a1, na);
      } else {
        if (!etaDistArgsFor(k, r, a1, na)) {
          out(r) += std::numeric_limits<double>::infinity();
          continue;
        }
        l = rxEtaDistLogD(etaDistFam(k), phiCols(r, c), a1);
      }
      out(r) += R_finite(l) ? -l : std::numeric_limits<double>::infinity();
    }
  }

  // The WEIGHTED PER-OBSERVATION log density for a declaration whose arguments
  // vary within the subject.
  //
  //   sum_r wt_r * log p(eta_i | args_r),   wt_r = 1/n_i
  //
  // A within-subject-varying covariate gives the declaration no single argument
  // set, so there is no one density to evaluate: the first record's answer
  // would be an arbitrary choice among n_i of them.  Weighting by 1/n_i keeps
  // the result on the scale of ONE density, so a subject with more records does
  // not thereby get a sharper prior -- which is what `rxEtaDistLoglikObj`'s
  // `wt` argument has documented as `1/n_i` all along.
  //
  // Any record whose arguments cannot be built makes the whole subject
  // non-finite, the same treatment a single unbuildable set gets: rejected,
  // never averaged around.
  // The same weighting for a copula-linked PAIR, scored JOINTLY at each record,
  // so a covariate on either member moves the joint density for that record.
  double etaDistWeightedPairLogD(int k, int j, unsigned int r,
                                 double xk, double xj, double rho,
                                 double *a1, double *a2, int na) const {
    unsigned int subj = (N > 0) ? (r % (unsigned int)N) : 0;
    int nrec = etaDistNrec(subj);
    if (nrec <= 0) return std::numeric_limits<double>::quiet_NaN();
    double acc = 0.0, w = 1.0/(double)nrec;
    for (int rec = 0; rec < nrec; ++rec) {
      if (!etaDistArgsAtRec(k, subj, rec, a1, na) ||
          !etaDistArgsAtRec(j, subj, rec, a2, na)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      double lr = rxEtaDistPairLogD(etaDistFam(k), xk, a1,
                                    etaDistFam(j), xj, a2, rho);
      if (!R_finite(lr)) return std::numeric_limits<double>::quiet_NaN();
      acc += w*lr;
    }
    return acc;
  }

  double etaDistWeightedLogD(int k, unsigned int r, double x,
                             double *a, int na) const {
    unsigned int subj = (N > 0) ? (r % (unsigned int)N) : 0;
    int nrec = etaDistNrec(subj);
    if (nrec <= 0) {
      if (!etaDistArgsFor(k, r, a, na)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      return rxEtaDistLogD(etaDistFam(k), x, a);
    }
    double acc = 0.0, w = 1.0/(double)nrec;
    for (int rec = 0; rec < nrec; ++rec) {
      if (!etaDistArgsAtRec(k, subj, rec, a, na)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
      double l = rxEtaDistLogD(etaDistFam(k), x, a);
      if (!R_finite(l)) return std::numeric_limits<double>::quiet_NaN();
      acc += w*l;
    }
    return acc;
  }

  // A copula-linked PAIR's contribution, scored jointly.  Both members are
  // resolved for the same row, so a covariate on either moves the joint density
  // for that subject alone.
  void etaDistDirectUPair(const mat &phiCols, int k, unsigned int c,
                          int j, unsigned int cj, double rho, int na,
                          double *a1, double *a2, arma::vec &out) const {
    bool varies = etaDistVaries(k) || etaDistVaries(j);
    for (unsigned int r = 0; r < phiCols.n_rows; ++r) {
      double l;
      if (varies) {
        l = etaDistWeightedPairLogD(k, j, r, phiCols(r, c), phiCols(r, cj),
                                    rho, a1, a2, na);
      } else {
        if (!etaDistArgsFor(k, r, a1, na) || !etaDistArgsFor(j, r, a2, na)) {
          out(r) += std::numeric_limits<double>::infinity();
          continue;
        }
        l = rxEtaDistPairLogD(etaDistFam(k), phiCols(r, c), a1,
                              etaDistFam(j), phiCols(r, cj), a2, rho);
      }
      out(r) += R_finite(l) ? -l : std::numeric_limits<double>::infinity();
    }
  }

  // The Gaussian half, with the declared columns removed.
  //
  // Zeroing dphi's declared columns rather than slicing IGamma2_phi down to the
  // undeclared ones: those agree only when omega is block diagonal between the
  // two sets, because a sub-block of an INVERSE is not the inverse of the
  // sub-block.  rxEtaDistExpand(param="direct") refuses to build a model where
  // a declared eta shares an omega block with an ordinary one, so the two are
  // block diagonal wherever this runs and the cheap form is the exact one.
  arma::vec etaDistDirectSplitU(const mat &dphi, const mat &iG,
                                const std::vector<int> &declOf) const {
    mat d = dphi;
    for (size_t c = 0; c < declOf.size(); ++c) {
      if (declOf[c] >= 0) d.col((unsigned int)c).zeros();
    }
    return 0.5*sum(d % (d*iG), 1);
  }

  // The CURRENT point's prior for a block, which is what U_phi is seeded with
  // before each random-walk kernel.  One call site per block instead of the
  // quadratic written out six times: the seed and the candidate MUST use the
  // same prior, and six copies of one of them is how they drift apart.
  arma::vec etaDistSeedUphi(const mcmcphi &mphi, const mat &phiM) const {
    mat dphi = phiM.cols(mphi.i) - mphi.mprior_phiM;
    std::vector<int> declOf;
    etaDistDirectLocal(mphi.i, declOf);
    return etaDistDirectPriorU(phiM.cols(mphi.i), dphi, mphi.IGamma2_phi, declOf);
  }

  // The whole prior for a block: family for the declared columns, Gaussian for
  // the rest.  Identical to the plain quadratic when nothing is declared.
  arma::vec etaDistDirectPriorU(const mat &phiCols, const mat &dphi,
                                const mat &iG,
                                const std::vector<int> &declOf) const {
    if (!etaDistAnyLocal(declOf)) return 0.5*sum(dphi % (dphi*iG), 1);
    return etaDistDirectSplitU(dphi, iG, declOf) + etaDistDirectU(phiCols, declOf);
  }

  // Draw a block's declared columns FROM THE PRIOR, which is what makes kernel
  // 1 an independence sampler whose proposal density cancels out of the
  // acceptance ratio.  A correlated pair is drawn through the Gaussian copula:
  // two correlated standard normals decoded by each marginal's quantile
  // function, which is what the copula IS -- so kernel 1 on the direct route
  // and the cdf construction draw the same law, as they must.
  // `_saemUE` (the uninformed-eta mask) is deliberately NOT applied here.
  //
  // On the Gaussian path masking a coordinate leaves it at mprior_phi, which is
  // that coordinate's prior MEAN -- the right answer for a subject whose data
  // say nothing about it.  Here mprior is 0 and 0 is outside a gamma's support,
  // so the same masking would put the row at zero prior density and keep it
  // there.  A draw from the prior is what "uninformed" means; it is what the
  // Gaussian mask is approximating, and on this route it is available exactly.
  // `edZ` holds each row's normals, two per phi column, drawn from its row seed.
  void etaDistDirectDraw(mat &phiCols, const std::vector<int> &declOf, const mat &edZ) {
    if (!etaDistDirectOn() || !etaDistAnyLocal(declOf)) return;
    const int na = (int)etaDistArgs.n_cols;
    const unsigned int nr = phiCols.n_rows;
    std::vector<char> seen(declOf.size(), 0);
    std::vector<int> colOf((size_t)etaDistNdist, -1);
    for (size_t c = 0; c < declOf.size(); ++c) {
      if (declOf[c] >= 0) colOf[(size_t)declOf[c]] = (int)c;
    }
    std::vector<double> a1((size_t)na), a2((size_t)na);
    for (size_t c = 0; c < declOf.size(); ++c) {
      int k = declOf[c];
      if (k < 0 || seen[c]) continue;
      double rho = 0.0;
      int j = etaDistPartnerOf(k, &rho);
      int cj = (j >= 0 && j < etaDistNdist) ? colOf[(size_t)j] : -1;
      const arma::vec z1 = edZ.col((unsigned int)(2 * c));
      // per-row arguments, for the same reason etaDistDirectU() resolves them
      // per row: a covariate makes each subject's prior its own distribution.
      // A row whose arguments cannot be built is LEFT ALONE rather than drawn
      // from a fabricated one.
      if (cj >= 0 && rho != 0.0) {
        const arma::vec z2 = edZ.col((unsigned int)(2 * c + 1));
        double sr = std::sqrt(1.0 - rho*rho);
        seen[c] = 1; seen[(size_t)cj] = 1;
        for (unsigned int r = 0; r < nr; ++r) {
          if (!etaDistArgsFor(k, r, &a1[0], na) ||
              !etaDistArgsFor(j, r, &a2[0], na)) continue;
          double w2 = rho*z1(r) + sr*z2(r);
          phiCols(r, (unsigned int)c)  = rxEtaDistQ(etaDistFam(k), etaDistNormCdf(z1(r)), &a1[0]);
          phiCols(r, (unsigned int)cj) = rxEtaDistQ(etaDistFam(j), etaDistNormCdf(w2),    &a2[0]);
        }
      } else {
        seen[c] = 1;
        for (unsigned int r = 0; r < nr; ++r) {
          if (!etaDistArgsFor(k, r, &a1[0], na)) continue;
          phiCols(r, (unsigned int)c) = rxEtaDistQ(etaDistFam(k), etaDistNormCdf(z1(r)), &a1[0]);
        }
      }
    }
  }

  // A bijected random-walk step for the declared columns, and the log-Jacobian
  // the acceptance ratio then owes.
  //
  // A declared eta is usually bounded (gamma is positive), so an unbijected
  // Gaussian step proposes outside the support and the chain simply rejects --
  // which is not a bias but is an efficiency floor that gets worse the closer
  // the chain sits to the bound.  Walking on u = log(x - lo) instead proposes
  // inside the support always; the target in u carries |dx/du|, so the
  // acceptance needs log|dx/du|(new) - log|dx/du|(cur).  Returned rather than
  // folded in, because kernels 2 and 3 add it at different points.
  //
  // `scaleCol` is PER COLUMN and `lamRow` per row (iacceptPerId).  Per column
  // and not per sweep: kernel 2 moves every coordinate at once, so its outer
  // index k1 is the sweep counter and has nothing to do with which column this
  // is -- scaling a declared column by another column's step size is a silent
  // mis-tuning that only shows up as a bad acceptance rate.
  arma::vec etaDistDirectRwCols(mat &phiCols, const mat &phiCur,
                                const std::vector<int> &declOf,
                                const arma::mat &noise,
                                const arma::vec &scaleCol,
                                const arma::vec &lamRow) {
    arma::vec logJ(phiCols.n_rows, arma::fill::zeros);
    if (!etaDistDirectOn() || !etaDistAnyLocal(declOf)) return logJ;
    const int na = (int)etaDistArgs.n_cols;
    std::vector<double> a((size_t)na);
    for (size_t c = 0; c < declOf.size(); ++c) {
      int k = declOf[c];
      if (k < 0) continue;
      int fam = etaDistFam(k);
      for (unsigned int r = 0; r < phiCols.n_rows; ++r) {
        // per-row: the BIJECTOR depends on the support, and with a covariate on
        // a bound-carrying argument the support is this subject's own
        double xc = phiCur(r, (unsigned int)c);
        double xn = xc;
        double dj = 0.0;
        if (etaDistArgsFor(k, r, &a[0], na)) {
          dj = etaDistRwOne(fam, &a[0], xc,
                            scaleCol((unsigned int)c)*lamRow(r)*noise(r, (unsigned int)c),
                            &xn);
        }
        phiCols(r, (unsigned int)c) = xn;
        logJ(r) += dj;
      }
    }
    return logJ;
  }

  // One coordinate's random walk, taken in BIJECTOR space: map to u, step,
  // map back, and return the log-Jacobian difference the acceptance ratio
  // needs.  Any non-finite intermediate leaves the coordinate where it was and
  // contributes nothing -- a step that cannot be described is not taken, rather
  // than taken and then corrected.
  double etaDistRwOne(int fam, const double *a, double xc, double step,
                      double *xOut) const {
    *xOut = xc;
    double uc = rxEtaDistToU(fam, xc, a);
    if (!R_finite(uc)) return 0.0;
    double un = uc + step;
    double xn = rxEtaDistFromU(fam, un, a);
    if (!R_finite(xn)) return 0.0;
    *xOut = xn;
    double jn = rxEtaDistLogJac(fam, un, a);
    double jc = rxEtaDistLogJac(fam, uc, a);
    if (!R_finite(jn) || !R_finite(jc)) return 0.0;
    return jn - jc;
  }

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
               int mixIdx = 0,
               vec *rwScale = nullptr,
               vec *rwLam = nullptr,
               vec *rwScale3 = nullptr,
               vec *rwLam3 = nullptr) {
    mat fcMat;
    // Each kernel carries its OWN acceptance scale: they are different
    // proposals with different optimal acceptance rates, and one shared vector
    // adapted toward two targets is what froze the chain.
    vec *rwK = (method == 3 && rwScale3 != nullptr) ? rwScale3 : rwScale;
    // Per-subject random-walk scale (saemControl(iacceptPerId=)), replicated
    // across the nmc chains so it lines up with phiM's row layout.  phiM
    // stacks the chains vertically (mprior_phiM = repmat(., nmc, 1)), so row r
    // belongs to subject r % N.  All ones when the option is off, which makes
    // every expression below bit-identical to the pooled path.
    // kernel 3 uses its own per-subject scale (see adaptRwPerIdIter)
    vec *rwLamK = (method == 3 && rwLam3 != nullptr) ? rwLam3 : rwLam;
    const bool perId = (rwLamK != nullptr) && ((int)rwLamK->n_elem == N);
    vec lamRow;
    if (perId) lamRow = repmat(*rwLamK, nmc, 1);
    else lamRow = ones<vec>(mx.nM);
    vec fc, fs, Uc_y, Uc_phi, deltu;
    uvec ind;
    arma::vec accU(mx.nM);   // per-block threefry acceptance uniforms

    uvec i=mphi.i;
    // Which of this block's columns are declared random effects on the DIRECT
    // route, and therefore scored by their family rather than by omega.  All
    // -1, and every branch below inert, on the cdf route.
    std::vector<int> etaDistDecl;
    etaDistDirectLocal(i, etaDistDecl);
    const bool edDirect = etaDistDirectOn() && etaDistAnyLocal(etaDistDecl);
    // the log-Jacobian a bijected proposal owes the acceptance ratio
    arma::vec edLogJ;
    double double_xmin = 1.0e-200;                               //FIXME hard-coded xmin, also in neldermean.hpp
    double xmax = 1e300;
    for (int u=0; u<nu; u++)
      for (int k1=0; k1<mphi.nphi; k1++) {
        mat phiMc=phiM;
        edLogJ = arma::vec(mx.nM, arma::fill::zeros);
        // proposal noise and acceptance uniforms, drawn before the solve
        const uint64_t seedOff = _seedLayout.step(kiter, mixIdx > 0 ? mixIdx - 1 : 0,
                                                  mphi.block, method, u, k1);
        switch (method) {
        case 1: {
          mat noise(mx.nM, mphi.nphi);
          // the declared-family draw's normals, from the same row seeds
          mat edZ(mx.nM, edDirect ? 2 * mphi.nphi : 0);
          _saemDrawRows(saemSeed, seedOff, noise, accU, &edZ);
          phiMc.cols(i)=noise*mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          // The declared columns are OVERWRITTEN with a draw from the declared
          // family, so kernel 1 stays an independence sampler whose proposal is
          // the prior -- which is exactly what lets the acceptance below keep
          // using deltu = Uc_y - U_y with no prior term at all.  Left as the
          // Gaussian draw it would still be a valid proposal, but the ratio
          // would then need a correction that the shared code path does not
          // apply, and the chain would target the wrong distribution.
          if (edDirect) {
            mat pc = phiMc.cols(i);
            etaDistDirectDraw(pc, etaDistDecl, edZ);
            phiMc.cols(i) = pc;
          }
          break;
        }
        case 2: {
          mat noise(mx.nM, mphi.nphi); _saemDrawRows(saemSeed, seedOff, noise, accU);
          // rwOmega: NONMEM's Z = lambda*Omega (eq. 1.139).  Otherwise the
          // historical diagonal, which is also what saemix uses.
          mat step = rwOmega ? (noise * mphi.Gfull_phi)
                             : (noise * mphi.Gdiag_phi);
          // kernel 2's own acceptance scale, per coordinate
          if (rwScale != nullptr && rwScale->n_elem == (unsigned int)mphi.nphi) {
            for (int c = 0; c < mphi.nphi; ++c) step.col(c) *= (*rwScale)(c);
          }
          if (perId) step.each_col() %= lamRow;
          phiMc.cols(i)=phiM.cols(i) + step % current_saem_state->_saemUE.cols(i);
          // ...then REPLACE the declared columns with a bijected step, since an
          // additive one proposes outside a bounded support and is rejected on
          // sight.  The step size is the one this kernel already adapted, read
          // per row so iacceptPerId still applies.
          if (edDirect) {
            // the same per-coordinate step this kernel just used additively:
            // Gdiag is diagonal, so (noise*Gdiag)(r,c) = noise(r,c)*Gdiag(c,c),
            // and the bijected move reproduces its magnitude exactly.  rwOmega's
            // full-matrix step has no bijected analogue -- a correlated move
            // between a bounded and an unbounded coordinate is not one step --
            // so a declared column takes the diagonal step either way.
            arma::vec scCol((unsigned int)mphi.nphi);
            for (int c = 0; c < mphi.nphi; ++c) {
              double sc2 = (rwScale != nullptr &&
                            rwScale->n_elem == (unsigned int)mphi.nphi)
                ? (*rwScale)(c) : 1.0;
              scCol((unsigned int)c) = mphi.Gdiag_phi(c, c)*sc2;
            }
            mat pc = phiMc.cols(i);
            edLogJ += etaDistDirectRwCols(pc, phiM.cols(i), etaDistDecl, noise,
                                          scCol, lamRow);
            phiMc.cols(i) = pc;
          }
          break;
        }
        case 3: {
          mat noiseM(mx.nM, 1); _saemDrawRows(saemSeed, seedOff, noiseM, accU);
          vec noise = noiseM.col(0);
          double s3 = (rwScale3 != nullptr && rwScale3->n_elem == (unsigned int)mphi.nphi)
            ? (*rwScale3)(k1) : 1.0;
          vec step = noise*(mphi.Gdiag_phi(k1,k1) * s3);
          if (perId) step %= lamRow;
          phiMc.col(i(k1))=phiM.col(i(k1))+
            step % current_saem_state->_saemUE.col(i(k1));
          // kernel 3 moves ONE coordinate, so it only bijects when that
          // coordinate is the declared one
          if (edDirect && etaDistDecl[(size_t)k1] >= 0) {
            std::vector<int> one(etaDistDecl.size(), -1);
            one[(size_t)k1] = etaDistDecl[(size_t)k1];
            arma::mat n1(mx.nM, mphi.nphi, arma::fill::zeros);
            n1.col((unsigned int)k1) = noise;
            arma::vec scCol((unsigned int)mphi.nphi, arma::fill::zeros);
            scCol((unsigned int)k1) = mphi.Gdiag_phi(k1, k1)*s3;
            mat pc = phiMc.cols(i);
            edLogJ += etaDistDirectRwCols(pc, phiM.cols(i), one, n1, scCol, lamRow);
            phiMc.cols(i) = pc;
          }
          break;
        }
        case 4: {
          // NONMEM mode 1B: independence proposal from each subject's own
          // accumulated conditional mean/variance (see buildMode1B()).
          mat noise(mx.nM, mphi.nphi); _saemDrawRows(saemSeed, seedOff, noise, accU);
          mat step = m1bFull ? mode1BNoise(noise) : (noise % m1bSd);
          // UE masks the NOISE, not the mean -- matching mode 1, so a masked
          // coordinate stays at its proposal centre rather than collapsing to 0
          phiMc.cols(i) = m1bMean + (step % current_saem_state->_saemUE.cols(i));
          break;
        }
        }

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
          // proposal IS the prior, so the prior terms cancel out of the ratio
          deltu=Uc_y-U_y;
          // ...which holds only while the CURRENT point has positive prior
          // density.  On the direct route it need not: phiM is initialized on
          // the Gaussian scale, so a bounded family starts some rows outside
          // its own support, where the prior is zero and the cancelled ratio is
          // not the ratio at all.  Such a row then accepts on likelihood alone,
          // and if the likelihood does not object it stays outside forever.
          //
          // Measured on dunif(1, 12): the fitted etas ran to 0.6272, below the
          // declared lower bound.  Gamma hid this because a negative clearance
          // wrecks the likelihood; for a uniform an eta of 0.627 is a perfectly
          // good clearance, so nothing pushed it back in.
          //
          // A point with zero prior density must accept any candidate with
          // positive density -- the same Metropolis limit the U_y rescue below
          // rests on.
          if (edDirect) {
            mat dphiC = phiM.cols(i) - mphi.mprior_phiM;
            vec uCur = etaDistDirectPriorU(phiM.cols(i), dphiC,
                                           mphi.IGamma2_phi, etaDistDecl);
            mat dphiP = phiMc.cols(i) - mphi.mprior_phiM;
            vec uPro = etaDistDirectPriorU(phiMc.cols(i), dphiP,
                                           mphi.IGamma2_phi, etaDistDecl);
            for (unsigned int _q = 0; _q < deltu.n_elem && _q < uCur.n_elem; ++_q) {
              if (!std::isfinite(uCur(_q)) && std::isfinite(uPro(_q))) {
                deltu(_q) = R_NegInf;   // accept: the current point is impossible
              }
            }
          }
        }
        else if (method==4) {
          // Mode 1B is an INDEPENDENCE sampler whose proposal is not the
          // prior, so the Metropolis-Hastings ratio keeps both the prior term
          // AND a proposal-density correction:
          //   log alpha = [logpi(new)-logpi(cur)] + [log k(cur)-log k(new)]
          // With U = -logpi and Q(x) = 0.5*sum((x-m)^2/v) (the -0.5*sum(log v)
          // normalizer is identical for both points and cancels), this is
          //   deltu = (Uc_y-U_y) + (Uc_phi-U_phi) - (Q_new - Q_cur)
          // and the existing `deltu < -log(u)` test applies unchanged.
          mat dphic=phiMc.cols(i)-mphi.mprior_phiM;
          // the same prior split as the random-walk kernels; the mode-1B
          // proposal correction (Qn - Qc) is a separate term and unaffected
          Uc_phi = etaDistDirectPriorU(phiMc.cols(i), dphic,
                                       mphi.IGamma2_phi, etaDistDecl);
          vec Qn = mode1BQ(phiMc.cols(i));
          vec Qc = mode1BQ(phiM.cols(i));
          deltu=Uc_y-U_y+Uc_phi-U_phi-(Qn-Qc);
        }
        else {
          mat dphic=phiMc.cols(i)-mphi.mprior_phiM;
          if (edDirect) {
            // Prior on the DIRECT route: the family for the declared columns,
            // the Gaussian quadratic for the rest.  U_phi comes in already
            // computed the same way (see the etaDistDirectPriorU() calls that
            // seed it before each random-walk kernel), so both sides of the
            // difference are on the same footing.
            //
            // `- edLogJ` and not `+`: the proposal walks on u and the target in
            // u is p(x(u))|dx/du|, so U_u = -log p - log|dx/du|.  U here is a
            // NEGATIVE log density, hence the Jacobian enters with the opposite
            // sign to the one it has in the log-density form.
            Uc_phi = etaDistDirectPriorU(phiMc.cols(i), dphic,
                                         mphi.IGamma2_phi, etaDistDecl);
            deltu = Uc_y - U_y + Uc_phi - U_phi - edLogJ;
          } else {
            Uc_phi=0.5*sum(dphic%(dphic*mphi.IGamma2_phi),1);
            deltu=Uc_y-U_y+Uc_phi-U_phi;
          }
        }

        // Accept, with a rescue for a chain that has latched.
        //
        // deltu = Uc_y - U_y (+ the prior terms).  Once a row's CURRENT
        // objective U_y is non-finite, deltu is NaN, `deltu < threshold` is
        // false for every proposal, and that row can never accept again --
        // including from kernel 1, whose proposal has nothing to do with any
        // step size.  That is why a scale runaway showed up as acceptance
        // exactly 0.000 on ALL THREE kernels rather than just the random-walk
        // ones: the freeze is permanent and kernel-independent.
        //
        // A row in that state accepts any FINITE candidate, which is the
        // correct Metropolis limit (the current point has zero density) and
        // restores the chain instead of leaving it dead for the rest of the
        // fit.
        arma::uvec accHit = (deltu < -log(accU));
        for (unsigned int _q = 0; _q < U_y.n_elem && _q < accHit.n_elem; ++_q) {
          if (!std::isfinite(U_y(_q)) && std::isfinite(Uc_y(_q))) accHit(_q) = 1;
        }
        // The same rescue for a non-finite PRIOR, which only the direct route
        // can produce.  phiM is initialized on the Gaussian scale, so a
        // positive-support family (gamma) starts some rows at or below 0, where
        // -log p is +Inf; deltu is then NaN, `deltu < threshold` is false, and
        // that row cannot accept from a random-walk kernel at all.  Kernel 1
        // repairs it -- its proposal IS the prior, so it carries no prior term
        // and accepts on the likelihood alone -- but only on the next sweep,
        // and observed at the first seed some rows were still at phi = -1.39.
        //
        // A row whose current point has zero prior density accepts any
        // candidate with positive density.  That is the Metropolis limit, not a
        // fudge, and it is the same argument the U_y rescue above makes.
        if (edDirect && Uc_phi.n_elem == accHit.n_elem &&
            U_phi.n_elem == accHit.n_elem) {
          for (unsigned int _q = 0; _q < accHit.n_elem; ++_q) {
            if (!std::isfinite(U_phi(_q)) && std::isfinite(Uc_phi(_q)) &&
                std::isfinite(Uc_y(_q))) accHit(_q) = 1;
          }
        }
        ind = find(accHit);
        mcmcRecordAccept(method, ind, mx.nM);
        // acceptance-rate adaptation of the random-walk scale (methods 2/3).
        // Method 1 draws from the prior, so its acceptance is not a function
        // of any step size -- NONMEM does not adapt its mode 1 either.
        if (method > 1 && mx.nM > 0) {
          double target = (method == 2) ? iaccept       // multivariate: ~0.234
                                        : iacceptSingle; // one-at-a-time: ~0.44
          if (perId) {
            // NONMEM tunes lambda for EACH SUBJECT so that subject's own
            // acceptance rate approaches IACCEPT (eq. 1.139 and the text after
            // it).  The update itself happens ONCE PER ITERATION, in
            // mcmcCloseIter(), against acceptances pooled over every
            // random-walk block -- NOT here, per block.  Adapting per block
            // estimates a subject's rate from only `nmc` trials, which is far
            // too noisy to drive a multiplicative update: compounded over the
            // ~1+nphi blocks each iteration runs it random-walks lambda into
            // its own clamps, and a subject whose lambda hits the ceiling
            // proposes steps so large it never accepts again.  Measured on
            // Bauer's gamma model, the per-block version left 83.4% of
            // subjects frozen against 30.7% for the pooled default -- the
            // exact opposite of what the option is for.
            (void)target;
          } else if (method == 2) {
            adaptRw(rwScale, arma::regspace<uvec>(0, mphi.nphi - 1),
                    (double)ind.n_elem / (double)mx.nM, iaccept);
          } else {
            uvec one(1); one(0) = (arma::uword)k1;
            adaptRw(rwK, one, (double)ind.n_elem / (double)mx.nM, iacceptSingle);
          }
        }
        phiM(ind,i)=phiMc(ind,i);
        U_y(ind)=Uc_y(ind);
        if (method>1) {
          U_phi(ind)=Uc_phi(ind);
        }
        ind = getObsIdx(ix_idM.rows(ind));
        cur_fsave(ind)=fs(ind);
        // only kernel 3 walks the columns one at a time; every other kernel's
        // proposal already covers all of them, so one pass IS the sweep
        if (method!=3) {
          break;
        }
      }
    setRxThreadId(-1);
  }

  // MSAEM (Lavielle & Mbogning 2014): mixture-weighted (log-sum-exp) observation loss for a
  // candidate phi, evaluated under every component. mix() only lets the component affect which
  // observation columns are read, so only this loss needs mixture treatment; the prior penalty
  // (in do_mcmc_msaem()) stays a standard single-Gaussian quadratic form.
  // ---- saemix's compute.Uy ------------------------------------------------
  //
  // The observation -log-likelihood at a candidate phi, with every column the
  // caller does not overwrite held at its currently sampled value.  This is
  // saemix's compute.Uy (R/func_aux.R:355), and the quantity NONMEM's
  // technical guide differentiates in eqs. 1.47-1.52 for a theta that is not
  // reachable through mu.  It is the M-step objective for a theta whose
  // Omega^-1-weighted GLS update cannot move it -- see zeroOmegaDirectStep().
  //
  // THREAD SAFETY.  Declared const, and every buffer it writes is
  // function-local, so it CANNOT touch the _scratch_ members do_mcmc() and
  // mixObsLoss() share for this same computation -- the compiler enforces
  // that, which is the point of the const.  The per-observation transform is
  // an OpenMP region over pure elementwise math (_powerD/handleF); it touches
  // no rxode2 solve state, so unlike phi1Objective()'s region it deliberately
  // does NOT call setRxThreadId(), and it makes no R API call.  What is not
  // reentrant is the ODE solve inside user_fn(), which drives the
  // process-global _rx: the parallelism lives INSIDE this function, and
  // concurrent CALLS to it must be serialized by the caller.
  double computeUy(const mat &phiCand) const {
    const double double_xmin = 1.0e-200, xmax = 1e300;
    rx_solving_options *op = getSolvingOptions(_rx);
    int cores = getOpCores(op);
    // the transform loop is pure math, so thread-safety of the ODE method is
    // irrelevant to it; only the per-element work has to be worth splitting
    bool doParallel = (cores > 1) && (ntotal >= 2048);
    mat fcMat = nmRngGuard([&]{ return user_fn(phiCand, evt, optM); });
    const vec fc = fcMat.col(0);
    const vec curCens = fcMat.col(1);
    const vec curLimit = fcMat.col(2);
    double total = 0.0;
    if (distribution == 1) {
      vec yt;
      if (hasFixedObsTransform) {
        yt = yTrans;
      } else {
        yt = y;
        for (int i = ntotal; i--;) {
          int cur = (int)ix_endpnt(i);
          yt(i) = _powerD(y(i), lambda(cur), yj(cur), low(cur), hi(cur));
        }
      }
      for (int k = 0; k < nmc; k++) {
        const int obs0 = k * ntotal;
        const vec fsk = fc.subvec(obs0, obs0 + ntotal - 1);
        const vec limitk = curLimit.subvec(obs0, obs0 + ntotal - 1);
        const vec censk = curCens.subvec(obs0, obs0 + ntotal - 1);
        vec ft(ntotal), limT(ntotal), ftT(ntotal), g(ntotal);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) if(doParallel)
#endif
        for (int i = 0; i < ntotal; i++) {
          int cur = (int)ix_endpnt(i);
          limT(i) = _powerD(limitk(i), lambda(cur), yj(cur), low(cur), hi(cur));
          double fi = fsk(i);
          double fti = _powerD(fi, lambda(cur), yj(cur), low(cur), hi(cur));
          ft(i) = fti;
          ftT(i) = handleF((int)propT(cur), fti, fi, false, true);
        }
        saemFormG(g, vecares, vecbres, ftT, veccres, vecaddProp);
        g.elem(find(g == 0.0)).fill(1.0);
        g.elem(find(g < double_xmin)).fill(double_xmin);
        g.elem(find(g > xmax)).fill(xmax);
        vec ftAr, gAr;
        vec dyf = arDYFinto(yt, ft, g, ftAr, gAr);
        // same translation applyCensLoss() performs: doCensNormal1 speaks
        // log-likelihood/variance, the chain carries a loss and an SD, and an
        // uncensored row comes back untouched
        for (int i = ntotal; i--;) {
          dyf(i) = -doCensNormal1(censk(i), yt(i), limT(i), -dyf(i),
                                  ftAr(i), gAr(i)*gAr(i), 0);
        }
        total += accu(dyf);
      }
    } else if (distribution == 2) {
      for (int k = 0; k < nmc; k++) {
        const vec fck = fc.subvec(k*ntotal, (k+1)*ntotal - 1);
        total += accu(-y % log(fck) + fck);
      }
    } else if (distribution == 3) {
      for (int k = 0; k < nmc; k++) {
        const vec fck = fc.subvec(k*ntotal, (k+1)*ntotal - 1);
        total += accu(-y % log(fck) - (1 - y) % log(1 - fck));
      }
    } else if (distribution == 4) {
      // general log-likelihood: the prediction column IS the per-obs loglik
      total = -accu(fc);
    }
    if (!std::isfinite(total)) return 1e300;
    return total;
  }

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

  // Model-aware naive classification, used for mixProbMethod="regress" membership and for
  // MSAEM's stratified init (both run once before iteration 1): for each
  // subject/hypothesis, shift phi1's owned columns toward/away from that hypothesis
  // (pertSd BSV-SD) and evaluate the compiled model's actual fit.  Returns the
  // per-subject argmin-hypothesis classification (1-indexed).
  uvec mixNaiveClassify(double pertSd) {
    double double_xmin = 1.0e-200;
    double xmax = 1e300;
    mat lossByHyp(N, nMix, fill::zeros);
    // phi0 columns (fixed-effect-only parameters) carry a search variance of 1 in
    // phiM, not a real BSV, so a subject's draw sits e-fold off the population value
    // and swamps the between-component signal the classification is measuring.  Judge
    // every hypothesis at the population value instead (#1058); the phi1 (eta) columns
    // stay at their draws, which are genuine subject-level values.
    mat basePhi = phiM;
    if (nphi0 > 0) basePhi.cols(i0) = repmat(mprior_phi0, nmc, 1);
    for (int mHyp = 0; mHyp < nMix; mHyp++) {
      mat cand = basePhi;
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
                      int kiter,
                      vec *rwScale = nullptr,
                      vec *rwLam = nullptr,
                      vec *rwScale3 = nullptr,
                      vec *rwLam3 = nullptr) {
    mat phiMc;
    vec Uc_y, Uc_phi, deltu;
    uvec ind;
    uvec i = mphi.i;
    arma::vec accU(mx.nM);
    // see do_mcmc(): all ones, and so inert, when iacceptPerId is off
    // kernel 3 uses its own per-subject scale (see adaptRwPerIdIter)
    vec *rwLamK = (method == 3 && rwLam3 != nullptr) ? rwLam3 : rwLam;
    const bool perId = (rwLamK != nullptr) && ((int)rwLamK->n_elem == N);
    vec lamRow;
    if (perId) lamRow = repmat(*rwLamK, nmc, 1);
    else lamRow = ones<vec>(mx.nM);

    for (int u = 0; u < nu; u++)
      for (int k1 = 0; k1 < mphi.nphi; k1++) {
        phiMc = phiM;
        const uint64_t seedOff = _seedLayout.step(kiter, 0, mphi.block, method, u, k1);
        switch (method) {
        case 1: {
          mat noise(mx.nM, mphi.nphi); _saemDrawRows(saemSeed, seedOff, noise, accU);
          phiMc.cols(i) = noise * mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          break;
        }
        case 2: {
          mat noise(mx.nM, mphi.nphi); _saemDrawRows(saemSeed, seedOff, noise, accU);
          // rwOmega: NONMEM's Z = lambda*Omega (eq. 1.139).  Otherwise the
          // historical diagonal, which is also what saemix uses.
          mat step = rwOmega ? (noise * mphi.Gfull_phi)
                             : (noise * mphi.Gdiag_phi);
          if (rwScale != nullptr && rwScale->n_elem == (unsigned int)mphi.nphi) {
            for (int c = 0; c < mphi.nphi; ++c) step.col(c) *= (*rwScale)(c);
          }
          if (perId) step.each_col() %= lamRow;
          phiMc.cols(i) = phiM.cols(i) +
            step % current_saem_state->_saemUE.cols(i);
          break;
        }
        case 3: {
          mat noiseM(mx.nM, 1); _saemDrawRows(saemSeed, seedOff, noiseM, accU);
          vec noise = noiseM.col(0);
          double s3 = (rwScale3 != nullptr && rwScale3->n_elem == (unsigned int)mphi.nphi)
            ? (*rwScale3)(k1) : 1.0;
          vec step = noise * (mphi.Gdiag_phi(k1, k1) * s3);
          if (perId) step %= lamRow;
          phiMc.col(i(k1)) = phiM.col(i(k1)) +
            step % current_saem_state->_saemUE.col(i(k1));
          break;
        }
        }
        Uc_y = mixObsLoss(phiMc, mx);

        if (method == 1) {
          deltu = Uc_y - U_y;
        } else {
          mat dphic = phiMc.cols(i) - mphi.mprior_phiM;
          Uc_phi = 0.5 * sum(dphic % (dphic * mphi.IGamma2_phi), 1);
          deltu = Uc_y - U_y + Uc_phi - U_phi;
        }

        arma::uvec accHit = (deltu < -log(accU));   // see do_mcmc() for the rescue
        for (unsigned int _q = 0; _q < U_y.n_elem && _q < accHit.n_elem; ++_q) {
          if (!std::isfinite(U_y(_q)) && std::isfinite(Uc_y(_q))) accHit(_q) = 1;
        }
        ind = find(accHit);
        mcmcRecordAccept(method, ind, mx.nM);
        if (method > 1 && mx.nM > 0) {
          double accRate = (double)ind.n_elem / (double)mx.nM;
          double target = (method == 2) ? iaccept : iacceptSingle;
          if (perId) {
            (void)target;   // adapted once per iteration; see do_mcmc()
          } else if (method == 2) {
            // multidimensional symmetric random walk: optimal ~0.234
            adaptRw(rwScale, arma::regspace<uvec>(0, mphi.nphi - 1), accRate, iaccept);
          } else {
            // one-at-a-time (Metropolis-within-Gibbs): optimal ~0.44, and its
            // OWN scale vector -- see do_mcmc()
            uvec one(1); one(0) = (arma::uword)k1;
            adaptRw(rwScale3 != nullptr ? rwScale3 : rwScale, one, accRate, iacceptSingle);
          }
        }
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



// Shi-difference fallback trampoline: the population prediction vector at the
// candidate free phi0 values `t`.
//
// Solves through user_fn with _saemSolveCompleteOnce left at 0, so the solve
// goes to _saemOwnSolveSlot (odeSlotPred) rather than the sensitivity peer --
// the original model without theta gradients, which is the whole point of the
// fallback.  The caller has already cleared _saemFreezeOde and snapshotted
// phi0, and restores both afterwards.
static arma::vec gShiPredFn(arma::vec &t, int id) {
  (void)id;
  if (gShiSelf == nullptr) return arma::vec();
  return gShiSelf->shiPredAt(t, gShiFreeIx);
}

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

// `n` consecutive seed offsets starting at `first`.
static void seedLayoutPush(uint64_t first, uint64_t n, std::vector<double> &out) {
  for (uint64_t i = 0; i < n; ++i) out.push_back((double)(first + i));
}

// One phi block's MCMC step offsets, in the kernel's order.
static void seedLayoutPushBlock(const saemSeedLayout &L, int kiter, int comp, int block,
                                std::vector<double> &out) {
  const uint64_t f = kiter == 0 ? 20u : 1u;
  const uint64_t count[4] = {f * L.nu[0], f * L.nu[1], f * L.nu[2], L.nu1B};
  for (int m = 1; m <= 4; ++m) {
    const int nk1 = m == 3 ? (int)L.nphi[block] : 1;
    for (uint64_t u = 0; u < count[m - 1]; ++u) {
      for (int k1 = 0; k1 < nk1; ++k1) {
        seedLayoutPush(L.step(kiter, comp, block, m, (int)u, k1), L.nM, out);
      }
    }
  }
}

// One iteration's censored-value offsets, in the kernel's order.
static void seedLayoutPushCens(const saemSeedLayout &L, int kiter, std::vector<double> &out) {
  for (uint64_t c = 0; c < L.nComp; ++c) {
    for (uint64_t k = 0; k < L.nmc; ++k) {
      seedLayoutPush(L.cens(kiter, (int)c, (int)k), L.ntotal, out);
    }
  }
}

// Test hook: every draw's seed offset, in the kernel's draw order.
//[[Rcpp::export]]
Rcpp::NumericVector saemSeedLayoutTest_(Rcpp::IntegerVector nu, int nu1B, int nphi1,
                                        int nphi0, int nMix, int nM, int nmc,
                                        int ntotal, int niter) {
  saemSeedLayout L;
  for (int j = 0; j < 3; ++j) L.nu[j] = (uint64_t)nu[j];
  L.nu1B = (uint64_t)nu1B;
  L.nphi[0] = (uint64_t)nphi1;
  L.nphi[1] = (uint64_t)nphi0;
  L.nComp = (uint64_t)std::max(nMix, 1);
  L.nM = (uint64_t)nM;
  L.nmc = (uint64_t)nmc;
  L.ntotal = (uint64_t)ntotal;
  std::vector<double> out;
  for (int kiter = 0; kiter < niter; ++kiter) {
    for (int c = 0; c < (int)L.nComp; ++c) {
      seedLayoutPushBlock(L, kiter, c, 0, out);
      seedLayoutPushBlock(L, kiter, c, 1, out);
    }
    seedLayoutPushCens(L, kiter, out);
  }
  return Rcpp::wrap(out);
}

//[[Rcpp::export]]
long saemEtaDistN_() { return _saemEtaDistN; }

//[[Rcpp::export]]
int saemEtaDistOn_() { return _saemEtaDistOn; }


// Not Rcpp-exported: only saem_fit_ below reads it, and an export would mean
// regenerating RcppExports AND hand-editing src/init.c's .Call table for a
// value nothing in R asks for.
static int saemEtaDistObsLik_() { return _saemEtaDistObsLik; }

// Exposed for testing.  The copula M-step is one line of arithmetic that was
// wrong in a way no fit-level assertion would localize: the product-moment it
// used is the constrained MLE only when the draws have unit variance, so any
// departure -- which every burn-in has -- inflated it into its clamp, and a
// clamped correlation makes both declared etas share one latent.
// The copula update the M-steps share (rxEtaDistCorFromRz), exposed so the
// formula can be pinned directly instead of inferred from a fit.  `rho` is the
// current correlation, `rz` the CORRELATION of the raw latents.
//[[Rcpp::export]]
double rxEtaDistCorFromRz_(double rho, double rz) {
  return rxEtaDistCorFromRz(rho, rz);
}

//[[Rcpp::export]]
double rxEtaDistCorTest_(Rcpp::NumericVector z1, Rcpp::NumericVector z2) {
  std::vector<double> a(z1.begin(), z1.end()), b(z2.begin(), z2.end());
  return rxEtaDistCorMle(a, b);
}

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

static double gZeroOmObj(const double *p) {
  if (gZeroOmEvalN >= gZeroOmEvalMax) {
    // out of budget: report worse than anything seen so the simplex collapses
    // onto the best point instead of wandering
    return gZeroOmBestF + 1.0e10;
  }
  std::vector<double> q(gZeroOmIx.size());
  for (size_t i = 0; i < gZeroOmIx.size(); ++i) {
    double v = p[i];
    if (v < gZeroOmLo((arma::uword)i)) v = gZeroOmLo((arma::uword)i);
    if (v > gZeroOmHi((arma::uword)i)) v = gZeroOmHi((arma::uword)i);
    q[i] = v;
  }
  double f = gZeroOmSelf->zeroOmegaObjective(q.data());
  gZeroOmEvalN++;
  if (gZeroOmEvalN == 1 || f < gZeroOmBestF) {
    gZeroOmBestF = f;
    for (size_t i = 0; i < gZeroOmIx.size(); ++i) gZeroOmBest((arma::uword)i) = q[i];
  }
  return f;
}

static void gZeroOmNmFn(double *p, double *fx) { *fx = gZeroOmObj(p); }

// n1qn1's simul contract: ind 2|4 -> objective, 3|4 -> gradient.  Both come
// from ONE solve here, so the two branches share it rather than solving twice.
static void gPhi0N1Cost(int *ind, int *n, double *x, double *f, double *g,
                        int *ti, float *tr, double *td, int *id) {
  (void)ti; (void)tr; (void)td; (void)id;
  if (gPhi0Self == nullptr || gPhi0N1Bad) return;
  std::vector<double> gg((size_t)*n, 0.0);
  double fv = 0.0;
  vec pasDummy;
  if (!gPhi0Self->nonMuGradPhi0(0, pasDummy, x, &fv, gg.data())) {
    gPhi0N1Bad = 1;
    return;
  }
  gPhi0N1Evals++;
  if (*ind == 2 || *ind == 4) *f = fv;
  if (*ind == 3 || *ind == 4) for (int i = 0; i < *n; ++i) g[i] = gg[(size_t)i];
}

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
// True when the sensitivity peer has no ODE to integrate (a linCmt() model is
// solved analytically).  Relaxing solver tolerances cannot change such a solve,
// so the retry ladder is skipped for it -- exactly as focei skips it.
bool _saemThetaSensAnalytic = false;
int _saemEtaDistCppMap = 0;
int _saemEtaDistRMap = 0;
bool _saemThetaSensActive = false;
arma::ivec _saemThetaSensPhi0Col;
arma::ivec _saemThetaSensTheta;
int _saemNonMuGradEvery = 1;
int _saemThetaSensPredOffset = -1;
int _saemThetaSensROffset = -1;
arma::ivec _saemThetaSensThetaKind;
arma::ivec _saemThetaSensThetaCol;
arma::vec _saemThetaSensThetaFixedVal;
arma::ivec _saemThetaSensEtaCol;
int _saemThetaSensDvCol = -1;
int _saemThetaSensSensOffset = -1;
arma::ivec _saemThetaSensSensIx;
int _saemThetaSensNlhs = 0;

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
    for (int k = 0; k < nEta; ++k) {
      // A mu-referenced parameter's whole combined phi value went into its
      // THETA[] above, so its ETA[] is 0.  A nonMuEta has no THETA[] at all --
      // the parameter IS the eta -- so its value belongs here instead, or the
      // model is evaluated at a latent eta of zero for every subject.
      double v = 0.0;
      if (k < (int)_saemPhi1EtaNonMu.n_elem && _saemPhi1EtaNonMu(k) != 0) {
        v = _phi(i, _saemPhi1I1(_saemPhi1H2EtaCol(k)));
      }
      setIndParPtr(ind, nH2Theta + k, v);
    }
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
// Set for individual i when its solve threw (see saemNoThrow).  Read by
// saemReadRowsPooled, which turns it into the same 1e99/hasNan the NaN path
// produces -- the solve and the read are separate passes, so the flag has to
// outlive the solve loop.
static std::vector<int> _saemPooledThrew;

// Which slot SAEM's OWN per-iteration solve goes through, and where rx_pred_
// sits in it.
//
// When the theta-sensitivity peer is live this is odeSlotThetaSens -- the
// COMPLETE system.  That model emits rx_pred_ alongside its d(f)/d(theta)
// columns, so one solve serves both SAEM's likelihood read and the non-mu
// gradient.  Solving odeSlotPred for the likelihood and then odeSlotThetaSens
// for the gradient, as this used to, solved the same system twice at identical
// parameters and doubled the per-iteration cost (imp already avoids exactly
// this -- impThetaSensCollect's reuseSolve path, src/inner.cpp).
//
// -1 means no peer: the ordinary bulk par_solve path, unchanged.
int _saemOwnSolveSlot = -1;
int _saemOwnPredOffset = -1;
// One-shot request for the COMPLETE system (odeSlotThetaSens) on the NEXT solve
// only.  Making that the standing choice is a net loss: refinePhi0Lik's search
// spends nonMuThetaMaxEval solves per refinement on the objective, and
// upgrading all of them from the 2-state predNoLhs to the 10-state sensitivity
// model costs more than sharing one solve saves (measured on Bauer's gamma
// model, 70 iterations: 164.0s -> 259.8s).  So the complete system is requested
// for exactly the establish-states solve the gradient then reads, and every
// other solve keeps the cheap peer.
int _saemSolveCompleteOnce = 0;
// Which slot the LAST pooled solve actually used, so the read pass and the
// gradient follow the solve rather than a standing preference.
int _saemReadSlot = -1;

// Resolve the own-solve slot from the odeSwap REGISTRY rather than from
// whatever offsets happen to have been ingested yet -- setupRx runs before the
// control ingestion that fills _saemPhi1PredOffset, so reading that here would
// silently pick the wrong slot depending on call order.
static void saemPickOwnSolveSlot() {
  _saemOwnSolveSlot = -1;
  _saemOwnPredOffset = -1;
  // The STANDING slot is the cheap prediction peer.  The sensitivity model is
  // requested per call through _saemSolveCompleteOnce, for the single solve the
  // gradient reads -- see that flag for the measurement showing why making it
  // standing is slower, not faster.
  if (_saemPhi1PoolActive && odeSwapLoaded(odeSlotPred)) {
    int o = odeSwapLhsIndex(odeSlotPred, "rx_pred_");
    if (o >= 0) { _saemOwnSolveSlot = odeSlotPred; _saemOwnPredOffset = o; }
  }
}

static void saemSolveIndividualsPooled(int nInd, int slot) {
  rx_solving_options *op = getSolvingOptions(_rx);
  _saemReadSlot = slot;
  _saemPooledThrew.assign((size_t)nInd, 0);
  int cores = getOpCores(op);
  bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  // The peer being solved may carry no event sensitivities of its own, but the
  // ES shape is
  // still a process global: without a batch here, a shape some earlier solve
  // left installed (this fit's own odeSlotHess2 call, or a prior FOCEi fit's
  // fit-wide inner load) stays live and handle_evid injects jumps sized for
  // that OTHER model into this (smaller, compacted) solve's dosing events.
  // OdeSwapEsBatch's own constructor is what deactivates it for a "no ES" slot
  // -- see its contract in src/odeSwap.cpp.  Must be constructed outside the
  // OpenMP region below, matching every other peer solve's own batch.
  OdeSwapEsBatch predEsBatch(slot);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int i = 0; i < nInd; ++i) {
#ifdef _OPENMP
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    rx_solving_options_ind *ind = getSolvingOptionsInd(_rx, i);
    OdeSwapScope neqGuard(slot, ind, op);
    OdeSwapCmtScope cmtGuard(slot, op, ind);
    setIndSolve(ind, -1);
    if (!saemNoThrow([&]{ odeSwapSolveInd(slot, i); })) {
      _saemPooledThrew[(size_t)i] = 1;
    }
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
    OdeSwapScope neqGuard(_saemReadSlot, ind, op);
    OdeSwapCmtScope cmtGuard(_saemReadSlot, op, ind);
    int nAll = getIndNallTimes(ind);
    std::vector<double> &obs = rowObs[(size_t)i];
    obs.reserve((size_t)nAll);
    bool threw = (i < (int)_saemPooledThrew.size()) && _saemPooledThrew[(size_t)i];
    if (threw || odeSwapIndBadSolveSlot(op, ind, _saemReadSlot)) {
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
      if (!saemNoThrow([&]{
            if (_saemReadSlot == odeSlotThetaSens) {
              rxThetaSens.calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
            } else {
              rxPred.calc_lhs(i, curT, getOpIndSolve(op, ind, j), lhs);
            } })) {
        obs.push_back(1.0e99);
        rowNan[i] = 1;
        continue;
      }
      // Only i < N writes (the first chain), so the writer is unique per
      // subject and this is safe under the OpenMP loop above.
      if (gAnchorSelf != nullptr && gAnchorAtState) {
        gAnchorSelf->etaDistAnchorHarvest(i, lhs);
      }
      double cur = lhs[(_saemReadSlot == odeSlotThetaSens)
                       ? _saemThetaSensPredOffset : _saemOwnPredOffset];
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
  int _slotNow = _saemOwnSolveSlot;
  if (_saemSolveCompleteOnce && _saemThetaSensActive &&
      odeSwapLoaded(odeSlotThetaSens)) {
    _slotNow = odeSlotThetaSens;
  }
  if (_slotNow >= 0) {
    saemSolveIndividualsPooled(_Nnlmixr2, _slotNow);
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
    if (_slotNow >= 0) {
      saemSolveIndividualsPooled(_Nnlmixr2, _slotNow);
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
  // The declared-distribution argument anchors are ordinary model lhs, computed
  // per observation by the compiled model.  Harvest them from this pass rather
  // than evaluating the argument expressions again elsewhere.
  // OFF by default.  Two reasons, both measured, and both have to be fixed
  // before this can be switched on:
  //
  //   * saem's model does not receive the declaration's thetas correctly.  With
  //     `rxEdA.eta.cl.shape=exp(-lclrv)` emitted, the model computed 0.698676
  //     where lclrv = -0.223144 gives exp(0.223144) = 1.25.  Nothing read those
  //     params before, so the defect was invisible.
  //   * `gAnchorSelf` is a file-static set during one fit's control ingestion
  //     and is not cleared when that SAEM object dies, so a later fit on a
  //     model with no declaration used a dangling pointer -- measured as
  //     "double free or corruption (out)" part-way through a test file, while a
  //     single-fit script was fine.
  bool _edAtState = (gAnchorSelf != nullptr) && gAnchorSelf->etaDistAnchorAtState(_phi);
  gAnchorAtState = _edAtState;
  if (gAnchorSelf != nullptr && _edAtState) {
    gAnchorSelf->etaDistAnchorReset();
  }
  if (_saemOwnSolveSlot >= 0) {
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
        if (gAnchorSelf != nullptr && _edAtState) {
          gAnchorSelf->etaDistAnchorHarvest((int)id, lhs);
        }
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
  // The theta-sensitivity peer for the non-mu (phi0) gradient refinement.
  // Declared for ANY model shape that produced one -- unlike the phi1 peers
  // below it is not tied to a general-likelihood fit.  Declaring it here also
  // means odeSwapPlan() sizes the shared pool for it, so the solve is threaded
  // across subjects like every other pooled solve.
  _saemThetaSensActive = opt.containsElementNamed("saemThetaSens") &&
    !Rf_isNull(opt["saemThetaSens"]);
  if (_saemThetaSensActive) {
    // Declared here, ONCE, ahead of either branch's sizing solve -- both the
    // pooled solve below and the ordinary one further down need odeSwapPlan()
    // to have already seen this peer.  Registration (which rxDynLoad's) waits
    // until after that solve, in whichever branch runs.
    if (!odeSwapDeclare(odeSlotThetaSens, "thetaSens", opt["saemThetaSens"])) {
      _saemThetaSensActive = false;
    }
  }
  if (_saemThetaSensActive) {
    _saemThetaSensPhi0Col = as<arma::ivec>(opt["saemThetaSensPhi0Col"]);
    _saemThetaSensTheta = as<arma::ivec>(opt["saemThetaSensTheta"]);
    _saemThetaSensThetaKind = as<arma::ivec>(opt["saemThetaSensThetaKind"]);
    _saemThetaSensThetaCol = as<arma::ivec>(opt["saemThetaSensThetaCol"]);
    _saemThetaSensThetaFixedVal = as<arma::vec>(opt["saemThetaSensThetaFixedVal"]);
    _saemThetaSensEtaCol = as<arma::ivec>(opt["saemThetaSensEtaCol"]);
    _saemThetaSensDvCol = as<int>(opt["saemThetaSensDvCol"]);
    _saemNonMuGradEvery = opt.containsElementNamed("nonMuThetaGradEvery") ?
      as<int>(opt["nonMuThetaGradEvery"]) : 1;
    if (_saemNonMuGradEvery < 1) _saemNonMuGradEvery = 1;
  } else {
    _saemThetaSensPhi0Col.reset();
    _saemThetaSensTheta.reset();
    _saemThetaSensThetaKind.reset();
    _saemThetaSensThetaCol.reset();
    _saemThetaSensThetaFixedVal.reset();
    _saemThetaSensEtaCol.reset();
    _saemThetaSensDvCol = -1;
  }
  _saemPhi1PoolActive = opt.containsElementNamed("saemPhi1Pred") &&
    !Rf_isNull(opt["saemPhi1Pred"]);
  // The pool is built ONCE, sized for whichever declared peer has the most ODE
  // states, and every peer must be declared before that sizing solve.  So this
  // block runs when EITHER family of peers is present: the phi1 pair (a
  // general-likelihood fit) or the theta-sensitivity model (any fit that wants
  // the exact-gradient non-mu refinement), or both.
  //
  // Registering before the pool exists rebinds rxode2's event-sensitivity
  // globals and corrupts the solve, and a stale registry from a PRIOR fit has
  // already been seen to segfault a later one -- hence odeSwapClearAll() above
  // and the strict declare -> size -> register order below.
  if (_saemPhi1PoolActive) {
    bool haveHess2 = _saemPhi1PoolActive &&
      opt.containsElementNamed("saemPhi1Hess2") &&
      !Rf_isNull(opt["saemPhi1Hess2"]);
    if (haveHess2) odeSwapDeclare(odeSlotHess2, "hess2", opt["saemPhi1Hess2"]);
    if (_saemPhi1PoolActive) odeSwapDeclare(odeSlotPred, "pred", opt["saemPhi1Pred"]);
    // the sensitivity peer is declared at the top of setupRx, ahead of this
    // sizing solve, so odeSwapPlan() already accounts for its neq and lhs
    // rxSolve_ on whichever peer has the most states (innerHess2's extra
    // eta-sensitivity states when it built, else predNoLhs) -- matches the
    // number-of-ODEs sizing rule odeSwap already uses for neq; both share
    // the identical THETA[]/ETA[]/DV parameter declaration (verified by
    // .saemPhi1TargetMap), so either works as the params-matrix source.
    // Widest = most ODE states, which is odeSwap's own neq sizing rule.  With
    // the sensitivity peer in play the answer is no longer "hess2 if it built":
    // its extra d(state)/d(theta) equations can exceed the eta-sensitivity ones.
    RObject widePar = R_NilValue;
    {
      int wideN = -1;
      auto consider = [&](RObject cand) {
        if (Rf_isNull(cand)) return;
        List mvC = _rxode2_rxModelVars_(cand);
        CharacterVector st = mvC[RxMv_state];
        if ((int)st.size() > wideN) { wideN = (int)st.size(); widePar = cand; }
      };
      if (haveHess2) consider(opt["saemPhi1Hess2"]);
      if (_saemPhi1PoolActive) consider(opt["saemPhi1Pred"]);
      if (_saemThetaSensActive) consider(opt["saemThetaSens"]);
    }
    if (Rf_isNull(widePar)) {
      // nothing usable to size against; fall through to the ordinary
      // (non-pooled) setup rather than solving a null model
      _saemThetaSensActive = false;
      _saemPhi1PoolActive = false;
      odeSwapClearAll();
    } else {
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
    if (_saemPhi1PoolActive) odeSwapRegister(odeSlotPred, "pred", opt["saemPhi1Pred"], &rxPred);
    if (_saemThetaSensActive) {
      if (!odeSwapRegister(odeSlotThetaSens, "thetaSens", opt["saemThetaSens"],
                           &rxThetaSens)) {
        // could not bind it -- run without the gradient step rather than solve
        // a slot whose entry points were never resolved
        _saemThetaSensActive = false;
      } else {
        _saemThetaSensPredOffset = odeSwapLhsIndex(odeSlotThetaSens, "rx_pred_");
        _saemThetaSensROffset = odeSwapLhsIndex(odeSlotThetaSens, "rx_r_");
        _saemThetaSensNlhs = odeSwapNlhs(odeSlotThetaSens);
        // The peer emits one rx__sens_rx_pred__BY_THETA_j___ per estimated
        // theta, ascending and contiguous, so the first one's lhs index plus
        // the output's position is the whole map (imp reads them the same way,
        // src/inner.cpp).
        _saemThetaSensSensOffset = -1;
        _saemThetaSensSensIx.set_size(_saemThetaSensTheta.n_elem);
        _saemThetaSensSensIx.fill(-1);
        for (unsigned int q = 0; q < _saemThetaSensTheta.n_elem; ++q) {
          std::string fq = "rx__sens_rx_pred__BY_THETA_" +
            std::to_string(_saemThetaSensTheta(q)) + "___";
          int ix = odeSwapLhsIndex(odeSlotThetaSens, fq.c_str());
          _saemThetaSensSensIx(q) = ix;
          if (q == 0) _saemThetaSensSensOffset = ix;
          if (getenv("NLMIXR2_SAEM_GRADCHECK") != NULL)
            Rprintf("  sens out %u: %s -> lhs %d\n", q, fq.c_str(), ix);
        }
        if (_saemThetaSensPredOffset < 0 || _saemThetaSensSensOffset < 0 ||
            _saemThetaSensNlhs <= 0) _saemThetaSensActive = false;
      }
    }
    }
    saemPickOwnSolveSlot();
    return;
  }

  saemPickOwnSolveSlot();
  rxUpdateFuns(mv["trans"], &rxInner);
  // Non-pooled fit (a plain normal model, which is most of them): carry the
  // theta-sensitivity peer here instead.  Declared BEFORE the sizing solve so
  // odeSwapPlan() accounts for it, and only registered after -- rxDynLoad-ing a
  // sensitivity model before rxSolve_ has built the pool rebinds rxode2's
  // event-sensitivity globals and the inner solve then frees a buffer sized for
  // the wrong neq.
  //
  // The pool is still sized by SAEM's OWN model here, not the peer.  If the
  // peer needs more states it is refused at solve time with
  // odeDenyPoolNotSized -- a loud refusal, not a corrupted solve -- and the
  // refinement falls back to the search alone.  That keeps the blast radius of
  // this on ordinary fits at zero, which matters more than covering every model
  // shape on the first pass.
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
    if (_saemThetaSensActive) {
      // Which model can size this pool is decided by MEASUREMENT, not by
      // assuming SAEM's own model is the biggest.  It usually is not: the
      // sensitivity peer expands linCmt() into explicit compartments and adds a
      // d(state)/d(theta) equation per estimated theta.  Measured on Bauer's
      // gamma model:
      //
      //     SAEM own   neq=2   nlhs=1   npars=9
      //     thetaSens  neq=10  nlhs=14  npars=12
      //
      // The peer is larger on every axis, so this pool cannot hold it and
      // solving it here reads and writes past the per-thread slice.
      // odeDenyPoolNotSized does NOT catch that (measured: it segfaulted in the
      // accumulate loop).
      //
      // Sizing the pool by the PEER instead -- solving it here the way the
      // pooled branch above solves its widest peer -- does not work on this
      // path either, and the reason is worth recording so it is not retried:
      // rxSolve_ lays out ONE parameter vector, the solved model's.  The peer
      // declares THETA[k]/ETA[k]; SAEM's own model declares native names
      // (lclm, rxz.eta.cl, ...).  Registering SAEM's own model as a peer and
      // switching to it then makes it read the peer's placeholder slots as its
      // own parameters -- measured: an immediate segfault.  The pooled branch
      // gets away with switching precisely because ALL of its peers share one
      // THETA[]/ETA[]/DV declaration.
      //
      // So the real fix is to give SAEM's own likelihood read that shared
      // declaration on this path too (what predNoLhs already does for a
      // general-likelihood fit), not to reorder the sizing.  Until then,
      // measure and decline: the refinement runs the derivative-free search
      // alone, exactly as it did before any of this.
      List mvTs = _rxode2_rxModelVars_(opt["saemThetaSens"]);
      CharacterVector stTs = mvTs[RxMv_state], stOwn = mv[RxMv_state];
      CharacterVector lhTs = mvTs[RxMv_lhs], lhOwn = mv[RxMv_lhs];
      CharacterVector pTs = mvTs[RxMv_params], pOwn = mv[RxMv_params];
      {
        IntegerVector fl = mvTs[RxMv_flags];
        _saemThetaSensAnalytic = (fl.size() > RxMvFlag_linCmt) &&
          (fl[RxMvFlag_linCmt] != 0);
      }
      if (getenv("NLMIXR2_SAEM_GRADCHECK") != NULL) {
        Rprintf("pool sizing (measured): own neq=%d nlhs=%d npars=%d | "
                "thetaSens neq=%d nlhs=%d npars=%d\n",
                (int)stOwn.size(), (int)lhOwn.size(), (int)pOwn.size(),
                (int)stTs.size(), (int)lhTs.size(), (int)pTs.size());
      }
      if (stTs.size() > stOwn.size() || lhTs.size() > lhOwn.size()) {
        _saemThetaSensActive = false;
      }
    }
    if (_saemThetaSensActive) {
      if (!odeSwapRegister(odeSlotThetaSens, "thetaSens", opt["saemThetaSens"],
                           &rxThetaSens)) {
        _saemThetaSensActive = false;
      } else {
        _saemThetaSensPredOffset = odeSwapLhsIndex(odeSlotThetaSens, "rx_pred_");
        _saemThetaSensROffset = odeSwapLhsIndex(odeSlotThetaSens, "rx_r_");
        _saemThetaSensNlhs = odeSwapNlhs(odeSlotThetaSens);
        // The peer emits one rx__sens_rx_pred__BY_THETA_j___ per estimated
        // theta, ascending and contiguous, so the first one's lhs index plus
        // the output's position is the whole map (imp reads them the same way,
        // src/inner.cpp).
        _saemThetaSensSensOffset = -1;
        _saemThetaSensSensIx.set_size(_saemThetaSensTheta.n_elem);
        _saemThetaSensSensIx.fill(-1);
        for (unsigned int q = 0; q < _saemThetaSensTheta.n_elem; ++q) {
          std::string fq = "rx__sens_rx_pred__BY_THETA_" +
            std::to_string(_saemThetaSensTheta(q)) + "___";
          int ix = odeSwapLhsIndex(odeSlotThetaSens, fq.c_str());
          _saemThetaSensSensIx(q) = ix;
          if (q == 0) _saemThetaSensSensOffset = ix;
          if (getenv("NLMIXR2_SAEM_GRADCHECK") != NULL)
            Rprintf("  sens out %u: %s -> lhs %d\n", q, fq.c_str(), ix);
        }
        if (_saemThetaSensPredOffset < 0 || _saemThetaSensSensOffset < 0 ||
            _saemThetaSensNlhs <= 0) _saemThetaSensActive = false;
      }
    }
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
  _saemEtaDistCorNotEst = saem.get_etaDistCorNotEst();
  // Rf_warningcall() is the C++ route onto the fit's $runInfo -- collected in
  // nlmixr2Est.R and printed under "Information about run", the same path
  // src/npde.cpp uses.  Worth a warning rather than a note: a correlation the
  // data do not identify still comes back with a number beside it, and that
  // number is whatever the search drifted to.  All the information about a
  // copula correlation is in the departure of S_z from the identity; when the
  // off-diagonal sits inside the sampling noise of zero there is none.
  if (_saemEtaDistCorNotEst > 0) {
    Rf_warningcall(R_NilValue,
                   "saem: %d declared copula correlation(s) are not identified by "
                   "the data -- the latent second-moment matrix is within sampling "
                   "noise of the identity, so the reported correlation reflects the "
                   "starting value and the search rather than information in the "
                   "data.  Treat it as unestimated.",
                   _saemEtaDistCorNotEst);
  }

  // etaDistMstep=TRUE with the step never firing is indistinguishable from the
  // option being off -- the fit converges and looks entirely normal.  The
  // spread guard can legitimately hold it back for a whole fit (a chain that
  // never settles below the unit prior), so say so rather than let it pass as a
  // silent no-op.  Same report imp makes for the same reason.
  // Not in the observation-likelihood mode: there the family M-step standing
  // down for the declared families is the POINT, not a no-op to report.  The
  // thetas were estimated -- by refinePhi0Lik against the observation
  // likelihood -- so saying "never ran ... had no effect" would be wrong twice.
  if (saemEtaDistOn_() && saemEtaDistN_() == 0 && !saemEtaDistObsLik_()) {
    int _mf = saem.get_etaDistMapFail();
    if (_mf > 0) {
      // Distinguish the two reasons.  This one is not a chain that needs longer
      // to settle -- it is a declaration this M-step cannot represent, and it
      // will not fix itself with more iterations.
      RSprintf("saem: the declared-distribution M-step never ran -- its fitted family parameters could not be mapped back to thetas (%d attempts).\n", _mf);
      RSprintf("      This is what happens when a distribution parameter depends on a COVARIATE: the map solves for one\n");
      RSprintf("      population-level set of native parameters, and a covariate gives every subject their own.  The family\n");
      RSprintf("      parameters were estimated by the rest of saem, not by this step; set etaDistMstep=FALSE to silence this.\n");
    } else {
      RSprintf("saem: the declared-distribution M-step never ran (etaDistMstep had no effect).\n");
      RSprintf("      The pooled latent spread never settled: it kept changing by more than etaDistSdTol\n");
      RSprintf("      between attempts, or left the [etaDistSdLo, etaDistSdHi] divergence cap.  The usual\n");
      RSprintf("      cause is too short a run -- the step needs at least two attempts after the chain\n");
      RSprintf("      has equilibrated, so raise nBurn or lower etaDistEvery.\n");
    }
  }

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
    Named("mixWeights") = wrap(saem.get_mixWeights()),
    Named("mcmcAccept") = saem.get_mcmcAccTrace(),
    Named("mcmcStuck") = saem.get_mcmcStuckTrace(),
    Named("mcmcPhiSd") = saem.get_phiSdTrace(),
    Named("mcmcPhiAcf") = saem.get_phiAcfTrace(),
    Named("etaDistRho") = wrap(saem.get_etaDistRho()),
    Named("etaDistCorWith") = wrap(saem.get_etaDistCorWith())
  );
  current_saem_state = nullptr;
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
}

// Expose the declared-distribution dispatch so every family's quantile and log
// density can be pinned against its R counterpart (the same reason
// saemFormGTest exists).  22 hand-written densities are exactly the kind of
// thing that is silently wrong otherwise.
//[[Rcpp::export]]
SEXP rxEtaDistTest_(int fam, SEXP inU, SEXP inArgs) {
  NumericVector u(inU), ar(inArgs);
  int na = rxEtaDistNarg(fam);
  if (na < 0) return R_NilValue;
  double a[4] = {0,0,0,0};
  for (int i = 0; i < na && i < ar.size(); ++i) a[i] = ar[i];
  NumericVector q(u.size()), ld(u.size());
  for (int i = 0; i < u.size(); ++i) {
    q[i] = rxEtaDistQ(fam, u[i], a);
    ld[i] = rxEtaDistLogD(fam, q[i], a);
  }
  return List::create(_["q"] = q, _["logd"] = ld, _["narg"] = na);
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
