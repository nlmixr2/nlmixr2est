#define STRICT_R_HEADER
#include <stdio.h>
#include <stdarg.h>
#include <thread>
#include <chrono>
#include <vector>
#include <R_ext/Rdynload.h>
#include <RcppArmadillo.h>
#include <rxode2ptr.h>
#include "utilc.h"
#include "censEst.h"
#include "nearPD.h"
#include "inner.h"

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
      g = ab0*ab0 + ab1*ab1*pow(fa, 2*pw);
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


// class def starts
class SAEM {
  typedef mat (*user_funct) (const mat&, const mat&, const List&);

public:

  SAEM() {
    user_fn = NULL;
  }

  ~SAEM() {}

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
  // (reduces to the independent 0.5*((yt-ft)/g)^2 + log(g) when no AR).
  vec arDYFhyp(const vec &yt, const vec &ft, const vec &g) {
    vec e = yt - ft;
    vec gg = g;
    arWhiten(e, gg);
    return 0.5*(e/gg)%(e/gg) + log(gg);
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
        double w = mixWeights(i_subj, weightCol);
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
        resk += mixWeights(i_subj, weightCol) * r_ji_sq;
      }
    }
    return resk;
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
    nb_correl = as<int>(x["nb_correl"]);
    nb_fixOmega = as<int>(x["nb_fixOmega"]);
    nb_fixResid = as<int>(x["nb_fixResid"]);
    resValue = as<vec>(x["resValue"]);
    resFixed = as<uvec>(x["resFixed"]);
    resKeep = find(resFixed==0);
    niter_phi0 = as<int>(x["niter_phi0"]);
    coef_phi0 = as<double>(x["coef_phi0"]);
    nb_sa = as<int>(x["nb_sa"]);
    coef_sa = as<double>(x["coef_sa"]);
    rmcmc = as<double>(x["rmcmc"]);
    pas = as<vec>(x["pas"]);
    pash = as<vec>(x["pash"]);
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
    nb_param = nphi1 + nlambda + 1;
    nphi = nphi1+nphi0;
    Plambda.zeros(nlambda);
    ilambda1 = as<uvec>(x["ilambda1"]);
    ilambda0 = as<uvec>(x["ilambda0"]);

    DYF = zeros<mat>(mlen, nM);
    phi.set_size(N, nphi, nmc);

    //FIXME
    nendpnt=as<int>(x["nendpnt"]);
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
    lambda = as<vec>(x["lambda"]);
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
    // Pre-allocate per-chain scratch buffers for the distribution==1 hot loops
    _scratch_ft.set_size(ntotal);
    _scratch_limitT.set_size(ntotal);
    _scratch_ftT.set_size(ntotal);
    _scratch_g.set_size(ntotal);
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
    int nprint = parHistThetaKeep.n_elem + parHistOmegaKeep.n_elem + resKeep.n_elem + (nMix > 1 ? nMix - 1 : 0);
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

    distribution=as<int>(x["distribution"]);
    DEBUG=as<int>(x["DEBUG"]);
    phiMFile=as<std::vector< std::string > >(x["phiMFile"]);
    //Rcout << phiMFile[0];

  }

  void saem_fit() {
    //arma_rng::set_seed(99);
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
    for (unsigned int kiter=0; kiter<(unsigned int)(niter); kiter++) {
      gamma2_phi1=Gamma2_phi1.diag();
      IGamma2_phi1=inv_sympd(Gamma2_phi1);
      D1Gamma21=LCOV1*IGamma2_phi1;
      D2Gamma21=D1Gamma21*LCOV1.t();
      CGamma21=COV21%D2Gamma21;

      gamma2_phi0=Gamma2_phi0.diag();
      IGamma2_phi0=inv_sympd(Gamma2_phi0);
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
      vec resy(nmc);
      vec fsM;
      fsM.set_size(0);

      if (nMix > 1 && mixSampleMethod == 1) {
        // MSAEM (Lavielle & Mbogning 2014): S-step uses a single MCMC trajectory per
        // subject (phiM) targeting the marginal mixture density; no label is simulated/masked.
        vec U_y = mixObsLoss(phiM, mx);
        if (nphi1 > 0) {
          vec U_phi;
          do_mcmc_msaem(1, nu1, mx, mphi1, phiM, U_y, U_phi);
          mat dphi = phiM.cols(i1) - mphi1.mprior_phiM;
          U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi1), 1);
          do_mcmc_msaem(2, nu2, mx, mphi1, phiM, U_y, U_phi);
          do_mcmc_msaem(3, nu3, mx, mphi1, phiM, U_y, U_phi);
        }
        if (nphi0 > 0) {
          vec U_phi;
          do_mcmc_msaem(1, nu1, mx, mphi0, phiM, U_y, U_phi);
          mat dphi = phiM.cols(i0) - mphi0.mprior_phiM;
          U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi0), 1);
          do_mcmc_msaem(2, nu2, mx, mphi0, phiM, U_y, U_phi);
          do_mcmc_msaem(3, nu3, mx, mphi0, phiM, U_y, U_phi);
        }
        if (DEBUG > 0) Rcout << "mcmc successful (msaem)\n";
        phiFile << phiM;

        // E-step: posterior responsibility gamma_{i,m} (softmax, see mixWeights below) from
        // the one simulated phi, plus per-hypothesis predictions for the residual term below.
        field<vec> fsave_hyp(nMix);
        mat Ly(N, nMix, fill::zeros);
        for (int mHyp = 0; mHyp < nMix; mHyp++) {
          current_saem_state->_saemMixest = mHyp + 1;
          mat hypMat = user_fn(phiM, evt, optM);
          fsave_hyp(mHyp) = hypMat.col(0);
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
              _scratch_g = vecares + vecbres % abs(_scratch_ftT);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = indio + (arma::uword)k * stride;
              DYFhyp(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              for (int j = ntotal; j--;) {
                DYFhyp(_scratch_indio(j)) = doCensNormal1(censk[j], y[j], _scratch_limitT[j],
                                                       DYFhyp(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
              }
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
          if (sumW > 0.0) {
            mixWeights.row(i) = w_i / sumW;
          } else {
            mixWeights.row(i) = mixProb.t();
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
              resk += arResk(b, f_cur, y_cur, mHyp);
            }
            statr[b] += resk;
            resy(k) = resk;
          }

          mat dphi1k = phi1k - mprior_phi1;
          mat dphi0k = phi0k - mprior_phi0;
          vec sdg1 = sum(dphi1k % dphi1k, 0).t() / gamma2_phi1;
          mat Md1 = (IGamma2_phi1 * (dphi1k.t() * Mcovariables)).t();
          mat Md0 = (IGamma2_phi0 * (dphi0k.t() * Mcovariables)).t();
          vec d1_mu_phi1 = Md1(ind_cov1);
          vec d1_mu_phi0 = Md0(ind_cov0);
          vec d1_loggamma2_phi1 = 0.5 * sdg1 - 0.5 * N;
          vec d1_logsigma2(1);
          d1_logsigma2[0] = 0.5 * resy(k) / sigma2[0] - 0.5 * ntotal;
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
          d2logk(nb_param - 1, nb_param - 1) = -0.5 * resy(k) / sigma2[0];
          D2 = D2 + d2logk;
        }
      } else if (nMix > 1) {
        field<mat> DYF_mix(nMix);
        field<vec> U_y_mix(nMix);

        for (int jMix = 0; jMix < nMix; jMix++) {
          current_saem_state->_saemMixest = jMix + 1;

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
              _scratch_g = vecares + vecbres % abs(_scratch_ftT);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = indio + (arma::uword)k * stride;
              cur_DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              for (int j = ntotal; j--;) {
                cur_DYF(_scratch_indio(j)) = doCensNormal1(censk[j], y[j], _scratch_limitT[j],
                                                       cur_DYF(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
              }
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
            do_mcmc(1, nu1, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
            mat dphi = cur_phiM.cols(i1) - mphi1.mprior_phiM;
            U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi1), 1);
            do_mcmc(2, nu2, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
            do_mcmc(3, nu3, mx, mphi1, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
          }
          if (nphi0 > 0) {
            vec U_phi;
            do_mcmc(1, nu1, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
            mat dphi = cur_phiM.cols(i0) - mphi0.mprior_phiM;
            U_phi = 0.5 * sum(dphi % (dphi * IGamma2_phi0), 1);
            do_mcmc(2, nu2, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
            do_mcmc(3, nu3, mx, mphi0, cur_DYF, cur_phiM, U_y, U_phi, cur_fsave, cur_cens, cur_limit);
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
          current_saem_state->_saemMixest = jMix + 1;
          cur_fsave = user_fn(cur_phiM, evt, optM).col(0);
        }
        current_saem_state->_saemMixest = 0;

        // Compute posterior weights a_ji
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

        if (DEBUG > 0) Rcout << "mcmc successful (mixture)\n";

        mat phiM_weighted(N * nmc, nphi, fill::zeros);
        for (int jMix = 0; jMix < nMix; jMix++) {
          for (int row = 0; row < N * nmc; row++) {
            int i_subj = row % N;
            phiM_weighted.row(row) += mixWeights(i_subj, jMix) * phiM_mix(jMix).row(row);
          }
        }
        phiFile << phiM_weighted;
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
              resk += arResk(b, f_cur, y_cur, jMix);
            }
            statr[b] += resk;
            resy(k) = resk;
          }

          vec sdg1 = sdg1_w / gamma2_phi1;
          mat Md1 = (IGamma2_phi1 * (dphi1k_w.t() * Mcovariables)).t();
          mat Md0 = (IGamma2_phi0 * (dphi0k_w.t() * Mcovariables)).t();
          vec d1_mu_phi1 = Md1(ind_cov1);
          vec d1_mu_phi0 = Md0(ind_cov0);
          vec d1_loggamma2_phi1 = 0.5 * sdg1 - 0.5 * N;
          vec d1_logsigma2(1);
          d1_logsigma2[0] = 0.5 * resy(k) / sigma2[0] - 0.5 * ntotal;
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
          d2logk(nb_param - 1, nb_param - 1) = -0.5 * resy(k) / sigma2[0];
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
            _scratch_g = vecares + vecbres % abs(_scratch_ftT);
            _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
            _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
            _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
            _scratch_indio = indio + (arma::uword)k * stride;
            DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
            for (int j = ntotal; j--;) {
              DYF(_scratch_indio(j)) = doCensNormal1(censk[j], y[j], _scratch_limitT[j],
                                                     DYF(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
            }
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
        }
        else {
          RSprintf("unknown distribution (id=%d)\n", distribution);
          return;
        }
        //U_y is a vec of subject llik; summed over obs for each subject
        vec U_y=sum(DYF, 0).t();

        if(nphi1>0) {
          vec U_phi;
          do_mcmc(1, nu1, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit);
          mat dphi = phiM.cols(i1)-mphi1.mprior_phiM;
          U_phi    = 0.5*sum(dphi%(dphi*IGamma2_phi1),1);
          do_mcmc(2, nu2, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit);
          do_mcmc(3, nu3, mx, mphi1, DYF, phiM, U_y, U_phi, fsave, cens, limit);
        }
        if(nphi0>0) {
          vec U_phi;
          do_mcmc(1, nu1, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit);
          mat dphi = phiM.cols(i0)-mphi0.mprior_phiM;
          U_phi    = 0.5*sum(dphi%(dphi*IGamma2_phi0),1);
          do_mcmc(2, nu2, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit);
          do_mcmc(3, nu3, mx, mphi0, DYF, phiM, U_y, U_phi, fsave, cens, limit);
        }
        if (DEBUG>0) Rcout << "mcmc successful\n";
        phiFile << phiM;

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
          vec gk, y_cur, f_cur;
          double ft, fa;
          //loop thru endpoints here
          for(int b=0; b<nendpnt; ++b) {
            if (hasFixedObsTransform) {
              y_cur = ysTrans(span(y_offset(b), y_offset(b+1)-1));
            } else {
              y_cur = ys(span(y_offset(b), y_offset(b+1)-1));
            }
            f_cur = fk(span(y_offset(b), y_offset(b+1)-1));
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

            if (res_mod(b) <= rmProp) {
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
            resy(k) = resk;                                          //FIXME: resy(b,k)?
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
          vec d1_logsigma2(1);
          d1_logsigma2[0] =  0.5*resy(k)/sigma2[0]-0.5*ntotal; //FIXME: sigma2[0], sigma2[b] instead?
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
          d2logk(nb_param-1,nb_param-1)=-0.5*resy(k)/sigma2[0];      //FIXME: sigma2[0], sigma2[b] instead?
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
      mprior_phi1=COV1*MCOV1;
      mprior_phi0=COV0*MCOV0;
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

      if (kiter<=(unsigned int)(nb_sa)) {
        Gamma2_phi1=max(Gamma2_phi1*coef_sa, diagmat(G1));
      } else {
        Gamma2_phi1=G1;
      }
      Gamma2_phi1=Gamma2_phi1%covstruct1;
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
          }
          break;
        }
        sigma2[b] = sig2;                                          //CHK: sigma2[] use
        if (sigma2[b]>1.0e99) sigma2[b] = 1.0e99;
        if (std::isnan(sigma2[b])) sigma2[b] = 1.0e99;
      }
      vecares = ares(ix_endpnt);
      vecbres = bres(ix_endpnt);
      if (DEBUG>0) Rcout << "par update successful\n";

      //    Fisher information
      DDa=(D1/nmc)*(D1/nmc).t()-D11/nmc-D2/nmc;
      DDb=-D11/nmc-D2/nmc;
      L=L+pash(kiter)*(D1/nmc-L);
      Ha=Ha+pash(kiter)*(DDa- Ha);
      Hb=Hb+pash(kiter)*(DDb- Hb);
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
      if (nMix > 1) {
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
      g2 = vcsig2.elem(resKeep);
      pl = join_cols(pl, g2);
      if (nMix > 1) {
        vec mixP = mixProb.head(nMix - 1);
        pl = join_cols(pl, mixP);
      }
      par_hist.row(kiter) = pl.t();
      // saem has no per-iteration objective function (scale.showOfv=0, so `f` is ignored here);
      // scalePrintFun gates printing, checks for user interrupt, and flushes solve warnings.
      scalePrintFun(&scale, pl.memptr(), NA_REAL);
    }//kiter
    phiFile.close();
  }


private:

  user_funct user_fn;

  uvec nu;
  int niter;
  int nb_sa;
  int nb_correl;
  int nb_fixOmega;
  int nb_fixResid;
  int niter_phi0;
  double coef_phi0;
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

  mat statphi01, statphi02, statphi11, statphi12;
  // Per-component, unblended sufficient statistic (never mixed across components); used to fix
  // tcl1/tcl2-style split-ETA fixed effects via weighted regression (see omegaShareSubpop block).
  field<mat> statphi11_mix;
  double statrese[MAXENDPNT];
  double sigma2[MAXENDPNT];
  vec ares, bres, cres, lres, lambda, low, hi;
  vec vecares, vecbres, veccres, veclres;
  uvec res_mod, yj, propT, addProp;

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

  int distribution;

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

  // Per-chain scratch buffers pre-allocated in inits() to avoid repeated heap
  // allocation in the hot distribution==1 loops in saem_fit() and do_mcmc().
  vec _scratch_ft;      // transformed-f per chain (replaces ftk/fck)
  vec _scratch_limitT;  // transformed-limit per chain (replaces limitTk)
  vec _scratch_ftT;     // handleF output per chain (replaces ftTk/fcTk)
  vec _scratch_g;       // residual SD per chain (replaces gk/gck)
  uvec _scratch_indio;  // DYF row indices per chain (replaces indio_k)

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

  static inline void doCens(mat &DYF, vec &cens, vec &limit, vec &fc, vec &r, const vec &dv) {
    for (int j = (int)cens.size(); j--;) {
      DYF(j) = doCensNormal1(cens[j], dv[j], limit[j], DYF(j), fc[j], r[j], 0);
    }
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
               vec &cur_limit) {
    mat fcMat;
    vec fc, fs, Uc_y, Uc_phi, deltu;
    uvec ind;

    uvec i=mphi.i;
    double double_xmin = 1.0e-200;                               //FIXME hard-coded xmin, also in neldermean.hpp
    double xmax = 1e300;
    for (int u=0; u<nu; u++)
      for (int k1=0; k1<mphi.nphi; k1++) {
        mat phiMc=phiM;
        switch (method) {
        case 1:
          phiMc.cols(i)=randn<mat>(mx.nM,mphi.nphi)*mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          break;
        case 2:
          phiMc.cols(i)=phiM.cols(i) +
            randn<mat>(mx.nM,mphi.nphi)*mphi.Gdiag_phi % current_saem_state->_saemUE.cols(i);
          break;
        case 3:
          phiMc.col(i(k1))=phiM.col(i(k1))+
            randn<vec>(mx.nM)*mphi.Gdiag_phi(k1,k1) % current_saem_state->_saemUE.col(i(k1));
          // Rcpp::print(Rcpp::wrap(phiM.cols(i(k1))));
          break;
        }

        fcMat = user_fn(phiMc, mx.evtM, mx.optM);
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
              _scratch_g = vecares + vecbres % abs(_scratch_ftT);
              _scratch_g.elem(find(_scratch_g == 0.0)).fill(1);
              _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
              _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
              _scratch_indio = mx.indio + (arma::uword)k * stride;
              DYF(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
              for (int j = ntotal; j--;) {
                DYF(_scratch_indio(j)) = doCensNormal1(censk[j], mx.y[j], _scratch_limitT[j],
                                                       DYF(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
              }
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

        ind=find( deltu < -log(randu<vec>(mx.nM)) );
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
      mat fcMat = user_fn(phiC, mx.evtM, mx.optM);
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
              _scratch_ftT(i) = handleF(propT(cur), fsk(i), _scratch_ft(i), false, true);
            }
            _scratch_g = vecares + vecbres % abs(_scratch_ftT);
            _scratch_g.elem(find(_scratch_g == 0.0)).fill(1);
            _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
            _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
            _scratch_indio = mx.indio + (arma::uword)k * stride;
            DYFm(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
            for (int j = ntotal; j--;) {
              DYFm(_scratch_indio(j)) = doCensNormal1(censk[j], mx.y[j], _scratch_limitT[j],
                                                     DYFm(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
            }
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
            _scratch_ftT(i) = handleF(propT(cur), fk(i), _scratch_ft(i), false, true);
          }
          _scratch_g = vecares + vecbres % abs(_scratch_ftT);
          _scratch_g.elem(find(_scratch_g == 0.0)).fill(1.0);
          _scratch_g.elem(find(_scratch_g < double_xmin)).fill(double_xmin);
          _scratch_g.elem(find(_scratch_g > xmax)).fill(xmax);
          _scratch_indio = indio + (arma::uword)k * stride;
          DYFhyp(_scratch_indio) = arDYFhyp(yt, _scratch_ft, _scratch_g);
          for (int j = ntotal; j--;) {
            DYFhyp(_scratch_indio(j)) = doCensNormal1(censk[j], y[j], _scratch_limitT[j],
                                                   DYFhyp(_scratch_indio(j)), _scratch_ft[j], _scratch_g[j], 0);
          }
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
                      vec &U_phi) {
    mat phiMc;
    vec Uc_y, Uc_phi, deltu;
    uvec ind;
    uvec i = mphi.i;

    for (int u = 0; u < nu; u++)
      for (int k1 = 0; k1 < mphi.nphi; k1++) {
        phiMc = phiM;
        switch (method) {
        case 1:
          phiMc.cols(i) = randn<mat>(mx.nM, mphi.nphi) * mphi.Gamma_phi % current_saem_state->_saemUE.cols(i) +
            mphi.mprior_phiM;
          break;
        case 2:
          phiMc.cols(i) = phiM.cols(i) +
            randn<mat>(mx.nM, mphi.nphi) * mphi.Gdiag_phi % current_saem_state->_saemUE.cols(i);
          break;
        case 3:
          phiMc.col(i(k1)) = phiM.col(i(k1)) +
            randn<vec>(mx.nM) * mphi.Gdiag_phi(k1, k1) % current_saem_state->_saemUE.col(i(k1));
          break;
        }

        Uc_y = mixObsLoss(phiMc, mx);

        if (method == 1) {
          deltu = Uc_y - U_y;
        } else {
          mat dphic = phiMc.cols(i) - mphi.mprior_phiM;
          Uc_phi = 0.5 * sum(dphic % (dphic * mphi.IGamma2_phi), 1);
          deltu = Uc_y - U_y + Uc_phi - U_phi;
        }

        ind = find(deltu < -log(randu<vec>(mx.nM)));
        phiM(ind, i) = phiMc(ind, i);
        U_y(ind) = Uc_y(ind);
        if (method > 1) {
          U_phi(ind) = Uc_phi(ind);
        }
        if (method < 3) {
          break;
        }
      }
  }
};


// closing for #ifndef __SAEM_CLASS_RCPP_HPP__
#endif

using namespace Rcpp;



t_calc_lhs saem_lhs = NULL;
t_update_inis saem_inis = NULL;


rx_solve* _rx = NULL;

RObject mat2NumMat(const mat &m) {
  RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
  x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
  return x;
}

CharacterVector parNames;

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


  // Fill in subject parameter information
  for (int _i = 0; _i < _Nnlmixr2; ++_i) {
    ind = getSolvingOptionsInd(_rx, _i);
    setIndSolve(ind, -1);
    if (current_saem_state != nullptr && current_saem_state->_saemMixest != 0) {
      setIndMixest(ind, current_saem_state->_saemMixest);
    } else if (indMixest != nullptr) {
      setIndMixest(ind, indMixest[_i]);
    }
    int k=0;
    for (int _j = 0; _j < nPar; _j++){
      if (doParam[_j] == 1) {
        setIndParPtr(ind, _j, _phi(_i, k++));
      }
    }
  }
  resetRxBadSolve(_rx);
  par_solve(_rx); // Solve the complete system (possibly in parallel)
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
    par_solve(_rx);
    j++;
  }
  if (!current_saem_state->_saemIndTolRelax && j != 0) {
    // Reset all subjects' tolFactor after the non-selective retry.
    for (int _i = 0; _i < _Nnlmixr2; _i++) {
      setIndTolFactor(getSolvingOptionsInd(_rx, _i), 1.0);
    }
  }
  // indTolRelax=TRUE: stiff subjects retain their loosened tolFactor across iterations.
  mat g(getRxNsim(_rx) * getRxNobs2(_rx), 3); // nobs across all chains
  int elt=0;
  bool hasNan = false;
  for (int id = 0; id < _Nnlmixr2; ++id) {
    ind = getSolvingOptionsInd(_rx, id);
    iniSubjectE(getOpNeq(op), 1, ind, op, _rx, saem_inis);
    for (int j = 0; j < getIndNallTimes(ind); ++j) {
      setIndIdx(ind, j);
      int kk = getIndIx(ind, getIndIdx(ind));
      double curT = getTime(kk, ind);
      double *lhs = getIndLhs(ind);
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
  rxUpdateFuns(mv["trans"], &rxInner);
  parNames = mv[RxMv_params];

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

//[[Rcpp::export]]
SEXP saem_do_pred(SEXP in_phi, SEXP in_evt, SEXP in_opt) {
  List opt = List(in_opt);
  mat phi = as<mat>(in_phi);
  setupRx(opt, in_evt, 1, (int)phi.n_rows);
  saem_lhs = rxInner.calc_lhs;
  saem_inis = rxInner.update_inis;
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
  saem_lhs = rxInner.calc_lhs;
  saem_inis = rxInner.update_inis;
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
