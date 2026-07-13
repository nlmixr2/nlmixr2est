// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#define iniRxodePtrs0 _nlmixr2est_iniRxodePtrs0
#include "armahead.h"
#include "utilc.h"
#include <lbfgsb3ptr.h>
#include "censEst.h"
#include "nearPD.h"
#include "shi21.h"
#include "inner.h"
#include "imp.h"
#include "rxomp.h"
#include "solveWarnHelper.h"
#include "vaeEncoder.h"

// scale.h uses `_("...")` for translatable strings; provide the trivial
// passthrough macro before including it (matches saem.cpp's usage).
#ifndef _
#define _(String) (String)
#endif
#include "scale.h"
#include <n1qn1c.h>
#include <Rinternals.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <atomic>
#include <chrono>
#include <memory>

// Set while inside the parallel inner optimization region: R-API calls and
// running-mean accumulation are deferred to the post-parallel phase. Atomic
// to avoid a TSan-flagged data race on worker threads' reads.
static std::atomic<int> _innerParallel{0};

extern "C" {
#define iniLbfgsb3ptr _nlmixr2est_iniLbfgsb3ptr
  iniLbfgsb3
#define iniRxodePtrs _nlmixr2est_iniRxodePtrs
  iniRxode2ptr
#define iniN1qn1cPtrs _nlmixr2est_iniN1qn1cPtrs
  iniN1qn1c
}

#define _(String) (String)

#define PHI(x) 0.5*(1.0+erf((x)/M_SQRT2))

using focei_wall_clock = std::chrono::steady_clock;

static inline double foceiElapsedSeconds(const focei_wall_clock::time_point& start) {
  return std::chrono::duration<double>(focei_wall_clock::now() - start).count();
}

void saveIntoEnvrionment(Environment e);
void restoreFromEnvrionment(Environment e);

#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define innerOde(id) ind_solve(rx, getRxId(id), rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)
#define predOde(id) ind_solve(rx, getRxId(id), rxPred.dydt_liblsoda, rxPred.dydt_lsoda_dum, rxPred.jdum_lsoda, rxPred.dydt, rxPred.update_inis, rxPred.global_jt)
#define thetaSensOde(id) ind_solve(rx, getRxId(id), rxThetaSens.dydt_liblsoda, rxThetaSens.dydt_lsoda_dum, rxThetaSens.jdum_lsoda, rxThetaSens.dydt, rxThetaSens.update_inis, rxThetaSens.global_jt)
#define getCholOmegaInv() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "chol.omegaInv", R_NilValue)))
#define getOmega() (as<NumericMatrix>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaMat() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omega", R_NilValue)))
#define getOmegaInv() (as<arma::mat>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "omegaInv", R_NilValue)))
#define getOmegaDet() (as<double>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "log.det.OMGAinv.5", R_NilValue)))
#define getOmegaN() as<int>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "ntheta", R_NilValue))
#define getOmegaTheta() as<NumericVector>(rxode2::rxSymInvCholEnvCalculate(_rxInv, "theta", R_NilValue));
#define setOmegaTheta(x) rxode2::rxSymInvCholEnvCalculate(_rxInv, "theta", x)
#define tbs(x) _powerD(x,    getIndLambda(ind), getIndLambdaYj(ind), getIndLogitLow(ind), getIndLogitHi(ind))
#define tbsL(x) _powerL(x,   getIndLambda(ind), getIndLambdaYj(ind), getIndLogitLow(ind), getIndLogitHi(ind))
#define tbsDL(x) _powerDL(x, getIndLambda(ind), getIndLambdaYj(ind), getIndLogitLow(ind), getIndLogitHi(ind))
#define tbsD(x) _powerDD(x,  getIndLambda(ind), getIndLambdaYj(ind), getIndLogitLow(ind), getIndLogitHi(ind))
#define _safe_log(a) (((a) <= 0.0) ? log(DBL_EPSILON) : log(a))
// #define _safe_log(a) log(a)
#define _safe_zero(a) ((a) == 0 ? DBL_EPSILON : (a))
//#define _safe_zero(a) (a)
#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))
//#define _safe_sqrt(a) sqrt(a)
#define _as_dbleps(a) (fabs(a) < sqrt(DBL_EPSILON) ? ((a) < 0 ? -sqrt(DBL_EPSILON)  : sqrt(DBL_EPSILON)) : a)

#define expit(alpha, low, high) _powerDi(alpha, 1.0, 4, low, high)
#define probitInv(alpha, low, high) _powerDi(alpha, 1.0, 6, low, high)

extern "C" {
  typedef double optimfn(int n, double *par, void *ex);

  typedef void optimgr(int n, double *par, double *gr, void *ex);

  void lbfgsbRX(int n, int lmm, double *x, double *lower,
                double *upper, int *nbd, double *Fmin, optimfn fn,
                optimgr gr, int *fail, void *ex, double factr,
                double pgtol, int *fncount, int *grcount,
                int maxit, char *msg, int trace, int nREPORT);

}

bool assignFn_ = false;

List _rxInv;

// These are focei inner options
struct focei_options {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  List mvi;
  double *etaUpper = NULL;
  double *etaLower = NULL;
  int *nbdInner = NULL;
  double *geta = NULL;
  double *getahf =NULL;
  double *getahr = NULL;
  double *getahh = NULL;
  double *goldEta = NULL;
  double *gtryEta = NULL;
  double *gsaveEta = NULL;
  double *gthetaGrad = NULL;
  bool mGthetaGrad = false;
  // n1qn1 specific vectors
  double *gZm = NULL;
  double *gZmH = NULL;
  double *gZmEta = NULL;
  double *gG = NULL;
  double *gVar = NULL;
  double *gX = NULL;

  double *glp = NULL;
  double *ga = NULL;
  double *gB = NULL;
  double *gc = NULL;
  // per-observation censored inner-Hessian coefficients (rho_ff/rho_fR/rho_RR): the
  // EXACT censored 2nd derivative for the Laplace determinant; normal obs carry the
  // Gauss-Newton 1/r, 0, 0.5/r^2 (so calcEtaHessian reduces bit-identically there).
  double *gcHff = NULL;
  double *gcHfr = NULL;
  double *gcHrr = NULL;
  double *gH = NULL;
  double *gVid = NULL;

  double *likSav = NULL;
  double *llikObsFull = NULL;

  // Integer of ETAs
  unsigned int gEtaGTransN;
  // Where likelihood is saved.

  int *etaTrans = NULL;
  int *etaFD = NULL;
  int *mixTrans = NULL;
  int predNeq;
  int eventType;

  // Index of rx_pred_ in the inner model's lhs.  Normally 0, but an AR(1)
  // endpoint emits lag()-referenced defs (the residual/time the lag needs as
  // real lhs) ahead of rx_pred_, so it shifts.  The d(f)/d(eta), rx_r_ and
  // d(r)/d(eta) columns follow rx_pred_ contiguously; located by name at setup.
  int predOffset;
  int predNoLhsOffset; // same, for the predNoLhs model used in the FD fallback
  int thetaSensOffset = -1;   // lhs offset of the first d(f)/d(theta) output (impmap)
  int thetaSensDvOffset = -1; // lhs offset of the first d(V)/d(theta) output (impmap)
  int thetaSensPredOffset = -1; // lhs offset of rx_pred_ (f) in the sensitivity model
  int thetaSensROffset = -1;    // lhs offset of rx_r_ (V) in the sensitivity model
  int thetaSensNeq = 0;       // ODE state count of the sensitivity model (impmap)
  int innerNeq = 0;           // inner model state count when the pool is sized larger (impmap)

  unsigned int neta;
  unsigned int ntheta;
  unsigned int npars;
  unsigned int thetan;
  unsigned int omegan;

  int calcGrad;
  int nF;
  int nF2;
  int nG;
  int derivMethod;
  int covDerivMethod;
  int covMethod;
  int derivMethodSwitch;
  double derivSwitchTol;
  double lastOfv;

  double *fullTheta = NULL;
  double *theta = NULL;
  double *thetaGrad = NULL;
  double *initPar = NULL;
  double *scaleC = NULL;
  // Adaptive Gaussian Quadrature
  //
  double *aqx;
  double *aqw;
  double aqLow;
  double aqHi;
  bool aqfirst = false;
  double scaleC0;
  int *xPar = NULL;
  // probitIdxArr (npars-sized): 0 = no probit transform, k>=1 = k-th probit
  // transform (bounds at probitThetaLow/Hi[k-1]). Kept separate from xPar's
  // logit encoding so it doesn't collide with xPar's omega scaling codes 2-5.
  int *probitIdxArr = NULL;
  NumericVector lowerIn;
  double *lower = NULL;
  NumericVector upperIn;
  double *upper = NULL;
  int *nbd = NULL;

  int *fixedTrans = NULL;
  int *thetaTrans = NULL;

  int scaleType;
  int normType;
  double scaleCmin;
  double scaleCmax;
  double c1;
  double c2;
  double scaleTo;
  double epsilon;

  int maxOuterIterations;
  int maxInnerIterations;

  double odeRecalcFactor;
  int maxOdeRecalc;
  int objfRecalN;
  int stickyRecalcN1;
  std::atomic<int> stickyRecalcN2{0};
  int stickyRecalcN;
  std::atomic<int> stickyTol{0};
  bool indTolRelax;

  int nsim;
  unsigned int nzm;
  int warm; // 1 = seed zm from calculated eta Hessian, 0 = classic behavior

  int imp;
  // int printInner;


  mat omega;
  mat omegaInv;
  mat cholOmegaInv;
  mat etaM;
  mat etaS;
  mat eta1SD;
  double n;
  double logDetOmegaInv5;

  int    *gillRetC = NULL;
  int    *gillRet = NULL;
  double *gillDf = NULL;
  double *gillDf2 = NULL;
  double *gillErr = NULL;
  double *rEps = NULL;
  double *aEps = NULL;
  double *rEpsC = NULL;
  double *aEpsC = NULL;

  //
  double factr;
  double pgtol;
  double abstol;
  double reltol;
  int lmm;
  int *skipCov = NULL;
  unsigned int skipCovN;

  int outerOpt;
  int eigen;
  int scaleObjective;
  double scaleObjectiveTo;
  int initObj;
  double initObjective;
  // Confidence Interval
  double ci;
  double sigdig;
  //
  clock_t t0 = clock();
  int cur = 0;
  int curTick = 0;
  int totTick = 100;
  double boundTol;
  int noabort;
  int interaction;
  int foceType; // FOCE residual-variance R: 0 = "nonmem" (eta=0 frozen R), 1 = "foce+" (live conditional R)
  int fast;     // analytic ("fast") outer gradient + Eq-48 eta extrapolation
  int curAnalytic = 0;      // this gradient came from the analytic ("fast") path
  int nAnalyticGrad = 0;    // # outer gradients from the analytic path
  int nFDGradFast = 0;      // # FD fallbacks while fast was requested
  int warnedAnalyticFallback = 0; // one-time FD-fallback warning latch
  double cholSEtol;
  double hessEps;
  double hessEpsLlik;
  double hessEpsInner;
  int shi21maxOuter;
  int shi21maxInner;
  int shi21maxInnerCov;
  int shi21maxFD;
  double cholAccept;
  double resetEtaSize;
  std::atomic<int> didEtaReset{0};
  double resetThetaSize = std::numeric_limits<double>::infinity();
  double resetThetaFinalSize = std::numeric_limits<double>::infinity();
  int checkTheta;
  int *muRef = NULL;
  unsigned int muRefN;
  // Per-eta opt-out (indexed like muRef) from eta-drift zero-reset and mu-ref
  // theta soft-shift; NULL (default) protects no etas.
  int *muRefEtaCovSkipReset = NULL;
  // mufocei/irlsfocei group structure: one group per mu-ref-covariate
  // population theta with an eta. Arrays NULL/0 when muModel="none".
  // muModel: 0=none, 1=lin (OLS), 2=irls (reweighted).
  int muModel = 0;
  unsigned int muGroupN = 0;
  int *muGroupTheta = NULL;      // [muGroupN] population theta index (0-based, ntheta space)
  int *muGroupEta = NULL;        // [muGroupN] eta index (0-based, neta space)
  int *muGroupCovStart = NULL;   // [muGroupN] start offset into the flattened covariate arrays
  int *muGroupCovCount = NULL;   // [muGroupN] number of covariates in this group
  unsigned int muGroupCovN = 0;  // total flattened covariate count across all groups
  int *muGroupCovTheta = NULL;   // [muGroupCovN] covariate-coefficient theta index (0-based, ntheta space)
  arma::mat muGroupCovData;      // nsub x muGroupCovN, baseline covariate values per subject
  // Per-theta skip flag (ntheta space): excludes group thetas from the outer
  // optimizer's free-parameter set. Built once at setup.
  int *muGroupThetaSkip = NULL;  // [ntheta]
  // User-fixed (ini(...~fix)) flattened covariate-coefficient thetas. Needed
  // because once a theta is in muGroupThetaSkip it's no longer in fixedTrans,
  // so isFixedTheta() can't tell "user fixed" from "mu-group excluded";
  // updateMuGroups() uses this to hold user-fixed coefficients out of the
  // regression design matrix while merely-excluded ones are estimated normally.
  int *muGroupCovUserFixed = NULL; // [muGroupCovN]
  // Flattened covariate-coefficient thetas with a finite bound
  // (ini(...~c(lower,est,upper))). Handled like user-fixed for the regression
  // (excluded from the design matrix, live gradient-updated offset subtracted)
  // but NOT added to muGroupThetaSkip, since it must stay in fixedTrans/npars
  // so the outer optimizer keeps optimizing it subject to its bound.
  int *muGroupCovBounded = NULL; // [muGroupCovN]
  // Display names for mu-group thetas, used only by printMuGroupThetaRow() to
  // label the extra print row; excluded from op_focei.scale's own column names.
  CharacterVector muGroupThetaNames;    // [muGroupN]
  CharacterVector muGroupCovThetaNames; // [muGroupCovN]
  // A single updateMuGroups() step per outer iteration isn't always enough:
  // the outer optimizer's convergence check ignores mu-group parameters, so
  // it can converge while the mu-group regress/re-optimize sub-problem is
  // still moving. innerOpt() instead cycles {re-optimize etas, updateMuGroups()}
  // up to muGroupMaxCycles times (or until movement < muGroupTol).
  double muGroupTol = 1e-3;
  int muGroupMaxCycles = 10;
  int resetHessianAndEta;
  std::atomic<int> didHessianReset{0};
  int cholSEOpt=0;
  int cholSECov;
  int fo;
  int covTryHarder;
  // Gill options
  int gillK;
  double gillStep;
  double gillFtol;
  double gillRtol;
  int gillKcov;
  int gillKcovLlik;
  double gillStepCov;
  double gillStepCovLlik;
  double gillFtolCov;
  double gillFtolCovLlik;
  double covSmall;
  int didGill;
  int smatNorm;
  int smatNormLlik;
  int rmatNorm;
  int rmatNormLlik;
  int covGillF;
  int optGillF;
  int mixDeriv;
  int covFull;       // covFull=TRUE: report the full theta+sigma+Omega covariance
  int covFdDirect;   // 1 while the FD-full cov perturbs fullTheta/Omega directly:
                     // updateTheta then skips the unscale + Omega rebuild (below)
  double gradTrim;
  double gradCalcCentralSmall;
  double gradCalcCentralLarge;
  double etaNudge;
  double etaNudge2;
  std::atomic<int> didEtaNudge{0};
  std::atomic<int> reducedTol{0};
  std::atomic<int> reducedTol2{0};
  int repeatGill;
  int repeatGillN;
  int repeatGillMax;
  int curGill;
  int printTop;
  double resetThetaCheckPer;
  int slow = 0;
  double gradProgressOfvTime;
  bool alloc=false;
  bool zeroGrad = false;
  int nfixed=0;
  NumericVector logitThetaHi;
  NumericVector logitThetaLow;
  NumericVector probitThetaHi;
  NumericVector probitThetaLow;
  double badSolveObjfAdj = 100;
  std::atomic<bool> didPredSolve{false};
  bool canDoFD  = false;
  bool adjLik = false;
  bool fallbackFD = false;
  bool needOptimHess = false;
  int optimHessType = 1;
  int optimHessCovType = 1;
  int censOption = 0;   // 0 gauss (historic Gauss-Newton, default) / 1 laplace (exact censored 2nd deriv)
  double smatPer;
  std::atomic<bool> didLikCalc{false};
  bool zeroGradFirstReset= false;
  bool zeroGradRunReset=false;
  bool zeroGradBobyqa=false;
  bool zeroGradBobyqaRun=false;
  int nEstOmega=0;
  int mceta= -1; // number of mc samples of ETA

  // Almquist Eq-48 warm-start extrapolation (fast=TRUE): per-subject d eta*/d(theta)
  // in the SCALED optimizer parameterization (etaP columns pre-multiplied by
  // dUnscaleParDx), the scaled theta snapshot at gradient time, and a validity flag.
  double *getaP = NULL;      // [neta * npars * nsub], indexed by id
  double *etaPTheta = NULL;  // [npars] scaled theta at the last analytic-gradient call
  int etaPValid = 0;
  size_t getaPn = 0;         // allocated element count (realloc guard)

  // Pre-drawn ETA samples for mceta >= 1 (neta x (mceta-1) x nsub); filled
  // serially before the parallel for-loop so workers avoid R API calls.
  arma::cube mcetaSamples;

  unsigned int mixIdxN = 0;
  int *mixIdx = NULL;
  double *mixProb = NULL;
  double *mixProbGrad = NULL;

  // Iteration-print formatting shared with saem/nlm via src/scale.h; populated
  // by foceiSetupScale(). scale.save=0 since focei records history separately.
  scaling scale;
  bool isSaem = false;
  bool isNlm = false;   // nlm-family outer optimizer (censOption is inert: FD outer Hessian)
  bool isImpmap = false; // importance-sampling EM (est="impmap"/"imp"/"qrpem"); outer runs impOuter
  bool isImp = false;    // est="imp": no per-iteration MAP search (proposal at the conditional mean)
  bool isAdvi = false;   // est="advi": reuses the theta-sensitivity model for the outer ADVI gradient
  bool isQrpem = false;  // est="qrpem": the impmap kernel labeled as QRPEM (qr+sir sugar)
  int impIsample = 300;  // importance samples drawn per subject per iteration
  double impGamma = 1.0; // proposal-variance inflation factor: cov = gamma * H^-1
  int impNiter = 100;    // maximum EM iterations
  double impIaccept = 0.4;   // target importance-sampling effective-sample fraction (adapts gamma)
  double impIscaleMin = 0.1; // lower bound for adapted gamma
  double impIscaleMax = 10.0;// upper bound for adapted gamma
  double impCtol = -1.0;     // windowed-convergence tolerance on the objective (<0: derive from sigdig)
  int impNconvWindow = 10;   // trailing-iteration window for the convergence check
  bool impCov = false;       // experimental: compute the MC observed-information theta covariance
  bool impQr = false;        // quasi-random (Sobol) importance samples (QRPEM)
  bool impQrShift = true;    // Cranley-Patterson random shift of the Sobol points
  bool impQrRefresh = true;  // redraw the shift each iteration (false: one shift/subject)
  bool impSir = false;       // SIR-accelerated non-mu/sigma M-step
  int impSirSample = 30;     // SIR resampled points per subject
  int impSeed = 42;          // base seed for the per-(iter,subject) draw streams
  std::string impDiagXform = "sqrt"; // Omega diagonal parameterization for the EM Omega update
  IntegerVector impMuThetaIdx; // 0-based theta indices of simple mu intercepts (no covariates)
  IntegerVector impMuEtaIdx;   // corresponding 0-based eta indices
  IntegerVector impThetaSensIdx; // 0-based theta indices with a d(f)/d(theta) sensitivity output
  IntegerVector impOmegaFixedEta; // 0-based eta indices whose Omega diagonal is fixed
};

focei_options op_focei;

static inline size_t getRxNsubAndMix(rx_solve* rx) {
  return (size_t)getRxNsub(rx) * (op_focei.mixIdxN + 1);
}

static inline size_t getRxNallAndMix(rx_solve* rx) {
  return (size_t)getRxNall(rx) * (op_focei.mixIdxN + 1);
}

static inline int getRxId(int id) {
  return id % getRxNsub(rx);
}

// Is eta j excluded from the eta-drift zero-reset / mu-ref theta soft-shift
// because it's driven by the mufocei/irlsfocei restart-loop's linear-model
// step instead? False when NULL/unset (muModel="none").
static inline bool isMuRefCovProtected(unsigned int j) {
  return op_focei.muRefEtaCovSkipReset != NULL && j < op_focei.muRefN &&
    op_focei.muRefEtaCovSkipReset[j] != 0;
}

static inline int getRxMixFromId(int id) {
  return std::floor(id / getRxNsub(rx)) + 1;
}

unsigned int _aqn = 0;
int _nagq = 0;

struct focei_ind {
  int nInnerF;
  int nInnerG;
  double lik[3]; // lik[0] = liklihood; For central difference: lik[1] = lower lik[2] = upper
  double *eta; // Eta includes the ID number for the patient
  double *etahf;
  double *etahr;
  double *etahh;
  //
  double *thetaGrad; // Theta gradient; Calculated on the individual level for S matrix calculation
  double thVal[2]; // thVal[0] = lower; thVal[2] = upper
  //
  // F and varaibility
  unsigned int nobs;
  unsigned int setup;

  double *saveEta; // Saved when lik[0] is saved.
  double *oldEta;
  double *tryEta;

  // Likilihood gradient
  double llik;
  double *a;
  double *B;
  double *c;
  double *cHff;   // per-obs censored inner-Hessian coeffs (rho_ff/rho_fR/rho_RR)
  double *cHfr;
  double *cHrr;
  double *lp;// = mat(neta,1);

  double *g;
  double *Vid;

  double *llikObs;

  double tbsLik;

  int mode; // 1 = dont use zm, 2 = use zm.
  double *zm;
  // warm="calc": packed lower-triangle eta Hessian + the eta it was computed at
  double *zmH;
  double *zmEta;
  int zmValid;
  double *var;
  double *x;
  unsigned int uzm;
  int doChol=1;
  int doFD=0;
  int doEtaNudge;
  int badSolve=0;
  // Per-subject tolerance-loosening retry count in likInner0's inner-retry loop;
  // replaces the old shared atomic, which raced across threads. Capped at
  // op_focei.stickyRecalcN; past that, tolFactor stays loose and stickyTol latches.
  int stickyRecalcN2=0;
  int parWarnBadHess=0;  // deferred "non-positive definite Hessian" warning
  int parErrorNoEta=0;   // deferred "could not find best eta" error
  double curF;
  double curT;
  double *curS;
  int nNonNormal = 0;
  int nObs=0;
  // Mixture options
  int    *mixest;
  double *mixProb;
  double *mixProbGrad;
};

// Zero every eta except protected (mu-ref-covariate) ones (isMuRefCovProtected).
// With no protected etas (default) this is a plain std::fill to 0.
static inline void resetEtaSelective(focei_ind *fInd, int neta) {
  if (op_focei.muRefEtaCovSkipReset == NULL) {
    std::fill(&fInd->eta[0], &fInd->eta[0] + neta, 0.0);
    return;
  }
  for (int j = 0; j < neta; ++j) {
    if (!isMuRefCovProtected((unsigned int)j)) {
      fInd->eta[j] = 0.0;
    }
  }
}

// Standardized-eta "p-value" bound test (the same criterion that drives the inner
// eta reset): eta is in bound when |chol(Omega^-1) eta|_j < resetEtaSize AND
// |eta/SD_j| < resetEtaSize for every non-mu-ref-protected j.  Used by the Eq-48
// warm-start extrapolation to decide whether to accept the extrapolated eta.
static inline bool etaInBound(double *eta) {
  if (!R_FINITE(op_focei.resetEtaSize) || op_focei.resetEtaSize <= 0) return true;
  arma::mat em(op_focei.neta, 1);
  for (int j = 0; j < op_focei.neta; j++) em(j, 0) = eta[j];
  arma::mat r1 = op_focei.cholOmegaInv * em;
  for (unsigned int j = 0; j < r1.n_rows; j++) {
    if (isMuRefCovProtected(j)) continue;
    if (std::fabs(r1(j, 0)) >= op_focei.resetEtaSize) return false;
  }
  arma::mat r2 = op_focei.eta1SD % em;
  for (unsigned int j = 0; j < r2.n_rows; j++) {
    if (isMuRefCovProtected(j)) continue;
    if (std::fabs(r2(j, 0)) >= op_focei.resetEtaSize) return false;
  }
  return true;
}

focei_ind *inds_focei = NULL;

// FOCE eta=0 population-R cache.  rPop depends only on theta, so it is constant
// across the whole inner optimization and only needs recomputing when theta
// changes.  Keyed by subject id (getRxNsubAndMix).  The generation counter is
// bumped in updateTheta() (single-threaded, before the parallel inner region)
// and only read inside that region, so no locking is needed.
static std::vector<arma::vec> _foceRPopCache;
static std::vector<long> _foceRPopGen;
static long _foceRPopCurGen = 0;

// Parameter table
std::vector<int> niter;
std::vector<int> iterType;
std::vector<double> vPar;
std::vector<double> vGrad;
std::vector<int> niterGrad;
std::vector<int> gradType;

static void releaseCovSolveArgs_(); // defined with covSolveArgs_ below; teardown backstop

extern "C" void rxOptionsFreeFocei() {
  releaseCovSolveArgs_(); // release any preserved covType="analytic" solve args
  if (op_focei.etaTrans != NULL) R_Free(op_focei.etaTrans);
  op_focei.etaTrans=NULL;

  if (op_focei.fullTheta != NULL) R_Free(op_focei.fullTheta);
  op_focei.fullTheta = NULL;

  if (op_focei.etaUpper != NULL) R_Free(op_focei.etaUpper);
  op_focei.etaUpper = NULL;

  if (op_focei.gillRet != NULL) R_Free(op_focei.gillRet);
  op_focei.gillRet = NULL;

  if (op_focei.gillDf != NULL) R_Free(op_focei.gillDf);
  op_focei.gillDf = NULL;

  if (op_focei.gthetaGrad != NULL && op_focei.mGthetaGrad) R_Free(op_focei.gthetaGrad);
  op_focei.gthetaGrad = NULL;
  op_focei.mGthetaGrad = false;

  // Eq-48 warm-start buffers (standalone, not part of the contiguous per-subject block)
  if (op_focei.getaP != NULL) R_Free(op_focei.getaP);
  op_focei.getaP = NULL; op_focei.getaPn = 0;
  if (op_focei.etaPTheta != NULL) R_Free(op_focei.etaPTheta);
  op_focei.etaPTheta = NULL;
  op_focei.etaPValid = 0;

  if (op_focei.muGroupTheta != NULL) R_Free(op_focei.muGroupTheta);
  op_focei.muGroupTheta = NULL; op_focei.muGroupEta = NULL;
  op_focei.muGroupCovStart = NULL; op_focei.muGroupCovCount = NULL;
  op_focei.muGroupCovTheta = NULL; op_focei.muGroupCovUserFixed = NULL;
  op_focei.muGroupCovBounded = NULL;
  op_focei.muGroupThetaSkip = NULL;

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei=NULL;

  op_focei.alloc = false;
  op_focei.didPredSolve = false;

  // Placement-new reset (copy-assignment is deleted by std::atomic members).
  op_focei.~focei_options();
  new (&op_focei) focei_options();

  vGrad.clear();
  vPar.clear();
  iterType.clear();
  gradType.clear();
  niter.clear();
  niterGrad.clear();
}

//[[Rcpp::export]]
void freeFocei(){
  rxOptionsFreeFocei();
}


rxSolveF rxInner;
rxSolveF rxPred;
rxSolveF rxThetaSens; // est="impmap": d(f)/d(theta) model (peer of rxInner/rxPred)

void rxUpdateFuns(SEXP trans, rxSolveF *inner){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda,
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda;
  lib = CHAR(STRING_ELT(trans, 0));
  s_dydt = CHAR(STRING_ELT(trans, 3));
  s_calc_jac = CHAR(STRING_ELT(trans, 4));
  s_calc_lhs = CHAR(STRING_ELT(trans, 5));
  s_inis = CHAR(STRING_ELT(trans, 8));
  s_dydt_lsoda_dum = CHAR(STRING_ELT(trans, 9));
  s_dydt_jdum_lsoda = CHAR(STRING_ELT(trans, 10));
  s_ode_solver_solvedata = CHAR(STRING_ELT(trans, 11));
  s_ode_solver_get_solvedata = CHAR(STRING_ELT(trans, 12));
  s_dydt_liblsoda = CHAR(STRING_ELT(trans, 13));
  inner->global_jt = 2;
  inner->global_mf = 22;
  inner->global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    inner->global_jt = 1;
    inner->global_mf = 21;
  } else {
    inner->global_jt = 2;
    inner->global_mf = 22;
  }
  inner->calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  inner->dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  inner->calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  inner->update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  inner->dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  inner->jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  inner->set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  inner->get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  inner->dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
}

void rxClearFuns(rxSolveF *inner){
  inner->calc_lhs              = NULL;
  inner->dydt                  = NULL;
  inner->calc_jac              = NULL;
  inner->update_inis           = NULL;
  inner->dydt_lsoda_dum        = NULL;
  inner->jdum_lsoda            = NULL;
  inner->set_solve             = NULL;
  inner->get_solve             = NULL;
  inner->dydt_liblsoda         = NULL;
}

rx_solve* rx;

////////////////////////////////////////////////////////////////////////////////
// n1qn1 functions
uvec lowerTri(mat H, bool diag = false){
  unsigned int d = H.n_rows;
  mat o(d, d, fill::ones);
  if (!diag){
    return find(trimatl(o,-1));
  } else {
    return find(trimatl(o));
  }
}

void updateZm(focei_ind *indF){
  std::fill(&indF->zm[0], &indF->zm[0]+op_focei.nzm,0.0);
  if (!indF->uzm){
    // Udate the curvature to Hessian to restart n1qn1
    int n = op_focei.neta;
    mat L = eye(n, n);
    mat D = mat(n, n, fill::zeros);
    mat H = mat(n, n);
    unsigned int l_n = n * (n + 1)/2;
    vec zmV(l_n);
    std::copy(&indF->zm[0], &indF->zm[0]+l_n, zmV.begin());
    H.elem(lowerTri(H, true)) = zmV;
    if (n == 1) H = D;
    else{
      L.elem(lowerTri(H,false)) = H.elem(lowerTri(H,0));
      D.diag() = H.diag();
      H = L*D*L.t();
    }
    // Hessian -> c.hess
    vec hessV = H.elem(lowerTri(H, true));
    std::copy(hessV.begin(),hessV.end(),&indF->zm[0]);
    indF->uzm = 1;
    indF->mode=2;
  }
}

static inline double getScaleC(int i){
  if (ISNA(op_focei.scaleC[i])) {
    switch (op_focei.xPar[i]){
    case 1: // log
      op_focei.scaleC[i]=1.0;
      break;
    case 2: // diag^2
      op_focei.scaleC[i]=1.0/fabs(op_focei.initPar[i]);
      break;
    case 3: // exp(diag)
      op_focei.scaleC[i] = 1.0/2.0;
      break;
    case 4: // Identity diagonal chol(Omega ^-1)
    case 5: // off diagonal chol(Omega^-1)
      op_focei.scaleC[i] = 1.0/(2.0*fabs(op_focei.initPar[i]));
      break;
    default:
      op_focei.scaleC[i]= 1.0/(fabs(op_focei.initPar[i]));
      break;
    }
  }
  return min2(max2(op_focei.scaleC[i], op_focei.scaleCmin),op_focei.scaleCmax);
}

////////////////////////////////////////////////////////////////////////////////
// Likelihood for inner functions
static inline double unscalePar(double *x, int i) {
  double scaleTo = op_focei.scaleTo, C=getScaleC(i);
  switch(op_focei.scaleType){
  case 1: // normalized
    return x[i]*op_focei.c2+op_focei.c1;
    break;
  case 2: // log vs linear scales and/or ranges
    if (op_focei.normType <= 5){
      scaleTo = (op_focei.initPar[i]-op_focei.c1)/op_focei.c2;
    } else if (scaleTo == 0){
      scaleTo=op_focei.initPar[i];
    }
    return (x[i]-scaleTo)*C + op_focei.initPar[i];
    break;
  case 3: // simple multiplicative scaling
    if (op_focei.scaleTo != 0){
      return x[i]*op_focei.initPar[i]/scaleTo;
    } else {
      return x[i];
    }
    break;
  case 4: // log non-log multiplicative scaling
    if (op_focei.scaleTo > 0){
      switch (op_focei.xPar[i]){
      case 1:
        return (x[i]-scaleTo) + op_focei.initPar[i];
      default:
        return x[i]*op_focei.initPar[i]/scaleTo;
      }
    } else {
      return x[i];
    }
  default:
    if (op_focei.scaleTo > 0){
      return (x[i]-scaleTo)*1 + op_focei.initPar[i];
    } else {
      return x[i];
    }
  }
  return 0;
}

static inline double scalePar(double *x, int i){
  double scaleTo = op_focei.scaleTo, C=getScaleC(i);
  switch(op_focei.scaleType){
  case 1:
    return (x[i]-op_focei.c1)/op_focei.c2;
  case 2:
    if (op_focei.normType <= 5){
      scaleTo = (op_focei.initPar[i]-op_focei.c1)/op_focei.c2;
    } else if (scaleTo == 0){
      scaleTo=op_focei.initPar[i];
    }
    return (x[i]-op_focei.initPar[i])/C + scaleTo;
    break;
  case 3: // simple multiplicative scaling
    if (op_focei.scaleTo > 0){
      return x[i]/op_focei.initPar[i]*op_focei.scaleTo;
    } else {
      return x[i];
    }
    break;
  case 4: // log non-log multiplicative scaling
    if (op_focei.scaleTo > 0){
      switch (op_focei.xPar[i]){
      case 1:
        return (x[i]-op_focei.initPar[i]) + op_focei.scaleTo;
      default:
        return x[i]/op_focei.initPar[i]*op_focei.scaleTo;
      }
    } else {
      return x[i];
    }
  default:
    if (op_focei.scaleTo > 0){
      return (x[i]-op_focei.initPar[i]) + op_focei.scaleTo;
    } else {
      return x[i];
    }
  }
  return 0;
}


void updateTheta(double *theta){
  // Theta is the acutal theta
  unsigned int j, k;
  // covFdDirect: the FD-full covariance has already written the natural-scale
  // op_focei.fullTheta (and op_focei.omegaInv/cholOmegaInv/logDetOmegaInv5 for the
  // Omega block) itself, so skip the unscale (and the Omega rebuild below) and just
  // propagate fullTheta into the per-subject solve.
  if (!op_focei.covFdDirect) {
    for (k = op_focei.npars; k--;){
      j=op_focei.fixedTrans[k];
      op_focei.fullTheta[j] = unscalePar(theta, k);
    }
  }
  // Update theta parameters in each individual
  rx = getRxSolve_();
  // Update theta parameters
  for (int id = getRxNsub(rx); id--;){
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    for (j = op_focei.ntheta; j--;){
      setIndParPtr(ind, op_focei.thetaTrans[j], op_focei.fullTheta[j]);
    }
  }
  // Update the mixture probabilities
  if (op_focei.mixIdxN != 0) {
    NumericVector curTheta(op_focei.ntheta);
    std::copy(&op_focei.fullTheta[0],
              &op_focei.fullTheta[0] + op_focei.ntheta,
              curTheta.begin());
    IntegerVector mixIdx(op_focei.mixIdxN);
    std::copy(&op_focei.mixIdx[0],
              &op_focei.mixIdx[0] + op_focei.mixIdxN,
              mixIdx.begin());
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    Environment nlmixr2 = loadNamespace("nlmixr2est");
    Function f = as<Function>(nlmixr2[".getMixFromLog"]);
    // Get the mix probabilities
    NumericVector mixProbs = f(curTheta, mixIdx);
    std::copy(mixProbs.begin(), mixProbs.end(), &op_focei.mixProb[0]);
    f = as<Function>(nlmixr2[".getMixJacFromLog"]);
    NumericVector mixJac = f(curTheta, mixIdx);
    std::copy(mixJac.begin(), mixJac.end(), &op_focei.mixProbGrad[0]);
  }
  // Update setOmegaTheta
  if (op_focei.neta > 0 && !op_focei.covFdDirect) {
    NumericVector omegaTheta(op_focei.omegan);
    std::copy(&op_focei.fullTheta[0] + op_focei.ntheta,
              &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
              omegaTheta.begin());
    setOmegaTheta(omegaTheta);
    if (op_focei.fo == 1){
      op_focei.omega = getOmegaMat();
    } else {
      op_focei.omegaInv = getOmegaInv();
      op_focei.cholOmegaInv = getCholOmegaInv();
      op_focei.logDetOmegaInv5 = getOmegaDet();
    }
  }
  //Now Setup Last theta
  if (!op_focei.calcGrad){
    // op_focei.estStr=sc + un + ex;
    std::copy(&theta[0], &theta[0] + op_focei.npars, &op_focei.theta[0]);
  }
  // Theta moved -> invalidate the FOCE eta=0 population-R cache (rPop is a
  // function of theta only).  Bumped even for FOCEI (cache simply unused there).
  _foceRPopCurGen++;
}

arma::mat cholSE__(arma::mat A, double tol);

typedef void (*gill83fn_type)(double *fp, double *theta, int id, int foceiGill);

void gill83fnF(double *fp, double *theta, int, int foceiGill);
int gill83(double *hf, double *hphif, double *df, double *df2, double *ef,
           double *theta, int cpar, double epsR, int K, double gillStep,
           double fTol, int cid, gill83fn_type gill83fn, int foceiGill, double gillF);

gill83fn_type gill83fnG = &gill83fnF;


void updateEta(double *eta, int cid) {
  rx_solving_options_ind *ind =  getSolvingOptionsInd(rx, getRxId(cid));
  for (int i = op_focei.neta; i--;) {
    setIndParPtr(ind, op_focei.etaTrans[i], eta[i]);
  }
}

// RAII guard: snapshots a per-individual ETA vector, restores it on
// destruction unless disarm() is called. Protects against ETA-leak when an
// ODE solve throws/returns NA mid finite-difference perturbation, which
// would otherwise leave a perturbed ETA in ind->par_ptr and corrupt later calls.
struct EtaRestoreGuard {
  rx_solving_options_ind *ind;
  arma::vec saved;
  bool armed;
  EtaRestoreGuard(int cid) : armed(true) {
    ind = getSolvingOptionsInd(rx, getRxId(cid));
    saved.set_size(op_focei.neta);
    for (int i = op_focei.neta; i--;) {
      saved[i] = getIndParPtr(ind, op_focei.etaTrans[i]);
    }
  }
  ~EtaRestoreGuard() {
    if (!armed) return;
    for (int i = op_focei.neta; i--;) {
      setIndParPtr(ind, op_focei.etaTrans[i], saved[i]);
    }
  }
  void disarm() { armed = false; }
};

// Per-subject "did this solve fail" check, replacing the shared
// hasOpBadSolve(op) poll: op->badSolve is set globally on ANY thread's
// failure, so reading it in a retry loop is racy across threads. Each
// subject's own ind->solve buffer holds NA_REAL on failure, so scanning it
// answers "did THIS subject fail" without touching shared state.
static inline bool indHasBadSolve(rx_solving_options *op,
                                  rx_solving_options_ind *ind) {
  int neq = getOpNeq(op);
  if (neq <= 0) return false;
  double *solve = getIndSolve(ind);
  int n = neq * getIndNallTimes(ind);
  for (int i = 0; i < n; ++i) {
    if (ISNA(solve[i]) || std::isnan(solve[i]) || std::isinf(solve[i])) {
      return true;
    }
  }
  return false;
}

// RAII guard for ind->neqOverride: switches one subject's effective neq for a
// predOde finite-difference pass without mutating shared op->neq (rxode2
// honors it via rxEffNeq(ind, op)), so other threads keep a stable state count.
struct IndNeqOverrideGuard {
  rx_solving_options_ind *ind;
  int saved;
  bool armed;
  IndNeqOverrideGuard(rx_solving_options_ind *_ind, int newOverride)
      : ind(_ind), armed(true) {
    saved = getIndNeqOverride(ind);
    setIndNeqOverride(ind, newOverride);
  }
  ~IndNeqOverrideGuard() {
    if (armed) setIndNeqOverride(ind, saved);
  }
  void disarm() { armed = false; }
};

arma::vec getCurEta(int cid) {
  rx_solving_options_ind *ind =  getSolvingOptionsInd(rx, getRxId(cid));
  arma::vec eta(op_focei.neta);
  for (int i = op_focei.neta; i--;) {
    eta[i] = getIndParPtr(ind, op_focei.etaTrans[i]);
  }
  return eta;
}

arma::mat grabRFmatFromInner(int id, bool predSolve) {
  int _rxId = getRxId(id); // base subject index for rxode2
  rx_solving_options_ind *ind =  getSolvingOptionsInd(rx, _rxId);
  focei_ind *fInd = &(inds_focei[id]);
  arma::vec retF(getIndNallTimes(ind));
  arma::vec retR(getIndNallTimes(ind));
  // this assumes the inner problem has been solved
  fInd->nObs = 0;
  rx_solving_options *op = getSolvingOptions(rx);
  int kk, k=0;
  double curT;
  if (predSolve) {
    iniSubjectE(_rxId, 1, ind, op, rx, rxPred.update_inis);
  } else {
    iniSubjectE(_rxId, 1, ind, op, rx, rxInner.update_inis);
  }
  iniSubjectE(_rxId, 1, ind, op, rx, rxPred.update_inis);
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      if (predSolve) {
        rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      } else {
        rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      }
      continue;
    }
    fInd->nObs++;
    if (predSolve) {
      rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      retF(k) = lhs[op_focei.predNoLhsOffset];
      retR(k) = lhs[op_focei.predNoLhsOffset + 1];
    } else {
      rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      retF(k) = lhs[op_focei.predOffset];
      retR(k) = lhs[op_focei.predOffset + op_focei.neta + 1];
    }
    k++;
    if (k >= getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind)) {
      // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
      break;
    }
  }
  arma::mat ret(fInd->nObs, 2);
  ret.col(0) = retF(span(0, fInd->nObs-1));
  ret.col(1) = retR(span(0, fInd->nObs-1));
  return ret;
}

// This is needed for shi21 h optimization
arma::vec shi21EtaGeneral(arma::vec &eta, int id, int w) {
  EtaRestoreGuard etaGuard(id); // restores ind->par_ptr on any exit path
  updateEta(eta.memptr(), id);
  focei_ind *fInd = &(inds_focei[id]);
  arma::vec ret(fInd->nObs);
  int _rxId = getRxId(id); // base subject index for rxode2 (only nSub subjects)
  rx_solving_options_ind *ind =  getSolvingOptionsInd(rx, _rxId);
  rx_solving_options *op = getSolvingOptions(rx);
  IndNeqOverrideGuard neqGuard(ind, op_focei.predNeq); // switches this subject's neq to predNeq
  predOde(_rxId); // Assumes same order of parameters; use base subject index
  int kk, k = 0;
  iniSubjectE(_rxId, 1, ind, op, rx, rxPred.update_inis);
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    }
    rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
    ret(k) = lhs[w];
    k++;
    if (k >= getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind)) {
      // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
      break;
    }
  }
  // neqGuard restores ind->neqOverride; etaGuard restores ind->par_ptr.
  return ret;
}

arma::vec shi21EtaF(arma::vec &eta, int id) {
  return shi21EtaGeneral(eta, id, 0);
}

arma::vec shi21EtaR(arma::vec &eta, int id) {
  return shi21EtaGeneral(eta, id, 1);
}

// est="impmap" gradient FD fallback: perturb the THETAS (holding the current eta,
// set by the caller) and solve the pred model, returning per-observation f (w=0,
// rx_pred_) or V (w=1, rx_r_).  Mirrors shi21EtaGeneral but over thetas; the
// caller restores fullTheta afterward.  Used by shi21Central()/shi21Forward().
arma::vec shi21ThetaGeneral(arma::vec &theta, int id, int w) {
  focei_ind *fInd = &(inds_focei[id]);
  int _rxId = getRxId(id);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, _rxId);
  rx_solving_options *op = getSolvingOptions(rx);
  for (int t = 0; t < (int)op_focei.ntheta; ++t) {
    setIndParPtr(ind, op_focei.thetaTrans[t], theta[t]);
  }
  IndNeqOverrideGuard neqGuard(ind, op_focei.predNeq);
  setIndSolve(ind, -1);
  predOde(_rxId);
  int nObs = getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind);
  fInd->nObs = nObs; // keep consistent for callers
  arma::vec ret(nObs);
  iniSubjectE(_rxId, 1, ind, op, rx, rxPred.update_inis);
  int kk, k = 0;
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    }
    rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
    if (k < nObs) ret(k) = lhs[op_focei.predNoLhsOffset + w];
    k++;
    if (k >= nObs) break;
  }
  return ret;
}
arma::vec shi21ThetaF(arma::vec &theta, int id) { return shi21ThetaGeneral(theta, id, 0); }
arma::vec shi21ThetaR(arma::vec &theta, int id) { return shi21ThetaGeneral(theta, id, 1); }

arma::vec calcGradForward(arma::vec &f0,
                          arma::vec &grPH,  double h) {
  if (grPH.is_finite()) {
    // forward
    return (grPH - f0)/h;
  }
  arma::vec ret(grPH.size());
  ret.zeros();
  return ret;
}

arma::vec calcGradCentral(arma::vec &grMH, arma::vec &f0,
                          arma::vec &grPH,  double h) {
  if (grMH.is_finite() && grPH.is_finite()) {
    return (grPH-grMH)/(2.0*h);
  } else if (grPH.is_finite()) {
    // forward
    return (grPH - f0)/(h);
  } else if (grMH.is_finite()) {
    // backward
    return (f0 - grMH)/(h);
  }
  arma::vec ret(grPH.size());
  ret.zeros();
  return ret;
}
// FOCE (interaction==0): R must enter the inner likelihood at the eta=0
// population prediction, held constant (the truncated Sheiner-Beal gradient
// drops dR/deta, so eta-dependent R destabilizes the optimizer).  Solve the
// inner model at eta=0 and read rx_r_ (lhs[neta+1]) per obs into rPop in the
// likInner0 observation k-order.  Call BEFORE the inner solve (this overwrites
// ind->solve; the inner solve re-establishes it).
static void getPopR(int id, arma::vec &rPop) {
  int _rxId = getRxId(id); // base subject index for rxode2
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, _rxId);
  rx_solving_options *op = getSolvingOptions(rx);
  int nObsMax = getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind);
  if ((int)rPop.n_elem < nObsMax) rPop.set_size(nObsMax);
  EtaRestoreGuard etaGuard(id); // restore the trial eta on any exit path
  arma::vec zeroEta(op_focei.neta, arma::fill::zeros);
  updateEta(zeroEta.memptr(), id); // eta = 0 -> population prediction
  // innerOde (not predOde): the pred model shares a global linCmt solver
  // pointer with the inner model, so solving it here would corrupt the next
  // inner linCmt gradient.  At eta=0 the inner model's states are population,
  // so its rx_r_ is the genuine eta=0 R for both ODE and linCmt.
  setIndSolve(ind, -1);
  innerOde(_rxId); // solve the inner model at eta=0
  iniSubjectE(_rxId, 1, ind, op, rx, rxInner.update_inis);
  int kk, k = 0;
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
      rPop(k) = lhs[op_focei.neta + 1]; // inner-model rx_r_ at eta=0
      k++;
      if (k >= nObsMax) break;
    }
  }
}

double likInner0(double *eta, int id) {
  rx = getRxSolve_();
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
  rx_solving_options *op = getSolvingOptions(rx);
  int i, j;
  bool recalc = false;
  focei_ind *fInd= &(inds_focei[id]);
  op_focei.didLikCalc.store(true, std::memory_order_relaxed);
  double *solve = getIndSolve(ind);
  if (op_focei.neta > 0){
    if (!fInd->setup){
      recalc = true;
      fInd->setup = 1;
    } else {
      // Check to see if old ETA matches.
      for (j = op_focei.neta; j--;){
        if (fInd->oldEta[j] != eta[j]){
          recalc = true;
          break;
        }
      }
    }
  } else {
    recalc = true;
  }
  arma::vec rPopVec;  // FOCE: per-obs eta=0 population R (empty for FOCEI)
  if (recalc){
    if (op_focei.mixIdxN != 0) {
      setIndMixest(ind, getRxMixFromId(id));
    }
    for (j = op_focei.neta; j--;){
      setIndParPtr(ind, op_focei.etaTrans[j], eta[j]);
    }
    // FOCE: capture eta=0 population R before the inner solve overwrites
    // ind->solve.  rPop is a function of theta only, so it is cached across inner
    // iterations and recomputed only when updateTheta() bumps the generation
    // counter -- keeping FOCE at one solve per inner iteration (like FOCEI) plus
    // one eta=0 solve per subject per outer iteration, not two solves every time.
    // Only "nonmem" FOCE freezes R at the eta=0 population value; "foce+"
    // (foceType==1) keeps the live conditional R and needs no eta=0 solve.
    if (op_focei.interaction == 0 && op_focei.neta > 0 && op_focei.fo == 0 &&
        op_focei.foceType == 0) {
      if (id >= 0 && id < (int)_foceRPopGen.size()) {
        if (_foceRPopGen[id] != _foceRPopCurGen) {
          getPopR(id, _foceRPopCache[id]);
          _foceRPopGen[id] = _foceRPopCurGen;
        }
        rPopVec = _foceRPopCache[id];
      } else {
        getPopR(id, rPopVec); // cache not sized (defensive); recompute directly
      }
    }
    // Reset the sticky-recalc counter only if this subject hasn't exhausted its
    // retry budget yet; otherwise subsequent calls keep short-circuiting the
    // retry loop since tolFactor is already permanently loose.
    if (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN) {
      fInd->stickyRecalcN2 = 0;
    }
    setIndSolve(ind, -1);
    // op->badSolve is racy under cores>1; our retry loop reads per-subject
    // indHasBadSolve() instead, this reset is just a courtesy for diagnostics.
    resetOpBadSolve(op);
    bool predSolve = false;
    // Guard for ind->neqOverride; lazily allocated in the doFD branch, lives
    // until likInner0 returns so all reads of ind->solve see the same predNeq
    // stride. Stays nullptr (no-op) on the common innerOde path.
    std::unique_ptr<IndNeqOverrideGuard> neqGuard;
    // Mixture subjects (id >= nSub) use the base subject index for rxode2 calls.
    int _rxId = getRxId(id);
    if (fInd->doFD == 0) {
      double prevTol = getIndTolFactor(ind);
      innerOde(_rxId);
      j = 0;
      while (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN
             && indHasBadSolve(op, ind) && j < op_focei.maxOdeRecalc) {
        fInd->stickyRecalcN2++;
        op_focei.reducedTol.store(1, std::memory_order_relaxed);
        op_focei.reducedTol2.store(1, std::memory_order_relaxed);
        atolRtolFactor_(op_focei.odeRecalcFactor);
        setIndSolve(ind, -1);
        resetOpBadSolve(op);
        innerOde(_rxId);
        j++;
      }
      if (j != 0) {
        if (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN) {
          setIndTolFactor(ind, prevTol); // succeeded within budget; un-stick
        } else {
          // Hit sticky threshold: tolFactor stays loose, latch global stickyTol.
          op_focei.stickyTol.store(1, std::memory_order_relaxed);
        }
      }
    } else {
      // Inner sensitivity solve failed; fall back to perturbing the simpler
      // prediction model, keeping the neq override alive through likInner0's
      // remaining reads of ind->solve.
      neqGuard.reset(new IndNeqOverrideGuard(ind, op_focei.predNeq));
      predOde(_rxId);
      predSolve=true;
      op_focei.didPredSolve.store(true, std::memory_order_relaxed);
    }
    bool isBadSolve = false;
    // predSolve lays the buffer out at predNeq stride; scan only those slots.
    int effNeq = predSolve ? op_focei.predNeq : getOpNeq(op);
    int nsolve = (effNeq + getOpNlin(op))*getIndNallTimes(ind);
    if (effNeq > 0) {
      for (int ns = 0; ns < nsolve; ++ns) {
        if (ISNA(solve[ns]) || std::isnan(solve[ns]) ||
            std::isinf(solve[ns])) {
          isBadSolve = true;
          break;
        }
      }
    }
    if (isBadSolve){
      return NA_REAL;
      //throw std::runtime_error("bad solve");
    } else {
      // Update eta.
      arma::mat lp(fInd->lp, op_focei.neta, 1, false, true);
      lp.zeros();
      arma::mat a(fInd->a, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind), op_focei.neta, false, true);
      arma::mat B(fInd->B, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind), 1, false, true);
      arma::mat c(fInd->c, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind),
                  op_focei.neta, false, true);
      arma::mat Vid(fInd->Vid, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind),
                    getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind), false, true);
      // a.zeros(); already reset
      // B.zeros(); already reset
      // c.zeros(); already reset
      // Vid.zeros(); // already reset when needed
      // Check to see if finite difference step size needs to be optimized
      bool finiteDiffNeeded = predSolve;
      for (int ii = 0; ii < op_focei.neta; ++ii) {
        if (op_focei.etaFD[ii]==1) {
          finiteDiffNeeded = true;
        }
        if (finiteDiffNeeded) break;
      }
      arma::mat etaGradF(fInd->a, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind),
                         op_focei.neta, false, true);
      etaGradF.zeros();
      arma::mat etaGradR(fInd->c, getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind),
                         op_focei.neta, false, true);
      etaGradR.zeros();
      if (finiteDiffNeeded) {
        // need to optimize finite difference
        // First get the f0 for F and R based on current solve
        arma::mat rf0mat = grabRFmatFromInner(id, predSolve);
        // now save the prior solve
        arma::vec solveSave(nsolve);
        std::copy(solve, solve + nsolve, solveSave.memptr());
        arma::vec f0 = rf0mat.col(0);
        arma::vec r0 = rf0mat.col(1);
        arma::vec curEta = getCurEta(id);
        arma::vec hEta(curEta.size());
        arma::vec grETA(fInd->nObs);

        arma::vec grPH(fInd->nObs);
        arma::vec grMH(fInd->nObs);

        for (int ii = 0; ii < op_focei.neta; ++ii) {
          if (predSolve || op_focei.etaFD[ii]==1) {
            if (fInd->etahf[ii] == 0.0) {
              double h = 0;
              //     .eventTypeIdx <- c("stencil" = 1L, "central" = 2L, "forward" = 3L)
              switch(op_focei.eventType) {
              case 2: // central
                fInd->etahf[ii] = shi21Central(shi21EtaF, curEta, h,
                                               f0, grETA, id, ii,
                                               op_focei.hessEpsInner, // ef,
                                               1.5,//double rl = 1.5,
                                               4.5,//double ru = 4.5,
                                               3.0,//double nu = 8.0);
                                               op_focei.shi21maxFD); // maxiter
                break;
              case 3: // forward
                fInd->etahf[ii] = shi21Forward(shi21EtaF, curEta, h,
                                               f0, grETA, id, ii,
                                               op_focei.hessEpsInner, // ef,
                                               1.5,  //double rl = 1.5,
                                               6.0,  //double ru = 6.0);;
                                               op_focei.shi21maxFD); // maxiter
              }
              etaGradF.col(ii) = grETA;
              if (op_focei.interaction == 1) {
                switch(op_focei.eventType) {
                case 2: //central
                  fInd->etahr[ii] = shi21Central(shi21EtaR, curEta, h,
                                                 r0, grETA, id, ii,
                                                 op_focei.hessEpsInner, // ef,
                                                 1.5,//double rl = 1.5,
                                                 4.5,//double ru = 4.5,
                                                 3.0,//double nu = 8.0);
                                                 op_focei.shi21maxFD); // maxiter
                  break;
                case 3: // forward
                  fInd->etahr[ii] = shi21Forward(shi21EtaR, curEta, h,
                                                 r0, grETA, id, ii,
                                                 op_focei.hessEpsInner, // ef,
                                                 1.5,  //double rl = 1.5,
                                                 6.0,  //double ru = 6.0);;
                                                 op_focei.shi21maxFD); // maxiter
                  break;
                }
                etaGradR.col(ii) = grETA;
              }
            } else {
              // etaGradF
              hEta = curEta;
              hEta[ii] += fInd->etahf[ii];
              grPH = shi21EtaF(hEta, id);
              bool useForward = false;
              if (op_focei.eventType == 3) {
                // if this isn't true try backward
                if (grPH.is_finite()) {
                  useForward = true;
                  etaGradF.col(ii) = calcGradForward(f0, grPH,  fInd->etahf[ii]);
                }
              }
              if (!useForward) {
                // stencil or central
                hEta = curEta;
                hEta[ii] -= fInd->etahf[ii];
                grMH = shi21EtaF(hEta, id);
                // central
                etaGradF.col(ii) = calcGradCentral(grMH, f0, grPH,  fInd->etahf[ii]);
              }
              if (op_focei.interaction == 1) {
                // etaGradR
                hEta = curEta;
                hEta[ii] += fInd->etahr[ii];
                grPH = shi21EtaR(hEta, id);
                useForward = false;
                if (op_focei.eventType == 3) {
                  if (grPH.is_finite()) {
                    useForward = true;
                    etaGradR.col(ii) = calcGradForward(r0, grPH,  fInd->etahr[ii]);
                  }
                }
                if (!useForward) {
                  hEta = curEta;
                  hEta[ii] -= fInd->etahr[ii];
                  grMH = shi21EtaR(hEta, id);
                  etaGradR.col(ii) = calcGradCentral(grMH, r0, grPH,  fInd->etahr[ii]);
                }
              }
            }
          }
        }
        // restore the prior solve
        std::copy(solveSave.begin(), solveSave.end(), solve);
      }
      if (op_focei.fo == 1){
        Vid.zeros();
      }
      // RSprintf("ID: %d; Solve #2: %f\n", id, solve[2]);
      // Calculate matricies
      int k = 0, kk=0;//getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind) - 1;
      fInd->llik=0.0;
      fInd->nNonNormal = 0;
      fInd->nObs = 0;
      fInd->tbsLik=0.0;
      double f, err, r, fpm, rp = 0,lnr, limit, dv,dv0, curT;
      int cens = 0;
      if (predSolve) {
        iniSubjectE(_rxId, 1, ind, op, rx, rxPred.update_inis);
      } else {
        iniSubjectE(_rxId, 1, ind, op, rx, rxInner.update_inis);
      }
      int dist=0, yj0=0, yj = 0;
      double *llikObs = fInd->llikObs;
      for (j = 0; j < getIndNallTimes(ind); ++j) {
        setIndIdx(ind, j);
        kk = getIndIx(ind, j);
        curT = getTime(kk, ind);
        dv0 = getIndDv(ind, kk);
        yj = getIndYj(ind);
        _splitYj(&yj, &dist,  &yj0);
        double *lhs = getIndLhs(ind);
        if (isDose(getIndEvid(ind, kk))) {
          llikObs[kk] = NA_REAL;
          // Need to calculate for advan sensitivities
          if (predSolve) {
            rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
            // Normalize the predNoLhs layout into the inner offset layout so the
            // shared reads below use op_focei.predOffset uniformly.
            double _pf = lhs[op_focei.predNoLhsOffset];
            double _pr = lhs[op_focei.predNoLhsOffset + 1];
            lhs[op_focei.predOffset] = _pf;
            lhs[op_focei.predOffset + op_focei.neta + 1] = _pr;
          }
          else {
            rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
          }
        } else if (getIndEvid(ind, kk) == 0) {
          if (predSolve) {
            rxPred.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
            // Normalize the predNoLhs layout into the inner offset layout.
            double _pf = lhs[op_focei.predNoLhsOffset];
            double _pr = lhs[op_focei.predNoLhsOffset + 1];
            lhs[op_focei.predOffset] = _pf;
            lhs[op_focei.predOffset + op_focei.neta + 1] = _pr;
          } else {
            rxInner.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, j), lhs);
          }

          f = lhs[op_focei.predOffset]; // TBS is performed in the rxode2 rx_pred_ statement. This allows derivatives of TBS to be propagated
          dv = tbs(dv0);
          if (ISNA(f) || std::isnan(f) || std::isinf(f)) {
            return NA_REAL;
            //throw std::runtime_error("bad solve");
          }
          // fInd->f(k, 0) = lhs[0];
          err = f - dv;
          limit = R_NegInf;
          if (hasRxLimit(rx)) {
            limit = getIndLimit(ind, kk);
            if (ISNA(limit)) {
              limit = R_NegInf;
            } else if (R_FINITE(limit)) {
              limit = tbs(limit);
            }
          }
          cens = 0;
          if (hasRxCens(rx)) cens = getIndCens(ind, kk);
          fInd->tbsLik+=tbsL(dv0);
          // fInd->err(k, 0) = lhs[0] - getIndDv(ind, k); // pred-dv
          if (ISNA(lhs[op_focei.predOffset + op_focei.neta + 1])){
            return NA_REAL;
            //throw std::runtime_error("bad solve");
          }
          if (dist == rxDistributionNorm) {
            r = lhs[op_focei.predOffset + op_focei.neta + 1];
            // "nonmem" FOCE: use the eta=0 population R (FOCEI and "foce+" keep
            // the live inner rx_r_ evaluated at the current eta)
            if (op_focei.interaction == 0 && op_focei.neta > 0 && op_focei.fo == 0 &&
                op_focei.foceType == 0) {
              r = rPopVec(k);
            }
            if (r <= sqrt(std::numeric_limits<double>::epsilon())) {
              r = 1.0;
            }
          } else {
            r = 1.0;
          }
          if (op_focei.neta == 0) {
            if (dist == rxDistributionNorm) {
              double ll = err/r;
              //r = variance
              ll =  -0.5 * ll * err - 0.5*log(r);
              ll = doCensNormal1((double)cens, dv, limit, ll, f, r,
                                 (int)op_focei.adjLik);
              llikObs[kk] = ll;
              fInd->llik += ll;
              fInd->nObs++;
            } else {
              llikObs[kk] = f;
              fInd->llik += f;
              fInd->nNonNormal++;
              fInd->nObs++;
            }
          } else if (op_focei.fo == 1) {
            // FO
            B(k, 0) = err; // res
            Vid(k, k) = r;
            for (i = op_focei.neta; i--; ) {
              if (predSolve || op_focei.etaFD[i]==1) {
                a(k, i) = etaGradF(k, i);
              } else {
                a(k, i) = lhs[op_focei.predOffset + i + 1];
              }
            }
            // Ci = fpm %*% omega %*% t(fpm) + Vi; Vi=diag(r)
          } else {
            // Use `a` (=fpm) for the gradient so dose-based etas share the same
            // approach for normal and non-normal log likelihoods; err/r are garbage here.
            // FOCE: log(R) from `r`, which already holds the mode-appropriate R
            // (eta=0 frozen for "nonmem", live conditional for "foce+"); FOCEI keeps lhs
            if (dist == rxDistributionNorm) {
              if (op_focei.interaction == 0 && op_focei.neta > 0 && op_focei.fo == 0) {
                lnr = _safe_log(r);
              } else {
                lnr = _safe_log(lhs[op_focei.predOffset + op_focei.neta + 1]);
              }
            }
            else lnr = 0;
            // fInd->r(k, 0) = lhs[op_focei.neta+1];
            // B(k, 0) = 2.0/lhs[op_focei.neta+1];
            // lhs 0 = F
            // lhs 1-eta = df/deta
            // FIXME faster initialization via copy or elm
            // RSprintf("id: %d k: %d j: %d\n", id, k, j);
            B(k, 0) = 2.0/_safe_zero(r);
            // per-obs censored inner-Hessian coefficients (rho_ff/rho_fR/rho_RR) for the
            // exact censored Laplace determinant; normal obs get the Gauss-Newton values
            // (1/r, 0, 0.5/r^2) so calcEtaHessian reduces bit-identically to the old form.
            {
              int isCensObs = (cens != 0) || (R_FINITE(limit) && !ISNA(limit));
              if (isCensObs && dist == rxDistributionNorm && op_focei.censOption == 1) {
                double _cp[9]; for (int _i = 0; _i < 9; _i++) _cp[_i] = 0.0;
                censNormalPartials((double)cens, dv, limit, f, r, 2, _cp);
                fInd->cHff[k] = _cp[2]; fInd->cHfr[k] = _cp[3]; fInd->cHrr[k] = _cp[4];
              } else {
                // gauss (historic) or a normal obs: uncensored Gauss-Newton curvature
                fInd->cHff[k] = 1.0/_safe_zero(r); fInd->cHfr[k] = 0.0;
                fInd->cHrr[k] = 0.5/_safe_zero(r*r);
              }
            }
            if (op_focei.interaction == 1) {
              for (i = op_focei.neta; i--; ) {
                if (predSolve || op_focei.etaFD[i]==1) {
                  fpm = a(k, i) = etaGradF(k, i);
                  rp = 0.0;
                  if (dist == rxDistributionNorm) {
                    rp = etaGradR(k, i);
                  }
                } else {
                  fpm = a(k, i) = lhs[op_focei.predOffset + i + 1]; // Almquist uses different a (see eq #15)
                  rp  = (dist == rxDistributionNorm)*lhs[op_focei.predOffset + i + op_focei.neta + 2];
                }
                if (fpm == 0.0) {
                  a(k, i) = fpm = sqrt(DBL_EPSILON);
                }
                if (rp == 0.0) {
                  rp = sqrt(DBL_EPSILON);
                }
                c(k, i) = rp/_safe_zero(r);
                //lp is eq 12 in Almquist 2015
                // .5*apply(eps*fp*B + .5*eps^2*B*c - c, 2, sum) - OMGAinv %*% ETA
                if (dist == rxDistributionNorm) {
                  double lpCur = 0.25 * err * err * B(k, 0) * c(k, i) -
                    0.5 * c(k, i) - 0.5 * err * fpm * B(k, 0);
                  lp(i, 0) += dCensNormal1((double)cens, dv, limit, lpCur, f, r, fpm, rp);
                } else {
                  lp(i, 0) += fpm;
                }
              }
              // Eq #10
              //llik <- -0.5 * sum(err ^ 2 / R + log(R));
               if (dist == rxDistributionNorm) {
                 double ll = -0.5 * err * err/_safe_zero(r) - 0.5 * lnr;
                 ll = doCensNormal1((double)cens, dv, limit, ll, f, r, (int) op_focei.adjLik);
                 llikObs[kk] = ll;
                 fInd->llik +=  ll;
                 fInd->nObs++;
               } else {
                 llikObs[kk] = f;
                 fInd->llik += f;
                 fInd->nNonNormal++;
                 fInd->nObs++;
               }
            } else if (op_focei.interaction == 0) {
              for (i = op_focei.neta; i--; ){
                if (predSolve || op_focei.etaFD[i]==1) {
                  a(k, i) = fpm = etaGradF(k, i);
                } else {
                  a(k, i) = fpm = lhs[op_focei.predOffset + i + 1];
                }
                if (dist == rxDistributionNorm) {
                  double lpCur = -0.5 * err * fpm * B(k, 0);
                  lp(i, 0) += dCensNormal1((double)cens, dv, limit, lpCur, f, r, fpm, rp);
                } else {
                  lp(i, 0) += fpm;
                }
              }
              // Eq #10
              //llik <- -0.5 * sum(err ^ 2 / R + log(R));
              if (dist == rxDistributionNorm) {
                double ll = -0.5 * err * err/_safe_zero(r) -0.5 * lnr;
                ll =doCensNormal1((double)cens, dv, limit, ll, f, r, (int) op_focei.adjLik);
                llikObs[kk] = ll;
                fInd->llik +=  ll;
                fInd->nObs++;
              } else {
                llikObs[kk] = f;
                fInd->llik += f;
                fInd->nNonNormal++;
                fInd->nObs++;
              }
            }
          }
          // k--;
          k++;
          if (k >= getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind)) {
            // With moving doses this may be at the very end, so drop out now if all the observations were accounted for
            break;
          }
        }
      }
      if (op_focei.neta == 0) {
        if (fInd->nNonNormal && op_focei.adjLik) {
          fInd->llik -= fInd->nNonNormal*M_LN_SQRT_2PI;
        }
      } else if (op_focei.fo == 1) {
        if (cens != 0) stop("FO censoring not supported.");
        if (dist != rxDistributionNorm) stop("Generalized llik for FO is not supported");
        mat Ci = a * op_focei.omega * trans(a) + Vid;
        mat cholCi = cholSE__(Ci, op_focei.cholSEtol);
        mat CiInv;
        bool success  = inv(CiInv, trimatu(cholCi));
        if (!success){
          CiInv = pinv(trimatu(cholCi));
        }
        CiInv = CiInv * CiInv.t();
        double lik = 0;
        // 2*sum(log(diag(chol(Ci))))
        for (unsigned int j = cholCi.n_rows; j--;){
          lik += 2*_safe_log(cholCi(j,j));
        }
        // + t(.$Ri) %*% solve(Ci) %*% .$Ri
        mat rest =trans(B) * CiInv * B;
        lik += rest(0,0);
        // lik = -2*ll
        fInd->llik = -0.5*lik;
      } else {
        // Now finalize lp
        mat etam = arma::mat(op_focei.neta, 1);
        std::copy(&eta[0], &eta[0] + op_focei.neta, etam.begin()); // fill in etam
        // print(wrap(lp));
        // print(wrap(etam));
        // Finalize eq. #12
        lp = -(lp - op_focei.omegaInv * etam);
        // Partially finalize #10
        fInd->llik = -trace(fInd->llik - 0.5*(etam.t() * op_focei.omegaInv * etam));
        if (fInd->nNonNormal && op_focei.adjLik) {
          fInd->llik -= fInd->nNonNormal*M_LN_SQRT_2PI;
        }
        std::copy(&eta[0], &eta[0] + op_focei.neta, &fInd->oldEta[0]);
      }
    }
  }
  return fInd->llik;
}

double *lpInner(double *eta, double *g, int id){
  focei_ind *fInd = &(inds_focei[id]);
  likInner0(eta, id);
  std::copy(&fInd->lp[0], &fInd->lp[0] + op_focei.neta,
            &g[0]);
  return &g[0];
}

//[[Rcpp::export]]
NumericVector foceiInnerLp(NumericVector eta, int id = 1){
  double *etad = new double[eta.size()];
  std::copy(eta.begin(),eta.end(),&etad[0]);
  NumericVector lp(eta.size());
  int curId = id - 1;
  lpInner(etad,&lp[0], curId);
  delete[] etad;
  return lp;
}

//[[Rcpp::export]]
double likInner(NumericVector eta, int id = 1){
  double *etad = new double[eta.size()];
  std::copy(eta.begin(),eta.end(),&etad[0]);
  int curId = id - 1;
  double llik = likInner0(etad, curId);
  delete[] etad;
  return llik;
}

// This is needed for shi21 h optimization
arma::vec getGradForOptimHess(arma::vec &t, int id) {
  arma::vec ret(op_focei.neta);
  lpInner(&t[0], &ret[0], id);
  return ret;
}

bool _finalObfCalc = false;

bool calcEtaHessian(double *eta, int likId, int id,
                    focei_ind *fInd,
                    rx_solving_options_ind *ind,
                    mat &H, mat &H0) {
  H.zeros();
  int k, l;
  // This is actually -H
  if (op_focei.needOptimHess) {
    arma::vec gr0(op_focei.neta, fill::zeros);
    std::copy(&fInd->lp[0], &fInd->lp[0] + op_focei.neta, &gr0[0]);

    arma::vec grPH(op_focei.neta, fill::zeros);
    arma::vec grMH(op_focei.neta, fill::zeros);

    double h = 0;

    for (k = op_focei.neta; k--;) {
      h = fInd->etahh[k];
      if (op_focei.optimHessType == 3 && h <= 0) {
        arma::vec t(eta, op_focei.neta);
        fInd->etahh[k] = shi21Forward(getGradForOptimHess, t, h,
                                      gr0, grPH, id, k,
                                      op_focei.hessEpsInner, //double ef = 7e-7,
                                      1.5,  //double rl = 1.5,
                                      6.0,  //double ru = 6.0);;
                                      op_focei.shi21maxInner);  //maxiter=15
        H.col(k) = grPH;
        continue;
      }
      if (op_focei.optimHessType == 1 && h <= 0) {
        // Central
        arma::vec t(eta, op_focei.neta);
        fInd->etahh[k] = shi21Central(getGradForOptimHess, t, h,
                                      gr0, grPH, id, k,
                                      op_focei.hessEpsInner, // ef,
                                      1.5,//double rl = 1.5,
                                      4.5,//double ru = 4.5,
                                      3.0,//double nu = 8.0);
                                      op_focei.shi21maxInner); // maxiter
        H.col(k) = grPH;
        continue;
      }
      // x + h
      eta[k] += h;
      lpInner(eta, &grPH[0], id);
      bool forwardFinite =  grPH.is_finite();
      if (op_focei.optimHessType == 3 && forwardFinite) { // forward
        H.col(k) = (grPH-gr0)/h;
        eta[k] -= h;
        continue;
      }

      // x - h
      eta[k] -= 2*h;
      lpInner(eta, &grMH[0], id);
      bool backwardFinite = grMH.is_finite();
      if (op_focei.optimHessType == 1 &&
          forwardFinite && backwardFinite) {
        // central
        eta[k] += h;
        H.col(k) = (grPH-grMH)/(2.0*h);
        continue;
      }
      if (forwardFinite && !backwardFinite) {
        // forward difference
        H.col(k) = (grPH-gr0)/h;
        eta[k] += h;
        continue;
      }
      if (!forwardFinite && backwardFinite) {
        // backward difference
        H.col(k) = (gr0-grMH)/h;
        eta[k] += h;
        continue;
      }
    }
    // symmetrize
    H = 0.5*(H + H.t());
    // Note that since the gradient includes omegaInv*etam,
    // op_focei.omegaInv(k, l) shouldn't be added.
  } else if (op_focei.interaction) {
    int nO = getIndNallTimes(ind) - getIndNdoses(ind) - getIndNevid2(ind);
    arma::mat a(fInd->a, nO, op_focei.neta, false, true);
    arma::mat B(fInd->B, nO, 1, false, true);
    arma::mat c(fInd->c, nO, op_focei.neta, false, true);
    // per-obs exact-Laplace inner-Hessian coefficients (rho_ff/rho_fR/rho_RR); normal
    // obs carry the Gauss-Newton 1/r, 0, 0.5/r^2 so this reduces bit-identically to the
    // old 0.5*(a*B*a + c*c) form.  rp = dr/deta = c*r = c*(2/B).
    arma::vec cHff(fInd->cHff, nO, false, true), cHfr(fInd->cHfr, nO, false, true), cHrr(fInd->cHrr, nO, false, true);
    arma::vec r2 = 2.0 / B.col(0);   // r per obs (= 2/B)
    // A non-observation row counted in nO but never filled by the per-obs loop
    // (a/B/cHff/cHfr/cHrr all left at their zero-initialized default -- e.g. a
    // DDE model's extra time-zero bookkeeping row) has B==0, so 2/B is +Inf;
    // that Inf then turns a legitimately-zero cHfr/cHrr*rp term into 0*Inf=NaN,
    // corrupting the whole sum even though the row should contribute nothing.
    // Zero it explicitly -- r2 is only ever multiplied by other per-row
    // coefficients that are themselves 0 for such a row, so this changes
    // nothing for any row with a real (nonzero) B.
    r2.elem(arma::find_nonfinite(r2)).zeros();
    for (k = op_focei.neta; k--;){
      arma::vec rpk = c.col(k) % r2;
      for (l = k+1; l--;){
        arma::vec rpl = c.col(l) % r2;
        H(k, l) = sum(cHff % a.col(l) % a.col(k) +
                      cHfr % (a.col(l) % rpk + rpl % a.col(k)) +
                      cHrr % rpl % rpk) +
          op_focei.omegaInv(k, l);
        if (!R_finite(H(k, l))) {
          return false;
        }
        H(l, k) = H(k, l);
      }
    }
  } else {
    arma::mat a(fInd->a, fInd->nObs, op_focei.neta, false, true);
    // FOCE (no interaction): frozen variance -> only the prediction curvature enters
    // the inner Hessian; censored obs use rho_ff^cens (cHff), normal use 1/r.
    arma::vec cHff(fInd->cHff, fInd->nObs, false, true);
    for (k = op_focei.neta; k--;){
      for (l = k+1; l--;) {
        H(k, l) = sum(cHff % a.col(l) % a.col(k)) +
          op_focei.omegaInv(k, l);
        if (!R_finite(H(k, l))) {
            return false;
        }
        H(l, k) = H(k, l);
      }
    }
  }
  k=0;
  if (!H.is_sympd()) {
    arma::mat H2;
    if (nmNearPD(H2, H)) {
      H=H2;
    }
  }
  if (fInd->doChol) {
    arma::mat Hout, Hin = H;
    bool success = chol_sym(Hout, Hin);
    if (!success) {
      return false;
    }
    H0 = Hout;
  } else {
    H0=cholSE__(H, op_focei.cholSEtol);
  }
  H = H0.t() * H0;
  return true;
}

double LikInner2(double *eta, int likId, int id) {
  focei_ind *fInd = &(inds_focei[id]);
  double lik=0;
  if (op_focei.neta == 0) {
    lik = fInd->llik;
  } else if (op_focei.fo == 1){
    // Already almost completely calculated.
    lik = fInd->llik;
  } else {
    // print(wrap(op_focei.logDetOmegaInv5));
    lik = -likInner0(eta, id);
    // print(wrap(lik));
    rx = getRxSolve_();
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    rx_solving_options *op = getSolvingOptions(rx);
    double *solve = getIndSolve(ind);
    if (getOpNeq(op) > 0 && ISNA(solve[0])){
      //return 1e300;
      return NA_REAL;
    }
    // Calculate lik first to calculate components for Hessian
    // Hessian
    mat H(op_focei.neta, op_focei.neta, fill::zeros);
    mat H0(op_focei.neta, op_focei.neta, fill::zeros);

    if (!calcEtaHessian(eta, likId, id, fInd, ind, H, H0)) {
      return NA_REAL;
    }
    if (_finalObfCalc) {
      std::copy(H.begin(), H.end(),
                op_focei.gH + id*op_focei.neta*op_focei.neta);
    }
    if (op_focei.warm == 1 && fInd->zmH != NULL) {
      // Save the calculated Hessian to warm-start the next n1qn1 inner problem
      vec hPack = H.elem(lowerTri(H, true));
      std::copy(hPack.begin(), hPack.end(), fInd->zmH);
      std::copy(&eta[0], &eta[0] + op_focei.neta, fInd->zmEta);
      fInd->zmValid = 1;
    }
    // - sum(log(H.diag()));
    double logH0diag = 0.0;
    for (unsigned int j = H0.n_rows; j--;){
      logH0diag -= _safe_log(H0(j,j));
    }
    if (_aqn == 0) {
      lik += logH0diag + op_focei.logDetOmegaInv5;
    } else {
      // Adaptive Gaussian Hermite Quadrature expansion around the EBE mode;
      // the SE of the EBE estimate (inverse cholesky of the Hessian) drives
      // the expansion (not needed for the Laplace case).
      // https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature
      // Get x and w from the AQ
      arma::mat aqx(op_focei.aqx, _aqn, op_focei.neta, false, true);
      arma::mat aqw(op_focei.aqw, _aqn, op_focei.neta, false, true);

      // llik = -likInner0(eta, id);
      int curi = 0;
      double slik = 0;
      if (op_focei.aqfirst) {
        // Already calculated:
        // lik = -likInner0(eta, id);
        arma::vec w = aqw.row(curi).t();
        curi=1;
        lik += sum(log(w)); // x % x  = x^2; here x=0
        // can be factored out
        //lik += op_focei.logDetOmegaInv5;
        lik = max2(lik, op_focei.aqLow);
        lik = min2(lik, op_focei.aqHi);
        slik = exp(lik);
      }
      arma::mat Ginv_5(op_focei.neta, op_focei.neta, fill::zeros);
      if (_aqn > 1) {
        // The Ginv_5 only needs to be calculated for the
        // non-Laplace case
        bool success;
        try {
          success = inv(Ginv_5, trimatu(H0));
          if (!success) {
            success = inv(Ginv_5, H0);
            if (!success) {
              return NA_REAL;
            }
          }
        } catch (...) {
          success = inv(Ginv_5, H0);
          if (!success) {
            return NA_REAL;
          }
        }
        arma::vec etahat(eta, op_focei.neta);
        for (; curi < _aqn; curi++) {
          // Get the x and w for the current iteration
          arma::vec x = aqx.row(curi).t();
          arma::vec w = aqw.row(curi).t();
          arma::vec etaCur = etahat +  Ginv_5 * x;
          lik  = -likInner0(etaCur.memptr(), id);
          lik += sum(log(w) +  0.5 * x % x); // x % x  = x^2
          // Can be factored out
          //lik += op_focei.logDetOmegaInv5;
          lik = max2(lik, op_focei.aqLow);
          lik = min2(lik, op_focei.aqHi);
          slik += exp(lik);
        }
      }
      //lik = 0.5*op_focei.neta * M_LN2 + det_Ginv_5 + log(slik);
      lik = log(slik) + logH0diag + op_focei.logDetOmegaInv5;
    }
  }
  lik += fInd->tbsLik;
  if (likId == 0){
    fInd->lik[0] = lik;
    if (op_focei.neta != 0) std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, &fInd->saveEta[0]);
  } else {
    // Use Objective function for Numeric Gradient.
    fInd->lik[likId] = -2*lik;
  }
  return lik;
}

// warm="calc": seed n1qn1's zm with the eta Hessian calculated at the starting
// eta; when no Hessian was calculated at this eta, calculate it here.
void warmZm(focei_ind *fInd, int id) {
  std::fill(&fInd->zm[0], &fInd->zm[0] + op_focei.nzm, 0.0);
  int neta = op_focei.neta;
  bool etaMatch = fInd->zmValid == 1;
  if (etaMatch) {
    for (int j = neta; j--;) {
      if (fInd->zmEta[j] != fInd->eta[j]) {
        etaMatch = false;
        break;
      }
    }
  }
  // For FD Hessians (needOptimHess) this optimizes the shi21 step (etahh) at
  // the starting eta (e.g. eta=0) instead of the EBE mode, shifting downstream
  // Hessians by FD error; the warm-start speedup is worth that small shift.
  if (!etaMatch) {
    double f = likInner0(fInd->eta, id);
    if (!ISNA(f)) {
      rx = getRxSolve_();
      rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
      mat H(neta, neta, fill::zeros);
      mat H0(neta, neta, fill::zeros);
      if (calcEtaHessian(fInd->eta, 0, id, fInd, ind, H, H0)) {
        vec hPack = H.elem(lowerTri(H, true));
        std::copy(hPack.begin(), hPack.end(), fInd->zmH);
        std::copy(&fInd->eta[0], &fInd->eta[0] + neta, fInd->zmEta);
        fInd->zmValid = 1;
        etaMatch = true;
      }
    }
  }
  if (etaMatch) {
    // n1qn1 mode=2 reads the packed lower-triangle Hessian from zm
    std::copy(&fInd->zmH[0], &fInd->zmH[0] + neta*(neta+1)/2, &fInd->zm[0]);
    fInd->mode = 2;
    fInd->uzm = 1;
  } else {
    fInd->mode = 1;
    fInd->uzm = 1;
  }
}

extern "C" double innerOptimF(int n, double *x, void *ex){
  int *id = (int*)ex;
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve == 1) return NA_REAL;
  double f = likInner0(x, *id);
  if (ISNA(f)) {
    fInd->badSolve = 1;
  }
  fInd->nInnerF++;
  return f;
}

extern "C" void innerOptimG(int n, double *x, double *g, void *ex) {
  int *id = (int*)ex;
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve == 1) return;
  lpInner(x, g, *id);
  fInd->nInnerG++;
}

// Scli-lab style cost function for inner
void innerCost(int *ind, int *n, double *x, double *f, double *g, int *ti, float *tr, double *td, int *id){
  rx = getRxSolve_();
  focei_ind *fInd = &(inds_focei[*id]);
  if (fInd->badSolve==1) {
    return;
  }
  if (*ind==2 || *ind==4) {
    // Function
    // Make sure ID remains installed
    *f = likInner0(x, *id);
    if (ISNA(*f))  {
      fInd->badSolve=1;
    }
    fInd->nInnerF++;
    // if (op_focei.printInner != 0 && fInd->nInnerF % op_focei.printInner == 0){
    // RSprintf("%03d: ", *id);
    // for (int i = 0; i < *n; i++) RSprintf(" %#10g", x[i]);
    // RSprintf(" (nG: %d)\n", fInd->nInnerG);
    // }
  }
  if (*ind==3 || *ind==4) {
    // Gradient

    lpInner(x, g, *id);
    fInd->nInnerG++;
  }
}

static inline int innerEval(int id){
  focei_ind *fInd = &(inds_focei[id]);
  // Use eta
  double lik0 = likInner0(fInd->eta, id);
  if (ISNA(lik0)) return 0;
  lik0 = LikInner2(fInd->eta, 0, id);
  if (ISNA(lik0)) return 0;
  return 1;
}

static inline int innerOpt1(int id, int likId) {
  focei_ind *fInd = &(inds_focei[id]);
  focei_options *fop = &op_focei;
  if (op_focei.neta == 0) {
    double lik = likInner0(NULL, id);
    if (ISNA(lik)) return 0;
    lik = LikInner2(NULL, 0, id);
    if (ISNA(lik)) return 0;
    return 1;
  }
  fInd->nInnerF=0;
  fInd->nInnerG=0;
  // Zero when eta has changed from the last computed value.
  if (fInd->setup) {
    bool _etaChanged = false;
    for (int _j = 0; _j < op_focei.neta; _j++) {
      if (fInd->eta[_j] != fInd->oldEta[_j]) { _etaChanged = true; break; }
    }
    if (_etaChanged) {
      rx = getRxSolve_();
      rx_solving_options_ind *_zInd = getSolvingOptionsInd(rx, getRxId(id));
      int _nobs = getIndNallTimes(_zInd) - getIndNdoses(_zInd) - getIndNevid2(_zInd);
      std::fill_n(fInd->lp, op_focei.neta, 0.0);
      std::fill_n(fInd->a,   (size_t)_nobs * op_focei.neta, 0.0);
      std::fill_n(fInd->B,   _nobs, 0.0);
      std::fill_n(fInd->c,   (size_t)_nobs * op_focei.neta, 0.0);
      std::fill_n(fInd->Vid, (size_t)_nobs * _nobs, 0.0);
      fInd->setup = 0;
    }
  }
  bool n1qn1Inner = true;
  // Use eta
  // Convert Zm to Hessian, if applicable.
  mat etaMat(fop->neta, 1, fill::zeros);
  if (op_focei.mceta == -2 || op_focei.mceta == -1) {
    // Almquist Eq-48 warm-start: extrapolate the next starting eta from the last
    // analytic gradient's EBE sensitivity, eta^0 = eta*_s + (d eta*/d theta)
    // (theta_now - theta_grad) (scaled space).  Only when the analytic gradient
    // supplied etaP (fast=TRUE) and this is a real inner solve (not a gradient
    // perturbation or the cov/linearization step); otherwise fall through keeping
    // the last eta (the legacy mceta=-1 behavior).
    if (!op_focei.calcGrad && op_focei.etaPValid && op_focei.getaP != NULL &&
        op_focei.maxInnerIterations > 0) {
      double *etaP_i = &op_focei.getaP[(size_t)id * op_focei.neta * op_focei.npars];
      std::vector<double> etaNew(op_focei.neta);
      for (int j = 0; j < op_focei.neta; j++) etaNew[j] = fInd->eta[j];
      for (int k = 0; k < (int)op_focei.npars; k++) {
        double dth = op_focei.theta[k] - op_focei.etaPTheta[k];
        if (dth == 0.0) continue;
        for (int j = 0; j < op_focei.neta; j++)
          etaNew[j] += etaP_i[(size_t)j + (size_t)op_focei.neta * k] * dth;
      }
      if (op_focei.mceta == -2) {
        // etaNew in bound -> start there (warmZm re-seeds the Hessian at etaNew);
        // else keep eta*_s -- the standardized-eta reset below zeros it when eta*_s
        // is itself out of bound (the "both out of bound -> reset to zero" case).
        if (etaInBound(&etaNew[0])) {
          std::copy(etaNew.begin(), etaNew.end(), &fInd->eta[0]);
        }
      } else {
        // mceta == -1: jump from etaNew to zero -- keep whichever gives the better
        // (lower -logLik) inner objective.
        double fNew = likInner0(&etaNew[0], id);
        std::fill(&fInd->tryEta[0], &fInd->tryEta[0] + op_focei.neta, 0.0);
        double fZero = likInner0(fInd->tryEta, id);
        if (R_FINITE(fNew) && (!R_FINITE(fZero) || fNew <= fZero)) {
          std::copy(etaNew.begin(), etaNew.end(), &fInd->eta[0]);
        } else {
          std::fill(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, 0.0);
        }
      }
    }
  } else if (op_focei.mceta == 0) {
    // always reset to zero
    std::fill(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, 0.0);
  } else if (op_focei.mceta >= 1 &&
             static_cast<arma::uword>(id) < op_focei.mcetaSamples.n_slices) {
    // mceta sampling: ETA samples pre-drawn serially in innerOpt() (mcetaSamples
    // cube); guard skips subjects with no slice (maxInnerIterations == 0, e.g.
    // covariance/linearization step) to avoid an out-of-bounds Cube::slice().
    int nmc = op_focei.mceta-1;
    double fcur = likInner0(fInd->eta, id); // last eta
    std::fill(&fInd->tryEta[0], &fInd->tryEta[0] + op_focei.neta, 0.0);
    double ftry = likInner0(fInd->tryEta, id); // zero eta
    int sampCol = 0;
    while (true) {
      if (ftry < fcur) {
        std::copy(&fInd->tryEta[0], &fInd->tryEta[0] + op_focei.neta, &fInd->eta[0]);
        fcur = ftry;
      }
      if (nmc <= 0) break;
      nmc--;
      // Read the next pre-drawn sample for this individual.
      arma::vec samp = op_focei.mcetaSamples.slice(id).col(sampCol);
      sampCol++;
      std::copy(samp.begin(), samp.end(), &fInd->tryEta[0]);
      ftry = likInner0(fInd->tryEta, id); // sampled eta
    }
  }
  if (!op_focei.calcGrad) {
    if (op_focei.resetEtaSize <= 0) {
      if (op_focei.resetHessianAndEta){
        fInd->mode = 1;
        fInd->uzm = 1;
        if (n1qn1Inner) op_focei.didHessianReset.store(1, std::memory_order_relaxed);
      }
      resetEtaSelective(fInd, op_focei.neta);
      op_focei.didEtaReset.store(1, std::memory_order_relaxed);
    } else if (R_FINITE(op_focei.resetEtaSize)) {
      std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
      // Standardized ETAs
      // chol(omega^-1) %*% eta
      mat etaRes = op_focei.cholOmegaInv * etaMat;
      bool doBreak = false;
      for (unsigned int j = etaRes.n_rows; j--;){
        if (isMuRefCovProtected(j)) continue; // mu-ref-covariate etas never trigger a reset
        if (std::fabs(etaRes(j, 0)) >= op_focei.resetEtaSize){
          if (op_focei.resetHessianAndEta){
            fInd->mode = 1;
            fInd->uzm = 1;
            if (n1qn1Inner) op_focei.didHessianReset.store(1, std::memory_order_relaxed);
          }
          resetEtaSelective(fInd, op_focei.neta);
          op_focei.didEtaReset.store(1, std::memory_order_relaxed);
          doBreak=true;
          break;
        }
      }
      if (!doBreak){
        etaRes = op_focei.eta1SD % etaMat;
        for (unsigned int j = etaRes.n_rows; j--;){
          if (isMuRefCovProtected(j)) continue; // mu-ref-covariate etas never trigger a reset
          if (std::fabs(etaRes(j, 0)) >= op_focei.resetEtaSize){
            if (op_focei.resetHessianAndEta){
              fInd->mode = 1;
              fInd->uzm = 1;
              if (n1qn1Inner) op_focei.didHessianReset.store(1, std::memory_order_relaxed);
            }
            resetEtaSelective(fInd, op_focei.neta);
            op_focei.didEtaReset.store(1, std::memory_order_relaxed);
            break;
          }
        }
      }
    }
  }
  if (n1qn1Inner) {
    if (op_focei.warm == 1) warmZm(fInd, id);
    else updateZm(fInd);
    std::fill_n(&fInd->var[0], fop->neta, 0.1);
  }
  int npar = fop->neta;
  std::copy(&fInd->eta[0], &fInd->eta[0]+fop->neta, fInd->x);
  double f = 0.0, epsilon = max2(fop->epsilon, sqrt(DBL_EPSILON));

  // Since these are pointers, without reassignment they are modified.
  int mode = fInd->mode, maxInnerIterations=fop->maxInnerIterations,
    nsim=fop->nsim, imp=fop->imp;
  int izs = 0; float rzs = 0.0f; double dzs = 0.0;

  int nF = fInd->nInnerF;
  if (n1qn1Inner) {
    fInd->badSolve = 0;
    n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
           fInd->var, &epsilon,
           &mode, &maxInnerIterations, &nsim,
           &imp, fInd->zm, &izs, &rzs, &dzs, &id);
    if (ISNA(f)) return 0;
    nF = fInd->nInnerF-nF;
    // REprintf("innerCost id: %d, fInd->nInnerF: %d", id, fInd->nInnerF);
    // If stays at zero try another point?
    if (fInd->doEtaNudge == 1 && op_focei.etaNudge != 0.0){
      bool tryAgain=false;
      // if (nF <= 3) tryAgain = true;
      if (!tryAgain){
        tryAgain = true;
        for (int i = fop->neta; i--;){
          if (fInd->x[i] != 0.0){
            tryAgain=false;
            break;
          }
        }
      }
      if (tryAgain) {
        // an ETA actually stayed at zero, so a real nudge is performed; only
        // then flag it (mirrors didHessianReset below) so the "initial ETAs
        // were nudged" warning is not raised when the inner optimization moved
        // the ETAs off zero on its own and no nudge was needed
        op_focei.didEtaNudge.store(1, std::memory_order_relaxed);
        fInd->mode = 1;
        fInd->uzm = 1;
        op_focei.didHessianReset.store(1, std::memory_order_relaxed);
        if (op_focei.warm == 1) mode = 1;
        std::fill_n(fInd->x, fop->neta, op_focei.etaNudge);
        //nF = fInd->nInnerF;
        fInd->badSolve = 0;
        n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
               fInd->var, &epsilon,
               &mode, &maxInnerIterations, &nsim,
               &imp, fInd->zm,
               &izs, &rzs, &dzs, &id);
        if (ISNA(f)) return 0;
        // nF = fInd->nInnerF - nF;
        // if (nF > 3) tryAgain = false;
        if (!tryAgain) {
          tryAgain = true;
          for (int i = fop->neta; i--;){
            if (fInd->x[i] != op_focei.etaNudge){
              tryAgain=false;
              break;
            }
          }
        }
        if (tryAgain) {
          fInd->mode = 1;
          fInd->uzm = 1;
          op_focei.didHessianReset.store(1, std::memory_order_relaxed);
          if (op_focei.warm == 1) mode = 1;
          std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge);
          nF = fInd->nInnerF;
          fInd->badSolve = 0;
          n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
                 fInd->var, &epsilon,
                 &mode, &maxInnerIterations, &nsim,
                 &imp, fInd->zm, &izs, &rzs, &dzs, &id);
          if (ISNA(f)) return 0;
          // nF = fInd->nInnerF - nF;
          // if (nF > 3) tryAgain = false;
          if (!tryAgain){
            tryAgain = true;
            for (int i = fop->neta; i--;){
              if (fInd->x[i] != -op_focei.etaNudge){
                tryAgain=false;
                break;
              }
            }
          }
          if (tryAgain){
            fInd->mode = 1;
            fInd->uzm = 1;
            op_focei.didHessianReset.store(1, std::memory_order_relaxed);
            if (op_focei.warm == 1) mode = 1;
            std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge2);
            nF = fInd->nInnerF;
            fInd->badSolve = 0;
            n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
                   fInd->var, &epsilon,
                   &mode, &maxInnerIterations, &nsim,
                   &imp, fInd->zm, &izs, &rzs, &dzs, &id);
            if (ISNA(f)) return 0;
            // nF = fInd->nInnerF - nF;
            // if (nF > 3) tryAgain = false;
            if (!tryAgain){
              tryAgain = true;
              for (int i = fop->neta; i--;){
                if (fInd->x[i] != -op_focei.etaNudge2){
                  tryAgain=false;
                  break;
                }
              }
            }
            if (tryAgain) {
              fInd->mode = 1;
              fInd->uzm = 1;
              op_focei.didHessianReset.store(1, std::memory_order_relaxed);
              if (op_focei.warm == 1) mode = 1;
              std::fill_n(fInd->x, fop->neta, +op_focei.etaNudge2);
              nF = fInd->nInnerF;
              fInd->badSolve = 0;
              n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
                     fInd->var, &epsilon,
                     &mode, &maxInnerIterations, &nsim,
                     &imp, fInd->zm, &izs, &rzs, &dzs, &id);
              if (ISNA(f)) return 0;
              // nF = fInd->nInnerF - nF;
              // if (nF > 3) tryAgain = false;
              if (!tryAgain){
                tryAgain = true;
                for (int i = fop->neta; i--;){
                  if (fInd->x[i] != +op_focei.etaNudge2){
                    tryAgain=false;
                    break;
                  }
                }
              }
              if (tryAgain) {
                std::fill_n(fInd->x, fop->neta, 0);
                std::fill_n(&fInd->var[0], fop->neta, 0.2);
                nF = fInd->nInnerF;
                fInd->badSolve = 0;
                n1qn1_(innerCost, &npar, fInd->x, &f, fInd->g,
                       fInd->var, &epsilon,
                       &mode, &maxInnerIterations, &nsim,
                       &imp, fInd->zm,
                       &izs, &rzs, &dzs, &id);
                if (ISNA(f)) return 0;
                //nF = fInd->nInnerF-nF;
                // if (nF > 3) tryAgain = false;
                if (!tryAgain){
                  tryAgain = true;
                  for (int i = fop->neta; i--;){
                    if (fInd->x[i] != 0.0){
                      tryAgain=false;
                      break;
                    }
                  }
                }
                if (tryAgain) {
                  std::fill_n(fInd->x, fop->neta, 0);
                }
                std::fill_n(&fInd->var[0], fop->neta, 0.1);
              }
            }
          }
        }
      }
    }
  } else {
    int fail=0, fncount=0, grcount=0;
    char msg[100];
    fInd->badSolve = 0;
    lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
             op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
             &fail, (void*)(&id), op_focei.factr,
             op_focei.pgtol, &fncount, &grcount,
             op_focei.maxInnerIterations, msg, 0, -1,
             op_focei.abstol, op_focei.reltol, fInd->g);
    if (ISNA(f)) return 0;
    // if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    //   // did not converge
    //   if (fInd->doEtaNudge == 1 && op_focei.etaNudge != 0.0){
    //  std::fill_n(fInd->x, fop->neta, op_focei.etaNudge);
    //  fail=0;
    //  lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
    //       op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
    //       &fail, (void*)(&id), op_focei.factr,
    //       op_focei.pgtol, &fncount, &grcount,
    //       op_focei.maxInnerIterations, msg, 0, -1,
    //       op_focei.abstol, op_focei.reltol, fInd->g);
    //  if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    //    std::fill_n(fInd->x, fop->neta, -op_focei.etaNudge);
    //    lbfgsb3C(npar, op_focei.lmm, fInd->x, op_focei.etaLower,
    //       op_focei.etaUpper, op_focei.nbdInner, &f, innerOptimF, innerOptimG,
    //       &fail, (void*)(&id), op_focei.factr,
    //       op_focei.pgtol, &fncount, &grcount,
    //       op_focei.maxInnerIterations, msg, 0, -1,
    //       op_focei.abstol, op_focei.reltol, fInd->g);
    //    if (fail != 6 && fail != 7 && fail != 8 && fail != 27){
    //      std::fill_n(fInd->x, fop->neta, 0);
    //    }
    //  }
    //   }
    // }
  }
  // only nudge once
  fInd->doEtaNudge=0;

  std::copy(&fInd->x[0],&fInd->x[0]+fop->neta,&fInd->eta[0]);
  if (_innerParallel.load(std::memory_order_acquire) == 0) {
    // Serial path: update running mean/variance of ETAs (Welford's algorithm)
    // Update variances
    std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
    op_focei.n = op_focei.n + 1.0;
    mat oldM = op_focei.etaM;
    op_focei.etaM = op_focei.etaM + (etaMat - op_focei.etaM)/op_focei.n;
    op_focei.etaS = op_focei.etaS + (etaMat - op_focei.etaM) %  (etaMat - oldM);
  }
  // In parallel mode: etaM/etaS accumulated after the parallel region
  fInd->llik = f;
  // Use saved Hessian on next opimization.
  fInd->mode=2;
  fInd->uzm =0;
  if (ISNA(LikInner2(fInd->eta, likId, id))) return 0;
  return 1;
}

void parHistData(Environment e, bool focei);

void foceiPrintInfo() {
  arma::Row<int> etaTrans(op_focei.etaTrans, op_focei.neta);
  arma::Row<int> nbdInner(op_focei.nbdInner, op_focei.neta);
  arma::Row<int> xPar(op_focei.xPar, op_focei.ntheta + op_focei.omegan);
  arma::Row<int> thetaTrans(op_focei.thetaTrans, op_focei.ntheta + op_focei.omegan);
  arma::Row<int> fixedTrans(op_focei.fixedTrans, op_focei.ntheta + op_focei.omegan);
  arma::Row<int> etaFD(op_focei.etaFD, op_focei.neta);

  arma::rowvec fullTheta(op_focei.fullTheta, op_focei.ntheta+op_focei.omegan);
  arma::rowvec theta(op_focei.theta, op_focei.ntheta+op_focei.omegan);
  arma::rowvec initPar(op_focei.initPar, op_focei.ntheta+op_focei.omegan);
  arma::rowvec scaleC(op_focei.scaleC, op_focei.ntheta+op_focei.omegan);

  REprintf("etaTrans\n");
  print(wrap(etaTrans));

  REprintf("nbdInner\n");
  print(wrap(nbdInner));

  REprintf("xPar\n");
  print(wrap(xPar));

  REprintf("thetaTrans\n");
  print(wrap(thetaTrans));

  REprintf("fixedTrans\n");
  print(wrap(fixedTrans));

  REprintf("etaFD\n");
  print(wrap(etaFD));

  REprintf("fullTheta\n");
  print(wrap(fullTheta));

  REprintf("theta\n");
  print(wrap(theta));

  REprintf("initPar\n");
  print(wrap(initPar));

  REprintf("scaleC\n");
  print(wrap(scaleC));
}


static inline void thetaReset00(NumericVector &thetaIni, NumericVector &omegaTheta, arma::mat &etaMat) {
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Environment thetaReset = nlmixr2[".thetaReset"];
  focei_options *fop = &op_focei;
  thetaReset["maxInnerIterations"]=fop->maxInnerIterations;
  thetaReset["etaMat"] = wrap(etaMat);
  thetaReset["thetaIni"]= thetaIni;
  thetaReset["omegaTheta"] = omegaTheta;
  thetaReset["nF"] = op_focei.nF+op_focei.nF2;
  // Save gill info to skip recalc.
  IntegerVector gillRetC(op_focei.npars);
  std::copy(&op_focei.gillRetC[0], &op_focei.gillRetC[0]+op_focei.npars, gillRetC.begin());
  thetaReset["gillRetC"] = gillRetC;
  IntegerVector gillRet(op_focei.npars);
  std::copy(&op_focei.gillRet[0], &op_focei.gillRet[0]+op_focei.npars, gillRet.begin());
  thetaReset["gillRet"] = gillRet;
  NumericVector gillDf(op_focei.npars);
  std::copy(&op_focei.gillDf[0], &op_focei.gillDf[0]+op_focei.npars, gillDf.begin());
  thetaReset["gillDf"] = gillDf;
  NumericVector gillDf2(op_focei.npars);
  std::copy(&op_focei.gillDf2[0], &op_focei.gillDf2[0]+op_focei.npars, gillDf2.begin());
  thetaReset["gillDf2"] = gillDf2;
  NumericVector gillErr(op_focei.npars);
  std::copy(&op_focei.gillErr[0], &op_focei.gillErr[0]+op_focei.npars, gillErr.begin());
  thetaReset["gillErr"] = gillErr;
  NumericVector rEps(op_focei.npars);
  std::copy(&op_focei.rEps[0], &op_focei.rEps[0]+op_focei.npars, rEps.begin());
  thetaReset["rEps"] = rEps;
  NumericVector aEps(op_focei.npars);
  std::copy(&op_focei.aEps[0], &op_focei.aEps[0]+op_focei.npars, aEps.begin());
  thetaReset["aEps"] = aEps;
  NumericVector rEpsC(op_focei.npars);
  std::copy(&op_focei.rEpsC[0], &op_focei.rEpsC[0]+op_focei.npars, rEpsC.begin());
  thetaReset["rEpsC"] = rEpsC;
  NumericVector aEpsC(op_focei.npars);
  std::copy(&op_focei.aEpsC[0], &op_focei.aEpsC[0]+op_focei.npars, aEpsC.begin());
  thetaReset["aEpsC"] = aEpsC;
  thetaReset["c1"] = op_focei.c1;
  thetaReset["c2"] = op_focei.c2;
  //foceiPrintInfo();
  parHistData(thetaReset, true);
  saveIntoEnvrionment(thetaReset);
}

static inline bool isFixedTheta(int m) {
  unsigned int j;
  for (unsigned int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    if (m == (int)j) return false; // here the parameter is estimated
  }
  return true; // here the parameter is fixed
}

static inline bool thetaReset0(bool forceReset = false) {
  if (op_focei.isSaem) return false;
  NumericVector thetaIni(op_focei.ntheta);
  NumericVector thetaUp(op_focei.ntheta);
  NumericVector thetaDown(op_focei.ntheta);
  LogicalVector adjustEta(op_focei.muRefN);
  bool doAdjust = false;
  for (int ii = (int)op_focei.ntheta; ii--;) {
    thetaIni[ii] = unscalePar(op_focei.fullTheta, ii);
    if (R_FINITE(op_focei.lower[ii])) {
      thetaDown[ii] = unscalePar(op_focei.lower, ii);
    } else {
      thetaDown[ii] = R_NegInf;
    }
    if (R_FINITE(op_focei.upper[ii])) {
      thetaUp[ii]= unscalePar(op_focei.upper, ii);
    } else {
      thetaUp[ii] = std::numeric_limits<double>::infinity();
    }
  }
  double ref=0;
  int ij = 0;
  for (unsigned int ii = op_focei.muRefN; ii--;){
    if (op_focei.muRef[ii] != -1 && op_focei.muRef[ii] < (int)op_focei.ntheta) {
      ij = op_focei.muRef[ii];
      if (isFixedTheta(ij) || isMuRefCovProtected(ii)) {
        // mu-ref-covariate thetas are only ever updated by the
        // mufocei/irlsfocei-family restart-loop's linear-model step, never
        // by this soft mu-shift.
        adjustEta[ii] = false;
      }  else {
        ref = thetaIni[ij] + op_focei.etaM(ii,0);
        if (thetaDown[ij] < ref && thetaUp[ij] > ref) {
          thetaIni[ij] = ref;
          adjustEta[ii] = true;
          doAdjust = true;
        } else {
          adjustEta[ii] = false;
        }
      }
    } else {
      adjustEta[ii] = false;
    }
  }
  if (!doAdjust && !forceReset) {
    return false;
  }

  arma::mat etaMat(getRxNsubAndMix(rx), op_focei.neta);

  for (int ii = getRxNsubAndMix(rx); ii--;) {
    focei_ind *fInd = &(inds_focei[ii]);
    for (int jj = op_focei.neta; jj--; ) {
      if (op_focei.muRef[jj] != -1  && op_focei.muRef[jj] < (int)op_focei.ntheta &&
          adjustEta[jj]) {
        etaMat(ii, jj) = fInd->eta[jj]-op_focei.etaM(jj,0);
      } else {
        etaMat(ii, jj) = fInd->eta[jj];
      }
    }
  }
  // Update omega estimates
  NumericVector omegaTheta(op_focei.omegan);

  std::copy(&op_focei.fullTheta[0] + op_focei.ntheta,
            &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
            omegaTheta.begin());
  thetaReset00(thetaIni, omegaTheta, etaMat);
  return true;
}

void thetaReset(double size){
  if (op_focei.isSaem) return;
  if (op_focei.maxOuterIterations <= 0) return;
  if (std::isinf(size)) return;
  mat etaRes =  op_focei.eta1SD % op_focei.etaM; //op_focei.cholOmegaInv * etaMat;
  double res=0;
  for (unsigned int j = etaRes.n_rows; j--;) {
    if (isMuRefCovProtected(j)) continue; // mu-ref-covariate etas never trigger a theta reset
    res = etaRes(j, 0);
    res = res < 0 ? -res : res;
    if (res >= size) { // Says reset;
      if (thetaReset0()) {
        if (op_focei.didEtaReset==1) {
          warning(_("mu-referenced Thetas were reset during optimization; (Can control by foceiControl(resetThetaP=.,resetThetaCheckPer=.,resetThetaFinalP=.))"));
        }
        stop("theta reset");
      }
    }
  }
}

void thetaResetZero() {
  if (op_focei.isSaem) return;
  if (op_focei.maxOuterIterations <= 0) return;
  thetaReset0(true);
  warning(_("thetas were reset during optimization because of a zero gradient"));
  stop("theta reset0");
}

void thetaResetObj(Environment e) {
  if (op_focei.isSaem) return;
  // Check to see if the last objective function is the lowest, otherwise reset to the lowest thetaResetObj
  if (op_focei.maxOuterIterations > 0) {
    parHistData(e, true);
    if (e.exists("parHistData")) {
      if (TYPEOF(e["parHistData"]) == VECSXP) {
        List parHistData = e["parHistData"];
        IntegerVector iter = parHistData["iter"];
        IntegerVector type = parHistData["type"];
        NumericVector obj = parHistData["objf"];
        int maxiter=-1, minObjId=-1;
        double minObj = std::numeric_limits<double>::infinity();
        for (int i = obj.size(); i--;) {
          if (type[i] == 5) {
            if (!ISNA(obj[i])) {
              if (obj[i] < minObj) {
                minObj = obj[i];
                minObjId = i;
              }
            }
            maxiter = max2(maxiter, iter[i]);
          }
        }
        if (iter[minObjId] != maxiter) {
          warning(_("last objective function was not at minimum, possible problems in optimization"));
          // stop("theta resetZ");
        }
      }
    }
  }
}

static inline int didInnerResetPointFail(focei_ind *indF, int& id, double point) {
  indF->mode = 1;
  indF->uzm = 1;
  op_focei.didHessianReset=1;
  if (point == 0.0) {
    std::fill(&indF->eta[0], &indF->eta[0] + op_focei.neta, point);
  } else {
    // chol(omega^-1) %*% eta = point for each eta
    //
    //mat etaRes = op_focei.cholOmegaInv * etaMat
    for (int j = op_focei.neta; j--; ) {
      indF->eta[j] =  point/op_focei.cholOmegaInv(j, j);
    }
  }
  return !innerOpt1(id, 0);
}

static inline int didInnerResetFailTryOne(focei_ind *indF, int& id) {
  if (didInnerResetPointFail(indF, id, 0.0)) {
    if (op_focei.etaNudge != 0.0) {
      if (didInnerResetPointFail(indF, id, op_focei.etaNudge) &&
          didInnerResetPointFail(indF, id, -op_focei.etaNudge)) {
        if (op_focei.etaNudge2 != 0.0) {
          if (didInnerResetPointFail(indF, id, op_focei.etaNudge2) &&
              didInnerResetPointFail(indF, id, -op_focei.etaNudge2)) {
            return 1; // failed
          } else {
            return 0; // reset on one of the eta nudge
          }
        } else {
          return 1; // failed reset
        }
      } else {
        return 0; // reset
      }
    } else {
      return 1; // failed reset
    }
  }
  return 0;
}

static inline int didInnerResetFail(focei_ind *indF, int& id) {
  if (didInnerResetFailTryOne(indF, id)) {
    if (op_focei.canDoFD && op_focei.fallbackFD) {
      indF->doFD = 1;
      if (didInnerResetFailTryOne(indF, id)) {
        indF->doFD = 0;
        return 1; // failed reset
      } else {
        indF->doFD = 0;
        return 0;
      }
    } else {
      return 1;
    }
  }
  return 0;
}

static inline void resetToZeroWithoutOpt(focei_ind *indF, int& id) {
  // Just use ETA=0
  std::fill(&indF->eta[0], &indF->eta[0] + op_focei.neta, 0.0);
  if (!innerEval(id)) {
    indF->lik[0] = NA_REAL;
    if (_innerParallel.load(std::memory_order_acquire) == 0) {
      warning(_("bad solve during optimization"));
    }
  }
}

static inline void innerOptId(int id) {
  focei_ind *indF = &(inds_focei[id]);
  indF->parWarnBadHess = 0;
  indF->parErrorNoEta = 0;
  if (!innerOpt1(id, 0)) {
    // First try resetting ETA
    if (didInnerResetFail(indF, id)) {
      if(!op_focei.noabort){
        if (_innerParallel.load(std::memory_order_acquire) == 0) {
          stop("Could not find the best eta even hessian reset and eta reset for ID %d.", id+1);
        } else {
          indF->parErrorNoEta = 1; // deferred to post-parallel
        }
      } else if (indF->doChol == 1){
        indF->doChol = 0; // Use generalized cholesky decomposition
        if (didInnerResetFail(indF, id)) {
          resetToZeroWithoutOpt(indF, id);
        }
        indF->doChol = 1; // Use cholesky again.
      } else {
        resetToZeroWithoutOpt(indF, id);
      }
    }
  }
}



// mufocei/irlsfocei regression update, called once per outer iteration after
// every subject's eta is finalized. Per group: build phi_i = fullTheta[popIdx]
// + covariate_effect_i + eta_i, regress on the free covariates (OLS "lin" or
// curvature-weighted "irls"), write new thetas into fullTheta and residuals
// into each subject's eta via setIndParPtr() (like updateTheta(), no R round-trip).
// Returns the max absolute theta change, used by innerOpt() to decide whether
// another {re-optimize etas, regress} cycle is needed (muGroupMaxCycles/Tol).
static inline double updateMuGroups() {
  if (op_focei.muModel == 0 || op_focei.muGroupN == 0) return 0.0;
  rx = getRxSolve_();
  int nsubAll = (int)getRxNsubAndMix(rx);
  double maxDelta = 0.0;

  for (unsigned int g = 0; g < op_focei.muGroupN; g++) {
    int popIdx = op_focei.muGroupTheta[g];
    int etaIdx = op_focei.muGroupEta[g];
    unsigned int covStart = (unsigned int)op_focei.muGroupCovStart[g];
    unsigned int covCount = (unsigned int)op_focei.muGroupCovCount[g];

    // Split covariates into known-offset (user-fixed or bounded, not estimated
    // here) vs free (estimated by the regression). Both read fullTheta[covTheta]
    // fresh each call, so bounded ones automatically pick up the outer
    // optimizer's latest value with no extra logic.
    std::vector<unsigned int> freeCovCols;
    std::vector<int> freeCovTheta;
    arma::vec offset(nsubAll, arma::fill::zeros);
    for (unsigned int c = 0; c < covCount; c++) {
      unsigned int flatIdx = covStart + c;
      int covTheta = op_focei.muGroupCovTheta[flatIdx];
      bool userFixed = op_focei.muGroupCovUserFixed != NULL &&
        op_focei.muGroupCovUserFixed[flatIdx] != 0;
      bool bounded = op_focei.muGroupCovBounded != NULL &&
        op_focei.muGroupCovBounded[flatIdx] != 0;
      if (userFixed || bounded) {
        double coefVal = op_focei.fullTheta[covTheta];
        for (int id = 0; id < nsubAll; id++) {
          int rid = getRxId(id);
          offset(id) += coefVal * op_focei.muGroupCovData(rid, flatIdx);
        }
      } else {
        freeCovCols.push_back(flatIdx);
        freeCovTheta.push_back(covTheta);
      }
    }

    int nFree = (int)freeCovCols.size();
    arma::vec y(nsubAll);
    arma::mat X(nsubAll, nFree + 1, arma::fill::ones); // column 0 = intercept
    arma::vec w(nsubAll, arma::fill::ones);
    for (int id = 0; id < nsubAll; id++) {
      focei_ind *fInd = &(inds_focei[id]);
      int rid = getRxId(id);
      double phi = op_focei.fullTheta[popIdx] + offset(id) + fInd->eta[etaIdx];
      for (int c = 0; c < nFree; c++) {
        double covVal = op_focei.muGroupCovData(rid, freeCovCols[c]);
        phi += covVal * op_focei.fullTheta[freeCovTheta[c]];
        X(id, c + 1) = covVal;
      }
      y(id) = phi;
      if (op_focei.muModel == 2) {
        // "irls": weight by each subject's n1qn1 per-eta precision estimate
        // (fInd->var) -- a cheap default, not independently validated.
        double v = fInd->var[etaIdx];
        w(id) = (v > 0 && R_FINITE(v)) ? (1.0 / v) : 1.0;
      }
    }

    arma::vec beta;
    bool solveOk = true;
    try {
      if (op_focei.muModel == 2) {
        arma::mat Wd = arma::diagmat(w);
        beta = arma::solve(X.t() * Wd * X, X.t() * Wd * y);
      } else {
        beta = arma::solve(X, y);
      }
    } catch (...) {
      solveOk = false;
    }
    if (!solveOk || beta.n_elem != (unsigned int)(nFree + 1) ||
        !beta.is_finite()) {
      // Singular/ill-conditioned design matrix: leave this group's thetas/etas
      // untouched this iteration rather than propagating garbage.
      continue;
    }

    maxDelta = std::max(maxDelta, std::fabs(beta(0) - op_focei.fullTheta[popIdx]));
    op_focei.fullTheta[popIdx] = beta(0);
    for (int c = 0; c < nFree; c++) {
      maxDelta = std::max(maxDelta, std::fabs(beta(c + 1) - op_focei.fullTheta[freeCovTheta[c]]));
      op_focei.fullTheta[freeCovTheta[c]] = beta(c + 1);
    }

    arma::vec fitted = X * beta;
    for (int id = 0; id < nsubAll; id++) {
      focei_ind *fInd = &(inds_focei[id]);
      fInd->eta[etaIdx] = y(id) - fitted(id);
    }

    // Propagate the updated theta(s)/eta to every subject's solve state,
    // exactly like updateTheta()'s setIndParPtr() loop.
    for (int id = 0; id < nsubAll; id++) {
      rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
      focei_ind *fInd = &(inds_focei[id]);
      setIndParPtr(ind, op_focei.thetaTrans[popIdx], op_focei.fullTheta[popIdx]);
      for (int c = 0; c < nFree; c++) {
        setIndParPtr(ind, op_focei.thetaTrans[freeCovTheta[c]],
                     op_focei.fullTheta[freeCovTheta[c]]);
      }
      setIndParPtr(ind, op_focei.etaTrans[etaIdx], fInd->eta[etaIdx]);
    }
  }
  return maxDelta;
}

void innerOpt() {
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  int cores = getOpCores(op);
  if (op_focei.neta > 0 && !op_focei.covFdDirect) {   // covFdDirect: keep the FD-perturbed Omega
    op_focei.omegaInv=getOmegaInv();
    op_focei.logDetOmegaInv5 = getOmegaDet();
  }
  // Pre-draw per-subject ETA samples serially before the parallel for-loop so
  // workers only do memory access (no R API calls), making mceta safe under cores > 1.
  if (op_focei.mceta >= 1 && op_focei.maxInnerIterations > 0) {
    int nsubAll = (int)getRxNsubAndMix(rx);
    int nmc = op_focei.mceta - 1;
    if (nmc > 0 && op_focei.neta > 0) {
      op_focei.mcetaSamples.set_size(op_focei.neta, nmc, nsubAll);
      NumericMatrix omega = getOmega();
      Function loadNamespace("loadNamespace", R_BaseNamespace);
      Environment nlmixr2 = loadNamespace("nlmixr2est");
      Function fSample = as<Function>(nlmixr2[".sampleOmega"]);
      for (int id = 0; id < nsubAll; ++id) {
        for (int k = 0; k < nmc; ++k) {
          NumericMatrix samp = fSample(omega);
          std::copy(samp.begin(), samp.end(),
                    op_focei.mcetaSamples.slice(id).colptr(k));
        }
      }
    } else {
      op_focei.mcetaSamples.reset();
    }
  } else {
    op_focei.mcetaSamples.reset();
  }
  if (op_focei.maxInnerIterations <= 0){
    std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
    for (int id = 0; id < getRxNsubAndMix(rx); id++){
      focei_ind *indF = &(inds_focei[id]);
      indF->doChol = 1;
      if (!innerEval(id)) {
        indF->doChol = 0; // Use generalized cholesky decomposition
        innerEval(id);
        // Not thread safe
        warning(_("non-positive definite individual Hessian at solution(ID=%d); FOCEi objective functions may not be comparable"),id);
        indF->doChol = 1; // Cholesky again.
      }
    }
  } else {
    int nsub_orig = getRxNsub(rx);
    int nMix = op_focei.mixIdxN + 1;
    int nsub = nsub_orig * nMix;
    bool _doParallel = (cores > 1) && solveMethodThreadSafe(op);

    // A single {re-optimize etas, updateMuGroups()} pass isn't always enough:
    // the regression's new thetas can shift each subject's conditional-mode
    // eta, changing what the next regression would produce, while the outer
    // optimizer's gradient check (which ignores mu-group params) could declare
    // convergence early. So cycle until mu-group thetas move less than
    // muGroupTol or muGroupMaxCycles is hit; reduces to one pass when there
    // are no mu-groups or during gradient/Hessian FD perturbations.
    int muCycles = (op_focei.muModel != 0 && op_focei.muGroupN > 0 &&
                     !op_focei.calcGrad) ? op_focei.muGroupMaxCycles : 1;
    for (int muCyc = 0; muCyc < muCycles; muCyc++) {
    if (_doParallel) {
      sortIds(rx, 2); // only initialize if rx->ordId == NULL
      _innerParallel.store(1, std::memory_order_release);
    }
    // Cross-DLL OpenMP thread-id fix: rxode2 and nlmixr2est each statically
    // link their own libgomp on Windows, so rxode2's omp_get_thread_num()
    // inside ind_solve returns 0 for every worker, collapsing its per-thread
    // solve buffers onto slot 0 and corrupting the heap. Hand rxode2 our real
    // thread id around each per-subject solve (resolved once on main thread).
    //
    // Mixture models (nMix > 1): components share the physical subjects' data
    // and solving structures (the memory-saving default from the saem/focei
    // setup), so no two mixture components may ever be solved concurrently.
    // Solve component-by-component (serial jMix OUTSIDE), parallelizing over
    // physical subjects INSIDE each component -- partial parallelism, but
    // mixture-safe. Non-mixture models (nMix == 1) keep the single fully
    // parallel subject loop.
    for (int jMix = 0; jMix < nMix; jMix++) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(_doParallel)
#endif
      for (int i = 0; i < nsub_orig; i++) {
        int _id0 = _doParallel ? (getOrdId(rx, i) - 1) : i;
        int _id = _id0 + jMix * nsub_orig;
#ifdef _OPENMP
        if (_doParallel) {
          setRxThreadId(omp_get_thread_num());
          try {
            innerOptId(_id);
          } catch (...) {
            inds_focei[_id].parErrorNoEta = 1;
          }
          setRxThreadId(-1);
        } else {
#endif
          innerOptId(_id);
#ifdef _OPENMP
        }
#endif
      }
    }
    _innerParallel.store(0, std::memory_order_release);
    if (_doParallel) {
      sortIds(rx, 0);
    }

    if (_doParallel) {
      // --- Post-parallel: deferred error handling ---
      for (int id = 0; id < nsub; id++) {
        focei_ind *indF = &(inds_focei[id]);
        if (indF->parErrorNoEta) {
          stop("Could not find the best eta even hessian reset and eta reset for ID %d.", id+1);
        }
      }
      // --- Post-parallel: recompute etaM/etaS from per-subject etas ---
      if (op_focei.neta > 0) {
        std::fill(op_focei.etaM.begin(), op_focei.etaM.end(), 0.0);
        std::fill(op_focei.etaS.begin(), op_focei.etaS.end(), 0.0);
        op_focei.n = 0.0;
        mat etaMat(op_focei.neta, 1);
        for (int id = 0; id < nsub; id++) {
          focei_ind *fInd = &(inds_focei[id]);
          std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
          op_focei.n = op_focei.n + 1.0;
          mat oldM = op_focei.etaM;
          op_focei.etaM = op_focei.etaM + (etaMat - op_focei.etaM)/op_focei.n;
          op_focei.etaS = op_focei.etaS + (etaMat - op_focei.etaM) % (etaMat - oldM);
        }
      }
    }

    // Gated by !calcGrad: this site is also reached during gradient/Hessian FD
    // perturbations (foceiCalcR()), where overwriting fullTheta would corrupt
    // the derivative being computed.
    if (!op_focei.calcGrad) {
      double muDelta = updateMuGroups();
      if (muDelta <= op_focei.muGroupTol) break;
    } else {
      break;
    }
    }

    // Reset ETA variances for next step
    if (op_focei.neta > 0){
      if (op_focei.zeroGrad) {
        thetaResetZero();
      }
      op_focei.eta1SD = 1/sqrt(op_focei.etaS);
      if (!op_focei.calcGrad && op_focei.maxOuterIterations > 0 &&
          (!op_focei.initObj || op_focei.checkTheta==1) &&
          R_FINITE(op_focei.resetThetaSize)){
        // Not thread safe...
        thetaReset(op_focei.resetThetaSize);
      }
      std::fill(op_focei.etaM.begin(),op_focei.etaM.end(), 0.0);
      std::fill(op_focei.etaS.begin(),op_focei.etaS.end(), 0.0);
      op_focei.n = 0.0;
    }
  }
  Rcpp::checkUserInterrupt();
}

static inline double foceiLik0(double *theta) {
  updateTheta(theta);
  innerOpt();
  double lik = 0.0;
  double cur;

  if (op_focei.mixIdxN == 0) {
    for (int id=getRxNsub(rx); id--;){
      focei_ind *fInd = &(inds_focei[id]);
      cur = fInd->lik[0];
      if (ISNA(cur) || std::isinf(cur) || std::isnan(cur)) {
        cur = -op_focei.badSolveObjfAdj;
      }
      lik += cur;
    }
  } else {
    focei_ind *fInd;
    for (int id=getRxNsub(rx); id--;){
      double tot = 0;
      bool allZero = true;
      double explast = 0.0;
      double mixprob = 0.0;
      for (unsigned int mn = 0; mn < op_focei.mixIdxN + 1; mn++) {
        fInd = &(inds_focei[id + mn*getRxNsub(rx)]);
        cur = fInd->lik[0];
        double ecur = exp(cur);
        if (ISNA(cur) || std::isinf(cur) || std::isnan(cur)) {
          fInd->mixProb[mn] = 0.0;
        } else {
          allZero = false;
          fInd->mixProb[mn] = ecur*op_focei.mixProb[mn];
          tot += fInd->mixProb[mn];
        }
        if (mn == op_focei.mixIdxN) {
          explast = ecur;
        } else {
          fInd->mixProbGrad[mn] = ecur;
        }
      }
      for (unsigned int mn = 0; mn < op_focei.mixIdxN + 1; mn++) {
        fInd->mixProb[mn] = fInd->mixProb[mn]/tot;
        if (mn != op_focei.mixIdxN) {
          // Finish Calculating the gradient based on the probability
          fInd->mixProbGrad[mn] = (fInd->mixProbGrad[mn]-explast)/tot;
        }
        if (mixprob < fInd->mixProb[mn]) {
          mixprob = fInd->mixProb[mn];
          fInd->mixest[0] = mn + 1;
        }
      }
      if (allZero) {
        cur = -op_focei.badSolveObjfAdj;
      } else {
        cur = log(tot);
      }
      lik += cur;
    }
  }
  // Now reset the saved ETAs
  if (op_focei.neta !=0) {
    // All etas = -42;  Unlikely if normal
    std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0);
  }
  return lik;
}


static inline double foceiOfv0(double *theta){
  if (op_focei.objfRecalN != 0 && !op_focei.calcGrad) {
    op_focei.stickyRecalcN1++;
    if (op_focei.indTolRelax) {
      // Per-individual tolFactor is sticky: stiff subjects keep their loosened
      // tolerance across optimizer calls; no tolerance reset.
      if (op_focei.stickyRecalcN1 > op_focei.stickyRecalcN) {
        op_focei.stickyTol=1;
      }
    } else {
      if (op_focei.stickyRecalcN1 <= op_focei.stickyRecalcN) {
        // Reset all subjects' tolerances to prepare for the next evaluation.
        double _resetFactor = pow(op_focei.odeRecalcFactor, -op_focei.objfRecalN);
        int _nsub = (int)getRxNsubAndMix(rx);
        for (int _i = 0; _i < _nsub; _i++) {
          rx_solving_options_ind *_indI = getSolvingOptionsInd(rx, _i);
          setIndTolFactor(_indI, getIndTolFactor(_indI) * _resetFactor);
        }
      } else {
        op_focei.stickyTol=1;
      }
    }
  }
  double ret = -2*foceiLik0(theta);
  rx_solving_options *_op0 = getSolvingOptions(rx);
  while (!op_focei.calcGrad && op_focei.stickyRecalcN1 <= op_focei.stickyRecalcN &&
         (std::isnan(ret) || std::isinf(ret)) &&
         op_focei.objfRecalN < op_focei.maxOdeRecalc){
    op_focei.reducedTol=1;
    if (op_focei.indTolRelax) {
      // Only loosen subjects whose ODE solve produced NaN/Inf; stiff subjects
      // accumulate loosening across calls; non-stiff subjects are unaffected.
      if (getOpNeq(_op0) > 0) {
        int _nsub = (int)getRxNsubAndMix(rx);
        for (int _i = 0; _i < _nsub; _i++) {
          rx_solving_options_ind *_indI = getSolvingOptionsInd(rx, _i);
          double *_solveI = getIndSolve(_indI);
          int _nsolveI = getOpNeq(_op0) * getIndNallTimes(_indI);
          for (int _ns = 0; _ns < _nsolveI; _ns++) {
            if (ISNA(_solveI[_ns]) || std::isnan(_solveI[_ns]) || std::isinf(_solveI[_ns])) {
              setIndTolFactor(_indI, getIndTolFactor(_indI) * op_focei.odeRecalcFactor);
              break;
            }
          }
        }
      }
    } else {
      // Loosen all subjects uniformly.
      int _nsub = (int)getRxNsubAndMix(rx);
      for (int _i = 0; _i < _nsub; _i++) {
        rx_solving_options_ind *_indI = getSolvingOptionsInd(rx, _i);
        setIndTolFactor(_indI, getIndTolFactor(_indI) * op_focei.odeRecalcFactor);
      }
    }
    ret = -2*foceiLik0(theta);
    op_focei.objfRecalN++;
  }
  if (!op_focei.initObj){
    op_focei.initObj=1;
    op_focei.initObjective=std::fabs(ret);
    if (std::isnan(ret)){
      stop(_("NaN while evaluating initial objective function"));
    } else if (std::isinf(ret)) {
      stop(_("infinite while evaluating initial objective function"));
    }
    if (op_focei.scaleObjective == 1) op_focei.scaleObjective=2;
  } else {
    if (std::isnan(ret) || std::isinf(ret)){
      ret=5e100;
    }
  }
  if (op_focei.scaleObjective == 2){
    ret = ret / op_focei.initObjective * op_focei.scaleObjectiveTo;
  }
  if (!op_focei.calcGrad){
    if (op_focei.derivMethodSwitch) {
      double diff = std::fabs(op_focei.lastOfv-ret);
      if (op_focei.derivMethod==0 && diff <= op_focei.derivSwitchTol){
        op_focei.derivMethod=1;
      } else if (op_focei.derivMethod==1 && diff > op_focei.derivSwitchTol){
        op_focei.derivMethod=0;
      }
    }
    if (fabs((op_focei.lastOfv-ret)/max2(op_focei.lastOfv, ret))*100 < op_focei.resetThetaCheckPer){
      op_focei.checkTheta=1;
    } else {
      op_focei.checkTheta=0;
    }
    op_focei.lastOfv = ret;
  }
  return ret;
}

//[[Rcpp::export]]
double foceiLik(NumericVector theta){
  return foceiLik0(&theta[0]);
}

//[[Rcpp::export]]
double foceiOfv(NumericVector theta){
  return foceiOfv0(&theta[0]);
}

void foceiPhiOne(Environment e, List &retC, List &retH, int mixest) {
  if (e.exists("idLvl")) {
    RObject idl = e["idLvl"];
    retH.attr("names") = idl;
    retC.attr("names") = idl;
  }
  bool doDimNames = false;
  List dimn(2);
  if (e.exists("etaNames")) {
    doDimNames=true;
    dimn[0] = e["etaNames"];
    dimn[1] = e["etaNames"];
  }
  for (int j=getRxNsub(rx); j--;) {
    arma::mat H(op_focei.gH +
                (j+mixest*getRxNsub(rx))*op_focei.neta*op_focei.neta,
                op_focei.neta, op_focei.neta, false, true);
    RObject cur = wrap(H);
    if (doDimNames) cur.attr("dimnames") = dimn;
    retH[j] = cur;
    arma::mat cov;
    bool success  = inv(cov, H);
    if (!success) {
      success = pinv(cov, H);
      if (!success) {
        retC[j] = NA_REAL;
      } else {
        cur = wrap(cov);
        if (doDimNames) cur.attr("dimnames") = dimn;
        retC[j] = cur;
      }
    } else {
      cur = wrap(cov);
      if (doDimNames) cur.attr("dimnames") = dimn;
      retC[j] = cur;
    }
  }
}

void foceiPhi(Environment e) {
  if (op_focei.neta==0) return;
  List retH(getRxNsub(rx));
  List retC(getRxNsub(rx));
  foceiPhiOne(e, retC, retH, 0);
  e["phiH"] = retH;
  e["phiC"] = retC;
}

SEXP foceiEtas(Environment e, bool bestMixEst=false) {
  if (op_focei.neta==0) return R_NilValue;
  int mixest = 0;
  if (op_focei.mixIdxN != 0) {
    mixest = 1;
  }
  List ret(op_focei.neta+2+mixest);
  CharacterVector nm(op_focei.neta+2+mixest);
  rx = getRxSolve_();
  IntegerVector ids(bestMixEst ? getRxNsub(rx) : getRxNsubAndMix(rx));
  NumericVector ofv(bestMixEst ? getRxNsub(rx) : getRxNsubAndMix(rx));
  IntegerVector mixesti(bestMixEst ? getRxNsub(rx) : getRxNsubAndMix(rx));
  int j,eta;
  for (j = op_focei.neta; j--;) {
    ret[j+1+mixest]=NumericVector(bestMixEst ? getRxNsub(rx) : getRxNsubAndMix(rx));
    nm[j+1+mixest] = "ETA[" + std::to_string(j+1) + "]";
  }
  NumericVector tmp;
  for (j=(bestMixEst ? (int)getRxNsub(rx) : (int)getRxNsubAndMix(rx)); j--;) {
    ids[j] = getRxId(j)+1;
    focei_ind *fInd = &(inds_focei[j]);
    // Update based on the best mix estimate when requested
    if (bestMixEst && op_focei.mixIdxN != 0) {
      int mixId = fInd->mixest[0]-1;
      if (mixId < 0) mixId = 0;  // guard: if inner OFV not run, default to mix 0
      fInd = &(inds_focei[getRxNsub(rx)*mixId + getRxId(j)]);
    }
    ofv[j] = -2*fInd->lik[0];
    if (mixest == 1) {
      // For bestMixEst: use the stored 1-indexed best mixture assignment
      // For full listing: getRxMixFromId already returns 1-indexed
      mixesti[j] = (bestMixEst ? inds_focei[j].mixest[0] : getRxMixFromId(j));
    }
    for (eta = op_focei.neta; eta--;) {
      tmp = ret[eta+1+mixest];
      // Save eta is what the ETAs are saved
      tmp[j] = fInd->saveEta[eta];
    }
  }
  if (e.exists("idLvl")) {
    RObject idl = e["idLvl"];
    if (idl.sexp_type() == STRSXP) {
      ids.attr("class") = "factor";
      ids.attr("levels") = idl;
    }
  }
  ret[0] = ids;
  nm[0] = "ID";
  if (mixest == 1) {
    ret[1] = mixesti;
    nm[1] = "MIXEST";
  }
  ret[op_focei.neta+1+mixest]=ofv;
  nm[op_focei.neta+1+mixest] = "OBJI";
  ret.attr("names") = nm;
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = IntegerVector::create(NA_INTEGER,
             bestMixEst ? -(int)getRxNsub(rx) : -(int)getRxNsubAndMix(rx));
  return(wrap(ret));
}


// R style optimfn
extern "C" double outerLikOpim(int n, double *par, void *ex){
  return(foceiOfv0(par));
}

// Gill 1983 Chat
static inline double Chat(double phi, double h, double epsA){
  if (phi == 0) return 2*epsA/(h*h);
  return 2*epsA/(h*fabs(phi));
}

static inline double ChatP(double phi, double h, double epsA){
  if (phi == 0) return 4*epsA/(h*h*h);
  return 4*epsA/(h*h*fabs(phi));
}

static inline double Phi(double fp, double f, double fn, double h){
  return (fp-2*f+fn)/(h*h);
}
static inline double phiC(double fp, double fn, double h){
  return (fp-fn)/(2*h);
}
static inline double phiF(double f, double fp, double h){
  return (fp-f)/h;
}
static inline double phiB(double f, double fn, double h){
  return (f-fn)/h;
}

// Call R function for gill83
double gillF=NA_REAL;
int gillThetaN=0;
Environment gillRfnE_;
// Fit environment + gate for the analytic ("fast") outer gradient.  The gate is
// TRUE only while the outer optimizer's gradient callback is live (set in
// foceiOuter around the optimizer switch), so foceiS's own numericGrad use and
// the exported foceiNumericGrad are unaffected.
Environment op_foceiFitEnv;
bool op_foceiFitEnvSet = false;
bool op_foceiUseAnalyticGrad = false;
Environment baseEnv = Environment::base_env();
Function doCall = baseEnv["do.call"];
Function gillRfn_ = baseEnv["invisible"];
int gillPar = 0;
double gillLong = false;
double gillRfn(double *theta){
  List par(1);
  NumericVector par0(gillThetaN);
  std::copy(&theta[0], &theta[0]+gillThetaN,par0.begin());
  par[0] = par0;
  NumericVector ret = as<NumericVector>(doCall(_["what"] = gillRfn_, _["args"]=par, _["envir"]=gillRfnE_));
  if (ret.size() == 1){
    return(ret[0]);
  } else {
    return(ret[gillPar]);
  }
}

static inline void gill83tickStep(int &k, int &K, int foceiGill) {
  if (foceiGill == 1 && op_focei.slow) {
    if (k < K) {
      op_focei.cur += (K-k);
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
    }
  }
}

// This is needed for shi21 h optimization
arma::vec shi21fnF(arma::vec &theta, int id) {
  updateTheta(theta.memptr());
  arma::vec ret(1);
  ret(0) = foceiOfv0(theta.memptr());
  if (op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
  return ret;
}

void gill83fnF(double *fp, double *theta, int, int foceiGill) {
  if (foceiGill == 1) {
    updateTheta(theta);
    *fp = foceiOfv0(theta);
    if (op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
  } else {
    *fp = gillRfn(theta);
  }
}

// *hf is the forward difference final estimate
// *hphif is central difference final estimate (when switching from forward to central differences)
// *df is the derivative estimate
// *df2 is the 2nd derivative estimate, useful for pre-conditioning.
// *ef is the err of the final estimate.
// *theta is the theta vector
// cpar is the parameter we are considering
// epsR is the relative error for the problem
// K is the maximum number of iterations before giving up on searching for the best interval.
// Returns 1 -- Success
//         2 -- Large error; Derivative estimate error 50% or more of the derivative
//         3 -- Function constant or nearly constant for this parameter
//         4 -- Function odd or nearly linear, df = K, df2 ~ 0
//         5 -- df2 increases rapidly as h decreases
int gill83(double *hf, double *hphif, double *df, double *df2, double *ef,
           double *theta, int cpar, double epsR, int K, double gillStep,
           double fTol, int cid, gill83fn_type gill83fn, int foceiGill,
           double gillF) {
  if (foceiGill == 1) op_focei.calcGrad=1;
  double f , x, hbar, h0, fp, fn=NA_REAL, phif, phib, phic, phicc = 0, phi, Chf, Chb,
    Ch, hs, hphi, hk, tmp, ehat, lasth, lastht=NA_REAL, lastfpt=NA_REAL, phict=NA_REAL;

  f = gillF;
  int k = 0;
  // Relative error should be given by the tolerances, I believe.
  double epsA=std::fabs(f)*epsR;
  x = theta[cpar];
  // FD1: // Initialization
  hbar = 2*(1+std::fabs(x))*_safe_sqrt(epsA/(1+std::fabs(f)));
  h0 = gillStep*hbar;
  lasth=h0;
  theta[cpar] = x + h0;
  gill83fn(&fp, theta, cid, foceiGill);
  theta[cpar] = x-h0;

  gill83fn(&fn, theta, cid, foceiGill);
  phif = phiF(f, fp, h0);
  phib = phiB(f, fn, h0);
  phic = phiC(fp, fn, h0);
  phi = Phi(fp, f, fn, h0);

  Chf = Chat(phif, h0, epsA);
  Chb = Chat(phib, h0, epsA);
  Ch  = ChatP(phi, h0, epsA);
  hs  = -1;
  hphi=hbar; // Not defined in Gill, but used for central difference switch if there are problems
  // FD2:  // Decide if to accept the interval
  hk = h0;
  if (max2(Chf, Chb) <= 0.1){
    hs=h0;
  }
  if (0.001 <= Ch && Ch <= 0.1){
    phicc=phic;
    hphi=h0;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      phict=phic;
      lastht  = lasth;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (Ch < 0.001){
    goto FD4;
  }
 FD3: // Increase h
  k++;
  hk=hk*gillStep;
  lasth=hk;
  // Compute the associated finite difference estimates and their
  // relative condition errors.
  theta[cpar] = x + hk;
  gill83fn(&fp, theta, cid, foceiGill);
  theta[cpar] = x-hk;
  gill83fn(&fn, theta, cid, foceiGill);
  phif = phiF(f, fp, hk);
  phib = phiB(f, fn, hk);
  phic = phiC(fp, fn, hk);
  phi = Phi(fp, f, fn, hk);
  Chf = Chat(phif, hk, epsA);
  Chb = Chat(phib, hk, epsA);
  Ch = ChatP(phi, hk, epsA);
  if (hs < 0 && max2(Chf, Chb) <= 0.1){
    hs = hk;
  }
  if (Ch <= 0.1){
    phicc=phic;
    hphi = hk;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (k == K) goto FD6;
  goto FD3;
 FD4: // Decrease h
  k++;
  hk=hk/gillStep;
  lasth=hk;
  // Compute the associated finite difference estimates and their
  // relative condition errors.
  theta[cpar] = x + hk;
  gill83fn(&fp, theta, cid, foceiGill);
  theta[cpar] = x-hk;
  gill83fn(&fn, theta, cid, foceiGill);
  phif = phiF(f, fp, hk);
  phib = phiB(f, fn, hk);
  tmp=phic;
  phic = phiC(fp, fn, hk);
  phi = Phi(fp, f, fn, hk);
  Chf = Chat(phif, hk, epsA);
  Chb = Chat(phib, hk, epsA);
  Ch = ChatP(phi, hk, epsA);
  if (Ch > .1){
    phicc=tmp;
    hphi=hk*gillStep; // hphi = h_k-1
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (max2(Chf, Chb) <= 0.1){
    hs = hk;
  }
  if (0.001 <= Ch && Ch <= 1){
    hphi = hk;
    if (fTol != 0 && fabs(phif) < fTol){
      lastfpt = fp;
      lastht  = lasth;
      phict=phic;
    }
    goto FD5;
  }
  if (fTol != 0 && fabs(phif) < fTol){
    lastfpt = fp;
    lastht  = lasth;
    phict=phic;
  }
  if (k == K) goto FD6;
  goto FD4;
 FD5: // Compute the estimate of the optimal interval
  *df2 = phi;
  *hf = 2*_safe_sqrt(epsA/fabs(phi));
  theta[cpar] = x + *hf;
  gill83fn(&fp, theta, cid, foceiGill);
  // Restore theta
  theta[cpar] = x;
  if (foceiGill == 1) updateTheta(theta);
  *df = phiF(f, fp, *hf);
  *ef = (*hf)*fabs(phi)/2+2*epsA/(*hf);
  *hphif=hphi;
  ehat = fabs(*df-phicc);
  if (max2(*ef, ehat) <= 0.5*(*df)){
    gill83tickStep(k, K, foceiGill);
    return 1;
  } else {
    // warning("The finite difference derivative err more than 50%% of the slope; Consider a different starting point.");
    if (!ISNA(lastht)){
      // Could be used;  Stick with the last below Ftol
      // *hf = lasth;
      // fp = lastfp;
      // *df = phiF(f, fp, *hf);
      // *df2=0;
      // // *df = 0.0; // Doesn't move.
      // *hphif=2*(*hf);
      // } else {
      *hf = lastht;
      fp = lastfpt;
      *df = phiF(f, fp, *hf);
      *df2=phic;
      // *df = 0.0; // Doesn't move.
      *hphif=phict;
    }
    gill83tickStep(k, K, foceiGill);
    return 2;
  }
 FD6: // Check unsatisfactory cases
  if (hs < 0){
    // F nearly constant.
    // Use sqrt(h0) as a last ditch effort.
    *hf = pow(DBL_EPSILON, 0.25);//hbar;
    // *df=phic;
    theta[cpar] = x + *hf;
    gill83fn(&fp, theta, cid, foceiGill);
    *df = phiF(f, fp, *hf);
    *df2=0;
    // *df = 0.0; // Doesn't move.
    *hphif=_safe_sqrt(h0);
    // warning("The surface around the initial estimate is nearly constant in one parameter grad=0.  Consider a different starting point.");
    gill83tickStep(k, K, foceiGill);
    return 3;
  }
  if (Ch > 0.1){ // Odd or nearly linear.
    *hf = h0;
    *df = phic;
    *df2 = 0;
    *ef = 2*epsA/(*hf);
    *hphif=hphi;
    // warning("The surface odd or nearly linear for one parameter; Check your function.");
    gill83tickStep(k, K, foceiGill);
    return 4;
  }
  // f'' is increasing rapidly as h decreases
  *hf = h0;
  *df = phic;
  *df2 = phi;
  *hphif=hphi;
  *ef = (*hf)*fabs(phi)/2+2*epsA/(*hf);
  // warning("The surface around the initial estimate is highly irregular in at least one parameter.  Consider a different starting point.");
  gill83tickStep(k, K, foceiGill);
  return 5;
}

// Calculate the mixture parameter gradient
//
// This notes that the mixture gradient does not need to be numerically, but
// can be calculated directly from the mixture probabilities, the translation from
// the by the mexpit, and the scaling factors.
//
// @param theta The parameter vector
//
// @param g The gradient vector to fill in
//
// @param cpar The parameter index to test/calculate the gradient for.
//
// @return 0 if the gradient was not calculated, 1 if it was.
//
int mixGrad(double *theta, double *g, int cpar) {
  if (op_focei.mixTrans == NULL) return 0;
  if (op_focei.mixTrans[cpar] != -1) {
    // This is a mixture grad
    int mi = op_focei.mixTrans[cpar];
    // First add the gradients from each individual contribution
    g[cpar] = 0.0;
    for (int i = 0; i < getRxNsub(rx); ++i) {
      focei_ind *fInd = &(inds_focei[i]);
      g[cpar] += fInd->mixProbGrad[mi];
    }
    // Next multiple the gradient from the mexpit() transformation
    // (from chain rule)
    g[cpar] *= op_focei.mixProbGrad[mi];
    // FIXME: Last apply the scaling gradient changes from chain rule
    double scaleTo = op_focei.scaleTo, C=getScaleC(cpar);
    switch (op_focei.scaleType){
    case 1: // normalized
      g[cpar] *= op_focei.c2;
      return 1;
      break;
    case 2: // log vs linear scales and/or ranges
      g[cpar] *= C;
      return 1;
      break;
    case 3: // simple multiplicative scaling
      if (op_focei.scaleTo != 0){
        g[cpar] *= op_focei.initPar[cpar]/scaleTo;
        return 1;
      } else {
        return 1;
      }
      break;
    case 4: // log non-log multiplicative scaling
      if (op_focei.scaleTo > 0){
        switch (op_focei.xPar[cpar]){
        case 1:
          return 1;
        default:
          g[cpar] *= op_focei.initPar[cpar]/scaleTo;
          return 1;
        }
      } else {
        return 1;
      }
    default:
      return 1;
    }
    return 0;

  }
  return 0;
}


// d(unscalePar)/d(x_i): the finite-difference outer gradient is d(OFV)/d(scaled
// par), so an analytic d(OFV)/d(theta) must be multiplied by this factor to land
// in the same optimizer scale.  Mirrors the linear coefficient of unscalePar().
static inline double dUnscaleParDx(int i) {
  double scaleTo = op_focei.scaleTo, C = getScaleC(i);
  switch (op_focei.scaleType) {
  case 1: return op_focei.c2;
  case 2: return C;
  case 3: return (op_focei.scaleTo != 0) ? op_focei.initPar[i]/scaleTo : 1.0;
  case 4:
    if (op_focei.scaleTo > 0) {
      return (op_focei.xPar[i] == 1) ? 1.0 : op_focei.initPar[i]/scaleTo;
    }
    return 1.0;
  default: return 1.0;
  }
}

static bool restoreFitSolve_();   // defined below (with covSolveArgs_)

// Read the per-subject etaP (neta x npars x nsub) the R analytic gradient stashed
// on the fit env (.foceiGradEtaP) into op_focei.getaP, converting each column to
// the SCALED optimizer parameterization (x dUnscaleParDx) so the Eq-48 step
// (op_focei.theta - etaPTheta) is in the same space.  Snapshots the scaled theta
// and sets etaPValid.  No-op (etaPValid=0) unless mceta is an Eq-48 mode.
static void loadAnalyticEtaP(double *theta) {
  op_focei.etaPValid = 0;
  if (op_focei.mceta != -2 && op_focei.mceta != -1) return;
  if (!op_foceiFitEnv.exists(".foceiGradEtaP")) return;
  RObject _epO = op_foceiFitEnv[".foceiGradEtaP"];
  if (Rf_isNull(_epO)) return;
  NumericVector _ep(_epO);
  int neta = op_focei.neta, npars = (int)op_focei.npars;
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  size_t need = (size_t)neta * npars * nsub;
  if (need == 0 || (size_t)_ep.size() != need) return;
  if (op_focei.getaP == NULL || op_focei.getaPn != need) {
    if (op_focei.getaP != NULL) R_Free(op_focei.getaP);
    op_focei.getaP = R_Calloc(need, double);
    op_focei.getaPn = need;
  }
  if (op_focei.etaPTheta == NULL) op_focei.etaPTheta = R_Calloc(npars, double);
  for (int k = 0; k < npars; k++) {
    double sc = dUnscaleParDx(k);
    for (int i = 0; i < nsub; i++) {
      for (int j = 0; j < neta; j++) {
        size_t idx = (size_t)j + (size_t)neta * ((size_t)k + (size_t)npars * i);
        op_focei.getaP[idx] = _ep[idx] * sc;
      }
    }
  }
  std::copy(theta, theta + npars, op_focei.etaPTheta);
  op_focei.etaPValid = 1;
}

// Analytic ("fast") outer gradient hook: fill g from .foceiCalcGradAnalytic when
// the gate is live.  Returns true on success (g filled), false to fall back to
// the finite-difference gradient -- ONLY when the sensitivity system could not
// be solved (a size mismatch is a wiring bug and stops).  theta is the current
// scaled optimizer point.
static bool analyticOuterGrad(double *theta, double *g) {
  if (!op_foceiUseAnalyticGrad || !op_foceiFitEnvSet) return false;
  op_focei.calcGrad = 1;
  // Ensure the inner solutions (eta*) and omega are current at this theta.
  foceiOfv0(theta);
  // Refresh the live state the R gradient reads; omega/theta/etaObf are
  // otherwise only written into the fit env at finalize.
  NumericVector _th(op_focei.ntheta);
  std::copy(&op_focei.fullTheta[0], &op_focei.fullTheta[0] + op_focei.ntheta,
            _th.begin());
  op_foceiFitEnv[".gradTheta"] = _th;
  op_foceiFitEnv["omega"] = getOmega();
  if (op_focei.mixIdxN != 0) {
    op_foceiFitEnv["etaObf"] = foceiEtas(op_foceiFitEnv, true);
  } else {
    op_foceiFitEnv["etaObf"] = foceiEtas(op_foceiFitEnv);
  }
  Environment _nlmixr2est = Environment::namespace_env("nlmixr2est");
  Function _agf = as<Function>(_nlmixr2est[".foceiCalcGradAnalytic"]);
  RObject _res = _agf(op_foceiFitEnv);
  // The augmented-sensitivity solves replaced the fit's global solve; restore it
  // so the next foceiOfv0 (objective/inner solve) reads the fit, not the last
  // augmented subject.  A failed restore -> fall back to FD (which re-solves).
  bool _restored = restoreFitSolve_();
  if (Rf_isNull(_res) || !_restored) {
    if (!op_focei.warnedAnalyticFallback) {
      op_focei.warnedAnalyticFallback = 1;
      Rf_warning("fast=TRUE: the analytic outer gradient could not be solved at this point; using the finite-difference gradient for the affected iteration(s)");
    }
    return false;
  }
  NumericVector _gv = as<NumericVector>(_res);
  if ((int)_gv.size() != (int)op_focei.npars) {
    // never silently degrade to FD on a length mismatch -- it means the R-side
    // parameter mapping disagrees with the optimizer's free-parameter set
    stop("fast=TRUE: analytic gradient length (%d) does not match the number of outer parameters (%d)",
         (int)_gv.size(), (int)op_focei.npars);
  }
  for (int i = 0; i < (int)op_focei.npars; i++) {
    if (!R_finite(_gv[i])) {
      // non-finite element = the sensitivity system did not solve cleanly
      if (!op_focei.warnedAnalyticFallback) {
        op_focei.warnedAnalyticFallback = 1;
        Rf_warning("fast=TRUE: the analytic outer gradient could not be solved at this point; using the finite-difference gradient for the affected iteration(s)");
      }
      return false;
    }
    g[i] = _gv[i] * dUnscaleParDx(i);
  }
  // Stash the per-subject EBE sensitivity etaP (Almquist Eq 46) that the R gradient
  // computed, for the Eq-48 warm-start extrapolation at the next outer step.
  loadAnalyticEtaP(theta);
  return true;
}

void numericGrad(double *theta, double *g){
  op_focei.mixDeriv=0;
  op_focei.reducedTol2=0;
  op_focei.curGill=0;
  op_focei.curAnalytic=0;
  if (analyticOuterGrad(theta, g)) {
    op_focei.curAnalytic=1;
    op_focei.nAnalyticGrad++;
    return;
  }
  if (op_foceiUseAnalyticGrad) op_focei.nFDGradFast++;
  if (op_focei.shi21maxOuter != 0 && op_focei.nF == 1) {
    clock_t t = clock() - op_focei.t0;
    int finalSlow = (op_focei.scale.every == 1) &&
      ((double)t)/CLOCKS_PER_SEC >= op_focei.gradProgressOfvTime;
    int maxiter = op_focei.shi21maxOuter;
    op_focei.totTick = op_focei.npars * maxiter * 2;
    op_focei.slow = (op_focei.scale.every == 1) &&
      ((double)t)/CLOCKS_PER_SEC*op_focei.totTick >= op_focei.gradProgressOfvTime;
    // Use Shi Difference
    if (op_focei.slow) {
      RSprintf(_("calculate Shi21 Difference and optimize forward difference step size:\n"));
      op_focei.t0 = clock();
      op_focei.cur=0;
      op_focei.curTick=0;
    }
    arma::vec grFinal(1);
    arma::vec f0(1);
    f0(0) = op_focei.lastOfv;
    arma::vec armaTheta(op_focei.npars);
    std::copy(theta, theta+op_focei.npars, armaTheta.begin());
    double h = 0;
    for (int cpar = (int)op_focei.npars; cpar--;) {
      if (mixGrad(theta, g, cpar) == 1) {
        continue;
      } else {
        op_focei.calcGrad=1;
        op_focei.aEps[cpar] = shi21Forward(shi21fnF, armaTheta, h,
                                           f0, grFinal, 0, cpar,
                                           op_focei.hessEpsInner, //double ef = 7e-7,
                                           1.5,  //double rl = 1.5,
                                           6.0,  //double ru = 6.0);;
                                           maxiter);  //maxiter=15
        op_focei.aEpsC[cpar] = op_focei.aEps[cpar];
        g[cpar] = grFinal(0);
      }
    }
    if (op_focei.slow) {
      op_focei.cur=op_focei.totTick;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      RSprintf("\n");
    }
    op_focei.calcGrad=0;
    op_focei.curGill=2;
    op_focei.slow = finalSlow;
  } else if ((op_focei.repeatGill == 1 || op_focei.nF == 1) && op_focei.gillK > 0) {
    clock_t t = clock() - op_focei.t0;
    int finalSlow = (op_focei.scale.every == 1) &&
      ((double)t)/CLOCKS_PER_SEC >= op_focei.gradProgressOfvTime;

    op_focei.repeatGill=0;
    op_focei.reducedTol2=0;
    double hf, hphif, err;
    op_focei.slow = (op_focei.scale.every == 1) &&
      ((double)t)/CLOCKS_PER_SEC*op_focei.npars * op_focei.gillK >= op_focei.gradProgressOfvTime;
    if (op_focei.slow){
      op_focei.cur = 0;
      op_focei.totTick = op_focei.npars * op_focei.gillK;
      op_focei.t0 = clock();
      if (op_focei.repeatGillN != 0){
        RSprintf(_("repeat %d Gill diff/forward difference step size:\n"),
                 op_focei.repeatGillN);
      } else {
        RSprintf(_("calculate Gill Difference and optimize forward difference step size:\n"));
      }
    }
    for (int cpar = (int)op_focei.npars; cpar--;) {
      if (mixGrad(theta, g, cpar) == 1) {
        continue;
      } else {
        op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
                                        theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
                                        -1, gill83fnG, 1, op_focei.lastOfv);
        err = 1/(std::fabs(theta[cpar])+1);
        if (op_focei.gillDf[cpar] == 0){
          op_focei.scaleC[cpar]=op_focei.scaleC0;
          op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
                                          theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
                                          -1, gill83fnG, 1, op_focei.lastOfv);
          if (op_focei.gillDf[cpar] == 0){
            op_focei.scaleC[cpar]=1/op_focei.scaleC0;
            op_focei.gillRet[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
                                            theta, cpar, op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
                                            -1, gill83fnG, 1, op_focei.lastOfv);
          }
        }
        // h=aEps*(|x|+1)/sqrt(1+fabs(f));
        // h*sqrt(1+fabs(f))/(|x|+1) = aEps
        // let err=2*sqrt(epsA/(1+f))
        // err*(aEps+|x|rEps) = h
        // Let aEps = rEps (could be a different ratio)
        // h/err = aEps(1+|x|)
        // aEps=h/err/(1+|x|)
        //
        op_focei.aEps[cpar]  = hf*err;
        op_focei.rEps[cpar]  = hf*err;
        if(op_focei.optGillF){
          op_focei.aEpsC[cpar] = hf*err;
          op_focei.rEpsC[cpar] = hf*err;
        } else {
          op_focei.aEpsC[cpar] = hphif*err;
          op_focei.rEpsC[cpar] = hphif*err;
        }
        g[cpar] = op_focei.gillDf[cpar];
      }
    }
    if(op_focei.slow){
      op_focei.cur=op_focei.totTick;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      RSprintf("\n");
    }
    op_focei.didGill=1;
    if (op_focei.reducedTol2 && op_focei.repeatGillN < op_focei.repeatGillMax){
      op_focei.repeatGill=1;
      op_focei.repeatGillN++;
      op_focei.reducedTol2=0;
    }
    op_focei.curGill=1;
    op_focei.slow = finalSlow;
  } else {
    if(op_focei.slow){
      op_focei.t0 = clock();
      op_focei.cur=0;
      op_focei.curTick=0;
      op_focei.totTick = op_focei.npars * 2;
    }
    op_focei.calcGrad=1;
    rx = getRxSolve_();
    int npars = (int)op_focei.npars;
    int cpar;
    double cur, delta, tmp, tmp0=NA_REAL;
    double f=0;
    // Do Forward difference if the OBJF for *theta has already been calculated.
    bool doForward=false;
    if (op_focei.derivMethod == 0) {
      doForward=true;
      // If the first derivative wasn't calculated, then calculate it.
      for (cpar = npars; cpar--;){
        if (theta[cpar] != op_focei.theta[cpar]){
          doForward=false;
          break;
        }
      }
      if (doForward){
        // Fill in lik0
        f=op_focei.lastOfv;
      } else {
        op_focei.calcGrad=0; // Save OBF and theta
        f = foceiOfv0(theta);
        if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        op_focei.calcGrad=1;
        doForward=true;
      }
    }
    for (cpar = npars; cpar--;) {
      if (mixGrad(theta, g, cpar) == 1) {
        continue;
      } else {
        if (doForward){
          delta = (std::fabs(theta[cpar])*op_focei.rEps[cpar] + op_focei.aEps[cpar]);
        } else {
          delta = (std::fabs(theta[cpar])*op_focei.rEpsC[cpar] + op_focei.aEpsC[cpar]);
        }
        cur = theta[cpar];
        theta[cpar] = cur + delta;
        if (doForward) {
          tmp = foceiOfv0(theta);
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          g[cpar] = (tmp-f)/delta;
        } else {
          tmp0 = foceiOfv0(theta);
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          theta[cpar] = cur - delta;
          tmp = foceiOfv0(theta);
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          g[cpar] = (tmp0-tmp)/(2*delta);
        }
        if (doForward && fabs(g[cpar]) > op_focei.gradCalcCentralLarge){
          doForward = false;
          theta[cpar] = cur - delta;
          g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          op_focei.mixDeriv=1;
        }
        if (std::isnan(g[cpar]) ||  ISNA(g[cpar]) || !R_FINITE(g[cpar])){
          if (doForward){
            // Switch to Backward difference method
            op_focei.mixDeriv=1;
            theta[cpar] = cur - delta;
            g[cpar] = (f-foceiOfv0(theta))/(delta);
            if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            if (R_FINITE(op_focei.gradTrim)){
              if (g[cpar] > op_focei.gradTrim){
                g[cpar]=op_focei.gradTrim;
              } else if (g[cpar] < op_focei.gradTrim){
                g[cpar]=-op_focei.gradTrim;
              }
            }
          } else {
            // We are using the central difference AND there is an NA in one of the terms
            // g[cpar] = (tmp0-tmp)/(2*delta);
            op_focei.mixDeriv=1;
            if (std::isnan(tmp0) || ISNA(tmp0) || !R_FINITE(tmp0)){
              // Backward
              g[cpar] = (f-tmp)/delta;
            } else {
              // Forward
              g[cpar] = (tmp0-f)/delta;
            }
            if (R_FINITE(op_focei.gradTrim)){
              if (g[cpar] > op_focei.gradTrim){
                g[cpar]=op_focei.gradTrim;
              } else if (g[cpar] < op_focei.gradTrim){
                g[cpar]=-op_focei.gradTrim;
              }
            }
          }
        }
        if (R_FINITE(op_focei.gradTrim) && g[cpar] > op_focei.gradTrim){
          if (doForward){
            op_focei.mixDeriv=1;
            theta[cpar] = cur - delta;
            g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
            if(op_focei.slow)  op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            if (g[cpar] > op_focei.gradTrim){
              g[cpar]=op_focei.gradTrim;
            } else if (g[cpar] < op_focei.gradTrim){
              g[cpar]=-op_focei.gradTrim;
            }
          } else {
            g[cpar]=op_focei.gradTrim;
          }
        } else if (R_FINITE(op_focei.gradTrim) && g[cpar] < -op_focei.gradTrim){
          if (doForward){
            op_focei.mixDeriv=1;
            theta[cpar] = cur - delta;
            g[cpar] = (tmp-foceiOfv0(theta))/(2*delta);
            if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            if (g[cpar] > op_focei.gradTrim){
              g[cpar]=op_focei.gradTrim;
            } else if (g[cpar] < op_focei.gradTrim){
              g[cpar]=-op_focei.gradTrim;
            }
          } else {
            g[cpar]=-op_focei.gradTrim;
          }
        } else if (doForward && fabs(g[cpar]) < op_focei.gradCalcCentralSmall){
          op_focei.mixDeriv = 1;
          theta[cpar]       = cur - delta;
          tmp = g[cpar];
          g[cpar]           = (tmp-foceiOfv0(theta))/(2*delta);
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          if (fabs(tmp) > fabs(g[cpar])) g[cpar] = tmp;
        } else if (doForward) {
          if(op_focei.slow) op_focei.curTick = par_progress(op_focei.cur++, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        }
        theta[cpar] = cur;
        // zero gradients mean reset
        if (g[cpar] == 0.0)  {
          op_focei.zeroGrad = 1;
          break;
        }
        if (std::isnan(g[cpar]) ||  ISNA(g[cpar]) || !R_FINITE(g[cpar])){
          op_focei.zeroGrad = 1;
          break;
        }
      }
    }
    if(op_focei.slow) {
      op_focei.cur=op_focei.totTick;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      RSprintf("\n");
    }
    op_focei.calcGrad=0;
  }
}

//[[Rcpp::export]]
NumericVector foceiNumericGrad(NumericVector theta){
  NumericVector ret(theta.size());
  numericGrad(&theta[0], &ret[0]);
  return ret;
}

////////////////////////////////////////////////////////////////////////////////
// Setup FOCEi functions
CharacterVector rxParams_(const RObject &obj);


static inline void foceiSetupTrans_(CharacterVector pars){
  unsigned int k, j,  ps = pars.size();
  k=ps;
  std::string thetaS;
  std::string etaS;
  std::string cur;
  if (op_focei.etaTrans != NULL) {
    R_Free(op_focei.etaTrans);
  }
  // One contiguous int buffer for the per-parameter index arrays.
  // Layout: etaTrans(neta) | nbdInner(neta) | xPar(N) | probitIdxArr(N)
  //         | thetaTrans(N) | fixedTrans(N) | etaFD(neta)
  //         | (gap of omegan when mixIdxN) | mixTrans(N)
  // where N = ntheta+omegan.  The (mixIdxN ? 4 : 3) ntheta+omegan
  // slots become (mixIdxN ? 5 : 4) once probitIdxArr is included.
  op_focei.etaTrans    = R_Calloc(op_focei.neta*3 +
                                  (4+(op_focei.mixIdxN != 0))*(op_focei.ntheta + op_focei.omegan), int); //[neta]
  op_focei.nbdInner    = op_focei.etaTrans + op_focei.neta;
  op_focei.xPar        = op_focei.nbdInner + op_focei.neta; // [ntheta+nomega]
  op_focei.probitIdxArr = op_focei.xPar + op_focei.ntheta + op_focei.omegan; // [ntheta+nomega]
  op_focei.thetaTrans  = op_focei.probitIdxArr + op_focei.ntheta + op_focei.omegan; // [ntheta+nomega]
  op_focei.fixedTrans  = op_focei.thetaTrans + op_focei.ntheta + op_focei.omegan; // [ntheta + nomega]
  op_focei.etaFD       = op_focei.fixedTrans + op_focei.ntheta + op_focei.omegan; // [neta]
  if (op_focei.mixIdxN) {
    op_focei.mixTrans    = op_focei.etaFD + op_focei.neta; // [ntheta+omegan]
  } else {
    op_focei.mixTrans    = NULL;
  }

  if (op_focei.fullTheta != NULL) R_Free(op_focei.fullTheta);
  op_focei.fullTheta   = R_Calloc(4*(op_focei.ntheta+op_focei.omegan) +
                                  2*(_aqn*op_focei.neta), double); // [ntheta+omegan]
  op_focei.theta       = op_focei.fullTheta+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]
  op_focei.initPar     = op_focei.theta+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]
  op_focei.scaleC      = op_focei.initPar+op_focei.ntheta+op_focei.omegan; // [ntheta + omegan]
  op_focei.aqx         = op_focei.scaleC + op_focei.ntheta+op_focei.omegan; // [aqn*neta]
  op_focei.aqw         = op_focei.aqx + _aqn*op_focei.neta; // [aqn*neta]

  unsigned int neta = 0;
  unsigned int ntheta = 0;
  for (;k--;){
    for (j = ps; j--;) {
      // Compare ETAS first since they are smaller strings.
      cur = as<std::string>(pars[k]);
      etaS = "ETA[" + std::to_string(j+1) + "]";
      if (cur == etaS) {
        op_focei.etaTrans[j] = k;
        neta++;
        break;
      } else {
        thetaS = "THETA[" + std::to_string(j+1) + "]";
        if (cur == thetaS){
          op_focei.thetaTrans[j] = k;
          ntheta++;
          break;
        }
      }
    }
  }
  if (op_focei.ntheta != ntheta) {
    stop("theta mismatch op_focei.ntheta %d, ntheta: %d\n", op_focei.ntheta, ntheta);
  }
  if (op_focei.neta != neta) {
    stop("eta mismatch op_focei.neta %d, neta: %d\n", op_focei.neta, neta);
  }
  op_focei.nzm = (op_focei.neta + 1) * (op_focei.neta + 2) / 2 + (op_focei.neta + 1)*6+1;
}

static inline void foceiSetupMixTrans(int k, int j) {
  if (!op_focei.mixIdxN) return;
  op_focei.mixTrans[k] = -1;
  // op_focei.mixIdx is allocated later in foceiSetup_; skip the lookup
  // when it hasn't been filled yet (the caller will re-run after it is set).
  if (op_focei.mixIdx == NULL) return;
  for (unsigned int m = 0; m < op_focei.mixIdxN; ++m) {
    if (op_focei.mixIdx[m] - 1 == j) {
      op_focei.mixTrans[k] = m;
      break;
    }
  }
}

static inline void foceiSetupTheta_(List mvi,
                                    NumericVector theta,
                                    Nullable<LogicalVector> thetaFixed,
                                    double scaleTo,
                                    bool alloc){
  // Get the fixed thetas
  // fixedTrans gives the theta->full theta translation
  // initPar is the initial parameters used for parameter scaling.
  op_focei.scaleTo=scaleTo;
  int thetan = theta.size();
  int omegan;
  NumericVector omegaTheta;
  if (op_focei.neta == 0) {
    omegan=0;
  } else {
    omegan = getOmegaN();
    omegaTheta = getOmegaTheta();
  }
  int fixedn = 0;
  int j;
  LogicalVector thetaFixed2;
  if (!thetaFixed.isNull()){
    thetaFixed2 = as<LogicalVector>(thetaFixed);
    for (j = thetaFixed2.size(); j--;){
      if (thetaFixed2[j]) fixedn++;
    }
  } else {
    thetaFixed2 =LogicalVector(0);
  }
  // Mu-group thetas are excluded from the outer optimizer's free-parameter set
  // separately from the user's own ini(...~fix) (see includeMuGroupThetas for
  // the final-Hessian-pass unlock counterpart); skip already-fixed thetas to
  // avoid double-counting.
  auto isMuGroupSkip = [&](int jj) -> bool {
    return op_focei.muModel != 0 && op_focei.muGroupThetaSkip != NULL &&
      jj < thetan && op_focei.muGroupThetaSkip[jj] != 0;
  };
  for (j = thetan; j--;){
    if (isMuGroupSkip(j) && !(j < thetaFixed2.size() && thetaFixed2[j])) fixedn++;
  }
  int npars = thetan+omegan-fixedn;
  if (alloc){
    rxUpdateFuns(as<SEXP>(mvi["trans"]), &rxInner);
    foceiSetupTrans_(as<CharacterVector>(mvi["params"]));
    // Locate rx_pred_ in the inner lhs (AR(1) lag defs may precede it).
    op_focei.predOffset = 0;
    CharacterVector innerLhs = as<CharacterVector>(mvi["lhs"]);
    for (int il = 0; il < innerLhs.size(); ++il) {
      if (as<std::string>(innerLhs[il]) == "rx_pred_") { op_focei.predOffset = il; break; }
    }
  } else if (!op_focei.alloc){
    stop("FOCEi problem not allocated\nThis can happen when symengine<->nlmixr2 interaction is not working correctly.");
  }
  std::copy(theta.begin(), theta.end(), &op_focei.fullTheta[0]);
  std::copy(omegaTheta.begin(), omegaTheta.end(), &op_focei.fullTheta[0]+thetan);
  if ((int)op_focei.ntheta != (int)thetan){
    rxOptionsFreeFocei();
    stop("op_focei.ntheta(%d) != thetan(%d)", op_focei.ntheta, thetan);
  }
  op_focei.ntheta = (unsigned int)thetan;
  op_focei.omegan = (unsigned int)omegan;
  int k = 0;
  for (j = 0; j < thetan+omegan; j++){
    if (isMuGroupSkip(j)) continue;
    if (j < thetaFixed2.size() && !thetaFixed2[j]){
      if (j < theta.size()){
        op_focei.initPar[k] = theta[j];
      } else if (j < theta.size() + omegan){
        op_focei.initPar[k] = omegaTheta[j-theta.size()];
      }
      foceiSetupMixTrans(k, j);
      op_focei.fixedTrans[k++] = j;
    } else if (j >= thetaFixed2.size()){
      if (j < theta.size()){
        foceiSetupMixTrans(k, j);
        op_focei.initPar[k] = theta[j];
        op_focei.fixedTrans[k++] = j;
      } else if (j < theta.size() + omegan){
        foceiSetupMixTrans(k, j);
        op_focei.initPar[k] = omegaTheta[j-theta.size()];
        op_focei.fixedTrans[k++] = j;
      }
    }
  }
  op_focei.npars  = (unsigned int)npars;
}

static inline void foceiSetupNoEta_(){

  // Mixtures only work in population only models;
  rx = getRxSolve_();

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei = R_Calloc(getRxNsub(rx), focei_ind);
  op_focei.gEtaGTransN=(op_focei.neta)*getRxNsub(rx);

  if (op_focei.gthetaGrad != NULL && op_focei.mGthetaGrad) R_Free(op_focei.gthetaGrad);
  op_focei.gthetaGrad = R_Calloc((size_t)op_focei.gEtaGTransN + (size_t)getRxNall(rx), double);
  op_focei.llikObsFull = op_focei.gthetaGrad + op_focei.gEtaGTransN; // [getRxNall(rx)]
  std::fill_n(op_focei.llikObsFull, getRxNall(rx), NA_REAL);
  op_focei.mGthetaGrad = true;
  focei_ind *fInd;
  int jj = 0, iLO=0;
  for (int i = getRxNsub(rx); i--;){
    fInd = &(inds_focei[i]);
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, i);
    fInd->doChol=!(op_focei.cholSEOpt);
    fInd->doFD=0;
    // ETA ini
    fInd->eta = NULL;
    fInd->oldEta = NULL;
    fInd->tryEta = NULL;
    fInd->saveEta = NULL;
    fInd->g = NULL;
    fInd->x = NULL;
    fInd->var = NULL;
    fInd->lp = NULL;
    fInd->Vid = NULL;
    fInd->a = NULL;
    fInd->c = NULL;
    fInd->B = NULL;
    fInd->zm = NULL;
    fInd->zmH = NULL;
    fInd->zmEta = NULL;
    fInd->zmValid = 0;
    fInd->thetaGrad = &op_focei.gthetaGrad[jj];
    jj+= op_focei.npars;
    fInd->mode = 1;
    fInd->uzm = 1;
    fInd->doEtaNudge=0;
    // llikObs
    fInd->llikObs = &op_focei.llikObsFull[iLO];
    iLO += getIndNallTimes(ind);
  }
  op_focei.alloc=true;
}

static inline void foceiSetupEta_(NumericMatrix etaMat0){
  rx = getRxSolve_();
  // Guard: nall_mix^2 is used in the etaUpper allocation below. When
  // nall_mix > 65535 the squared term exceeds the range of a 32-bit integer,
  // and on 32-bit size_t platforms the product would overflow.
  {
    size_t nall_mix_chk = getRxNallAndMix(rx);
    if (nall_mix_chk > 65535) { // nocov
      stop("focei: dataset too large for this mixture model configuration " // nocov
               "(getRxNall * (mixIdxN+1) = %zu would produce an infeasibly large " // nocov
               "allocation). Reduce the number of observations or mixture components.", // nocov
               nall_mix_chk); // nocov
    } // nocov
  }

  if (inds_focei != NULL) R_Free(inds_focei);
  inds_focei = R_Calloc(getRxNsubAndMix(rx), focei_ind);
  // Size the FOCE eta=0 population-R cache alongside inds_focei (single-threaded
  // here); gen = -1 forces a compute on each subject's first inner call.
  _foceRPopCache.assign(getRxNsubAndMix(rx), arma::vec());
  _foceRPopGen.assign(getRxNsubAndMix(rx), -1L);
  _foceRPopCurGen = 0;
  RObject etaMat0s = transpose(etaMat0);
  double *etaMat0d = REAL(etaMat0s);
  {
    size_t _gEtaGTransN = (op_focei.neta + 1) * getRxNsubAndMix(rx);
    if (_gEtaGTransN > UINT_MAX) { // nocov
      stop("focei: neta * nsub_mix too large for gEtaGTransN"); // nocov
    } // nocov
    op_focei.gEtaGTransN = (unsigned int)_gEtaGTransN;
  }
  size_t nz = ((op_focei.neta+1)*(op_focei.neta+2)/2 + 6*(op_focei.neta+1) + 1) *
               getRxNsubAndMix(rx);

  if (op_focei.etaUpper != NULL) R_Free(op_focei.etaUpper);

  {
    size_t nall_mix = getRxNallAndMix(rx);
    size_t nsub_mix = getRxNsubAndMix(rx);
    op_focei.etaUpper = R_Calloc(
      (size_t)op_focei.gEtaGTransN * 10 +
      op_focei.npars * (nsub_mix + 1) + nz +
      2 * op_focei.neta * nall_mix +
      nall_mix +
      nall_mix * nall_mix +
      op_focei.neta * 6 +
      2 * op_focei.neta * op_focei.neta * nsub_mix +
      nall_mix +
      // gZmH (packed lower-tri Hessian) + gZmEta per subject for warm="calc"
      nsub_mix * ((size_t)op_focei.neta * (op_focei.neta + 1) / 2 + op_focei.neta) +
      // per-obs censored inner-Hessian coefficients gcHff/gcHfr/gcHrr
      3 * nall_mix,
      double);
  }
  op_focei.etaLower =  op_focei.etaUpper + op_focei.neta;
  op_focei.geta     = op_focei.etaLower + op_focei.neta;
  op_focei.gtryEta  = op_focei.geta + op_focei.neta;
  op_focei.goldEta  = op_focei.gtryEta + op_focei.gEtaGTransN;
  op_focei.getahf   = op_focei.goldEta + op_focei.gEtaGTransN;
  op_focei.getahr   = op_focei.getahf + op_focei.gEtaGTransN;
  op_focei.getahh   = op_focei.getahr + op_focei.gEtaGTransN;
  op_focei.gsaveEta = op_focei.getahh + op_focei.gEtaGTransN;
  op_focei.gG       = op_focei.gsaveEta + op_focei.gEtaGTransN;
  op_focei.gVar     = op_focei.gG + op_focei.gEtaGTransN;
  op_focei.gX       = op_focei.gVar + op_focei.gEtaGTransN;
  op_focei.glp      = op_focei.gX + op_focei.gEtaGTransN;
  op_focei.gthetaGrad = op_focei.glp + op_focei.gEtaGTransN;  // op_focei.npars*(getRxNsub(rx) + 1)
  op_focei.gZm      = op_focei.gthetaGrad + op_focei.npars*(getRxNsubAndMix(rx) + 1); // nz
  op_focei.ga       = op_focei.gZm + nz;//[op_focei.neta * getRxNall(rx)]
  op_focei.gc       = op_focei.ga + op_focei.neta * getRxNallAndMix(rx);//[op_focei.neta * getRxNall(rx)]
  op_focei.gB       = op_focei.gc + op_focei.neta * getRxNallAndMix(rx);//[getRxNall(rx)]
  op_focei.gH       = op_focei.gB + getRxNallAndMix(rx); //[op_focei.neta*op_focei.neta*getRxNsub(rx)]
  op_focei.llikObsFull = op_focei.gH + op_focei.neta*op_focei.neta*getRxNsubAndMix(rx); // [getRxNall(rx)]
  op_focei.gVid     = op_focei.llikObsFull + getRxNallAndMix(rx);
  op_focei.gZmH     = op_focei.gVid + (size_t)getRxNallAndMix(rx)*getRxNallAndMix(rx); //[neta*(neta+1)/2 * nsub]
  op_focei.gZmEta   = op_focei.gZmH + (size_t)(op_focei.neta*(op_focei.neta+1)/2)*getRxNsubAndMix(rx); //[neta*nsub]
  op_focei.gcHff    = op_focei.gZmEta + (size_t)op_focei.neta*getRxNsubAndMix(rx); //[nall_mix]
  op_focei.gcHfr    = op_focei.gcHff + getRxNallAndMix(rx); //[nall_mix]
  op_focei.gcHrr    = op_focei.gcHfr + getRxNallAndMix(rx); //[nall_mix]
  // Could use .zeros() but since I used Calloc, they are already zero.
  // Yet not doing it causes the theta reset error.
  op_focei.etaM     = mat(op_focei.neta, 1, arma::fill::zeros);
  op_focei.etaS     = mat(op_focei.neta, 1, arma::fill::zeros);
  op_focei.eta1SD   = mat(op_focei.neta, 1, arma::fill::zeros);
  op_focei.n        = 1.0;

  // Prefill to 0.1 or 10%
  std::fill_n(&op_focei.gVar[0], op_focei.gEtaGTransN, 0.1);
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal


  unsigned int i, j = 0, k = 0, ii=0, jj = 0, iA=0, iB=0, iH=0, iVid=0, iLO=0;
  focei_ind *fInd;
  for (i = getRxNsubAndMix(rx); i--;){
    fInd = &(inds_focei[i]);
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(i));

    // Guard against null-pointer arithmetic (UBSan): mixIdx/mixProb are
    // allocated later in foceiSetup_; when mixIdxN==0 these pointers are
    // never dereferenced, so NULL is the correct sentinel value.
    if (op_focei.mixIdx != NULL) {
      fInd->mixest = op_focei.mixIdx + op_focei.mixIdxN + getRxId(i);
    } else {
      fInd->mixest = NULL;
    }
    if (op_focei.mixProb != NULL) {
      fInd->mixProb    = op_focei.mixProb +
        (getRxId(i) + 1)*(op_focei.mixIdxN + 1);
      fInd->mixProbGrad = op_focei.mixProbGrad +
        (getRxId(i) + 1)*(op_focei.mixIdxN);
    } else {
      fInd->mixProb    = NULL;
      fInd->mixProbGrad = NULL;
    }

    fInd->doChol=!(op_focei.cholSEOpt);
    fInd->doFD = 0;
    // ETA ini
    fInd->eta = &op_focei.geta[j];
    fInd->etahf = &op_focei.getahf[j];
    fInd->etahr = &op_focei.getahr[j];
    fInd->etahh = &op_focei.getahh[j];
    fInd->oldEta = &op_focei.goldEta[j];
    fInd->tryEta = &op_focei.gtryEta[j];
    fInd->saveEta = &op_focei.gsaveEta[j];
    fInd->g = &op_focei.gG[j];
    fInd->x = &op_focei.gX[j];
    fInd->var = &op_focei.gVar[j];
    fInd->lp = &op_focei.glp[j];
    fInd->Vid = &op_focei.gVid[iVid];
    iH += op_focei.neta*op_focei.neta;
    iVid += (getIndNallTimes(ind) -
             getIndNdoses(ind) -
             getIndNevid2(ind)
             )*(getIndNallTimes(ind) -
                getIndNdoses(ind) -
                getIndNevid2(ind));
    fInd->llikObs = &op_focei.llikObsFull[iLO];
    iLO += getIndNallTimes(ind);

    // Copy in etaMat0 to the inital eta stored (0 if unspecified)
    // std::copy(&etaMat0[i*op_focei.neta], &etaMat0[(i+1)*op_focei.neta], &fInd->saveEta[0]);
    if (etaMat0.nrow() == (int)getRxNsubAndMix(rx)) {
      std::copy(&etaMat0d[i*op_focei.neta], &etaMat0d[(i+1)*op_focei.neta], &fInd->eta[0]);
    } else {
      int origId = getRxId(i);
      std::copy(&etaMat0d[origId*op_focei.neta], &etaMat0d[(origId+1)*op_focei.neta], &fInd->eta[0]);
    }

    fInd->eta[op_focei.neta] = i;
    fInd->saveEta[op_focei.neta] = i;
    fInd->oldEta[op_focei.neta] = i;
    fInd->tryEta[op_focei.neta] = i;

    j+=op_focei.neta+1;

    k+=op_focei.neta;

    fInd->a = &op_focei.ga[iA];
    fInd->c = &op_focei.gc[iA];
    iA += op_focei.neta * (getIndNallTimes(ind) -
                           getIndNdoses(ind) -
                           getIndNevid2(ind));

    fInd->B = &op_focei.gB[iB];
    fInd->cHff = &op_focei.gcHff[iB];
    fInd->cHfr = &op_focei.gcHfr[iB];
    fInd->cHrr = &op_focei.gcHrr[iB];
    iB += (getIndNallTimes(ind) -
           getIndNdoses(ind) -
           getIndNevid2(ind));

    fInd->zm = &op_focei.gZm[ii];
    ii+= (op_focei.neta+1) * (op_focei.neta + 2) / 2 +
      6*(op_focei.neta + 1)+1;

    fInd->zmH = op_focei.gZmH + (size_t)i*(op_focei.neta*(op_focei.neta+1)/2);
    fInd->zmEta = op_focei.gZmEta + (size_t)i*op_focei.neta;
    fInd->zmValid = 0;

    fInd->thetaGrad = &op_focei.gthetaGrad[jj];
    jj+= op_focei.npars;

    fInd->mode = 1;
    fInd->uzm = 1;
    fInd->doEtaNudge=1;
  }
  op_focei.thetaGrad = &op_focei.gthetaGrad[jj];
  op_focei.alloc=true;
}

extern "C" double foceiOfvOptim(int n, double *x, void *ex);
extern "C" void outerGradNumOptim(int n, double *par, double *gr, void *ex);

// covType="analytic" FD fallback: the augmented cov solves run through
// rxode2::rxSolve, which calls rxSolveFree() and frees the fit's global solve.
// Stash foceiSetup_'s rxSolve_ setup args so foceiCalcR can re-run them and restore
// the fit solve before the finite-difference Hessian.
static SEXP covSolveArgs_ = R_NilValue;
// est="impmap": theta-sensitivity model that sizes the shared pool (the largest
// structure), or R_NilValue to size for the inner model as usual.
SEXP _impPoolModel = R_NilValue;
static void storeCovSolveArgs_(SEXP obj, SEXP rxControl, SEXP params, SEXP data) {
  List L = List::create(obj, rxControl, params, data);
  if (covSolveArgs_ != R_NilValue) R_ReleaseObject(covSolveArgs_);
  R_PreserveObject(L);
  covSolveArgs_ = L;
}
// release the preserved model+rxControl+params+data once the cov step is done
// so an analytic fit doesn't pin its whole dataset until the next analytic fit
static void releaseCovSolveArgs_() {
  if (covSolveArgs_ != R_NilValue) {
    R_ReleaseObject(covSolveArgs_);
    covSolveArgs_ = R_NilValue;
  }
}
static bool restoreFitSolve_() {
  if (covSolveArgs_ == R_NilValue) return false;
  try {
    List args = as<List>(covSolveArgs_);
    RObject obj = args[0]; List rxControl = as<List>(args[1]);
    RObject params = args[2]; RObject data = args[3];
    rxode2::rxSolve_(obj, rxControl, R_NilValue, R_NilValue, params, data, R_NilValue, 1);
    rx = getRxSolve_();
    return true;
  } catch (Rcpp::internal::InterruptedException&) {
    throw;   // Ctrl-C: Rcpp requires the interrupt to propagate, never swallow it
  } catch (Rcpp::LongjumpException&) {
    throw;   // Rcpp longjump protocol must be rethrown intact, not turned into false
  } catch (...) {
    return false;   // genuine solve error -> report a failed restore
  }
}

// RAII: tighten the ODE solve tolerances to covSolveTol for the covariance-step
// finite-difference solves, restoring the fit's tolerances on exit.  No-op unless the
// user set foceiControl(covSolveTol=); the analytic R-matrix applies covSolveTol on its
// own augmented solves (.foceiAnalyticSolveTol).
struct CovSolveTolGuard {
  bool active = false;
  double savAtol = NA_REAL, savRtol = NA_REAL;
  CovSolveTolGuard(Environment e) {
    if (!e.exists("control")) return;
    List ctl = as<List>(e["control"]);
    if (!ctl.containsElementNamed("covSolveTol")) return;
    RObject cst = ctl["covSolveTol"];
    if (cst.isNULL() || Rf_length(cst) < 1) return;
    double tol = as<double>(cst);
    if (!R_FINITE(tol) || tol <= 0) return;
    rxGetSolveAtolRtol(&savAtol, &savRtol);
    if (!R_FINITE(savAtol) || !R_FINITE(savRtol)) return;   // no live solve to retune
    rxSetSolveAtolRtol(tol, tol);
    active = true;
  }
  ~CovSolveTolGuard() {
    if (active) rxSetSolveAtolRtol(savAtol, savRtol);
  }
};

// [[Rcpp::export]]
NumericVector foceiSetup_(const RObject &obj,
                          const RObject &data,
                          NumericVector theta,
                          IntegerVector mixIdx,
                          Nullable<LogicalVector> thetaFixed = R_NilValue,
                          Nullable<LogicalVector> skipCov    = R_NilValue,
                          RObject rxInv                      = R_NilValue,
                          Nullable<NumericVector> lower      = R_NilValue,
                          Nullable<NumericVector> upper      = R_NilValue,
                          Nullable<NumericMatrix> etaMat     = R_NilValue,
                          Nullable<List> control             = R_NilValue){
  if (control.isNull()){
    stop("ODE options must be specified.");
  }
  List foceiO = as<List>(control);
  List rxControl = as<List>(foceiO["rxControl"]);
  // This fills in op_focei.neta
  List mvi;
  if (!Rf_isNull(obj)){
    if (!rxode2::rxDynLoad(obj)){
      stop("Cannot load rxode2 dlls for this model.");
    }
    mvi = rxode2::rxModelVars_(obj);
  }
  rxOptionsFreeFocei();
  op_focei.mvi = mvi;
  op_focei.mixIdxN = (unsigned int)mixIdx.size();
  op_focei.adjLik = as<bool>(foceiO["adjLik"]);
  op_focei.badSolveObjfAdj=fabs(as<double>(foceiO["badSolveObjfAdj"]));
  if (foceiO.containsElementNamed("est") && TYPEOF(foceiO["est"]) == STRSXP) {
    std::string estStr = as<std::string>(foceiO["est"]);
    op_focei.isSaem = (estStr == "saem");
    op_focei.isNlm = (estStr == "nlm" || estStr == "nlminb" || estStr == "bobyqa" ||
                      estStr == "newuoa" || estStr == "n1qn1" || estStr == "lbfgsb3c" ||
                      estStr == "optim" || estStr == "uobyqa" || estStr == "nls");
    op_focei.isImpmap = (estStr == "impmap" || estStr == "imp" || estStr == "qrpem");
    op_focei.isImp = (estStr == "imp");
    op_focei.isAdvi = (estStr == "advi");
    op_focei.isQrpem = (estStr == "qrpem");
  } else {
    op_focei.isSaem = false;
    op_focei.isNlm = false;
    op_focei.isImpmap = false;
    op_focei.isImp = false;
    op_focei.isAdvi = false;
    op_focei.isQrpem = false;
  }
  if (op_focei.isImpmap) {
    if (foceiO.containsElementNamed("isample")) op_focei.impIsample = as<int>(foceiO["isample"]);
    if (foceiO.containsElementNamed("gamma")) op_focei.impGamma = as<double>(foceiO["gamma"]);
    if (foceiO.containsElementNamed("nIter")) op_focei.impNiter = as<int>(foceiO["nIter"]);
    if (foceiO.containsElementNamed("iaccept")) op_focei.impIaccept = as<double>(foceiO["iaccept"]);
    if (foceiO.containsElementNamed("iscaleMin")) op_focei.impIscaleMin = as<double>(foceiO["iscaleMin"]);
    if (foceiO.containsElementNamed("iscaleMax")) op_focei.impIscaleMax = as<double>(foceiO["iscaleMax"]);
    if (foceiO.containsElementNamed("ctol") && !Rf_isNull(foceiO["ctol"]))
      op_focei.impCtol = as<double>(foceiO["ctol"]);
    if (foceiO.containsElementNamed("nConvWindow")) op_focei.impNconvWindow = as<int>(foceiO["nConvWindow"]);
    if (foceiO.containsElementNamed("impCov")) op_focei.impCov = as<bool>(foceiO["impCov"]);
    if (foceiO.containsElementNamed("qr")) op_focei.impQr = as<bool>(foceiO["qr"]);
    if (foceiO.containsElementNamed("qrShift")) op_focei.impQrShift = as<bool>(foceiO["qrShift"]);
    if (foceiO.containsElementNamed("qrRefresh")) op_focei.impQrRefresh = as<bool>(foceiO["qrRefresh"]);
    if (foceiO.containsElementNamed("sir")) op_focei.impSir = as<bool>(foceiO["sir"]);
    if (foceiO.containsElementNamed("sirSample")) op_focei.impSirSample = as<int>(foceiO["sirSample"]);
    if (foceiO.containsElementNamed("impSeed")) op_focei.impSeed = as<int>(foceiO["impSeed"]);
    if (foceiO.containsElementNamed("diagXform") && TYPEOF(foceiO["diagXform"]) == STRSXP)
      op_focei.impDiagXform = as<std::string>(foceiO["diagXform"]);
    if (foceiO.containsElementNamed("impMuThetaIdx"))
      op_focei.impMuThetaIdx = as<IntegerVector>(foceiO["impMuThetaIdx"]);
    if (foceiO.containsElementNamed("impMuEtaIdx"))
      op_focei.impMuEtaIdx = as<IntegerVector>(foceiO["impMuEtaIdx"]);
    if (foceiO.containsElementNamed("impThetaSensIdx"))
      op_focei.impThetaSensIdx = as<IntegerVector>(foceiO["impThetaSensIdx"]);
    if (foceiO.containsElementNamed("impOmegaFixedEta"))
      op_focei.impOmegaFixedEta = as<IntegerVector>(foceiO["impOmegaFixedEta"]);
  }
  // est="advi" reuses the theta-sensitivity model (impThetaSensIdx) for the outer
  // population gradient, but is not isImpmap; load the index here too.
  if (op_focei.isAdvi && foceiO.containsElementNamed("impThetaSensIdx")) {
    op_focei.impThetaSensIdx = as<IntegerVector>(foceiO["impThetaSensIdx"]);
  }

  op_focei.zeroGrad = false;
  op_focei.resetThetaCheckPer = as<double>(foceiO["resetThetaCheckPer"]);
  op_focei.printTop = as<int>(foceiO["printTop"]);
  op_focei.nF2 = as<int>(foceiO["nF"]);
  if (op_focei.nF2 == 0){
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
  }
  op_focei.didLikCalc = false;
  op_focei.maxOuterIterations = as<int>(foceiO["maxOuterIterations"]);
  op_focei.maxInnerIterations = as<int>(foceiO["maxInnerIterations"]);
  op_focei.mceta = as<int>(foceiO["mceta"]);
  op_focei.warm = foceiO.containsElementNamed("warm") ? as<int>(foceiO["warm"]) : 0;
  op_focei.maxOdeRecalc = as<int>(foceiO["maxOdeRecalc"]);
  op_focei.objfRecalN=0;
  op_focei.odeRecalcFactor = as<double>(foceiO["odeRecalcFactor"]);
  op_focei.indTolRelax = foceiO.containsElementNamed("indTolRelax") ? as<bool>(foceiO["indTolRelax"]) : true;
  op_focei.reducedTol = 0;
  op_focei.repeatGill=0;
  op_focei.repeatGillN=0;
  op_focei.repeatGillMax=as<int>(foceiO["repeatGillMax"]);
  op_focei.stickyRecalcN=as<int>(foceiO["stickyRecalcN"]);
  op_focei.neta = (unsigned int)as<int>(foceiO["neta"]);
  // Zero gradient reset setup: NA -> zeroGradReset=FALSE, zeroGradBobyqa=FALSE,
  // zeroGradBobyqaRun tracks the zeroGradBobyqa option (last reset calls bobyqa
  // instead); TRUE/FALSE -> zeroGradBobyqaRun matches zeroGradBobyqa directly;
  // zeroGradBobyqa=NA -> zeroGradBobyqa=true, zeroGradBobyqaRun=false.
  int curI = as<int>(foceiO["zeroGradFirstReset"]);
  if (curI == 1) {
    op_focei.zeroGradFirstReset = true;
    curI = as<int>(foceiO["zeroGradBobyqa"]);
    if (curI == 1) {
      op_focei.zeroGradBobyqaRun = op_focei.zeroGradBobyqa = true;
    } else if (curI == 0) {
      op_focei.zeroGradBobyqaRun = op_focei.zeroGradBobyqa = false;
    } else {
      op_focei.zeroGradBobyqaRun = false;
      op_focei.zeroGradBobyqa = true;
    }
  } else {
    op_focei.zeroGradFirstReset = false;
    if (curI == NA_INTEGER) {
      op_focei.zeroGradBobyqa = false;
      op_focei.zeroGradBobyqaRun = as<bool>(foceiO["zeroGradBobyqa"]);
    } else {
      curI = as<int>(foceiO["zeroGradBobyqa"]);
      if (curI == 1) {
        op_focei.zeroGradBobyqaRun = op_focei.zeroGradBobyqa = true;
      } else if (curI == 0)  {
        op_focei.zeroGradBobyqaRun = op_focei.zeroGradBobyqa = false;
      } else {
        op_focei.zeroGradBobyqaRun = false;
        op_focei.zeroGradBobyqa = true;
      }
    }
  }
  op_focei.zeroGradRunReset = as<bool>(foceiO["zeroGradRunReset"]);
  if (op_focei.neta != 0) {
    if (!rxode2::rxIs(rxInv, "rxSymInvCholEnv")){
      stop(_("Omega isn't in the proper format"));
    } else {
      _rxInv = as<List>(rxInv);
    }
  }
  op_focei.omegan = (unsigned int)as<int>(foceiO["nomega"]);
  op_focei.ntheta = (unsigned int)as<int>(foceiO["ntheta"]);
  // op_focei.ntheta = op_focei.ntheta;
  op_focei.nfixed = as<int>(foceiO["nfixed"]);
  int* tempMixIdx = NULL;
  if (op_focei.mixIdxN > 0) {
    tempMixIdx = R_Calloc(op_focei.mixIdxN, int);
    std::copy(mixIdx.begin(), mixIdx.end(), tempMixIdx);
    op_focei.mixIdx = tempMixIdx;
  }

  // Mu-group setup must happen before foceiSetupTheta_() below, which reads
  // muGroupThetaSkip while building fixedTrans/npars. Kept as its own R_Calloc
  // block since the gillRet combined buffer's size depends on npars, which
  // is only known after foceiSetupTheta_() runs.
  if (op_focei.muGroupTheta != NULL) R_Free(op_focei.muGroupTheta);
  op_focei.muGroupTheta = NULL; op_focei.muGroupEta = NULL;
  op_focei.muGroupCovStart = NULL; op_focei.muGroupCovCount = NULL;
  op_focei.muGroupCovTheta = NULL; op_focei.muGroupCovUserFixed = NULL;
  op_focei.muGroupCovBounded = NULL;
  op_focei.muGroupThetaSkip = NULL;
  op_focei.muModel = foceiO.containsElementNamed("foceiMuModel") ?
    as<int>(foceiO["foceiMuModel"]) : 0;
  IntegerVector muGroupTheta, muGroupEta, muGroupCovStart, muGroupCovCount,
    muGroupCovTheta, muGroupCovUserFixed, muGroupCovBounded;
  if (op_focei.muModel != 0 && foceiO.containsElementNamed("foceiMuGroupTheta")) {
    muGroupTheta = as<IntegerVector>(foceiO["foceiMuGroupTheta"]);
    muGroupEta = as<IntegerVector>(foceiO["foceiMuGroupEta"]);
    muGroupCovStart = as<IntegerVector>(foceiO["foceiMuGroupCovStart"]);
    muGroupCovCount = as<IntegerVector>(foceiO["foceiMuGroupCovCount"]);
    muGroupCovTheta = as<IntegerVector>(foceiO["foceiMuGroupCovTheta"]);
    muGroupCovUserFixed = as<IntegerVector>(foceiO["foceiMuGroupCovUserFixed"]);
    muGroupCovBounded = foceiO.containsElementNamed("foceiMuGroupCovBounded") ?
      as<IntegerVector>(foceiO["foceiMuGroupCovBounded"]) : IntegerVector(muGroupCovTheta.size());
  }
  op_focei.muGroupN = (unsigned int)muGroupTheta.size();
  op_focei.muGroupCovN = (unsigned int)muGroupCovTheta.size();
  if (op_focei.muGroupN > 0) {
    size_t totalInt = 4*(size_t)op_focei.muGroupN + 3*(size_t)op_focei.muGroupCovN +
      (size_t)op_focei.ntheta;
    int *buf = R_Calloc(totalInt, int);
    op_focei.muGroupTheta = buf;
    op_focei.muGroupEta = buf + op_focei.muGroupN;
    op_focei.muGroupCovStart = op_focei.muGroupEta + op_focei.muGroupN;
    op_focei.muGroupCovCount = op_focei.muGroupCovStart + op_focei.muGroupN;
    op_focei.muGroupCovTheta = op_focei.muGroupCovCount + op_focei.muGroupN;
    op_focei.muGroupCovUserFixed = op_focei.muGroupCovTheta + op_focei.muGroupCovN;
    op_focei.muGroupCovBounded = op_focei.muGroupCovUserFixed + op_focei.muGroupCovN;
    op_focei.muGroupThetaSkip = op_focei.muGroupCovBounded + op_focei.muGroupCovN;
    std::copy(muGroupTheta.begin(), muGroupTheta.end(), op_focei.muGroupTheta);
    std::copy(muGroupEta.begin(), muGroupEta.end(), op_focei.muGroupEta);
    std::copy(muGroupCovStart.begin(), muGroupCovStart.end(), op_focei.muGroupCovStart);
    std::copy(muGroupCovCount.begin(), muGroupCovCount.end(), op_focei.muGroupCovCount);
    if (op_focei.muGroupCovN > 0) {
      std::copy(muGroupCovTheta.begin(), muGroupCovTheta.end(), op_focei.muGroupCovTheta);
      std::copy(muGroupCovUserFixed.begin(), muGroupCovUserFixed.end(), op_focei.muGroupCovUserFixed);
      std::copy(muGroupCovBounded.begin(), muGroupCovBounded.end(), op_focei.muGroupCovBounded);
    }
    for (unsigned int g = 0; g < op_focei.muGroupN; g++) {
      op_focei.muGroupThetaSkip[op_focei.muGroupTheta[g]] = 1;
    }
    for (unsigned int c = 0; c < op_focei.muGroupCovN; c++) {
      // Bounded covariates stay OUT of muGroupThetaSkip: they remain ordinary
      // outer-optimized free parameters (updateMuGroups() still excludes them
      // from the regression design matrix).
      if (op_focei.muGroupCovBounded[c] == 0) {
        op_focei.muGroupThetaSkip[op_focei.muGroupCovTheta[c]] = 1;
      }
    }
    if (foceiO.containsElementNamed("foceiMuGroupCovData")) {
      NumericMatrix covData = as<NumericMatrix>(foceiO["foceiMuGroupCovData"]);
      op_focei.muGroupCovData = as<arma::mat>(covData);
    }
    op_focei.muGroupTol = foceiO.containsElementNamed("foceiMuGroupTol") ?
      as<double>(foceiO["foceiMuGroupTol"]) : 1e-3;
    op_focei.muGroupMaxCycles = foceiO.containsElementNamed("foceiMuGroupMaxCycles") ?
      as<int>(foceiO["foceiMuGroupMaxCycles"]) : 10;
  }

  if (op_focei.maxOuterIterations <= 0){
    // No scaling.
    foceiSetupTheta_(mvi, theta, thetaFixed, 0.0, !Rf_isNull(obj));
    op_focei.scaleObjective=0;
  } else {
    foceiSetupTheta_(mvi, theta, thetaFixed, as<double>(foceiO["scaleTo"]), !Rf_isNull(obj));
    op_focei.scaleObjectiveTo=as<double>(foceiO["scaleObjective"]);
    if (op_focei.scaleObjectiveTo <= 0){
      op_focei.scaleObjective=0;
    } else {
      op_focei.scaleObjective=1;
    }
  }
  if (tempMixIdx != NULL) {
    R_Free(tempMixIdx);
    op_focei.mixIdx = NULL;
  }
  // First see if etaMat is null.
  NumericMatrix etaMat0;
  unsigned int nsub=0;
  // Find the number of IDs to create an etaMat
  List df = as<List>(data);
  CharacterVector dfN = df.names();
  int idn = -1;
  std::string cur;
  for (unsigned int j = dfN.size(); j--;){
    cur = as<std::string>(dfN[j]);
    if (cur == "ID" || cur == "id" || cur == "Id" || cur == "iD"){
      idn=j;
      break;
    }
  }
  if  (idn == -1){
    stop("Can't find ID in dataset.");
  }
  IntegerVector ids = as<IntegerVector>(df[idn]);
  int last = ids[ids.size()-1]-1;
  for (unsigned int j = ids.size(); j--;){
    if (last != ids[j]){
      last = ids[j];
      nsub++;
    }
  }
  int expected_nsub = nsub * (op_focei.mixIdxN + 1);
  if (op_focei.neta > 0) {
    if (!etaMat.isNull()){
      if (TYPEOF(wrap(etaMat)) != REALSXP) {
        // Rcpp::print("etaMat is not a numeric matrix");
        Rcpp::print(etaMat);
        stop("etaMat must be a numeric matrix");
      }
      NumericMatrix etaMat1 = NumericMatrix(etaMat);
      if (etaMat1.nrow() == expected_nsub) {
        // ok
      } else if (etaMat1.nrow() != (int)nsub){
        print(etaMat1);
        stop("The etaMat must have the same number of ETAs (rows) as subjects.");
      }
      if (etaMat1.ncol() != op_focei.neta){
        print(etaMat1);
        stop("The etaMat must have the same number of ETAs (cols) as the model.");
      }
      etaMat0 = NumericMatrix(expected_nsub, op_focei.neta);
      if (etaMat1.nrow() == expected_nsub) {
        std::copy(etaMat1.begin(),etaMat1.end(),etaMat0.begin());
      } else {
        int nMix = op_focei.mixIdxN + 1;
        for (int m = 0; m < nMix; m++) {
          for (int r = 0; r < (int)nsub; r++) {
            for (int c = 0; c < (int)op_focei.neta; c++) {
              etaMat0(m * nsub + r, c) = etaMat1(r, c);
            }
          }
        }
      }
    } else {
      etaMat0 = NumericMatrix(expected_nsub, op_focei.neta);
      std::fill(etaMat0.begin(), etaMat0.end(), 0.0);
    }
  }
  List params(theta.size()+op_focei.neta);
  CharacterVector paramsNames(theta.size()+op_focei.neta);
  unsigned int j;
  for (j = theta.size();j--;){
    params[j] = NumericVector(expected_nsub,theta[j]);
    if (theta.hasAttribute("names")){
      paramsNames[j] = (as<CharacterVector>(theta.names()))[j];
    } else {
      paramsNames[j] = "THETA[" + std::to_string(j + 1) + "]";
    }
  }
  if (op_focei.neta > 0) {
    bool hasDimn = etaMat0.hasAttribute("dimnames");
    CharacterVector dims;
    if (hasDimn){
      List diml = etaMat0.attr("dimnames");
      if (!Rf_isNull(as<RObject>(diml[1]))){
        dims = as<CharacterVector>(diml[1]);
      } else {
        hasDimn=false;
      }
    }
    for (j=op_focei.neta; j--;){
      params[j+theta.size()]= etaMat0(_, j);
      if (hasDimn){
        paramsNames[j+theta.size()] = dims[j];
      } else {
        paramsNames[j+theta.size()] = "ETA[" + std::to_string(j + 1) + "]";
      }
    }
  }
  params.names() = paramsNames;
  params.attr("class") = "data.frame";
  params.attr("row.names") = IntegerVector::create(NA_INTEGER,-expected_nsub);
  // Now pre-fill parameters.
  if (!Rf_isNull(obj)) {
    // Don't force cores=1: the user-requested core count must propagate to
    // rxSolve_ so per-thread buffers are sized correctly (parallelization
    // itself is managed in innerOpt()).
    // stash the solve args for covType="analytic" and fast (both replace the fit's
    // global solve with augmented-sensitivity solves and need restoreFitSolve_);
    // avoids retaining the model/dataset for the session on every plain FD fit
    bool _needSolveArgs =
      (foceiO.containsElementNamed("covType") &&
       as<std::string>(foceiO["covType"]) == "analytic") ||
      (foceiO.containsElementNamed("fast") && as<int>(foceiO["fast"]) != 0);
    if (_needSolveArgs) {
      storeCovSolveArgs_(obj, rxControl, params, data);
    }
    // est="impmap": size the shared solve pool for the theta-sensitivity model (the
    // largest structure) when set; the inner MAP then uses ind->neqOverride.
    RObject _poolObj = (_impPoolModel != R_NilValue) ? RObject(_impPoolModel) : obj;
    rxode2::rxSolve_(_poolObj, rxControl,
                     R_NilValue,//const Nullable<CharacterVector> &specParams =
                     R_NilValue,//const Nullable<List> &extraArgs =
                     as<RObject>(params),//const RObject &params =
                     data,//const RObject &events =
                     R_NilValue, // inits
                     1);//const int setupOnly = 0
    rx = getRxSolve_();
    if (op_focei.neta == 0) foceiSetupNoEta_();
    else foceiSetupEta_(etaMat0);
  }
  op_focei.epsilon=as<double>(foceiO["epsilon"]);
  op_focei.nsim=as<int>(foceiO["n1qn1nsim"]);
  op_focei.imp=0;
  op_focei.resetThetaSize = std::numeric_limits<double>::infinity();
  // op_focei.printInner=as<int>(foceiO["printInner"]);
  // if (op_focei.printInner < 0) op_focei.printInner = -op_focei.printInner;
  unsigned int totN=op_focei.ntheta + op_focei.omegan;
  NumericVector cEps=foceiO["derivEps"];
  if (cEps.size() != 2){
    stop("derivEps must be 2 elements for determining forward difference step size.");
  }
  NumericVector covDerivEps=foceiO["centralDerivEps"];
  if (covDerivEps.size() != 2){
    stop("centralDerivEps must be 2 elements for determining central difference step size.");
  }
  op_focei.derivMethod = as<int>(foceiO["derivMethod"]);
  if (op_focei.derivMethod == 3){
    op_focei.derivMethod=0;
    op_focei.derivSwitchTol=as<double>(foceiO["derivSwitchTol"]);
    op_focei.derivMethodSwitch=1;
  }

  IntegerVector muRef;
  if (foceiO.containsElementNamed("foceiMuRef")){
    if (rxode2::rxIs(foceiO["foceiMuRef"], "integer")) {
      muRef = as<IntegerVector>(foceiO["foceiMuRef"]);
    }
  }
  if (muRef.size() == 0){
    op_focei.resetThetaSize = std::numeric_limits<double>::infinity();
    op_focei.muRefN=0;
  } else{
    op_focei.muRefN=(unsigned int)muRef.size();
  }

  // Per-eta opt-out from eta-drift zero-reset / mu-ref theta soft-shift,
  // same length/ordering as foceiMuRef; unused (stays zeroed) when muModel="none".
  IntegerVector muCovEta;
  bool haveMuCovEta = false;
  if (foceiO.containsElementNamed("foceiMuCovEta")){
    if (rxode2::rxIs(foceiO["foceiMuCovEta"], "integer")) {
      muCovEta = as<IntegerVector>(foceiO["foceiMuCovEta"]);
      haveMuCovEta = ((unsigned int)muCovEta.size() == op_focei.muRefN);
    }
  }

  IntegerVector skipCov1;
  if (skipCov.isNull()){
    op_focei.skipCovN = 0;
  } else {
    skipCov1 = as<IntegerVector>(skipCov);
    op_focei.skipCovN = (unsigned int)skipCov1.size();
  }

  if (op_focei.gillRet != NULL) R_Free(op_focei.gillRet);
  op_focei.gillRet = R_Calloc(2*totN + op_focei.npars +
                              op_focei.muRefN + op_focei.muRefN + op_focei.skipCovN +
                              op_focei.mixIdxN +
                              (size_t)getRxNsub(rx),
                              int);
  op_focei.gillRetC= op_focei.gillRet + totN;
  op_focei.nbd     = op_focei.gillRetC + totN;//[op_focei.npars]

  op_focei.muRef   = op_focei.nbd + op_focei.npars; //[op_focei.muRefN]
  if (op_focei.muRefN) {
    std::copy(&op_focei.muRef[0], &op_focei.muRef[0]+op_focei.muRefN, muRef.begin());
  }

  op_focei.muRefEtaCovSkipReset = op_focei.muRef + op_focei.muRefN; //[op_focei.muRefN]
  if (op_focei.muRefN && haveMuCovEta) {
    std::copy(muCovEta.begin(), muCovEta.end(), &op_focei.muRefEtaCovSkipReset[0]);
  }

  op_focei.mixIdx = op_focei.muRefEtaCovSkipReset + op_focei.muRefN; // [op_focei.mixIdxN  + getRxNsub(rx)]
  if (op_focei.mixIdxN) {
    std::copy(mixIdx.begin(), mixIdx.end(), op_focei.mixIdx);
    // Re-run mixTrans now that mixIdx is filled (earlier call only set -1 placeholders).
    for (unsigned int _k = op_focei.npars; _k--;) {
      foceiSetupMixTrans((int)_k, op_focei.fixedTrans[_k]);
    }
  }

  op_focei.skipCov   = op_focei.mixIdx + op_focei.mixIdxN + getRxNsub(rx); //[op_focei.skipCovN]
  if (op_focei.skipCovN) {
    std::copy(skipCov1.begin(),skipCov1.end(),op_focei.skipCov); //
  }

  if (op_focei.gillDf != NULL) R_Free(op_focei.gillDf);

  op_focei.gillDf = R_Calloc(7*totN + 2*op_focei.npars +
                             (op_focei.mixIdxN + 1)*((size_t)getRxNsub(rx)+1) +
                             op_focei.mixIdxN*((size_t)getRxNsub(rx)+1) +
                             (size_t)getRxNsub(rx), double);
  op_focei.mixProb = op_focei.gillDf+totN; // [op_focei.mixIdN+1 + (op_focei.mixIdN+1)*getRxNsub(rx)]
  op_focei.mixProbGrad = op_focei.mixProb + (op_focei.mixIdxN+1)*(getRxNsub(rx)+1);
  op_focei.gillDf2 = op_focei.mixProbGrad + (op_focei.mixIdxN)*(getRxNsub(rx)+1);
  // Patch per-individual pointers left NULL during earlier inds_focei init.
  if (op_focei.mixIdxN != 0) {
    for (size_t _mi = getRxNsubAndMix(rx); _mi--;) {
      focei_ind *_fI = &(inds_focei[_mi]);
      _fI->mixest     = op_focei.mixIdx + op_focei.mixIdxN + getRxId(_mi);
      _fI->mixProb    = op_focei.mixProb + (getRxId(_mi) + 1)*(op_focei.mixIdxN + 1);
      _fI->mixProbGrad= op_focei.mixProbGrad + (getRxId(_mi) + 1)*(op_focei.mixIdxN);
    }
  }
  op_focei.gillErr = op_focei.gillDf2+totN;
  op_focei.rEps=op_focei.gillErr + totN;
  op_focei.aEps = op_focei.rEps + totN;
  op_focei.rEpsC = op_focei.aEps + totN;
  op_focei.aEpsC = op_focei.rEpsC + totN;
  op_focei.lower = op_focei.aEpsC + totN;
  op_focei.upper = op_focei.lower +op_focei.npars;
  op_focei.likSav      = op_focei.upper + op_focei.npars;//[getRxNsub(rx)]

  if (op_focei.mixIdxN > 0) {
    for (unsigned int i = 0; i < getRxNsubAndMix(rx); i++) {
      focei_ind *fInd = &(inds_focei[i]);
      fInd->mixest = op_focei.mixIdx + op_focei.mixIdxN + getRxId(i);
      fInd->mixProb = op_focei.mixProb + (getRxId(i) + 1)*(op_focei.mixIdxN + 1);
      fInd->mixProbGrad = op_focei.mixProbGrad + (getRxId(i) + 1)*(op_focei.mixIdxN);
    }
  }

  if (op_focei.derivMethod){
    std::fill_n(&op_focei.rEps[0], totN, std::fabs(cEps[0])/2.0);
    std::fill_n(&op_focei.aEps[0], totN, std::fabs(cEps[1])/2.0);
    std::fill_n(&op_focei.rEpsC[0], totN, std::fabs(covDerivEps[0])/2.0);
    std::fill_n(&op_focei.aEpsC[0], totN, std::fabs(covDerivEps[1])/2.0);
  } else {
    std::fill_n(&op_focei.rEps[0], totN, std::fabs(cEps[0]));
    std::fill_n(&op_focei.aEps[0], totN, std::fabs(cEps[1]));
    std::fill_n(&op_focei.rEpsC[0], totN, std::fabs(covDerivEps[0]));
    std::fill_n(&op_focei.aEpsC[0], totN, std::fabs(covDerivEps[1]));
  }
  NumericVector lowerIn(totN);
  NumericVector upperIn(totN);
  if (lower.isNull()){
    std::fill_n(lowerIn.begin(), totN, R_NegInf);
  } else {
    NumericVector lower1=as<NumericVector>(lower);
    if (lower1.size() == 1){
      std::fill_n(lowerIn.begin(), totN, lower1[0]);
    } else if (lower1.size() < totN){
      std::copy(lower1.begin(), lower1.end(), lowerIn.begin());
      std::fill_n(lowerIn.begin()+lower1.size(), totN - lower1.size(), R_NegInf);
    } else if (lower1.size() > totN){
      warning(_("lower bound is larger than the number of parameters being estimated"));
      std::copy(lower1.begin(), lower1.begin()+totN, lowerIn.begin());
    } else {
      lowerIn = lower1;
    }
  }
  if (upper.isNull()){
    std::fill_n(upperIn.begin(), totN, std::numeric_limits<double>::infinity());
  } else {
    NumericVector upper1=as<NumericVector>(upper);
    if (upper1.size() == 1){
      std::fill_n(upperIn.begin(), totN, upper1[0]);
    } else if (upper1.size() < totN){
      std::copy(upper1.begin(), upper1.end(), upperIn.begin());
      std::fill_n(upperIn.begin()+upper1.size(), totN - upper1.size(), std::numeric_limits<double>::infinity());
    } else if (upper1.size() > totN){
      warning(_("upper bound is larger than the number of parameters being estimated"));
      std::copy(upper1.begin(), upper1.begin()+totN, upperIn.begin());
    } else {
      upperIn = upper1;
    }
  }
  op_focei.lowerIn =lowerIn;
  op_focei.upperIn =upperIn;

  std::fill_n(&op_focei.nbd[0], op_focei.npars, 0);

  NumericVector ret(op_focei.npars, op_focei.scaleTo);
  op_focei.calcGrad=0;

  // Outer options
  op_focei.outerOpt =as<int>(foceiO["outerOpt"]);
  // lbfgsb options
  op_focei.factr    = as<double>(foceiO["lbfgsFactr"]);
  op_focei.pgtol    = as<double>(foceiO["lbfgsPgtol"]);
  op_focei.lmm      = as<int>(foceiO["lbfgsLmm"]);
  op_focei.covDerivMethod = as<int>(foceiO["covDerivMethod"]);
  int type = TYPEOF(foceiO["covMethod"]);
  if (type == INTSXP || type == REALSXP) {
    op_focei.covMethod = as<int>(foceiO["covMethod"]);
  } else {
    Rcpp::stop("'covMethod' needs to be an integer in lower level inner");
  }
  op_focei.eigen = as<int>(foceiO["eigen"]);
  op_focei.ci=as<double>(foceiO["ci"]);
  op_focei.sigdig=as<double>(foceiO["sigdigTable"]);
  op_focei.boundTol=as<double>(foceiO["boundTol"]);
  scaleApplyIterPrintControl(&op_focei.scale,
                             as<List>(foceiO["iterPrintControl"]));
  op_focei.noabort=as<int>(foceiO["noAbort"]);
  op_focei.interaction=as<int>(foceiO["interaction"]);
  op_focei.foceType=as<int>(foceiO["foceType"]);
  op_focei.fast=foceiO.containsElementNamed("fast") ? as<int>(foceiO["fast"]) : 0;
  op_focei.cholSEtol=as<double>(foceiO["cholSEtol"]);
  op_focei.hessEps=as<double>(foceiO["hessEps"]);
  op_focei.hessEpsLlik=as<double>(foceiO["hessEpsLlik"]);
  op_focei.hessEpsInner=as<double>(rxControl[Rxc_atolSens]);
  op_focei.shi21maxOuter = as<int>(foceiO["shi21maxOuter"]);
  op_focei.shi21maxInner = as<int>(foceiO["shi21maxInner"]);
  op_focei.shi21maxInnerCov = as<int>(foceiO["shi21maxInnerCov"]);
  op_focei.optimHessCovType=as<int>(foceiO["optimHessCovType"]);
  op_focei.shi21maxFD = as<int>(foceiO["shi21maxFD"]);
  op_focei.optimHessType=as<int>(foceiO["optimHessType"]);
  // censOption: 0 "gauss" (historic uncensored Gauss-Newton, default) / 1 "laplace" (exact
  // censored 2nd deriv).  Tolerate an older control missing the field (-> gauss).
  op_focei.censOption = foceiO.containsElementNamed("censOption") ? as<int>(foceiO["censOption"]) : 0;
  op_focei.cholAccept=as<double>(foceiO["cholAccept"]);
  op_focei.resetEtaSize=as<double>(foceiO["resetEtaSize"]);
  op_focei.resetThetaSize=as<double>(foceiO["resetThetaSize"]);
  op_focei.resetThetaFinalSize = as<double>(foceiO["resetThetaFinalSize"]);
  op_focei.needOptimHess = as<bool>(foceiO["needOptimHess"]);

  op_focei.cholSEOpt=as<double>(foceiO["cholSEOpt"]);
  op_focei.cholSECov=as<double>(foceiO["cholSECov"]);
  op_focei.fo=as<int>(foceiO["fo"]);
  if (op_focei.neta == 0) op_focei.fo = 0;
  if (op_focei.fo) op_focei.maxInnerIterations=0;

  op_focei.covTryHarder=as<int>(foceiO["covTryHarder"]);
  op_focei.resetHessianAndEta=as<int>(foceiO["resetHessianAndEta"]);
  op_focei.gillK = as<int>(foceiO["gillK"]);
  op_focei.gillKcov = as<int>(foceiO["gillKcov"]);
  op_focei.gillKcovLlik = as<int>(foceiO["gillKcovLlik"]);
  op_focei.gillStep    = fabs(as<double>(foceiO["gillStep"]));
  op_focei.gillStepCov = fabs(as<double>(foceiO["gillStepCov"]));
  op_focei.gillStepCovLlik=fabs(as<double>(foceiO["gillStepCovLlik"]));
  op_focei.gillFtol    = fabs(as<double>(foceiO["gillFtol"]));
  op_focei.gillFtolCov = fabs(as<double>(foceiO["gillFtolCov"]));
  op_focei.gillFtolCovLlik = fabs(as<double>(foceiO["gillFtolCovLlik"]));
  op_focei.didGill = 0;
  op_focei.gillRtol = as<double>(foceiO["gillRtol"]);
  op_focei.scaleType = as<int>(foceiO["scaleType"]);
  op_focei.normType = as<int>(foceiO["normType"]);
  op_focei.scaleC0=as<double>(foceiO["scaleC0"]);
  op_focei.scaleCmin=as<double>(foceiO["scaleCmin"]);
  op_focei.scaleCmax=as<double>(foceiO["scaleCmax"]);
  op_focei.abstol=as<double>(foceiO["abstol"]);
  op_focei.reltol=as<double>(foceiO["reltol"]);
  op_focei.smatNorm=as<int>(foceiO["smatNorm"]);
  op_focei.smatNormLlik=as<int>(foceiO["smatNormLlik"]);
  op_focei.rmatNorm=as<int>(foceiO["rmatNorm"]);
  op_focei.rmatNormLlik=as<int>(foceiO["rmatNormLlik"]);
  op_focei.covGillF=as<int>(foceiO["covGillF"]);
  op_focei.covFull = foceiO.containsElementNamed("covFull") ? as<int>(foceiO["covFull"]) : 0;
  op_focei.covFdDirect = 0;
  op_focei.optGillF=as<int>(foceiO["optGillF"]);
  op_focei.covSmall = as<double>(foceiO["covSmall"]);
  op_focei.gradTrim = as<double>(foceiO["gradTrim"]);
  op_focei.gradCalcCentralLarge = as<double>(foceiO["gradCalcCentralLarge"]);
  op_focei.gradCalcCentralSmall = as<double>(foceiO["gradCalcCentralSmall"]);
  op_focei.etaNudge = as<double>(foceiO["etaNudge"]);
  op_focei.etaNudge2 = as<double>(foceiO["etaNudge2"]);
  op_focei.eventType = as<int>(foceiO["eventType"]);
  op_focei.predNeq = as<int>(foceiO["predNeq"]);
  op_focei.gradProgressOfvTime = as<double>(foceiO["gradProgressOfvTime"]);
  op_focei.fallbackFD = as<int>(foceiO["fallbackFD"]);
  op_focei.smatPer = as<double>(foceiO["smatPer"]);
  op_focei.initObj=0;
  op_focei.lastOfv=std::numeric_limits<double>::max();
  for (unsigned int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    ret[k] = op_focei.fullTheta[j];
    if (R_FINITE(lowerIn[j])){
      op_focei.lower[k] = lowerIn[j];
      op_focei.lower[k] += 2*(op_focei.lower[k]*op_focei.rEps[k] + op_focei.aEps[k]);
      // lower bound only = 1
      op_focei.nbd[k]=1;
    } else {
      op_focei.lower[k] = R_NegInf;//std::numeric_limits<double>::lowest();
    }
    if (R_FINITE(upperIn[j])){
      op_focei.upper[k] = upperIn[j];
      op_focei.upper[k] -= 2*(op_focei.upper[k]*op_focei.rEps[k] - op_focei.aEps[k]);
      // Upper bound only = 3
      // Upper and lower bound = 2
      op_focei.nbd[k]= 3 - op_focei.nbd[j];
    } else {
      op_focei.upper[k] = std::numeric_limits<double>::infinity();//std::numeric_limits<double>::max();
    }
  }
  double mn = (op_focei.npars > 0) ? op_focei.initPar[op_focei.npars-1] : 0.0;
  double mx = mn, mean=0, oN=0, oM=0, s=0;
  double len=0;
  unsigned int k;
  if (op_focei.nF2 > 0 && foceiO.containsElementNamed("c1") && foceiO.containsElementNamed("c2")){
    op_focei.c1 = foceiO["c1"];
    op_focei.c2 = foceiO["c2"];
  } else {
    switch (op_focei.normType){
    case 1:
      // OptdesX
      // http://apmonitor.com/me575/uploads/Main/optimization_book.pdf
      for (k = op_focei.npars-1; k--;){
        mn = min2(op_focei.initPar[k],mn);
        mx = max2(op_focei.initPar[k],mx);
      }
      if (mx == mn) {
        warning(_("all parameters are the same value, switch to length normType"));
        for (unsigned int k = op_focei.npars-1; k--;){
          len += op_focei.initPar[k]*op_focei.initPar[k];
        }
        op_focei.c1 = 0;
        op_focei.c2 = _safe_sqrt(len);
        op_focei.normType = 5;
      } else {
        op_focei.c1 = (mx+mn)/2;
        op_focei.c2 = (mx-mn)/2;
      }
      break;
    case 2: // Rescaling (min-max normalization)
      for (k = op_focei.npars-1; k--;){
        mn = min2(op_focei.initPar[k],mn);
        mx = max2(op_focei.initPar[k],mx);
      }
      if (mx == mn) {
        warning(_("all parameters are the same value, switch to length normType"));
        for (unsigned int k = op_focei.npars-1; k--;){
          len += op_focei.initPar[k]*op_focei.initPar[k];
        }
        op_focei.c1 = 0;
        op_focei.c2 = _safe_sqrt(len);
        op_focei.normType = 5;
      } else {
        op_focei.c1 = mn;
        op_focei.c2 = (mx-mn);
      }
      break;
    case 3: // Mean normalization
      for (k = op_focei.npars-1; k--;){
        mn = min2(op_focei.initPar[k],mn);
        mx = max2(op_focei.initPar[k],mx);
        oN++;
        mean += (op_focei.initPar[k]-mean)/oN;
      }
      if (mx == mn) {
        warning(_("all parameters are the same value, switch to length normType"));
        for (unsigned int k = op_focei.npars-1; k--;){
          len += op_focei.initPar[k]*op_focei.initPar[k];
        }
        op_focei.c1 = 0;
        op_focei.c2 = _safe_sqrt(len);
        op_focei.normType = 5;
      } else {
        op_focei.c1 = mean;
        op_focei.c2 = (mx-mn);
      }
      break;
    case 4: // Standardization
      for (k = op_focei.npars-1; k--;){
        mn = min2(op_focei.initPar[k],mn);
        mx = max2(op_focei.initPar[k],mx);
        oM= mean;
        oN++;
        mean += (op_focei.initPar[k]-mean)/oN;
        s += (op_focei.initPar[k]-mean)*(op_focei.initPar[k]-oM);
      }
      if (mx == mn) {
        warning("all parameters are the same value, switch to length norm type");
        for (unsigned int k = op_focei.npars-1; k--;){
          len += op_focei.initPar[k]*op_focei.initPar[k];
        }
        op_focei.c1 = 0;
        op_focei.c2 = _safe_sqrt(len);
        op_focei.normType = 5;
      } else {
        op_focei.c1 = mean;
        op_focei.c2 = _safe_sqrt(s/(oN-1));
      }
      break;
    case 5: // Normalize to length.
      for (unsigned int k = op_focei.npars-1; k--;){
        len += op_focei.initPar[k]*op_focei.initPar[k];
      }
      op_focei.c1 = 0;
      op_focei.c2 = _safe_sqrt(len);
      break;
    case 6:
      // No Normalization
      op_focei.c1 = 0;
      op_focei.c2 = 1;
      break;
    default:
      stop("unrecognized normalization (normType=%d)",op_focei.normType);
    }
  }
  return ret;
}

LogicalVector nlmixr2EnvSetup(Environment e, double fmin){
  if (e.exists("theta") && rxode2::rxIs(e["theta"], "data.frame") &&
      e.exists("omega") && e.exists("etaObf")) {
    int nobs2=0;
    if (e.exists("nobs2")) {
      nobs2=as<int>(e["nobs2"]);
    } else {
      nobs2 = getRxNobs2(rx);
    }
    if (op_focei.scaleObjective) {
      fmin = fmin * op_focei.initObjective / op_focei.scaleObjectiveTo;
    }

    bool doAdj = false;
    bool doObf = false;
    if (!e.exists("objective")){
      if (op_focei.adjLik){
        doAdj = true;
      } else {
        fmin -= 0.5*(nobs2)*log(M_2PI);
        doAdj=true;
      }
      e["objective"] = fmin;
    } else {
      fmin = as<double>(e["objective"]);
      if (op_focei.adjLik){
        doObf=true;
      }
    }
    e["OBJF"] = fmin;
    e["objf"] = fmin;
    NumericVector logLik(1);
    double adj= 0;
    if (doAdj){
      adj=nobs2*log(2*M_PI)/2;
    }
    e["adj"]=adj;
    logLik[0]=-fmin/2-adj;
    logLik.attr("df") = op_focei.npars;
    if (e.exists("nobs")){
      logLik.attr("nobs") = e["nobs"];
      e["BIC"] = fmin+2*adj + log(as<double>(e["nobs"]))*op_focei.npars;
    } else {
      logLik.attr("nobs") = nobs2;
      e["BIC"] = fmin + 2*adj + log((double)nobs2)*op_focei.npars;
      e["nobs"] = getRxNobs(rx);
    }
    if (!e.exists("nsub")) {
      e["nsub"] = getRxNsub(rx);
    }
    logLik.attr("class") = "logLik";
    e["logLik"] = logLik;

    e["AIC"] = fmin+2*adj+2*op_focei.npars;
    if (doObf){
      // -2 * object$logLik - object$dim$N * log(2 * pi)
      adj = -2*as<double>(logLik) - (nobs2)*log(M_2PI);
      e["OBJF"] = adj;
      e["objf"] = adj;
      e["objective"] = adj;
    }
    return true;
  } else {
    stop("Not Setup right\nNeeds: theta = data.frame for theta information \n omega= matrix\n etaObf=eta objective function");
    return false;
  }
}

void foceiOuterFinal(double *x, Environment e){
  // reset the optimal step size to have reproducible likelihoods with
  // generalized log-likelihood focei
  std::fill_n(op_focei.getahh, op_focei.gEtaGTransN, 0.0);
  // for (int i = op_focei.gEtaGTransN; i--;){
  //   op_focei.getahh[i] = -op_focei.getahh[i];
  // }
  // This will give reproducible likelihoods with fd events
  std::fill_n(op_focei.getahf, op_focei.gEtaGTransN, 0.0);
  std::fill_n(op_focei.getahr, op_focei.gEtaGTransN, 0.0);
  op_focei.optimHessType = op_focei.optimHessCovType;
  op_focei.shi21maxInner = op_focei.shi21maxInnerCov;
  _finalObfCalc = true;
  double fmin = foceiOfv0(x);
  _finalObfCalc = false;
  NumericVector theta(op_focei.ntheta);
  std::copy(&op_focei.fullTheta[0],  &op_focei.fullTheta[0] + op_focei.ntheta,
            theta.begin());

  NumericVector fullTheta(op_focei.ntheta+op_focei.omegan);
  std::copy(&op_focei.fullTheta[0],  &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
            fullTheta.begin());
  LogicalVector thetaFixed(op_focei.ntheta);
  std::fill_n(thetaFixed.begin(),op_focei.ntheta, true);
  int j;
  for (int k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    if (j < thetaFixed.size()) thetaFixed[j]=false;
  }
  // std::copy(&op_focei.thetaFixed[0],  &op_focei.thetaFixed[0] + op_focei.ntheta,
  //           thetaFixed.begin());
  NumericVector lowerIn(op_focei.ntheta);
  NumericVector upperIn(op_focei.ntheta);
  std::copy(&op_focei.lowerIn[0],  &op_focei.lowerIn[0] + op_focei.ntheta,
            lowerIn.begin());
  std::copy(&op_focei.upperIn[0],  &op_focei.upperIn[0] + op_focei.ntheta,
            upperIn.begin());
  e["theta"] = DataFrame::create(_["lower"]=lowerIn, _["theta"]=theta, _["upper"]=upperIn,
                                 _["fixed"]=thetaFixed);
  e["fullTheta"] = fullTheta;
  if (op_focei.neta == 0){
    e["omega"] = R_NilValue;
    e["etaObf"] = R_NilValue;
  } else {
    e["omega"] = getOmega();
    if (op_focei.mixIdxN != 0) {
      e["etaObfFull"] = foceiEtas(e);
      e["etaObf"] = foceiEtas(e, true);
    } else {
      e["etaObf"] = foceiEtas(e);
    }
    // create phi object for standard errors
    foceiPhi(e);
  }
  nlmixr2EnvSetup(e, fmin);
}

static inline void foceiPrintLine(int ncol){
  RSprintf("|-----+---------------+");
  for (int i = 0; i < ncol; i++){
    if (i == ncol-1)
      RSprintf("-----------|");
    else
      RSprintf("-----------+");
  }
  RSprintf("\n");
}

// Mu-group thetas are excluded from op_focei.scale's parameter vector, so
// scalePrintFun/scalePrintGrad's table never shows them; this adds one extra
// self-labeled row after each call (current regression value, or NA for the
// gradient print since these thetas have no FD gradient). Kept out of the
// shared scale.h to avoid touching its column layout; uses the same
// scale->every/cn throttle so it no-ops when printing is off or muModel="none".
static inline void printMuGroupThetaRow(bool isGrad) {
  if (op_focei.muModel == 0 || op_focei.muGroupN == 0) return;
  scaling *scale = &op_focei.scale;
  if (scale->every == 0 || scale->cn % scale->every != 0) return;
  RSprintf(isGrad ? "|   mu|   (no gradient - regression update)  |" :
                     "|   mu|                                     |");
  for (unsigned int g = 0; g < op_focei.muGroupN; g++) {
    std::string nm = (g < (unsigned int)op_focei.muGroupThetaNames.size() &&
                       op_focei.muGroupThetaNames[g] != NA_STRING) ?
      as<std::string>(op_focei.muGroupThetaNames[g]) : "?";
    if (isGrad) {
      RSprintf(" %8s:         NA |", nm.c_str());
    } else {
      RSprintf(" %8s: %#10.4g |", nm.c_str(), op_focei.fullTheta[op_focei.muGroupTheta[g]]);
    }
  }
  for (unsigned int c = 0; c < op_focei.muGroupCovN; c++) {
    std::string nm = (c < (unsigned int)op_focei.muGroupCovThetaNames.size() &&
                       op_focei.muGroupCovThetaNames[c] != NA_STRING) ?
      as<std::string>(op_focei.muGroupCovThetaNames[c]) : "?";
    if (isGrad) {
      RSprintf(" %8s:         NA |", nm.c_str());
    } else {
      RSprintf(" %8s: %#10.4g |", nm.c_str(), op_focei.fullTheta[op_focei.muGroupCovTheta[c]]);
    }
  }
  RSprintf("\n");
}

////////////////////////////////////////////////////////////////////////////////
// Outer l-BFGS-b from R
extern "C" double foceiOfvOptim(int n, double *x, void *ex){
  double ret = foceiOfv0(x);
  niter.push_back(op_focei.nF2+(++op_focei.nF));
  // Scaled
  vPar.push_back(ret);
  iterType.push_back(5);
  int finalize = 0, i = 0;
  for (i = 0; i < n; i++){
    vPar.push_back(x[i]);
  }
  // Unscaled
  iterType.push_back(6);
  niter.push_back(niter.back());
  if (op_focei.scaleObjective){
    vPar.push_back(op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
  } else {
    vPar.push_back(ret);
  }
  for (i = 0; i < n; i++){
    vPar.push_back(unscalePar(x, i));
  }
  // Back-transformed (7)
  iterType.push_back(7);
  niter.push_back(niter.back());
  if (op_focei.scaleObjective){
    vPar.push_back(op_focei.initObjective * ret / op_focei.scaleObjectiveTo);
  } else {
    vPar.push_back(ret);
  }
  for (i = 0; i < n; i++){
    int probitCode = (op_focei.scale.probitIdx != NULL) ? op_focei.scale.probitIdx[i] : 0;
    vPar.push_back(scaleBackTransform(unscalePar(x, i),
                                      op_focei.xPar[i], probitCode,
                                      op_focei.scale.logitThetaLow,
                                      op_focei.scale.logitThetaHi,
                                      op_focei.scale.probitThetaLow,
                                      op_focei.scale.probitThetaHi));
  }
  // Emit the per-iteration #/U/X rows via scale.h; when scaleObjective is on,
  // rescale back to the unscaled OFV so all rows show the number the user expects.
  double displayedOfv = op_focei.scaleObjective
    ? op_focei.initObjective * ret / op_focei.scaleObjectiveTo
    : ret;
  scalePrintFun(&op_focei.scale, x, displayedOfv);
  printMuGroupThetaRow(false);
  // One summary line per printed iteration for any ODE-solve warnings (e.g.
  // intdy window-misses) accumulated since the last flush; without this a
  // difficult model floods the console with identical lines each iteration.
  nmFlushRxSolveWarn(5);
  return ret;
}

//[[Rcpp::export]]
double foceiOuterF(NumericVector &theta){
  int n = theta.size();
  void *ex = NULL;
  return foceiOfvOptim(n, theta.begin(), ex);
}

extern "C" void outerGradNumOptim(int n, double *par, double *gr, void *ex){
  numericGrad(par, gr);
  op_focei.nG++;
  int finalize=0, i = 0;
  niterGrad.push_back(niter.back());
  if (op_focei.curAnalytic){
    gradType.push_back(9);
  } else if (op_focei.derivMethod == 0){
    if (op_focei.curGill == 1){
      gradType.push_back(1);
    } else if (op_focei.curGill == 2){
      gradType.push_back(5);
    } else if (op_focei.mixDeriv){
      gradType.push_back(2);
    } else{
      gradType.push_back(3);
    }
  } else {
    gradType.push_back(4);
  }
  // gradType convention: 1=Gill, 2=Mixed, 3=Forward, 4=Central, 5=Shi21,
  // 8=nlm forward sensitivity, 9=analytic outer gradient (fast=TRUE).
  scalePrintGrad(&op_focei.scale, gr, gradType.back());
  printMuGroupThetaRow(true);
  vGrad.push_back(NA_REAL); // Gradient doesn't record objf
  for (i = 0; i < n; i++){
    if (gr[i] == 0){
      if (op_focei.nF+op_focei.nF2 == 1) {
        if (op_focei.zeroGradFirstReset) {
          op_focei.zeroGrad=true;
          gr[i]=sqrt(DBL_EPSILON);
        } else {
          if (op_focei.zeroGradBobyqa) {
            stop("On initial gradient evaluation, one or more parameters have a zero gradient\ntrying outerOpt=\"bobyqa\")");
          }  else {
            stop("On initial gradient evaluation, one or more parameters have a zero gradient");
          }
        }
      } else {
        if (op_focei.zeroGradRunReset) {
          op_focei.zeroGrad=true;
          gr[i]=sqrt(DBL_EPSILON);
        } else {
          if (op_focei.zeroGradBobyqaRun) {
            stop("Zero gradient while searching, trying outerOpt=\"bobyqa\"");
          } else {
            stop("Zero gradient while searching");
          }
        }
      }
    }
    vGrad.push_back(gr[i]);
  }
}

//[[Rcpp::export]]
NumericVector foceiOuterG(NumericVector &theta){
  int n = theta.size();
  void *ex = NULL;
  NumericVector gr(n);
  outerGradNumOptim(n, theta.begin(), gr.begin(), ex);
  return gr;
}

void foceiLbfgsb3(Environment e){
  void *ex = NULL;
  double Fmin;
  int fail, fncount=0, grcount=0;
  NumericVector x(op_focei.npars);
  NumericVector g(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  char msg[100];
  std::fill_n(msg, 100, 0);
  lbfgsb3C(op_focei.npars, op_focei.lmm, x.begin(), op_focei.lower,
           op_focei.upper, op_focei.nbd, &Fmin, foceiOfvOptim,
           outerGradNumOptim, &fail, ex, op_focei.factr,
           op_focei.pgtol, &fncount, &grcount,
           op_focei.maxOuterIterations, msg, 0, -1,
           op_focei.abstol, op_focei.reltol, g.begin());
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = fail;
  e["message"] = msg;
  e["lastGrad"] = g;
}

void foceiLbfgsb(Environment e){
  void *ex = NULL;
  double Fmin;
  int fail, fncount=0, grcount=0;
  NumericVector x(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  char msg[100];
  lbfgsbRX(op_focei.npars, op_focei.lmm, x.begin(), op_focei.lower,
           op_focei.upper, op_focei.nbd, &Fmin, foceiOfvOptim,
           outerGradNumOptim, &fail, ex, op_focei.factr,
           op_focei.pgtol, &fncount, &grcount,
           op_focei.maxOuterIterations, msg, 0, op_focei.maxOuterIterations+1);
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = fail;
  e["message"] = msg;
}

void foceiCustomFun(Environment e){
  NumericVector x(op_focei.npars);
  NumericVector lower(op_focei.npars);
  NumericVector upper(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;){
    x[k]=scalePar(op_focei.initPar, k);
  }
  std::copy(&op_focei.upper[0], &op_focei.upper[0]+op_focei.npars, &upper[0]);
  std::copy(&op_focei.lower[0], &op_focei.lower[0]+op_focei.npars, &lower[0]);
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Function f = as<Function>(nlmixr2["foceiOuterF"]);
  Function g = as<Function>(nlmixr2["foceiOuterG"]);
  List ctl = e["control"];
  Function opt = as<Function>(ctl["outerOptFun"]);
  //.bobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...)
  List ret = as<List>(opt(_["par"]=x, _["fn"]=f, _["gr"]=g, _["lower"]=lower,
                          _["upper"]=upper,_["control"]=ctl));
  x = ret["x"];
  // Recalculate OFV in case the last calculated OFV isn't at the minimum....
  // Otherwise ETAs may be off
  if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
  // Finalize environment
  foceiOuterFinal(x.begin(), e);
  e["convergence"] = ret["convergence"];
  e["message"] = ret["message"];
  e["optReturn"] = ret;
}


////////////////////////////////////////////////////////////////////////////////
// Overall Outer Problem

//[[Rcpp::export]]
Environment foceiOuter(Environment e){
  op_focei.nF=0;
  op_focei.nG=0;
  op_focei.curAnalytic=0;
  op_focei.nAnalyticGrad=0;
  op_focei.nFDGradFast=0;
  op_focei.warnedAnalyticFallback=0;
  if (op_focei.maxOuterIterations > 0){
    for (unsigned int k = op_focei.npars; k--;){
      if (R_FINITE(op_focei.lower[k])){
        op_focei.lower[k]=scalePar(op_focei.lower, k);
      }
      if (R_FINITE(op_focei.upper[k])) {
        op_focei.upper[k]=scalePar(op_focei.upper,k);
      }
    }

    // Enable the analytic outer gradient only for the duration of the outer
    // optimizer's gradient callbacks; foceiS's own numericGrad use runs later
    // (in foceiFitCpp_) with the gate cleared.
    op_foceiFitEnv = e;
    op_foceiFitEnvSet = true;
    op_foceiUseAnalyticGrad = (op_focei.fast != 0);
    switch(op_focei.outerOpt){
    case 0:
      foceiLbfgsb(e);
      break;
    case 1:
      foceiLbfgsb3(e);
      break;
    case -1:
      foceiCustomFun(e);
    }
    op_foceiUseAnalyticGrad = false;
  } else {
    NumericVector x(op_focei.npars);
    for (unsigned int k = op_focei.npars; k--;){
      x[k]=scalePar(op_focei.initPar, k);
    }
    foceiOuterFinal(x.begin(), e);
    if (op_focei.maxInnerIterations == 0){
      e["fail"] = NA_INTEGER;
      if (!e.exists("message")) {
        e["message"] = "Likelihood evaluation with provided ETAs";
      }
    } else {
      e["fail"] = 0;
      if (!e.exists("message")) {
        e["message"] = "Posthoc prediction with provided THETAs";
      }
    }
  }
  return e;
}

//[[Rcpp::export]]
List nlmixr2Gill83_(Function what, NumericVector args, Environment envir,
                    LogicalVector which,
                    double gillRtol, int gillK=10, double gillStep=2, double gillFtol=0, bool optGillF=true){
  if (args.size()!=which.size()) stop("'args' must have same size as 'which'");
  gillRfn_=what;
  gillThetaN=args.size();
  gillRfnE_=envir;
  double *theta;
  theta = &args[0];
  NumericVector hfN(args.size());
  NumericVector hphifN(args.size());
  NumericVector gillDfN(args.size());
  NumericVector gillDf2N(args.size());
  NumericVector gillErrN(args.size());
  NumericVector aEps(args.size());
  NumericVector rEps(args.size());
  NumericVector aEpsC(args.size());
  NumericVector rEpsC(args.size());
  IntegerVector retN(args.size());
  gillLong=false;
  NumericVector fN(args.size());
  double gillF;
  for (int i = args.size(); i--;){
    if (which[i]){
      gillPar=i;
      if (i == args.size()-1 || gillLong){
        gillF = gillRfn(theta);
      }
      fN[i] = gillF;
      retN[i] = gill83(&hfN[i], &hphifN[i], &gillDfN[i], &gillDf2N[i], &gillErrN[i],
                       theta, i, gillRtol, gillK, gillStep,
                       gillFtol,
                       -1, gill83fnG, 0, gillF) + 1;
      double err=1/(std::fabs(theta[i])+1);
      aEps[i]  = hfN[i]*err;
      rEps[i]  = hfN[i]*err;
      if(optGillF){
        aEpsC[i] = hfN[i]*err;
        rEpsC[i] = hfN[i]*err;
      } else {
        aEpsC[i] = hphifN[i]*err;
        rEpsC[i] = hphifN[i]*err;
      }
    } else {
      retN[i] = 1;
      hfN[i] = NA_REAL;
      hphifN[i] = NA_REAL;
      gillDfN[i] = NA_REAL;
      gillDf2N[i] = NA_REAL;
      gillErrN[i] = NA_REAL;
      aEps[i]  = NA_REAL;
      rEps[i]  = NA_REAL;
      aEpsC[i]  = NA_REAL;
      rEpsC[i]  = NA_REAL;
    }
  }
  List df(11);
  retN.attr("levels") = CharacterVector::create("Not Assessed","Good","High Grad Error",
                                                "Constant Grad","Odd/Linear Grad",
                                                "Grad changes quickly");
  retN.attr("class") = "factor";
  df[0] = retN;
  df[1] = hfN;
  df[2] = hphifN;
  df[3] = gillDfN;
  df[4] = gillDf2N;
  df[5] = gillErrN;
  df[6] = aEps;
  df[7] = rEps;
  df[8] = aEpsC;
  df[9] = rEpsC;
  df[10] = fN;
  df.attr("names") = CharacterVector::create("info","hf","hphi","df","df2","err","aEps","rEps","aEpsC","rEpsC","f");
  if (args.hasAttribute("names")){
    df.attr("row.names") = args.attr("names");
  } else {
    df.attr("row.names") = IntegerVector::create(NA_INTEGER, -args.size());
  }
  CharacterVector cls = CharacterVector::create("nlmixr2Gill83", "data.frame");
  List info = List::create(_["which"]=which,
                           _["gillRtol"]=gillRtol,
                           _["gillK"]=gillK,
                           _["gillStep"]=gillStep,
                           _["gillFtol"]=gillFtol,
                           _["optGillF"]=optGillF);
  info.attr("class") = CharacterVector::create("nlmixr2LstSilent");
  cls.attr(".nlmixr2Gill") = info;
  df.attr("class") = cls;
  return df;
}
//' @rdname nlmixr2GradFun
//' @export
//[[Rcpp::export]]
double nlmixr2Eval_(NumericVector theta, std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Environment gradInfo = nlmixr2[".nlmixr2GradInfo"];
  std::string EF = md5 + ".f";
  std::string EE = md5 + ".e";
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  if (!gradInfo.exists(EW)){
    LogicalVector tmp(theta.size());
    for (int i = theta.size(); i--;){
      tmp[i] = true;
    }
    gradInfo[EW] = tmp;
  }
  lEW=gradInfo[EW];
  if (lEW.size() != theta.size()) stop("invalid theta size");
  Function cFun = as<Function>(gradInfo[EF]);
  Environment cEnvir = as<Environment>(gradInfo[EE]);
  double f0;
  List par(1);
  par[0] = theta;
  f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  std::string f0s = md5 + ".fc";
  std::string f0t = md5 + ".ft";
  std::string cns = md5 + ".n";
  f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  gradInfo[f0s] = f0;
  gradInfo[f0t] = theta;
  int cn = gradInfo[cns]; cn++;
  gradInfo[cns] = cn;
  bool useColor = as<bool>(gradInfo["useColor"]);
  int printNcol=as<int>(gradInfo["printNcol"]);
  int printN=as<int>(gradInfo["print"]);
  int i, finalize=0, n=theta.size();
  bool isRstudio=as<bool>(gradInfo["isRstudio"]);
  if (cn == 1){
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
    if (printN != 0){
      foceiPrintLine(min2(n, printNcol));
      if (gradInfo.exists("thetaNames")){
        CharacterVector tn;
        tn = gradInfo["thetaNames"];
        if (tn.size()!=lEW.size()){
          CharacterVector tn2(lEW.size());
          for (int i = 0; i < lEW.size(); i++){
            tn2[i] = "t" + std::to_string(i+1);
          }
          gradInfo["thetaNames"]=tn2;
        }
      } else {
        CharacterVector tn(lEW.size());
        for (int i = 0; i < lEW.size(); i++){
          tn[i] = "t" + std::to_string(i+1);
        }
        gradInfo["thetaNames"]=tn;
      }
      CharacterVector thetaNames = gradInfo["thetaNames"];
      RSprintf("|    #| Objective Fun |");
      int i=0, finalize=0;
      std::string tmpS;
      for (i = 0; i < n; i++){
        tmpS = thetaNames[i];
        RSprintf("%#10s |", tmpS.c_str());
        if ((i + 1) != n && (i + 1) % printNcol == 0){
          if (useColor && printNcol + i  >= n){
            RSprintf("\n\033[4m|.....................|");
          } else {
            RSprintf("\n|.....................|");
          }
          finalize=1;
        }
      }
      if (finalize){
        while(true){
          if ((i++) % printNcol == 0){
            if (useColor) RSprintf("\033[0m");
            RSprintf("\n");
            break;
          } else {
            RSprintf("...........|");
          }
        }
      } else {
        RSprintf("\n");
      }
    }
  }
  bool doUnscaled = false;
  std::string unscaledPar = md5 + ".uPar";
  NumericVector thetaU;
  niter.push_back(cn);
  // Scaled
  vPar.push_back(f0);
  if (gradInfo.exists(unscaledPar)){
    thetaU=as<NumericVector>(gradInfo[unscaledPar]);
    if (thetaU.size() != theta.size()){
      iterType.push_back(6);
    } else {
      doUnscaled=true;
      iterType.push_back(5);
    }
  } else {
    // Actually unscaled
    iterType.push_back(6);
  }
  for (i = 0; i < n; i++){
    vPar.push_back(theta[i]);
  }
  if (printN != 0 && cn % printN == 0){
    if (useColor && isRstudio)
      RSprintf("|\033[1m%5d\033[0m|%#14.8g |", cn, f0);
    else
      RSprintf("|%5d|%#14.8g |", cn, f0);
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g |", theta[i]);
      if ((i + 1) != n && (i + 1) % printNcol == 0){
        if (useColor && printNcol + i  > n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % printNcol == 0){
          if (useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
  }
  if (doUnscaled){
    iterType.push_back(6);
    niter.push_back(niter.back());
    finalize=0;
    // No obj scaling currently
    vPar.push_back(f0);
    for (i = 0; i < n; i++){
      vPar.push_back(thetaU[i]);
    }
    if (printN != 0 && cn % printN == 0){
      if (useColor && isRstudio)
        RSprintf("|    U|%#14.8g |", f0);
      else
        RSprintf("|    U|%#14.8g |", f0);
      for (i = 0; i < n; i++){
        RSprintf("%#10.4g |", thetaU[i]);
        if ((i + 1) != n && (i + 1) % printNcol == 0){
          if (useColor && printNcol + i  > n){
            RSprintf("\n\033[4m|.....................|");
          } else {
            RSprintf("\n|.....................|");
          }
          finalize=1;
        }
      }
      if (finalize){
        while(true){
          if ((i++) % printNcol == 0){
            if (useColor) RSprintf("\033[0m");
            RSprintf("\n");
            break;
          } else {
            RSprintf("...........|");
          }
        }
      } else {
        RSprintf("\n");
      }
    }
  }
  return f0;
}

void nlmixr2GradPrint(NumericVector gr, int gradType, int cn, bool useColor,
                      int printNcol, int printN, bool isRstudio){
  int n = gr.size(), finalize=0, i;
  if (printN != 0 && cn % printN == 0){
    if (useColor && printNcol >= n){
      switch(gradType){
      case 1:
        RSprintf("|\033[4m    G|    Gill Diff. |");
        break;
      case 2:
        RSprintf("|\033[4m    M|   Mixed Diff. |");
        break;
      case 3:
        RSprintf("|\033[4m    F| Forward Diff. |");
        break;
      case 4:
        RSprintf("|\033[4m    C| Central Diff. |");
        break;
      }
    } else {
      switch(gradType){
      case 1:
        RSprintf("|    G|    Gill Diff. |");
        break;
      case 2:
        RSprintf("|    M|   Mixed Diff. |");
        break;
      case 3:
        RSprintf("|    F| Forward Diff. |");
        break;
      case 4:
        RSprintf("|    C| Central Diff. |");
        break;
      }
    }
    for (i = 0; i < n; i++){
      RSprintf("%#10.4g ", gr[i]);
      if (useColor && printNcol >= n && i == n-1){
        RSprintf("\033[0m");
      }
      RSprintf("|");
      if ((i + 1) != n && (i + 1) % printNcol == 0){
        if (useColor && printNcol + i  >= n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % printNcol == 0){
          if (useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
    if (!useColor){
      foceiPrintLine(min2(n, printNcol));
    }
  }
}

//' @rdname nlmixr2GradFun
//' @export
//[[Rcpp::export]]
RObject nlmixr2Unscaled_(NumericVector theta, std::string md5){
  // Unscaled
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Environment gradInfo = nlmixr2[".nlmixr2GradInfo"];
  std::string unscaledPar = md5 + ".uPar";
  gradInfo[unscaledPar] = theta;
  return R_NilValue;
}

//' @rdname nlmixr2GradFun
//' @export
//[[Rcpp::export]]
NumericVector nlmixr2Grad_(NumericVector theta, std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Environment gradInfo = nlmixr2[".nlmixr2GradInfo"];

  std::string Egill = md5 + ".g";
  std::string EF = md5 + ".f";
  std::string EE = md5 + ".e";

  Function cFun = as<Function>(gradInfo[EF]);
  Environment cEnvir = as<Environment>(gradInfo[EE]);
  List Lgill;
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  bool useColor = as<bool>(gradInfo["useColor"]);
  int printNcol=as<int>(gradInfo["printNcol"]);
  int printN=as<int>(gradInfo["print"]);
  bool isRstudio=as<bool>(gradInfo["isRstudio"]);
  if (gradInfo.exists(Egill)){
    lEW=gradInfo[EW];
    if (lEW.size() != theta.size()){
      stop("Invalid theta size");
    }
    Lgill = gradInfo[Egill];
  } else {
    if (!gradInfo.exists(EW)){
      LogicalVector tmp(theta.size());
      for (int i = theta.size(); i--;){
        tmp[i] = true;
      }
      gradInfo[EW] = tmp;
    }
    lEW=gradInfo[EW];
    if (lEW.size() != theta.size()){
      stop("Invalid theta size (or which size)");
    }
    std::string Ertol = md5 + ".rtol";
    std::string EK = md5 + ".k";
    std::string Estep = md5 + ".s";
    std::string EFtol = md5 + ".ftol";
    Lgill = nlmixr2Gill83_(cFun, theta, gradInfo[EE],
                           lEW, gradInfo[Ertol],
                           gradInfo[EK], gradInfo[Estep],
                           gradInfo[EFtol]);
    gradInfo[Egill]=Lgill;
    niterGrad.push_back(niter.back());
    gradType.push_back(1);
    vGrad.push_back(NA_REAL); // Gradient doesn't record objf
    NumericVector gr = as<NumericVector>(Lgill["df"]);
    for (int i = 0; i < gr.size(); i++){
      if (gr[i] == 0){
        stop("On initial gradient evaluation, one or more parameters have a zero gradient\nChange model, try different initial estimates or try derivative free optimization)");
      }
      vGrad.push_back(gr[i]);
    }
    nlmixr2GradPrint(gr, gradType.back(), niter.back(), useColor,
                     printNcol, printN, isRstudio);
    return gr;
  }
  NumericVector aEps = as<NumericVector>(Lgill["aEps"]);
  NumericVector rEps = as<NumericVector>(Lgill["rEps"]);
  NumericVector aEpsC = as<NumericVector>(Lgill["aEpsC"]);
  NumericVector rEpsC = as<NumericVector>(Lgill["rEpsC"]);
  NumericVector g(theta.size());
  double f0, delta, cur, tmp=0, tmp0;
  bool doForward=true;
  // FIXME
  List par(1);
  par[0] = theta;
  std::string f0s = md5 + ".fc";
  std::string f0t = md5 + ".ft";
  bool reEval = true;
  if (gradInfo.exists(f0s)){
    NumericVector thetaL = gradInfo[f0t];
    if (thetaL.size() == theta.size()){
      reEval=false;
      for (int i = theta.size(); i--;){
        if (thetaL[i] != theta[i]){
          reEval=true;
          break;
        }
      }
      if (!reEval){
        NumericVector tmp  = gradInfo[f0s];
        if (tmp.size() == 1){
          f0 = tmp[0];
        } else {
          reEval=true;
        }
      }
    }
  }
  if (reEval){
    f0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
  }
  niterGrad.push_back(niter.back());
  vGrad.push_back(NA_REAL); // Gradient doesn't record objf
  bool isMixed=false;
  for (int i = theta.size(); i--;){
    cur = theta[i];
    if (doForward){
      delta = (std::fabs(theta[i])*rEps[i] + aEps[i]);
      theta[i] = cur + delta;
      par[0] = theta;
      tmp = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      g[i] = (tmp-f0)/delta;
      theta[i] = cur;
    } else {
      delta = (std::fabs(theta[i])*rEpsC[i] + aEpsC[i]);
      theta[i] = cur + delta;
      par[0] = theta;
      tmp0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      theta[i] = cur - delta;
      par[0] = theta;
      tmp = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
      g[i] = (tmp0-tmp)/(2*delta);
      theta[i] = cur;
    }

    // Check for bad grad
    if (std::isnan(g[i]) ||  ISNA(g[i]) || !R_FINITE(g[i])){
      if (doForward){
        // Switch to Backward difference method
        // op_focei.mixDeriv=1;
        theta[i] = cur - delta;
        par[0] = theta;
        tmp0 = as<double>(doCall(_["what"] = cFun, _["args"]=par, _["envir"]=cEnvir));
        g[i] = (f0-tmp0)/(delta);
        isMixed=true;
      } else {
        // We are using the central difference AND there is an NA in one of the terms
        // g[cpar] = (tmp0-tmp)/(2*delta);
        // op_focei.mixDeriv=1;
        isMixed=true;
        if (std::isnan(tmp0) || ISNA(tmp0) || !R_FINITE(tmp0)){
          // Backward
          g[i] = (f0-tmp)/delta;
        } else {
          // Forward
          g[i] = (tmp0-f0)/delta;
        }
      }
    }
  }
  for (int i = 0; i < theta.size(); i++){
    vGrad.push_back(g[i]);
  }
  if (isMixed){
    gradType.push_back(2);
  } else if (doForward) {
    gradType.push_back(3);
  } else {
    gradType.push_back(4);
  }
  nlmixr2GradPrint(g, gradType.back(), niter.back(), useColor,
                   printNcol, printN, isRstudio);
  return g;
}
//' @rdname nlmixr2GradFun
//' @export
//[[Rcpp::export]]
RObject nlmixr2ParHist_(std::string md5){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Environment gradInfo = nlmixr2[".nlmixr2GradInfo"];
  std::string EW = md5 + ".w";
  LogicalVector lEW;
  if (!gradInfo.exists(EW)){
    LogicalVector tmp(lEW.size());
    for (int i = lEW.size(); i--;){
      tmp[i] = true;
    }
    gradInfo[EW] = tmp;
  }
  lEW=gradInfo[EW];
  if (gradInfo.exists("thetaNames")){
    CharacterVector tn;
    tn = gradInfo["thetaNames"];
    if (tn.size()!=lEW.size()){
      CharacterVector tn2(lEW.size());
      for (int i = 0; i < lEW.size(); i++){
        tn2[i] = "t" + std::to_string(i+1);
      }
      gradInfo["thetaNames"]=tn2;
    }
  } else {
    CharacterVector tn(lEW.size());
    for (int i = 0; i < lEW.size(); i++){
      tn[i] = "t" + std::to_string(i+1);
    }
    gradInfo["thetaNames"]=tn;
  }
  std::string cns = md5 + ".n";
  gradInfo[cns] = 0;
  parHistData(gradInfo, false);
  return gradInfo["parHistData"];
}

//[[Rcpp::export]]
RObject nlmixr2Hess_(RObject thetaT, RObject fT, RObject e,
                     RObject gillInfoT){
  List par(1);
  NumericVector theta = as<NumericVector>(thetaT);
  Function f = as<Function>(fT);
  List gillInfo = as<List>(gillInfoT);
  arma::mat H(theta.size(), theta.size(), fill::zeros);
  double epsI, epsJ;
  NumericVector rEpsC = as<NumericVector>(gillInfo["rEpsC"]);
  NumericVector aEpsC = as<NumericVector>(gillInfo["aEpsC"]);
  NumericVector nF = as<NumericVector>(gillInfo["f"]);
  double lastOfv=nF[0];
  int n = theta.size();
  double f1,f2,f3,f4;
  double ti, tj;
  int i, j;
  int totTick= 4*n +2*n*(n-1);
  int cur = 0, curTick=0;
  clock_t t0=clock();
  for (i=n; i--;){
    epsI = (std::fabs(theta[i])*rEpsC[i] + aEpsC[i]);
    ti = theta[i];
    theta[i] = ti + 2*epsI;
    par[0]=theta;
    f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti + epsI;
    par[0]=theta;
    f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - epsI;
    par[0]=theta;
    f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti - 2*epsI;
    par[0]=theta;
    f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
    cur++;
    curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
    theta[i] = ti;
    H(i,i)=(-f1+16*f2-30*lastOfv+16*f3-f4)/(12*epsI*epsI);
    for (j = i; j--;){
      epsJ = (std::fabs(theta[j])*rEpsC[j] + aEpsC[j]);
      // eps = sqrt(epsI*epsJ);// 0.5*epsI+0.5*epsJ;
      // epsI = eps;
      // epsJ = eps;
      tj = theta[j];
      theta[i] = ti + epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f1 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti + epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f2 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj + epsJ;
      par[0]=theta;
      f3 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      theta[i] = ti - epsI;
      theta[j] = tj - epsJ;
      par[0]=theta;
      f4 = as<double>(doCall(_["what"] = f, _["args"]=par, _["envir"]=e));
      cur++;
      curTick = par_progress(cur, totTick, curTick, 1, t0, 0);
      H(i,j)= (f1-f2-f3+f4)/(4*epsI*epsJ);
      H(j,i) = H(i,j);
      theta[i] = ti;
      theta[j] = tj;
    }
  }
  par_progress(totTick, totTick, cur, 1, t0, 0);
  if (isRstudio){
    RSprintf("\n");
  } else {
    RSprintf("\r                                                                                \r");
  }
  return wrap(H);
}

////////////////////////////////////////////////////////////////////////////////
// Covariance functions

// R-callable bridge to shi21Central (analytic-cov 3rd-order tensor): `f` returns
// the vector to difference (NULL/short -> NaN, which shi21Central tolerates).
// Serial use only (single static holder); the covariance step is not parallel.
static Rcpp::Function *shiWrapFn_ = NULL;
static unsigned int shiWrapN_ = 0;
static arma::vec shiWrapCall_(arma::vec &t, int id) {
  arma::vec bad(shiWrapN_); bad.fill(NA_REAL);
  RObject res = (*shiWrapFn_)(NumericVector(t.begin(), t.end()));
  if (Rf_isNull(res)) return bad;
  NumericVector rv(res);
  if ((unsigned int)rv.size() != shiWrapN_) return bad;
  return arma::vec(rv.begin(), rv.size());
}
struct ShiWrapGuard { ~ShiWrapGuard() { shiWrapFn_ = NULL; } };
//[[Rcpp::export]]
Rcpp::List shi21CentralWrap(Rcpp::Function f, arma::vec t, arma::vec f0, int idx, double ef) {
  ShiWrapGuard _g;                                  // clear the holder even on an R error
  shiWrapFn_ = &f; shiWrapN_ = f0.n_elem;
  arma::vec gr(f0.n_elem, arma::fill::zeros);
  double h = 0.0;
  h = shi21Central(shiWrapCall_, t, h, f0, gr, 0, idx - 1, ef);
  return Rcpp::List::create(_["h"] = h, _["gr"] = gr);
}

int foceiCalcR(Environment e){
  rx = getRxSolve_();
  // covType="analytic": the exact analytic observed-information R-matrix while the
  // optimizer is live; `.foceiCalcRanalytic` returns the npars x npars R or NULL to
  // fall through to the finite-difference Hessian below.
  {
    List _ctl = as<List>(e["control"]);
    std::string _covType = _ctl.containsElementNamed("covType") ?
      as<std::string>(_ctl["covType"]) : "fd";
    if (_covType == "analytic") {
      Environment nlmixr2est = Environment::namespace_env("nlmixr2est");
      Function af = as<Function>(nlmixr2est[".foceiCalcRanalytic"]);
      RObject res = af(e);
      if (!Rf_isNull(res)) {
        arma::mat H0 = as<arma::mat>(res);
        e["R.0"] = wrap(H0);
        arma::mat cholR0, RE0;
        bool rpd0 = cholSE0(cholR0, RE0, H0, op_focei.cholSEtol);
        e["R.pd"] = wrap(rpd0);
        e["R.E"]  = wrap(RE0);
        e["cholR"] = wrap(cholR0);
        // the augmented sensitivity solves replaced the fit's global solve; restore it
        // so foceiFinalizeTables (llikObs, tolFactor) reads the fit, not the last
        // subject.  If the restore fails the global solve is unusable -> abort the cov
        // step (foceiCalcCov honors covMethod==0 as a clean skip) rather than let
        // finalize read a dangling solve.
        if (!restoreFitSolve_()) {
          op_focei.covMethod = 0;
          return 0;
        }
        return 1;
      }
      // analytic declined -> FD Hessian below.  Restore the freed global solve BEFORE
      // the warning (options(warn=2) longjmps out); a failed restore -> failed cov,
      // don't run FD against a freed solve.
      if (e.exists(".analyticStarted") && as<bool>(e[".analyticStarted"])) {
        if (!restoreFitSolve_()) {
          op_focei.covMethod = 0;
          return 0;
        }
      }
      // tell the user (not silent).  RSprintf is the visible channel (like "Could not
      // calculate covariance matrix"); the warning condition is for programmatic capture.
      RSprintf("\rcovType=\"analytic\" not available for this model (out of scope, or "
               "the augmented model would not build/solve); using the finite-difference "
               "sandwich (\"r,s\") covariance.\n");
      Rf_warning("covType=\"analytic\": the analytic covariance is not available for "
                 "this model; used the finite-difference sandwich (\"r,s\") covariance instead.");
      // analytic requested but unavailable -> fall back to the finite-difference SANDWICH
      // ("r,s"), not the R-matrix alone ("r"): covMethod=2 was only the internal slot the
      // "analytic" token maps to.  The S block in foceiCalcCov fires once this is 1.
      if (op_focei.covMethod == 2) op_focei.covMethod = 1;
    }
  }
  arma::mat H(op_focei.npars, op_focei.npars);
  arma::vec theta(op_focei.npars);
  unsigned int i, j, k;
  for (k = op_focei.npars; k--;){
    j=op_focei.fixedTrans[k];
    theta[k] = op_focei.fullTheta[j];
  }

  // arma::vec df1(op_focei.npars);
  // arma::vec df2(op_focei.npars);
  double epsI, epsJ;

  bool doForward=false;
  if (op_focei.derivMethod == 0){
    doForward=true;
  }
  double f1,f2,f3,f4;
  double ti, tj;
  if (doForward){
    stop("Not implemented for finite differences.");
  } else {
    //https://pdfs.semanticscholar.org/presentation/e7d5/aff49eb17fd155e75725c295859d983cfda4.pdf
    // https://v8doc.sas.com/sashtml/ormp/chap5/sect28.htm
    double fnscale = 1.0;
    if (op_focei.scaleObjective == 2){
      fnscale = op_focei.initObjective / op_focei.scaleObjectiveTo;
    }
    double parScaleI=1.0, parScaleJ=1.0;
    for (i=op_focei.npars; i--;){
      epsI = (std::fabs(theta[i])*op_focei.rEpsC[i] + op_focei.aEpsC[i]);
      // > (45/12)^(1/5)
      // [1] 1.302585542348676073132
      // based on central error to stencil error which is closer to optimal for this difference
      // epsI = pow(epsI, 3.0/5.0)*1.302585542348676073132;
      ti = theta[i];
      theta[i] = ti + 2*epsI;
      updateTheta(theta.begin());
      f1 = foceiOfv0(theta.begin());
      if (ISNA(f1)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti + epsI;
      updateTheta(theta.begin());
      f2 = foceiOfv0(theta.begin());
      if (ISNA(f2)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti - epsI;
      updateTheta(theta.begin());
      f3 = foceiOfv0(theta.begin());
      if (ISNA(f3)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti - 2*epsI;
      updateTheta(theta.begin());
      f4 = foceiOfv0(theta.begin());
      if (ISNA(f4)) return 0;
      op_focei.cur++;
      op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
      theta[i] = ti;
      // RSprintf("-- i:%d, i: %d\n", i, i);
      // print(NumericVector::create(f1,f2,f3,f4,op_focei.lastOfv));
      H(i,i)=fnscale*(-f1+16*f2-30*op_focei.lastOfv+16*f3-f4)/(12*epsI*epsI*parScaleI*parScaleI);
      for (j = i; j--;){
        epsJ = (std::fabs(theta[j])*op_focei.rEpsC[j] + op_focei.aEpsC[j]);
        epsI = (std::fabs(theta[i])*op_focei.rEpsC[i] + op_focei.aEpsC[i]);
        // eps = sqrt(epsI*epsJ);// 0.5*epsI+0.5*epsJ;
        // epsI = eps;
        // epsJ = eps;
        tj = theta[j];
        theta[i] = ti + epsI;
        theta[j] = tj + epsJ;
        updateTheta(theta.begin());
        f1 = foceiOfv0(theta.begin());
        if (ISNA(f1)) return 0;
        op_focei.cur++;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        theta[i] = ti + epsI;
        theta[j] = tj - epsJ;
        updateTheta(theta.begin());
        f2 = foceiOfv0(theta.begin());
        if (ISNA(f2)) return 0;
        op_focei.cur++;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        theta[i] = ti - epsI;
        theta[j] = tj + epsJ;
        updateTheta(theta.begin());
        f3 = foceiOfv0(theta.begin());
        if (ISNA(f3)) return 0;
        op_focei.cur++;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        theta[i] = ti - epsI;
        theta[j] = tj - epsJ;
        updateTheta(theta.begin());
        f4 = foceiOfv0(theta.begin());
        if (ISNA(f4)) return 0;
        op_focei.cur++;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        // RSprintf("-- i:%d, j: %d\n", i, j);
        // print(NumericVector::create(f1,f2,f3,f4));
        H(i,j)= fnscale*(f1-f2-f3+f4)/(4*epsI*epsJ*parScaleI*parScaleJ);
        H(j,i) = H(i,j);
        theta[i] = ti;
        theta[j] = tj;
      }
    }
  }
  // R matrix = Hessian/2
  H = H*0.5;
  // https://github.com/cran/nmw/blob/59478fcc91f368bb3bbc23e55d8d1d5d53726a4b/R/CovStep.R
  // H = 0.25*H + 0.25*H.t();
  if (e.exists("R.1")){
    // This is the 2nd attempt
    arma::mat H2 = 0.5*H + 0.5*as<arma::mat>(e["R.1"]);
    arma::mat cholR;
    arma::mat RE;
    bool rpd = cholSE0(cholR, RE, H2, op_focei.cholSEtol);
    if (rpd){
      e["R.pd"] =  rpd;
      e["R.E"] =  wrap(RE);
      e["cholR"] = wrap(cholR);
    } else {
      e["R.pd2"] = false;
      e["R.2"] = H2;
      e["R.E2"] = wrap(RE);
      e["cholR2"] = wrap(cholR);
      e["R.pd"] = cholSE0(cholR, RE, H, op_focei.cholSEtol);
      e["R.E"] =  wrap(RE);
      e["cholR"] = wrap(cholR);
    }
  } else {
    e["R.0"] = H;
    arma::mat cholR;
    arma::mat RE;
    bool rpd = cholSE0(cholR, RE, H, op_focei.cholSEtol);
    e["R.pd"] =  wrap(rpd);
    e["R.E"] =  wrap(RE);
    e["cholR"] = wrap(cholR);
  }
  return 1;
}

// Necessary for S-matrix calculation
int foceiS(double *theta, Environment e, bool &hasZero){
  int npars = op_focei.npars;
  int oldCalcGrad = op_focei.calcGrad;
  op_focei.calcGrad = 1;
  arma::vec gfull(npars);
  numericGrad(theta, gfull.memptr());
  op_focei.calcGrad = oldCalcGrad;
  hasZero = false;
  rx = getRxSolve_();
  op_focei.calcGrad=1;
  int cpar, gid;
  double cur, delta;
  focei_ind *fInd;
  // Do Forward difference if the OBJF for *theta has already been calculated.
  bool doForward=false;
  if (op_focei.derivMethod == 0){
    doForward=true;
    // If the first derivative wasn't calculated, then calculate it.
    for (cpar = npars; cpar--;){
      if (theta[cpar] != op_focei.theta[cpar]){
        doForward=false;
        break;
      }
    }
    if (doForward){
      // Fill in lik0
      for (gid = getRxNsub(rx); gid--;){
        fInd = &(inds_focei[gid]);
        op_focei.likSav[gid] = -2*fInd->lik[0];
      }
    }
  }
  bool smatNorm = op_focei.smatNorm;
  if (op_focei.needOptimHess) {
    smatNorm = op_focei.smatNormLlik;
  }
  double sInfoPer = npars * getRxNsub(rx);
  for (cpar = npars; cpar--;){
    double rEps = op_focei.rEps[cpar];
    double rEpsC = op_focei.rEpsC[cpar];
    if (smatNorm){
      if (doForward){
        delta = (std::fabs(theta[cpar])*rEps + op_focei.aEps[cpar])/_safe_sqrt(1+std::fabs(min2(op_focei.initObjective, op_focei.lastOfv)));
      } else {
        delta = (std::fabs(theta[cpar])*rEpsC + op_focei.aEpsC[cpar])/_safe_sqrt(1+std::fabs(min2(op_focei.initObjective, op_focei.lastOfv)));
      }
    } else {
      if (doForward){
        delta = std::fabs(theta[cpar])*rEps + op_focei.aEps[cpar];
      } else {
        delta = std::fabs(theta[cpar])*rEpsC + op_focei.aEpsC[cpar];
      }
    }
    if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0); // All etas = -42;  Unlikely if normal
    cur = theta[cpar];
    theta[cpar] = cur + delta;
    updateTheta(theta);
    {
      int _nsub = (int)getRxNsub(rx);
      rx_solving_options *_op = getSolvingOptions(rx);
      int _cores = getOpCores(_op);
      bool _doParallel = (_cores > 1) && solveMethodThreadSafe(_op);
      std::vector<int> _opt1Res(_nsub, 0);
      if (_doParallel) sortIds(rx, 2);
      _innerParallel.store(1, std::memory_order_release);
#ifdef _OPENMP
#pragma omp parallel for num_threads(_cores) schedule(dynamic) if(_doParallel)
#endif
      for (int _i = 0; _i < _nsub; _i++) {
        int _gid = _doParallel ? (getOrdId(rx, _i) - 1) : _i;
        focei_ind *fIndL = &(inds_focei[_gid]);
        fIndL->thetaGrad[cpar] = NA_REAL;
        // Set thread id for windows
        setRxThreadId(omp_get_thread_num());
        _opt1Res[_gid] = innerOpt1(_gid, 2);
        setRxThreadId(-1);
        if (_opt1Res[_gid] && doForward) {
          fIndL->thetaGrad[cpar] = (fIndL->lik[2] - op_focei.likSav[_gid]) / delta;
        }
      }
      _innerParallel.store(0, std::memory_order_release);
      if (_doParallel) sortIds(rx, 0);
      // Serial fallback for subjects where innerOpt1(gid, 2) failed
      for (int _gid = 0; _gid < _nsub; _gid++) {
        if (!_opt1Res[_gid]) {
          fInd = &(inds_focei[_gid]);
          if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0);
          theta[cpar] = cur - delta;
          updateTheta(theta);
          if (!innerOpt1(_gid, 2)) {
            hasZero = true;
            sInfoPer -= 1.0;
            fInd->thetaGrad[cpar] = gfull[cpar];
          }
          theta[cpar] = cur + delta;
          updateTheta(theta);
          fInd->thetaGrad[cpar] = (op_focei.likSav[_gid] - fInd->lik[2]) / delta;
        }
      }
    }
    if (!doForward){
      if (op_focei.neta != 0) std::fill_n(&op_focei.goldEta[0], op_focei.gEtaGTransN, -42.0);
      theta[cpar] = cur - delta;
      updateTheta(theta);
      // Second inner loop: run innerOpt1(gid, 1) over subjects in parallel.
      {
        int _nsub = (int)getRxNsub(rx);
        rx_solving_options *_op = getSolvingOptions(rx);
        int _cores = getOpCores(_op);
        bool _doParallel = (_cores > 1) && solveMethodThreadSafe(_op);
        if (_doParallel) sortIds(rx, 2);
        _innerParallel.store(1, std::memory_order_release);
#ifdef _OPENMP
#pragma omp parallel for num_threads(_cores) schedule(dynamic) if(_doParallel)
#endif
        for (int _i = 0; _i < _nsub; _i++) {
          int _gid = _doParallel ? (getOrdId(rx, _i) - 1) : _i;
          focei_ind *fIndL = &(inds_focei[_gid]);
          if (ISNA(fIndL->thetaGrad[cpar])) {
            setRxThreadId(omp_get_thread_num());
            if (!innerOpt1(_gid, 1)) {
              // forward only
              fIndL->thetaGrad[cpar] = (fIndL->lik[2] - op_focei.likSav[_gid]) / delta;
            } else {
              // central
              fIndL->thetaGrad[cpar] = (fIndL->lik[2] - fIndL->lik[1]) / (2*delta);
            }
            setRxThreadId(-1);
          }
        }
        _innerParallel.store(0, std::memory_order_release);
        if (_doParallel) sortIds(rx, 0);
      }
    }
    theta[cpar] = cur;
    op_focei.cur++;
    op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
  }
  op_focei.calcGrad=0;
  // Now calculate S matrix
  arma::mat m1(1, op_focei.npars), S(op_focei.npars, op_focei.npars, fill::zeros), s1(1, op_focei.npars,fill::ones);
  for (gid = getRxNsub(rx); gid--;){
    fInd = &(inds_focei[gid]);
    std::copy(&fInd->thetaGrad[0],&fInd->thetaGrad[0]+op_focei.npars,&m1[0]);
    S = S + m1.t() * m1;
  }
  // S matrix = S/4
  // According to https://github.com/cran/nmw/blob/59478fcc91f368bb3bbc23e55d8d1d5d53726a4b/R/Objs.R
  S=S*0.25;
  e["S0"] = wrap(S);
  sInfoPer = sInfoPer / (npars * getRxNsub(rx));
  e["Sper"] = sInfoPer;
  // fixme hard coded
  if (sInfoPer < op_focei.smatPer) {
    return 0;
  }
  arma::mat cholS;
  arma::mat SE;
  e["S.pd"] =  cholSE0(cholS, SE, S, op_focei.cholSEtol);
  e["S.E"] =  wrap(SE);
  e["cholS"] = wrap(cholS);
  return 1;
}
//' Return the square root of general square matrix A
//'
//' @param m Matrix to take the square root of.
//'
//' @return A square root general square matrix of m
//'
//' @export
//[[Rcpp::export]]
NumericMatrix sqrtm(NumericMatrix m){
  arma::cx_mat ret = sqrtmat(as<arma::mat>(m));
  mat im = arma::imag(ret);
  mat re = arma::real(ret);
  if (arma::any(arma::any(im,0))){
    stop("Some components of sqrtm are imaginary.");
  }
  return wrap(re);
}


void setupAq0_(Environment e) {
  if (e.exists("aqn")) {
    _aqn = (unsigned int)as<int>(e["aqn"]);
    _nagq = as<int>(e["nAGQ"]);
  } else {
    _aqn = 0;
    _nagq = 0;
  }
}


//[[Rcpp::export]]
NumericMatrix foceiCalcCov(Environment e){
  std::string boundStr = "";
  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  // covType="analytic" stashes the fit's solve args (foceiSetup_); restoreFitSolve_
  // needs them for the whole cov step (incl. foceiCalcR below), so release only on
  // exit.  RAII covers every early return, the catch, and the covMethod=="" no-op.
  struct CovSolveArgsRelease { ~CovSolveArgsRelease() { releaseCovSolveArgs_(); } } _covSolveArgsRelease;
  try {
    if (op_focei.covMethod) {
      // Mu-referenced-FOCEI-family (muModel = lin/irls): the covariance must be
      // computed on the FULL corresponding focei/foce/focep model, NOT the mu->phi
      // reduced model used during estimation (the reduced parameterization gives
      // wrong SEs on the mu-referenced/linear parameters).  Bail here and recompute
      // the covariance at the R level with muModel="none" (.foceiRecomputeMuCov).
      if (op_focei.muModel != 0) {
        op_focei.cur = op_focei.totTick;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        e["covMethod"] = CharacterVector::create("");
        NumericMatrix ret;
        return ret;
      }
      op_focei.derivMethodSwitch=0;
      // Check boundaries
      unsigned int j, k;
      double cur;
      bool boundary=false;
      bool checkLowerBound=false;
      bool checkUpperBound=false;
      rx = getRxSolve_();
      if (op_focei.neta == 0) op_focei.covMethod = 2; // Always use hessian for NLS
      for (unsigned int k = op_focei.npars; k--;){
        if (R_FINITE(op_focei.lower[k])){
          op_focei.lower[k]=unscalePar(op_focei.lower,k);
        }
        if (R_FINITE(op_focei.upper[k])) {
          op_focei.upper[k]=unscalePar(op_focei.upper,k);
        }
      }
      if (op_focei.boundTol > 0){
        // Subtract nEstOmega so that Omega boundaries are not counted.
        for (k = op_focei.npars-op_focei.nEstOmega; k--;){
          if (op_focei.nbd[k] != 0){
            // bounds
            j=op_focei.fixedTrans[k];
            cur = op_focei.fullTheta[j];
            if (op_focei.nbd[k] == 1){
              // Lower only
              checkLowerBound = true;
            } else if (op_focei.nbd[k] == 2) {
              // Upper and lower
              checkLowerBound = true;
              checkUpperBound = true;
            } else {
              // Upper only
              checkUpperBound = true;
            }
            if (checkLowerBound &&
                j < thetaNames.size() &&
                (std::fabs((cur-op_focei.lower[k])/cur) < op_focei.boundTol)) {
              boundary = true;
              boundStr += "\"" + thetaNames[j] + "\" ";
            }
            if (checkUpperBound &&
                j < thetaNames.size() &&
                (std::fabs((op_focei.upper[k]-cur)/cur) < op_focei.boundTol)) {
              boundary = true;
              boundStr += "\"" + thetaNames[j] + "\" ";
            }
          }
        }
      }
      for (unsigned int j = getRxNsub(rx); j--;){
        focei_ind *fInd = &(inds_focei[j]);
        fInd->doChol=!(op_focei.cholSECov);
        fInd->doFD = 0;
      }
      op_focei.resetEtaSize = std::numeric_limits<double>::infinity(); // Dont reset ETAs
      op_focei.resetEtaSize=0; // Always reset ETAs.
      if (!e.exists("fullTheta")) {
        stop("focei environment requires 'fullTheta'");
      }
      NumericVector fullT = e["fullTheta"];
      NumericVector fullT2(op_focei.ntheta);
      std::copy(fullT.begin(), fullT.begin()+fullT2.size(), fullT2.begin());
      LogicalVector skipCov(op_focei.ntheta+op_focei.omegan);//skipCovN
      if (op_focei.skipCovN == 0){
        std::fill_n(skipCov.begin(), op_focei.ntheta, false);
        std::fill_n(skipCov.begin()+op_focei.ntheta, skipCov.size() - op_focei.ntheta, true);
      } else {
        std::copy(&op_focei.skipCov[0],&op_focei.skipCov[0]+op_focei.skipCovN,skipCov.begin());
        std::fill_n(skipCov.begin()+op_focei.skipCovN,skipCov.size()-op_focei.skipCovN,true);
      }
      e["skipCov"] = skipCov;
      // Unscaled objective and parameters.
      if (op_focei.scaleObjective){
        op_focei.scaleObjective=0;
        op_focei.lastOfv = op_focei.lastOfv * op_focei.initObjective / op_focei.scaleObjectiveTo;
      }
      // foceiSetupTheta_(op_focei.mvi, fullT2, skipCov, op_focei.scaleTo, false);
      setupAq0_(e);
      // muModel is always 0 here: ordinary methods never set it, and mu-referenced
      // (lin/irls) families bail above and recompute the covariance on the full model.
      foceiSetupTheta_(op_focei.mvi, fullT2, skipCov, 0, false);
      op_focei.scaleType=10;
      if (op_focei.covMethod && !boundary) {
        rx = getRxSolve_();
        op_focei.t0 = clock();
        op_focei.totTick=0;
        op_focei.cur=0;
        op_focei.curTick=0;
        RSprintf(_("calculating covariance matrix\n"));
        // Change options to covariance options
        // op_focei.scaleObjective = 0;
        op_focei.derivMethod = op_focei.covDerivMethod;

        arma::mat Rinv;
        op_focei.totTick=1;
        if (op_focei.covMethod == 1 || op_focei.covMethod == 3){
          op_focei.totTick+=op_focei.npars;
        }
        if (op_focei.covMethod == 1 || op_focei.covMethod == 2){
          op_focei.totTick += 2*op_focei.npars +2*(op_focei.npars*op_focei.npars);
        }
        op_focei.totTick += op_focei.npars;
        double hf, hphif, err;
        unsigned int j, k;
        arma::vec theta(op_focei.npars);
        for (k = op_focei.npars; k--;){
          j=op_focei.fixedTrans[k];
          theta[k] = op_focei.fullTheta[j];
        }
        std::copy(&theta[0], &theta[0] + op_focei.npars, &op_focei.theta[0]);
        arma::vec f0(1);
        f0(0) = op_focei.lastOfv;
        arma::vec grf(1);
        grf(0) = 0;

        arma::vec armaTheta(op_focei.npars);
        std::copy(&theta[0], &theta[0] + op_focei.npars, armaTheta.begin());
        double h = 0;
        int gillKcov;
        double gillStepCov;
        double gillFtolCov;
        double hessEps;
        for (int cpar = (int)op_focei.npars; cpar--;){
          if (op_focei.needOptimHess) {
            err = op_focei.rmatNormLlik ? 1/(std::fabs(theta[cpar])+1) : 1;
            gillKcov = op_focei.gillKcovLlik;
            gillStepCov=op_focei.gillStepCovLlik;
            gillFtolCov=op_focei.gillFtolCovLlik;
            hessEps = op_focei.hessEpsLlik;
          } else {
            err = op_focei.rmatNorm ? 1/(std::fabs(theta[cpar])+1) : 1;
            gillKcov = op_focei.gillKcov;
            gillStepCov=op_focei.gillStepCov;
            gillFtolCov=op_focei.gillFtolCov;
            hessEps = op_focei.hessEps;
          }
          if (op_focei.shi21maxOuter != 0) {
            op_focei.calcGrad=1;
            h = op_focei.aEps[cpar];
            op_focei.aEpsC[cpar] = shi21Central(shi21fnF, armaTheta, h,
                                                f0, grf, 0, cpar,
                                                op_focei.hessEpsInner, //double ef = 7e-7,
                                                1.5,  //double rl = 1.5,
                                                4.0,  //double ru = 6.0);;
                                                3.0, // nu
                                                op_focei.shi21maxOuter);  //maxiter=15
          } if (op_focei.gillKcov != 0){
            op_focei.gillRetC[cpar] = gill83(&hf, &hphif, &op_focei.gillDf[cpar], &op_focei.gillDf2[cpar], &op_focei.gillErr[cpar],
                                             &theta[0], cpar, hessEps, gillKcov, gillStepCov,
                                             gillFtolCov, -1, gill83fnG, 1, op_focei.lastOfv);
            // h=aEps*(|x|+1)/sqrt(1+fabs(f));
            // h*sqrt(1+fabs(f))/(|x|+1) = aEps
            // let err=2*sqrt(epsA/(1+f))
            // err*(aEps+|x|rEps) = h
            // Let aEps = rEps (could be a different ratio)
            // h/err = aEps(1+|x|)
            // aEps=h/err/(1+|x|)
            //
            op_focei.aEps[cpar]  = hf*err;
            op_focei.rEps[cpar]  = hf*err;
            if (op_focei.covGillF){
              op_focei.aEpsC[cpar] = hf*err;
              op_focei.rEpsC[cpar] = hf*err;
            } else {
              op_focei.aEpsC[cpar] = hphif*err;
              op_focei.rEpsC[cpar] = hphif*err;
            }
          } else {
            hf = hessEps;
            op_focei.aEps[cpar]  = hf*err;
            op_focei.rEps[cpar]  = hf*err;
            op_focei.aEpsC[cpar] = hf*err;
            op_focei.rEpsC[cpar] = hf*err;
          }
          op_focei.cur++;
          op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        }
        op_focei.didGill+=1;

        bool isPd;
        std::string rstr = "r";
        bool checkSandwich = false, checkSandwich2 = false;
        if (op_focei.covMethod == 1 || op_focei.covMethod == 2) {
          // R matrix based covariance
          arma::mat cholR;
          if (!e.exists("cholR")){
            foceiCalcR(e);
            // covType="analytic" signals an unrecoverable abort (the fit's global
            // solve could not be restored after the augmented sensitivity solves) by
            // zeroing covMethod.  Honor it as a clean cov-skip; do NOT fall through to
            // the S-matrix path below, which would run foceiS against a freed solve.
            if (op_focei.covMethod == 0) {
              warning(_("covariance step failed"));
              e["covMethod"] = CharacterVector::create("failed");
              NumericMatrix ret;
              return ret;
            }
          } else {
            op_focei.cur += op_focei.npars*2;
            op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          }
          if (e.exists("cholR")) {
            isPd = as<bool>(e["R.pd"]);
            if (!isPd){
              isPd = true;
              arma::vec E = as<arma::vec>(e["R.E"]);
              for (int j = E.size(); j--;){
                if (E[j] > op_focei.cholAccept){
                  isPd=false;
                  break;
                }
              }
              if (isPd){
                rstr = "r+";
                checkSandwich = true;
              }
            }
            if (!isPd){
              // Suggted by https://www.tandfonline.com/doi/pdf/10.1198/106186005X78800
              mat H0 = as<arma::mat>(e["R.0"]);
              H0 = H0*H0;
              cx_mat H1;
              bool success = sqrtmat(H1,H0);
              if (success){
                mat im = arma::imag(H1);
                mat re = arma::real(H1);
                if (!arma::any(arma::any(im,0))){
                  success= chol(H0, re);
                  if (success){
                    e["cholR"] = wrap(H0);
                    rstr = "|r|";
                    checkSandwich = true;
                    isPd = true;
                  }
                }
              }
            }
            op_focei.cur += op_focei.npars*2;
            op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            if (!isPd){
              warning(_("R matrix non-positive definite"));
              e["R"] = wrap(e["R.0"]);
              op_focei.covMethod = 3;
              op_focei.cur += op_focei.npars*2;
              op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            } else {
              cholR = as<arma::mat>(e["cholR"]);
              e["R"] = wrap(trans(cholR) * cholR);
              if (!e.exists("Rinv")){
                bool success  = inv(Rinv, trimatu(cholR));
                if (!success){
                  warning(_("Hessian (R) matrix seems singular; Using pseudo-inverse"));
                  Rinv = pinv(trimatu(cholR));
                  checkSandwich = true;
                }
                Rinv = Rinv * Rinv.t();
                e["Rinv"] = wrap(Rinv);
              } else {
                Rinv = as<arma::mat>(e["Rinv"]);
              }
              op_focei.cur++;
              op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, getRxCores(rx), op_focei.t0, 0);
              if (!e.exists("covR")){
                // Issue #666: the R-matrix covariance is the inverse observed
                // information Rinv = R^{-1} (R = 0.5*Hessian(-2LL) here).  This
                // was 2*Rinv, which made every covMethod="r" SE sqrt(2) too
                // large vs NONMEM $COV MATRIX=R and vs the (correct) sandwich
                // "r,s".  covR feeds both the final "r" output and the
                // sandwich-selection heuristic below, so it is fixed at source.
                e["covR"] = wrap(Rinv);
              }
              if (op_focei.covMethod == 2){
                e["cov"] = as<NumericMatrix>(e["covR"]);
              }
            }
          } else {
            RSprintf("\rR matrix calculation failed; Switch to S-matrix covariance.\n");
            op_focei.covMethod = 3;
            op_focei.cur += op_focei.npars*2;
            op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          }
        }
        arma::mat cholS;
        int origCov = op_focei.covMethod;
        std::string sstr="s";
        bool sHasZero = false;
        if (op_focei.covMethod == 1 || op_focei.covMethod == 3) {
          arma::vec theta(op_focei.npars);
          unsigned int j, k;
          for (k = op_focei.npars; k--;){
            j=op_focei.fixedTrans[k];
            theta[k] = op_focei.fullTheta[j];
          }
          if (!e.exists("cholS")){
            foceiS(&theta[0], e, sHasZero);
          } else {
            op_focei.cur += op_focei.npars;
            op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
          }
          if (e.exists("cholS")) {
            isPd = as<bool>(e["S.pd"]);
            if (!isPd){
              isPd=true;
              arma::vec E = as<arma::vec>(e["S.E"]);
              for (int j = E.size(); j--;){
                if (E[j] > op_focei.cholAccept){
                  isPd=false;
                  break;
                }
              }
              if (isPd){
                sstr="s+";
                checkSandwich = true;
              }
            }
            if (!isPd){
              // Suggted by https://www.tandfonline.com/doi/pdf/10.1198/106186005X78800
              mat H0 = as<arma::mat>(e["S0"]);
              H0 = H0*H0;
              cx_mat H1;
              bool success = sqrtmat(H1,H0);
              if (success){
                mat im = arma::imag(H1);
                mat re = arma::real(H1);
                if (!arma::any(arma::any(im,0))){
                  success= chol(H0,re);
                  if (success){
                    e["cholS"] = wrap(H0);
                    sstr = "|s|";
                    checkSandwich = true;
                    isPd = true;
                  }
                }
              }
            }
            if (!isPd){
              warning(_("S matrix non-positive definite"));
              if (op_focei.covMethod == 1){
                e["cov"] = as<NumericMatrix>(e["covR"]);
                op_focei.covMethod = 2;
              } else {
                warning(_("cannot calculate covariance"));
              }
              op_focei.cur += op_focei.npars*2;
              op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            } else {
              cholS = as<arma::mat>(e["cholS"]);
              arma::mat S;
              if (e.exists("S")){
                S = as<arma::mat>(e["S"]);
              } else {
                S = trans(cholS) * cholS;
                e["S"] = wrap(S);
              }
              if (op_focei.covMethod == 1){
                e["covRS"] = Rinv * S *Rinv;
                arma::mat covRS = as<arma::mat>(e["covRS"]);
                mat Sinv;
                bool success;
                success = inv(Sinv, trimatu(cholS));
                if (!success){
                  warning(_("S matrix seems singular; Using pseudo-inverse"));
                  Sinv = pinv(trimatu(cholS));
                }
                Sinv = Sinv * Sinv.t();
                // Issue #666: the S-matrix (OPG/cross-product) covariance is
                // Sinv = S^{-1}.  This was 4*Sinv, which made every covMethod="s"
                // SE 2x too large.  Confirmed against the empirical sampling
                // covariance of a known data-generating model: r, the sandwich,
                // and S^{-1} all match the true Cov(theta_hat), while 4*S^{-1}
                // is ~2x.  covS also feeds the selection heuristic below, so it
                // is fixed at source alongside covR.
                e["covS"]= Sinv;
                if (!checkSandwich) {
                  bool covRSsmall = arma::any(abs(covRS.diag()) < op_focei.covSmall);
                  if (covRSsmall) {
                    checkSandwich2 = true;
                  }
                }
                if (checkSandwich || checkSandwich2){
                  if (!checkSandwich2 && rstr == "r"){
                    // Use covR
                    e["cov"] = as<NumericMatrix>(e["covR"]);
                    op_focei.covMethod=2;
                  } else if (!checkSandwich2 && sstr == "s"){
                    // use covS
                    e["cov"] = as<NumericMatrix>(e["covS"]);
                    op_focei.covMethod=3;
                  } else {
                    // Now check sandwich matrix against R and S methods.
                    // issue #666 rescaled covR (2*Rinv->Rinv) and covS (4*Sinv->Sinv)
                    // but left the sandwich covRS unchanged.  This selector's covSmall
                    // diagonal floors and magnitude ordering were calibrated to the OLD
                    // covR/covS scale, so compare in that scale (covR*2, covS*4) to keep
                    // the estimator choice invariant to the rescale; the *installed*
                    // covR/covS/covRS (below) stay on the corrected #666 scale.
                    bool covRSsmall = arma::any(abs(covRS.diag()) < op_focei.covSmall);
                    double covRSd= sum(covRS.diag());
                    arma::mat covR = as<arma::mat>(e["covR"]);
                    bool covRsmall = arma::any(abs(2.0*covR.diag()) < op_focei.covSmall);
                    double covRd= sum(2.0*covR.diag());
                    arma::mat covS = as<arma::mat>(e["covS"]);
                    bool covSsmall = arma::any(abs(4.0*covS.diag()) < op_focei.covSmall);
                    double  covSd= sum(4.0*covS.diag());
                    if ((covRSsmall && covSsmall && covRsmall)){
                      e["cov"] = covRS;
                    } else if (covRSsmall && covSsmall && !covRsmall) {
                      e["cov"] = covR;
                      op_focei.covMethod=2;
                    } else if (covRSsmall && !covSsmall && covRsmall) {
                      e["cov"] = covS;
                      op_focei.covMethod=3;
                    } else if (covRSd > covRd){
                      // SE(RS) > SE(R)
                      if (covRd > covSd){
                        // SE(R) > SE(S)
                        e["cov"] = covS;
                        op_focei.covMethod=3;
                      } else {
                        e["cov"] = covR;
                        op_focei.covMethod=2;
                      }
                    } else if (covRSd > covSd){
                      e["cov"] = covS;
                      op_focei.covMethod=3;
                    } else {
                      e["cov"] = covRS;
                    }
                  }
                } else {
                  e["cov"] = covRS;
                }
              } else {
                mat Sinv;
                bool success;
                success = inv(Sinv, trimatu(cholS));
                if (!success){
                  warning(_("S matrix seems singular; Using pseudo-inverse"));
                  Sinv = pinv(trimatu(cholS));
                }
                Sinv = Sinv * Sinv.t();
                op_focei.cur++;
                op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
                e["cov"]= Sinv;   // issue #666: S-matrix covariance is S^{-1}, not 4*S^{-1}
              }
            }
          } else {
            if (op_focei.covMethod == 1){
              RSprintf("\rS matrix calculation failed; Switch to R-matrix covariance.\n");
              e["cov"] = wrap(e["covR"]);
              op_focei.covMethod = 2;
            } else {
              op_focei.covMethod=0;
              RSprintf("\rCould not calculate covariance matrix.\n");
              warning(_("cannot calculate covariance"));
              op_focei.cur++;
              op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
            }
          }
        }
        op_focei.cur=op_focei.totTick;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        if (e.exists("cov")){
          arma::mat cov = as<arma::mat>(e["cov"]);
          if (all(cov.diag() < 1e-7)){
            warning(_("The variance of all elements are unreasonably small, <1e-7"));
            op_focei.covMethod=0;
            e.remove("cov");
          }
        }
        if (op_focei.covMethod==0){
          warning(_("covariance step failed"));
          e["covMethod"] = CharacterVector::create("failed");
          NumericMatrix ret;
          return ret;
        } else {
          if (sHasZero) {
            warning(_("S matrix had problems solving for some subject and parameters"));
          }
          if (op_focei.covMethod == 1){
            bool doWarn=false;
            if (rstr == "|r|"){
              warning(_("R matrix non-positive definite but corrected by R = sqrtm(R%%*%%R)"));
              doWarn=true;
            } else if (rstr == "r+"){
              warning(_("R matrix non-positive definite but corrected (because of cholAccept)"));
              doWarn=true;
            }
            if (sstr == "|s|"){
              warning(_("S matrix non-positive definite but corrected by S = sqrtm(S%%*%%S)"));
              doWarn=true;
            } else if (sstr == "s+"){
              warning(_("S matrix non-positive definite but corrected (because of cholAccept)"));
              doWarn=true;
            }
            if (doWarn){
              warning(_("since sandwich matrix is corrected, you may compare to $covR or $covS if you wish"));
            }
            rstr =  rstr + "," + sstr;
            e["covMethod"] = wrap(rstr);
          } else if (op_focei.covMethod == 2){
            if (rstr == "|r|"){
              warning(_("R matrix non-positive definite but corrected by R = sqrtm(R%%*%%R)"));
            } else if (rstr == "r+"){
              warning(_("R matrix non-positive definite but corrected (because of cholAccept)"));
            }
            e["covMethod"] = wrap(rstr);
            if (origCov != 2){
              if (checkSandwich){
                warning(_("using R matrix to calculate covariance, can check sandwich or S matrix with $covRS and $covS"));
              } else {
                warning(_("using R matrix to calculate covariance"));
              }
            }
          } else if (op_focei.covMethod == 3){
            if (sHasZero) {
              warning(_("S matrix had problems solving for some subject and parameters"));
            }
            e["covMethod"] = wrap(sstr);
            if (origCov != 2){
              if (checkSandwich){
                warning(_("using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR"));
              } else {
                warning(_("using S matrix to calculate covariance"));
              }
            }
          }
          if (e.exists("cov")) {
            RObject covRO = e["cov"];
            if (covRO.sexp_type() == REALSXP &&
                Rf_isMatrix(covRO)) {
              return as<NumericMatrix>(covRO);
            }
          }
          NumericMatrix ret;
          return ret;
        }
      } else {
        if (op_focei.covMethod && boundary){
          warning(_("parameter estimate near boundary; covariance not calculated:\n   ") + boundStr +
                  _("\n use 'getVarCov' to calculate anyway"));
          e["covMethod"] = "Boundary issue; Get SEs with `getVarCov()`: " + boundStr;
        }
        op_focei.cur=op_focei.totTick;
        op_focei.curTick = par_progress(op_focei.cur, op_focei.totTick, op_focei.curTick, 1, op_focei.t0, 0);
        NumericMatrix ret;
        return ret;
      }
    }
  } catch (...) {
    warning(_("covariance step failed"));
    e["covMethod"] = CharacterVector::create("failed");
    NumericMatrix ret;
    return ret;
  }
  NumericMatrix ret;
  return ret;
}

// --- covType="fd", covFull=TRUE: full theta+sigma+Omega covariance by central finite
// differences of the objective, at the "normal fd seam" (the objective is the one
// foceiCalcR differences).  Structural + residual thetas are perturbed on
// op_focei.fullTheta (natural scale); the Omega block on the variance-covariance scale
// directly (op_focei.omegaInv/cholOmegaInv/logDetOmegaInv5 set from the perturbed Omega,
// op_focei.covFdDirect makes updateTheta skip the unscale + Cholesky rebuild) -- so the
// Hessian is natural-scale with no Jacobian (the rxOmegaVarCovDeriv parameterization).
// Split into small helpers below; the orchestrator is foceiCalcRFdFull. ---
struct FdFullCtx { IntegerVector thPos, omA, omB; arma::mat Om0; int nth = 0, nom = 0; };

// set op_focei Omega state from a variance-covariance Omega; false if it is not PD.
static bool foceiFdSetOmega(const arma::mat &Om) {
  arma::mat OmInv, ch;
  if (!arma::inv_sympd(OmInv, Om) || !arma::chol(ch, OmInv)) return false;
  double ld, sgn; arma::log_det(ld, sgn, Om);
  op_focei.omegaInv = OmInv; op_focei.cholOmegaInv = ch;
  op_focei.logDetOmegaInv5 = -0.5 * ld;       // 0.5*log|OmInv| = -0.5*log|Om|
  return true;
}

// objective at a full natural parameter vector x (thetas on fullTheta, Omega on the
// variance-covariance scale); NA_REAL if the perturbed Omega is not PD.
static double foceiFdObjAt(const FdFullCtx &c, const std::vector<double> &x) {
  for (int i = 0; i < c.nth; ++i) op_focei.fullTheta[c.thPos[i]] = x[i];
  arma::mat Om = c.Om0;
  for (int q = 0; q < c.nom; ++q) {
    Om(c.omA[q]-1, c.omB[q]-1) = x[c.nth+q]; Om(c.omB[q]-1, c.omA[q]-1) = x[c.nth+q];
  }
  if (!foceiFdSetOmega(Om)) return NA_REAL;
  return foceiOfv0(&op_focei.theta[0]);
}

// gill83 callback (plain function pointer -> file-static context): the objective at the
// perturbed natural parameter vector.  Same role as gill83fnF for the theta gradient.
static const FdFullCtx *g_fdGillCtx = nullptr;
static void foceiFdGill83fn(double *fp, double *theta, int, int) {
  const FdFullCtx &c = *g_fdGillCtx;
  std::vector<double> x(theta, theta + c.nth + c.nom);
  *fp = foceiFdObjAt(c, x);
}

// Per-parameter finite-difference step for coordinate i via the Gill-Murray-Saunders-
// Wright (1983) optimal-interval routine gill83 -- the SAME infrastructure foceiCalcR uses
// for the theta gradient/Hessian steps (it grows the step until the 2nd-difference
// condition error sits in [0.001, 0.1]).  Returns the accepted central step `hphif`, or
// NA_REAL if gill83 fails (caller falls back to a step-doubling search).
static double foceiFdGillStep(const FdFullCtx &c, int i, const std::vector<double> &x0, double f0) {
  g_fdGillCtx = &c;
  std::vector<double> xg = x0;                       // gill83 perturbs then restores this
  double hf, hphif, df, df2, ef;
  int gret = gill83(&hf, &hphif, &df, &df2, &ef, xg.data(), i,
                    op_focei.gillRtol, op_focei.gillK, op_focei.gillStep, op_focei.gillFtol,
                    -1, foceiFdGill83fn, 0, f0);
  return (gret == 1 && R_FINITE(hphif) && hphif > 0) ? hphif : NA_REAL;
}

// 5-point central diagonal 2nd-difference at step `e` (the foceiCalcR stencil).
static double foceiFdDiag5(const FdFullCtx &c, int i, const std::vector<double> &x0, double f0, double e) {
  std::vector<double> xp2=x0,xp1=x0,xm1=x0,xm2=x0;
  xp2[i]+=2*e; xp1[i]+=e; xm1[i]-=e; xm2[i]-=2*e;
  double f1=foceiFdObjAt(c,xp2), f2=foceiFdObjAt(c,xp1), f3=foceiFdObjAt(c,xm1), f4=foceiFdObjAt(c,xm2);
  if (!R_FINITE(f1)||!R_FINITE(f2)||!R_FINITE(f3)||!R_FINITE(f4)) return NA_REAL;
  return (-f1 + 16*f2 - 30*f0 + 16*f3 - f4) / (12*e*e);
}

// Gill-style adaptive diagonal 2nd-difference for coordinate i: grow the step until the
// 2nd-difference stabilizes (below that it is round-off dominated, which makes the
// off-diagonal Hessian -- and the covariance -- indefinite).  Returns it (or NA_REAL);
// writes the accepted step to *hOut for reuse on the off-diagonals.
static double foceiFdDiag(const FdFullCtx &c, int i, const std::vector<double> &x0,
                          double f0, double *hOut) {
  double base = std::max(std::fabs(x0[i]), 1e-3);
  double d2prev = NA_REAL, hi = base * 5e-4, d2Use = NA_REAL;
  *hOut = base * 1.6e-2;
  for (int s = 0; s < 6; ++s) {
    std::vector<double> xp = x0, xm = x0; xp[i] += hi; xm[i] -= hi;
    double fp = foceiFdObjAt(c, xp), fm = foceiFdObjAt(c, xm);
    if (R_FINITE(fp) && R_FINITE(fm)) {
      double d2 = (fp - 2*f0 + fm) / (hi*hi);
      *hOut = hi; d2Use = d2;
      if (R_FINITE(d2prev) && std::fabs(d2 - d2prev) <= 0.01 * std::fabs(d2)) break;
      d2prev = d2;
    }
    hi *= 2.0;
  }
  return d2Use;
}

// off-diagonal (i,j) 4-point mixed 2nd-difference with the accepted per-parameter steps h.
static double foceiFdOffDiag(const FdFullCtx &c, int i, int j, const std::vector<double> &x0,
                             const std::vector<double> &h) {
  std::vector<double> xpp=x0,xpm=x0,xmp=x0,xmm=x0;
  xpp[i]+=h[i]; xpp[j]+=h[j]; xpm[i]+=h[i]; xpm[j]-=h[j];
  xmp[i]-=h[i]; xmp[j]+=h[j]; xmm[i]-=h[i]; xmm[j]-=h[j];
  double a=foceiFdObjAt(c,xpp), b=foceiFdObjAt(c,xpm), cc=foceiFdObjAt(c,xmp), d=foceiFdObjAt(c,xmm);
  if (!R_FINITE(a)||!R_FINITE(b)||!R_FINITE(cc)||!R_FINITE(d)) return NA_REAL;
  return (a - b - cc + d) / (4*h[i]*h[j]);
}

// full natural-scale FD Hessian at x0 (f0=obj(x0)), mirroring foceiCalcR: a Gill-optimal
// per-parameter step (gill83) with the 5-point diagonal and 4-point off-diagonal stencils.
// If gill83 fails for a coordinate, fall back to the step-doubling diagonal.  Returns false
// on any non-finite probe.
static bool foceiFdHessian(const FdFullCtx &c, const std::vector<double> &x0, double f0, arma::mat &H) {
  int np = c.nth + c.nom;
  H.zeros(np, np);
  if (!R_FINITE(f0)) return false;
  std::vector<double> h(np);
  for (int i = 0; i < np; ++i) {
    double e = foceiFdGillStep(c, i, x0, f0);
    double d2 = R_FINITE(e) ? foceiFdDiag5(c, i, x0, f0, e) : NA_REAL;
    if (R_FINITE(d2)) { h[i] = e; H(i, i) = d2; }
    else              { H(i, i) = foceiFdDiag(c, i, x0, f0, &h[i]); }   // step-doubling fallback
    if (!R_FINITE(H(i, i))) return false;
  }
  for (int i = 0; i < np - 1; ++i) for (int j = i+1; j < np; ++j) {
    double v = foceiFdOffDiag(c, i, j, x0, h);
    if (!R_FINITE(v)) return false;
    H(i, j) = H(j, i) = v;
  }
  return true;
}

// enumeration/names from the R helper .foceiFdFullParams(e); false if nothing to do.
static bool foceiFdParams(Environment e, FdFullCtx &c, CharacterVector &nm) {
  Environment nlmixr2 = Environment::namespace_env("nlmixr2est");
  if (!nlmixr2.exists(".foceiFdFullParams")) return false;
  Function pf = as<Function>(nlmixr2[".foceiFdFullParams"]);
  List pl;
  try {
    RObject plo = pf(e);
    if (Rf_isNull(plo)) return false;
    pl = as<List>(plo);
  } catch (Rcpp::internal::InterruptedException&) {
    throw;
  } catch (Rcpp::LongjumpException&) {
    throw;
  } catch (...) { return false; }
  c.thPos = pl["thPos"]; c.omA = pl["omA"]; c.omB = pl["omB"]; nm = pl["names"];
  c.nth = c.thPos.size(); c.nom = c.omA.size();
  return (c.nth + c.nom) > 0;
}

// Orchestrator: assemble the FD Hessian around the fit, restore live state, install the
// natural cov solve(0.5*H) in e[".fdFullCov"] (installed as fit$cov by .foceiInstallFdFullCov).
void foceiCalcRFdFull(Environment e) {
  if (op_focei.neta <= 0) return;
  FdFullCtx c; CharacterVector nm;
  if (!foceiFdParams(e, c, nm)) return;
  c.Om0 = as<arma::mat>(getOmega());
  int np = c.nth + c.nom;

  // save live state (restored on exit)
  arma::mat omegaInv0 = op_focei.omegaInv, cholOmegaInv0 = op_focei.cholOmegaInv;
  double logDet0 = op_focei.logDetOmegaInv5;
  std::vector<double> fth0(op_focei.ntheta);
  std::copy(&op_focei.fullTheta[0], &op_focei.fullTheta[0] + op_focei.ntheta, fth0.begin());
  std::vector<double> x0(np);
  for (int i = 0; i < c.nth; ++i) x0[i] = op_focei.fullTheta[c.thPos[i]];
  for (int q = 0; q < c.nom; ++q) x0[c.nth + q] = c.Om0(c.omA[q]-1, c.omB[q]-1);

  op_focei.covFdDirect = 1;
  arma::mat H;
  bool ok = foceiFdHessian(c, x0, foceiFdObjAt(c, x0), H);

  // restore live state
  std::copy(fth0.begin(), fth0.end(), &op_focei.fullTheta[0]);
  op_focei.omegaInv = omegaInv0; op_focei.cholOmegaInv = cholOmegaInv0; op_focei.logDetOmegaInv5 = logDet0;
  op_focei.covFdDirect = 0;
  if (!ok) return;

  arma::mat cov;
  if (!arma::inv_sympd(cov, 0.5 * H) && !arma::inv(cov, 0.5 * H)) return;   // cov = 2 H^{-1}
  NumericMatrix covR = wrap(cov);
  covR.attr("dimnames") = List::create(nm, nm);
  e[".fdFullCov"] = covR;
}

void addLlikObs(Environment e) {
  if (op_focei.didLikCalc) {
    rx = getRxSolve_();
    NumericVector llikObs(getRxNall(rx));
    std::copy(&op_focei.llikObsFull[0], &op_focei.llikObsFull[0] + getRxNall(rx), llikObs.begin());
    e["llikObs"] = llikObs;
  }
}

void parHistData(Environment e, bool focei){
  if (!e.exists("method") && iterType.size() > 0) {
    CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
    CharacterVector dfNames;
    if (focei){
      CharacterVector dfNames2(3+op_focei.npars);
      dfNames = dfNames2;
    } else {
      CharacterVector dfNames2(3+thetaNames.size());
      dfNames = dfNames2;
    }
    dfNames[0] = "iter";
    dfNames[1] = "type";
    dfNames[2] = "objf";
    int i, j, k=1;
    if (focei){
      for (i = 0; i < op_focei.npars; i++){
        j=op_focei.fixedTrans[i];
        if (j < thetaNames.size()){
          dfNames[i+3] = thetaNames[j];
        } else {
          dfNames[i+3] = "o" + std::to_string(k++);
        }
      }
    } else {
      for (i = 0; i < thetaNames.size(); i++){
        dfNames[i+3] = thetaNames[i];
      }
    }
    // iter type parameters
    List ret;
    if (focei){
      ret = List(3+op_focei.npars);
    } else {
      ret = List(3+thetaNames.size());
    }
    int sz = niter.size()+niterGrad.size();
    IntegerVector tmp;
    std::vector<int> iter;
    iter.reserve(sz);
    iter.insert(iter.end(), niter.begin(), niter.end());
    iter.insert(iter.end(), niterGrad.begin(), niterGrad.end());
    ret[0] = iter;
    tmp = IntegerVector(sz);
    std::vector<int> typ;
    typ.reserve(sz);
    typ.insert(typ.end(), iterType.begin(), iterType.end());
    typ.insert(typ.end(), gradType.begin(), gradType.end());
    tmp = typ;
    tmp.attr("levels") = CharacterVector::create("Gill83 Gradient", "Mixed Gradient",
                                                 "Forward Difference", "Central Difference",
                                                 "Scaled", "Unscaled", "Back-Transformed",
                                                 "Forward Sensitivity", "Analytic Gradient");
    tmp.attr("class") = "factor";
    ret[1] = tmp;
    arma::mat cPar(vPar.size()/iterType.size(), iterType.size());
    std::copy(vPar.begin(), vPar.end(), cPar.begin());
    arma::mat vals;
    if (vGrad.size() > 0){
      arma::mat cGrad(vGrad.size()/gradType.size(), gradType.size());
      std::copy(vGrad.begin(), vGrad.end(), cGrad.begin());
      cPar = cPar.t();
      cGrad = cGrad.t();
      vals = arma::join_cols(cPar, cGrad);
    } else {
      cPar = cPar.t();
      vals = cPar;
    }
    if (focei){
      for (i = 0; i < min2(op_focei.npars+1, vals.n_cols); i++){
        ret[i+2]= vals.col(i);
      }
    } else {
      for (i = 0; i < thetaNames.size()+1; i++){
        ret[i+2]= vals.col(i);
      }
    }
    vGrad.clear();
    vPar.clear();
    iterType.clear();
    gradType.clear();
    niter.clear();
    niterGrad.clear();
    ret.attr("names")=dfNames;
    ret.attr("class") = "data.frame";
    ret.attr("row.names")=IntegerVector::create(NA_INTEGER, -sz);
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    Environment nlmixr2 = loadNamespace("nlmixr2est");
    Environment thetaReset = nlmixr2[".thetaReset"];
    if (thetaReset.exists("parHistData")) {
      // rbind  data.
      if (TYPEOF(thetaReset["parHistData"]) == VECSXP) {
        Function loadNamespace("loadNamespace", R_BaseNamespace);
        Environment nlmixr2 = loadNamespace("nlmixr2est");
        Function rbind = nlmixr2[".rbindParHistory"];
        ret = rbind(thetaReset["parHistData"], ret);
      }
      thetaReset.remove("parHistData");
    }
    e["parHistData"] = ret;
  }
}


void foceiFinalizeTables(Environment e){
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Function preFinalParTableHooksRun = nlmixr2[".preFinalParTableHooksRun"];
  preFinalParTableHooksRun(e);

  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  {
    // censInformation text; when censoring is present, note the 2nd-derivative treatment
    // (censOption: "laplace" = exact censored Hessian, "gauss" = historic Gauss-Newton) so
    // it is clear which was used.  Only the FOCEI-family conditional methods use censOption --
    // SAEM and the NLM family bow out (censOption is inert there), so their text stays plain.
    RObject ciR = censEstGetFactor();
    IntegerVector ci = as<IntegerVector>(ciR);
    CharacterVector lvls = as<CharacterVector>(ci.attr("levels"));
    std::string ciStr = as<std::string>(lvls[ci[0] - 1]);
    if (ci[0] > 1 && !op_focei.isSaem && !op_focei.isNlm)
      ciStr += (op_focei.censOption == 1 ? " (laplace)" : " (gauss)");
    e["censInformation"] = ciStr;
  }
  resetCensFlag();
  arma::mat cov;
  bool covExists = e.exists("cov");
  if (covExists) {
    if (rxode2::rxIs(e["cov"], "matrix")){
      cov= as<arma::mat>(e["cov"]);
    } else {
      covExists = false;
    }
  }
  LogicalVector skipCov = e["skipCov"];

  if (covExists) {
    Function getCor = nlmixr2[".cov2cor"];
    e["fullCor"] = getCor(e["cov"]);
    arma::mat cor = as<arma::mat>(e["fullCor"]);
    cor.diag().ones();
    arma::vec eigval;
    arma::mat eigvec;
    eig_sym(eigval, eigvec, cor);
    e["eigenCor"] = eigval;
    e["eigenVecCor"] = eigvec;
    unsigned int k=0;
    if (eigval.size() > 0){
      double mx=std::fabs(eigval[0]), mn, cur;
      mn=mx;
      for (k = eigval.size(); k--;){
        cur = std::fabs(eigval[k]);
        if (cur > mx){
          mx=cur;
        }
        if (cur < mn){
          mn=cur;
        }
      }
      e["conditionNumberCor"] = mx/mn;
    } else {
      e["conditionNumberCor"] = NA_REAL;
    }
  }

  if (covExists && op_focei.eigen) {
    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, cov);
    e["eigenCov"] = eigval;
    e["eigenVecCov"] = eigvec;
    unsigned int k=0;
    if (eigval.size() > 0){
      double mx=std::fabs(eigval[0]), mn, cur;
      mn=mx;
      for (k = eigval.size(); k--;){
        cur = std::fabs(eigval[k]);
        if (cur > mx){
          mx=cur;
        }
        if (cur < mn){
          mn=cur;
        }
      }
      e["conditionNumberCov"] = mx/mn;
    } else {
      e["conditionNumberCov"] = NA_REAL;
    }
  }
  arma::vec se1;
  if (covExists){
    se1 = sqrt(cov.diag());
  }
  DataFrame thetaDf = as<DataFrame>(e["theta"]);
  arma::vec theta = as<arma::vec>(thetaDf["theta"]);
  NumericVector se(theta.size());
  NumericVector cv(theta.size());
  std::fill_n(&se[0], theta.size(), NA_REAL);
  std::fill_n(&cv[0], theta.size(), NA_REAL);
  int j=0, k=0, l=0;
  if (covExists){
    for (int k = 0; k < se.size(); k++){
      if (k >= skipCov.size()) break;
      if (!skipCov[k]){
        se[k] = se1[j++];
        cv[k] = std::fabs(se[k]/theta[k])*100;
      }
    }
  }
  e["se"] = se;
  List popDf = List::create(_["Estimate"]=thetaDf["theta"], _["SE"]=se,
                            _["%RSE"]=cv);
  popDf.attr("class") = "data.frame";
  popDf.attr("row.names") = IntegerVector::create(NA_INTEGER,-theta.size());
  e["popDf"] = popDf;

  e["fixef"]=thetaDf["theta"];
  NumericMatrix tmpNM;
  List tmpL;
  int i;
  if (!Rf_isNull(e["etaObf"])) {
    List etas = e["etaObf"];
    IntegerVector idx = seq_len(etas.length())-1;
    etas = etas[idx != etas.length()-1];
    if (e.exists("idLvl")) {
      RObject idl = e["idLvl"];
      RObject eta0 = etas[0];
      if (idl.sexp_type() == STRSXP) {
        eta0.attr("class") = "factor";
        eta0.attr("levels") = idl;
      }
    }
    e["ranef"]=etas;

    // Now put names on the objects
    ////////////////////////////////////////////////////////////////////////////////
    // Eta Names
    CharacterVector etaNames=as<CharacterVector>(e["etaNames"]);
    tmpNM = getOmega();
    tmpNM.attr("dimnames") = List::create(etaNames, etaNames);
    e["omega"] = tmpNM;

    tmpL  = as<List>(e["ranef"]);
    List tmpL2 = as<List>(e["etaObf"]);
    CharacterVector tmpN  = tmpL.attr("names");
    CharacterVector tmpN2 = tmpL2.attr("names");
    // For mixture models the columns are: ID, MIXEST, ETA[1], ..., ETA[n]
    // so eta names start at offset 2; for non-mixture: offset 1 (ID only)
    {
      int offset = 1;
      for (int k = 0; k < tmpN2.size(); k++) {
        if (tmpN2[k] == "MIXEST") {
          offset = 2;
          break;
        }
      }
      for (i = 0; i < etaNames.size(); i++){
        if (i + offset <  tmpN.size())  tmpN[i+offset] = etaNames[i];
        if (i + offset < tmpN2.size()) tmpN2[i+offset] = etaNames[i];
      }
    }
    ////////////////////////////////////////////////////////////////////////////////
    tmpL.attr("names") = tmpN;
    tmpL2.attr("names") = tmpN2;
    e["ranef"] = tmpL;
    e["etaObf"] = tmpL2;
  } else {
    e["ranef"] = R_NilValue;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Theta names
  //
  tmpL = as<List>(e["theta"]);
  tmpL.attr("row.names") = thetaNames;
  e["theta"] = tmpL;

  tmpL=e["popDf"];
  // Add a few columns

  NumericVector Estimate = tmpL["Estimate"];
  NumericVector SE = tmpL["SE"];
  NumericVector RSE = tmpL["%RSE"];
  NumericVector EstBT(Estimate.size());
  NumericVector EstLower(Estimate.size());
  NumericVector EstUpper(Estimate.size());

  CharacterVector EstS(Estimate.size());
  CharacterVector SeS(Estimate.size());
  CharacterVector rseS(Estimate.size());
  CharacterVector btCi(Estimate.size());
  // LogicalVector EstBT(Estimate.size());
  // Rf_pt(stat[7],(double)n1,1,0)
  // Per-theta transform codes from R's .iterPrintXParFromUi: xform$xPar[i] /
  // xform$probitIdx[i] hold the log/logit/probit code; bounds vectors are
  // indexed by abs(code)-1.
  List xform = as<List>(e["xform"]);
  IntegerVector thetaXPar      = as<IntegerVector>(xform["xPar"]);
  IntegerVector thetaProbitIdx = as<IntegerVector>(xform["probitIdx"]);
  NumericVector logitThetaLow  = as<NumericVector>(xform["logitThetaLow"]);
  NumericVector logitThetaHi   = as<NumericVector>(xform["logitThetaHi"]);
  NumericVector probitThetaLow = as<NumericVector>(xform["probitThetaLow"]);
  NumericVector probitThetaHi  = as<NumericVector>(xform["probitThetaHi"]);
  const double *logitLowP  = logitThetaLow.size()  ? &logitThetaLow[0]  : NULL;
  const double *logitHiP   = logitThetaHi.size()   ? &logitThetaHi[0]   : NULL;
  const double *probitLowP = probitThetaLow.size() ? &probitThetaLow[0] : NULL;
  const double *probitHiP  = probitThetaHi.size()  ? &probitThetaHi[0]  : NULL;
  double qn = Rf_qnorm5(1.0-(1-op_focei.ci)/2, 0.0, 1.0, 1, 0);
  std::string cur;
  char buff[100];
  LogicalVector thetaFixed = thetaDf["fixed"];
  for (i = Estimate.size(); i--;) {
    int xpc = (i < thetaXPar.size())       ? thetaXPar[i]      : 0;
    int ppc = (i < thetaProbitIdx.size())  ? thetaProbitIdx[i] : 0;
    snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, Estimate[i]);
    EstS[i] = buff;
    EstBT[i] = scaleBackTransform(Estimate[i], xpc, ppc,
                                  logitLowP, logitHiP, probitLowP, probitHiP);
    snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstBT[i]);
    cur = buff;
    if (ISNA(SE[i])) {
      EstLower[i] = NA_REAL;
      EstUpper[i] = NA_REAL;
      if (thetaFixed[i]) {
        SeS[i]  = "FIXED";
        rseS[i] = "FIXED";
      } else {
        SeS[i]  = "";
        rseS[i] = "";
      }
    } else {
      EstLower[i] = scaleBackTransform(Estimate[i] - SE[i]*qn, xpc, ppc,
                                       logitLowP, logitHiP, probitLowP, probitHiP);
      EstUpper[i] = scaleBackTransform(Estimate[i] + SE[i]*qn, xpc, ppc,
                                       logitLowP, logitHiP, probitLowP, probitHiP);
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, SE[i]);
      SeS[i] = buff;
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, RSE[i]);
      rseS[i] = buff;
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstLower[i]);
      cur = cur + " (" + buff + ", ";
      snprintf(buff, sizeof(buff), "%.*g", (int)op_focei.sigdig, EstUpper[i]);
      cur = cur + buff + ")";
    }
    btCi[i] = cur;
  }
  tmpL["Back-transformed"] = EstBT;
  tmpL["CI Lower"] = EstLower;
  tmpL["CI Upper"] = EstUpper;
  tmpL.attr("row.names") = thetaNames;
  tmpL.attr("class") = "data.frame";
  e["popDf"]=tmpL;
  std::string bt = "Back-transformed(" + std::to_string((int)(op_focei.ci*100)) + "%CI)";

  List popDfSig;
  if (e.exists("cov") && rxode2::rxIs(e["cov"], "matrix")){
    popDfSig = List::create(_["Est."]=EstS,
                            _["SE"]=SeS,
                            _["%RSE"]=rseS,
                            _[bt]=btCi);
  } else {
    popDfSig = List::create(_["Est."]=EstS,
                            _["Back-transformed"] = btCi);
  }

  popDfSig.attr("row.names") = thetaNames;
  popDfSig.attr("class") = "data.frame";
  e["popDfSig"]=popDfSig;

  NumericVector tmpNV = e["fixef"];
  tmpNV.names() = thetaNames;
  e["fixef"] = tmpNV;


  tmpNV = e["se"];
  tmpNV.names() = thetaNames;
  e["se"] = tmpNV;

  // Now get covariance names
  if (e.exists("cov")  && rxode2::rxIs(e["cov"], "matrix")){
    tmpNM = as<NumericMatrix>(e["cov"]);
    CharacterVector thetaCovN(tmpNM.nrow());
    LogicalVector skipCov = e["skipCov"];
    int j=0;
    for (unsigned int k = 0; k < thetaNames.size(); k++){
      if (k >= skipCov.size()) break;
      if (j >= thetaCovN.size()) break;
      if (!skipCov[k]){
        thetaCovN[j++] = thetaNames[k];
      }
    }
    List thetaDim = List::create(thetaCovN,thetaCovN);
    tmpNM.attr("dimnames") = thetaDim;
    e["cov"]=tmpNM;
    if (e.exists("Rinv") && rxode2::rxIs(e["Rinv"], "matrix")){
      tmpNM = as<NumericMatrix>(e["Rinv"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["Rinv"]=tmpNM;
    }
    if (e.exists("Sinv") && rxode2::rxIs(e["Sinv"], "matrix")){
      tmpNM = as<NumericMatrix>(e["Sinv"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["Sinv"]=tmpNM;
    }
    if (e.exists("S") && rxode2::rxIs(e["S"], "matrix")){
      tmpNM = as<NumericMatrix>(e["S"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["S"]=tmpNM;
    }
    if (e.exists("R") && rxode2::rxIs(e["R"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R"]=tmpNM;
    }
    if (e.exists("covR") && rxode2::rxIs(e["covR"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covR"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covR"]=tmpNM;
    }
    if (e.exists("covRS") && rxode2::rxIs(e["covRS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covRS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covRS"]=tmpNM;
    }
    if (e.exists("covS") && rxode2::rxIs(e["covS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["covS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["covS"]=tmpNM;
    }
    if (e.exists("R.1") && rxode2::rxIs(e["R.1"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R.1"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R.1"]=tmpNM;
    }
    if (e.exists("R.2") && rxode2::rxIs(e["R.2"], "matrix")){
      tmpNM = as<NumericMatrix>(e["R.2"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["R.2"]=tmpNM;
    }
    if (e.exists("cholR") && rxode2::rxIs(e["cholR"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholR"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholR"]=tmpNM;
    }
    if (e.exists("cholR2") && rxode2::rxIs(e["cholR2"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholR2"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholR2"]=tmpNM;
    }
    if (e.exists("cholS") && rxode2::rxIs(e["cholS"], "matrix")){
      tmpNM = as<NumericMatrix>(e["cholS"]);
      tmpNM.attr("dimnames") = thetaDim;
      e["cholS"]=tmpNM;
    }
  }
  List objDf;
  if (e.exists("conditionNumberCov")){
    objDf = List::create(_["OBJF"] = as<double>(e["objective"]), _["AIC"]=as<double>(e["AIC"]),
                         _["BIC"] = as<double>(e["BIC"]), _["Log-likelihood"]=as<double>(e["logLik"]),
                         _["Condition#(Cov)"]=as<double>(e["conditionNumberCov"]),
                         _["Condition#(Cor)"]=as<double>(e["conditionNumberCor"]));
  } else {
    objDf = List::create(_["OBJF"] = as<double>(e["objective"]), _["AIC"]=as<double>(e["AIC"]),
                         _["BIC"] = as<double>(e["BIC"]), _["Log-likelihood"]=as<double>(e["logLik"]));
  }
  if (op_focei.neta == 0){
    if (op_focei.needOptimHess) {
      objDf.attr("row.names") = CharacterVector::create("lPop");
    } else {
      objDf.attr("row.names") = CharacterVector::create("Pop");
    }
    addLlikObs(e);
    e["ofvType"] = "Pop";
  } else if (op_focei.fo == 1){
    objDf.attr("row.names") = CharacterVector::create("FO");
    e["ofvType"] = "fo";
  } else if (op_focei.interaction){
    objDf.attr("row.names") = CharacterVector::create("FOCEi");
    addLlikObs(e);
    e["ofvType"] = "focei";
  } else if (e.exists("ofvType")) {
    std::string ofvType = as<std::string>(e["ofvType"]);
    objDf.attr("row.names") = ofvType;
    e["ofvType"]= ofvType;
  } else if (_nagq > 0)  {
    if (_nagq == 1) {
      objDf.attr("row.names") = CharacterVector::create("Laplace");
    } else {
      objDf.attr("row.names") = CharacterVector::create("AGQ" + std::to_string(_nagq));
    }
    addLlikObs(e);
    e["ofvType"] = "agq";
  } else {
    if (op_focei.needOptimHess) {
      objDf.attr("row.names") = CharacterVector::create("lFOCEi");
    } else {
      objDf.attr("row.names") = CharacterVector::create("FOCE");
    }
    addLlikObs(e);
    e["ofvType"] = "foce";
  }
  objDf.attr("class") = "data.frame";
  e["objDf"]=objDf;
  if (!e.exists("method")){
    if (op_focei.isImpmap) {
      // Importance-sampling EM; the objDf row stays "FOCEi" because the final
      // objective is a FOCEi evaluation at the EM estimates.
      e["method"] = op_focei.isImp ? "imp" : (op_focei.isQrpem ? "qrpem" : "impmap");
    } else if (_aqn > 0) {
      e["method"] ="AGQ";
    } else if (op_focei.neta == 0){
      e["method"] = "Population Only";
    } else if (op_focei.fo == 1){
      e["method"] = "FO";
    } else {
      e["method"] = "FOCE";
    }
  }
  if (!e.exists("extra")){
    if (op_focei.isImpmap) {
      // The EM drives itself (no outer optimizer, no gradients, and the mu
      // update is intrinsic to the M-step), so there is nothing to append.
      e["extra"] = "";
    } else {
      if (op_focei.neta == 0){
        e["extra"] = "";
        e["skipTable"] = NA_LOGICAL;
      } else if (op_focei.fo == 1){
        e["extra"] = "";
        e["skipTable"] = LogicalVector::create(true);
      } else if (op_focei.interaction || op_focei.needOptimHess){
        if(op_focei.scale.useColor){
          e["extra"] = "\033[31;1mi\033[0m";
        } else {
          e["extra"] = "i";
        }
      } else {
        e["extra"] = "";
      }
      List ctl = e["control"];
      // What the outer optimizer actually consumed: the analytic ("fast")
      // gradient, the finite-difference gradient, or a mix (per-iteration solve
      // fallbacks); plus the mu-referenced regression variant when active.
      std::string _details = as<std::string>(ctl["outerOptTxt"]);
      if (op_focei.fast && op_focei.maxOuterIterations > 0) {
        if (op_focei.nAnalyticGrad > 0 && op_focei.nFDGradFast == 0) {
          _details += "; grad: analytic";
        } else if (op_focei.nAnalyticGrad > 0) {
          _details += "; grad: analytic+fd";
        } else if (op_focei.nG > 0 || op_focei.nFDGradFast > 0) {
          _details += "; grad: fd";
        }
      }
      if (op_focei.muModel == 1) {
        _details += "; mu: lin";
      } else if (op_focei.muModel == 2) {
        _details += "; mu: irls";
      }
      if (_aqn == 0) {
        e["extra"] = as<std::string>(e["extra"]) +
          " (outer: " + _details +
          ")";
      } else if (_aqn == 1) {
        e["extra"] = as<std::string>(e["extra"]) +
          " (outer: " + _details +
          "; Laplace)";
      } else {
        e["extra"] = as<std::string>(e["extra"]) +
          " (outer: " + _details +
          "; nAGQ=" + std::to_string(_nagq)  + ")";
      }
    }
  }
  // rxode2::rxSolveFree();
  {
    int _nsub = (int)getRxNsub(rx);
    NumericVector _tf(_nsub);
    for (int _i = 0; _i < _nsub; _i++) {
      _tf[_i] = getIndTolFactor(getSolvingOptionsInd(rx, _i));
    }
    e["tolFactor"] = _tf;
  }
  e.attr("class") = "nlmixr2FitCore";
}

void setupAq1_(Environment e) {
  if (_aqn >= 1) {
    double *tmp = REAL(e["qx"]);
    std::copy(tmp, tmp + _aqn*op_focei.neta, op_focei.aqx);
    tmp = REAL(e["qw"]);
    std::copy(tmp, tmp + _aqn*op_focei.neta, op_focei.aqw);
    op_focei.aqfirst = as<bool>(e["qfirst"]);
    op_focei.aqLow = as<double>(e["aqLow"]);
    op_focei.aqHi = as<double>(e["aqHi"]);
  } else {
    op_focei.aqfirst = false;
  }
}

// ---- Importance-sampling EM interface (imp.h); reuses the FOCEI inner state.
// These run on the live op_focei / inds_focei set up by foceiSetup_ and are
// called from the impOuter() driver (src/imp.cpp).
int impNsub() {
  rx = getRxSolve_();
  return getRxNsub(rx);
}

int impNeta() {
  return op_focei.neta;
}

int impNsample() {
  return op_focei.impIsample;
}

double impGammaProp() {
  return op_focei.impGamma;
}

int impCores() {
  rx = getRxSolve_();
  return getOpCores(getSolvingOptions(rx));
}

// Override the solve core count (returns the previous value).  Used to force the
// mixture EM single-threaded: the expanded pseudo-subject inner solves are not
// thread-safe for a mixture, so parallel E-step / innerOpt passes race.
int impSetSolveCores(int cores) {
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  int old = op->cores;
  op->cores = cores;
  return old;
}

// 0.5 * log|Omega^-1| ( = -0.5 * log|Omega| ), the per-iteration normalizer that
// enters each subject's importance-sampling objective L_i.
double impLogDetOmegaInv5() {
  return op_focei.logDetOmegaInv5;
}

int impNiter() {
  return op_focei.impNiter;
}

std::string impDiagXform() {
  return op_focei.impDiagXform;
}

double impIaccept() { return op_focei.impIaccept; }
double impIscaleMin() { return op_focei.impIscaleMin; }
double impIscaleMax() { return op_focei.impIscaleMax; }
int impNconvWindow() { return op_focei.impNconvWindow; }

// Windowed-convergence tolerance on the (relative) objective change; derived
// from the table sigdig (10^-sigdig) when the control leaves it unset (<0).
double impCtol() {
  if (op_focei.impCtol >= 0) return op_focei.impCtol;
  return std::pow(10.0, -op_focei.sigdig);
}

int impMuGroupN() {
  return (int)op_focei.muGroupN;
}

// ---- Monte-Carlo covariance support (imp.cpp orchestrates the sampling + FD) ----

int impNtheta() { return (int)op_focei.ntheta; }

bool impCovEnabled() { return op_focei.impCov; }

// ---- quasi-random (QRPEM) + SIR controls -----------------------------------
bool impQrEnabled() { return op_focei.impQr; }
bool impQrShiftEnabled() { return op_focei.impQrShift; }
bool impQrRefreshEnabled() { return op_focei.impQrRefresh; }
bool impSirEnabled() { return op_focei.impSir; }
int impSirN() { return op_focei.impSirSample; }
int impBaseSeed() { return op_focei.impSeed; }

// ---- mixture (sub-population) support -------------------------------------
// impmap computes its OWN importance-sampling mixture posterior + proportion
// update (NONMEM-style, separate from FOCEI's Laplace fInd->mixProb).  It only
// READS op_focei.mixProb (the population proportions = the parameter) and reuses
// the per-component MAP modes solved for the expanded pseudo-subjects.

int impNmix() { return (int)op_focei.mixIdxN + 1; }        // number of components (1 if none)
double impMixProb(int j) { return op_focei.mixProb[j]; }    // population proportion of component j (0-based)

// Install absolute $MIX thetas (the multinomial-logit of the target proportions)
// and recompute the population proportions / Jacobian -- used for the stable
// mean-posterior EM update (mirrors the Omega-theta path in updateTheta()).
void impSetMixThetas(const arma::vec& theta) {
  for (int m = 0; m < (int)op_focei.mixIdxN; ++m)
    op_focei.fullTheta[op_focei.mixIdx[m] - 1] = theta[m];
  NumericVector curTheta(op_focei.ntheta);
  std::copy(&op_focei.fullTheta[0], &op_focei.fullTheta[0] + op_focei.ntheta, curTheta.begin());
  IntegerVector mixIdx(op_focei.mixIdxN);
  std::copy(&op_focei.mixIdx[0], &op_focei.mixIdx[0] + op_focei.mixIdxN, mixIdx.begin());
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment nlmixr2 = loadNamespace("nlmixr2est");
  Function f = as<Function>(nlmixr2[".getMixFromLog"]);
  NumericVector mp = f(curTheta, mixIdx);
  std::copy(mp.begin(), mp.end(), &op_focei.mixProb[0]);
  f = as<Function>(nlmixr2[".getMixJacFromLog"]);
  NumericVector mj = f(curTheta, mixIdx);
  std::copy(mj.begin(), mj.end(), &op_focei.mixProbGrad[0]);
}

// fullTheta indices of the estimated (non-fixed) thetas, in free-parameter order.
void impGetEstThetaIdx(std::vector<int>& idx) {
  idx.clear();
  for (unsigned int k = 0; k < op_focei.npars; ++k) {
    int j = op_focei.fixedTrans[k];
    if (j >= 0 && j < (int)op_focei.ntheta) idx.push_back(j);
  }
}

// fullTheta index of every free (estimated) parameter, in the optimizer's
// free-parameter (fixedTrans) order -- the order the fit's covariance uses.
// idx[k] < ntheta is a theta; idx[k] >= ntheta is the Omega parameter idx-ntheta.
void impGetCovParList(std::vector<int>& idx) {
  idx.clear();
  for (unsigned int k = 0; k < op_focei.npars; ++k) idx.push_back(op_focei.fixedTrans[k]);
}

double impGetFullThetaVal(int idx) { return op_focei.fullTheta[idx]; }

// Set theta fullTheta[idx] on every subject's parameter pointer (FD perturbation).
void impSetThetaAll(int idx, double val) {
  rx = getRxSolve_();
  op_focei.fullTheta[idx] = val;
  int nsub = getRxNsub(rx);
  for (int id = 0; id < nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    setIndParPtr(ind, op_focei.thetaTrans[idx], val);
  }
}

// Force likInner0 to re-solve subject id on its next call (theta changed but its
// eta may repeat a cached value, which would otherwise skip the re-solve).
void impForceResolve(int id) { inds_focei[id].setup = 0; }

// ---- Omega block of the MC covariance ----

// Number of parameterized Omega free parameters (fullTheta[ntheta .. +omegan-1]).
int impOmegaN() { return (int)op_focei.omegan; }

double impGetOmegaThetaVal(int m) { return op_focei.fullTheta[op_focei.ntheta + m]; }

// Set Omega free parameter m (FD perturbation) and rebuild omegaInv / cholOmegaInv
// / logDetOmegaInv5 from the full Omega-parameter vector (the same path
// updateTheta() takes), so likInner0's prior term and the objective normalizer
// both reflect it.
void impSetOmegaThetaAll(int m, double val) {
  op_focei.fullTheta[op_focei.ntheta + m] = val;
  NumericVector omegaTheta(op_focei.omegan);
  std::copy(&op_focei.fullTheta[0] + op_focei.ntheta,
            &op_focei.fullTheta[0] + op_focei.ntheta + op_focei.omegan,
            omegaTheta.begin());
  setOmegaTheta(omegaTheta);
  op_focei.omegaInv = getOmegaInv();
  op_focei.cholOmegaInv = getCholOmegaInv();
  op_focei.logDetOmegaInv5 = getOmegaDet();
}

// M-step helpers (EM loop lives in impOuter, src/imp.cpp).

// Overwrite subject id's eta (used to seed updateMuGroups with the conditional mean).
void impSetEta(int id, const arma::vec& eta) {
  focei_ind *fInd = &(inds_focei[id]);
  std::copy(eta.begin(), eta.begin() + op_focei.neta, &fInd->eta[0]);
}

// Current Omega matrix (its zero pattern gives the estimated-element structure).
void impGetOmega(arma::mat& Om) {
  Om = getOmegaMat();
}

// est="imp": no per-iteration MAP search (the E-step proposal is centered at the
// running conditional mean with covariance gamma * V_i instead of gamma * H_i^-1).
bool impIsImp() { return op_focei.isImp; }

// Read subject id's current eta (after updateMuGroups this is the recentered residual).
void impGetEta(int id, arma::vec& eta) {
  focei_ind *fInd = &(inds_focei[id]);
  eta.set_size(op_focei.neta);
  std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, eta.begin());
}

// Mu-referenced population-parameter update: OLS-regress the (conditional-mean)
// etas onto the mu-group covariates, updating fullTheta and recentering etas.
double impUpdateMuThetas() {
  return updateMuGroups();
}

// EM update for the simple mu-referenced intercepts (thetas that are the
// population mean of an eta, with no covariates -- these are not handled by the
// covariate regression in updateMuGroups()).  Operates on the current per-subject
// etas (which the caller has set to the conditional means): shifts each theta by
// the mean eta, recenters that eta to a mean-zero residual, and propagates both
// to every subject's solve.
void impMuInterceptStep() {
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  if (nsub == 0) return;
  IntegerVector &thIdx = op_focei.impMuThetaIdx;
  IntegerVector &etIdx = op_focei.impMuEtaIdx;
  for (int g = 0; g < thIdx.size(); ++g) {
    int th = thIdx[g], et = etIdx[g];
    double s = 0.0;
    for (int id = 0; id < nsub; ++id) s += inds_focei[id].eta[et];
    double delta = s / (double)nsub;
    op_focei.fullTheta[th] += delta;
    for (int id = 0; id < nsub; ++id) {
      rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
      inds_focei[id].eta[et] -= delta;
      setIndParPtr(ind, op_focei.thetaTrans[th], op_focei.fullTheta[th]);
      setIndParPtr(ind, op_focei.etaTrans[et], inds_focei[id].eta[et]);
    }
  }
}

int impThetaSensN() { return op_focei.impThetaSensIdx.size(); }

// 0-based eta indices whose Omega diagonal is fixed (the EM Omega update restores
// their rows/columns to the starting Omega so fix()ed variances are held).
void impGetOmegaFixedEta(std::vector<int>& idx) {
  idx.assign(op_focei.impOmegaFixedEta.begin(), op_focei.impOmegaFixedEta.end());
}

// ---- iteration print + parameter-history (shared scale.h machinery) ----
// The main setup already populated op_focei.scale (column names, back-transform
// codes, iterPrintControl cadence) for the FOCEI free parameters in fixedTrans
// order -- the same order impGetEstPar() returns.  Reconfigure it for impmap's
// EM walk, which is directly on the estimation scale: scaleTypeNone drops the
// redundant "Unscaled" rows (identical to the iteration rows) while the
// "Back-Transformed" rows still apply the per-theta transform.  The focei
// gradient-legend keyExtra is replaced (the EM has no gradients).  Records every
// iteration (save=1); prints at the control's cadence (scale.every, already set
// from iterPrintControl -- 0 when print is off).
static std::vector<double> _impScaleC;
void impIterPrintStart() {
  scaling *s = &op_focei.scale;
  int np = (int)op_focei.npars;
  _impScaleC.assign(std::max(np, 1), 1.0);
  s->scaleC = _impScaleC.data();
  s->scaleType = scaleTypeNone;
  s->scaleTo = 0.0; s->scaleCmin = 0.0; s->scaleCmax = 0.0;
  s->save = 1; s->simple = 0; s->showOfv = 1;
  s->keyExtra =
    "Omegas=chol(solve(omega));\n"
    "Diagonals are transformed, as specified by impmapControl(diagXform=)\n";
  s->vPar.clear(); s->niter.clear(); s->iterType.clear();
  s->vGrad.clear(); s->gradType.clear(); s->niterGrad.clear();
  s->cn = 0; s->printCount = 0;
  if (s->every != 0) scalePrintHeader(s);
}

void impIterPrintRow(arma::vec& par, double obj) {
  scalePrintFun(&op_focei.scale, par.memptr(), obj);
}

// Emit the closing rule (if printing) and stash the recorded walk on the fit
// environment as parHistData (iter/type/objf + one column per parameter);
// scaleParHisDf clears the accumulators.
void impIterPrintGet(Environment e) {
  scaling *s = &op_focei.scale;
  if (s->every != 0) scalePrintLine(s, min2((int)s->npars, s->ncol));
  RObject ph = scaleParHisDf(s);
  if (!ph.isNULL()) e["parHistData"] = ph;
}

// Newton step on the non-mu structural thetas: add step[s] to theta
// impThetaSensIdx[s] in fullTheta and propagate to every subject's solve.
void impUpdateStructThetas(const arma::vec& step) {
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  IntegerVector &thIdx = op_focei.impThetaSensIdx;
  for (int s = 0; s < thIdx.size(); ++s) {
    int th = thIdx[s];
    op_focei.fullTheta[th] += step[s];
    for (int id = 0; id < nsub; ++id) {
      rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
      setIndParPtr(ind, op_focei.thetaTrans[th], op_focei.fullTheta[th]);
    }
  }
}

// Re-optimize every subject's conditional mode at the current parameters.
void impReMap() {
  innerOpt();
}

// Install a new Omega: rebuild the rxSymInvChol environment (reusing the rxode2
// matrix->parameterization machinery), refresh the cached inverse/Cholesky/log-
// determinant, and copy the new Omega thetas into fullTheta so the next MAP and
// the output see them.
void impSetOmega(const arma::mat& Omega, const std::string& diagXform) {
  Environment rxode2ns = Environment::namespace_env("rxode2");
  Function f = rxode2ns["rxSymInvCholCreate"];
  _rxInv = as<List>(f(Rcpp::Named("mat") = wrap(Omega),
                      Rcpp::Named("diag.xform") = diagXform));
  if (op_focei.fo == 1) {
    op_focei.omega = getOmegaMat();
  } else {
    op_focei.omegaInv = getOmegaInv();
    op_focei.cholOmegaInv = getCholOmegaInv();
    op_focei.logDetOmegaInv5 = getOmegaDet();
  }
  NumericVector omegaTheta = getOmegaTheta();
  std::copy(omegaTheta.begin(),
            omegaTheta.begin() + op_focei.omegan,
            &op_focei.fullTheta[0] + op_focei.ntheta);
}

// Sync the optimizer's reference point (initPar) to the current natural
// fullTheta so the final foceiOuterFinal populates the fit at the converged
// estimates rather than the initial ones.
void impSyncInitParToFullTheta() {
  for (unsigned int k = 0; k < op_focei.npars; k++) {
    op_focei.initPar[k] = op_focei.fullTheta[op_focei.fixedTrans[k]];
  }
}

// Current estimated (free) parameter vector -- the estimated thetas plus the
// parameterized Omega elements -- for the EM convergence check.
void impGetEstPar(arma::vec& par) {
  par.set_size(op_focei.npars);
  for (unsigned int k = 0; k < op_focei.npars; k++) {
    par[k] = op_focei.fullTheta[op_focei.fixedTrans[k]];
  }
}

void impMapPass(Environment e) {
  // Single MAP pass at the initial parameters -- the same posthoc path
  // foceiOuter() takes when maxOuterIterations == 0.
  NumericVector x(op_focei.npars);
  for (unsigned int k = op_focei.npars; k--;) {
    x[k] = scalePar(op_focei.initPar, k);
  }
  foceiOuterFinal(x.begin(), e);
}

void impGetMode(int id, arma::vec& mode) {
  focei_ind *fInd = &(inds_focei[id]);
  mode.set_size(op_focei.neta);
  std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, mode.begin());
}

double impGetIndLik(int id) {
  focei_ind *fInd = &(inds_focei[id]);
  return fInd->lik[0];
}

bool impGetHessian(int id, arma::mat& H) {
  focei_ind *fInd = &(inds_focei[id]);
  int neta = op_focei.neta;
  // Establish the inner solve at this subject's mode before the FD Hessian.
  double f = likInner0(fInd->eta, id);
  if (ISNA(f)) {
    // The MAP already evaluated this mode, so an NA here means the cached solve
    // state is stale (e.g. left by a prior fit of a different-sized model).
    // Force a fresh solve and retry at the same mode before giving up.
    rx = getRxSolve_();
    setIndSolve(getSolvingOptionsInd(rx, getRxId(id)), -1);
    resetOpBadSolve(getSolvingOptions(rx));
    f = likInner0(fInd->eta, id);
    if (ISNA(f)) return false;
  }
  rx = getRxSolve_();
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
  H.set_size(neta, neta);
  H.zeros();
  arma::mat H0(neta, neta, arma::fill::zeros);
  return calcEtaHessian(fInd->eta, 0, id, fInd, ind, H, H0);
}

double impEvalJointLik(const arma::vec& eta, int id) {
  std::vector<double> ev(eta.begin(), eta.end());
  return likInner0(ev.data(), id);
}

// Central-FD fallback for the theta-sensitivity gradient (mirrors the pred-model
// finite-difference fallback the inner problem uses when the sensitivity ODE will
// not solve).  Solves the pred model at the current (theta, eta) for per-obs f and
// V, then shi21-optimized central differences over each estimated theta for
// d(f)/d(theta) and d(V)/d(theta).  Fills fvec/Vvec (nobs) and dfmat/dVmat
// (nobs x nSens); restores fullTheta on exit.
static void impThetaSensFD(int id, arma::vec& curTheta,
                           arma::vec& fvec, arma::vec& Vvec,
                           arma::mat& dfmat, arma::mat& dVmat) {
  int nSens = op_focei.impThetaSensIdx.size();
  int _rxId = getRxId(id);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, _rxId);
  focei_ind *fInd = &(inds_focei[id]);
  fvec = shi21ThetaGeneral(curTheta, id, 0); // f0 (also sets fInd->nObs)
  Vvec = shi21ThetaGeneral(curTheta, id, 1); // V0
  int nobs = (int)fInd->nObs;
  dfmat.set_size(nobs, nSens); dfmat.zeros();
  dVmat.set_size(nobs, nSens); dVmat.zeros();
  for (int s = 0; s < nSens; ++s) {
    int tIdx = op_focei.impThetaSensIdx[s];
    arma::vec grF(nobs), grV(nobs);
    double h = 0.0;
    shi21Central(shi21ThetaF, curTheta, h, fvec, grF, id, tIdx,
                 op_focei.hessEpsInner, 1.5, 4.5, 3.0, op_focei.shi21maxFD);
    h = 0.0;
    shi21Central(shi21ThetaR, curTheta, h, Vvec, grV, id, tIdx,
                 op_focei.hessEpsInner, 1.5, 4.5, 3.0, op_focei.shi21maxFD);
    dfmat.col(s) = grF;
    dVmat.col(s) = grV;
  }
  for (int t = 0; t < (int)op_focei.ntheta; ++t) {
    setIndParPtr(ind, op_focei.thetaTrans[t], op_focei.fullTheta[t]);
  }
}

// Accumulate subject id's importance-sampling contribution to the score `g`
// (length nSens) and Fisher-information Hessian `H` (nSens x nSens) for the
// estimated non-mu thetas, from its samples `S` (nsamp x neta) and normalized
// weights `zk`.  The theta-sensitivity model outputs per-observation f (rx_pred_),
// V (rx_r_), d(f)/d(theta) and d(V)/d(theta) in one solve; a bad ODE solve triggers
// the FOCEI-style tolerance-relaxation retry and, if it still fails, the pred-model
// central-FD fallback (impThetaSensFD).  For a Gaussian endpoint the score is
//   sum_j [ -(f-dv)/V d(f)/d(theta) + 0.5 ((dv-f)^2/V^2 - 1/V) d(V)/d(theta) ]
// and the (mean+variance) Fisher information is
//   sum_j [ d(f)/d(theta) d(f)/d(theta)'/V + 0.5 d(V)/d(theta) d(V)/d(theta)'/V^2 ].
void impThetaScore(int id, const arma::mat& S, const arma::vec& zk,
                   arma::vec& g, arma::mat& H) {
  int nSens = op_focei.impThetaSensIdx.size();
  if (nSens == 0 || op_focei.thetaSensOffset < 0 || op_focei.thetaSensDvOffset < 0 ||
      op_focei.thetaSensPredOffset < 0 || op_focei.thetaSensROffset < 0 ||
      rxThetaSens.calc_lhs == NULL) return;
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  int _rxId = getRxId(id);
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, _rxId);
  focei_ind *fInd = &(inds_focei[id]);
  int nsamp = S.n_rows;
  // Mixture: when called with an expanded pseudo-subject id, pin the theta-sens
  // solve to that subject's mixture component (like likInner0), so the sensitivity
  // reflects the component's structural parameters.
  if (op_focei.mixIdxN != 0) setIndMixest(ind, getRxMixFromId(id));
  // The shared solve pool is sized for the theta-sensitivity model (the largest
  // structure); the inner MAP runs with ind->neqOverride = innerNeq, so switch it
  // to thetaSensNeq for these solves and restore on exit (mirrors the predNeq FD
  // pattern).
  IndNeqOverrideGuard neqGuard(ind, op_focei.thetaSensNeq);
  // Sync the current thetas into this subject before solving (etas set per sample).
  arma::vec curTheta((unsigned int)op_focei.ntheta);
  for (int t = 0; t < (int)op_focei.ntheta; ++t) {
    curTheta[t] = op_focei.fullTheta[t];
    setIndParPtr(ind, op_focei.thetaTrans[t], op_focei.fullTheta[t]);
  }
  // per-observation DV and censoring info (same across samples), on the
  // transformed scale -- matching how the inner problem reads them.
  std::vector<double> dvv, limv;
  std::vector<int> censv;
  { int kk;
    for (int jj = 0; jj < getIndNallTimes(ind); ++jj) {
      setIndIdx(ind, jj); kk = getIndIx(ind, jj);
      if (getIndEvid(ind, kk) == 0) {
        dvv.push_back(tbs(getIndDv(ind, kk)));
        double lim = R_NegInf;
        if (hasRxLimit(rx)) {
          lim = getIndLimit(ind, kk);
          if (ISNA(lim)) lim = R_NegInf;
          else if (R_FINITE(lim)) lim = tbs(lim);
        }
        limv.push_back(lim);
        censv.push_back(hasRxCens(rx) ? getIndCens(ind, kk) : 0);
      }
    }
  }
  int nobs = (int)dvv.size();
  if (nobs == 0) return;
  arma::vec fvec(nobs), Vvec(nobs);
  arma::mat dfmat(nobs, nSens), dVmat(nobs, nSens);
  for (int k = 0; k < nsamp; ++k) {
    for (int j = 0; j < (int)op_focei.neta; ++j) {
      setIndParPtr(ind, op_focei.etaTrans[j], S(k, j));
    }
    // Solve the sensitivity model, relaxing ODE tolerances on a bad solve (as the
    // inner problem does), up to maxOdeRecalc; then fall back to central FD.
    if (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN) fInd->stickyRecalcN2 = 0;
    double prevTol = getIndTolFactor(ind);
    setIndSolve(ind, -1);
    resetOpBadSolve(op);
    thetaSensOde(_rxId);
    int jr = 0;
    while (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN &&
           indHasBadSolve(op, ind) && jr < op_focei.maxOdeRecalc) {
      fInd->stickyRecalcN2++;
      op_focei.reducedTol.store(1, std::memory_order_relaxed);
      op_focei.reducedTol2.store(1, std::memory_order_relaxed);
      atolRtolFactor_(op_focei.odeRecalcFactor);
      setIndSolve(ind, -1);
      resetOpBadSolve(op);
      thetaSensOde(_rxId);
      jr++;
    }
    if (jr != 0) {
      if (fInd->stickyRecalcN2 <= op_focei.stickyRecalcN) setIndTolFactor(ind, prevTol);
      else op_focei.stickyTol.store(1, std::memory_order_relaxed);
    }
    if (!indHasBadSolve(op, ind)) {
      // read f, V, d(f)/d(theta), d(V)/d(theta) from the sensitivity model lhs
      int ko = 0, kk;
      double curT;
      for (int jj = 0; jj < getIndNallTimes(ind) && ko < nobs; ++jj) {
        setIndIdx(ind, jj); kk = getIndIx(ind, jj); curT = getTime(kk, ind);
        double *lhs = getIndLhs(ind);
        if (isDose(getIndEvid(ind, kk))) {
          rxThetaSens.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, jj), lhs);
          continue;
        } else if (getIndEvid(ind, kk) == 0) {
          rxThetaSens.calc_lhs(_rxId, curT, getOpIndSolve(op, ind, jj), lhs);
          fvec[ko] = lhs[op_focei.thetaSensPredOffset];
          Vvec[ko] = lhs[op_focei.thetaSensROffset];
          for (int s = 0; s < nSens; ++s) {
            dfmat(ko, s) = lhs[op_focei.thetaSensOffset + s];
            dVmat(ko, s) = lhs[op_focei.thetaSensDvOffset + s];
          }
          ko++;
        }
      }
    } else {
      // sensitivity ODE unsolvable: central FD on the pred model.  Flag it as an
      // optimization problem (same warning as the inner FD fallback).
      op_focei.didPredSolve.store(true, std::memory_order_relaxed);
      impThetaSensFD(id, curTheta, fvec, Vvec, dfmat, dVmat);
      // restore this sample's eta (impThetaSensFD only touched thetas)
      for (int j = 0; j < (int)op_focei.neta; ++j) {
        setIndParPtr(ind, op_focei.etaTrans[j], S(k, j));
      }
      if ((int)dfmat.n_rows != nobs) continue;
    }
    for (int jo = 0; jo < nobs; ++jo) {
      double f = fvec[jo], V = Vvec[jo];
      if (!R_finite(V) || V <= 0.0 || !R_finite(f)) continue;
      arma::rowvec df = dfmat.row(jo), dV = dVmat.row(jo);
      if (!df.is_finite() || !dV.is_finite()) continue;
      bool isCens = (censv[jo] != 0) || (R_FINITE(limv[jo]) && !ISNA(limv[jo]));
      if (isCens) {
        // Exact censored score/information from the analytic partials of
        // rho = -logLik (out[0]=rho_f, out[1]=rho_r, out[2..4]=rho_ff/rho_fr/rho_rr),
        // so the M-step gradient for BLQ/M2/M3/M4 points is correct (no FD).  The
        // log-likelihood score is -rho_f, -rho_r; the information is the rho 2nd
        // derivatives.  A normal obs reduces to the Gauss-Newton form below.
        double cp[9]; for (int _i = 0; _i < 9; ++_i) cp[_i] = 0.0;
        censNormalPartials((double)censv[jo], dvv[jo], limv[jo], f, V, 2, cp);
        g += zk[k] * (-cp[0] * df.t() - cp[1] * dV.t());
        H += zk[k] * (cp[2] * (df.t() * df) +
                      cp[3] * (df.t() * dV + dV.t() * df) +
                      cp[4] * (dV.t() * dV));
      } else {
        double err = f - dvv[jo];
        g += zk[k] * ((-err / V) * df.t() +
                      (0.5 * (err * err / (V * V) - 1.0 / V)) * dV.t());
        H += zk[k] * ((df.t() * df) / V + 0.5 * (dV.t() * dV) / (V * V));
      }
    }
  }
}

// est="impmap": the theta-sensitivity model is the largest structure, so it sizes
// the single shared solve pool (foceiSetup_'s rxSolve_ uses _impPoolModel when set,
// mirroring how the pool is sized for the inner model in plain FOCEI).  The inner
// MAP then runs with ind->neqOverride = op_focei.innerNeq (like the predNeq FD
// path), and the theta-sensitivity solve overrides to thetaSensNeq.
// (_impPoolModel is declared earlier, before foceiSetup_.)

// True when the solve pool is sized for the larger theta-sensitivity model
// (multi-endpoint pool-sizing).  The inner solves then run with ind->neqOverride
// against a pool sized differently from op->neq, whose per-thread work arrays are
// not thread-safe under this override -- so the EM must run serial (like the
// mixture path), else the parallel E-step non-deterministically rejects a
// subject's importance samples (neff collapse).
bool impPoolSizing() { return op_focei.innerNeq > 0; }

// Pin every subject's effective inner state count to op_focei.innerNeq, since the
// pool (and op->neq) is sized for the larger theta-sensitivity model.
void impSetInnerNeqOverride() {
  if (op_focei.innerNeq <= 0) return;
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  for (int id = 0; id < nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    setIndNeqOverride(ind, op_focei.innerNeq);
  }
}

// Clear the persistent inner neqOverride at the end of the fit so it does not
// bleed into a later fit whose op->neq differs (the shared rx_global is reused).
void impClearInnerNeqOverride() {
  if (op_focei.innerNeq <= 0) return;
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  for (int id = 0; id < nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    setIndNeqOverride(ind, -1);
  }
}

//' Fit/Evaluate FOCEi
//'
//' This shouldn't be called directly.
//'
//' @param e Environment
//'
//' @return A focei fit object
//'
//' @keywords internal
//' @export
//[[Rcpp::export]]
Environment foceiFitCpp_(Environment e){
  focei_wall_clock::time_point wallT0 = focei_wall_clock::now();
  List model = e["model"];
  bool doPredOnly = false;
  op_focei.canDoFD = false;
  if (e.exists("nEstOmega")) {
    op_focei.nEstOmega = as<int>(e["nEstOmega"]);
  } else {
    op_focei.nEstOmega = 0;
  }
  setupAq0_(e);
  if (model.containsElementNamed("inner")) {
    RObject inner;
    if (model.containsElementNamed("innerLlik")) {
      inner = model["innerLlik"];
    } else {
      inner = model["inner"];
    }
    if (rxode2::rxIs(inner, "rxode2")) {
      RObject _dataSav = as<RObject>(e["dataSav"]);
      NumericVector _thetaIni = as<NumericVector>(e["thetaIni"]);
      IntegerVector _mixIdx = as<IntegerVector>(e["mixIdx"]);
      Nullable<LogicalVector> _thetaFixed =  as<Nullable<LogicalVector>>(e["thetaFixed"]);
      Nullable<LogicalVector> _skipCov = as<Nullable<LogicalVector>>(e["skipCov"]);
      RObject _rxInv = e["rxInv"];
      Nullable<NumericVector> _lower = as<Nullable<NumericVector>> (e["lower"]);
      Nullable<NumericVector> _upper = as<Nullable<NumericVector>> (e["upper"]);
      Nullable<NumericMatrix> _etaMat = as<Nullable<NumericMatrix>>(e["etaMat"]);
      Nullable<List> _control = as<Nullable<List>>(e["control"]);
      setupAq0_(e);
      // est="impmap": if the theta-sensitivity model has more ODE states than the
      // inner model, use it to size the shared solve pool (foceiSetup_'s rxSolve_)
      // and record the inner state count so inner solves run with neqOverride.
      _impPoolModel = R_NilValue;
      op_focei.innerNeq = 0;
      if (model.containsElementNamed("thetaSens")) {
        RObject _tsm = model["thetaSens"];
        if (rxode2::rxIs(_tsm, "rxode2")) {
          int _tsNeq = as<CharacterVector>(rxode2::rxModelVars_(_tsm)["state"]).size();
          int _innNeq = as<CharacterVector>(rxode2::rxModelVars_(inner)["state"]).size();
          if (_tsNeq > _innNeq) { _impPoolModel = _tsm; op_focei.innerNeq = _innNeq; }
        }
      }
      foceiSetup_(inner, _dataSav, _thetaIni, _mixIdx, _thetaFixed, _skipCov,
                  _rxInv, _lower, _upper, _etaMat, _control);
      _impPoolModel = R_NilValue; // consumed by foceiSetup_
      if (model.containsElementNamed("predNoLhs")) {
        RObject noLhs;
        if (model.containsElementNamed("predNoLhsLlik")) {
          noLhs = model["predNoLhsLlik"];
        } else {
          noLhs = model["predNoLhs"];
        }
        if (rxode2::rxIs(noLhs, "rxode2")) {
          List mvp = rxode2::rxModelVars_(noLhs);
          rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);
          op_focei.canDoFD = true;
          // Locate rx_pred_ in the predNoLhs lhs (AR(1) lag defs may precede it).
          op_focei.predNoLhsOffset = 0;
          CharacterVector predLhs = as<CharacterVector>(mvp["lhs"]);
          for (int il = 0; il < predLhs.size(); ++il) {
            if (as<std::string>(predLhs[il]) == "rx_pred_") { op_focei.predNoLhsOffset = il; break; }
          }
        } else {
          stop(_("focei cannot be run without rxode2 'predNoLhs'"));
        }
      } else {
        stop(_("focei cannot be run without 'predNoLhs'"));
      }
      // est="impmap": load the theta-sensitivity model (peer of rxInner/rxPred)
      // and record its ODE state count + the lhs offset of the first
      // d(f)/d(theta) output rx__sens_rx_pred__BY_THETA_1___.
      op_focei.thetaSensOffset = -1;
      op_focei.thetaSensNeq = 0;
      if ((op_focei.isImpmap || op_focei.isAdvi) && model.containsElementNamed("thetaSens")) {
        RObject ts = model["thetaSens"];
        if (rxode2::rxIs(ts, "rxode2")) {
          List mvts = rxode2::rxModelVars_(ts);
          rxUpdateFuns(as<SEXP>(mvts["trans"]), &rxThetaSens);
          op_focei.thetaSensNeq = as<CharacterVector>(mvts["state"]).size();
          rxThetaSens.neq = op_focei.thetaSensNeq;
          // The model outputs rx__sens_rx_pred__BY_THETA_j___ (d(f)/d(theta)) and
          // rx__sens_rx_r__BY_THETA_j___ (d(V)/d(theta)) for each estimated non-mu
          // theta j (1-based), in two ascending, contiguous blocks.  Record the lhs
          // offsets of the first of each (impThetaSensIdx[0] + 1).
          CharacterVector tsLhs = as<CharacterVector>(mvts["lhs"]);
          op_focei.thetaSensPredOffset = -1;
          op_focei.thetaSensROffset = -1;
          if (op_focei.impThetaSensIdx.size() > 0) {
            std::string j0 = std::to_string(op_focei.impThetaSensIdx[0] + 1);
            std::string firstF = "rx__sens_rx_pred__BY_THETA_" + j0 + "___";
            std::string firstV = "rx__sens_rx_r__BY_THETA_" + j0 + "___";
            for (int il = 0; il < tsLhs.size(); ++il) {
              std::string nm = as<std::string>(tsLhs[il]);
              if (nm == firstF) op_focei.thetaSensOffset = il;
              else if (nm == firstV) op_focei.thetaSensDvOffset = il;
              else if (nm == "rx_pred_") op_focei.thetaSensPredOffset = il;
              else if (nm == "rx_r_") op_focei.thetaSensROffset = il;
            }
          }
        }
        // The pool is sized for the theta-sensitivity model (the largest); pin the
        // inner solves to the inner state count via ind->neqOverride.
        if (op_focei.thetaSensOffset >= 0 && op_focei.innerNeq > 0) {
          impSetInnerNeqOverride();
        }
      }
      // Now setup which ETAs need a finite difference
      if (model.containsElementNamed("eventEta")) {
        IntegerVector eventEta = model["eventEta"];
        std::copy(eventEta.begin(), eventEta.end(),&op_focei.etaFD[0]);
      }
    } else if (model.containsElementNamed("predOnly")){
      if (model.containsElementNamed("predOnlyLlik")){
        inner = model["predOnlyLlik"];
      } else {
        inner = model["predOnly"];
      }
      if (rxode2::rxIs(inner, "rxode2")){
        RObject _dataSav = as<RObject>(e["dataSav"]);
        NumericVector _thetaIni = as<NumericVector>(e["thetaIni"]);
        IntegerVector _mixIdx = as<IntegerVector>(e["mixIdx"]);
        Nullable<LogicalVector> _thetaFixed =  as<Nullable<LogicalVector>>(e["thetaFixed"]);
        Nullable<LogicalVector> _skipCov = as<Nullable<LogicalVector>>(e["skipCov"]);
        RObject _rxInv = e["rxInv"];
        Nullable<NumericVector> _lower = as<Nullable<NumericVector>> (e["lower"]);
        Nullable<NumericVector> _upper = as<Nullable<NumericVector>> (e["upper"]);
        Nullable<NumericMatrix> _etaMat = as<Nullable<NumericMatrix>>(e["etaMat"]);
        Nullable<List> _control = as<Nullable<List>>(e["control"]);
        setupAq0_(e);
        foceiSetup_(inner, _dataSav, _thetaIni,  _mixIdx, _thetaFixed, _skipCov,
                    _rxInv, _lower, _upper, _etaMat, _control);
        doPredOnly = true;
        if (op_focei.neta == 0) doPredOnly = false;
      } else {
        stop(_("cannot run this function"));
      }
    } else {
      doPredOnly=true;
      setupAq0_(e);
      foceiSetupTrans_(as<CharacterVector>(e[".params"]));
      RObject _dataSav = as<RObject>(e["dataSav"]);
      NumericVector _thetaIni = as<NumericVector>(e["thetaIni"]);
      IntegerVector _mixIdx = as<IntegerVector>(e["mixIdx"]);
      Nullable<LogicalVector> _thetaFixed =  as<Nullable<LogicalVector>>(e["thetaFixed"]);
      Nullable<LogicalVector> _skipCov = as<Nullable<LogicalVector>>(e["skipCov"]);
      RObject _rxInv = e["rxInv"];
      Nullable<NumericVector> _lower = as<Nullable<NumericVector>> (e["lower"]);
      Nullable<NumericVector> _upper = as<Nullable<NumericVector>> (e["upper"]);
      Nullable<NumericMatrix> _etaMat = as<Nullable<NumericMatrix>>(e["etaMat"]);
      Nullable<List> _control = as<Nullable<List>>(e["control"]);
      setupAq0_(e);
      foceiSetup_(inner, _dataSav, _thetaIni,  _mixIdx, _thetaFixed, _skipCov,
                  _rxInv, _lower, _upper, _etaMat, _control);
    }
  } else {
    doPredOnly=true;
    if (!model.containsElementNamed("predOnly")) {
      stop(_("with focei inner, 'model$predOnly' needs to be present in the environment"));
    }
    RObject inner;
    if (model.containsElementNamed("predOnlyLlik")) {
      inner = model["predOnlyLlik"];
    } else {
      inner = model["predOnly"];
    }
    // foceiSetupTrans_(as<CharacterVector>(e[".params"]));
    if (!e.exists("dataSav")) {
      stop(_("without focei inner setup, this needs a data frame 'dataSav' in the environment"));
    }
    if (!e.exists("thetaFixed")) {
      stop(_("without focei inner setup, this needs a logical vector 'thetaFixed' in the environment"));
    }
    if (!e.exists("skipCov")) {
      stop(_("without focei inner setup, this needs a logical vector 'skipCov' in the environment"));
    }
    if (!e.exists("rxInv")) {
      stop(_("without focei inner setup, this needs a data frame 'rxInv' in the environment"));
    }
    if (!e.exists("lower")) {
      stop(_("without focei inner setup, this needs a numeric vector 'lower' in the environment"));
    }
    if (!e.exists("upper")) {
      stop(_("without focei inner setup, this needs a numeric vector 'upper' in the environment"));
    }
    if (!e.exists("upper")) {
      stop(_("without focei inner setup, this needs a numeric vector 'etaMat' in the environment"));
    }
    if (!e.exists("control")) {
      stop(_("without focei inner setup, this needs a foceiControl in the 'control' in the environment"));
    }
    if (!e.exists("thetaNames")) {
      stop(_("without focei inner setup, this needs a character vector in the 'thetaNames' in the environment"));
    }
    if (rxode2::rxIs(inner, "rxode2")){
      RObject _dataSav = as<RObject>(e["dataSav"]);
      NumericVector _thetaIni = as<NumericVector>(e["thetaIni"]);
      IntegerVector _mixIdx = as<IntegerVector>(e["mixIdx"]);
      Nullable<LogicalVector> _thetaFixed =  as<Nullable<LogicalVector>>(e["thetaFixed"]);
      Nullable<LogicalVector> _skipCov = as<Nullable<LogicalVector>>(e["skipCov"]);
      RObject _rxInv = e["rxInv"];
      Nullable<NumericVector> _lower = as<Nullable<NumericVector>> (e["lower"]);
      Nullable<NumericVector> _upper = as<Nullable<NumericVector>> (e["upper"]);
      Nullable<NumericMatrix> _etaMat = as<Nullable<NumericMatrix>>(e["etaMat"]);
      Nullable<List> _control = as<Nullable<List>>(e["control"]);
      setupAq0_(e);
      foceiSetup_(inner, _dataSav, _thetaIni,  _mixIdx, _thetaFixed, _skipCov,
                  _rxInv, _lower, _upper, _etaMat, _control);
      doPredOnly = true;
      if (op_focei.neta == 0) doPredOnly = false;
    } else {
      stop(_("model$predOnly needs to be an rxode2 object"));
    }
  }
  setupAq1_(e);
  if (e.exists("setupTime")){
    e["setupTime"] = as<double>(e["setupTime"]) + foceiElapsedSeconds(wallT0);
  } else {
    e["setupTime"] = foceiElapsedSeconds(wallT0);
  }
  wallT0 = focei_wall_clock::now();
  CharacterVector thetaNames=as<CharacterVector>(e["thetaNames"]);
  // Capture mu-group theta display names here while the full thetaNames is
  // available (op_focei.scale's own is npars-sized and excludes them); used
  // only by printMuGroupThetaRow() below.
  if (op_focei.muModel != 0 && op_focei.muGroupN > 0) {
    CharacterVector muThetaNm(op_focei.muGroupN);
    for (unsigned int g = 0; g < op_focei.muGroupN; g++) {
      int jj = op_focei.muGroupTheta[g];
      muThetaNm[g] = (jj < thetaNames.size()) ? thetaNames[jj] : NA_STRING;
    }
    op_focei.muGroupThetaNames = muThetaNm;
    if (op_focei.muGroupCovN > 0) {
      CharacterVector muCovNm(op_focei.muGroupCovN);
      for (unsigned int c = 0; c < op_focei.muGroupCovN; c++) {
        int jj = op_focei.muGroupCovTheta[c];
        muCovNm[c] = (jj < thetaNames.size()) ? thetaNames[jj] : NA_STRING;
      }
      op_focei.muGroupCovThetaNames = muCovNm;
    }
  }
  IntegerVector xType = e["xType"];
  std::fill_n(&op_focei.scaleC[0], op_focei.ntheta+op_focei.omegan, NA_REAL);
  if (e.exists("scaleC")){
    arma::vec scaleC = as<arma::vec>(e["scaleC"]);
    std::copy(scaleC.begin(), scaleC.end(), &op_focei.scaleC[0]);
  }
  // Theta transforms (ntheta-indexed, from R's .iterPrintXParFromUi xform list,
  // length ntheta_total); re-indexed below via fixedTrans into the npars-sized
  // xPar/probitIdxArr.
  List xform = as<List>(e["xform"]);
  IntegerVector thetaXPar      = as<IntegerVector>(xform["xPar"]);
  IntegerVector thetaProbitIdx = as<IntegerVector>(xform["probitIdx"]);
  op_focei.logitThetaLow  = as<NumericVector>(xform["logitThetaLow"]);
  op_focei.logitThetaHi   = as<NumericVector>(xform["logitThetaHi"]);
  op_focei.probitThetaLow = as<NumericVector>(xform["probitThetaLow"]);
  op_focei.probitThetaHi  = as<NumericVector>(xform["probitThetaHi"]);
  // Theta slots (j < ntheta) get their code from thetaXPar[j]/probitIdx[j];
  // omega slots (j >= ntheta) use the existing xType scaling code (xPar 2-5),
  // and never get probit.
  for (unsigned int k = op_focei.npars; k--;){
    int j = op_focei.fixedTrans[k];
    op_focei.xPar[k] = 0;
    op_focei.probitIdxArr[k] = 0;
    if ((int)op_focei.ntheta < j){
      op_focei.xPar[k] = xType[j-op_focei.ntheta];
    } else {
      if (j < thetaXPar.size())       op_focei.xPar[k]         = thetaXPar[j];
      if (j < thetaProbitIdx.size())  op_focei.probitIdxArr[k] = thetaProbitIdx[j];
    }
  }
  std::string tmpS;
  if (op_focei.nF2) {
    Function loadNamespace("loadNamespace", R_BaseNamespace);
    Environment nlmixr2 = loadNamespace("nlmixr2est");
    Environment thetaReset = nlmixr2[".thetaReset"];
    restoreFromEnvrionment(thetaReset);
  }
  // Populate the shared scale.h struct unconditionally (even when print is
  // off) so scalePrintFun calls below can no-op via scale.every==0. scale.save=0
  // since focei keeps its own iteration history in module-level globals.
  {
    // Build a CharacterVector with the printed column names in the same
    // order scalePrintFun emits them (fixedTrans-indexed).
    CharacterVector scaleNames(op_focei.npars);
    int k = 1;
    for (unsigned int i = 0; i < op_focei.npars; i++) {
      int jj = op_focei.fixedTrans[i];
      if (jj < thetaNames.size()) {
        scaleNames[i] = thetaNames[jj];
      } else {
        scaleNames[i] = "o" + std::to_string(k++);
      }
    }
    op_focei.scale.npars         = op_focei.npars;
    op_focei.scale.initPar       = op_focei.initPar;
    op_focei.scale.scaleC        = op_focei.scaleC;
    op_focei.scale.xPar          = op_focei.xPar;
    op_focei.scale.logitThetaLow = op_focei.logitThetaLow.size() ? &op_focei.logitThetaLow[0] : NULL;
    op_focei.scale.logitThetaHi  = op_focei.logitThetaHi.size()  ? &op_focei.logitThetaHi[0]  : NULL;
    op_focei.scale.probitIdx       = op_focei.probitIdxArr;
    op_focei.scale.probitThetaLow  = op_focei.probitThetaLow.size() ? &op_focei.probitThetaLow[0] : NULL;
    op_focei.scale.probitThetaHi   = op_focei.probitThetaHi.size()  ? &op_focei.probitThetaHi[0]  : NULL;
    op_focei.scale.thetaNames    = scaleNames;
    op_focei.scale.normType      = op_focei.normType;
    op_focei.scale.scaleType     = op_focei.scaleType;
    op_focei.scale.scaleCmin     = op_focei.scaleCmin;
    op_focei.scale.scaleCmax     = op_focei.scaleCmax;
    op_focei.scale.scaleTo       = op_focei.scaleTo;
    op_focei.scale.c1            = op_focei.c1;
    op_focei.scale.c2            = op_focei.c2;
    op_focei.scale.simple        = 0;
    // focei has a meaningful per-iteration objective function, so show the
    // "Function Val." column.  op_focei is a zero-initialized global and
    // focei configures its scaling struct by hand here (it does not call
    // scaleSetup(), which is where showOfv would otherwise default to 1),
    // so this must be set explicitly — otherwise the column, its header,
    // and the gradient-row label that lives in the same slot are dropped
    // for every outer optimizer (e.g. outerOpt="bobyqa").
    op_focei.scale.showOfv       = 1;
    // focei's richer Key suffix — gradient-method legend and omega note.
    // Appended after "X: Back-transformed parameters; " by scalePrintHeader.
    op_focei.scale.keyExtra =
      "G: Gill difference gradient approximation\n"
      "F: Forward difference gradient approximation\n"
      "C: Central difference gradient approximation\n"
      "M: Mixed forward and central difference gradient approximation\n"
      "A: Analytic (forward sensitivity) gradient (fast=TRUE)\n"
      "Unscaled parameters for Omegas=chol(solve(omega));\n"
      "Diagonals are transformed, as specified by foceiControl(diagXform=)\n";
    op_focei.scale.printCount    = 0;
    op_focei.scale.save          = 0;  // focei records into module-level globals
    op_focei.scale.cn            = 0;
    op_focei.scale.showOfv       = 1;
  }
  if (op_focei.maxOuterIterations > 0 && op_focei.printTop == 1){
    op_focei.t0 = clock();
    scalePrintHeader(&op_focei.scale);
  }
  if (doPredOnly){
    if (e.exists("objective")){
      nlmixr2EnvSetup(e, as<double>(e["objective"]));
    } else {
      stop(_("not setup right, needs objective in the environment"));
    }
  } else {
    op_focei.didHessianReset=0;
    op_focei.didEtaNudge =0;
    op_focei.didEtaReset=0;
    op_focei.stickyRecalcN2=0;
    op_focei.stickyRecalcN1=0;
    if (op_focei.isImpmap) {
      // est="impmap": run the importance-sampling EM driver (src/imp.cpp) in
      // place of the FOCEI outer optimizer.  Module M1 does a single MAP pass.
      impOuter(e);
    } else {
      foceiOuter(e);
    }
    if (op_focei.didHessianReset==1){
      warning(_("Hessian reset during optimization; (can control by foceiControl(resetHessianAndEta=.))"));
    }
    if (op_focei.didEtaNudge==1){
      warning(_("initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=))"));
    }
    if (op_focei.didEtaReset==1){
      warning(_("ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.))"));
    }
    if (op_focei.repeatGillN > 0){
      warning(_("tolerances were reduced during Gill Gradient, so it was repeated %d/%d times\nYou can control this with foceiControl(repeatGillMax=.)"), op_focei.repeatGillN, op_focei.repeatGillMax);
    }
    if (op_focei.neta != 0 && op_focei.maxOuterIterations > 0  && R_FINITE(op_focei.resetThetaFinalSize)) {
      focei_options *fop = &op_focei;
      op_focei.etaM.zeros();
      op_focei.etaS.zeros();
      double n = 1.0;
      for (int id=getRxNsub(rx); id--;){
        focei_ind *fInd = &(inds_focei[id]);
        mat etaMat(fop->neta, 1);
        std::copy(&fInd->eta[0], &fInd->eta[0] + op_focei.neta, etaMat.begin());
        mat oldM = op_focei.etaM;
        op_focei.etaM = op_focei.etaM + (etaMat - op_focei.etaM)/n;
        op_focei.etaS = op_focei.etaS + (etaMat - op_focei.etaM) %  (etaMat - oldM);
        n += 1.0;
      }
      op_focei.eta1SD = 1/sqrt(op_focei.etaS);
      thetaReset(op_focei.resetThetaFinalSize);
    }
  }
  e["optimTime"] = foceiElapsedSeconds(wallT0);
  NumericVector scaleSave(op_focei.ntheta+op_focei.omegan);
  for (unsigned int i =op_focei.ntheta+op_focei.omegan;i--;){
    scaleSave[i] = getScaleC(i);
  }
  e["scaleC"] = scaleSave;
  parHistData(e, true); // Need to calculate before the parameter translations are mangled
  thetaResetObj(e);
  IntegerVector gillRet(op_focei.ntheta+op_focei.omegan);
  NumericVector gillAEps(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillREps(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillAEpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillREpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillCAEpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  NumericVector gillCREpsC(op_focei.ntheta+op_focei.omegan,NA_REAL);
  bool warnGill = false;
  int j = op_focei.npars;
  for (int i = op_focei.ntheta+op_focei.omegan; i--;){
    gillRet[i] = op_focei.gillRet[i]+1;
    if (gillRet[i] != 1) {
      gillAEps[i] = op_focei.aEps[--j];
      gillREps[i] = op_focei.rEps[j];
      gillAEpsC[i] = op_focei.aEpsC[j];
      gillREpsC[i] = op_focei.rEpsC[j];
    }
    if (gillRet[i] >= 3) warnGill=true;
  }
  CharacterVector gillLvl = CharacterVector::create("Not Assessed","Good","High Grad Error", "Constant Grad","Odd/Linear Grad",
                                                    "Grad changes quickly");
  gillRet.attr("levels") = gillLvl;
  gillRet.attr("class") = "factor";
  e["gillRet"] = gillRet;
  wallT0 = focei_wall_clock::now();
  {
    // covSolveTol tightens the finite-difference cov solves (R/S + full-cov FD)
    CovSolveTolGuard _covTolGuard(e);
    foceiCalcCov(e);
    // covType="fd" + covFull=TRUE: the full theta+sigma+Omega FD covariance (installed
    // by .foceiInstallFdFullCov).  covType="analytic" fills the full cov its own way.
    if (op_focei.covFull && op_focei.covMethod != 0 && e.exists("control")) {
      List _ctlF = as<List>(e["control"]);
      std::string _covTypeF = _ctlF.containsElementNamed("covType") ?
        as<std::string>(_ctlF["covType"]) : "fd";
      if (_covTypeF != "analytic") {
        try { foceiCalcRFdFull(e); }
        catch (Rcpp::internal::InterruptedException&) { throw; }
        catch (Rcpp::LongjumpException&) { throw; }
        catch (...) {}
      }
    }
  }
  if (op_focei.didPredSolve) {
    warning(_("numerical difficulties solving forward sensitivity inner problem, tried approximating with more inaccurate numeric differences"));
  }
  IntegerVector gillRetC(op_focei.ntheta+op_focei.omegan);
  bool warnGillC = false;
  j = op_focei.npars;
  for (int i = op_focei.ntheta+op_focei.omegan; i--;){
    gillRetC[i] = op_focei.gillRetC[i]+1;
    if (gillRetC[i] >= 3) warnGillC=true;
    if (gillRetC[i] != 1) {
      gillCAEpsC[i] = op_focei.aEpsC[--j];
      gillCREpsC[i] = op_focei.rEpsC[j];
    }
  }
  gillRetC.attr("levels") = gillLvl;
  gillRetC.attr("class") = "factor";
  e["gillRetC"] = gillRetC;
  e["covTime"] = foceiElapsedSeconds(wallT0);
  List timeDf = List::create(_["setup"]=as<double>(e["setupTime"]),
                             _["optimize"]=as<double>(e["optimTime"]),
                             _["covariance"]=as<double>(e["covTime"]));
  timeDf.attr("class") = "data.frame";
  timeDf.attr("row.names") = CharacterVector::create("elapsed");
  e["time"] = timeDf;
  List scaleInfo = List::create(as<NumericVector>(e["fullTheta"]),
                                as<NumericVector>(e["scaleC"]), gillRet,
                                gillAEps,
                                gillREps,
                                gillAEpsC,
                                gillREpsC,
                                gillRetC,
                                gillCAEpsC,
                                gillCREpsC);
  scaleInfo.attr("names") = CharacterVector::create("est","scaleC","Initial Gradient",
                                                    "Forward aEps","Forward rEps",
                                                    "Central aEps","Central rEps",
                                                    "Covariance Gradient",
                                                    "Covariance aEps","Covariance rEps");
  scaleInfo.attr("class") = "data.frame";
  scaleInfo.attr("row.names") = IntegerVector::create(NA_INTEGER,-gillRet.size());
  e["scaleInfo"] = scaleInfo;
  if (warnGillC && warnGill){
    warning(_("gradient problems with initial estimate and covariance; see $scaleInfo"));
  } else if (warnGill){
    warning(_("gradient problems with initial estimate; see $scaleInfo"));
  } else if (warnGillC){
    warning(_("gradient problems with covariance; see $scaleInfo"));
  }
  if (op_focei.reducedTol){
    if (op_focei.stickyTol){
      warning(_("tolerances (atol/rtol) were increased (after %d bad solves) for some difficult ODE solving during the optimization.\ncan control with foceiControl(stickyRecalcN=)\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"), op_focei.stickyRecalcN);
    } else {
      warning(_("tolerances (atol/rtol) were temporarily increased for some difficult ODE solving during the optimization.\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"));
    }
  }
  if (op_focei.zeroGrad){
    warning(_("zero gradient replaced with small number (%f)"), sqrt(DBL_EPSILON));
  }
  foceiFinalizeTables(e);
  // NumericVector scaleC(op_focei.ntheta+op_focei.omegan);
  // std::copy(&op_focei.scaleC[0], &op_focei.scaleC[0]+op_focei.ntheta+op_focei.omegan, scaleC.begin());
  // e["scaleC"]= scaleC;
  if (op_focei.maxOuterIterations){
    RSprintf(_("done\n"));
  }
  return e;
}

// ---------------------------------------------------------------------------
// VAE inner-likelihood interface: set up the FOCEi inner problem ONCE (no outer
// optimizer), then evaluate the inner objective/gradient per subject (and per
// mixture component) directly through likInner0/lpInner in a parallel OpenMP
// loop -- reusing the inner logic (mixtures via setIndMixest, multiple
// endpoints, error structures, log-likelihood, censoring) without the nlmixr2 R
// interface. `e` is the same env foceiFitCpp_ consumes (built by .vaeInnerSetup
// in R). updateTheta() populates op_focei's omega inverse + theta state that the
// inner solve needs (normally done inside the outer objective evaluation).
static NumericVector _vaeInitPar;

//[[Rcpp::export]]
RObject vaeInnerSetup_(Environment e) {
  op_focei.canDoFD = false;
  op_focei.nEstOmega = e.exists("nEstOmega") ? as<int>(e["nEstOmega"]) : 0;
  setupAq0_(e);
  List model = e["model"];
  RObject inner = model.containsElementNamed("innerLlik") ? model["innerLlik"] : model["inner"];
  if (!rxode2::rxIs(inner, "rxode2")) stop("vaeInnerSetup_ needs an rxode2 inner model");
  RObject _dataSav = as<RObject>(e["dataSav"]);
  NumericVector _thetaIni = as<NumericVector>(e["thetaIni"]);
  IntegerVector _mixIdx = as<IntegerVector>(e["mixIdx"]);
  Nullable<LogicalVector> _thetaFixed = as<Nullable<LogicalVector>>(e["thetaFixed"]);
  Nullable<LogicalVector> _skipCov = as<Nullable<LogicalVector>>(e["skipCov"]);
  RObject _rxInv = e["rxInv"];
  Nullable<NumericVector> _lower = as<Nullable<NumericVector>>(e["lower"]);
  Nullable<NumericVector> _upper = as<Nullable<NumericVector>>(e["upper"]);
  Nullable<NumericMatrix> _etaMat = as<Nullable<NumericMatrix>>(e["etaMat"]);
  Nullable<List> _control = as<Nullable<List>>(e["control"]);
  // est="advi": if the theta-sensitivity model has more ODE states than the inner
  // model, size the shared solve pool to it and record the inner state count so
  // inner solves run with neqOverride (mirrors the impmap full-fit path).
  _impPoolModel = R_NilValue;
  op_focei.innerNeq = 0;
  if (model.containsElementNamed("thetaSens")) {
    RObject _tsm = model["thetaSens"];
    if (rxode2::rxIs(_tsm, "rxode2")) {
      int _tsNeq = as<CharacterVector>(rxode2::rxModelVars_(_tsm)["state"]).size();
      int _innNeq = as<CharacterVector>(rxode2::rxModelVars_(inner)["state"]).size();
      if (_tsNeq > _innNeq) { _impPoolModel = _tsm; op_focei.innerNeq = _innNeq; }
    }
  }
  NumericVector initPar = foceiSetup_(inner, _dataSav, _thetaIni, _mixIdx, _thetaFixed, _skipCov,
                                      _rxInv, _lower, _upper, _etaMat, _control);
  _impPoolModel = R_NilValue; // consumed by foceiSetup_
  // predNoLhs -> finite-difference event sensitivities
  if (model.containsElementNamed("predNoLhs")) {
    RObject noLhs = model.containsElementNamed("predNoLhsLlik") ? model["predNoLhsLlik"] : model["predNoLhs"];
    if (rxode2::rxIs(noLhs, "rxode2")) {
      List mvp = rxode2::rxModelVars_(noLhs);
      rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxPred);
      op_focei.canDoFD = true;
      op_focei.predNoLhsOffset = 0;
      CharacterVector predLhs = as<CharacterVector>(mvp["lhs"]);
      for (int il = 0; il < predLhs.size(); ++il)
        if (as<std::string>(predLhs[il]) == "rx_pred_") { op_focei.predNoLhsOffset = il; break; }
    }
  }
  if (model.containsElementNamed("eventEta")) {
    IntegerVector eventEta = model["eventEta"];
    std::copy(eventEta.begin(), eventEta.end(), &op_focei.etaFD[0]);
  }
  // est="advi": wire the theta-sensitivity model offsets (rx_pred_, rx_r_, and
  // the first d(f)/d(theta) / d(V)/d(theta) columns) so impThetaScore can supply
  // the outer population gradient -- mirrors the impmap full-fit path (which sets
  // these in foceiFitCpp_, a code path vaeInnerSetup_ does not go through).
  op_focei.thetaSensOffset = -1;
  op_focei.thetaSensNeq = 0;
  if ((op_focei.isImpmap || op_focei.isAdvi) && model.containsElementNamed("thetaSens")) {
    RObject ts = model["thetaSens"];
    if (rxode2::rxIs(ts, "rxode2")) {
      List mvts = rxode2::rxModelVars_(ts);
      rxUpdateFuns(as<SEXP>(mvts["trans"]), &rxThetaSens);
      op_focei.thetaSensNeq = as<CharacterVector>(mvts["state"]).size();
      rxThetaSens.neq = op_focei.thetaSensNeq;
      CharacterVector tsLhs = as<CharacterVector>(mvts["lhs"]);
      op_focei.thetaSensPredOffset = -1;
      op_focei.thetaSensROffset = -1;
      if (op_focei.impThetaSensIdx.size() > 0) {
        std::string j0 = std::to_string(op_focei.impThetaSensIdx[0] + 1);
        std::string firstF = "rx__sens_rx_pred__BY_THETA_" + j0 + "___";
        std::string firstV = "rx__sens_rx_r__BY_THETA_" + j0 + "___";
        for (int il = 0; il < tsLhs.size(); ++il) {
          std::string nm = as<std::string>(tsLhs[il]);
          if (nm == firstF) op_focei.thetaSensOffset = il;
          else if (nm == firstV) op_focei.thetaSensDvOffset = il;
          else if (nm == "rx_pred_") op_focei.thetaSensPredOffset = il;
          else if (nm == "rx_r_") op_focei.thetaSensROffset = il;
        }
      }
    }
    if (op_focei.thetaSensOffset >= 0 && op_focei.innerNeq > 0) {
      impSetInnerNeqOverride();
    }
  }
  // populate the omega inverse + theta-dependent state for the inner solve
  _vaeInitPar = clone(initPar);
  updateTheta(&_vaeInitPar[0]);
  return initPar;
}

// Re-parameterize the already-set-up VAE inner problem at new natural-scale
// theta/omega values WITHOUT re-running foceiSetup_ -- the per-gradient-step
// fast path.  Each non-fixed parameter goes through scalePar() into the
// reduced par vector and updateTheta() rebuilds fullTheta, the per-id solve
// parameters and the omega inverse (exactly what focei's outer objective does
// per evaluation).  The omega block requires the setup's diagonal rxInv with
// diag.xform="sqrt", whose parameter for a diagonal omega is omega_kk^(-1/4)
// (the square root of the chol(Omega^-1) diagonal).
static void vaeInnerUpdateParCore(const arma::vec& thFull, const arma::vec& omegaDiag) {
  if ((int)_vaeInitPar.size() != (int)op_focei.npars)
    stop("vaeInnerUpdatePar_ requires vaeInnerSetup_ to be run first");
  if ((int)thFull.n_elem != (int)op_focei.ntheta)
    stop("vaeInnerUpdatePar_: theta size %d != ntheta %d",
         (int)thFull.n_elem, (int)op_focei.ntheta);
  if ((int)omegaDiag.n_elem != (int)op_focei.omegan)
    stop("vaeInnerUpdatePar_: omega size %d != omegan %d (setup needs a diagonal rxInv)",
         (int)omegaDiag.n_elem, (int)op_focei.omegan);
  // Reproduce vaeInnerSetup_'s tail EXACTLY without the full foceiSetup_:
  // foceiSetupTheta_ sets op_focei.initPar[k] to the natural parameter value and
  // then vaeInnerSetup_ calls updateTheta(initPar).  updateTheta()'s unscale is
  // relative to op_focei.initPar, so initPar MUST be refreshed to the new values
  // first (a stale initPar makes the unscaled fullTheta wrong).  Natural theta
  // for the theta block, omega_kk^(-1/4) for the diagonal "sqrt"-xform omega
  // block.  Fixed thetas are absent from fixedTrans and stay at their setup
  // values (which is what the VAE wants for a held theta).
  std::vector<double> par(op_focei.npars);
  for (unsigned int k = 0; k < op_focei.npars; ++k) {
    int j = op_focei.fixedTrans[k];
    par[k] = (j < (int)op_focei.ntheta) ? thFull[j] :
      pow(omegaDiag[j - op_focei.ntheta], -0.25);
    op_focei.initPar[k] = par[k];
  }
  updateTheta(par.data());
  // likInner0 short-circuits on an unchanged eta (fInd->oldEta), which is only
  // safe while theta/omega are fixed -- force a re-solve (like impForceResolve)
  rx = getRxSolve_();
  for (int id = getRxNsubAndMix(rx); id--;) inds_focei[id].setup = 0;
}

//[[Rcpp::export]]
RObject vaeInnerUpdatePar_(NumericVector thFull, NumericVector omegaDiag) {
  arma::vec th(thFull.begin(), thFull.size());
  arma::vec om(omegaDiag.begin(), omegaDiag.size());
  vaeInnerUpdateParCore(th, om);
  return R_NilValue;
}

// Per-id inner objective (and optional gradient) at the supplied etas, parallel
// over ids. id in [0, nSub) are the physical subjects; for mixtures the caller
// passes nSub*nMix rows (id = component*nSub + subject) and combines.  The
// arma-native core is shared with the C++ training loop (vaeTrainCpp_) so the
// per-step evaluations allocate no R objects.
static void vaeInnerLikCore(const arma::mat& etaMat, int cores, bool grad, bool preds,
                            arma::vec& obj, arma::mat& lp,
                            std::vector<std::vector<double> >& pf) {
  const int nid = (int)etaMat.n_rows;
  const int neta = (int)etaMat.n_cols;
  obj.set_size(nid);
  if (grad) lp.set_size(nid, neta);
  if (preds) { pf.clear(); pf.resize(nid); }
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  // rxode2 sizes every per-thread solve buffer by op->cores (1 for models it
  // flags not thread safe, e.g. linCmtB); more external threads than that
  // index past the pools and corrupt the heap.
  cores = min2(cores, getOpCores(op));
  // Mixture components (id = m*nsub + subject) share the physical subjects'
  // data/solving structures (the memory-saving default), so no two components
  // may ever be solved concurrently: solve component-by-component (serial m
  // OUTSIDE), parallelizing over physical subjects INSIDE each component.
  int nsub = (int)getRxNsub(rx);
  int nMix = (nsub > 0 && nid % nsub == 0) ? nid / nsub : 0;
  if (nMix == 0) { nsub = nid; nMix = 1; cores = 1; } // unexpected shape: serial
  const bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  if (doParallel) {
    sortIds(rx, 2);
    _innerParallel.store(1, std::memory_order_release);
  }
  for (int m = 0; m < nMix; ++m) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
    for (int i = 0; i < nsub; ++i) {
      int base = doParallel ? (getOrdId(rx, i) - 1) : i;
#ifdef _OPENMP
      if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
      int id = base + m * nsub;
      std::vector<double> eta(neta);
      for (int j = 0; j < neta; ++j) eta[j] = etaMat(id, j);
      if (grad) {
        std::vector<double> g(neta);
        lpInner(&eta[0], &g[0], id);
        for (int j = 0; j < neta; ++j) lp(id, j) = g[j];
      }
      obj[id] = likInner0(&eta[0], id);
      if (preds) {
        arma::mat rf = grabRFmatFromInner(id, false); // F,R for the solved component
        pf[id].assign(rf.colptr(0), rf.colptr(0) + rf.n_rows);
      }
#ifdef _OPENMP
      if (doParallel) setRxThreadId(-1);
#endif
    }
  }
  if (doParallel) {
    _innerParallel.store(0, std::memory_order_release);
    sortIds(rx, 0);
  }
}

//[[Rcpp::export]]
List vaeInnerLik(NumericMatrix etaMat, int cores, bool grad = false, bool preds = false) {
  const int nid = etaMat.nrow();
  const int neta = etaMat.ncol();
  arma::mat eta(etaMat.begin(), nid, neta, false, true);
  arma::vec obj;
  arma::mat lp;
  std::vector<std::vector<double> > pf;
  vaeInnerLikCore(eta, cores, grad, preds, obj, lp, pf);
  NumericVector objR(obj.begin(), obj.end());
  NumericMatrix lpR(grad ? nid : 1, grad ? neta : 1);
  if (grad) std::copy(lp.begin(), lp.end(), lpR.begin());
  List fl(preds ? nid : 0);
  if (preds) for (int id = 0; id < nid; ++id) fl[id] = NumericVector(pf[id].begin(), pf[id].end());
  return List::create(_["obj"] = objR, _["lp"] = lpR, _["f"] = fl);
}

// ===========================================================================
// ADVI (Kucukelbir et al. 2017) outer likelihood + gradient.
//
// The variational family lives in the unconstrained real coordinate space.  In
// the (mean-field, point-estimate) case each subject i carries variational
// params (mu_i, omega_i) with q(eta_i) = N(mu_i, diag(exp(2 omega_i))); the
// reparameterization eta_i = mu_i + exp(omega_i) .* eps_i pushes the gradient
// inside the expectation.  The per-subject log-joint and its eta-gradient are
// reused verbatim from the FOCEi inner problem: likInner0 returns
// -log p(y_i, eta_i) up to constants (the eta prior QUADRATIC is included via
// op_focei.omegaInv, but NOT the -0.5 logdet(2*pi*Omega) normalization), and its
// gradient d(likInner0)/d(eta) = fInd->lp.  So with L = likInner0, lp = fInd->lp:
//   grad_mu_ik     = -lp_ik
//   grad_omega_ik  = -lp_ik * exp(omega_ik) * eps_ik + 1        (entropy grad = 1)
// The population theta gradient (the OUTER gradient) uses the mu-ref identity
// d/dtheta = d/deta for mu-referenced thetas (data-term eta score
// Omega^-1 eta - lp), and impThetaScore (theta forward sensitivities) for the
// non-mu structural + residual-error thetas.  The population between-subject
// log-variances get the closed-form ELBO gradient
//   grad_logOmega_k = sum_i [0.5 (mu_ik^2 + var_ik)/w_k - 0.5]   (w_k = exp(logOmega_k))
// where the reparam draw enters -obj through omegaInv, and the -0.5 N logdet term
// supplies the -0.5.  All dropped additive terms are constant in every parameter,
// so the returned ELBO and gradients are mutually consistent (FD-checkable).
//
// Sets op_focei.fullTheta (natural scale) + propagates to each subject, and sets
// op_focei.omegaInv = diag(exp(-logPopOmega)) directly (diagonal population
// Omega, as est="vae" also estimates only the diagonal).  Serial for now; the
// optimization loop (adviLoop_) adds parallelism and the counter-based RNG.
// ===========================================================================

// VAE optimization printing + parameter-history capture. Reuses the SHARED
// iteration-print machinery in scale.h (scaleSetup / scalePrintHeader /
// scalePrintFun / scaleParHisDf) that saem, focei and the nlm family use, so the
// VAE prints the exact same iteration table and produces parHistData in the
// standard format -- no duplicated formatting. The VAE never scales its
// parameters, so scaleTypeNone drops the redundant "U" (Unscaled) rows; the
// R-side `xform` list (.iterPrintXParFromUi) drives the "X" back-transform row
// (e.g. exp() thetas), like focei/saem.
static scaling _vaeScale;
static std::vector<double> _vaeIpInit;
static std::vector<double> _vaeIpScaleC;
static std::vector<int> _vaeIpXPar;
static std::string _vaeIpPhase;

//[[Rcpp::export]]
RObject vaeIterPrintStart_(NumericVector initPar, CharacterVector names,
                           List iterPrintControl, RObject xform = R_NilValue) {
  int np = initPar.size();
  _vaeIpInit.assign(initPar.begin(), initPar.end());
  _vaeIpScaleC.assign(np, 1.0);
  scaleSetup(&_vaeScale, _vaeIpInit.data(), _vaeIpScaleC.data(), names,
             /*useColor*/ 0, /*printNcol*/ np, /*print*/ 1,
             /*normType*/ normTypeConstant, /*scaleType*/ scaleTypeNone,
             /*scaleCmin*/ 0.0, /*scaleCmax*/ 0.0, /*scaleTo*/ 0.0, np);
  if (!Rf_isNull(xform)) {
    scaleAttachXform(&_vaeScale, as<List>(xform));
  }
  if (_vaeScale.xPar == NULL) {
    // zero xform: natural-scale values, no back-transform. xPar must be
    // non-NULL because scalePrintFun indexes it unconditionally.
    _vaeIpXPar.assign(np, 0);
    _vaeScale.xPar = _vaeIpXPar.data();
    _vaeScale.probitIdx = NULL;
    _vaeScale.logitThetaLow = _vaeScale.logitThetaHi = NULL;
    _vaeScale.probitThetaLow = _vaeScale.probitThetaHi = NULL;
  }
  // phase legend for the labels shown in the Function-Val cell
  _vaeScale.keyExtra =
    "Burn in: encoder-only burn-in; KL anneal: KL-weight ramp;\n"
    "EM: main EM phase; Smooth: EMA-smoothing phase\n";
  scaleApplyIterPrintControl(&_vaeScale, iterPrintControl);
  scalePrintHeader(&_vaeScale);
  return R_NilValue;
}

//[[Rcpp::export]]
RObject vaeIterPrintRow_(NumericVector x, double f, std::string phase = "") {
  _vaeIpPhase = phase;
  _vaeScale.phaseLabel = _vaeIpPhase.empty() ? NULL : _vaeIpPhase.c_str();
  scalePrintFun(&_vaeScale, &x[0], f);
  return R_NilValue;
}

//[[Rcpp::export]]
RObject vaeIterPrintGet_(bool printLine = true) {
  _vaeScale.save = 0;
  _vaeScale.every = 0;
  if (printLine) scalePrintLine(&_vaeScale, min2(_vaeScale.npars, _vaeScale.ncol));
  return scaleParHisDf(&_vaeScale);
}

// ADVI iteration printing: same shared scale.h machinery/state as the vae
// helpers above (one loop runs at a time).  Rows are always captured into
// standard parHistData; iterPrintControl$every gates the console output.
// `it0` keeps a resumed run's iteration numbering global.
static void adviIterPrintStart(NumericVector initPar, CharacterVector names,
                               List iterPrintControl, RObject xform, int it0) {
  int np = initPar.size();
  _vaeIpInit.assign(initPar.begin(), initPar.end());
  _vaeIpScaleC.assign(np, 1.0);
  scaleSetup(&_vaeScale, _vaeIpInit.data(), _vaeIpScaleC.data(), names,
             /*useColor*/ 0, /*printNcol*/ np, /*print*/ 1,
             /*normType*/ normTypeConstant, /*scaleType*/ scaleTypeNone,
             /*scaleCmin*/ 0.0, /*scaleCmax*/ 0.0, /*scaleTo*/ 0.0, np);
  if (!Rf_isNull(xform)) {
    scaleAttachXform(&_vaeScale, as<List>(xform));
  }
  if (_vaeScale.xPar == NULL) {
    _vaeIpXPar.assign(np, 0);
    _vaeScale.xPar = _vaeIpXPar.data();
    _vaeScale.probitIdx = NULL;
    _vaeScale.logitThetaLow = _vaeScale.logitThetaHi = NULL;
    _vaeScale.probitThetaLow = _vaeScale.probitThetaHi = NULL;
  }
  _vaeScale.keyExtra =
    "Function-Val: negative Monte-Carlo ELBO (lower is better)\n"
    "srch <eta>: step-size-scale search run; SGA: main run\n";
  scaleApplyIterPrintControl(&_vaeScale, iterPrintControl);
  _vaeScale.cn = it0;
  scalePrintHeader(&_vaeScale);
}

// Print/record one ADVI iteration: the just-written parHist row + -ELBO.
static void adviIterPrintRow(NumericMatrix &parHist, int it, int np, double elbo,
                             const std::string &phase) {
  NumericVector r(np);
  for (int c = 0; c < np; ++c) r[c] = parHist(it, c);
  vaeIterPrintRow_(r, -elbo, phase);
}

// Set the natural-scale population thetas into op_focei.fullTheta and propagate
// to every subject's solve (mirrors updateTheta's propagation, but takes the
// natural theta directly rather than the scaled free-parameter vector).
static void adviSetTheta(NumericVector theta) {
  rx = getRxSolve_();
  int nsub = getRxNsub(rx);
  int ntheta = (int)op_focei.ntheta;
  int nt = min2((int)theta.size(), ntheta);
  for (int j = 0; j < nt; ++j) op_focei.fullTheta[j] = theta[j];
  for (int id = 0; id < nsub; ++id) {
    rx_solving_options_ind *ind = getSolvingOptionsInd(rx, getRxId(id));
    for (int j = 0; j < nt; ++j) setIndParPtr(ind, op_focei.thetaTrans[j], op_focei.fullTheta[j]);
  }
}

// Set the diagonal population Omega used by the inner prior term (likInner0 reads
// op_focei.omegaInv only; the logdet normalization is added in the ELBO).
static void adviSetPopOmega(NumericVector logPopOmega) {
  int neta = (int)op_focei.neta;
  arma::mat oi(neta, neta, arma::fill::zeros);
  for (int k = 0; k < neta; ++k) oi(k, k) = std::exp(-logPopOmega[k]);
  op_focei.omegaInv = oi;
}

// Diagnostic: report the wired theta-sensitivity offsets after setup.
//[[Rcpp::export]]
List adviThetaSensInfo_() {
  return List::create(_["nSens"] = op_focei.impThetaSensIdx.size(),
                      _["predOffset"] = op_focei.thetaSensPredOffset,
                      _["rOffset"] = op_focei.thetaSensROffset,
                      _["fOffset"] = op_focei.thetaSensOffset,
                      _["dvOffset"] = op_focei.thetaSensDvOffset,
                      _["neq"] = op_focei.thetaSensNeq,
                      _["calcLhsNull"] = (rxThetaSens.calc_lhs == NULL),
                      _["isAdvi"] = op_focei.isAdvi,
                      _["isImpmap"] = op_focei.isImpmap);
}

// Divergence test for the adaptEta step-size search: a non-finite ELBO or any
// recorded population parameter exceeding 1e4 in magnitude (same criterion the R
// scorer .adviAdaptEta applies to the full trace).  `row` is the just-written
// parHist iteration, `ncol` its width.
static bool adviDiverged(double elbo, NumericMatrix &parHist, int row, int ncol) {
  if (!R_finite(elbo)) return true;
  for (int c = 0; c < ncol; ++c) {
    double v = parHist(row, c);
    if (!R_finite(v) || std::fabs(v) > 1e4) return true;
  }
  return false;
}

// Core (mean-field, point-estimate) ELBO + gradient at a fixed eps draw.  Fills
// the pre-allocated gMu/gOmega/gTheta/gPopLogOmega (zeroed here) and returns the
// ELBO.  Shared by the R-facing adviElboGrad_ (FD test) and the C++ loop adviLoop_.
static double adviElboGradCore(NumericMatrix mu, NumericMatrix omega, NumericVector theta,
                               NumericVector logPopOmega, NumericMatrix eps,
                               IntegerVector muRefThetaIdx,
                               NumericMatrix &gMu, NumericMatrix &gOmega,
                               NumericVector &gTheta, NumericVector &gPopLogOmega,
                               int cores) {
  const int N = mu.nrow();
  const int neta = mu.ncol();
  const int ntheta = (int)op_focei.ntheta;
  adviSetTheta(theta);
  adviSetPopOmega(logPopOmega);
  arma::mat omegaInv = op_focei.omegaInv;
  // Make every ELBO evaluation a pure, deterministic function of its inputs: the
  // FOCEi bad-solve tol relaxation otherwise persists per-subject (and in the
  // global atomics) across evaluations, so a solve that loosened tolerance in an
  // earlier iteration/call would make the result depend on history -- breaking
  // bit-for-bit resume==fresh reproducibility.  Reset to base each call.
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  op_focei.reducedTol.store(0, std::memory_order_relaxed);
  op_focei.reducedTol2.store(0, std::memory_order_relaxed);
  op_focei.stickyTol.store(0, std::memory_order_relaxed);
  resetOpBadSolve(op);
  std::fill(gMu.begin(), gMu.end(), 0.0);
  std::fill(gOmega.begin(), gOmega.end(), 0.0);
  std::fill(gTheta.begin(), gTheta.end(), 0.0);
  std::fill(gPopLogOmega.begin(), gPopLogOmega.end(), 0.0);
  const int nSens = op_focei.impThetaSensIdx.size();
  // Per-subject buffers reduced serially (in id order) after the parallel loop so
  // the floating-point summation order -- hence the result -- is identical for any
  // thread count and identical to the serial code.  gMu/gOmega rows are already
  // per-subject, so they are written directly inside the loop.
  std::vector<double> objBuf(N, 0.0);            // -obj + mean-field entropy
  arma::mat gThetaBuf(N, ntheta, arma::fill::zeros);
  arma::mat gLpoBuf(N, neta, arma::fill::zeros);
  // rxode2 sizes every per-thread solve buffer by op->cores; more external threads
  // than that index past the pools and corrupt the heap (see vaeInnerLikCore).
  cores = min2(cores, getOpCores(op));
  int nsub = (int)getRxNsub(rx);
  if (nsub != N) { nsub = N; cores = 1; }        // unexpected shape: serial
  const bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  if (doParallel) { sortIds(rx, 2); _innerParallel.store(1, std::memory_order_release); }
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int ii = 0; ii < nsub; ++ii) {
    int i = doParallel ? (getOrdId(rx, ii) - 1) : ii;
#ifdef _OPENMP
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    // reparameterized draw eta_i = mu_i + exp(omega_i) .* eps_i
    std::vector<double> eta(neta);
    arma::vec etav(neta);
    for (int k = 0; k < neta; ++k) {
      eta[k] = mu(i, k) + std::exp(omega(i, k)) * eps(i, k);
      etav[k] = eta[k];
    }
    // force a full inner recompute: likInner0 caches on unchanged eta, but theta
    // and the population Omega also change between ELBO evaluations.  Reset this
    // subject's tolerance factor + solve state to base so the solve is history-
    // independent (deterministic).
    inds_focei[i].setup = 0;
    { rx_solving_options_ind *indI = getSolvingOptionsInd(rx, getRxId(i));
      setIndTolFactor(indI, 1.0);
      setIndSolve(indI, -1); }
    double obj = likInner0(&eta[0], i);          // -log p(y_i, eta_i) + const
    focei_ind *fInd = &(inds_focei[i]);
    arma::vec lp(neta);
    for (int k = 0; k < neta; ++k) lp[k] = fInd->lp[k];

    // ELBO: data+prior term (-obj) + mean-field entropy (sum omega, + const)
    double e = -obj;
    for (int k = 0; k < neta; ++k) e += omega(i, k);
    objBuf[i] = e;

    // variational-parameter gradients
    for (int k = 0; k < neta; ++k) {
      gMu(i, k) = -lp[k];
      gOmega(i, k) = -lp[k] * std::exp(omega(i, k)) * eps(i, k) + 1.0;
    }

    // population between-subject log-variance gradient (per-subject):
    //   0.5 (eta_ik^2)/w_k - 0.5   with var absorbed via the reparam draw
    for (int k = 0; k < neta; ++k) {
      double wk = std::exp(logPopOmega[k]);
      gLpoBuf(i, k) = 0.5 * (eta[k] * eta[k]) / wk - 0.5;
    }

    // OUTER population theta gradient
    // mu-referenced thetas: data-term eta score (Omega^-1 eta - lp)
    arma::vec dataScore = omegaInv * etav - lp;
    for (int k = 0; k < neta; ++k) {
      int p = muRefThetaIdx[k];                  // 1-based ntheta, NA_INTEGER if free
      if (p != NA_INTEGER && p >= 1 && p <= ntheta) gThetaBuf(i, p - 1) += dataScore[k];
    }
    // non-mu structural + residual-error thetas: theta forward-sensitivity score
    if (nSens > 0) {
      arma::mat S(1, neta);
      for (int k = 0; k < neta; ++k) S(0, k) = eta[k];
      arma::vec zk(1); zk[0] = 1.0;
      arma::vec g(nSens, arma::fill::zeros);
      arma::mat H(nSens, nSens, arma::fill::zeros);
      impThetaScore(i, S, zk, g, H);
      for (int s = 0; s < nSens; ++s) {
        int th = op_focei.impThetaSensIdx[s];    // 0-based
        if (th >= 0 && th < ntheta) gThetaBuf(i, th) += g[s];
      }
    }
#ifdef _OPENMP
    if (doParallel) setRxThreadId(-1);
#endif
  }
  if (doParallel) { _innerParallel.store(0, std::memory_order_release); sortIds(rx, 0); }
  // serial reduction in id order (thread-count invariant, matches serial code)
  double elbo = 0.0;
  for (int i = 0; i < N; ++i) {
    elbo += objBuf[i];
    for (int p = 0; p < ntheta; ++p) gTheta[p] += gThetaBuf(i, p);
    for (int k = 0; k < neta; ++k) gPopLogOmega[k] += gLpoBuf(i, k);
  }
  // population prior logdet: -0.5 N sum_k logPopOmega_k
  for (int k = 0; k < neta; ++k) elbo += -0.5 * N * logPopOmega[k];
  return elbo;
}

//[[Rcpp::export]]
List adviElboGrad_(NumericMatrix mu, NumericMatrix omega, NumericVector theta,
                   NumericVector logPopOmega, NumericMatrix eps,
                   IntegerVector muRefThetaIdx) {
  const int N = mu.nrow(), neta = mu.ncol(), ntheta = (int)op_focei.ntheta;
  NumericMatrix gMu(N, neta), gOmega(N, neta);
  NumericVector gTheta(ntheta), gPopLogOmega(neta);
  double elbo = adviElboGradCore(mu, omega, theta, logPopOmega, eps, muRefThetaIdx,
                                 gMu, gOmega, gTheta, gPopLogOmega, 1);
  return List::create(_["elbo"] = elbo, _["gMu"] = gMu, _["gOmega"] = gOmega,
                      _["gTheta"] = gTheta, _["gPopLogOmega"] = gPopLogOmega);
}

// ADVI stochastic-gradient-ascent loop (mean-field, point-estimate), 100% in
// C++.  Each iteration draws eps from the counter-based threefry stream keyed by
// (seed, global-iteration, mc-sample) -- so a shorter run is a bit-for-bit prefix
// of a longer one and results are independent of thread count -- averages the
// ELBO gradient over nMc samples, then updates every parameter (variational
// mu/omega and the point-estimate theta/logPopOmega) with the paper's adaptive
// step-size (Eqs 10-11): s = alpha g^2 + (1-alpha) s (s^(1)=g^2),
// rho = etaScale * i^(-1/2+eps) / (tau + sqrt(s)).  Fixed thetas / omega diagonal
// entries are held.  `it0` and the passed-in s-accumulators support warm resume
// (continue a finished fit): global iteration index = it0 + local.
//[[Rcpp::export]]
List adviLoop_(NumericMatrix mu0, NumericMatrix omega0, NumericVector theta0,
               NumericVector logPopOmega0, IntegerVector muRefThetaIdx,
               IntegerVector thetaMuRefEta,
               LogicalVector thetaFix, LogicalVector omegaFix,
               int iters, double seed, double etaScale, double tau, double alpha,
               int nMc, int it0,
               NumericMatrix sMu0, NumericMatrix sOmega0,
               NumericVector sTheta0, NumericVector sLpo0, int cores, int divergeStop,
               CharacterVector parNames, RObject iterPrintControl, RObject xform,
               std::string ipPhase, int ipStart, int ipEnd) {
  const int N = mu0.nrow(), neta = mu0.ncol(), ntheta = theta0.size();
  NumericMatrix mu = clone(mu0), omega = clone(omega0);
  NumericVector theta = clone(theta0), logPopOmega = clone(logPopOmega0);
  NumericMatrix sMu = clone(sMu0), sOmega = clone(sOmega0);
  NumericVector sTheta = clone(sTheta0), sLpo = clone(sLpo0);
  // gradient buffers reused each iteration
  NumericMatrix gMu(N, neta), gOmega(N, neta);
  NumericVector gTheta(ntheta), gPopLogOmega(neta);
  NumericMatrix gMuAcc(N, neta), gOmAcc(N, neta);
  NumericVector gThAcc(ntheta), gLpoAcc(neta);
  NumericMatrix eps(N, neta);
  NumericVector elboTrace(iters);
  NumericMatrix parHist(iters, ntheta + neta);   // theta ... , popOmega diag ...
  const double eeps = 1e-16, invM = 1.0 / nMc;
  int itRun = iters;                              // may shrink if divergeStop fires
  const bool doIp = !Rf_isNull(iterPrintControl);
  if (doIp && ipStart) {
    NumericVector row0(ntheta + neta);
    for (int p = 0; p < ntheta; ++p) row0[p] = theta[p];
    for (int k = 0; k < neta; ++k) row0[ntheta + k] = std::exp(logPopOmega[k]);
    adviIterPrintStart(row0, parNames, as<List>(iterPrintControl), xform, it0);
  }

  for (int it = 0; it < iters; ++it) {
    int gi = it0 + it;                            // global iteration index
    std::fill(gMuAcc.begin(), gMuAcc.end(), 0.0);
    std::fill(gOmAcc.begin(), gOmAcc.end(), 0.0);
    std::fill(gThAcc.begin(), gThAcc.end(), 0.0);
    std::fill(gLpoAcc.begin(), gLpoAcc.end(), 0.0);
    double elboAcc = 0.0;
    for (int m = 0; m < nMc; ++m) {
      // counter-based standard-normal draw, stream keyed by (seed, gi, m)
      setSeedEng1((uint32_t)seed + (uint32_t)((gi * nMc + m) & 0x7fffffff));
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < neta; ++k) eps(i, k) = rxNormEng(0.0, 1.0);
      elboAcc += adviElboGradCore(mu, omega, theta, logPopOmega, eps, muRefThetaIdx,
                                  gMu, gOmega, gTheta, gPopLogOmega, cores);
      for (int j = 0; j < N * neta; ++j) { gMuAcc[j] += gMu[j]; gOmAcc[j] += gOmega[j]; }
      for (int p = 0; p < ntheta; ++p) gThAcc[p] += gTheta[p];
      for (int k = 0; k < neta; ++k) gLpoAcc[k] += gPopLogOmega[k];
    }
    elboAcc *= invM;
    for (int j = 0; j < N * neta; ++j) { gMuAcc[j] *= invM; gOmAcc[j] *= invM; }
    for (int p = 0; p < ntheta; ++p) gThAcc[p] *= invM;
    for (int k = 0; k < neta; ++k) gLpoAcc[k] *= invM;
    elboTrace[it] = elboAcc;

    // paper adaptive step-size (Eq 10-11); global-iteration decay so resume keeps
    // the same schedule as a fresh longer run.
    double idecay = std::pow((double)(gi + 1), -0.5 + eeps);
    for (int j = 0; j < N * neta; ++j) {
      double g = gMuAcc[j];
      sMu[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sMu[j];
      mu[j] += etaScale * idecay / (tau + std::sqrt(sMu[j])) * g;
      g = gOmAcc[j];
      sOmega[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sOmega[j];
      omega[j] += etaScale * idecay / (tau + std::sqrt(sOmega[j])) * g;
    }
    for (int p = 0; p < ntheta; ++p) {
      // mu-referenced intercepts are not gradient-updated: theta and eta share a
      // flat direction (data constrains only theta+eta), which makes joint SGA
      // drift/diverge.  They are updated by the recentering M-step below instead.
      if (thetaFix[p] || thetaMuRefEta[p] >= 0) continue;
      double g = gThAcc[p];
      sTheta[p] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sTheta[p];
      theta[p] += etaScale * idecay / (tau + std::sqrt(sTheta[p])) * g;
    }
    // Recentering M-step for mu-referenced intercepts: shift the mean of each
    // mu-referenced eta's variational means into its typical-value theta
    // (theta+eta is invariant, so the data fit is unchanged; centering the etas
    // lowers the prior penalty -- an ELBO-ascent move -- and removes the flat
    // direction that otherwise diverges).
    for (int p = 0; p < ntheta; ++p) {
      int k = thetaMuRefEta[p];
      if (k < 0 || thetaFix[p]) continue;
      double mbar = 0.0;
      for (int i = 0; i < N; ++i) mbar += mu(i, k);
      mbar /= N;
      theta[p] += mbar;
      for (int i = 0; i < N; ++i) mu(i, k) -= mbar;
    }
    // Population between-subject variance M-step, EMA-DAMPED.  The undamped
    // ELBO-maximizing value is w_k = mean_i(mu_ik^2 + var_ik) (var_ik =
    // exp(2 omega_ik)) -- the vae/SAEM omega update -- but applying it every
    // iteration makes the prior track the variational scale exactly, which
    // cancels the restoring force on the variational log-sd (exp(2 omega)/w -> 1)
    // and the posterior inflates without bound.  A lagging (damped) w stays
    // tighter than the variational scale during transients, preserving the
    // restoring force; it converges to the same fixed point once the E-step
    // settles.  gamma anneals from a small value toward 1.
    double gam = 1.0 / (10.0 + (double)gi);
    for (int k = 0; k < neta; ++k) {
      if (omegaFix[k]) continue;
      double s = 0.0;
      for (int i = 0; i < N; ++i) s += mu(i, k) * mu(i, k) + std::exp(2.0 * omega(i, k));
      double wTarget = s / N, wCur = std::exp(logPopOmega[k]);
      logPopOmega[k] = std::log((1.0 - gam) * wCur + gam * wTarget);
    }
    for (int p = 0; p < ntheta; ++p) parHist(it, p) = theta[p];
    for (int k = 0; k < neta; ++k) parHist(it, ntheta + k) = std::exp(logPopOmega[k]);
    if (doIp) adviIterPrintRow(parHist, it, ntheta + neta, elboAcc, ipPhase);
    // adaptEta candidates only (divergeStop != 0): stop paying for a run that has
    // clearly blown up.  The offending iteration is recorded first, so the
    // (zero-padded) trace still trips R's divergence test; the main run passes
    // divergeStop == 0 and is bit-for-bit unchanged.
    if (divergeStop != 0 && adviDiverged(elboAcc, parHist, it, ntheta + neta)) {
      itRun = it + 1; break;
    }
  }
  RObject parHistData = (doIp && ipEnd) ? vaeIterPrintGet_(_vaeScale.every != 0) : RObject(R_NilValue);
  return List::create(_["mu"] = mu, _["omega"] = omega, _["theta"] = theta,
                      _["logPopOmega"] = logPopOmega, _["elbo"] = elboTrace,
                      _["parHist"] = parHist, _["it0"] = it0 + itRun, _["itRun"] = itRun,
                      _["sMu"] = sMu, _["sOmega"] = sOmega,
                      _["sTheta"] = sTheta, _["sLpo"] = sLpo,
                      _["parHistData"] = parHistData);
}

// ===========================================================================
// Block full-rank ADVI: q(eta_i) = N(mu_i, L_i L_i^T) with L_i lower-triangular
// (neta x neta), stored packed row-major lower-tri (index (i,j), i>=j, at
// i*(i+1)/2 + j).  The reparameterization is eta_i = mu_i + L_i eps_i.  The
// entropy is 0.5 log det(L_i L_i^T) = sum_k log|L_i(k,k)| -- it depends ONLY on
// the diagonal, so the entropy gradient is 1/L_kk on the diagonal and 0 off it
// (no matrix inverse needed).  The data+prior gradient is
//   grad_L(i,j) = E[(d/deta log p)_i eps_j] = E[(-lp)_i eps_j]   (i >= j),
// giving grad_L(i,j) = (-lp)_i eps_j + [i==j] / L_ii.  Everything else (theta
// gradient, population omega M-step) matches the mean-field path with the
// per-subject variance var_ik = (L_i L_i^T)_kk = sum_{j<=k} L_i(k,j)^2.
// ===========================================================================

static double adviElboGradCoreFR(NumericMatrix mu, NumericMatrix Lpack, NumericVector theta,
                                 NumericVector logPopOmega, NumericMatrix eps,
                                 IntegerVector muRefThetaIdx,
                                 NumericMatrix &gMu, NumericMatrix &gL,
                                 NumericVector &gTheta, NumericVector &gPopLogOmega,
                                 int cores) {
  const int N = mu.nrow();
  const int neta = mu.ncol();
  const int ntheta = (int)op_focei.ntheta;
  adviSetTheta(theta);
  adviSetPopOmega(logPopOmega);
  arma::mat omegaInv = op_focei.omegaInv;
  rx = getRxSolve_();
  rx_solving_options *op = getSolvingOptions(rx);
  op_focei.reducedTol.store(0, std::memory_order_relaxed);
  op_focei.reducedTol2.store(0, std::memory_order_relaxed);
  op_focei.stickyTol.store(0, std::memory_order_relaxed);
  resetOpBadSolve(op);
  std::fill(gMu.begin(), gMu.end(), 0.0);
  std::fill(gL.begin(), gL.end(), 0.0);
  std::fill(gTheta.begin(), gTheta.end(), 0.0);
  std::fill(gPopLogOmega.begin(), gPopLogOmega.end(), 0.0);
  const int nSens = op_focei.impThetaSensIdx.size();
  // Per-subject buffers reduced serially (id order) after the parallel loop so the
  // result is thread-count invariant and identical to the serial code (see the
  // mean-field core).  gMu/gL rows are per-subject, written directly in the loop.
  std::vector<double> objBuf(N, 0.0);
  arma::mat gThetaBuf(N, ntheta, arma::fill::zeros);
  arma::mat gLpoBuf(N, neta, arma::fill::zeros);
  cores = min2(cores, getOpCores(op));
  int nsub = (int)getRxNsub(rx);
  if (nsub != N) { nsub = N; cores = 1; }        // unexpected shape: serial
  const bool doParallel = (cores > 1) && solveMethodThreadSafe(op);
  if (doParallel) { sortIds(rx, 2); _innerParallel.store(1, std::memory_order_release); }
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic) if(doParallel)
#endif
  for (int ss = 0; ss < nsub; ++ss) {
    int s = doParallel ? (getOrdId(rx, ss) - 1) : ss;
#ifdef _OPENMP
    if (doParallel) setRxThreadId(omp_get_thread_num());
#endif
    // eta_i = mu_i + L_i eps_i (L_i lower-tri, packed)
    std::vector<double> eta(neta);
    arma::vec etav(neta);
    for (int i = 0; i < neta; ++i) {
      double e = mu(s, i);
      for (int j = 0; j <= i; ++j) e += Lpack(s, i * (i + 1) / 2 + j) * eps(s, j);
      eta[i] = e; etav[i] = e;
    }
    inds_focei[s].setup = 0;
    { rx_solving_options_ind *indI = getSolvingOptionsInd(rx, getRxId(s));
      setIndTolFactor(indI, 1.0); setIndSolve(indI, -1); }
    double obj = likInner0(&eta[0], s);
    focei_ind *fInd = &(inds_focei[s]);
    arma::vec lp(neta);
    for (int k = 0; k < neta; ++k) lp[k] = fInd->lp[k];

    double e = -obj;
    for (int k = 0; k < neta; ++k)
      e += std::log(std::fabs(Lpack(s, k * (k + 1) / 2 + k)));   // entropy: sum log|L_kk|
    objBuf[s] = e;

    for (int i = 0; i < neta; ++i) {
      gMu(s, i) = -lp[i];
      for (int j = 0; j <= i; ++j) {
        double g = -lp[i] * eps(s, j);
        if (i == j) g += 1.0 / Lpack(s, i * (i + 1) / 2 + i);   // diagonal entropy grad
        gL(s, i * (i + 1) / 2 + j) = g;
      }
    }
    for (int k = 0; k < neta; ++k) {
      double wk = std::exp(logPopOmega[k]);
      gLpoBuf(s, k) = 0.5 * (eta[k] * eta[k]) / wk - 0.5;
    }
    arma::vec dataScore = omegaInv * etav - lp;
    for (int k = 0; k < neta; ++k) {
      int p = muRefThetaIdx[k];
      if (p != NA_INTEGER && p >= 1 && p <= ntheta) gThetaBuf(s, p - 1) += dataScore[k];
    }
    if (nSens > 0) {
      arma::mat S(1, neta);
      for (int k = 0; k < neta; ++k) S(0, k) = eta[k];
      arma::vec zk(1); zk[0] = 1.0;
      arma::vec g(nSens, arma::fill::zeros);
      arma::mat H(nSens, nSens, arma::fill::zeros);
      impThetaScore(s, S, zk, g, H);
      for (int t = 0; t < nSens; ++t) {
        int th = op_focei.impThetaSensIdx[t];
        if (th >= 0 && th < ntheta) gThetaBuf(s, th) += g[t];
      }
    }
#ifdef _OPENMP
    if (doParallel) setRxThreadId(-1);
#endif
  }
  if (doParallel) { _innerParallel.store(0, std::memory_order_release); sortIds(rx, 0); }
  double elbo = 0.0;
  for (int i = 0; i < N; ++i) {
    elbo += objBuf[i];
    for (int p = 0; p < ntheta; ++p) gTheta[p] += gThetaBuf(i, p);
    for (int k = 0; k < neta; ++k) gPopLogOmega[k] += gLpoBuf(i, k);
  }
  for (int k = 0; k < neta; ++k) elbo += -0.5 * N * logPopOmega[k];
  return elbo;
}

//[[Rcpp::export]]
List adviElboGradFR_(NumericMatrix mu, NumericMatrix Lpack, NumericVector theta,
                     NumericVector logPopOmega, NumericMatrix eps,
                     IntegerVector muRefThetaIdx) {
  const int N = mu.nrow(), neta = mu.ncol(), ntheta = (int)op_focei.ntheta;
  const int nL = neta * (neta + 1) / 2;
  NumericMatrix gMu(N, neta), gL(N, nL);
  NumericVector gTheta(ntheta), gPopLogOmega(neta);
  double elbo = adviElboGradCoreFR(mu, Lpack, theta, logPopOmega, eps, muRefThetaIdx,
                                   gMu, gL, gTheta, gPopLogOmega, 1);
  return List::create(_["elbo"] = elbo, _["gMu"] = gMu, _["gL"] = gL,
                      _["gTheta"] = gTheta, _["gPopLogOmega"] = gPopLogOmega);
}

//[[Rcpp::export]]
List adviLoopFR_(NumericMatrix mu0, NumericMatrix Lpack0, NumericVector theta0,
                 NumericVector logPopOmega0, IntegerVector muRefThetaIdx,
                 IntegerVector thetaMuRefEta,
                 LogicalVector thetaFix, LogicalVector omegaFix,
                 int iters, double seed, double etaScale, double tau, double alpha,
                 int nMc, int it0,
                 NumericMatrix sMu0, NumericMatrix sL0,
                 NumericVector sTheta0, NumericVector sLpo0, int cores, int divergeStop,
                 CharacterVector parNames, RObject iterPrintControl, RObject xform,
                 std::string ipPhase, int ipStart, int ipEnd) {
  const int N = mu0.nrow(), neta = mu0.ncol(), ntheta = theta0.size();
  const int nL = neta * (neta + 1) / 2;
  NumericMatrix mu = clone(mu0), Lpack = clone(Lpack0);
  NumericVector theta = clone(theta0), logPopOmega = clone(logPopOmega0);
  NumericMatrix sMu = clone(sMu0), sL = clone(sL0);
  NumericVector sTheta = clone(sTheta0), sLpo = clone(sLpo0);
  NumericMatrix gMu(N, neta), gL(N, nL);
  NumericVector gTheta(ntheta), gPopLogOmega(neta);
  NumericMatrix gMuAcc(N, neta), gLAcc(N, nL);
  NumericVector gThAcc(ntheta), gLpoAcc(neta);
  NumericMatrix eps(N, neta);
  NumericVector elboTrace(iters);
  NumericMatrix parHist(iters, ntheta + neta);
  const double eeps = 1e-16, invM = 1.0 / nMc;
  int itRun = iters;
  const bool doIp = !Rf_isNull(iterPrintControl);
  if (doIp && ipStart) {
    NumericVector row0(ntheta + neta);
    for (int p = 0; p < ntheta; ++p) row0[p] = theta[p];
    for (int k = 0; k < neta; ++k) row0[ntheta + k] = std::exp(logPopOmega[k]);
    adviIterPrintStart(row0, parNames, as<List>(iterPrintControl), xform, it0);
  }

  for (int it = 0; it < iters; ++it) {
    int gi = it0 + it;
    std::fill(gMuAcc.begin(), gMuAcc.end(), 0.0);
    std::fill(gLAcc.begin(), gLAcc.end(), 0.0);
    std::fill(gThAcc.begin(), gThAcc.end(), 0.0);
    std::fill(gLpoAcc.begin(), gLpoAcc.end(), 0.0);
    double elboAcc = 0.0;
    for (int m = 0; m < nMc; ++m) {
      setSeedEng1((uint32_t)seed + (uint32_t)((gi * nMc + m) & 0x7fffffff));
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < neta; ++k) eps(i, k) = rxNormEng(0.0, 1.0);
      elboAcc += adviElboGradCoreFR(mu, Lpack, theta, logPopOmega, eps, muRefThetaIdx,
                                    gMu, gL, gTheta, gPopLogOmega, cores);
      for (int j = 0; j < N * neta; ++j) gMuAcc[j] += gMu[j];
      for (int j = 0; j < N * nL; ++j) gLAcc[j] += gL[j];
      for (int p = 0; p < ntheta; ++p) gThAcc[p] += gTheta[p];
      for (int k = 0; k < neta; ++k) gLpoAcc[k] += gPopLogOmega[k];
    }
    elboAcc *= invM;
    for (int j = 0; j < N * neta; ++j) gMuAcc[j] *= invM;
    for (int j = 0; j < N * nL; ++j) gLAcc[j] *= invM;
    for (int p = 0; p < ntheta; ++p) gThAcc[p] *= invM;
    for (int k = 0; k < neta; ++k) gLpoAcc[k] *= invM;
    elboTrace[it] = elboAcc;

    double idecay = std::pow((double)(gi + 1), -0.5 + eeps);
    for (int j = 0; j < N * neta; ++j) {
      double g = gMuAcc[j];
      sMu[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sMu[j];
      mu[j] += etaScale * idecay / (tau + std::sqrt(sMu[j])) * g;
    }
    for (int j = 0; j < N * nL; ++j) {
      double g = gLAcc[j];
      sL[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sL[j];
      Lpack[j] += etaScale * idecay / (tau + std::sqrt(sL[j])) * g;
    }
    for (int p = 0; p < ntheta; ++p) {
      if (thetaFix[p] || thetaMuRefEta[p] >= 0) continue;
      double g = gThAcc[p];
      sTheta[p] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sTheta[p];
      theta[p] += etaScale * idecay / (tau + std::sqrt(sTheta[p])) * g;
    }
    // recenter mu-referenced intercepts (same as mean-field)
    for (int p = 0; p < ntheta; ++p) {
      int k = thetaMuRefEta[p];
      if (k < 0 || thetaFix[p]) continue;
      double mbar = 0.0;
      for (int i = 0; i < N; ++i) mbar += mu(i, k);
      mbar /= N;
      theta[p] += mbar;
      for (int i = 0; i < N; ++i) mu(i, k) -= mbar;
    }
    // damped population variance M-step; var_ik = (L_i L_i^T)_kk = sum_{j<=k} L(k,j)^2
    double gam = 1.0 / (10.0 + (double)gi);
    for (int k = 0; k < neta; ++k) {
      if (omegaFix[k]) continue;
      double sAcc = 0.0;
      for (int i = 0; i < N; ++i) {
        double var = 0.0;
        for (int j = 0; j <= k; ++j) { double l = Lpack(i, k * (k + 1) / 2 + j); var += l * l; }
        sAcc += mu(i, k) * mu(i, k) + var;
      }
      double wTarget = sAcc / N, wCur = std::exp(logPopOmega[k]);
      logPopOmega[k] = std::log((1.0 - gam) * wCur + gam * wTarget);
    }
    for (int p = 0; p < ntheta; ++p) parHist(it, p) = theta[p];
    for (int k = 0; k < neta; ++k) parHist(it, ntheta + k) = std::exp(logPopOmega[k]);
    if (doIp) adviIterPrintRow(parHist, it, ntheta + neta, elboAcc, ipPhase);
    if (divergeStop != 0 && adviDiverged(elboAcc, parHist, it, ntheta + neta)) {
      itRun = it + 1; break;
    }
  }
  RObject parHistData = (doIp && ipEnd) ? vaeIterPrintGet_(_vaeScale.every != 0) : RObject(R_NilValue);
  return List::create(_["mu"] = mu, _["Lpack"] = Lpack, _["theta"] = theta,
                      _["logPopOmega"] = logPopOmega, _["elbo"] = elboTrace,
                      _["parHist"] = parHist, _["it0"] = it0 + itRun, _["itRun"] = itRun,
                      _["sMu"] = sMu, _["sL"] = sL,
                      _["sTheta"] = sTheta, _["sLpo"] = sLpo,
                      _["parHistData"] = parHistData);
}

// ===========================================================================
// FULL-BAYES ADVI: the population parameters phi = (free thetas, free population
// log-variances) get their OWN full-rank variational Gaussian q(phi)=N(mPop,
// Lpop Lpop^T) (flat priors), on top of the per-subject q(eta_i) (mean-field or
// block full-rank).  Each iteration draws phi = mPop + Lpop eps_pop, maps it into
// the theta / logPopOmega vectors, and evaluates the SAME per-subject ELBO core
// -- whose gTheta / gPopLogOmega are exactly d(ELBO)/d(phi) at the drawn phi.
// The population block is then a standard reparameterization update:
//   grad_mPop = grad_phi,  grad_Lpop(i,j) = grad_phi_i eps_pop_j + [i==j]/Lpop_ii.
// The mu-referenced flat direction is handled the same way as point-estimate: the
// per-subject eta means are recentered into the corresponding phi (theta) means.
// phiThetaIdx[j] / phiOmIdx[j] map phi component j to a 0-based theta / eta index
// (the other is -1); phiMuRef[j] is the 0-based eta a mu-ref-theta phi recenters.
// ===========================================================================

//[[Rcpp::export]]
List adviLoopFB_(NumericMatrix mu0, NumericMatrix scale0, NumericVector theta0,
                 NumericVector logPopOmega0, NumericVector mPop0, NumericVector LpopPack0,
                 IntegerVector phiThetaIdx, IntegerVector phiOmIdx, IntegerVector phiMuRef,
                 IntegerVector muRefThetaIdx, int fr,
                 int iters, double seed, double etaScale, double tau, double alpha,
                 int nMc, int it0,
                 NumericMatrix sMu0, NumericMatrix sScale0,
                 NumericVector smPop0, NumericVector sLpop0, int cores, int divergeStop,
                 CharacterVector parNames, RObject iterPrintControl, RObject xform,
                 std::string ipPhase, int ipStart, int ipEnd) {
  const int N = mu0.nrow(), neta = mu0.ncol(), ntheta = theta0.size();
  const int nScale = scale0.ncol();
  const int npop = mPop0.size(), nLpop = LpopPack0.size();
  NumericMatrix mu = clone(mu0), scale = clone(scale0);
  NumericVector mPop = clone(mPop0), Lpop = clone(LpopPack0);
  NumericVector theta = clone(theta0), logPopOmega = clone(logPopOmega0);
  NumericMatrix sMu = clone(sMu0), sScale = clone(sScale0);
  NumericVector smPop = clone(smPop0), sLpop = clone(sLpop0);
  // buffers from the per-subject core
  NumericMatrix gMu(N, neta), gScale(N, nScale);
  NumericVector gTheta(ntheta), gPopLogOmega(neta);
  // accumulators
  NumericMatrix gMuAcc(N, neta), gScaleAcc(N, nScale);
  NumericVector gmPopAcc(npop), gLpopAcc(nLpop), epsPop(npop);
  NumericMatrix eps(N, neta);
  // The mu-referenced theta data score (sum_i Omega^-1 eta - lp) IS needed here:
  // it provides the curvature that pins the population posterior VARIANCE of the
  // mu-referenced thetas (the recentering below only fixes their posterior mean).
  // So pass the real muRefThetaIdx to the core, not an all-NA vector.
  NumericVector elboTrace(iters);
  NumericMatrix parHist(iters, ntheta + neta);
  const double eeps = 1e-16, invM = 1.0 / nMc;
  int itRun = iters;
  const bool doIp = !Rf_isNull(iterPrintControl);
  if (doIp && ipStart) {
    NumericVector row0(ntheta + neta);
    for (int p = 0; p < ntheta; ++p) row0[p] = theta[p];
    for (int k = 0; k < neta; ++k) row0[ntheta + k] = std::exp(logPopOmega[k]);
    adviIterPrintStart(row0, parNames, as<List>(iterPrintControl), xform, it0);
  }

  for (int it = 0; it < iters; ++it) {
    int gi = it0 + it;
    std::fill(gMuAcc.begin(), gMuAcc.end(), 0.0);
    std::fill(gScaleAcc.begin(), gScaleAcc.end(), 0.0);
    std::fill(gmPopAcc.begin(), gmPopAcc.end(), 0.0);
    std::fill(gLpopAcc.begin(), gLpopAcc.end(), 0.0);
    double elboAcc = 0.0;
    for (int m = 0; m < nMc; ++m) {
      setSeedEng1((uint32_t)seed + (uint32_t)((gi * nMc + m) & 0x7fffffff));
      // population draw first, then the per-subject eta noise (one stream)
      for (int j = 0; j < npop; ++j) epsPop[j] = rxNormEng(0.0, 1.0);
      for (int i = 0; i < N; ++i)
        for (int k = 0; k < neta; ++k) eps(i, k) = rxNormEng(0.0, 1.0);
      // phi = mPop + Lpop eps_pop  -> map into theta / logPopOmega
      NumericVector phi(npop);
      for (int i = 0; i < npop; ++i) {
        double v = mPop[i];
        for (int j = 0; j <= i; ++j) v += Lpop[i * (i + 1) / 2 + j] * epsPop[j];
        phi[i] = v;
      }
      for (int j = 0; j < npop; ++j) {
        if (phiThetaIdx[j] >= 0) theta[phiThetaIdx[j]] = phi[j];
        else if (phiOmIdx[j] >= 0) logPopOmega[phiOmIdx[j]] = phi[j];
      }
      double e = fr
        ? adviElboGradCoreFR(mu, scale, theta, logPopOmega, eps, muRefThetaIdx, gMu, gScale, gTheta, gPopLogOmega, cores)
        : adviElboGradCore  (mu, scale, theta, logPopOmega, eps, muRefThetaIdx, gMu, gScale, gTheta, gPopLogOmega, cores);
      // population entropy: sum log|Lpop_kk|
      for (int k = 0; k < npop; ++k) e += std::log(std::fabs(Lpop[k * (k + 1) / 2 + k]));
      elboAcc += e;
      for (int j = 0; j < N * neta; ++j) gMuAcc[j] += gMu[j];
      for (int j = 0; j < N * nScale; ++j) gScaleAcc[j] += gScale[j];
      // grad_phi (component j) = gTheta[.] or gPopLogOmega[.]
      NumericVector gPhi(npop);
      for (int j = 0; j < npop; ++j)
        gPhi[j] = (phiThetaIdx[j] >= 0) ? gTheta[phiThetaIdx[j]] : gPopLogOmega[phiOmIdx[j]];
      for (int j = 0; j < npop; ++j) gmPopAcc[j] += gPhi[j];
      for (int i = 0; i < npop; ++i)
        for (int j = 0; j <= i; ++j) gLpopAcc[i * (i + 1) / 2 + j] += gPhi[i] * epsPop[j];
    }
    elboAcc *= invM;
    for (int j = 0; j < N * neta; ++j) gMuAcc[j] *= invM;
    for (int j = 0; j < N * nScale; ++j) gScaleAcc[j] *= invM;
    for (int j = 0; j < npop; ++j) gmPopAcc[j] *= invM;
    for (int j = 0; j < nLpop; ++j) gLpopAcc[j] *= invM;
    // population entropy gradient: +1/Lpop_kk on the diagonal
    for (int k = 0; k < npop; ++k) gLpopAcc[k * (k + 1) / 2 + k] += 1.0 / Lpop[k * (k + 1) / 2 + k];
    elboTrace[it] = elboAcc;

    double idecay = std::pow((double)(gi + 1), -0.5 + eeps);
    for (int j = 0; j < N * neta; ++j) {
      double g = gMuAcc[j];
      sMu[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sMu[j];
      mu[j] += etaScale * idecay / (tau + std::sqrt(sMu[j])) * g;
    }
    for (int j = 0; j < N * nScale; ++j) {
      double g = gScaleAcc[j];
      sScale[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sScale[j];
      scale[j] += etaScale * idecay / (tau + std::sqrt(sScale[j])) * g;
    }
    for (int j = 0; j < npop; ++j) {
      double g = gmPopAcc[j];
      smPop[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * smPop[j];
      mPop[j] += etaScale * idecay / (tau + std::sqrt(smPop[j])) * g;
    }
    for (int j = 0; j < nLpop; ++j) {
      double g = gLpopAcc[j];
      sLpop[j] = (gi == 0) ? g * g : alpha * g * g + (1.0 - alpha) * sLpop[j];
      Lpop[j] += etaScale * idecay / (tau + std::sqrt(sLpop[j])) * g;
    }
    // recenter mu-referenced etas into their phi (theta) posterior means
    for (int j = 0; j < npop; ++j) {
      int k = phiMuRef[j];
      if (k < 0) continue;
      double mbar = 0.0;
      for (int i = 0; i < N; ++i) mbar += mu(i, k);
      mbar /= N;
      mPop[j] += mbar;
      for (int i = 0; i < N; ++i) mu(i, k) -= mbar;
    }
    // record theta / popOmega at the population posterior mean mPop
    NumericVector thetaBar = clone(theta), lpoBar = clone(logPopOmega);
    for (int j = 0; j < npop; ++j) {
      if (phiThetaIdx[j] >= 0) thetaBar[phiThetaIdx[j]] = mPop[j];
      else if (phiOmIdx[j] >= 0) lpoBar[phiOmIdx[j]] = mPop[j];
    }
    for (int p = 0; p < ntheta; ++p) parHist(it, p) = thetaBar[p];
    for (int k = 0; k < neta; ++k) parHist(it, ntheta + k) = std::exp(lpoBar[k]);
    if (doIp) adviIterPrintRow(parHist, it, ntheta + neta, elboAcc, ipPhase);
    if (divergeStop != 0 && adviDiverged(elboAcc, parHist, it, ntheta + neta)) {
      itRun = it + 1; break;
    }
  }
  // final theta / logPopOmega at the population posterior mean
  for (int j = 0; j < npop; ++j) {
    if (phiThetaIdx[j] >= 0) theta[phiThetaIdx[j]] = mPop[j];
    else if (phiOmIdx[j] >= 0) logPopOmega[phiOmIdx[j]] = mPop[j];
  }
  RObject parHistData = (doIp && ipEnd) ? vaeIterPrintGet_(_vaeScale.every != 0) : RObject(R_NilValue);
  return List::create(_["mu"] = mu, _["scale"] = scale, _["theta"] = theta,
                      _["logPopOmega"] = logPopOmega, _["mPop"] = mPop, _["Lpop"] = Lpop,
                      _["elbo"] = elboTrace, _["parHist"] = parHist, _["it0"] = it0 + itRun,
                      _["itRun"] = itRun,
                      _["sMu"] = sMu, _["sScale"] = sScale, _["smPop"] = smPop, _["sLpop"] = sLpop,
                      _["parHistData"] = parHistData);
}

// Format eta like R's paste0("srch ", signif(x, 3)) for the stage label.
static std::string adviSignif3(double x) {
  if (x == 0.0) return "0";
  double m = std::pow(10.0, 2.0 - std::floor(std::log10(std::fabs(x))));
  double r = std::round(x * m) / m;
  char buf[32];
  snprintf(buf, sizeof(buf), "%g", r);
  return std::string(buf);
}

// Whole-optimization driver: the full port of R's .adviOptimize core.  Builds
// the initial variational/population state (or unpacks a resumed one), the
// mu-ref recentering map and the full-Bayes phi maps, chooses the step-size
// scale eta -- a resumed value, an adaptive search over the candidates (paper
// Sec 2.6: short deterministic runs from the initial state, best late-iteration
// mean ELBO, divergent candidates rejected; the old .adviAdaptEta) or the fixed
// fallback -- runs the main loop, and installs the derived result fields
// (popOmega, normalized scale/sScale, the full-Bayes adviCov).  One .Call for
// the whole ADVI optimization; the R wrapper keeps only data prep, the inner
// setup, and print-name metadata.  The search runs print as "srch <eta>" stages
// and the main run as "SGA" in one continuous table.
//[[Rcpp::export]]
List adviOptimize_(List args) {
  const bool pe = as<bool>(args["pointEstimate"]);
  const int frInt = as<int>(args["fr"]);
  const int N = as<int>(args["N"]);
  NumericVector thetaIni = as<NumericVector>(args["theta"]);
  NumericVector omegaIni = as<NumericVector>(args["omega"]);
  IntegerVector muRefThetaIdx = as<IntegerVector>(args["muRefThetaIdx"]);
  LogicalVector thetaFix = as<LogicalVector>(args["thetaFix"]);
  LogicalVector omegaFix = as<LogicalVector>(args["omegaFix"]);
  const int iters = as<int>(args["iters"]);
  const double tau = as<double>(args["tau"]);
  const double alpha = as<double>(args["alpha"]);
  const int nMc = as<int>(args["nMc"]);
  const int cores = as<int>(args["cores"]);
  const bool adapt = as<bool>(args["adaptEta"]);
  NumericVector cands = as<NumericVector>(args["etaCandidates"]);
  const int nAdapt = as<int>(args["nAdapt"]);
  CharacterVector parNames = as<CharacterVector>(args["parNames"]);
  RObject ipc = args["iterPrintControl"];
  RObject xform = args["xform"];
  RObject resumeR = args["resume"];
  const bool hasResume = !Rf_isNull(resumeR);
  List resume;
  if (hasResume) resume = as<List>(resumeR);
  const int ntheta = (int)thetaIni.size(), neta = (int)omegaIni.size();
  const int nL = neta * (neta + 1) / 2;

  // a resumed run continues its original counter-based stream and step-size scale
  double seed = as<double>(args["seed"]);
  double etaResume = NA_REAL;
  if (hasResume) {
    if (resume.containsElementNamed("seed") && !Rf_isNull(resume["seed"]))
      seed = as<double>(resume["seed"]);
    if (resume.containsElementNamed("etaScale") && !Rf_isNull(resume["etaScale"]))
      etaResume = as<double>(resume["etaScale"]);
  }

  // ---- initial variational/population state (fresh) or resume unpack ----
  NumericMatrix mu0, scale0, sMu0, sScale0;
  NumericVector theta0, logPopOmega0, sTheta0, sLpo0;
  int it0 = 0;
  if (!hasResume) {
    theta0 = clone(thetaIni);
    logPopOmega0 = NumericVector(neta);
    for (int k = 0; k < neta; ++k) logPopOmega0[k] = std::log(omegaIni[k]);
    mu0 = NumericMatrix(N, neta);
    sMu0 = NumericMatrix(N, neta);
    sTheta0 = NumericVector(ntheta);
    sLpo0 = NumericVector(neta);
    if (frInt) {
      // full-rank: L_i starts diagonal with L_kk = sqrt(popOmega_k) (packed)
      scale0 = NumericMatrix(N, nL);
      for (int k = 0; k < neta; ++k) {
        double d = std::exp(0.5 * logPopOmega0[k]);
        int col = (k + 1) * (k + 2) / 2 - 1;
        for (int i = 0; i < N; ++i) scale0(i, col) = d;
      }
      sScale0 = NumericMatrix(N, nL);
    } else {
      // mean-field: log-sd starts at the prior scale
      scale0 = NumericMatrix(N, neta);
      for (int k = 0; k < neta; ++k)
        for (int i = 0; i < N; ++i) scale0(i, k) = 0.5 * logPopOmega0[k];
      sScale0 = NumericMatrix(N, neta);
    }
  } else {
    mu0 = as<NumericMatrix>(resume["mu"]);
    theta0 = as<NumericVector>(resume["theta"]);
    logPopOmega0 = as<NumericVector>(resume["logPopOmega"]);
    it0 = as<int>(resume["it0"]);
    sMu0 = as<NumericMatrix>(resume["sMu"]);
    // full-Bayes stores a generic $scale; point-estimate stores $Lpack/$omega
    if (resume.containsElementNamed("scale") && !Rf_isNull(resume["scale"]))
      scale0 = as<NumericMatrix>(resume["scale"]);
    else if (frInt) scale0 = as<NumericMatrix>(resume["Lpack"]);
    else scale0 = as<NumericMatrix>(resume["omega"]);
    sScale0 = as<NumericMatrix>(resume["sScale"]);
    if (pe) {
      sTheta0 = as<NumericVector>(resume["sTheta"]);
      sLpo0 = as<NumericVector>(resume["sLpo"]);
    }
  }

  // per-theta recentering eta (0-based) for mu-referenced intercepts, else -1
  IntegerVector thetaMuRefEta(ntheta, -1);
  for (int k = 0; k < neta; ++k) {
    int p = muRefThetaIdx[k];
    if (p != NA_INTEGER) thetaMuRefEta[p - 1] = k;
  }

  // ---- full-Bayes: phi = c(theta[free], logPopOmega[free]) maps + state ----
  NumericVector mPop0, Lpop0, smPop0, sLpop0;
  IntegerVector phiThetaIdx, phiOmIdx, phiMuRef;
  int npop = 0;
  if (!pe) {
    std::vector<int> thFree, omFree;
    for (int p = 0; p < ntheta; ++p) if (!thetaFix[p]) thFree.push_back(p);
    for (int k = 0; k < neta; ++k) if (!omegaFix[k]) omFree.push_back(k);
    const int nth = (int)thFree.size();
    npop = nth + (int)omFree.size();
    phiThetaIdx = IntegerVector(npop, -1);
    phiOmIdx = IntegerVector(npop, -1);
    phiMuRef = IntegerVector(npop, -1);
    for (int j = 0; j < nth; ++j) {
      phiThetaIdx[j] = thFree[j];
      phiMuRef[j] = thetaMuRefEta[thFree[j]];
    }
    for (int j = 0; j < (int)omFree.size(); ++j) phiOmIdx[nth + j] = omFree[j];
    if (!hasResume) {
      const int nLpop = npop * (npop + 1) / 2;
      mPop0 = NumericVector(npop);
      for (int j = 0; j < nth; ++j) mPop0[j] = theta0[thFree[j]];
      for (int j = 0; j < (int)omFree.size(); ++j) mPop0[nth + j] = logPopOmega0[omFree[j]];
      Lpop0 = NumericVector(nLpop);
      for (int k = 1; k <= npop; ++k) Lpop0[k * (k + 1) / 2 - 1] = 0.1; // init pop posterior sd
      smPop0 = NumericVector(npop);
      sLpop0 = NumericVector(nLpop);
    } else {
      mPop0 = as<NumericVector>(resume["mPop"]);
      Lpop0 = as<NumericVector>(resume["Lpop"]);
      smPop0 = as<NumericVector>(resume["smPop"]);
      sLpop0 = as<NumericVector>(resume["sLpop"]);
    }
  }

  bool ipStarted = false;
  auto run1 = [&](double eta, int itersRun, int itRun0, int divergeStop,
                  const std::string &phase, int ipEnd) -> List {
    int ipStart = ipStarted ? 0 : 1;
    ipStarted = true;
    if (!pe) {
      return adviLoopFB_(mu0, scale0, theta0, logPopOmega0, mPop0, Lpop0,
                         phiThetaIdx, phiOmIdx, phiMuRef, muRefThetaIdx, frInt,
                         itersRun, seed, eta, tau, alpha, nMc, itRun0,
                         sMu0, sScale0, smPop0, sLpop0, cores, divergeStop,
                         parNames, ipc, xform, phase, ipStart, ipEnd);
    }
    if (frInt) {
      return adviLoopFR_(mu0, scale0, theta0, logPopOmega0, muRefThetaIdx,
                         thetaMuRefEta, thetaFix, omegaFix,
                         itersRun, seed, eta, tau, alpha, nMc, itRun0,
                         sMu0, sScale0, sTheta0, sLpo0, cores, divergeStop,
                         parNames, ipc, xform, phase, ipStart, ipEnd);
    }
    return adviLoop_(mu0, scale0, theta0, logPopOmega0, muRefThetaIdx,
                     thetaMuRefEta, thetaFix, omegaFix,
                     itersRun, seed, eta, tau, alpha, nMc, itRun0,
                     sMu0, sScale0, sTheta0, sLpo0, cores, divergeStop,
                     parNames, ipc, xform, phase, ipStart, ipEnd);
  };

  double eta;
  if (R_finite(etaResume)) {
    eta = etaResume;                          // resumed run keeps its step-size scale
  } else if (adapt && cands.size() > 1) {
    double best = R_NegInf, bestEta = cands[0];
    for (int c = 0; c < (int)cands.size(); ++c) {
      double e = cands[c];
      List r = run1(e, nAdapt, 0, 1, std::string("srch ") + adviSignif3(e), 0);
      NumericVector el = r["elbo"];
      NumericMatrix ph = r["parHist"];
      // reject a candidate that diverges (non-finite, or the recorded population
      // estimates blow up); an early-aborted run (itRun < requested) also diverged
      bool div = as<int>(r["itRun"]) < nAdapt;
      if (!div) {
        for (int i = 0; i < (int)el.size(); ++i)
          if (!R_finite(el[i])) { div = true; break; }
      }
      if (!div) {
        for (int i = 0; i < (int)ph.size(); ++i) {
          double v = ph[i];
          if (!R_finite(v) || std::fabs(v) > 1e4) { div = true; break; }
        }
      }
      double score;
      if (div) {
        score = R_NegInf;
      } else {
        // late-iteration mean ELBO, same window as the old R scorer
        int len = (int)el.size();
        int lo = len - nAdapt / 3;
        if (lo < 1) lo = 1;
        long double s = 0.0;
        for (int i = lo - 1; i < len; ++i) s += el[i];
        score = (double)(s / (len - lo + 1));
      }
      if (score > best) { best = score; bestEta = e; }
    }
    eta = bestEta;
  } else if (cands.size() == 1) {
    eta = cands[0];
  } else {
    eta = 0.1;
  }

  List res = run1(eta, iters, it0, 0, "SGA", 1);
  res["etaScale"] = eta;
  res["seed"] = seed;
  NumericVector lpo = res["logPopOmega"];
  NumericVector popOmega((int)lpo.size());
  for (int k = 0; k < (int)lpo.size(); ++k) popOmega[k] = std::exp(lpo[k]);
  res["popOmega"] = popOmega;
  if (pe) {
    res["pointEstimate"] = true;
    // normalize the per-subject scale field name (Lpack for full-rank, omega else)
    if (frInt) {
      res["scale"] = res["Lpack"];
      res["sScale"] = res["sL"];
    } else {
      res["scale"] = res["omega"];
      res["sScale"] = res["sOmega"];
    }
  } else {
    res["pointEstimate"] = false;
    res["phiThetaIdx"] = phiThetaIdx;
    res["phiOmIdx"] = phiOmIdx;
    // population variational covariance in phi space (Lpop Lpop^T)
    NumericVector Lpop = res["Lpop"];
    NumericMatrix cov(npop, npop);
    for (int i = 0; i < npop; ++i) {
      for (int j = 0; j < npop; ++j) {
        int m = i < j ? i : j;
        long double sAcc = 0.0;
        for (int k = 0; k <= m; ++k)
          sAcc += (long double)Lpop[i * (i + 1) / 2 + k] * (long double)Lpop[j * (j + 1) / 2 + k];
        cov(i, j) = (double)sAcc;
      }
    }
    res["adviCov"] = cov;
  }
  return res;
}

// f-SAEM (Karimi, Lavielle & Moulines 2020) proposal builder: for each physical
// subject, optimize the conditional MAP of the random effects and return that
// MAP together with the FOCEi inner information matrix H = Gamma_i^-1 (the
// proposal precision).  innerOpt1(id, 0) runs the n1qn1 inner optimizer and
// finalizes with LikInner2, which -- with _finalObfCalc set -- stores H into
// op_focei.gH.  H is the no-interaction Jacobian information (paper Eq 17) or
// the interaction/Laplace Hessian (Eq 13), selected by the inner control the
// caller set up.  The proposal lives in eta-space; the SAEM chain adds the
// prior mean (phi = mprior + eta) so the covariance transfers unchanged.
// Serial by design here (the SAEM solve is not the active rx during this call,
// and the validation path wants determinism); id-parallelism is added later.
//[[Rcpp::export]]
List fsaemInnerMap_(int cores) {
  (void)cores;
  const int neta = op_focei.neta;
  if (neta == 0) stop("fsaemInnerMap_ requires a model with random effects");
  rx = getRxSolve_();
  const int nsub = (int)getRxNsub(rx);
  NumericMatrix etaHat(nsub, neta);
  NumericMatrix hess(nsub, neta*neta); // row id = vectorized neta x neta H (col-major)
  IntegerVector ok(nsub);
  bool savedFinal = _finalObfCalc;
  _finalObfCalc = true; // make LikInner2 stash H into op_focei.gH
  for (int id = 0; id < nsub; ++id) {
    ok[id] = innerOpt1(id, 0);
    focei_ind *fInd = &(inds_focei[id]);
    for (int j = 0; j < neta; ++j) etaHat(id, j) = fInd->eta[j];
    if (ok[id]) {
      double *Hid = op_focei.gH + (size_t)id*neta*neta;
      for (int j = 0; j < neta*neta; ++j) hess(id, j) = Hid[j];
    } else {
      for (int j = 0; j < neta*neta; ++j) hess(id, j) = NA_REAL;
    }
  }
  _finalObfCalc = savedFinal;
  return List::create(_["eta"] = etaHat, _["hess"] = hess, _["ok"] = ok);
}

// f-SAEM independent Metropolis-Hastings kernel (paper Alg 2 / Eq 23).  Runs on
// the FOCEi inner allocation: the Gaussian proposal N(eta_hat_i, Gamma_i) is
// independent of the chain state, and the joint target p(y_i, eta_i) is the
// FOCEi inner objective (exp(-likInner0)).  For each chain c and subject id a
// candidate is drawn, accepted/rejected by the exact IMH ratio, and the updated
// state returned.  `etaCur` is chain-major ((nchain*nsub) x neta, row = c*nsub +
// id); `cholGamma` row id is the vectorized lower-triangular L (col-major) with
// Gamma_i = L L'.  The proposal-density normalizers (log|L|, (p/2)log 2pi)
// cancel in the ratio, so only the quadratic forms enter.
// Sign: likInner0 = -log p(y_i, eta_i) + C (the minimized inner objective), so
// log alpha = (fCur - fProp) + (0.5||z||^2 - 0.5||L^-1(eta_cur - eta_hat)||^2).
// Proposal draws use rxode2's threefry engine (setSeedEng1 + rxNormEng), seeded
// per (call, chain, subject) via `streamBase` so the sequence is reproducible
// regardless of thread count -- like est="imp"/"impmap" -- instead of R's global
// RNG.  `mprior` (nsub x neta) is the per-subject prior mean, so the proposed
// individual parameter is phi = mprior + eta; when a bounded parameter (nbd/
// lower/upper on the phi scale, from the theta iniDf, matching the L-BFGS-B
// bounds) is violated the normal draw is repeated up to `nRetry` times, then the
// value is clamped to the boundary it last violated.
//[[Rcpp::export]]
List fsaemImhKernel_(NumericMatrix etaCur, NumericMatrix etaHat,
                     NumericMatrix cholGamma, int nchain, int cores,
                     NumericMatrix mprior, NumericVector lower, NumericVector upper,
                     IntegerVector nbd, double streamBase, int nRetry) {
  const int neta = op_focei.neta;
  if (neta == 0) stop("fsaemImhKernel_ requires a model with random effects");
  const int nsub = etaHat.nrow();
  NumericMatrix etaOut = clone(etaCur);
  IntegerVector nAcc(nsub);
  std::vector<arma::mat> L(nsub), Linv(nsub);
  std::vector<bool> good(nsub, false);
  for (int id = 0; id < nsub; ++id) {
    arma::mat Lid(neta, neta);
    for (int j = 0; j < neta*neta; ++j) Lid(j) = cholGamma(id, j);
    if (Lid.is_finite()) {
      arma::mat Li;
      if (arma::inv(Li, arma::trimatl(Lid))) { L[id] = Lid; Linv[id] = Li; good[id] = true; }
    }
  }
  const bool hasBounds = ((int)nbd.size() == neta);
  const uint32_t base = (uint32_t)streamBase;
  // Pin thread 0 (this kernel is serial) so setSeedEng1/rxNormEng use a fixed
  // engine.  The engine array is already allocated + seeded by rxSetSeed() at the
  // start of the fit (saem_fit.R), so we do NOT call seedEng here -- re-seeding it
  // would disrupt the shared rxode2 engine other fits rely on.
  (void)cores;
  setRxThreadId(0);
  for (int c = 0; c < nchain; ++c) {
    for (int id = 0; id < nsub; ++id) {
      int row = c*nsub + id;
      if (!good[id]) continue; // no valid proposal -> keep current state
      // reproducible per-(call, chain, subject) threefry stream
      setSeedEng1(base + (uint32_t)(c*nsub + id));
      arma::vec ecur(neta), ehat(neta), z(neta);
      for (int j = 0; j < neta; ++j) { ecur(j) = etaOut(row, j); ehat(j) = etaHat(id, j); }
      arma::vec eprop(neta);
      // draw the Gaussian proposal, repeating on bound violation up to nRetry
      bool inBounds = true;
      for (int attempt = 0; attempt <= nRetry; ++attempt) {
        for (int j = 0; j < neta; ++j) z(j) = rxNormEng(0.0, 1.0);
        eprop = ehat + L[id]*z;
        if (!hasBounds) break;
        inBounds = true;
        for (int j = 0; j < neta; ++j) {
          double phi = mprior(id, j) + eprop(j);
          if ((nbd[j] == 1 || nbd[j] == 2) && phi < lower[j]) inBounds = false;
          if ((nbd[j] == 3 || nbd[j] == 2) && phi > upper[j]) inBounds = false;
        }
        if (inBounds) break;
      }
      if (hasBounds && !inBounds) {
        // retries exhausted: clamp each violated component to its boundary
        for (int j = 0; j < neta; ++j) {
          double phi = mprior(id, j) + eprop(j);
          if ((nbd[j] == 1 || nbd[j] == 2) && phi < lower[j]) eprop(j) = lower[j] - mprior(id, j);
          if ((nbd[j] == 3 || nbd[j] == 2) && phi > upper[j]) eprop(j) = upper[j] - mprior(id, j);
        }
      }
      std::vector<double> ec(neta), ep(neta);
      for (int j = 0; j < neta; ++j) { ec[j] = ecur(j); ep[j] = eprop(j); }
      double fCur = likInner0(ec.data(), id);
      double fProp = likInner0(ep.data(), id);
      if (!R_FINITE(fProp)) continue; // reject un-solvable candidate
      arma::vec dCur = Linv[id]*(ecur - ehat);
      double logAlpha = (fCur - fProp) + 0.5*(arma::dot(z, z) - arma::dot(dCur, dCur));
      // uniform(0,1) from the same threefry stream: PHI(standard normal)
      double u = rxUnifEng(0.0, 1.0);
      if (R_FINITE(fCur) ? (std::log(u) < logAlpha) : true) {
        for (int j = 0; j < neta; ++j) etaOut(row, j) = eprop(j);
        nAcc[id]++;
      }
    }
  }
  setRxThreadId(-1);
  return List::create(_["eta"] = etaOut, _["nAcc"] = nAcc, _["nchain"] = nchain);
}

// MAP + proposal covariance (Gamma_i = H_i^-1, lower Cholesky) + iteration-indexed
// IMH sweeps, for an inner that is ALREADY re-parameterized.  Shared by the
// no-covariate step (fsaemStepCpp_) and the covariate step (fsaemMapImhCpp_,
// where R rebuilds the mprior-as-data inner first).
static NumericMatrix fsaemMapImh(NumericMatrix mprior, NumericMatrix etaCur, int nchain,
                                 int nsweep, int cores, NumericVector lower, NumericVector upper,
                                 IntegerVector nbd, double seed, int nRetry, int kiter) {
  List mapL = fsaemInnerMap_(cores);
  NumericMatrix etaHat = mapL["eta"];
  NumericMatrix hess = mapL["hess"];
  IntegerVector ok = mapL["ok"];
  const int neta = op_focei.neta;
  const int nsub = etaHat.nrow();
  // cholGamma row id = vec(lower Cholesky of Gamma_i = H_i^-1); NA where the MAP
  // failed.  R did solve(H) then t(chol(.)); chol reads the upper triangle, so
  // symmatu() reproduces that input and "lower" == t(chol()).
  NumericMatrix cholGamma(nsub, neta * neta);
  for (int id = 0; id < nsub; ++id) {
    bool bad = (ok[id] == 0);
    arma::mat L(neta, neta, arma::fill::zeros);
    if (!bad) {
      arma::mat H(neta, neta);
      for (int j = 0; j < neta * neta; ++j) H(j) = hess(id, j);
      arma::mat G;
      if (!arma::inv(G, H)) bad = true;
      else if (!arma::chol(L, arma::symmatu(G), "lower")) bad = true;
    }
    for (int j = 0; j < neta * neta; ++j) cholGamma(id, j) = bad ? NA_REAL : L(j);
  }
  NumericMatrix eta = clone(etaCur);
  for (int s = 0; s < nsweep; ++s) {
    double streamBase = seed + ((double)kiter * nsweep + s) * ((double)nchain * nsub);
    List r = fsaemImhKernel_(eta, etaHat, cholGamma, nchain, cores,
                             mprior, lower, upper, nbd, streamBase, nRetry);
    eta = as<NumericMatrix>(r["eta"]);
  }
  return eta;
}

// C++-native f-SAEM step: the whole no-covariate per-iteration orchestration that
// was the R closure (.fsaemInnerUpdate + .fsaemInnerMap + .fsaemImh).  Keeping it
// in C++ lets the SAEM loop drive it directly (no per-iteration R round-trip),
// which is what makes it safe to change a phase's step count dynamically.
//[[Rcpp::export]]
NumericMatrix fsaemStepCpp_(Environment env, NumericVector theta, NumericVector omega,
                            NumericMatrix mprior, NumericMatrix etaCur, int nchain,
                            int nsweep, int cores, NumericVector lower, NumericVector upper,
                            IntegerVector nbd, double seed, int nRetry, int kiter) {
  // ---- re-parameterize the inner (mirror .fsaemInnerUpdate) ----
  int nth = theta.size();
  NumericVector th = clone(theta);
  CharacterVector thNames(nth);
  for (int i = 0; i < nth; ++i) thNames[i] = "THETA[" + std::to_string(i + 1) + "]";
  th.attr("names") = thNames;
  env["thetaIni"] = th;
  int no = omega.size();
  NumericMatrix om(no, no);
  for (int i = 0; i < no; ++i) om(i, i) = omega[i];
  Environment rxns = Environment::namespace_env("rxode2");
  Function symInvCreate = rxns["rxSymInvCholCreate"];
  env["rxInv"] = symInvCreate(_["mat"] = om, _["diag.xform"] = "sqrt");
  env["etaMat"] = NumericMatrix(mprior.nrow(), no);   // zeros
  vaeInnerSetup_(env);
  return fsaemMapImh(mprior, etaCur, nchain, nsweep, cores, lower, upper, nbd, seed, nRetry, kiter);
}

// Covariate path: R re-parameterizes the mprior-as-data inner (data rebuild +
// full setup), then this does the C++ MAP + IMH.
//[[Rcpp::export]]
NumericMatrix fsaemMapImhCpp_(NumericMatrix mprior, NumericMatrix etaCur, int nchain,
                              int nsweep, int cores, NumericVector lower, NumericVector upper,
                              IntegerVector nbd, double seed, int nRetry, int kiter) {
  return fsaemMapImh(mprior, etaCur, nchain, nsweep, cores, lower, upper, nbd, seed, nRetry, kiter);
}

//[[Rcpp::export]]
RObject vaeInnerFree_() {
  rxOptionsFreeFocei();
  return R_NilValue;
}


// ---------------------------------------------------------------------------
// VAE training loop (est="vae"), fully in C++.  This is a straight port of the
// former R orchestration (.vaeTrain / .vaeElboStepInner / .vaeMStep* /
// .vaeUpdateErr / .vaeAdamStep in R/vaeFit.R + R/vaeInner.R): burn-in
// (encoder-only, tiny KL) -> main EM (KL anneal + closed-form M-step, optional
// BICc-ELBO covariate selection) -> EMA smoothing.  The heavy pieces it calls
// are already C++: the LSTM encoder fwd/bwd (vaeEncoderFwdBwdCore), the parallel
// FOCEi inner likelihood (vaeInnerLikCore) and its per-step reparameterization
// (vaeInnerUpdateParCore), and the shared iteration-print/parHist machinery
// (vaeIterPrint*_).  Keeping the loop in C++ removes the per-gradient-step R
// round-trip AND the per-step full inner re-setup that made est="vae" too slow.

// Per-block Adam optimizer state for the 6 encoder parameter blocks.
struct VaeAdamBlk { arma::mat m, v; };

// R IntegerVector -> arma::ivec (arma sword is 64-bit; no direct memory ctor).
static inline arma::ivec vaeToIvec(IntegerVector v) {
  arma::ivec o(v.size());
  for (int i = 0; i < v.size(); ++i) o[i] = v[i];
  return o;
}

// One reparameterization-noise draw, matching R's matrix(rxnorm(N*zDim),N,zDim)
// column-major fill.  Uses rxode2's threefry engine seeded per (phase,it,l) so a
// step is reproducible regardless of the total iteration count (like the R
// rxWithSeed(.es, ..., rxseed=.es) scheme, though the engine stream itself
// differs from R's rxnorm -- the fit is validated against simulation truth).
static inline void vaeDrawEps(arma::mat& eps, uint32_t es) {
  setSeedEng1(es);
  double *p = eps.memptr();
  const arma::uword n = eps.n_elem;
  for (arma::uword i = 0; i < n; ++i) p[i] = rxNormEng(0.0, 1.0);
}

static inline arma::vec vaeClampVec(arma::vec v, const arma::vec& lo, const arma::vec& hi) {
  for (arma::uword k = 0; k < v.n_elem; ++k) {
    if (R_FINITE(lo[k]) && v[k] < lo[k]) v[k] = lo[k];
    if (R_FINITE(hi[k]) && v[k] > hi[k]) v[k] = hi[k];
  }
  return v;
}

// Assemble the full theta vector (ntheta order) from the current structural
// population means (zPop, on the mu-referenced thetas) and residual-error a.
// zPopThetaIdx0[k] is the 0-based theta index for eta k, or -1 (free eta:
// theta forced to 0, structure carried elsewhere).  errThetaIdx0 likewise.
static inline arma::vec vaeBuildTh(const arma::vec& th, const arma::ivec& zPopThetaIdx0,
                                   const arma::vec& baseline, const arma::ivec& errThetaIdx0,
                                   const arma::vec& a) {
  arma::vec out = th;
  for (arma::uword k = 0; k < zPopThetaIdx0.n_elem; ++k)
    if (zPopThetaIdx0[k] >= 0) out[zPopThetaIdx0[k]] = baseline[k];
  for (arma::uword e = 0; e < errThetaIdx0.n_elem; ++e)
    out[errThetaIdx0[e]] = a[e];
  return out;
}

// Closed-form error-parameter M-step (add/prop/combined), robust to non-finite
// predictions.  errTypeCode: 0=add, 1=prop, 2=other (kept at current value).
static arma::vec vaeUpdateErr(const std::vector<std::vector<double> >& preds,
                              const std::vector<std::vector<double> >& yList,
                              const arma::ivec& errTypeCode, const arma::vec& a,
                              const arma::vec& errLower, const arma::vec& errUpper) {
  if (a.n_elem == 0) return a;
  std::vector<double> res, f;
  for (size_t i = 0; i < preds.size(); ++i) {
    const std::vector<double>& fi = preds[i];
    const std::vector<double>& yi = yList[i];
    size_t n = std::min(fi.size(), yi.size());
    for (size_t o = 0; o < n; ++o) {
      double r = yi[o] - fi[o], ff = fi[o];
      if (R_FINITE(r) && R_FINITE(ff)) { res.push_back(r); f.push_back(ff); }
    }
  }
  arma::vec aNew = a;
  if (res.empty()) return a;
  arma::uword nAdd = 0, nProp = 0; arma::uword iAdd = 0, iProp = 0;
  for (arma::uword e = 0; e < errTypeCode.n_elem; ++e) {
    if (errTypeCode[e] == 0 && nAdd == 0) { iAdd = e; nAdd = 1; }
    if (errTypeCode[e] == 1 && nProp == 0) { iProp = e; nProp = 1; }
  }
  const size_t m = res.size();
  if (nAdd && nProp) {
    // combined: res^2 ~ a^2 + b^2 f^2 via 2-col least squares (normal equations)
    double s0 = m, s1 = 0, s11 = 0, b0 = 0, b1 = 0;
    for (size_t o = 0; o < m; ++o) {
      double x1 = f[o] * f[o], y = res[o] * res[o];
      s1 += x1; s11 += x1 * x1; b0 += y; b1 += x1 * y;
    }
    double det = s0 * s11 - s1 * s1;
    if (std::fabs(det) > 0) {
      double v0 = (s11 * b0 - s1 * b1) / det;
      double v1 = (s0 * b1 - s1 * b0) / det;
      v0 = std::max(v0, DBL_EPSILON); v1 = std::max(v1, DBL_EPSILON);
      if (R_FINITE(v0)) aNew[iAdd] = std::sqrt(v0);
      if (R_FINITE(v1)) aNew[iProp] = std::sqrt(v1);
    }
  } else if (nAdd) {
    double ss = 0; for (double r : res) ss += r * r;
    aNew[iAdd] = std::sqrt(ss / m);
  } else if (nProp) {
    double ss = 0; size_t nn = 0;
    for (size_t o = 0; o < m; ++o) if (std::fabs(f[o]) > 1e-8) { double q = res[o] / f[o]; ss += q * q; nn++; }
    if (nn) aNew[iProp] = std::sqrt(ss / nn);
  }
  for (arma::uword e = 0; e < aNew.n_elem; ++e) if (!R_FINITE(aNew[e])) aNew[e] = a[e];
  return vaeClampVec(aNew, errLower, errUpper);
}

// Result of one ELBO evaluation.
struct VaeStepOut {
  bool ok;
  double pxz, DKL;
  arma::mat mu, z;      // [N, zDim]
  arma::cube L;         // [zDim, zDim, N]
  arma::mat lp;         // [N, zDim] decoder eta-gradient (best mixture component)
  std::vector<std::vector<double> > preds;  // per-subject predictions f
  arma::ivec mixnum;    // [N] selected mixture component (1-based)
  // encoder parameter gradients
  arma::mat gWih, gWhh, gFcW; arma::vec gbih, gbhh, gFcB;
};

// One ELBO evaluation using the FOCEi inner likelihood (port of
// .vaeElboStepInner).  zPopMat [N,zDim] is the (possibly subject-specific) KL
// center; baseline [zDim] is the inner-problem THETA baseline.
static VaeStepOut vaeElboStepCpp(const arma::mat& Wih, const arma::mat& Whh,
                                 const arma::vec& bih, const arma::vec& bhh,
                                 const arma::mat& fcW, const arma::vec& fcB,
                                 const arma::cube& dataIn, const arma::ivec& lengths,
                                 const arma::mat& covIn, const arma::mat& eps, int zDim,
                                 const arma::vec& th, const arma::ivec& zPopThetaIdx0,
                                 const arma::ivec& errThetaIdx0,
                                 const arma::mat& zPopMat, const arma::vec& baseline,
                                 const arma::vec& omega, const arma::vec& a, double alphaKL,
                                 int nMix, const arma::vec& mixProb, int cores,
                                 bool withGrad = true) {
  VaeStepOut S; S.ok = true;
  const int N = dataIn.n_rows;
  // encoder forward (no backward yet -- need z to form gZ first)
  arma::mat mu(N, zDim), logSigma(N, zDim), zOut(N, zDim);
  arma::cube Lout(zDim, zDim, N);
  arma::mat dummyGz(N, zDim, arma::fill::zeros), dummyGls(N, zDim, arma::fill::zeros);
  arma::mat gWih, gWhh, gFcW; arma::vec gbih, gbhh, gFcB;
  vaeEncoderFwdBwdCore(dataIn, lengths, covIn, eps, Wih, Whh, bih, bhh, fcW, fcB, zDim,
                       dummyGz, dummyGls, false, mu, logSigma, Lout, zOut,
                       gWih, gWhh, gbih, gbhh, gFcW, gFcB);
  arma::mat eta = zOut;
  eta.each_row() -= baseline.t();                       // z - baseline
  // inner-problem re-parameterization + evaluation
  arma::vec thv = vaeBuildTh(th, zPopThetaIdx0, baseline, errThetaIdx0, a);
  vaeInnerUpdateParCore(thv, omega);
  arma::mat etaEval = eta;
  if (nMix > 1) { etaEval.set_size(nMix * N, zDim); for (int m = 0; m < nMix; ++m) etaEval.rows(m * N, m * N + N - 1) = eta; }
  arma::vec obj; arma::mat lp; std::vector<std::vector<double> > pf;
  vaeInnerLikCore(etaEval, cores, true, true, obj, lp, pf);

  const double ln2pi = std::log(2 * M_PI);
  arma::vec pzI(N);
  for (int i = 0; i < N; ++i) {
    double s = 0; for (int k = 0; k < zDim; ++k) s += eta(i, k) * eta(i, k) / omega[k] + std::log(omega[k]) + ln2pi;
    pzI[i] = 0.5 * s;
  }
  double jointTot;
  arma::mat lpBest(N, zDim);
  S.preds.resize(N);
  S.mixnum.set_size(N); S.mixnum.fill(1);
  if (nMix > 1) {
    jointTot = 0;
    for (int i = 0; i < N; ++i) {
      double mmax = -std::numeric_limits<double>::infinity(); int best = 0;
      arma::vec ll(nMix);
      for (int m = 0; m < nMix; ++m) {
        double v = std::log(mixProb[m]) - 0.5 * obj[m * N + i];
        if (!R_FINITE(v)) v = -std::numeric_limits<double>::infinity();
        ll[m] = v; if (v > mmax) { mmax = v; best = m; }
      }
      if (R_FINITE(mmax)) {
        double se = 0; for (int m = 0; m < nMix; ++m) se += std::exp(ll[m] - mmax);
        jointTot += -2 * (mmax + std::log(se));
      }
      lpBest.row(i) = lp.row(best * N + i);
      S.preds[i] = pf[best * N + i];
      S.mixnum[i] = best + 1;
    }
  } else {
    jointTot = arma::accu(obj);
    lpBest = lp;
    for (int i = 0; i < N; ++i) S.preds[i] = pf[i];
  }
  double pxz = jointTot - arma::accu(pzI);
  double pz = 0, qz = 0;
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < zDim; ++k) {
      double d = zOut(i, k) - zPopMat(i, k);
      pz += 0.5 * (d * d / omega[k] + std::log(omega[k]) + ln2pi);
      qz += 0.5 * (eps(i, k) * eps(i, k) + ln2pi + 2 * logSigma(i, k));
    }
  }
  double DKL = pz - qz;
  if (withGrad) {
    // encoder upstream: d(pxz)/dz = lp - eta/omega; KL adds alphaKL*(z-zPopMat)/omega
    arma::mat gZ = lpBest;
    for (int i = 0; i < N; ++i)
      for (int k = 0; k < zDim; ++k) {
        double g = gZ(i, k) - eta(i, k) / omega[k] + alphaKL * (zOut(i, k) - zPopMat(i, k)) / omega[k];
        gZ(i, k) = R_FINITE(g) ? g : 0.0;
      }
    arma::mat gLS(N, zDim); gLS.fill(-alphaKL);
    arma::mat mu2(N, zDim), ls2(N, zDim), z2(N, zDim); arma::cube L2(zDim, zDim, N);
    vaeEncoderFwdBwdCore(dataIn, lengths, covIn, eps, Wih, Whh, bih, bhh, fcW, fcB, zDim,
                         gZ, gLS, true, mu2, ls2, L2, z2,
                         S.gWih, S.gWhh, S.gbih, S.gbhh, S.gFcW, S.gFcB);
  }
  S.pxz = pxz; S.DKL = DKL; S.mu = mu; S.z = zOut; S.L = Lout; S.lp = lpBest;
  if (!R_FINITE(pxz)) S.ok = true; // a failed subject is zeroed in gZ, not fatal
  return S;
}

// R-facing wrapper for one ELBO step (the vaeTrainCpp_ loop calls vaeElboStepCpp
// directly; this export lets R unit-test the same core).  The FOCEi inner problem
// must ALREADY be set up (.vaeInnerSetup / vaeInnerSetup_) -- this reuses the
// active op_focei allocation, exactly as the training loop does.  `prep` is the
// .vaeDataPrep output; `zPop` is either a length-zDim vector (shared KL center) or
// an N x zDim matrix (subject-specific covariate centers).  Returns the same list
// the former R .vaeElboStepInner produced.
//[[Rcpp::export]]
List vaeElboStepCpp_(List params, List prep, RObject zPopR, NumericVector omegaR,
                     NumericVector aR, double alphaKL, NumericMatrix epsR,
                     int nMix, NumericVector mixProbR, int cores, bool withGrad = true) {
  arma::mat Wih = as<arma::mat>(params["Wih"]);
  arma::mat Whh = as<arma::mat>(params["Whh"]);
  arma::vec bih = as<arma::vec>(params["bih"]);
  arma::vec bhh = as<arma::vec>(params["bhh"]);
  arma::mat fcW = as<arma::mat>(params["fcW"]);
  arma::vec fcB = as<arma::vec>(params["fcB"]);
  const int N = as<int>(prep["N"]);
  const int zDim = as<int>(prep["zDim"]);
  arma::cube dataIn = as<arma::cube>(prep["dataIn"]);
  arma::ivec lengths = vaeToIvec(prep["lengths"]);
  arma::mat covIn = as<arma::mat>(prep["covIn"]);
  arma::vec th = as<arma::vec>(prep["th"]);
  // .vaeDataPrep stores 1-based theta indices with NA for a free/mixture eta; map
  // to the 0-based (-1 = free) form the core uses.
  IntegerVector zpi = prep["zPopThetaIdx"];
  arma::ivec zPopThetaIdx0(zpi.size());
  for (int k = 0; k < zpi.size(); ++k)
    zPopThetaIdx0[k] = IntegerVector::is_na(zpi[k]) ? -1 : (int)zpi[k] - 1;
  IntegerVector eti = prep["errThetaIdx"];
  arma::ivec errThetaIdx0(eti.size());
  for (int e = 0; e < eti.size(); ++e) errThetaIdx0[e] = (int)eti[e] - 1;
  arma::vec omega(omegaR.begin(), omegaR.size());
  arma::vec a(aR.begin(), aR.size());
  arma::mat eps(epsR.begin(), N, zDim, false, true);
  arma::vec mixProb(mixProbR.begin(), mixProbR.size());
  // zPop -> KL center matrix + inner-problem baseline (mirrors the R is.matrix branch)
  arma::mat zPopMat(N, zDim); arma::vec baseline;
  if (Rf_isMatrix(zPopR)) {
    zPopMat = as<arma::mat>(zPopR);
    baseline = arma::mean(zPopMat, 0).t();
  } else {
    arma::vec zp = as<arma::vec>(zPopR);
    zPopMat.each_row() = zp.t();
    baseline = zp;
  }
  VaeStepOut S = vaeElboStepCpp(Wih, Whh, bih, bhh, fcW, fcB, dataIn, lengths, covIn,
                                eps, zDim, th, zPopThetaIdx0, errThetaIdx0, zPopMat, baseline,
                                omega, a, alphaKL, nMix, mixProb, cores, withGrad);
  RObject grads = R_NilValue;
  if (withGrad) grads = List::create(_["Wih"] = S.gWih, _["Whh"] = S.gWhh, _["bih"] = S.gbih,
                                     _["bhh"] = S.gbhh, _["fcW"] = S.gFcW, _["fcB"] = S.gFcB);
  List preds(N);
  for (int i = 0; i < N; ++i) preds[i] = NumericVector(S.preds[i].begin(), S.preds[i].end());
  IntegerVector mixnum(N);
  for (int i = 0; i < N; ++i) mixnum[i] = S.mixnum[i];
  return List::create(_["loss"] = S.pxz + alphaKL * S.DKL, _["pxz"] = S.pxz, _["DKL"] = S.DKL,
                      _["grads"] = grads, _["mu"] = S.mu, _["L"] = S.L, _["z"] = S.z,
                      _["preds"] = preds, _["mixnum"] = mixnum);
}

// One Adam update over a block (m,v are the block's moment estimates).
static inline void vaeAdam(arma::mat& p, const arma::mat& g, VaeAdamBlk& st,
                           double lr, int t, double b1 = 0.9, double b2 = 0.999, double eps = 1e-8) {
  st.m = b1 * st.m + (1 - b1) * g;
  st.v = b2 * st.v + (1 - b2) * (g % g);
  arma::mat mhat = st.m / (1 - std::pow(b1, t));
  arma::mat vhat = st.v / (1 - std::pow(b2, t));
  p -= lr * (mhat / (arma::sqrt(vhat) + eps));
}

//[[Rcpp::export]]
List vaeTrainCpp_(List params, List prep, List control, int nMix, NumericVector mixProbR,
                  int cores, NumericVector row0, CharacterVector parNames,
                  List iterPrintControl, RObject xform, IntegerVector structIdx0) {
  // ---- unpack encoder parameters (mutated in place by Adam) ----
  arma::mat Wih = as<arma::mat>(params["Wih"]);
  arma::mat Whh = as<arma::mat>(params["Whh"]);
  arma::vec bih = as<arma::vec>(params["bih"]);
  arma::vec bhh = as<arma::vec>(params["bhh"]);
  arma::mat fcW = as<arma::mat>(params["fcW"]);
  arma::vec fcB = as<arma::vec>(params["fcB"]);
  // ---- unpack prep ----
  const int N = as<int>(prep["N"]);
  const int zDim = as<int>(prep["zDim"]);
  arma::cube dataIn = as<arma::cube>(prep["dataIn"]);
  arma::ivec lengths = vaeToIvec(prep["lengths"]);
  arma::mat covIn = as<arma::mat>(prep["covIn"]);
  arma::mat covMat = as<arma::mat>(prep["covMat"]);
  arma::vec th = as<arma::vec>(prep["th"]);
  arma::ivec zPopThetaIdx0 = vaeToIvec(prep["zPopThetaIdx0"]);
  arma::ivec errThetaIdx0 = vaeToIvec(prep["errThetaIdx0"]);
  LogicalVector isFreeR = prep["isFree"];
  LogicalVector omegaFixR = prep["omegaFix"];
  arma::vec zPop = as<arma::vec>(prep["zPop"]);
  arma::vec omega = as<arma::vec>(prep["omega"]);
  arma::vec a = as<arma::vec>(prep["a"]);
  arma::vec zPopLower = as<arma::vec>(prep["zPopLower"]);
  arma::vec zPopUpper = as<arma::vec>(prep["zPopUpper"]);
  arma::vec errLower = as<arma::vec>(prep["errLower"]);
  arma::vec errUpper = as<arma::vec>(prep["errUpper"]);
  arma::ivec errTypeCode = vaeToIvec(prep["errTypeCode"]);
  List yListR = prep["yList"];
  std::vector<std::vector<double> > yList(N);
  for (int i = 0; i < N; ++i) { NumericVector yi = yListR[i]; yList[i].assign(yi.begin(), yi.end()); }
  // ---- control ----
  const int seed = as<int>(control["seed"]);
  const int itersBurnIn = as<int>(control["itersBurnIn"]);
  const int klWarmup = as<int>(control["klWarmup"]);
  const int gammaIter = as<int>(control["gammaIter"]);
  const int iters = as<int>(control["iters"]);
  const int Lg = as<int>(control["nGradStep"]);
  const double learningRate = as<double>(control["learningRate"]);
  const double burnInLearningRate = as<double>(control["burnInLearningRate"]);
  const bool covariateSelection = as<bool>(control["covariateSelection"]);
  const int printCtl = as<int>(control["print"]);
  arma::vec mixProb(mixProbR.begin(), mixProbR.size());
  const int nCov = covMat.n_cols;

  // ---- Adam state (zeros, per block) ----
  VaeAdamBlk aWih{arma::zeros(Wih.n_rows, Wih.n_cols), arma::zeros(Wih.n_rows, Wih.n_cols)};
  VaeAdamBlk aWhh{arma::zeros(Whh.n_rows, Whh.n_cols), arma::zeros(Whh.n_rows, Whh.n_cols)};
  VaeAdamBlk aBih{arma::zeros(bih.n_elem, 1), arma::zeros(bih.n_elem, 1)};
  VaeAdamBlk aBhh{arma::zeros(bhh.n_elem, 1), arma::zeros(bhh.n_elem, 1)};
  VaeAdamBlk aFcW{arma::zeros(fcW.n_rows, fcW.n_cols), arma::zeros(fcW.n_rows, fcW.n_cols)};
  VaeAdamBlk aFcB{arma::zeros(fcB.n_elem, 1), arma::zeros(fcB.n_elem, 1)};
  arma::mat bihM(bih), bhhM(bhh), fcBM(fcB); // matrix views for the vector Adam blocks

  // ---- parameter-history row helper ----
  arma::ivec structIdx = vaeToIvec(structIdx0);
  auto parRow = [&](const arma::vec& zP, const arma::vec& om, const arma::vec& av) {
    NumericVector r(structIdx.n_elem + zDim + av.n_elem);
    int p = 0;
    for (arma::uword s = 0; s < structIdx.n_elem; ++s) r[p++] = zP[structIdx[s]];
    for (int k = 0; k < zDim; ++k) r[p++] = om[k];
    for (arma::uword e = 0; e < av.n_elem; ++e) r[p++] = av[e];
    return r;
  };

  vaeIterPrintStart_(row0, parNames, iterPrintControl, xform);

  setRxThreadId(0);
  int tstep = 0;
  VaeStepOut last;

  // ---- burn-in: encoder-only, tiny fixed KL weight ----
  for (int it = 1; it <= itersBurnIn; ++it) {
    for (int l = 1; l <= Lg; ++l) {
      arma::mat eps(N, zDim);
      vaeDrawEps(eps, (uint32_t)(seed + (it - 1) * Lg + l));
      arma::mat zPopMat(N, zDim); zPopMat.each_row() = zPop.t();
      VaeStepOut st = vaeElboStepCpp(Wih, Whh, bih, bhh, fcW, fcB, dataIn, lengths, covIn,
                                     eps, zDim, th, zPopThetaIdx0, errThetaIdx0, zPopMat, zPop,
                                     omega, a, 0.001, nMix, mixProb, cores);
      tstep++;
      bihM = bih; bhhM = bhh; fcBM = fcB;
      vaeAdam(Wih, st.gWih, aWih, burnInLearningRate, tstep);
      vaeAdam(Whh, st.gWhh, aWhh, burnInLearningRate, tstep);
      vaeAdam(bihM, st.gbih, aBih, burnInLearningRate, tstep); bih = bihM.col(0);
      vaeAdam(bhhM, st.gbhh, aBhh, burnInLearningRate, tstep); bhh = bhhM.col(0);
      vaeAdam(fcW, st.gFcW, aFcW, burnInLearningRate, tstep);
      vaeAdam(fcBM, st.gFcB, aFcB, burnInLearningRate, tstep); fcB = fcBM.col(0);
      // M-step (gamma = 1)
      arma::vec zPopCur = arma::mean(st.mu, 0).t();
      arma::vec omegaCur(zDim);
      for (int k = 0; k < zDim; ++k) {
        double v = 0;
        for (int i = 0; i < N; ++i) { double d = st.mu(i, k) - zPopCur[k]; v += d * d + arma::accu(arma::square(st.L.slice(i).row(k))); }
        omegaCur[k] = v / N;
      }
      arma::vec aCur = vaeUpdateErr(st.preds, yList, errTypeCode, a, errLower, errUpper);
      for (int k = 0; k < zDim; ++k) { if (isFreeR[k]) zPopCur[k] = 0; if (omegaFixR[k]) omegaCur[k] = omega[k]; }
      if (!zPopCur.is_finite()) zPopCur = zPop;
      if (!omegaCur.is_finite()) omegaCur = omega;
      zPopCur = vaeClampVec(zPopCur, zPopLower, zPopUpper);
      zPop = zPopCur; omega = omegaCur; a = aCur;  // gamma = 1
      last = st;
    }
    vaeIterPrintRow_(parRow(zPop, omega, a), last.pxz + last.DKL, "Burn in");
    Rcpp::checkUserInterrupt();
  }

  // ---- main EM (optionally BICc-ELBO covariate selection) ----
  const bool doCov = covariateSelection && nCov > 0;
  arma::vec intercept = zPop;
  arma::mat beta(zDim, nCov, arma::fill::zeros);
  arma::umat selected(zDim, nCov, arma::fill::zeros);
  arma::mat zPopArg(N, zDim); zPopArg.each_row() = zPop.t();
  bool isCovStep = false;
  arma::vec elboTrace(iters, arma::fill::zeros);
  const double logN = std::log((double)N);

  for (int it = 1; it <= iters; ++it) {
    double gamma = (it <= gammaIter) ? 1.0 : 1.0 / (1.0 + it - gammaIter);
    arma::vec baseline;
    if (doCov) {
      // BICc-ELBO covariate M-step
      arma::mat X(N, 1 + nCov); X.col(0).ones(); if (nCov > 0) X.cols(1, nCov) = covMat;
      arma::mat zPopMat(N, zDim, arma::fill::zeros);
      intercept.zeros(); beta.zeros(); selected.zeros();
      for (int k = 0; k < zDim; ++k) {
        if (isFreeR[k]) continue;
        arma::vec yk = last.mu.col(k);
        double bestScore = std::numeric_limits<double>::infinity();
        arma::uvec bestCols; arma::vec bestCoef; arma::uvec bestS;
        for (unsigned int mask = 0; mask < (1u << nCov); ++mask) {
          std::vector<unsigned> Sset; for (int b = 0; b < nCov; ++b) if (mask & (1u << b)) Sset.push_back(b);
          arma::uvec cols(1 + Sset.size()); cols[0] = 0; for (size_t s = 0; s < Sset.size(); ++s) cols[s + 1] = 1 + Sset[s];
          arma::mat Xs = X.cols(cols);
          arma::vec coef;
          if (!arma::solve(coef, Xs, yk)) { if (!arma::solve(coef, Xs, yk, arma::solve_opts::force_approx)) continue; }
          double rss = arma::accu(arma::square(yk - Xs * coef));
          double score = rss / omega[k] + logN * (double)Sset.size();
          if (score < bestScore) { bestScore = score; bestCols = cols; bestCoef = coef; bestS = arma::uvec(Sset.size()); for (size_t s = 0; s < Sset.size(); ++s) bestS[s] = Sset[s]; }
        }
        double ic = bestCoef[0];
        if (R_FINITE(zPopLower[k]) && ic < zPopLower[k]) ic = zPopLower[k];
        if (R_FINITE(zPopUpper[k]) && ic > zPopUpper[k]) ic = zPopUpper[k];
        intercept[k] = ic; bestCoef[0] = ic;
        for (arma::uword s = 0; s < bestS.n_elem; ++s) { beta(k, bestS[s]) = bestCoef[s + 1]; selected(k, bestS[s]) = 1; }
        zPopMat.col(k) = X.cols(bestCols) * bestCoef;
      }
      arma::vec omegaCur(zDim);
      for (int k = 0; k < zDim; ++k) {
        double v = 0;
        for (int i = 0; i < N; ++i) { double d = last.mu(i, k) - zPopMat(i, k); v += d * d + arma::accu(arma::square(last.L.slice(i).row(k))); }
        omegaCur[k] = v / N;
      }
      arma::vec aCur = vaeUpdateErr(last.preds, yList, errTypeCode, a, errLower, errUpper);
      for (int k = 0; k < zDim; ++k) if (omegaFixR[k]) omegaCur[k] = omega[k];
      if (!omegaCur.is_finite()) omegaCur = omega;
      zPop = intercept;
      zPopArg = zPopMat;
      omega = omega + gamma * (omegaCur - omega);
      a = a + gamma * (aCur - a);
      isCovStep = true;
      baseline = arma::mean(zPopMat, 0).t();
    } else {
      // plain closed-form M-step (EMA gamma)
      arma::vec zPopCur = arma::mean(last.mu, 0).t();
      arma::vec omegaCur(zDim);
      for (int k = 0; k < zDim; ++k) {
        double v = 0;
        for (int i = 0; i < N; ++i) { double d = last.mu(i, k) - zPopCur[k]; v += d * d + arma::accu(arma::square(last.L.slice(i).row(k))); }
        omegaCur[k] = v / N;
      }
      arma::vec aCur = vaeUpdateErr(last.preds, yList, errTypeCode, a, errLower, errUpper);
      for (int k = 0; k < zDim; ++k) { if (isFreeR[k]) zPopCur[k] = 0; if (omegaFixR[k]) omegaCur[k] = omega[k]; }
      if (!zPopCur.is_finite()) zPopCur = zPop;
      if (!omegaCur.is_finite()) omegaCur = omega;
      zPopCur = vaeClampVec(zPopCur, zPopLower, zPopUpper);
      zPop = zPop + gamma * (zPopCur - zPop);
      omega = omega + gamma * (omegaCur - omega);
      a = a + gamma * (aCur - a);
      zPopArg.each_row() = zPop.t();
      isCovStep = false;
      baseline = zPop;
    }
    double alphaKL = (it <= klWarmup) ? 0.01 + 0.99 * (it - 1) / std::max(1, klWarmup - 1) : 1.0;
    double esum = 0;
    for (int l = 1; l <= Lg; ++l) {
      arma::mat eps(N, zDim);
      vaeDrawEps(eps, (uint32_t)(seed + 1000003 + (it - 1) * Lg + l));
      VaeStepOut st = vaeElboStepCpp(Wih, Whh, bih, bhh, fcW, fcB, dataIn, lengths, covIn,
                                     eps, zDim, th, zPopThetaIdx0, errThetaIdx0, zPopArg, baseline,
                                     omega, a, alphaKL, nMix, mixProb, cores);
      tstep++;
      bihM = bih; bhhM = bhh; fcBM = fcB;
      vaeAdam(Wih, st.gWih, aWih, learningRate, tstep);
      vaeAdam(Whh, st.gWhh, aWhh, learningRate, tstep);
      vaeAdam(bihM, st.gbih, aBih, learningRate, tstep); bih = bihM.col(0);
      vaeAdam(bhhM, st.gbhh, aBhh, learningRate, tstep); bhh = bhhM.col(0);
      vaeAdam(fcW, st.gFcW, aFcW, learningRate, tstep);
      vaeAdam(fcBM, st.gFcB, aFcB, learningRate, tstep); fcB = fcBM.col(0);
      esum += st.pxz + st.DKL;
      last = st;
    }
    elboTrace[it - 1] = esum / Lg;
    const char* phase = (it <= klWarmup) ? "KL anneal" : (it <= gammaIter) ? "EM" : "Smooth";
    vaeIterPrintRow_(parRow(zPop, omega, a), elboTrace[it - 1], phase);
    Rcpp::checkUserInterrupt();
  }
  setRxThreadId(-1);

  RObject parHist = vaeIterPrintGet_(printCtl >= 1);
  arma::mat zPopMatOut(N, zDim);
  if (isCovStep) zPopMatOut = zPopArg; else zPopMatOut.each_row() = zPop.t();

  List paramsOut = List::create(_["Wih"] = Wih, _["Whh"] = Whh, _["bih"] = bih,
                                _["bhh"] = bhh, _["fcW"] = fcW, _["fcB"] = fcB);
  IntegerVector mixnumOut(N);
  for (int i = 0; i < N; ++i) mixnumOut[i] = last.mixnum[i];
  return List::create(_["params"] = paramsOut, _["zPop"] = zPop, _["omega"] = omega,
                      _["a"] = a, _["intercept"] = intercept, _["beta"] = beta,
                      _["selected"] = selected, _["elboTrace"] = elboTrace,
                      _["parHist"] = parHist, _["mu"] = last.mu, _["zPopMat"] = zPopMatOut,
                      _["mixnum"] = mixnumOut);
}

//[[Rcpp::export]]
NumericVector boxCox_(NumericVector x = 1, double lambda=1, int yj = 0){
  NumericVector ret(x.size());
  for (unsigned int i = x.size(); i--;){
    ret[i] = _powerD(x[i], lambda, yj, 0.0, 1.0);
  }
  return ret;
}

//[[Rcpp::export]]
NumericVector iBoxCox_(NumericVector x = 1, double lambda=1, int yj = 0){
  NumericVector ret(x.size());
  for (unsigned int i = x.size(); i--;){
    ret[i] = _powerDi(x[i], lambda, yj, 0.0, 1.0);
  }
  return ret;
}

void saveIntoEnvrionment(Environment e) {
  unsigned int totN=op_focei.ntheta + op_focei.omegan;
  arma::Col<int> etaTrans(op_focei.etaTrans, op_focei.neta*3 + 3*(op_focei.ntheta + op_focei.omegan));
  e[".etaTrans"] = etaTrans;
  arma::vec fullTheta(op_focei.fullTheta, 4*(op_focei.ntheta+op_focei.omegan));
  e[".fullTheta"] = fullTheta;
  // no eta
  if (op_focei.neta == 0) {
    arma::vec gthetaGrad(op_focei.fullTheta, 4*(op_focei.ntheta+op_focei.omegan));
    e[".gthetaGrad"] = gthetaGrad;
  } else {
    size_t _nz = ((op_focei.neta+1)*(op_focei.neta+2)/2 + 6*(op_focei.neta+1) + 1) *
                  (size_t)getRxNsub(rx);
    {
      size_t _nall = (size_t)getRxNall(rx);
      size_t _nsub = (size_t)getRxNsub(rx);
      arma::vec etaUpper(op_focei.etaUpper,
                         (size_t)op_focei.gEtaGTransN*10 + op_focei.npars*(_nsub + 1) + _nz +
                         2*op_focei.neta * _nall + _nall + _nall*_nall +
                         op_focei.neta*5 + 2*op_focei.neta*op_focei.neta*_nsub + _nall);
      e[".etaUpper"] = etaUpper;
    }
  }
  arma::Col<int> gillRet(op_focei.gillRet,
                     2*totN+op_focei.npars+
                              op_focei.muRefN + op_focei.skipCovN);
  e[".gillRet"] = gillRet;
  arma::vec gillDf(op_focei.gillDf,7*totN + 2*op_focei.npars + getRxNsub(rx));
  e[".gillDf"] = gillDf;
}

void restoreFromEnvrionment(Environment e) {
  arma::Col<int> etaTrans = e[".etaTrans"];
  std::copy(etaTrans.begin(), etaTrans.end(), op_focei.etaTrans);
  arma::vec fullTheta = e[".fullTheta"];
  std::copy(fullTheta.begin(), fullTheta.end(), op_focei.fullTheta);
  // no eta
  if (op_focei.neta == 0) {
    arma::vec gthetaGrad = e[".gthetaGrad"];
    std::copy(gthetaGrad.begin(), gthetaGrad.end(), op_focei.fullTheta);
  } else {
    arma::vec etaUpper = e[".etaUpper"];
    std::copy(etaUpper.begin(), etaUpper.end(), op_focei.etaUpper);
  }
  arma::Col<int> gillRet = e[".gillRet"];
  std::copy(gillRet.begin(), gillRet.end(), op_focei.gillRet);
  arma::vec gillDf = e[".gillDf"];
  std::copy(gillDf.begin(), gillDf.end(), op_focei.gillDf);
}
