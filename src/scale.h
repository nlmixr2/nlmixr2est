#include <math.h>
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define expit(alpha, low, high) _powerDi(alpha, 1.0, 4, low, high)
#define probitInv(alpha, low, high) _powerDi(alpha, 1.0, 6, low, high)

#define scaleTypeNorm    1
#define scaleTypeNlmixr2 2
#define scaleTypeMult    3
#define scaleTypeMultAdd 4
#define scaleTypeNone    5

//  .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)

#define normTypeRescale2 1
#define normTypeRescale  2
#define normTypeMean     3
#define normTypeStd      4
#define normTypeLen      5
#define normTypeConstant 6

struct scaling {
  int npars; // number of parameters
  int scaleType; // scaling type
  int normType; // normalization type
  double scaleTo; // internal constant to say what we scale to
  double c1; // internal scaling constant
  double c2; // internal scaling constant
  double scaleCmin; // Cmin scaling constant
  double scaleCmax; // Cmax scaling constant
  // Iteration-print formatting, populated from the user's iterPrintControl()
  // sub-list via scaleApplyIterPrintControl().  Field names mirror the R
  // argument names of iterPrintControl() for code-search consistency.
  int useColor;
  int ncol;        // iterPrintControl$ncol         — columns per row
  int every;       // iterPrintControl$every        — print row every N evals (0=silent)
  int simple;      // iterPrintControl$simple       — single-row mode (skip U/X)
  int headerEvery; // iterPrintControl$headerEvery  — re-emit header every N prints (0=once)
  // showOfv: 1 = emit the "Function Val." column (the OFV column shown
  //              after the iteration counter).  Default.
  //          0 = skip that column entirely.  Used by estimators (like
  //              saem) that have no per-iteration objective function, so
  //              that NaN doesn't appear in the OFV slot.  The header,
  //              separator, and per-row prefix all shrink accordingly.
  int showOfv;
  // keyExtra: estimator-specific text appended to the "Key:" legend line
  // (between "X: Back-transformed parameters; " and the trailing newline).
  // NULL emits just the standard Key.  Used by focei to inject its
  // G/F/C/M gradient-method legend and omega note.  Setting it is the
  // estimator's responsibility; iterPrintControl() does not expose it.
  const char *keyExtra;
  // printCount: internal — count of parameter-print events emitted so far
  // (not iterations).  Gates periodic header re-prints when headerEvery > 0.
  int printCount;
  int save;
  int cn;
  double *initPar; // initial parameter estimates before scaling
  double *scaleC; // scaling C vector
  // xPar is an indicator for transformation type:
  // - xPar = -m-1 is logit transformation, ie logit(a, lower, upper); m
  //          matches the lower/upper bounds index
  // - xPar = 1 is a log transformation ie log(a)
  // - xPar = 0  exponential indicator.  Is the variable exp(x) (=1) or something else (=0)
  int *xPar;
  double *logitThetaLow;
  double *logitThetaHi;
  // Probit transform marker + bounds, parallel to logitThetaLow/Hi.
  // probitIdx[i] is 0 when parameter i has no probit transform, else
  // the 1-based index into probitThetaLow / probitThetaHi.  Stored
  // separately from xPar because positive xPar values 2-5 are already
  // reserved by scaleGetScaleC for omega scaling.  probitIdx may be
  // NULL when no parameter uses probit; scaleBackTransform treats that
  // as "no probit anywhere".
  int *probitIdx;
  double *probitThetaLow;
  double *probitThetaHi;
  // Backing storage for the six transform pointers above, populated
  // by scaleAttachXform() from the R-side xform list returned by
  // .iterPrintXParFromUi().  Estimators that own their own buffers
  // (focei merges per-theta xPar with per-omega xType, so the npars-
  // sized buffer lives on op_focei) leave these empty and assign the
  // raw pointers to their own arrays directly.
  std::vector<int>    xParStorage;
  std::vector<int>    probitIdxStorage;
  std::vector<double> logitLowStorage;
  std::vector<double> logitHiStorage;
  std::vector<double> probitLowStorage;
  std::vector<double> probitHiStorage;
  CharacterVector thetaNames;
  std::vector<int> niter;
  std::vector<int> iterType;
  std::vector<double> vPar;
  std::vector<double> vGrad;
  std::vector<int> niterGrad;
  std::vector<int> gradType;
#define iterTypeGill 1
#define iterTypeMixed 2
#define iterTypeForward 3
#define iterTypeCentral 4
#define iterTypeScaled 5
#define iterTypeUnscaled 6
#define iterTypeBack 7
#define iterTypeSens 8
};

static inline void scaleNone(scaling *scale) {
  scale->scaleType= scaleTypeNone;
  scale->normType=normTypeConstant;
  scale->scaleTo = 0.0;
  scale->every = 0;
  scale->save = 0;
}

static inline void scaleSetup(scaling *scale,
                              double *initPar,
                              double *scaleC,
                              int *xPar,
                              double *logitThetaLow,
                              double *logitThetaHi,
                              CharacterVector thetaNames,
                              int useColor,
                              int printNcol,
                              int print,
                              int normType,
                              int scaleType,
                              double scaleCmin,
                              double scaleCmax,
                              double scaleTo,
                              int npars) {
  scale->npars   = npars;
  scale->initPar = initPar;
  scale->scaleC = scaleC;
  scale->normType = normType;
  scale->scaleType = scaleType;
  scale->xPar = xPar;
  scale->logitThetaLow =logitThetaLow;
  scale->logitThetaHi =logitThetaHi;
  // Probit fields default to "no probit"; callers that wire it up
  // (saem/focei when those estimators gain probit support) set these
  // directly on the struct after scaleSetup returns.
  scale->probitIdx = NULL;
  scale->probitThetaLow = NULL;
  scale->probitThetaHi = NULL;
  scale->thetaNames = thetaNames;

  scale->useColor = useColor;
  scale->ncol = printNcol;
  scale->every = print;
  scale->simple = 0;
  scale->showOfv = 1;
  scale->keyExtra = NULL;
  scale->headerEvery = 10;
  scale->printCount = 0;
  scale->save = 1;

  scale->scaleCmin = scaleCmin;
  scale->scaleCmax = scaleCmax;
  scale->scaleTo = scaleTo;

  scale->vGrad.clear();
  scale->vPar.clear();
  scale->iterType.clear();
  scale->gradType.clear();
  scale->niter.clear();
  scale->niterGrad.clear();

  scale->cn = 0;

  double mn = scale->initPar[scale->npars-1],
    mx=scale->initPar[scale->npars-1],
    mean=0, oN=0, oM=0,s=0;
  double len=0;
  switch (scale->normType){
  case normTypeRescale2:
    // OptdesX
    // http://apmonitor.com/me575/uploads/Main/optimization_book.pdf
    for (unsigned int k = scale->npars-1; k--;){
      mn = min2(scale->initPar[k],mn);
      mx = max2(scale->initPar[k],mx);
      scale->scaleC[k] = NA_REAL;
    }
    if (fabs(mx-mn) < DBL_EPSILON) {
      for (unsigned int k = scale->npars-1; k--;){
        len += scale->initPar[k]*scale->initPar[k];
      }
      if (len < DBL_EPSILON){
        warning(_("all parameters are zero, cannot scale, run unscaled"));
        scale->c1 = 0;
        scale->c2 = 1;
        scale->normType = normTypeConstant;
      } else {
        warning(_("all parameters are the same value, switch to length normType"));
        scale->c1 = 0;
        scale->c2 = sqrt(len);
        scale->normType = 5;
      }
    } else {
      scale->c1 = (mx+mn)/2;
      scale->c2 = (mx-mn)/2;
    }
    break;
  case normTypeRescale: // Rescaling (min-max normalization)
    for (unsigned int k = scale->npars-1; k--;){
      mn = min2(scale->initPar[k],mn);
      mx = max2(scale->initPar[k],mx);
      scale->scaleC[k] = NA_REAL;
    }
    if (fabs(mx-mn) < DBL_EPSILON) {
      for (unsigned int k = scale->npars-1; k--;){
        len += scale->initPar[k]*scale->initPar[k];
      }
      if (len < DBL_EPSILON){
        warning(_("all parameters are zero, cannot scale, run unscaled"));
        scale->c1 = 0;
        scale->c2 = 1;
        scale->normType = normTypeConstant;
      } else {
        warning(_("all parameters are the same value, switch to length normType"));
        scale->c1 = 0;
        scale->c2 = sqrt(len);
        scale->normType = 5;
      }
    } else {
      scale->c1 = mn;
      scale->c2 = (mx-mn);
    }
    break;
  case normTypeMean: // Mean normalization
    for (unsigned int k = scale->npars-1; k--;){
      mn = min2(scale->initPar[k],mn);
      mx = max2(scale->initPar[k],mx);
      oN++;
      mean += (scale->initPar[k]-mean)/oN;
      scale->scaleC[k] = NA_REAL;
    }
    if (fabs(mx-mn) < DBL_EPSILON) {
      for (unsigned int k = scale->npars-1; k--;){
        len += scale->initPar[k]*scale->initPar[k];
      }
      if (len < DBL_EPSILON){
        warning(_("all parameters are zero, cannot scale, run unscaled"));
        scale->c1 = 0;
        scale->c2 = 1;
        scale->normType = normTypeConstant;
      } else {
        warning(_("all parameters are the same value, switch to length normType"));
        scale->c1 = 0;
        scale->c2 = sqrt(len);
        scale->normType = 5;
      }
    } else {
      scale->c1 = mean;
      scale->c2 = (mx-mn);
    }
    break;
  case normTypeStd: // Standardization
    for (unsigned int k = scale->npars-1; k--;){
      mn = min2(scale->initPar[k],mn);
      mx = max2(scale->initPar[k],mx);
      oM= mean;
      oN++;
      mean += (scale->initPar[k]-mean)/oN;
      s += (scale->initPar[k]-mean)*(scale->initPar[k]-oM);
      scale->scaleC[k] = NA_REAL;
    }
    if (fabs(mx-mn) < DBL_EPSILON) {
      for (unsigned int k = scale->npars-1; k--;){
        len += scale->initPar[k]*scale->initPar[k];
      }
      if (len < DBL_EPSILON){
        warning(_("all parameters are zero, cannot scale, run unscaled"));
        scale->c1 = 0;
        scale->c2 = 1;
        scale->normType = normTypeConstant;
      } else {
        warning(_("all parameters are the same value, switch to length normType"));
        scale->c1 = 0;
        scale->c2 = sqrt(len);
        scale->normType = 5;
      }
    } else {
      scale->c1 = mean;
      scale->c2 = sqrt(s/(oN-1));
    }
    break;
  case normTypeLen: // Normalize to length.
    for (unsigned int k = scale->npars-1; k--;){
      len += scale->initPar[k]*scale->initPar[k];
      scale->scaleC[k] = NA_REAL;
    }
    if (len < DBL_EPSILON){
      warning(_("all parameters are zero, cannot scale, run unscaled"));
      scale->c1 = 0;
      scale->c2 = 1;
      scale->normType = normTypeConstant;
    } else {
      scale->c1 = 0;
      scale->c2 = sqrt(len);
    }
    break;
  case normTypeConstant:
    // No Normalization
    for (unsigned int k = scale->npars-1; k--;){
      scale->scaleC[k] = NA_REAL;
    }
    scale->c1 = 0;
    scale->c2 = 1;
    break;
  default:
    stop("unrecognized normalization (normType=%d)", scale->normType);
  }
}


static inline double scaleGetScaleC(scaling *scale, int i){
  if (ISNA(scale->scaleC[i]) || isnan(scale->scaleC[i])) {
    switch (scale->xPar[i]){
    case 1: // log
      scale->scaleC[i]=1.0;
      break;
    case 2: // diag^2
      scale->scaleC[i]=1.0/max2(fabs(scale->initPar[i]), scale->scaleCmin);
      break;
    case 3: // exp(diag)
      scale->scaleC[i] = 1.0/2.0;
      break;
    case 4: // Identity diagonal chol(Omega ^-1)
    case 5: // off diagonal chol(Omega^-1)
      scale->scaleC[i] = 1.0/max2(2.0*fabs(scale->initPar[i]), scale->scaleCmin);
      break;
    default:
      scale->scaleC[i]= 1.0/max2(fabs(scale->initPar[i]), scale->scaleCmin);
      break;
    }
  }
  return min2(max2(scale->scaleC[i], scale->scaleCmin), scale->scaleCmax);
}

static inline double scaleAdjustGradScale(scaling *scale, double grad, double *x, int i) {
  // Here we have f(unscalePar(x)) hence the derivative is df/du*du/dx
  // (by chain rule); df/du is provided by the routine, du/dx is provided here
  double scaleTo = scale->scaleTo, C=scaleGetScaleC(scale, i);
  switch(scale->scaleType) {
  case scaleTypeNorm: // normalized
    // here du/dx = scale->c2
    return grad*scale->c2;
    break;
  case scaleTypeNlmixr2: // log vs linear scales and/or ranges
    // here du/dx = C*x[i]
    return grad*C;
    break;
  case scaleTypeMult: // simple multiplicative scaling
    if (scale->scaleTo != 0){
      return grad*scale->initPar[i]/scaleTo;
    } else {
      return grad;
    }
    break;
  case scaleTypeMultAdd: // log non-log multiplicative scaling
    if (scale->scaleTo > 0){
      switch (scale->xPar[i]){
      case 1:
        return grad;
      default:
        return grad*scale->initPar[i]/scaleTo;
      }
    } else {
      return grad;
    }
  case scaleTypeNone:
    // no scaling
    return grad;
  default:
    return grad;
  }
  return 0;
}

static inline double scaleUnscalePar(scaling *scale, double *x, int i){
  double scaleTo = scale->scaleTo, C=scaleGetScaleC(scale, i);
  switch(scale->scaleType) {
  case scaleTypeNorm: // normalized
    return x[i]*scale->c2+scale->c1;
    break;
  case scaleTypeNlmixr2: // log vs linear scales and/or ranges
    if (scale->normType != normTypeConstant){
      scaleTo = (scale->initPar[i] - scale->c1)/scale->c2;
    } else if (scaleTo == 0){
      scaleTo=scale->initPar[i];
    }
    return (x[i]-scaleTo)*C + scale->initPar[i];
    break;
  case scaleTypeMult: // simple multiplicative scaling
    if (scale->scaleTo != 0){
      return x[i]*scale->initPar[i]/scaleTo;
    } else {
      return x[i];
    }
    break;
  case scaleTypeMultAdd: // log non-log multiplicative scaling
    if (scale->scaleTo > 0){
      switch (scale->xPar[i]){
      case 1:
        return (x[i]-scaleTo) + scale->initPar[i];
      default:
        return x[i]*scale->initPar[i]/scaleTo;
      }
    } else {
      return x[i];
    }
  case scaleTypeNone: // no scaling
    return x[i];
  default:
    if (scale->scaleTo > 0){
      return (x[i]-scaleTo)*1 + scale->initPar[i];
    } else {
      return x[i];
    }
  }
  return 0;
}

static inline double scaleScalePar(scaling *scale, double *x, int i){
  double scaleTo = scale->scaleTo, C=scaleGetScaleC(scale, i);
  switch(scale->scaleType) {
  case scaleTypeNorm:
    return (x[i]-scale->c1)/scale->c2;
  case scaleTypeNlmixr2:
    if (scale->normType <= 5){
      scaleTo = (scale->initPar[i]-scale->c1)/scale->c2;
    } else if (scaleTo == 0){
      scaleTo=scale->initPar[i];
    }
    return (x[i]-scale->initPar[i])/C + scaleTo;
    break;
  case scaleTypeMult: // simple multiplicative scaling
    if (scale->scaleTo > 0){
      return x[i]/scale->initPar[i]*scale->scaleTo;
    } else {
      return x[i];
    }
    break;
  case scaleTypeMultAdd: // log non-log multiplicative scaling
    if (scale->scaleTo > 0){
      switch (scale->xPar[i]){
      case 1:
        return (x[i]-scale->initPar[i]) + scale->scaleTo;
      default:
        return x[i]/scale->initPar[i] * scale->scaleTo;
      }
    } else {
      return x[i];
    }
  case scaleTypeNone: // no scaling
    return x[i];
  default:
    if (scale->scaleTo > 0){
      return (x[i]-scale->initPar[i]) + scale->scaleTo;
    } else {
      return x[i];
    }
  }
  return 0;
}

// Read the user-facing iterPrintControl() sub-list and populate the
// matching iteration-print fields on a scaling struct (every, ncol,
// headerEvery, useColor, simple).  Each estimator's setup code calls
// this exactly once to wire its R-side configuration into the shared
// printer.  The keyExtra pointer is estimator-internal and remains
// the caller's responsibility to set.
static inline void scaleApplyIterPrintControl(scaling *scale,
                                              const Rcpp::List &ipc) {
  scale->every       = Rcpp::as<int>(ipc["every"]);
  scale->ncol        = Rcpp::as<int>(ipc["ncol"]);
  scale->headerEvery = Rcpp::as<int>(ipc["headerEvery"]);
  scale->useColor    = Rcpp::as<int>(ipc["useColor"]);
  scale->simple      = Rcpp::as<int>(ipc["simple"]);
}

// Wire every per-parameter transform array onto a scaling struct from
// a single R-side xform sub-list (output of .iterPrintXParFromUi()).
// Required list elements:
//
//   xPar            integer (length scale->npars).  1=log, -m=m-th
//                    logit, 0=no log/logit transform.
//   probitIdx       integer (same length).  k=k-th probit transform,
//                    0=no probit transform.
//   logitThetaLow   numeric.  Bounds indexed by -xPar[i]-1 when xPar<0.
//   logitThetaHi    numeric.  Same length as logitThetaLow.
//   probitThetaLow  numeric.  Bounds indexed by probitIdx[i]-1 when >0.
//   probitThetaHi   numeric.  Same length as probitThetaLow.
//
// The struct's own std::vector backing storage absorbs the data so
// the pointers survive the caller's stack frame.  Empty arrays leave
// the corresponding pointer NULL (scaleBackTransform treats that as
// "no transform of that type anywhere").
//
// This is the single point of wiring for log/logit/probit back-
// transforms — every estimator that calls it inherits identical
// behavior across the iteration X row, parHistData, and the final-
// fit-summary parFixed column.  Estimators with method-specific
// per-parameter codes (focei mixes per-theta xPar with per-omega
// xType in the same npars-sized array) skip this helper and assign
// the raw pointers directly.
static inline void scaleAttachXform(scaling *scale,
                                    const Rcpp::List &xform) {
  IntegerVector x  = Rcpp::as<IntegerVector>(xform["xPar"]);
  IntegerVector pi = Rcpp::as<IntegerVector>(xform["probitIdx"]);
  NumericVector ll = Rcpp::as<NumericVector>(xform["logitThetaLow"]);
  NumericVector lh = Rcpp::as<NumericVector>(xform["logitThetaHi"]);
  NumericVector pl = Rcpp::as<NumericVector>(xform["probitThetaLow"]);
  NumericVector ph = Rcpp::as<NumericVector>(xform["probitThetaHi"]);
  scale->xParStorage     .assign(x.begin(),  x.end());
  scale->probitIdxStorage.assign(pi.begin(), pi.end());
  scale->logitLowStorage .assign(ll.begin(), ll.end());
  scale->logitHiStorage  .assign(lh.begin(), lh.end());
  scale->probitLowStorage.assign(pl.begin(), pl.end());
  scale->probitHiStorage .assign(ph.begin(), ph.end());
  scale->xPar           = scale->xParStorage     .empty() ? NULL : scale->xParStorage     .data();
  scale->probitIdx      = scale->probitIdxStorage.empty() ? NULL : scale->probitIdxStorage.data();
  scale->logitThetaLow  = scale->logitLowStorage .empty() ? NULL : scale->logitLowStorage .data();
  scale->logitThetaHi   = scale->logitHiStorage  .empty() ? NULL : scale->logitHiStorage  .data();
  scale->probitThetaLow = scale->probitLowStorage.empty() ? NULL : scale->probitLowStorage.data();
  scale->probitThetaHi  = scale->probitHiStorage .empty() ? NULL : scale->probitHiStorage .data();
}

// Wrap-continuation marker emitted at the start of a wrapped row when a
// row has more parameter columns than fit in the chosen `ncol` width.
// Width matches the left-hand "label" prefix on the data rows so the
// wrapped parameter columns line up under the leading params:
//
//   showOfv = 1 → "|.....................|"  (23 chars; covers
//                  "|    #| Function Val. |" or "|    U|               |")
//   showOfv = 0 → "|.....|"                  (7 chars; covers "|    #|")
//
// `colored` chooses the ANSI-bracketed variant used when the row is
// emitted with `useColor` and the wrap reaches the final visible column.
static inline const char *scaleWrapMarker(scaling *scale, int colored) {
  if (colored) {
    return scale->showOfv ? "\n\033[4m|.....................|"
                          : "\n\033[4m|.....|";
  }
  return scale->showOfv ? "\n|.....................|" : "\n|.....|";
}

// Emit the separator line under the column header / between iteration
// blocks.  When `showOfv` is 0 (saem and similar — no per-iteration
// objective function) the Function-Val segment is skipped so the
// separator's column count matches the header above and the rows below.
static inline void scalePrintLine(scaling *scale, int ncol){
  if (scale->showOfv) {
    RSprintf("|-----+---------------+");
  } else {
    RSprintf("|-----+");
  }
  for (int i = 0; i < ncol; i++){
    if (i == ncol-1)
      RSprintf("-----------|");
    else
      RSprintf("-----------+");
  }
  RSprintf("\n");
}

// withKey: 1 = emit the "Key:" legend block (U/X transform legend plus any
//              estimator-specific keyExtra) above the column-label header.
//          0 = emit only the column-label header and its separator line.
// The legend explains the row codes and only needs to be shown once, so the
// startup header (printed by each estimator at fit start) passes withKey=1
// while the periodic re-emits every `headerEvery` prints pass withKey=0 —
// repeating just the compact column labels for readability without
// reprinting the multi-line explanation on every header refresh.
static inline void scalePrintHeader(scaling *scale, int withKey = 1) {
  if (scale->every != 0) {
    if (withKey) {
      // Match scalePrintFun's auto-skip rules so the Key text only
      // mentions rows that will actually appear in iteration output.
      int skipU = (scale->scaleType == scaleTypeNone);
      int anyXform = 0;
      for (int k = 0; k < scale->npars; k++) {
        if (scale->xPar[k] != 0) { anyXform = 1; break; }
        if (scale->probitIdx != NULL && scale->probitIdx[k] != 0) { anyXform = 1; break; }
      }
      int skipX = skipU && !anyXform;
      if (!scale->simple && (!skipU || !skipX || scale->keyExtra != NULL)) {
        if (scale->useColor)
          RSprintf("\033[1mKey:\033[0m ");
        else
          RSprintf("Key: ");
        if (!skipU) RSprintf("U: Unscaled Parameters; ");
        if (!skipX) RSprintf("X: Back-transformed parameters; ");
        // Estimator-specific Key suffix (e.g. focei's G/F/C/M gradient
        // legend and omega note).  When NULL the standard "Key:" line is
        // closed with a newline so the column header follows on a fresh row.
        if (scale->keyExtra != NULL) {
          RSprintf("%s", scale->keyExtra);
        } else {
          RSprintf("\n");
        }
      }
    }
    int i, finalize=0, n=scale->thetaNames.size();
    if (scale->showOfv) {
      RSprintf("\n|    #| Function Val. |");
    } else {
      RSprintf("\n|    #|");
    }
    std::string tmpS;
    for (i = 0; i < n; i++){
      tmpS = scale->thetaNames[i];
      RSprintf("%#10s |", tmpS.c_str());
      if ((i + 1) != n && (i + 1) % scale->ncol == 0){
        if (scale->useColor && scale->ncol + i  >= n){
          RSprintf("%s", scaleWrapMarker(scale, 1));
        } else {
          RSprintf("%s", scaleWrapMarker(scale, 0));
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->ncol == 0){
          if (scale->useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
    scalePrintLine(scale, min2(scale->npars, scale->ncol));
  }
}

// Apply the back-transform encoded by xParCode/probitCode to a single
// (unscaled, model-scale) estimate.  Used by both the iteration-print
// X row (scalePrintFun, below) and the final-fit-summary parFixed
// loop in src/inner.cpp so the two paths stay in sync.
//
//   xParCode == 1              -> log:    exp(est)
//   xParCode <= -1             -> logit:  expit(est, logitLow[m], logitHi[m])
//                                          where m = -xParCode-1
//   probitCode >= 1            -> probit: probitInv(est, probitLow[m], probitHi[m])
//                                          where m = probitCode-1
//   anything else              -> identity (no transform)
//
// Bounds pointers may be NULL when the corresponding code is never
// triggered (xParCode never < 0, or probitCode always 0).  The caller
// is responsible for ensuring the bounds arrays are large enough for
// any encoded index that does appear.
static inline double scaleBackTransform(double est, int xParCode, int probitCode,
                                        const double *logitLow, const double *logitHi,
                                        const double *probitLow, const double *probitHi) {
  if (xParCode == 1) {
    return exp(est);
  }
  if (xParCode <= -1) {
    int m = -xParCode - 1;
    return expit(est, logitLow[m], logitHi[m]);
  }
  if (probitCode >= 1) {
    int m = probitCode - 1;
    return probitInv(est, probitLow[m], probitHi[m]);
  }
  return est;
}

static inline void scalePrintFun(scaling *scale, double *x, double f) {
  // Scaled
  int finalize = 0, i = 0;
  scale->cn = scale->cn+1;
  // Auto-skip degenerate rows.  When scaleType == None, scaleUnscalePar
  // returns x[i] unchanged so the U row is identical to # and is skipped.
  // When U is skipped AND no xPar entry asks for a back-transform, the
  // X row is also identical to # and is skipped too — methods with no
  // scaling and no log/logit-transformed parameters collapse to a single
  // # row.  The user's iterPrintControl(simple=TRUE) override still
  // forces both rows off regardless of auto-detection.
  int skipU = (scale->scaleType == scaleTypeNone);
  int anyXform = 0;
  for (i = 0; i < scale->npars; i++){
    if (scale->xPar[i] != 0) { anyXform = 1; break; }
    if (scale->probitIdx != NULL && scale->probitIdx[i] != 0) { anyXform = 1; break; }
  }
  int skipX = skipU && !anyXform;
  if (scale->save) {
    scale->niter.push_back(scale->cn);

    scale->vPar.push_back(f);
    scale->iterType.push_back(iterTypeScaled);
    for (i = 0; i < scale->npars; i++){
      scale->vPar.push_back(x[i]);
    }
    if (!scale->simple && !skipU) {
      // Unscaled
      scale->iterType.push_back(iterTypeUnscaled);
      scale->niter.push_back(scale->niter.back());
      scale->vPar.push_back(f);
      for (i = 0; i < scale->npars; i++){
        scale->vPar.push_back(scaleUnscalePar(scale, x, i));
      }
    }
    if (!scale->simple && !skipX) {
      // Back-transformed (7)
      scale->iterType.push_back(iterTypeBack);
      scale->niter.push_back(scale->niter.back());
      scale->vPar.push_back(f);
      for (i = 0; i < scale->npars; i++){
        int probitCode = (scale->probitIdx != NULL) ? scale->probitIdx[i] : 0;
        scale->vPar.push_back(scaleBackTransform(scaleUnscalePar(scale, x, i),
                                                 scale->xPar[i], probitCode,
                                                 scale->logitThetaLow, scale->logitThetaHi,
                                                 scale->probitThetaLow, scale->probitThetaHi));
      }
    }
  }
  if (scale->every != 0 &&
      scale->cn % scale->every == 0){
    // Count this parameter-print event and, when configured, re-emit the
    // column header every `headerEvery` events (event 1 already has the
    // startup header printed elsewhere, so re-prints happen at 1+N, 1+2N, ...).
    scale->printCount++;
    if (scale->headerEvery > 0 &&
        scale->printCount > 1 &&
        ((scale->printCount - 1) % scale->headerEvery == 0)) {
      // Periodic refresh: repeat only the column labels, not the legend.
      scalePrintHeader(scale, 0);
    }
    if (scale->showOfv) {
      if (scale->useColor && !isRstudio())
        RSprintf("|\033[1m%5d\033[0m|%#14.8g |", scale->cn, f);
      else
        RSprintf("|%5d|%#14.8g |", scale->cn, f);
    } else {
      if (scale->useColor && !isRstudio())
        RSprintf("|\033[1m%5d\033[0m|", scale->cn);
      else
        RSprintf("|%5d|", scale->cn);
    }
    for (i = 0; i < scale->npars; i++){
      RSprintf("%#10.4g |", x[i]);
      if ((i + 1) != scale->npars && (i + 1) % scale->ncol == 0){
        if (scale->useColor && scale->ncol + i  > scale->npars){
          RSprintf("%s", scaleWrapMarker(scale, 1));
        } else {
          RSprintf("%s", scaleWrapMarker(scale, 0));
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->ncol == 0){
          if (scale->useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
    if (!scale->simple && !skipU) {
      if (scale->showOfv) RSprintf("|    U|               |");
      else                RSprintf("|    U|");
      for (i = 0; i < scale->npars; i++){
        RSprintf("%#10.4g |", scaleUnscalePar(scale, x, i));
        if ((i + 1) != scale->npars && (i + 1) % scale->ncol == 0){
          if (scale->useColor && scale->ncol + i  > scale->npars){
            RSprintf("%s", scaleWrapMarker(scale, 1));
          } else {
            RSprintf("%s", scaleWrapMarker(scale, 0));
          }
        }
      }
      if (finalize){
        while(true){
          if ((i++) % scale->ncol == 0){
            if (scale->useColor) RSprintf("\033[0m");
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
    if (!scale->simple && !skipX) {
      if (scale->showOfv) RSprintf("|    X|               |");
      else                RSprintf("|    X|");
      for (i = 0; i < scale->npars; i++){
        int probitCode = (scale->probitIdx != NULL) ? scale->probitIdx[i] : 0;
        RSprintf("%#10.4g |",
                 scaleBackTransform(scaleUnscalePar(scale, x, i),
                                    scale->xPar[i], probitCode,
                                    scale->logitThetaLow, scale->logitThetaHi,
                                    scale->probitThetaLow, scale->probitThetaHi));
        if ((i + 1) != scale->npars && (i + 1) % scale->ncol == 0){
          if (scale->useColor && scale->ncol + i >= scale->npars){
            RSprintf("%s", scaleWrapMarker(scale, 1));
          } else {
            RSprintf("%s", scaleWrapMarker(scale, 0));
          }
        }
      }
      if (finalize){
        while(true){
          if ((i++) % scale->ncol == 0){
            if (scale->useColor) RSprintf("\033[0m");
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
  // Universal user-interrupt check at each per-iteration print event.  Doing
  // this here (rather than in each estimator's outer loop) ensures every
  // method routed through scalePrintFun — saem, nlm, optim, nls, nlminb — is
  // interruptible without each one needing its own Rcpp::checkUserInterrupt()
  // call site.
  Rcpp::checkUserInterrupt();
}

static inline void scalePrintGrad(scaling *scale, double *gr, int type) {
  int finalize = 0, i = 0;
  // if (op_focei.derivMethod == 0){
  //   if (op_focei.curGill == 1){
  //     gradType.push_back(1);
  //   } else if (op_focei.curGill == 2){
  //     gradType.push_back(5);
  //   } else if (op_focei.mixDeriv){
  //     gradType.push_back(2);
  //   } else{
  //     gradType.push_back(3);
  //   }
  // } else {
  //   gradType.push_back(4);
  // }
  if (scale->save) {
    scale->niterGrad.push_back(scale->niter.back());
    scale->gradType.push_back(type);
  }
  if (scale->every != 0 &&
      scale->cn % scale->every == 0){
    // Method-specific label for the gradient row, keyed by `type`:
    //   1=Gill, 2=Mixed, 3=Forward, 4=Central, 5=Shi21.
    // Any other code (e.g. iterTypeSens=8 from nlm/optim) falls through to a
    // generic "Gradient" label.
    const char *label = NULL;
    switch (type) {
    case 1:  label = "    G|    Gill Diff. |"; break;  // Gill
    case 2:  label = "    M|   Mixed Diff. |"; break;  // Mixed
    case 3:  label = "    F| Forward Diff. |"; break;  // Forward
    case 4:  label = "    C| Central Diff. |"; break;  // Central
    case 5:  label = "    S|   Shi21 Diff. |"; break;  // Shi21
    default: label = "    G|    Gradient   |"; break;
    }
    if (scale->useColor && scale->ncol >= scale->npars){
      RSprintf("|\033[4m%s", label);
    } else {
      RSprintf("|%s", label);
    }
    for (i = 0; i < scale->npars; i++){
      RSprintf("%#10.4g ", gr[i]);
      if (scale->useColor && scale->ncol >= scale->npars && i == scale->npars-1){
        RSprintf("\033[0m");
      }
      RSprintf("|");
      if ((i + 1) != scale->npars && (i + 1) % scale->ncol == 0){
        if (scale->useColor && scale->ncol + i  >= scale->npars){
          RSprintf("%s", scaleWrapMarker(scale, 1));
        } else {
          RSprintf("%s", scaleWrapMarker(scale, 0));
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->ncol == 0){
          if (scale->useColor) RSprintf("\033[0m");
          RSprintf("\n");
          break;
        } else {
          RSprintf("...........|");
        }
      }
    } else {
      RSprintf("\n");
    }
    if (!scale->useColor){
      scalePrintLine(scale, min2(scale->npars, scale->ncol));
    }
  }
  if (scale->save) {
    scale->vGrad.push_back(NA_REAL); // Gradient doesn't record objf
    for (i = 0; i < scale->npars; i++){
      scale->vGrad.push_back(gr[i]);
    }
  }
}

static inline RObject scaleParHisDf(scaling *scale) {
  if (scale->iterType.size() == 0)  return R_NilValue;
  CharacterVector dfNames(3+scale->thetaNames.size());
  dfNames[0] = "iter";
  dfNames[1] = "type";
  dfNames[2] = "objf";
  int i;
  for (i = 0; i < scale->thetaNames.size(); i++){
    dfNames[i+3] = scale->thetaNames[i];
  }
  List ret(3+scale->thetaNames.size());
  int sz = scale->niter.size()+scale->niterGrad.size();
  std::vector<int> iter;
  iter.reserve(sz);
  iter.insert(iter.end(), scale->niter.begin(), scale->niter.end());
  iter.insert(iter.end(), scale->niterGrad.begin(), scale->niterGrad.end());
  ret[0] = iter;
  IntegerVector tmp;
  tmp = IntegerVector(sz);
  std::vector<int> typ;
  typ.reserve(sz);
  typ.insert(typ.end(), scale->iterType.begin(), scale->iterType.end());
  typ.insert(typ.end(), scale->gradType.begin(), scale->gradType.end());
  tmp = typ;
  tmp.attr("levels") = CharacterVector::create("Gill83 Gradient", "Mixed Gradient",
                                               "Forward Difference", "Central Difference",
                                               "Scaled", "Unscaled", "Back-Transformed",
                                               "Forward Sensitivity");
  tmp.attr("class") = "factor";
  ret[1] = tmp;
  arma::mat cPar(scale->vPar.size()/scale->iterType.size(), scale->iterType.size());
  std::copy(scale->vPar.begin(), scale->vPar.end(), cPar.begin());
  arma::mat vals;
  if (scale->vGrad.size() > 0){
    arma::mat cGrad(scale->vGrad.size()/scale->gradType.size(), scale->gradType.size());
    std::copy(scale->vGrad.begin(), scale->vGrad.end(), cGrad.begin());
    cPar = cPar.t();
    cGrad = cGrad.t();
    vals = arma::join_cols(cPar, cGrad);
  } else {
    cPar = cPar.t();
    vals = cPar;
  }
  for (i = 0; i < scale->thetaNames.size()+1; i++){
    ret[i+2]= vals.col(i);
  }
  scale->vGrad.clear();
  scale->vPar.clear();
  scale->iterType.clear();
  scale->gradType.clear();
  scale->niter.clear();
  scale->niterGrad.clear();
  ret.attr("names")=dfNames;
  ret.attr("class") = "data.frame";
  ret.attr("row.names")=IntegerVector::create(NA_INTEGER, -sz);
  return ret;
}
