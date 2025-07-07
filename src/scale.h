#include <math.h>
#define min2( a , b )  ( (a) < (b) ? (a) : (b) )
#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define expit(alpha, low, high) _powerDi(alpha, 1.0, 4, low, high)

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
  int useColor;
  int printNcol;
  int print;
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
  scale->print = 0;
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
  scale->thetaNames = thetaNames;

  scale->useColor = useColor;
  scale->printNcol = printNcol;
  scale->print = print;
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
  default:
    if (scale->scaleTo > 0){
      return (x[i]-scale->initPar[i]) + scale->scaleTo;
    } else {
      return x[i];
    }
  }
  return 0;
}

static inline void scalePrintLine(int ncol){
  RSprintf("|-----+---------------+");
  for (int i = 0; i < ncol; i++){
    if (i == ncol-1)
      RSprintf("-----------|");
    else
      RSprintf("-----------+");
  }
  RSprintf("\n");
}

void scalePrintHeader(scaling *scale) {
  if (scale->print != 0) {
    if (scale->useColor)
      RSprintf("\033[1mKey:\033[0m ");
    else
      RSprintf("Key: ");
    RSprintf("U: Unscaled Parameters; ");
    RSprintf("X: Back-transformed parameters; \n");
    int i, finalize=0, n=scale->thetaNames.size();
    RSprintf("\n|    #| Function Val. |");
    std::string tmpS;
    for (i = 0; i < n; i++){
      tmpS = scale->thetaNames[i];
      RSprintf("%#10s |", tmpS.c_str());
      if ((i + 1) != n && (i + 1) % scale->printNcol == 0){
        if (scale->useColor && scale->printNcol + i  >= n){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->printNcol == 0){
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
    scalePrintLine(min2(scale->npars, scale->printNcol));
  }
}

void scalePrintFun(scaling *scale, double *x, double f) {
  // Scaled
  int finalize = 0, i = 0;
  scale->cn = scale->cn+1;
  if (scale->save) {
    scale->niter.push_back(scale->cn);

    scale->vPar.push_back(f);
    scale->iterType.push_back(iterTypeScaled);
    for (i = 0; i < scale->npars; i++){
      scale->vPar.push_back(x[i]);
    }
    // Unscaled
    scale->iterType.push_back(iterTypeUnscaled);
    scale->niter.push_back(scale->niter.back());
    scale->vPar.push_back(f);
    for (i = 0; i < scale->npars; i++){
      scale->vPar.push_back(scaleUnscalePar(scale, x, i));
    }
    // Back-transformed (7)
    scale->iterType.push_back(iterTypeBack);
    scale->niter.push_back(scale->niter.back());
    scale->vPar.push_back(f);
    for (i = 0; i < scale->npars; i++){
      if (scale->xPar[i] == 1){
        scale->vPar.push_back(exp(scaleUnscalePar(scale, x, i)));
      } else if (scale->xPar[i] < 0){
        int m = -scale->xPar[i]-1;
        scale->vPar.push_back(expit(scaleUnscalePar(scale, x, i), scale->logitThetaLow[m], scale->logitThetaHi[m]));
      } else {
        scale->vPar.push_back(scaleUnscalePar(scale, x, i));
      }
    }
  }
  if (scale->print != 0 &&
      scale->cn % scale->print == 0){
    if (scale->useColor && !isRstudio())
      RSprintf("|\033[1m%5d\033[0m|%#14.8g |", scale->cn, f);
    else
      RSprintf("|%5d|%#14.8g |", scale->cn, f);
    for (i = 0; i < scale->npars; i++){
      RSprintf("%#10.4g |", x[i]);
      if ((i + 1) != scale->npars && (i + 1) % scale->printNcol == 0){
        if (scale->useColor && scale->printNcol + i  > scale->npars){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->printNcol == 0){
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
    RSprintf("|    U|               |");
    for (i = 0; i < scale->npars; i++){
      RSprintf("%#10.4g |", scaleUnscalePar(scale, x, i));
      if ((i + 1) != scale->npars && (i + 1) % scale->printNcol == 0){
        if (scale->useColor && scale->printNcol + i  > scale->npars){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->printNcol == 0){
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
    RSprintf("|    X|               |");
    for (i = 0; i < scale->npars; i++){
      if (scale->xPar[i] == 1){
        RSprintf("%#10.4g |", exp(scaleUnscalePar(scale, x, i)));
      } else if (scale->xPar[i] < 0){
        int m = -scale->xPar[i]-1;
        RSprintf("%#10.4g |", expit(scaleUnscalePar(scale, x, i), scale->logitThetaLow[m], scale->logitThetaHi[m]));
      } else {
        RSprintf("%#10.4g |", scaleUnscalePar(scale, x, i));
      }
      if ((i + 1) != scale->npars && (i + 1) % scale->printNcol == 0){
        if (scale->useColor && scale->printNcol + i >= scale->npars){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->printNcol == 0){
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

void scalePrintGrad(scaling *scale, double *gr, int type) {
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
  if (scale->print != 0 &&
      scale->cn % scale->print == 0){
    if (scale->useColor && scale->printNcol >= scale->npars){
      RSprintf("|\033[4m    G|   Gradient    |");
    } else {
      RSprintf("|    G|    Gradient   |");
    }
    for (i = 0; i < scale->npars; i++){
      RSprintf("%#10.4g ", gr[i]);
      if (scale->useColor && scale->printNcol >= scale->npars && i == scale->npars-1){
        RSprintf("\033[0m");
      }
      RSprintf("|");
      if ((i + 1) != scale->npars && (i + 1) % scale->printNcol == 0){
        if (scale->useColor && scale->printNcol + i  >= scale->npars){
          RSprintf("\n\033[4m|.....................|");
        } else {
          RSprintf("\n|.....................|");
        }
        finalize=1;
      }
    }
    if (finalize){
      while(true){
        if ((i++) % scale->printNcol == 0){
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
      scalePrintLine(min2(scale->npars, scale->printNcol));
    }
  }
  if (scale->save) {
    scale->vGrad.push_back(NA_REAL); // Gradient doesn't record objf
    for (i = 0; i < scale->npars; i++){
      scale->vGrad.push_back(gr[i]);
    }
  }
}

RObject scaleParHisDf(scaling *scale) {
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
