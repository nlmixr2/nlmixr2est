#define STRICT_R_HEADER
#include "censResid.h"

bool censTruncatedMvnReturnInterestingLimits(arma::vec& dv, arma::vec& dvt,
                                             arma::vec& ipred, arma::vec &ipredt,
                                             arma::vec& pred, arma::vec &predt,
                                             arma::Col<int> &cens, arma::vec &limit,
                                             arma::vec& lambda, arma::vec &yj, arma::vec& low, arma::vec& hi, 
                                             arma::vec &lowerLim, arma::vec &upperLim, arma::vec &ri,
                                             bool &doSim, int& censMethod,
                                             arma::uvec &normRelated, arma::uvec &normIdx,
                                             arma::uvec &nonNormIdx) {
  bool interestingLim = false;
  for (int i = dv.size(); i--;) {
    // ipredt is the transformation information
    int inYj = (int)yj[i];
    int cyj, dist;
    _splitYj(&inYj, &dist,  &cyj);
    if (dist == rxDistributionNorm ||
        dist == rxDistributionDnorm) {
      normRelated[i] = 1;
      ipred[i]= _powerDi(ipredt[i], lambda[i], (int)yj[i], low[i], hi[i]);
      pred[i]= _powerDi(predt[i], lambda[i], (int)yj[i], low[i], hi[i]);
      switch (cens[i]){
      case 1:
        // (limit, dv) ; limit could be -inf
        if (R_FINITE(limit[i])){
          double lim0TBS = _powerD(limit[i], lambda[i], (int)yj[i], low[i], hi[i]);
          double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
          double sd = sqrt(ri[i]);
          interestingLim=true;
          lowerLim[i] = limit[i];
          upperLim[i] = dv[i];
          if (doSim && (censMethod == CENS_TNORM || censMethod == CENS_IPRED || censMethod == CENS_PRED)) {
            if (censMethod == CENS_TNORM) {
              dvt[i] = truncnorm(ipredt[i], sd, lim0TBS, lim1TBS);
            } else if (censMethod == CENS_PRED) {
              dvt[i] = predt[i];
            } else {
              dvt[i] = ipredt[i];
            }
            if (!R_finite(ipredt[i])) {
              normRelated[i] = 0;
            }
            if (R_finite(dvt[i])) {
              dv[i] =_powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
            } else {
              normRelated[i] = 0;
            }
          } else if (R_finite(dv[i])) {
            dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
          } else {
            dvt[i] =NA_REAL;
            normRelated[i] = 0;
          }
        } else {
          // (-Inf, dv)
          double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
          double sd = sqrt(ri[i]);
          interestingLim=true;
          lowerLim[i] = R_NegInf;
          upperLim[i] = dv[i];
          if (doSim && (censMethod == CENS_TNORM || censMethod == CENS_IPRED || censMethod == CENS_PRED)) {
            if (censMethod == CENS_TNORM) {
              dvt[i] = truncnorm(ipredt[i], sd, R_NegInf, lim1TBS);
            } else if (censMethod == CENS_PRED) {
              dvt[i] = predt[i];
            } else {
              dvt[i] = ipredt[i];
            }
            if (!R_finite(ipredt[i])) {
              normRelated[i] = 0;
            }
            if (R_finite(dvt[i])) {
              dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
            } else {
              dv[i] = NA_REAL;
              normRelated[i] = 0;
            }
          } else {
            dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
          }
        }
        break;
      case -1:
        // (dv, limit); limit could be +inf
        if (R_FINITE(limit[i])){
          //(dv, limit)
          double lim1TBS = _powerD(limit[i], lambda[i], (int)yj[i], low[i], hi[i]);
          double lim0TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
          double sd = sqrt(ri[i]);
          interestingLim=true;
          lowerLim[i] = dv[i];
          upperLim[i] = limit[i];
          if (!R_finite(ipredt[i])) {
            normRelated[i] = 0;
          }
          if (doSim && (censMethod == CENS_TNORM || censMethod == CENS_IPRED || censMethod == CENS_PRED)) {
            if (censMethod == CENS_TNORM) {
              dvt[i] = truncnorm(ipredt[i], sd, lim0TBS, lim1TBS);
            } else if (censMethod == CENS_PRED) {
              dvt[i] = predt[i];
            } else {
              dvt[i] = ipredt[i];
            }
            if (R_finite(dvt[i])) {
              dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
            } else {
              normRelated[i] = 0;
            }
          } else if (R_finite(dv[i])) {
            dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
          } else {
            dvt[i] = NA_REAL;
            normRelated[i] = 0;
          }
        } else {
          // (dv, Inf)
          double lim1TBS = _powerD(dv[i], lambda[i], (int) yj[i], low[i], hi[i]);
          double sd = sqrt(ri[i]);
          interestingLim=true;
          lowerLim[i] = dv[i];
          upperLim[i] = R_PosInf;
          if (!R_finite(ipredt[i])) {
            normRelated[i] = 0;
          }
          if (doSim && (censMethod == CENS_TNORM || censMethod == CENS_IPRED || censMethod == CENS_PRED)) {
            if (censMethod == CENS_TNORM) {
              dvt[i] = truncnorm(ipredt[i], sd, lim1TBS, R_PosInf);
            } else if (censMethod == CENS_PRED) {
              dvt[i] = predt[i];
            } else {
              dvt[i] = ipredt[i];
            }
            if (R_finite(dvt[i])) {
              dv[i] = _powerDi(dvt[i], lambda[i], (int) yj[i], low[i], hi[i]);
            } else {
              dv[i] = NA_REAL;
              normRelated[i] = 0;
            }
          } else if (R_finite(dv[i])) {
            dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
          } else {
            dvt[i] = NA_REAL;
            normRelated[i] = 0;
          }
        }
        break;
      case 0:
        lowerLim[i] = NA_REAL;
        upperLim[i] = NA_REAL;
        if (R_finite(dv[i])) {
          dvt[i] = _powerD(dv[i], lambda[i], (int)yj[i], low[i], hi[i]);
        } else {
          dvt[i] = NA_REAL;
          normRelated[i] = 0;
        }
        break;
      }
    } else if (dist == rxDistributionT ||
               dist == rxDistributionCauchy) {
      ipred[i]= _powerDi(ipredt[i], lambda[i], (int)yj[i], low[i], hi[i]);
      pred[i]= _powerDi(predt[i], lambda[i], (int)yj[i], low[i], hi[i]);
      if (R_finite(ipredt[i])) {
              normRelated[i] = 1;
      } else {
        normRelated[i] = 0;
      }
    } else {
      ipred[i]= ipredt[i];
      pred[i]= predt[i];
      normRelated[i] = 0;
    }
  }
  normIdx = find(normRelated == 1);
  nonNormIdx = find(normRelated == 0);
  return interestingLim;
}
