#ifndef __CENSRESID_H__
#define __CENSRESID_H__
#include "armahead.h"
#define CENS_OMIT 1
#define CENS_CDF 2
#define CENS_TNORM 3
#define CENS_IPRED 4
#define CENS_PRED 5
#define CENS_EPRED 6

#if defined(__cplusplus)

extern "C" {
  typedef SEXP (*_rxode2random_rxRmvnSEXP_t)(SEXP nSSEXP, SEXP muSSEXP, SEXP sigmaSSEXP, SEXP lowerSSEXP, SEXP upperSSEXP, SEXP ncoresSSEXP, SEXP isCholSSEXP, SEXP keepNamesSSEXP, SEXP aSSEXP, SEXP tolSSEXP, SEXP nlTolSSEXP, SEXP nlMaxiterSSEXP);
  extern _rxode2random_rxRmvnSEXP_t rxRmvnSEXPnlmixrEst;
}

static inline double truncnorm(double mean, double sd, double low, double hi){
  if (R_finite(mean) && R_finite(sd)) {
    NumericMatrix sigma(1,1);
    sigma(0,0)=sd;
    SEXP ret = rxRmvnSEXPnlmixrEst(wrap(IntegerVector::create(1)),
                                       wrap(NumericVector::create(mean)),
                                       wrap(sigma),
                                       wrap(NumericVector::create(low)),
                                       wrap(NumericVector::create(hi)),
                                       wrap(IntegerVector::create(1)),
                                       wrap(LogicalVector::create(false)),
                                       wrap(LogicalVector::create(false)),
                                       wrap(NumericVector::create(0.4)),
                                       wrap(NumericVector::create(2.05)),
                                       wrap(NumericVector::create(1e-10)),
                                       wrap(IntegerVector::create(100)));
    return REAL(ret)[0];
  }
  return NA_REAL;
}

bool censTruncatedMvnReturnInterestingLimits(arma::vec& dv, arma::vec& dvt,
                                             arma::vec& ipred, arma::vec &ipredt,
                                             arma::vec& pred, arma::vec &predt,
                                             arma::Col<int> &cens, arma::vec &limit,
                                             arma::vec& lambda, arma::vec &yj, arma::vec& low, arma::vec& hi,
                                             arma::vec &lowerLim, arma::vec &upperLim, arma::vec &ri,
                                             bool &doSim, int& censMethod,
                                             arma::uvec &normRelated, arma::uvec &normIdx,
                                             arma::uvec &nonNormIdx);

#endif

#endif
