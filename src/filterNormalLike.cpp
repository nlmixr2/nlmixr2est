//#undef NDEBUG
#ifndef NDEBUG
#define NDEBUG // just in case
#endif
#define USE_FC_LEN_T
#define STRICT_R_HEADERS
#ifndef NDEBUG
#define NDEBUG // just in case
#endif
#include "armahead.h"
using namespace Rcpp;


bool normalLike(int curCmt, IntegerVector &inDistribution, IntegerVector& inDistCmt,
                int &nnorm, int &nlik, int &nother) {
  for (unsigned int i = inDistCmt.size(); i--;) {
    if (curCmt == inDistCmt[i]) {
      int dist = inDistribution[i];
      if (dist == rxDistributionNorm ||
          dist == rxDistributionDnorm ||
          dist == rxDistributionT ||
          dist == rxDistributionCauchy) {
        nnorm++;
        return true;
      } else {
        nlik++;
        return false;
      }
    }
  }
  nother++;
  return true;
}

//[[Rcpp::export]]
List filterNormalLikeAndDoses(IntegerVector& inCmt, IntegerVector& inDistribution,
                              IntegerVector& inDistCmt) {
  int nnorm = 0;
  int nlik = 0;
  int nother = 0;
  LogicalVector filter(inCmt.size());
  for (unsigned int i = filter.size(); i--;) {
    filter[i] = normalLike(inCmt[i], inDistribution, inDistCmt, nnorm, nlik, nother);
  }
  return List::create(_["filter"]=filter,
                      _["nnorm"]=nnorm,
                      _["nlik"]=nlik,
                      _["nother"]=nother);
}


//[[Rcpp::export]]
LogicalVector rxode2hasLlik() {
  return LogicalVector::create(rxHasFoceiLlik);
}
