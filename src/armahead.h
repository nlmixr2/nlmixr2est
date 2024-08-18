#ifndef __ARMAHEAD_H__
#define __ARMAHEAD_H__
#define ARMA_WARN_LEVEL 1
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <rxode2ptr.h>


#ifndef rxDistributionNorm
// backward compatability for new rxode2
#define rxDistributionNorm     1
#define rxDistributionPois     2
#define rxDistributionBinom    3
#define rxDistributionBeta     4
#define rxDistributionT        5
#define rxDistributionChisq    6
#define rxDistributionDexp     7
#define rxDistributionF        8
#define rxDistributionGeom     9
#define rxDistributionHyper   10
#define rxDistributionUnif    11
#define rxDistributionWeibull 12
#define rxDistributionCauchy  13
#define rxDistributionGamma   14
#define rxDistributionOrdinal 15
#define rxDistributionN2ll    16
#define rxDistributionDnorm   17
#define rxHasFoceiLlik false
static inline void _splitYj(int *yj, int *dist,  int *trans) {
  *dist  = 1;
}
#else
#define rxHasFoceiLlik true
#endif

arma::mat cholSE__(arma::mat A, double tol);
bool cholSE0(arma::mat &Ao, arma::mat &E, arma::mat A, double tol);

using namespace arma;
using namespace Rcpp;

#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))


#endif
