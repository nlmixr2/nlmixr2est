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
#include <rxode2.h>



arma::mat cholSE__(arma::mat A, double tol);
bool cholSE0(arma::mat &Ao, arma::mat &E, arma::mat A, double tol);

using namespace arma;
using namespace Rcpp;

#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))


#endif
