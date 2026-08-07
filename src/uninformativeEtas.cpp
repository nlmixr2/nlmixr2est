#define STRICT_R_HEADER
#include "armahead.h"


using namespace Rcpp;
using namespace arma;

extern "C" SEXP _nlmixr2est_uninformativeEta(SEXP rhoS) {
  Environment rho = as<Environment>(rhoS);
  int nid = as<int>(rho["nid"]);
  int neta = as<int>(rho["neta"]);
  double tol = as<double>(rho["tol"]);
  IntegerVector simId = as<IntegerVector>(rho["simId"]);
  IntegerVector id = as<IntegerVector>(rho["id"]);
  IntegerVector pm = as<IntegerVector>(rho["pm"]);
  IntegerVector w = as<IntegerVector>(rho["w"]);
  NumericVector val = as<NumericVector>(rho["val"]);
  NumericMatrix ret=NumericMatrix(nid, neta);
  std::fill(ret.begin(), ret.end(), 0.0);
  // size of the predictions `ret` is built from, to tell "the eta does not move the
  // prediction" apart from "there is no prediction to move".  Largest, not a sum: a sum
  // grows with the observation count, so the same underflowing subject would change
  // verdict on nothing but how often it happened to be sampled.  NaN never wins the
  // comparison, so a failed solve leaves this at zero.
  NumericMatrix scale=NumericMatrix(nid, neta);
  std::fill(scale.begin(), scale.end(), 0.0);
  IntegerMatrix retL=IntegerMatrix(nid, neta);
  for (int i = 0; i < simId.size(); i++) {
    int curid = id[i] - 1;
    int curpm = pm[i];
    int cureta = w[i] - 1;
    double curAbs = std::abs(val[i]);
    if (curAbs > scale(curid, cureta)) scale(curid, cureta) = curAbs;
    switch (curpm) {
    case -1:
    case 1:
      ret(curid, cureta) += val[i];
      break;
    case 0:
      ret(curid, cureta) -= 2.0*val[i];
    }
  }
  for (int i = 0; i < nid; ++i) {
    for (int j = 0; j < neta; ++j) {
      // This runs once, at the initial estimates.  If those are poor enough that the
      // prediction underflows (or fails outright) at every observed time, perturbing the
      // eta moves nothing for a reason unrelated to the eta -- do not freeze it for the
      // whole fit on that.  `ret` is checked too: a solve that failed on ONE arm leaves it
      // NaN while the other arms keep `scale` finite and over tol, and NaN > tol is false.
      bool havePred = R_finite(ret(i, j)) && R_finite(scale(i, j)) && scale(i, j) > tol;
      retL(i, j) = (std::abs(ret(i, j)) > tol) || !havePred;
    }
  }
  return as<SEXP>(retL);
}
