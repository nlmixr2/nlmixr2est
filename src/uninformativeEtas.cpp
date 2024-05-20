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
  IntegerMatrix retL=IntegerMatrix(nid, neta);
  for (int i = 0; i < simId.size(); i++) {
    int curid = id[i] - 1;
    int curpm = pm[i];
    int cureta = w[i] - 1;
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
      retL(i, j) = std::abs(ret(i, j)) > tol;
    }
  }
  return as<SEXP>(retL);
}
