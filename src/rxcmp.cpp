// [[Rcpp::plugins(openmp)]]
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"

int _rxode2ver4 = 0;
extern "C" SEXP _rxode2version4(SEXP v4) {
  // This is a work-round to allow rxode2 3 and rxode 4 to both work
  _rxode2ver4 = INTEGER(v4)[0];
  return R_NilValue;
}

//[[Rcpp::export]]
SEXP rxode2stateOde(SEXP inp) {
  // This is a work-round to allow rxode2 3 and rxode 4 to both work
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment rxode2 = loadNamespace("rxode2");
  Function stateOde = rxode2["rxState"];
  if (_rxode2ver4) {
    stateOde = rxode2["rxStateOde"];
  }
  return stateOde(inp);
}


extern "C" SEXP _rxode2rxFixRes(SEXP inp, SEXP n) {
  Function loadNamespace("loadNamespace", R_BaseNamespace);
  Environment rxode2 = loadNamespace("rxode2");
  if (_rxode2ver4) {
    Function rxFixRes = rxode2["rxFixRes"];
    return rxFixRes(inp, n);
  }
  return inp;
}
