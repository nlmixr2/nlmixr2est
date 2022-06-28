#define STRICT_R_HEADER
#include "censEst.h"

int globalCensFlag = 0;

SEXP censEstGetFactor() {
  Rcpp::IntegerVector ret(1);
  ret[0] = globalCensFlag + 1;
  ret.attr("class") = "factor";
  // M2 = 1
  // M3 = 2
  // M4 = 4
  ret.attr("levels") = Rcpp::CharacterVector::create("No censoring", //1
                                  "M2 censoring", //2
                                  "M3 censoring", //3
                                  "M2 and M3 censoring", //4
                                  "M4 censoring", // 5
                                  "M2 and M4 censoring", // 6
                                  "M3 and M4 censoring", // 7
                                  "M2, M3 and M4 censoring" // 8
                                  );
  return Rcpp::wrap(ret);
}
