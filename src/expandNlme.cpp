// -*- mode: c++; c-basic-offset: 2; tab-width: 2; indent-tabs-mode: t; -*-
#define USE_FC_LEN_T
#define STRICT_R_HEADERS
// [[Rcpp::interfaces(r, cpp)]]
//#undef NDEBUG
#include <Rcpp.h>
using namespace Rcpp;

#define _(String) (String)

std::string symengineRes(std::string val){
  if (val == "e" ||
      val == "E" ||
      val == "EulerGamma" ||
      val == "Catalan" ||
      val == "GoldenRatio" ||
      val == "I"){
    return "rx_SymPy_Res_" + val;
  }
  return val;
}


//' Expand Gradient for nlme
//'
//' @param state is the state to expand
//' @param params is the parameters to expand
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
List nlmixrExpandFdParNlme_(CharacterVector state, CharacterVector vars) {
  int neta = vars.size();
  CharacterVector fe(neta);
  CharacterVector calcS(neta);
  int nstate = state.size();
  for (int i = 0; i < neta; i++){
    std::string etaN = as<std::string>(vars[i]);
    std::string feta;
    std::string calc;
    feta = "rxD_" + etaN;
    calc = "assign(\"" + feta + "\",with(.s,-D(rx_pred_, " +
      etaN + ")";
    for (int j = nstate; j--;){
      calc += "-rx__sens_" + as<std::string>(state[j]) + "_BY_" + etaN +
        "__*D(rx_pred_,\""+ symengineRes(as<std::string>(state[j])) + "\")";
    }
    calc += "), envir=.s)";
    fe[i] = feta;
    calcS[i] = calc;
  }
  List out(2);
  out[0] = fe;
  out[1] = calcS;
  out.attr("names") = CharacterVector::create("dfe","calc");
  out.attr("class") = "data.frame";
  out.attr("row.names") = IntegerVector::create(NA_INTEGER, -neta);
  return out;
}
